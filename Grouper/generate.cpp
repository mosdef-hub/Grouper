#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>

#include <atomic>
#include <condition_variable>
#include <csignal>
#include <cstdlib>
#include <deque>
#include <exception>
#include <iostream>
#include <fstream>
#include <mutex>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <random>
#include <algorithm>

#include <libpq-fe.h>
#include <omp.h>
#include <stdio.h>

#include <pybind11/pybind11.h>

#include <nauty/nauty.h>

namespace py = pybind11;

#include "dataStructures.hpp"
#include "processColoredGraphs.hpp"
#include "generate.hpp"

#define MAX_EDGES 100

// Function to update and display progress bar
void update_progress(int current, int total) {
    float progress = static_cast<float>(current) / total;
    int bar_width = 100;

    std::cout << "[";
    int pos = bar_width * progress;
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << current << "/" << total << " (" << int(progress * 100.0) << "%)\r";
    std::cout.flush();
}

struct hash_vector {
    std::size_t operator()(const std::vector<setword>& v) const {
        std::size_t seed = 0;
        for (const auto& i : v) {
            boost::hash_combine(seed, i);
        }
        return seed;
    }
};

// Process-wide pointer used by the SIGINT handler to flip the per-run
// `interrupted` flag. We can't rely on Python's PyErr_CheckSignals here:
// with the GIL released, signal delivery is racy across the OMP master
// vs. our std::thread producer, and on macOS the signal often lands on
// a non-Python-aware OMP worker thread that doesn't wake up the main
// thread's flag. A direct C-level handler fires regardless of which
// thread receives the signal and just sets our atomic.
static std::atomic<std::atomic<bool>*> g_interrupt_flag{nullptr};
extern "C" void grouper_sigint_handler(int) {
    auto* p = g_interrupt_flag.load(std::memory_order_acquire);
    if (p) p->store(true, std::memory_order_release);
}

// RAII guard: install our SIGINT handler on entry, restore the previous
// (typically Python's) handler on exit. Saving and restoring is what
// keeps `KeyboardInterrupt` working in the surrounding Python session
// after the call returns.
class SigIntScope {
public:
    explicit SigIntScope(std::atomic<bool>& flag) {
        g_interrupt_flag.store(&flag, std::memory_order_release);
        struct sigaction new_act{};
        new_act.sa_handler = grouper_sigint_handler;
        sigemptyset(&new_act.sa_mask);
        new_act.sa_flags = 0; // no SA_RESTART — let blocking syscalls
                              // (e.g. fgets) return EINTR so the
                              // producer can notice and exit.
        ::sigaction(SIGINT, &new_act, &old_act_);
    }
    ~SigIntScope() {
        ::sigaction(SIGINT, &old_act_, nullptr);
        g_interrupt_flag.store(nullptr, std::memory_order_release);
    }
    SigIntScope(const SigIntScope&) = delete;
    SigIntScope& operator=(const SigIntScope&) = delete;
private:
    struct sigaction old_act_{};
};

// Bounded thread-safe queue used to stream vcolg output lines from the
// reader thread into the OMP worker pool. The cap exists so the producer
// blocks if consumers can't keep up — at large N (n=7+) the full vcolg
// output can run to GB; without the cap we'd buffer it all in memory.
namespace {
class BoundedLineQueue {
public:
    explicit BoundedLineQueue(std::size_t capacity) : capacity_(capacity) {}

    // Producer side. Blocks if the queue is full.
    void push(std::string line) {
        std::unique_lock<std::mutex> lock(mtx_);
        not_full_.wait(lock, [&] { return q_.size() < capacity_ || closed_; });
        if (closed_) return; // dropped — consumers gave up
        q_.push_back(std::move(line));
        not_empty_.notify_one();
    }

    // Consumer side. Returns false iff the queue is closed AND empty.
    bool pop(std::string& out) {
        std::unique_lock<std::mutex> lock(mtx_);
        not_empty_.wait(lock, [&] { return !q_.empty() || closed_; });
        if (q_.empty()) return false;
        out = std::move(q_.front());
        q_.pop_front();
        not_full_.notify_one();
        return true;
    }

    // Producer signals end of stream. Wakes any blocked consumers so
    // they can exit their pop loops.
    void close() {
        std::lock_guard<std::mutex> lock(mtx_);
        closed_ = true;
        not_empty_.notify_all();
        not_full_.notify_all();
    }

private:
    std::size_t capacity_;
    std::deque<std::string> q_;
    std::mutex mtx_;
    std::condition_variable not_empty_;
    std::condition_variable not_full_;
    bool closed_ = false;
};
} // anonymous namespace

std::unordered_map<std::string, std::string> parseConfig(const std::string& configFile) {
    std::unordered_map<std::string, std::string> configParams;
    std::ifstream config(configFile);
    std::string line;
    while (std::getline(config, line)) {
        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            configParams[key] = value;
        }
    }
    return configParams;
}

void executeQuery(PGconn* conn, const std::string& query, const std::vector<std::string>& params = {}) {
    const char* paramValues[params.size()];
    for (size_t i = 0; i < params.size(); ++i) {
        paramValues[i] = params[i].c_str();
    }
    PGresult* res = PQexecParams(conn, query.c_str(), params.size(), nullptr, paramValues, nullptr, nullptr, 0);
    if (PQresultStatus(res) != PGRES_COMMAND_OK) {
        std::cerr << "Query execution failed: " << PQerrorMessage(conn) << std::endl;
    }
    PQclear(res);
}

std::tuple< int, std::vector<int>, std::vector<std::pair<int, int>> > parse_nauty_graph_line_tmp(const std::string& line, const std::unordered_set<GroupGraph::Group>& node_defs) {
    // Split the line into node_description and edge_description
    size_t split_pos = line.find("  ");
    if (split_pos == std::string::npos) {
        throw std::runtime_error("Invalid nauty output line...");
    }

    std::string node_description_str = line.substr(0, split_pos);
    std::string edge_description_str = line.substr(split_pos + 2);

    std::istringstream node_desc_iss(node_description_str);
    std::vector<std::string> node_description(std::istream_iterator<std::string>{node_desc_iss},
                                              std::istream_iterator<std::string>());

    std::istringstream edge_desc_iss(edge_description_str);
    std::vector<std::string> edge_description(std::istream_iterator<std::string>{edge_desc_iss},
                                              std::istream_iterator<std::string>());

    std::vector<std::pair<int, int>> edge_list;
    for (size_t i = 0; i < edge_description.size(); i += 2) {
        edge_list.emplace_back(std::stoi(edge_description[i]), std::stoi(edge_description[i + 1]));
    }



    int n_vertices = std::stoi(node_description[0]);
    // int n_edges = std::stoi(node_description[1]);
    std::vector<int> colors;
    for (size_t i = 2; i < node_description.size(); ++i) {
        colors.push_back(std::stoi(node_description[i]));
    }

    if (node_defs.size() < static_cast<size_t>(*std::max_element(colors.begin(), colors.end()) + 1)) {
        throw std::runtime_error("Number of nodes in node_defs does not match the number of nodes in the nauty_output_file...");
    }

    return {n_vertices, colors, edge_list};
}

bool check_max_bond_not_exceeded_tmp(
    const std::vector<std::pair<int,int>>& edge_list,
    const std::vector<int>& colors,
    const std::unordered_map<std::string, std::vector<int>>& node_types,
    const std::unordered_map<int, std::string>& int_to_node_type) {

    std::vector<int> node_bond_count;
    for (std::size_t i = 0; i < colors.size(); ++i) {
        std::string node_type = int_to_node_type.at(colors[i]);
        int n_ports = node_types.at(node_type).size();
        node_bond_count.push_back(n_ports);
    }
    for (const auto& edge : edge_list) {
        int src = edge.first;
        int dst = edge.second;
        node_bond_count[src] -= 1;
        node_bond_count[dst] -= 1;
    }
    for (const auto& count : node_bond_count) {
        if (count < 0) {
            return false;
        }
    }
    return true;
}




void insertGraph(PGconn* conn, const GroupGraph& graph, const std::string& table_name) {
    std::string smiles = graph.toSmiles();
    std::string graph_data = graph.serialize();
    // Bind to a named local so the c_str() pointer remains valid through
    // PQexecParams. The previous std::to_string(n_nodes).c_str() returned
    // a pointer into a temporary destroyed at the semicolon, so the int
    // column was being read from freed memory.
    std::string n_nodes_str = std::to_string(graph.nodes.size());

    const char* paramValues[3];
    paramValues[0] = smiles.c_str();
    paramValues[1] = graph_data.c_str();
    paramValues[2] = n_nodes_str.c_str();

    std::string insert_query = "INSERT INTO " + table_name + " (smiles, graph_data, n_nodes) VALUES ($1, $2, $3) ON CONFLICT (smiles) DO NOTHING";
    PGresult* insertRes = PQexecParams(conn, insert_query.c_str(), 3, nullptr, paramValues, nullptr, nullptr, 0);

    if (PQresultStatus(insertRes) != PGRES_COMMAND_OK && PQresultStatus(insertRes) != PGRES_TUPLES_OK) {
    }
    PQclear(insertRes);
}

std::unordered_set<GroupGraph> exhaustiveGenerate(
    int n_nodes,
    std::unordered_set<GroupGraph::Group> node_defs,
    int num_procs = -1,
    std::string vcolg_output_file = "",
    std::unordered_map<std::string, int> positiveConstraints = {},
    std::unordered_set<std::string> negativeConstraints = {},
    std::string config_path = ""
) {


    // Error handling
    if (n_nodes < 1) {
        throw std::invalid_argument("Number of nodes must be greater than 0...");
    }
    if (node_defs.size() < 1) {
        throw std::invalid_argument("Group definitions must not be empty...");
    }
    if (num_procs <= -1){
        num_procs = omp_get_max_threads();
    }
    if (!positiveConstraints.empty()){
        for (const auto& constraint : positiveConstraints) {
            if (constraint.second < 0) {
                throw std::invalid_argument("Positive constraint value must be greater than or equal to 0...");
            }
        }
    }


    // Stream vcolg output through a bounded queue. A producer thread
    // reads lines from the popen'd `geng | vcolg` pipe (or the
    // user-provided file) and pushes them into the queue; the OMP
    // workers below consume from it. Memory is capped at the queue's
    // capacity, so even gigabyte-scale vcolg runs at large N don't
    // buffer their full output in RAM.
    BoundedLineQueue queue(/*capacity=*/8192);
    std::atomic<long long> producer_total{0};
    // Set to true by our SIGINT handler (installed below). Producer
    // and workers check this on each loop iteration and exit early;
    // after the OMP region we re-raise as KeyboardInterrupt so the
    // run doesn't silently complete after Ctrl-C.
    std::atomic<bool> interrupted{false};
    SigIntScope sigint_scope(interrupted);
    // vcolg/geng diagnostic lines (the >A header / >Z trailer lines)
    // captured via shell stderr redirection to a temp file, then
    // replayed to stderr after the bar finishes so they don't race
    // the bar on the user's terminal.
    std::vector<std::string> diag_lines;
    std::exception_ptr producer_exc = nullptr;
    const bool from_pipe = vcolg_output_file.empty();
    const std::string vcolg_input_path = vcolg_output_file;

    // Pre-count vcolg's output so the progress bar has a stable
    // denominator. The producer's queue is bounded (capacity 8192),
    // so under load it pushes at the *workers'* rate — meaning
    // producer_total tracks consumption, not production, and any
    // bar built on it would only show a meaningful percentage in
    // the last few % of the run. A separate `geng | vcolg | wc -l`
    // pass is single-threaded and finishes long before the workers
    // do, so the cost is small relative to the OMP processing
    // phase. For the file-input branch we count lines in the file.
    long long total_lines = 0;
    if (from_pipe) {
        std::string count_cmd = "geng " + std::to_string(n_nodes) +
                                " -c 2>/dev/null | vcolg -T -m" +
                                std::to_string(node_defs.size()) +
                                " 2>/dev/null | wc -l";
        FILE* cp = popen(count_cmd.c_str(), "r");
        if (cp) {
            if (std::fscanf(cp, "%lld", &total_lines) != 1) total_lines = 0;
            pclose(cp);
        }
    } else {
        std::ifstream count_in(vcolg_input_path);
        if (count_in.is_open()) {
            std::string l;
            while (std::getline(count_in, l)) ++total_lines;
        }
    }

    std::thread producer([&]() {
        try {
            auto push_line = [&](std::string&& line) {
                if (line.empty()) return;
                queue.push(std::move(line));
                producer_total.fetch_add(1, std::memory_order_relaxed);
            };
            if (from_pipe) {
                // Capture nauty's stderr (>A/>Z lines) to a temp file
                // rather than 2>&1-merging it: the merge is byte-level,
                // not line-level, so a partial stdout line + a stderr
                // write can concatenate into one corrupted line. With a
                // separate fd, line atomicity is preserved and the
                // diagnostics appear cleanly under the bar at the end.
                char err_template[] = "/tmp/grouper_nauty_err_XXXXXX";
                int err_fd = ::mkstemp(err_template);
                if (err_fd < 0) {
                    throw std::runtime_error("mkstemp failed for nauty stderr capture");
                }
                ::close(err_fd);
                const std::string err_path = err_template;
                std::string cmd = "(geng " + std::to_string(n_nodes) + " -c | vcolg -T -m" +
                                  std::to_string(node_defs.size()) + ") 2>" + err_path;
                FILE* pipe = popen(cmd.c_str(), "r");
                if (!pipe) {
                    ::unlink(err_path.c_str());
                    throw std::runtime_error("popen failed for: " + cmd);
                }
                char buf[4096];
                std::string accumulated;
                while (std::fgets(buf, sizeof(buf), pipe) != nullptr) {
                    if (interrupted.load(std::memory_order_acquire)) break;
                    accumulated.append(buf);
                    if (!accumulated.empty() && accumulated.back() == '\n') {
                        accumulated.pop_back();
                        push_line(std::move(accumulated));
                        accumulated.clear();
                    }
                }
                if (!accumulated.empty()) push_line(std::move(accumulated));
                int rc = pclose(pipe);
                {
                    std::ifstream err_in(err_path);
                    std::string err_line;
                    while (std::getline(err_in, err_line)) {
                        diag_lines.push_back(std::move(err_line));
                    }
                }
                ::unlink(err_path.c_str());
                if (rc != 0 && !interrupted.load(std::memory_order_acquire)) {
                    throw std::runtime_error(
                        "geng | vcolg failed (exit " + std::to_string(rc) +
                        ", is nauty on PATH? command: " + cmd + ")");
                }
            } else {
                std::ifstream input_file(vcolg_input_path);
                if (!input_file.is_open()) {
                    throw std::runtime_error("Error opening vcolg output file: " + vcolg_input_path);
                }
                std::string line;
                while (std::getline(input_file, line)) {
                    if (interrupted.load(std::memory_order_acquire)) break;
                    push_line(std::move(line));
                }
            }
        } catch (...) {
            producer_exc = std::current_exception();
        }
        queue.close();
    });
    // Make sure the producer is reaped on every exit path.
    struct ProducerJoiner {
        std::thread& t;
        ~ProducerJoiner() { if (t.joinable()) t.join(); }
    } joiner{producer};

    std::unordered_set<GroupGraph> global_basis;
    omp_set_num_threads(num_procs);
    std::atomic<long long> n_finished{0};
    // Progress bar with streaming: producer_total is the running count of
    // lines pushed by the producer. It races ahead of n_finished early
    // and stabilises at the true total once the producer is done — by
    // then the percentage is exact. Update every progress_step leaves to
    // keep the omp critical section out of the hot loop.
    constexpr long long progress_step = 256;
    std::cout << "Using " << num_procs << " processors (streaming)" << std::endl;
    // Detect whether stdout is an interactive terminal. The TTY path
    // uses ANSI escapes (\033[2K\r) for an in-place updating bar; the
    // non-TTY path (CI logs, piped output, redirected to a file) emits
    // a plain "N% (a/b)" line at every 10% milestone, since dumping
    // the escape codes and the ~750 progress prints of a long run
    // produces unreadable log files.
    const bool stdout_is_tty = ::isatty(::fileno(stdout));
    int last_pct_printed = -1;
    auto print_progress = [&](long long finished) {
        if (total_lines <= 0) {
            // Pre-count failed or wasn't run; show count only. In
            // non-TTY mode skip — without a denominator the count is
            // not very informative and would spam logs.
            if (stdout_is_tty) {
                std::cout << "\033[2K\rprocessed " << finished << " lines" << std::flush;
            }
            return;
        }
        if (finished > total_lines) finished = total_lines;
        int pct = static_cast<int>((100 * finished) / total_lines);
        if (!stdout_is_tty) {
            int rounded = (pct / 10) * 10;
            if (rounded == 0 || rounded <= last_pct_printed) return;
            last_pct_printed = rounded;
            std::cout << rounded << "% (" << finished << "/" << total_lines
                      << ")\n" << std::flush;
            return;
        }
        // \033[2K clears the current terminal line; pairs with the \r
        // before each write so we always start from a known-empty
        // line, even if a stray write to stdout/stderr (or a shorter
        // previous bar) left residue.
        constexpr int bar_width = 100;
        int pos = static_cast<int>((bar_width * finished) / total_lines);
        std::cout << "\033[2K\r[";
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << finished << "/" << total_lines << " (" << pct << "%)" << std::flush;
    };

    if (!config_path.empty()) {
        std::unordered_map<std::string, std::string> configParams = parseConfig(config_path);
        std::string table_name = configParams["table_name"];
        std::string conninfo = "dbname=" + configParams["dbname"] +
                       " user=" + configParams["user"] +
                       " password=" + configParams["password"] +
                       " hostaddr=" + configParams["hostaddr"] +
                       " port=" + configParams["port"];

        PGconn* main_conn = PQconnectdb(conninfo.c_str());
        if (PQstatus(main_conn) != CONNECTION_OK) {
            std::cerr << "Connection to database failed: " << PQerrorMessage(main_conn) << std::endl;
            PQfinish(main_conn);
            throw std::runtime_error("Database connection failed.");
        }
        std::string create_table_query =
            "CREATE TABLE IF NOT EXISTS " + table_name + " (" +
            "smiles TEXT PRIMARY KEY, " +
            "graph_data TEXT, " +
            "n_nodes INT" +
            ");";
        executeQuery(main_conn, create_table_query);
        PQfinish(main_conn);

        #pragma omp parallel
        {
            PGconn* conn = PQconnectdb(conninfo.c_str());
            if (PQstatus(conn) != CONNECTION_OK) {
                #pragma omp critical
                {
                    std::cerr << "Thread " << omp_get_thread_num() << " failed to connect to DB: " << PQerrorMessage(conn) << std::endl;
                }
                // libpq requires PQfinish even on failed connection
                // objects to release the conn struct's memory.
                PQfinish(conn);
            } else {
                LocalBasis local_basis;
                // n is the maximum number of nauty vertices we'll see in
                // this run — exactly n_nodes, since geng emits graphs of
                // size n_nodes. The previous hard-coded 20 silently
                // overflowed the preallocated buffers when callers asked
                // for n_nodes > 20.
                const int n = n_nodes;
                int m = SETWORDSNEEDED(n);
                std::vector<setword> g(m * n, 0);
                std::vector<int> lab(n, 0), ptn(n, 0), orbits(n, 0);
                DEFAULTOPTIONS_GRAPH(options);
                statsblk stats;

                std::string line;
                while (queue.pop(line)) {
                    if (interrupted.load(std::memory_order_acquire)) {
                        queue.close();
                        break;
                    }
                    local_basis.clear();
                    process_nauty_output(
                        line, node_defs, &local_basis,
                        positiveConstraints, negativeConstraints,
                        g.data(), lab.data(), ptn.data(), orbits.data(), &options, &stats
                    );
                    for (const auto& [canon, graph] : local_basis) {
                        insertGraph(conn, graph, table_name);
                    }
                    long long finished = ++n_finished;
                    if (finished % progress_step == 0) {
                        #pragma omp critical
                        { print_progress(finished); }
                    }
                }
                PQfinish(conn);
            }
        }
        // Producer is done at this point (queue is closed and drained).
        // Join explicitly so we can surface its exception, if any. The
        // joiner declared above will see joinable()==false and no-op.
        producer.join();
        if (producer_exc) std::rethrow_exception(producer_exc);
        print_progress(n_finished.load());
        // TTY bar leaves the cursor on the bar line with no \n; force
        // one so subsequent output (the diagnostics replay, "Number
        // of unique graphs", or the user's own prints) starts on a
        // clean row. In non-TTY mode the last milestone already ended
        // with \n, so an unconditional endl would just add a blank.
        if (stdout_is_tty) std::cout << std::endl;
        for (const auto& l : diag_lines) std::cerr << l << '\n';
        if (interrupted.load(std::memory_order_acquire)) {
            // Our SIGINT handler set the flag. Raise it to Python as a
            // KeyboardInterrupt so the user gets the standard interrupt
            // experience instead of a silent return after Ctrl-C.
            py::gil_scoped_acquire gil;
            PyErr_SetString(PyExc_KeyboardInterrupt, "");
            throw py::error_already_set();
        }
        std::cout << "Graphs saved to database." << std::endl;
        return {};
    } else {
        size_t previous_n_colorings = 0;
        int max_threads = num_procs;
        // Each thread keeps a canon-keyed map: same molecule produced by
        // different vcolg lines collapses into one entry, instead of
        // accumulating one duplicate copy per line. Memory now scales
        // with the unique-output count rather than total work performed.
        std::vector<LocalBasis> all_local_bases(max_threads);

    #pragma omp parallel
    {
            int tid = omp_get_thread_num();
            LocalBasis& local_basis = all_local_bases[tid];
            // n is the maximum number of nauty vertices we'll see —
            // matches geng's output size exactly.
            const int n = n_nodes;
            int m = SETWORDSNEEDED(n);
            std::vector<setword> g(m * n, 0);
            std::vector<int> lab(n, 0), ptn(n, 0), orbits(n, 0);
            DEFAULTOPTIONS_GRAPH(options);
            statsblk stats;

            std::string line;
            while (queue.pop(line)) {
                if (interrupted.load(std::memory_order_acquire)) {
                    queue.close();
                    break;
                }
                process_nauty_output(
                    line, node_defs, &local_basis,
                    positiveConstraints, negativeConstraints,
                    g.data(), lab.data(), ptn.data(), orbits.data(), &options, &stats
                );
                long long finished = ++n_finished;
                if (finished % progress_step == 0) {
                    #pragma omp critical
                    { print_progress(finished); }
                }
            }
        }

        // Producer must have finished (queue closed/drained) for the
        // OMP region above to have exited. Surface any producer error.
        producer.join();
        if (producer_exc) std::rethrow_exception(producer_exc);
        // Final 100% bar update so the visible state lands at N/N, not
        // the last 256-boundary it happened to cross.
        print_progress(n_finished.load());
        // TTY bar leaves the cursor on the bar line with no \n; force
        // one so subsequent output (the diagnostics replay, "Number
        // of unique graphs", or the user's own prints) starts on a
        // clean row. In non-TTY mode the last milestone already ended
        // with \n, so an unconditional endl would just add a blank.
        if (stdout_is_tty) std::cout << std::endl;
        // Replay nauty's >A/>Z diagnostics now that the bar is settled,
        // so they appear cleanly under the bar instead of racing it.
        for (const auto& l : diag_lines) std::cerr << l << '\n';
        if (interrupted.load(std::memory_order_acquire)) {
            // Our SIGINT handler set the flag. Raise it to Python as a
            // KeyboardInterrupt so the user gets the standard interrupt
            // experience instead of a silent return after Ctrl-C.
            py::gil_scoped_acquire gil;
            PyErr_SetString(PyExc_KeyboardInterrupt, "");
            throw py::error_already_set();
        }

        // Cross-thread merge. Each thread already deduplicated within
        // itself by canon, so we just need a single global canon set
        // here. No re-canonicalization — the canon is the map key.
        std::unordered_set<std::vector<setword>, hash_vector> canon_basis;
        for (auto& local_basis : all_local_bases) {
            for (auto& [canon, graph] : local_basis) {
                if (canon_basis.insert(canon).second) {
                    global_basis.insert(std::move(graph));
                }
            }
        }
        // Free the per-thread maps (and their canon keys) before the
        // py::set conversion in the binding allocates 4519 wrappers.
        // Without this, container destruction would happen *after* the
        // function returns, contributing to the post-print delay.
        all_local_bases.clear();
        canon_basis.clear();

        std::cout<< "Number of unique graphs: " << global_basis.size() << std::endl;

        return global_basis;
    }
}

// This is the randomGenerate that utilizes the nauty library
std::unordered_set<GroupGraph> randomGenerate(
    int n_nodes,
    const std::unordered_set<GroupGraph::Group>& node_defs,
    int num_graphs = 100,
    int num_procs = -1,
    const std::unordered_map<std::string, int>& positiveConstraints = {},
    const std::unordered_set<std::string>& negativeConstraints = {}
) {
    if (n_nodes < 1) {
        throw std::invalid_argument("Number of nodes must be greater than 0...");
    }
    if (node_defs.empty()) {
        throw std::invalid_argument("Group definitions must not be empty...");
    }
    if (num_graphs < 1) {
        throw std::invalid_argument("Number of graphs must be greater than 0...");
    }
    for (const auto& constraint : positiveConstraints) {
        if (constraint.second < 0) {
            throw std::invalid_argument("Positive constraint value must be >= 0...");
        }
    }
    if (num_procs <= -1) {
        num_procs = omp_get_max_threads();
    }

    // Pipe geng directly into vcolg via popen, same as exhaustiveGenerate.
    std::vector<std::string> lines;
    {
        std::string cmd = "geng " + std::to_string(n_nodes) + " -c | vcolg -T -m" +
                          std::to_string(node_defs.size());
        FILE* pipe = popen(cmd.c_str(), "r");
        if (!pipe) {
            throw std::runtime_error("popen failed for: " + cmd);
        }
        char buf[4096];
        std::string accumulated;
        while (std::fgets(buf, sizeof(buf), pipe) != nullptr) {
            accumulated.append(buf);
            if (!accumulated.empty() && accumulated.back() == '\n') {
                accumulated.pop_back();
                if (!accumulated.empty()) lines.push_back(std::move(accumulated));
                accumulated.clear();
            }
        }
        if (!accumulated.empty()) lines.push_back(std::move(accumulated));
        int rc = pclose(pipe);
        if (rc != 0) {
            throw std::runtime_error(
                "geng | vcolg failed (exit " + std::to_string(rc) +
                ", is nauty on PATH? command: " + cmd + ")");
        }
    }
    if (lines.empty()) {
        throw std::runtime_error("No valid graphs found from geng | vcolg.");
    }

    std::vector<std::tuple<int, std::vector<int>, std::vector<std::pair<int, int>>>> possible_node_colored_graphs;

    std::unordered_set<GroupGraph> global_basis;
    std::unordered_set<std::string> global_smiles_basis;
    // std::random_device rd;

    // Create necessary maps
    std::unordered_map<int, std::string> int_to_node_type;
    std::unordered_map<int, GroupGraph::Group> color_to_group;
    std::unordered_map<std::string, std::vector<int>> node_types;
    std::unordered_map<std::string, std::vector<int>> node_type_to_hub;
    std::unordered_map<int, std::string> int_to_pattern;
    std::unordered_map<std::string, std::string> type_to_pattern;
    for (const auto& node : node_defs) {
        node_types[node.ntype] = node.ports;
        type_to_pattern[node.ntype] = node.pattern;
    }
    for (const auto& node: node_defs) {
        for (const auto& h : node.hubs) {
            node_type_to_hub[node.ntype].push_back(h);
        }
    }
    int i = 0;
    for (const auto& [node_type, ports] : node_types) {
        int_to_node_type[i] = node_type;
        int_to_pattern[i] = type_to_pattern[node_type];
        for (auto& g : node_defs) {
            if (g.ntype == node_type) {
                color_to_group[i] = g;
                break;  // Exit the loop once a match is found
            }
        }

        i++;
    }

    // Process lines and filter out incompatible graphs
    std::vector<std::tuple<int, std::vector<int>, std::vector<std::pair<int, int>>>> valid_node_colored_graphs;
    for (int line_idx = 0; line_idx < lines.size(); ++line_idx) {
        const auto& l = lines[line_idx];
        auto [n_vertices, colors, edge_list] = parse_nauty_graph_line_tmp(l, node_defs);
        if (!check_max_bond_not_exceeded_tmp(edge_list, colors, node_types, int_to_node_type)) continue;

        bool is_valid_graph = true;
        std::unordered_map<int, int> node_degrees;
        for (const auto& [src, dst] : edge_list) {
            node_degrees[src]++;
            node_degrees[dst]++;
        }

        for (int i = 0; i < n_vertices; ++i) {
            int color = colors[i];
            int degree = node_degrees[i];
            if (color_to_group.at(color).hubs.size() < degree) {
                is_valid_graph = false;
                break;
            }
        }

        if (is_valid_graph) {
            valid_node_colored_graphs.push_back(std::make_tuple(line_idx, colors, edge_list));
        }
    }
    possible_node_colored_graphs = valid_node_colored_graphs;

    if (possible_node_colored_graphs.empty()) {
        throw std::runtime_error("No valid node-colored graphs found after filtering based on hub counts.");
    }

    std::uniform_int_distribution<> dist_graph(0, int(possible_node_colored_graphs.size()) - 1);


    // Find maximum degree of any node in the node colored graphs
    std::unordered_map<int, int> max_degree; // color -> max_degree
    for (const auto& [line_idx, colors, edge_list] : possible_node_colored_graphs) {
        std::unordered_map<int, int> index_to_degree;
        for (const auto& [src, dst]: edge_list) {
            index_to_degree[src]++;
            index_to_degree[dst]++;
        }
        for (const auto& [index, degree] : index_to_degree) {
            if (max_degree.find(colors[index]) == max_degree.end()) {
                max_degree[colors[index]] = degree;
            } else {
                max_degree[colors[index]] = std::max(max_degree[colors[index]], degree);
            }
        }
    }

    // Step 1: Precompute all possible attachments for each group
    std::unordered_map<GroupGraph::Group, std::unordered_map<int, std::vector<std::vector<int>>>> group_attachments;
    for(const auto& [color, degree] : max_degree) {
        for( int i = 1; i <= degree; i++) {
            group_attachments[color_to_group[color]][i] = color_to_group[color].getPossibleAttachments(i);
        }
    }

    bool reachedMaxAttempts = true;

    omp_set_num_threads(num_procs);
    #pragma omp parallel
    {
        std::unordered_set<GroupGraph> local_basis;
        // std::random_device rd;
        std::mt19937 gen(std::random_device{}());
        std::uniform_int_distribution<> dist_graph(0, possible_node_colored_graphs.size() - 1);

        #pragma omp for schedule(dynamic)
        for (int attempt = 0; attempt < num_graphs * 1000; ++attempt) {
            if (global_basis.size() >= num_graphs) { // Check if another thread has already reached the target
                continue; // Continue to next iteration, but don't generate new graphs
            }
            GroupGraph candidate_graph;
            auto [vcolg_line_idx, colors, edge_list] = possible_node_colored_graphs[dist_graph(gen)];
            for (int color : colors) {
                GroupGraph::Group group = color_to_group[color];
                candidate_graph.addNode(group.ntype, group.pattern, group.hubs, group.patternType);
            }
            std::unordered_map<int, int> node_degrees;
            for (const auto& [src, dst] : edge_list) {
                node_degrees[src]++;
                node_degrees[dst]++;
            }
            std::vector<std::vector<int>> choosen_attachments;
            int c = 0;
            for (int color : colors) {
                int node_degree = node_degrees[c];
                std::vector<std::vector<int>> possible_attachments = group_attachments[color_to_group[color]][node_degree];
                if (possible_attachments.empty()) {
                    continue; // Skip this graph if no valid attachments for this degree
                }
                std::uniform_int_distribution<> dist_attachment(0, possible_attachments.size() - 1);
                choosen_attachments.push_back(possible_attachments[dist_attachment(gen)]);
                c++;
            }

            // Shuffle attachments for each node once to ensure random port allocation
            for (auto& attachments : choosen_attachments) {
                std::shuffle(attachments.begin(), attachments.end(), gen);
            }

            std::unordered_map<int, int> current_degree;
            for (int index = 0; index < n_nodes; index++) {
                current_degree[index] = 0;
            }
            for (const auto& [src, dst] : edge_list) {
                candidate_graph.addEdge(
                    {src, choosen_attachments[src][current_degree[src]]},
                    {dst, choosen_attachments[dst][current_degree[dst]]}
                );
                current_degree[src]++;
                current_degree[dst]++;
            }
            std::string smiles = candidate_graph.toSmiles();
            #pragma omp critical
            {
                if (global_smiles_basis.insert(smiles).second) { // If SMILES is unique
                    global_basis.insert(candidate_graph); // Add graph to the result set
                    update_progress(global_basis.size(), num_graphs);
                    if (int(global_basis.size()) >= num_graphs) {
                        reachedMaxAttempts = false; // Target reached, so max attempts not reached
                    }
                }
            }
        }
    }

    if (reachedMaxAttempts) {
        std::cout << "Warning: Maximum number of attempts reached without generating desired number of graphs..." << std::endl;
    }

    return global_basis;
}
