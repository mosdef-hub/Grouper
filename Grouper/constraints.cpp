#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <stdexcept>
#include <exprtk.hpp>




std::unordered_map<std::string, double> extract_features(const GroupGraph& g) {
    std::unordered_map<std::string, double> vars;

    // Example group counts
    vars["count_A"] = count_group(g, "A");
    vars["count_B"] = count_group(g, "B");

    // Example connection presence
    vars["connected_A_B"] = is_connected(g, "A", "B") ? 1.0 : 0.0;
    vars["connected_X_Y"] = is_connected(g, "X", "Y") ? 1.0 : 0.0;

    // Example topology
    vars["node_count"] = g.num_nodes();
    
    return vars;
}




class ExprRule {
    public:
        ExprRule(const std::string& expression_str)
            : expression_str_(expression_str) {}
    
        bool evaluate(const GroupGraph& g) const {
            exprtk::symbol_table<double> symbol_table;
            std::unordered_map<std::string, double> vars = extract_features(g);
            std::vector<double> values;
    
            for (auto& [name, val] : vars) {
                values.push_back(val);
                symbol_table.add_variable(name, values.back());
            }
    
            exprtk::expression<double> expression;
            expression.register_symbol_table(symbol_table);
    
            exprtk::parser<double> parser;
            if (!parser.compile(expression_str_, expression)) {
                throw std::runtime_error("Failed to parse constraint: " + expression_str_);
            }
    
            return expression.value() != 0.0;
        }
    
    private:
        std::string expression_str_;
    };
    


class RuleEngine {
    
    std::vector<ExprRule> rules;

    public:
        RuleEngine(const std::vector<std::string>& expressions) {
            for (const auto& expr : expressions)
                rules.emplace_back(expr);
        }

        bool satisfies_all(const GroupGraph& g) const {
            for (const auto& rule : rules)
                if (!rule.evaluate(g)) return false;
            return true;
        }
    };
