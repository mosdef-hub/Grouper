# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Grouper is a hybrid C++/Python library for building and manipulating molecular **port graphs** (graphs where edges connect specific labeled ports on nodes — see https://doi.org/10.1017/S0960129518000270). The C++ core is exposed to Python via pybind11 and shipped as a single compiled extension `Grouper._Grouper`; the Python layer adds higher-level fragmentation, visualization, library/basis-set helpers, and PyG/NetworkX interop.

## Build & environment

The build is **not** plain pip — it requires a conda environment because `CMakeLists.txt` hard-fatals if `CONDA_PREFIX` is unset and pulls headers/libs (rdkit, boost, nauty, omp, cairo, libpq) from the conda prefix.

```sh
conda env create -f environment.yml   # creates `grouper-dev`
conda activate grouper-dev
pip install -e .                      # scikit-build-core invokes CMake
```

- Build system: `scikit-build-core` + CMake ≥ 3.15, C++20 (`CMAKE_CXX_STANDARD = 20`).
- Build artifacts go to `build/{wheel_tag}/`; the editable install is in "redirect" mode.
- Native module sources are listed in `CMakeLists.txt`: `dataStructures.cpp`, `binding.cpp`, `processColoredGraphs.cpp`, `generate.cpp`, `fragmentation.cpp` → compiled into `_Grouper`.
- Apple Silicon and x86_64 macOS branches set `-Xpreprocessor -fopenmp`; Linux defines `_Thread_local=thread_local`.
- After editing C++ sources, re-run `pip install -e .` (or `cmake --build build/...`) to rebuild — there is no incremental rebuild on import.

## Testing

```sh
python -m pytest Grouper/tests                  # full suite
python -m pytest Grouper/tests/test_generate.py # one file
python -m pytest Grouper/tests/test_data_structures.py::TestGroup::test_group_hash
python -m pytest --no-skip Grouper/tests        # also runs @skip / @skipif tests
```

- `pyproject.toml` adds `--cov=Grouper --cov-report=xml:coverage.xml` automatically.
- `Grouper/tests/conftest.py` registers the custom `--no-skip` flag.
- `Grouper/tests/base_test.py` provides shared fixtures (`empty_graph`, `basic_graph`, `single_node_graph`, `single_edge_graph`, `five_member_ring_graph`); test classes inherit from `BaseTest` to consume them.
- CI (`.github/workflows/run_test.yml`) runs the matrix `{macOS-latest, ubuntu-latest} × {Python 3.9, 3.12}` via micromamba.

## Lint & format

`pre-commit` is the source of truth (`.pre-commit-config.yaml`):
- `ruff` + `ruff-format` with `--line-length=80`
- `isort --profile=black --line-length=80`
- `Grouper/tests/**` and `Grouper/libraries/Libraries.py` are excluded from the formatter and trailing-whitespace hooks — don't reformat them.

## Architecture

### C++ core (exposed via `Grouper._Grouper`)
Defined in `Grouper/dataStructures.hpp`. The two central classes:

- **`AtomGraph`** — atom-level graph. Nodes are `Atom{ntype, valency}`; edges are `(src, dst, bondOrder)`. Knows how to round-trip SMILES/SMARTS via RDKit and canonize via nauty.
- **`GroupGraph`** — port graph. Nodes are `Group{ntype, pattern, hubs, ports, patternType}`. Edges are 5-tuples `(srcNodeID, srcPort, dstNodeID, dstPort, bondOrder)`. `to_smiles()` / `to_atomic_graph()` lower a `GroupGraph` to its underlying chemistry; `canonize()` calls into nauty.

Generation entrypoints (in `generate.cpp`, bound in `binding.cpp`):
- `exhaustive_generate(n_nodes, node_defs, ...)` — enumerate the full chemical space for a size; uses nauty's `geng` + `vcolg` (the conda env supplies the binaries; some call sites still accept a `nauty_path` for legacy reasons).
- `random_generate(n_nodes, node_defs, n_graphs, ...)` — sample non-isomorphic graphs.
- `exhaustive_fragment(...)` — graph-side fragmentation; the Python wrapper is `Grouper.fragmentation.fragment`.

Custom exceptions registered in pybind: `GrouperParseException`, `GrouperFragmentationError`, `GrouperNotYetImplementedException`.

### Python layer (`Grouper/`)
- `__init__.py` re-exports `Atom`, `AtomGraph`, `Group`, `GroupGraph`, `exhaustive_fragment`, `exhaustive_generate`, `random_generate` from the compiled `_Grouper`, plus the Python `fragment` wrapper.
- `fragmentation.py` — high-level `fragment(smiles, nodeDefs, returnHandler={"ideal","quick","exhaustive"}, ...)` that decomposes a SMILES into `GroupGraph`s using a set of `Group`/`GroupExtension` definitions. Read its docstring before changing matching semantics — the `returnHandler`, `nodeDefsSorter`, `incompleteGraphHandler`, and `matchHubs` flags interact subtly.
- `utils.py` — performance benchmarking helpers, `convert_to_nx`, and PyG (`torch_geometric.Data`) conversion.
- `counting.py` — Pólya-style enumeration helpers built on `igraph` + `sympy.combinatorics`.
- `io.py` — soft-import shim that raises `DelayImportError` (a `SkipTest` subclass) with install instructions when optional deps (`mbuild`, `torch`, `rdkit`) are missing. Use `import_("torch")` instead of bare `import torch` in code paths that should degrade gracefully in tests.
- `libraries/Libraries.py` — predefined group libraries: `Joback`, `Unifac`, `SaftGammaMie`, `BasisSet`, `GroupExtension`. **Excluded from formatters.**
- `visualization/` — graph drawing (`visualize_graph`), single-group rendering (`visualize_node`), and chemical-space network analytics (`visualize_space`: centrality, modularity, diversity, etc.).

### Data flow worth knowing
1. User defines `Group(ntype, pattern, hubs)` objects (pattern is SMILES by default; pass `pattern_type="SMARTS"` for SMARTS).
2. Graphs are constructed by `add_node` + `add_edge((nodeID, port), (nodeID, port), bondOrder)`. Ports are 0-indexed positions into `hubs`.
3. `to_atomic_graph()` → `AtomGraph`; `to_smiles()` → string. Both rely on RDKit at the C++ layer.
4. Enumeration funnels through nauty for canonicalization; that's why the conda env pinning the nauty headers/libs is non-negotiable.

## Notes for working in this repo

- Many top-level files in the repo root (`debug*.py`, `*.txt`, `*.png`, `*.pkl`, `publication_figures/`, `notebooks/`, untracked figure scripts) are research artefacts, not part of the package. The shipped package lives entirely under `Grouper/`. Don't refactor or delete root-level scratch files unless the user asks.
- When adding a new C++ symbol, expose it in `binding.cpp`, then re-export from `Grouper/__init__.py` so it shows up as `Grouper.<name>`.
- The repo ships prebuilt `.so` files in `Grouper/` (e.g. `_Grouper.cpython-312-darwin.so`) — these are stale local builds, not authoritative. Always rebuild for the active interpreter.
