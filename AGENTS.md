# AGENTS.md — OOCCuPY Development Guide

## 1. Project mission

OOCCuPY is a free, Python-based computational chemistry project being developed as a modular API and CLI for preparing, running, monitoring, parsing, and analyzing molecular simulations through reproducible scientific protocols.

The immediate objective is to prepare a stable, demonstrable release for the Escola de Modelagem Molecular de Sistemas Biológicos at LNCC. The following objectives are to harden the software, register it, and prepare a software/methodology article for the *Journal of Chemical Information and Modeling*.

pDynamo3 remains the central backend for the first stable release, especially for QM/MM, molecular dynamics, geometry optimization, relaxed scans, umbrella sampling, energy refinement, trajectory analysis, and electronic-result generation used by PRIMoRDiA. GROMACS, OpenMM, xTB, and CP2K must be added incrementally through isolated adapters and small, testable workflows.

## 2. Read this before changing code

Before proposing or implementing a change:

1. Inspect the repository tree and the source files involved.
2. Inspect existing tests, examples, fixtures, configuration files, and CLI entry points.
3. Do not assume that the README accurately reflects the current implementation.
4. Distinguish clearly between:
   - implemented and tested behavior;
   - implemented but untested behavior;
   - partial or placeholder behavior;
   - planned behavior.
5. Do not invent classes, functions, command options, file formats, or backend behavior.
6. When required information is absent, state the assumption and identify what must be inspected.
7. Prefer a small, verifiable change over a broad rewrite.

The currently observed repository snapshot contains important incomplete areas:

- `OOCCuPY.py` routes the `pDynamo`, `mdtools`, `QM_inputs`, and `config` commands.
- `Interface._run_subprocess_real_time()` is a placeholder.
- `Interface.MD_prep_handler()` is a placeholder.
- `Interface.QM_input_handler()` is a placeholder.
- the CLI redirects global `sys.stdout` through `Tee`;
- `pyproject.toml` currently contains only the build-system section;
- `requirements.txt` contains runtime packages but no development/test dependency groups;
- capabilities described in the README must be verified against source and tests before being treated as public functionality.

Treat these observations as a starting point, not as a complete architecture audit.

## 3. Release priorities

Use this order unless a task explicitly overrides it.

### P0 — release blockers and scientific correctness

1. Installation and editable installation work in a clean environment.
2. The CLI starts, shows help, returns meaningful exit codes, and reports errors without tracebacks for expected user mistakes.
3. Critical pDynamo3 workflows run reproducibly.
4. Existing scientifically validated examples are converted into regression or smoke tests.
5. Subprocess failures, missing executables, parsing failures, missing files, and invalid configurations are never silently ignored.
6. Changes that can affect numerical results are protected by regression tests.
7. Inputs, outputs, software versions, parameters, units, and working directories are recorded.

### P1 — first LNCC release scope

1. Stabilize the pDynamo3 automation layer.
2. Build a minimal GROMACS workflow:
   - preparation;
   - minimization;
   - equilibration;
   - production;
   - basic trajectory/energy analysis.
3. Build a minimal OpenMM workflow:
   - system construction from supported inputs;
   - minimization;
   - equilibration;
   - production;
   - state and trajectory reporters.
4. Define and implement a safe rapid-parameterization path for common organic molecules in biological simulations.
5. Add xTB support for:
   - single-point calculations;
   - geometry optimization;
   - molecular dynamics.
6. Add an initial CP2K adapter for:
   - input generation;
   - single-point calculations;
   - geometry optimization;
   - parsing of essential results.
7. Align the README and examples with verified functionality.

### P2 — post-release hardening

1. Standardize configuration and result models.
2. Separate protocol definitions from backend adapters.
3. Centralize local execution and process monitoring.
4. Expand real-backend integration tests.
5. Improve documentation, examples, changelog, citations, and known limitations.
6. Freeze and document a minimal public API.

### P3 — later scientific expansion

Do not include these in the first stable release unless explicitly requested and fully tested:

- Q-Force-based specialized parameterization;
- relaxed QM/MM scans with xTB;
- CP2K QM/MM metadynamics;
- advanced CP2K free-energy workflows;
- artificial common abstractions that hide scientifically important backend differences.

## 4. Architecture direction

Keep the following responsibilities separate:

- data models and validation;
- molecular structure handling;
- system preparation;
- scientific protocol configuration;
- backend-specific input generation;
- local or cluster execution;
- process monitoring;
- output parsing;
- normalized result models;
- scientific analysis;
- plots, images, and reports;
- CLI;
- Python API.

Recommended package direction:

```text
ooccupy/
├── cli/
├── config/
├── core/
│   ├── models/
│   ├── exceptions.py
│   ├── units.py
│   └── paths.py
├── execution/
│   ├── local.py
│   ├── process.py
│   └── records.py
├── protocols/
│   ├── optimization.py
│   ├── dynamics.py
│   ├── single_point.py
│   └── sampling.py
├── backends/
│   ├── pdynamo/
│   ├── gromacs/
│   ├── openmm/
│   ├── xtb/
│   └── cp2k/
├── preparation/
├── parsing/
├── analysis/
└── reporting/
```

Do not perform a repository-wide migration to this structure in one change. Introduce boundaries incrementally around working code.

### Protocols versus adapters

A scientific protocol describes intent and parameters, such as minimization, NVT equilibration, NPT equilibration, production MD, single point, or optimization.

A backend adapter translates that protocol into a specific engine invocation and parses its results.

Do not force incompatible concepts into one universal interface. Shared interfaces are appropriate for common lifecycle operations such as:

- validate configuration;
- prepare inputs;
- execute;
- parse;
- return a result object.

Backend-specific capabilities and options may remain explicit.

## 5. Scientific safety rules

Changes involving any item below are high risk:

- energies and energy decomposition;
- units or unit conversion;
- atom indexing or atom mapping;
- residue and atom selections;
- QM-region selection;
- total charge and spin multiplicity;
- force-field assignment;
- constraints and restraints;
- reaction coordinates;
- trajectory formats;
- periodic boxes;
- thermostat and barostat settings;
- integration timestep;
- random seeds;
- convergence thresholds;
- Molden generation;
- PySCF/pDynamo3 interface behavior.

For high-risk changes:

1. identify the scientific behavior being preserved or changed;
2. add a regression test before or with the change;
3. use justified numerical tolerances;
4. preserve a small reference input and expected result;
5. document any intentional change in output;
6. never silently convert units or infer charges;
7. never silently reorder atoms without preserving a mapping.

## 6. Configuration and result models

Prefer Pydantic models or focused dataclasses for new configuration and result objects.

Every simulation configuration should make relevant values explicit:

- backend and backend version;
- method;
- basis or parameter set where applicable;
- charge and multiplicity;
- ensemble;
- temperature and pressure;
- timestep and number of steps;
- thermostat/barostat;
- constraints/restraints;
- seed;
- resources;
- input structure;
- working directory;
- output names;
- restart policy.

Configuration should be serializable to YAML or JSON.

Results should not be represented only by loosely structured dictionaries. A result object should include, when applicable:

- success/failure status;
- backend;
- command;
- return code;
- start/end timestamps;
- input and output paths;
- parsed energies and units;
- trajectory/checkpoint paths;
- warnings;
- software version;
- metadata needed for reproduction.

Never store a numerical value without its unit or a documented unit convention.

## 7. External process execution

Centralize subprocess execution. Do not add new scattered `subprocess.run()` or `Popen()` calls inside scientific classes unless wrapping legacy code temporarily.

The common execution layer should support:

- command as a sequence of arguments;
- explicit working directory;
- controlled environment;
- real-time stdout/stderr streaming;
- complete log capture;
- timeout;
- return code validation;
- executable discovery;
- version probing;
- dry-run mode where useful;
- structured execution record;
- specific exceptions.

Never use `shell=True` unless there is a documented, unavoidable reason and inputs are fully controlled.

A missing executable should raise a backend-specific availability error. Tests must be able to replace real execution with a fake runner.

## 8. Logging and errors

Use `logging`, not dispersed `print()` calls, for library behavior.

CLI code may render concise user-facing messages, but internal code should emit structured logs.

Define specific exceptions, for example:

- `OOCCuPYError`;
- `ConfigurationError`;
- `BackendUnavailableError`;
- `InputGenerationError`;
- `ExecutionError`;
- `ParsingError`;
- `ScientificValidationError`.

Do not catch broad exceptions only to print a message and continue. Preserve the original exception with `raise ... from exc`.

Do not rely on `__del__` for critical cleanup. Prefer context managers and `try/finally`.

## 9. Testing strategy

Use `pytest`.

Organize tests into:

```text
tests/
├── unit/
├── integration/
├── regression/
├── cli/
├── smoke/
├── fixtures/
└── optional_backends/
```

### Unit tests

Cover:

- validators;
- parsers;
- path handling;
- configuration models;
- input writers;
- command construction;
- unit conversion;
- result normalization;
- error branches.

Unit tests must not require installed scientific executables.

### Integration tests

Cover interactions between:

- configuration and protocol;
- protocol and adapter;
- adapter and fake runner;
- runner and parser;
- CLI and public API.

### Regression tests

Convert scientifically validated examples into regression cases. Prefer small systems and store only essential fixtures.

Regression assertions may include:

- final energy within a justified tolerance;
- number and identity of selected QM atoms;
- reaction-coordinate values;
- expected output files;
- atom count and ordering;
- parsed metadata;
- Molden file presence and essential sections.

Do not compare large binary trajectories byte-for-byte.

### CLI tests

Use pytest facilities or a CLI testing library appropriate to the selected CLI framework.

Test:

- top-level help;
- subcommand help;
- missing arguments;
- invalid paths;
- missing backend;
- successful dry run;
- successful fake execution;
- nonzero exit codes;
- log creation;
- no-command behavior.

### External backends

Use three layers:

1. **fake runner tests** — always run;
2. **mocked executable tests** — always run where practical;
3. **real backend tests** — optional and marked.

Suggested markers:

```ini
slow
integration
regression
requires_pdynamo
requires_gromacs
requires_openmm
requires_xtb
requires_cp2k
```

A skipped optional backend test is acceptable only when the reason is explicit.

### Numerical assertions

Use `pytest.approx` or NumPy testing helpers. Set tolerances based on method stability, precision, platform variation, and the scientific purpose of the test. Document non-obvious tolerances.

## 10. Backend-specific first-release expectations

### pDynamo3

Prioritize existing, scientifically validated workflows. Before refactoring a workflow:

1. locate its example/input;
2. identify generated files and numerical outputs;
3. capture the baseline;
4. add a regression test;
5. refactor incrementally.

Keep compatibility with the modified pDynamo3 fork where required. Do not replace fork-specific behavior without inspecting the corresponding source changes.

### GROMACS

The minimal adapter should expose a protocol-driven sequence rather than one monolithic class.

Expected stages:

- topology/coordinate preparation from supported inputs;
- energy minimization;
- NVT equilibration;
- NPT equilibration when applicable;
- production;
- basic energy and trajectory outputs.

Record generated `.mdp` files, topology, coordinates, commands, logs, and GROMACS version.

Do not claim that arbitrary ligands are supported until the parameterization path and validation are implemented.

### OpenMM

Prefer its native Python API. Keep OpenMM-specific objects inside the adapter boundary.

The minimal workflow should support:

- platform selection;
- deterministic seed where supported;
- minimization;
- equilibration;
- production;
- checkpointing;
- trajectory and state reporters;
- reproducibility metadata.

Tests should run a tiny system on the Reference or CPU platform when OpenMM is installed.

### Rapid parameterization

The first version targets common organic molecules and routine biological simulation contexts.

Requirements:

- use established external APIs, toolkits, or accessible models;
- preserve provenance of parameters;
- validate charge, atom mapping, element support, and output completeness;
- expose unsupported chemistry clearly;
- avoid silent fallback to guessed atom types;
- cache external responses only with version/provenance metadata;
- make remote API use optional and explicit;
- provide a fully local test path with recorded fixtures or fakes.

Q-Force and specialized parameterization are future scope.

### xTB

First-release calculations:

- single point;
- geometry optimization;
- molecular dynamics.

The adapter should generate explicit commands/inputs, capture xTB version, parse final energy and convergence, preserve optimized coordinates, and expose trajectory/restart files for MD.

Do not assume that every xTB feature is available in every installed version. Probe or validate capabilities.

### CP2K

Initial scope:

- input generation;
- single point;
- geometry optimization;
- essential output parsing.

Keep the input model extensible but do not expose untested metadynamics or QM/MM free-energy workflows as stable API.

Preserve CP2K input files exactly as executed and record data-file dependencies, basis sets, potentials, MPI command, and executable version.

## 11. CLI and public API

The CLI should be a thin layer over the same service functions used by the Python API.

Avoid embedding scientific logic inside argument handlers.

Required release behavior:

- `ooccupy --help`;
- `ooccupy --version`;
- useful subcommand help;
- consistent option spelling;
- validated paths;
- stable exit codes;
- clear backend-unavailable messages;
- optional verbose/debug logging;
- no global `sys.stdout` replacement in library code.

When modifying existing command names, provide a compatibility path or document a breaking change.

## 12. Packaging and dependencies

Complete `pyproject.toml` before release with, at minimum:

- project metadata;
- package discovery;
- supported Python range;
- runtime dependencies;
- optional dependency groups;
- console script entry point;
- license;
- authors;
- repository/documentation URLs;
- package data rules.

Separate dependency groups conceptually:

- core;
- development;
- documentation;
- pDynamo;
- GROMACS helpers;
- OpenMM;
- xTB;
- CP2K;
- analysis/visualization.

Do not list Python standard-library packages such as `argparse` as runtime dependencies.

Avoid forcing all external backends into the base installation.

## 13. Code style

For new or substantially modified code:

- use Python 3 type hints;
- use `pathlib.Path`;
- add focused docstrings;
- keep functions small;
- avoid mutable default arguments;
- avoid wildcard imports;
- use explicit encodings for text files;
- validate user input at boundaries;
- return structured values;
- keep scientific constants named and documented.

Gradually introduce:

- Ruff for linting and formatting;
- pytest and pytest-cov;
- mypy or pyright for selected modules;
- pre-commit;
- CI.

Do not perform formatting-only changes across unrelated scientific modules in the same pull request as behavioral changes.

## 14. Change workflow for Codex

For any nontrivial task, follow this sequence.

### Step 1 — inspect

Report:

- files and symbols involved;
- current behavior;
- incomplete behavior;
- tests/examples that cover it;
- scientific risks.

### Step 2 — plan

Provide small steps with:

- files affected;
- purpose of each change;
- compatibility concerns;
- test to add;
- acceptance criterion.

### Step 3 — implement

Make the smallest coherent change. Include imports, typing, error handling, and documentation.

### Step 4 — test

Run the narrowest relevant tests first, then the broader suite.

Report exact commands and results. Never claim a test passed if it was not executed.

### Step 5 — review

Check:

- error paths;
- cleanup;
- units;
- atom ordering/mapping;
- output preservation;
- backwards compatibility;
- documentation impact.

## 15. Pull request and commit expectations

Each change should state:

- problem;
- scope;
- implementation;
- scientific impact;
- tests executed;
- limitations;
- follow-up work.

Keep commits focused. Do not mix:

- broad formatting;
- architecture migration;
- scientific behavior changes;
- dependency updates;
- documentation rewrites

unless they are inseparable.

## 16. Definition of done

A task is complete only when:

1. behavior is implemented;
2. validation and errors are handled;
3. tests cover success and important failure paths;
4. relevant tests were executed;
5. documentation or examples are updated;
6. no undocumented scientific behavior changed;
7. generated inputs and outputs are reproducible;
8. the task's acceptance criterion is demonstrably satisfied.

A backend workflow is not complete merely because it generates a file. It must validate inputs, execute or dry-run predictably, detect failure, parse essential outputs, and preserve provenance.

## 17. Release gates

Do not mark the LNCC release ready until:

- clean installation succeeds;
- CLI smoke tests pass;
- P0 pDynamo workflows pass;
- known critical scientific regressions pass;
- missing optional backends fail gracefully;
- minimal GROMACS and OpenMM demonstrations are reproducible;
- xTB and CP2K advertised features are tested;
- README matches actual behavior;
- version, license, citation, changelog, and known limitations are present;
- a release candidate has been tested from a clean checkout.

## 18. Communication style

Use clear technical Portuguese when discussing the project with the maintainer, unless English is requested.

In code reviews, list bugs, risks, and inconsistencies before improvements.

In planning, provide small tasks, dependencies, acceptance criteria, and tests.

In implementation, favor readable and maintainable code over clever abstractions.

Never remove or alter validated scientific functionality without explicit justification and regression evidence.
