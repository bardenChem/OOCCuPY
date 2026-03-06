# OOCCuPY

**Object-Oriented Computational Chemistry in Python**

[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![GitHub](https://img.shields.io/badge/GitHub-OOCCuPY-blue)](https://github.com/igorChem/OOCCuPY)

OOCCuPY is a comprehensive Python framework for computational chemistry that provides object-oriented interfaces to molecular dynamics, quantum mechanics, and molecular structure analysis tools. It features powerful wrappers for popular quantum chemistry packages (MOPAC, GAMESS, ORCA), molecular dynamics engines (AMBER, GROMACS), and the pDynamo library for hybrid QM/MM simulations.

## Table of Contents

- [Features](#features)
- [System Requirements](#system-requirements)
- [Installation](#installation)
  - [Quick Start (PyPI)](#quick-start-pypi)
  - [Installation from Source](#installation-from-source)
  - [Development Installation](#development-installation)
- [Configuration](#configuration)
- [Usage](#usage)
  - [Command-Line Interface](#command-line-interface)
  - [Python API](#python-api)
  - [pDynamo Wrapper Examples](#pdynamo-wrapper-examples)
  - [Quantum Chemistry Examples](#quantum-chemistry-examples)
  - [Molecular Dynamics Examples](#molecular-dynamics-examples)
- [Module Overview](#module-overview)
- [Examples and Tests](#examples-and-tests)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)
- [Citation](#citation)
- [Support](#support)

## Features

### Core Capabilities

- **Object-Oriented Design**: Clean, Pythonic interfaces to computational chemistry workflows
- **Hybrid QM/MM Simulations**: Seamless integration of quantum mechanics and molecular mechanics through pDynamo wrapper
- **Quantum Chemistry Support**: 
  - MOPAC (semi-empirical methods)
  - GAMESS (ab initio methods)
  - ORCA (modern DFT and wavefunction methods)
- **Molecular Dynamics Support**:
  - AMBER (molecular dynamics)
  - GROMACS (enhanced sampling, MD)
- **Structure Analysis**: Tools for PDB and XYZ file manipulation and analysis
- **Advanced Simulations**:
  - Molecular dynamics with multiple force fields
  - Energy minimization and geometry optimization
  - Umbrella sampling and enhanced sampling techniques
  - Relaxed scans and reaction coordinate analysis
  - QM/MM refinement procedures

### Simulation Capabilities

- Geometry optimization and initial structure preparation
- Molecular dynamics simulations (equilibration and production)
- Energy analysis and trajectory analysis
- Reaction coordinate studies
- Potential of mean force (PMF) calculations
- Quantum chemical input generation and output parsing

## System Requirements

### Python Requirements

- Python 3.8 or higher
- pip (Python package installer)

### Dependencies

```
numpy>=1.21.0
scipy>=1.7.0
matplotlib>=3.5.0
pandas>=1.3.0
pyyaml>=6.0
```

### External Software (Optional but Recommended)

For full functionality, install:

- **Quantum Chemistry**: MOPAC, GAMESS, or ORCA
- **Molecular Dynamics**: AMBER or GROMACS
- **Hybrid QM/MM**: pDynamo3 library
- **Visualization**: PyMOL (for structure visualization)

### System Requirements

- **OS**: Linux, macOS, or Windows (with appropriate MD/QC software)
- **RAM**: Minimum 4GB (8GB+ recommended for MD simulations)
- **Disk Space**: 2GB for installation + additional space for simulation outputs

## Installation

### Quick Start (PyPI)

The easiest way to install OOCCuPY is from PyPI:

```bash
pip install ooccupy
```

This installs OOCCuPY with all required Python dependencies and sets up user configuration directories.

### Installation from Source

If you want to install from the latest development version:

```bash
# Clone the repository
git clone https://github.com/igorChem/OOCCuPY.git
cd OOCCuPY

# Install dependencies
pip install -r requirements.txt
pip install pyyaml

# Install OOCCuPY
pip install .
```

### Development Installation

For development work and contributions:

```bash
# Clone the repository
git clone https://github.com/igorChem/OOCCuPY.git
cd OOCCuPY

# Install in development mode
pip install -e .

# Install additional development dependencies (optional)
pip install pytest pytest-cov black flake8
```

### Post-Installation Setup

After installation, OOCCuPY automatically sets up configuration directories:

```
~/.ooccupy/
├── config.yaml          # Configuration file
├── data/                # Data files directory
├── Examples/            # Example input files
├── Tests/               # Test files
├── logs/                # Log files
├── cache/               # Cache directory
└── temp/                # Temporary files
```

## Configuration

### Viewing Configuration

```bash
ooccupy config --show
```

### Custom Configuration

Edit the configuration file at `~/.ooccupy/config.yaml` to customize:

- Quantum chemistry method defaults
- Molecular dynamics parameters
- Computational resource limits
- Logging preferences

Example configuration:

```yaml
general:
  verbose: true
  log_level: INFO
  log_to_file: true
  use_cache: true
  cache_size: 1000

computational:
  default_method: DFT
  default_basis: '6-31G*'
  max_memory: 4GB

quantum_chemistry:
  mopac_method: RM1
  gamess_contrl: ''
  orca_keywords: 'B3LYP 6-31G(d)'

molecular_dynamics:
  force_field: AMBER99SB
  temperature: 298.15
  pressure: 1.0
```

## Usage

### Command-Line Interface

#### Help

```bash
# Show all available commands
ooccupy --help

# Show configuration
ooccupy config --show
```

#### pDynamo Wrapper

```bash
# Run pDynamo wrapper with input file
ooccupy pdynamo --input example.input --output results.txt

# Run in specific project folder
ooccupy pdynamo --input example.input --proj-folder ./my_project

# Run pDynamo tests
ooccupy pdynamo --tests

# Run specific test number
ooccupy pdynamo --test 01
```

#### Other Commands

```bash
# Version information
ooccupy --version

# Verbose output
ooccupy --verbose pdynamo --input example.input
```

### Python API

#### Basic pDynamo Usage

```python
from pDynamoWrapper.pDynamoWrapper import Wrapper
from pDynamoWrapper.SimulationSystem import SimulationSystem

# Initialize wrapper
wrapper = Wrapper(project_folder="my_simulation")

# Create or load a system
system = SimulationSystem("protein.pdb")

# Perform geometry optimization
system.minimize(method="DFT", basis="6-31G*")

# Run molecular dynamics
system.run_md(temperature=298.15, pressure=1.0, steps=10000)

# Analyze results
system.analyze_trajectory("trajectory.dcd")
```

#### Working with Quantum Chemistry

```python
from QM_inputs.mopac_module import MOPAC_input
from QM_inputs.OrcaModule import ORCA_input

# Create MOPAC input
mopac = MOPAC_input("molecule.xyz", method="PM6")
mopac.set_keywords("PRECISE GRADIENTS")
mopac.write_input("mopac_input.mop")

# Create ORCA input
orca = ORCA_input("molecule.xyz", method="B3LYP", basis="6-31G*")
orca.write_input("orca_input.inp")
```

#### Molecular Dynamics Preparation

```python
from MD_tools.md_prep import MDSimulation

# Setup MD simulation
md = MDSimulation("protein.pdb")

# Configure force field (AMBER)
md.setup_amber(charge_method="gasteiger")

# Prepare system for simulation
md.minimize()  # Energy minimization
md.equilibrate(steps=1000)  # Equilibration
md.produce(steps=10000)  # Production
```

#### Structure Analysis

```python
from Structure.pdb_class import PDB

# Load and analyze structure
pdb = PDB("protein.pdb")

# Get information
print(f"Number of atoms: {len(pdb.atoms)}")
print(f"Chains: {pdb.get_chains()}")

# Calculate properties
rmsd = pdb.calculate_rmsd(reference_pdb)
contacts = pdb.calculate_contacts()
```

### pDynamo Wrapper Examples

OOCCuPY includes extensive examples for the pDynamo wrapper. Run them with:

```bash
# Run all pDynamo examples
ooccupy pdynamo --tests

# Run specific example (e.g., Example 1)
python Examples/pDynamoWrapper/Example_1.input
```

#### Example 1: Basic Geometry Optimization

```python
from pDynamoWrapper.pDynamoWrapper import Wrapper

wrapper = Wrapper.From_Input("Examples/pDynamoWrapper/Example_1.input")
```

#### Example 2: Molecular Dynamics Simulation

```python
# Build system and run MD simulation with pDynamo
# See Examples/pDynamoWrapper/Example_2.input
```

#### Example 3: QM/MM Calculations

```python
# Hybrid quantum mechanics / molecular mechanics
# See Examples/pDynamoWrapper/Example_3.input
```

#### Example 4-16: Advanced Simulations

- Umbrella sampling studies
- Relaxed scans
- Energy refinement
- Reaction coordinate analysis
- And more!

Browse the `Examples/pDynamoWrapper/` directory for complete input file examples.

### Quantum Chemistry Examples

#### MOPAC Calculations

```python
from QM_inputs.mopac_module import mopac_mod
from Structure.xyz_class import XYZ

# Load structure
xyz = XYZ("molecule.xyz")

# Setup MOPAC calculation
mopac = mopac_mod("molecule.xyz")
mopac.set_method("PM6")
mopac.add_keywords("PRECISE GRADIENTS")
mopac.run()

# Parse output
from QM_inputs.mopac_out import mopac_out
result = mopac_out("molecule.out")
energy = result.get_energy()
gradients = result.get_gradients()
```

#### GAMESS Calculations

```python
from QM_inputs.gms_inp import GAMESS_input

# Create GAMESS input
gamess = GAMESS_input("molecule.xyz")
gamess.set_runtype("OPTIMIZE")
gamess.set_method("RHF")
gamess.set_basis("STO-3G")
gamess.write_input("gamess_input.inp")
```

#### ORCA Calculations

```python
from QM_inputs.OrcaModule import ORCA_input

# Create ORCA input
orca = ORCA_input("molecule.xyz")
orca.set_method("B3LYP")
orca.set_basis("6-31G*")
orca.set_keywords("OPT FREQ")
orca.write_input("orca_input.inp")
```

### Molecular Dynamics Examples

#### AMBER Simulations

```python
from MD_tools.amber_module import amber_mod

# Setup AMBER simulation
amber = amber_mod("protein.pdb")

# Prepare system
amber.tleap_call()  # Prepare topology
amber.equilibration()  # Equilibrate
amber.production()  # Produce trajectory

# Analyze results
amber.analysis_MD()
```

#### GROMACS Simulations

```python
from MD_tools.gmx_module import min_prot

# Setup GROMACS minimization
gromacs = min_prot("protein.gro", "protein.top")
gromacs.run()  # Run minimization
```

## Module Overview

### pDynamoWrapper
- **Description**: Object-oriented wrapper for pDynamo3 molecular dynamics library
- **Key Classes**: `Wrapper`, `SimulationSystem`, `Simulation`
- **Features**: QM/MM simulations, geometry optimization, MD, analysis
- **Use Case**: Hybrid quantum mechanics/molecular mechanics calculations

### QM_inputs
- **Description**: Quantum chemistry input/output handling
- **Supported Software**: MOPAC, GAMESS, ORCA, FMO
- **Key Classes**: `mopac_out`, `GAMESS_input`, `ORCA_input`
- **Features**: Input generation, output parsing, property extraction
- **Use Case**: Quantum chemical calculations and electronic structure analysis

### MD_tools
- **Description**: Molecular dynamics preparation and analysis
- **Supported Software**: AMBER, GROMACS
- **Key Modules**: `amber_module`, `gmx_module`, `md_prep`, `md_analysis`
- **Features**: System setup, simulation protocols, trajectory analysis
- **Use Case**: Classical molecular dynamics simulations and analysis

### Structure
- **Description**: Molecular structure manipulation and analysis
- **Supported Formats**: PDB, XYZ
- **Key Classes**: `PDB`, `XYZ`
- **Features**: Structure parsing, geometry calculations, coordinate manipulation
- **Use Case**: Structure preparation and analysis

### Analysis
- **Description**: Comprehensive simulation analysis tools
- **Key Classes**: `Analysis`, `EnergyAnalysis`, `TrajectoryAnalysis`
- **Features**: Energy breakdown, trajectory parsing, property calculations
- **Use Case**: Post-simulation analysis and result extraction

### Simulation & SimulationSystem
- **Description**: Core simulation management classes
- **Features**: System setup, simulation execution, result collection
- **Use Case**: Orchestrating complex multi-step simulations

## Examples and Tests

### Running Examples

```bash
# Run all pDynamo examples
cd Examples/pDynamoWrapper
python Example_1.input
python Example_2.input
# ... and so on
```

### Running Tests

```bash
# Run all tests
ooccupy pdynamo --tests

# Run specific test
ooccupy pdynamo --test 01

# Or directly with Python
cd Tests/pDynamoWrapper
python test_01.py
```

Test files are located in:
- Development: `Tests/pDynamoWrapper/test_*.py`
- Installed: `~/.ooccupy/Tests/pDynamoWrapper/`

### Sample Data

Sample molecular structures are provided in the `data/` directory:
- `1atp_peptide.gro` / `1atp_peptide.top` - ATP-peptide GROMACS files
- `7tim.crd` / `7tim.top` - Triose phosphate isomerase AMBER files
- `cyclohexane_single_frame.xyz` - Cyclohexane XYZ file

## Troubleshooting

### Installation Issues

**Problem**: `ModuleNotFoundError: No module named 'pDynamoWrapper'`

**Solution**: Ensure OOCCuPY is properly installed:
```bash
pip install --force-reinstall ooccupy
```

**Problem**: Missing external dependencies (AMBER, GROMACS, etc.)

**Solution**: Install required external software separately and ensure they're in your PATH:
```bash
# Check if AMBER is available
which sander
which tleap

# Check if GROMACS is available
which gmx
```

### Configuration Issues

**Problem**: Configuration not found

**Solution**: Manually initialize configuration:
```python
from config import OOCCuPYConfig
config = OOCCuPYConfig()
```

**Problem**: External software not found

**Solution**: Update `~/.ooccupy/config.yaml` with correct paths:
```yaml
paths:
  amber_home: /path/to/amber
  gromacs_home: /path/to/gromacs
  mopac_exe: /path/to/mopac
```

### Simulation Issues

**Problem**: Simulation crashes with "pDynamo not found"

**Solution**: Ensure pDynamo3 is properly installed. OOCCuPY requires the pDynamo3 library to be available in your Python environment. Install it separately if needed.

**Problem**: Memory errors during large MD simulations

**Solution**: Reduce simulation size or increase available RAM. Modify `max_memory` in `~/.ooccupy/config.yaml`.

### Path Issues

**Problem**: Data files not found

**Solution**: Use the `find_data_file` utility:
```python
from config import find_data_file
data_file = find_data_file("1atp_peptide.gro")
```

## Contributing

We welcome contributions! To contribute:

1. Fork the repository on GitHub
2. Create a feature branch:
   ```bash
   git checkout -b feature/your-feature-name
   ```
3. Make your changes and add tests
4. Ensure code quality:
   ```bash
   flake8 your_module.py
   black your_module.py
   ```
5. Commit with clear messages:
   ```bash
   git commit -m "Add feature: description"
   ```
6. Push to your fork and submit a pull request

### Development Setup

```bash
# Clone and setup development environment
git clone https://github.com/igorChem/OOCCuPY.git
cd OOCCuPY
pip install -e .
pip install pytest black flake8

# Run tests
pytest Tests/

# Format code
black .

# Check code quality
flake8 .
```

## License

OOCCuPY is licensed under the **Mozilla Public License 2.0 (MPL-2.0)**. 

This is a permissive open-source license that allows you to:
- ✅ Use commercially
- ✅ Modify the software
- ✅ Distribute modified versions
- ⚠️ Must include a copy of the license
- ⚠️ Must document significant changes

See [LICENSE](LICENSE) for full details.

## Citation

If you use OOCCuPY in your research, please cite:

```bibtex
@software{ooccupy2024,
  title={OOCCuPY: Object-Oriented Computational Chemistry in Python},
  author={Grillo, Igor Barden},
  year={2024},
  url={https://github.com/igorChem/OOCCuPY}
}
```

## Support

### Documentation

- **Examples**: See `Examples/pDynamoWrapper/` for input file examples
- **Tests**: Run `Tests/pDynamoWrapper/` test files as references
- **Module Docstrings**: Use Python's `help()` function for module details:
  ```python
  from pDynamoWrapper import Wrapper
  help(Wrapper)
  ```

### Getting Help

- **GitHub Issues**: Report bugs at https://github.com/igorChem/OOCCuPY/issues
- **Documentation**: Check module docstrings and examples
- **Email**: Contact the developer at barden.igor@gmail.com

### Useful Resources

- [pDynamo Documentation](http://www.christophergronbeck.net/pDynamo/)
- [MOPAC Documentation](http://openmopac.net/)
- [ORCA Documentation](https://www.orcasoftware.de/)
- [AMBER Documentation](https://ambermd.org/)
- [GROMACS Documentation](https://manual.gromacs.org/)

---

**Version**: 1.0.0  
**Author**: Igor Barden Grillo  
**Repository**: https://github.com/igorChem/OOCCuPY  
**License**: Mozilla Public License 2.0 (MPL-2.0)
