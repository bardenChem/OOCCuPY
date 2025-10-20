"""
OOCCuPY - Object Oriented Computational Chemistry in Python

A comprehensive computational chemistry package with modules for:
- pDynamo wrappers
- Molecular dynamics tools  
- QM input generation
"""

__version__ = "1.0.0"
__author__ = "Igor Barden Grillo"
__email__ = "barden.igor@gmail.com"

# Import main classes/functions for easy access
from .pDynamoWrapper import Wrapper
from .mdtools import MDSimulation
from .QM_inputs import QMInputGenerator


# Define what gets imported with "from OOCCuPY import *"
__all__ = [
    'Wrapper',
    'MDSimulation', 
    'QMInputGenerator',
    'pDynamoWrapper',
    'mdtools',
    'QM_inputs'
]

# Package metadata
__package_name__ = "OOCCuPY"
__description__  = "Object Oriented Computational Chemistry in Python"
__url__          = "https://github.com/bardenChem/OOCCuPY"
__license__      = "MP2"