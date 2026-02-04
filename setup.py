#!/usr/bin/env python3
"""
OOCCuPY Installer - Flat Structure Version
"""

import os
import sys
import shutil
from pathlib import Path
from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop


class PostInstallMixin:
    """Mixin for post-installation tasks for flat structure"""
    
    def setup_user_environment(self):
        """Setup user environment after installation"""
        print("\n" + "="*50)
        print("OOCCuPY Post-Installation Setup")
        print("="*50)
        
        # Create user config directory
        user_config_dir = Path.home() / '.ooccupy'
        user_config_dir.mkdir(exist_ok=True)
        
        # Create subdirectories
        for subdir in ['Tests', 'Examples', 'data', 'logs', 'cache']:
            (user_config_dir / subdir).mkdir(exist_ok=True)
        
        # Initialize configuration using your config.py
        try:
            # Add current directory to path
            sys.path.insert(0, str(Path(__file__).parent))
            
            from config import OOCCuPYConfig
            config = OOCCuPYConfig(user_config_dir)
            
            print(f"✓ Configuration created at: {user_config_dir}")
            print(f"✓ Test data directory: {config.get_test_data_path()}")
            print(f"✓ Examples directory: {config.get_examples_path()}")
            
            # Copy package data if user directories are empty
            self.copy_package_data_to_user(user_config_dir)
            
        except ImportError as e:
            print(f"⚠ Could not initialize configuration: {e}")
        
        print("="*50)
    
    def copy_package_data_to_user(self, user_config_dir):
        """Copy Tests and Examples to user directory if empty"""
        package_root = Path(__file__).parent
        
        # Copy Tests if user Tests directory is empty
        user_tests_dir = user_config_dir / 'Tests'
        package_tests_dir = package_root / 'Tests'
        
        if package_tests_dir.exists() and self.is_dir_empty(user_tests_dir):
            print(f"\nCopying Tests to user directory...")
            self.copy_directory(package_tests_dir, user_tests_dir)
            print(f"✓ Copied Tests to: {user_tests_dir}")
        
        # Copy Examples if user Examples directory is empty
        user_examples_dir = user_config_dir / 'Examples'
        package_examples_dir = package_root / 'Examples'
        
        if package_examples_dir.exists() and self.is_dir_empty(user_examples_dir):
            print(f"\nCopying Examples to user directory...")
            self.copy_directory(package_examples_dir, user_examples_dir)
            print(f"✓ Copied Examples to: {user_examples_dir}")
    
    def is_dir_empty(self, dir_path):
        """Check if directory is empty"""
        if not dir_path.exists():
            return True
        return len(list(dir_path.iterdir())) == 0
    
    def copy_directory(self, src, dst):
        """Copy directory contents"""
        for item in src.iterdir():
            if item.is_dir():
                shutil.copytree(item, dst / item.name, dirs_exist_ok=True)
            else:
                shutil.copy2(item, dst / item.name)


class CustomInstallCommand(PostInstallMixin, install):
    """Custom installation command"""
    
    def run(self):
        install.run(self)
        self.setup_user_environment()
        
        print("\n" + "="*50)
        print("Installation Complete!")
        print("="*50)
        print("\nUsage:")
        print("  ooccupy --help                 # Show available commands")
        print("  ooccupy config --show          # Show configuration")
        print("\nYour data is available at:")
        print(f"  {Path.home() / '.ooccupy'}")
        print("="*50)


class CustomDevelopCommand(PostInstallMixin, develop):
    """Custom development installation command"""
    
    def run(self):
        develop.run(self)
        self.setup_user_environment()
        print("\n✓ Development installation complete")


# Read requirements
def read_requirements():
    try:
        with open("requirements.txt", "r") as f:
            return [line.strip() for line in f if line.strip() and not line.startswith('#')]
    except:
        return []

def read_long_description():
    try:
        with open("README.md", "r") as f:
            return f.read()
    except:
        return "Object Oriented Computational Chemistry in Python"


#==================================================================
if __name__ == "__main__":
    # Ensure pyyaml is in requirements
    requirements = read_requirements()
    if not any('yaml' in req.lower() for req in requirements):
        requirements.append('pyyaml>=6.0')
    
    setup(
        name="ooccupy",
        version="1.0.0",
        author="Igor Barden Grillo",
        author_email="barden.igor@gmail.com",
        description="Object Oriented Computational Chemistry in Python",
        long_description=read_long_description(),
        long_description_content_type="text/markdown",
        url="https://github.com/bardenChem/OOCCuPY",
        
        # For flat structure, treat root as package
        packages=find_packages(include=['*']),  # Find all packages
        py_modules=[
            'config',        # Include config.py as module
            'OOCCuPY',       # Include OOCCuPY.py as module  
        ],
        
        # Include data files
        package_data={
            '': [  # Empty string means root package
                'Tests/**/*',
                'Examples/**/*',
                'data/**/*',
                '*.py', '*.txt', '*.md', '*.yaml', '*.json',
            ],
        },
        
        # This is CRITICAL for flat structure
        include_package_data=True,
        
        install_requires=requirements,
        python_requires=">=3.8",
        
        # Entry point points to OOCCuPY.py in root
        entry_points={
            "console_scripts": [
                "ooccupy=OOCCuPY:main",  # Directly from OOCCuPY.py module
            ],
        },
        
        cmdclass={
            'install': CustomInstallCommand,
            'develop': CustomDevelopCommand,
        },
        
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Chemistry",
            "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
        ],
        
        license="MPL-2.0",
    )