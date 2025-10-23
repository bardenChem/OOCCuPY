#!/usr/bin/env python3
"""
OOCCuPY Installer
Object Oriented Computational Chemistry in Python
"""

import os
import sys
import platform
import subprocess
from pathlib import Path
from setuptools import setup, find_packages
from setuptools.command.install import install

class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        # Run standard installation
        install.run(self)
        
        # Run post-install tasks
        self.post_install()
    
    def post_install(self):
        """Post-installation tasks"""
        print("\n" + "="*50)
        print("OOCCuPY Post-Installation Setup")
        print("="*50)
        
        # Set environment variables
        self.set_environment_variables()
        
        # Create configuration directory
        self.create_config_dir()
        
        # Verify installation
        self.verify_installation()
        
        print("OOCCuPY installation completed successfully!")
        print("You can now use 'ooccupy' command from anywhere.")
        print("="*50)

    def set_environment_variables(self):
        """Set OOCCuPY environment variables"""
        print("\nSetting up environment variables...")
        
        # Get installation paths
        install_path = Path(__file__).parent.absolute()
        script_path = install_path / "OOCCuPY.py"
        
        # Platform-specific environment setup
        system = platform.system().lower()
        
        if system == "windows":
            self.set_windows_environment(install_path, script_path)
        else:  # Linux, macOS, etc.
            self.set_unix_environment(install_path, script_path)

    def set_windows_environment(self, install_path, script_path):
        """Set environment variables on Windows"""
        try:
            # Add to PATH
            subprocess.run([
                'setx', 'OOCCUPY_ROOT', str(install_path)
            ], check=True, shell=True)
            
            # Create batch file for easy access
            batch_content = f'''@echo off
python "{script_path}" %*
'''
            batch_file = install_path / "ooccupy.bat"
            with open(batch_file, 'w') as f:
                f.write(batch_content)
                
            print(f" Created batch file: {batch_file}")
            print(" Added OOCCUPY_ROOT to environment variables")
            
        except subprocess.CalledProcessError as e:
            print(f" Could not set Windows environment variables: {e}")

    def set_unix_environment(self, install_path, script_path):
        """Set environment variables on Unix-like systems"""
        home = Path.home()
        shell_rc = self.get_shell_rc()
        
        if shell_rc:
            env_setup = f'''
# OOCCuPY Environment Variables
export OOCCUPY_ROOT="{install_path}"
export PYTHONPATH="$OOCCUPY_ROOT:$PYTHONPATH"
alias ooccupy="python3 $OOCCUPY_ROOT/OOCCuPY.py"
'''
            
            # Append to shell rc file
            with open(shell_rc, 'a') as f:
                f.write(env_setup)
            
            print(f" Added environment variables to {shell_rc}")
            print(" Created 'ooccupy' alias")
            
            # Make main script executable
            script_path.chmod(0o755)
            print(f" Made {script_path} executable")
        
        # Create symlink in /usr/local/bin (requires sudo)
        try:
            symlink_path = "/usr/local/bin/ooccupy"
            if not Path(symlink_path).exists():
                subprocess.run([
                    'sudo', 'ln', '-sf', str(script_path), symlink_path
                ], check=True)
                print(f" Created symlink: {symlink_path}")
        except (subprocess.CalledProcessError, PermissionError):
            print(" Could not create system-wide symlink (run with sudo for system-wide installation)")

    def get_shell_rc(self):
        """Determine the appropriate shell rc file"""
        shell = os.environ.get('SHELL', '')
        home = Path.home()
        
        if 'zsh' in shell:
            return home / '.zshrc'
        elif 'bash' in shell:
            return home / '.bashrc'
        elif 'fish' in shell:
            return home / '.config/fish/config.fish'
        else:
            # Default to bashrc
            return home / '.bashrc'

    def create_config_dir(self):
        """Create configuration directory"""
        config_dir = Path.home() / '.ooccupy'
        config_dir.mkdir(exist_ok=True)
        
        # Create default config file
        config_file = config_dir / 'config.py'
        if not config_file.exists():
            default_config = '''# OOCCuPY Configuration File

# Default paths
PDYNAMO_PATH = "/path/to/your/pdynamo"  # Update this path
SCRATCH_DIR = "/tmp/ooccupy"

# Default computational settings
DEFAULT_METHOD = "DFT"
DEFAULT_BASIS = "6-31G*"

# Logging settings
VERBOSE = True
LOG_LEVEL = "INFO"
'''
            with open(config_file, 'w') as f:
                f.write(default_config)
            print(f"✓ Created default config: {config_file}")

    def verify_installation(self):
        """Verify the installation was successful"""
        print("\nVerifying installation...")
        try:
            # Test import
            sys.path.insert(0, str(Path(__file__).parent))
            from OOCCuPY import __version__
            print(f" OOCCuPY imported successfully (v{__version__})")
            
            # Test submodule imports
            from OOCCuPY.pDynamoWrapper import Wrapper
            print(" pDynamoWrapper imported successfully")
            
        except ImportError as e:
            print(f"✗ Import test failed: {e}")
            print("  Some features might not work correctly")

# Package metadata
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = fh.read().splitlines()

setup(
    name="ooccupy",
    version="1.0.0",
    author="Igor Barden Grillo",
    author_email="barden.igor@gmail.com",
    description="Object Oriented Computational Chemistry in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bardenChem/OOCCuPY",
    packages=find_packages(),
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
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "ooccupy=OOCCuPY.OOCCuPY:main",
        ],
    },
    license="MPL-2.0",
    # ... rest of setup config
)