#!/usr/bin/env python3
"""
Configuration management for OOCCuPY - Flat Structure Version
"""

import os
import sys
import json
import yaml
from pathlib import Path
from typing import Dict, Any, Optional

class OOCCuPYConfig:
    """Configuration manager for OOCCuPY - Flat Structure"""
    
    def __init__(self, config_dir: Optional[Path] = None):
        if config_dir is None:
            self.config_dir = Path.home() / ".ooccupy"
        else:
            self.config_dir = Path(config_dir)
        
        # Create config directory
        self.config_dir.mkdir(exist_ok=True)
        
        # Define subdirectories - match your setup
        self.dirs = {
            'data': self.config_dir / "data",
            'Tests': self.config_dir / "Tests",      # Capital T to match your structure
            'Examples': self.config_dir / "Examples", # Capital E to match
            'logs': self.config_dir / "logs",
            'cache': self.config_dir / "cache",
            'temp': self.config_dir / "temp",
        }
        
        # Configuration file
        self.config_file = self.config_dir / "config.yaml"
        
        # Default configuration
        self.default_config = {
            'general': {
                'verbose': True,
                'log_level': 'INFO',
                'log_to_file': True,
                'use_cache': True,
                'cache_size': 1000,
            },
            'computational': {
                'default_method': 'DFT',
                'default_basis': '6-31G*',
                'max_memory': '4GB',
                'n_cores': 'auto',
                'scratch_dir': '/tmp/ooccupy',
            },
            'paths': {
                'pdynamo_path': os.environ.get('PDYNAMO_PATH', '/path/to/your/pdynamo'),
                'test_data_path': str(self.dirs['Tests']),
                'examples_path': str(self.dirs['Examples']),
                'user_config_dir': str(self.config_dir),
            },
            'modules': {
                'pDynamoWrapper': {
                    'enabled': True,
                    'test_data_path': str(self.dirs['Tests'] / 'pDynamoWrapper'),
                    'examples_path': str(self.dirs['Examples'] / 'pDynamoWrapper'),
                },
                'QM_inputs': {
                    'enabled': True,
                },
                'MD_tools': {
                    'enabled': True,
                },
                'Structure': {
                    'enabled': True,
                },
            }
        }
        
        # Initialize
        self._create_directories()
        self._load_config()
        self._discover_package_paths()
    
    def _create_directories(self):
        """Create all necessary directories"""
        for dir_name, dir_path in self.dirs.items():
            dir_path.mkdir(exist_ok=True)
            print(f"Created directory: {dir_path}")
    
    def _load_config(self):
        """Load or create configuration"""
        if self.config_file.exists():
            try:
                with open(self.config_file, 'r') as f:
                    self.config = yaml.safe_load(f)
                print(f"Loaded configuration from {self.config_file}")
            except Exception as e:
                print(f"Error loading config: {e}. Using defaults.")
                self.config = self.default_config.copy()
        else:
            self.config = self.default_config.copy()
            self._save_config()
            print(f"Created default configuration at {self.config_file}")
    
    def _save_config(self):
        """Save configuration"""
        try:
            with open(self.config_file, 'w') as f:
                yaml.dump(self.config, f, default_flow_style=False)
        except Exception as e:
            print(f"Error saving config: {e}")
    
    def _discover_package_paths(self):
        """Discover package installation paths"""
        try:
            # Try to find package root
            # Method 1: Check if OOCCuPY is installed as package
            try:
                import OOCCuPY
                package_root = Path(OOCCuPY.__file__).parent
            except ImportError:
                # Method 2: Check common installation locations
                import site
                for site_dir in site.getsitepackages():
                    potential_path = Path(site_dir) / 'OOCCuPY'
                    if potential_path.exists():
                        package_root = potential_path
                        break
                else:
                    # Method 3: Use current directory (for development)
                    package_root = Path(__file__).parent
            
            # Store package root in config
            self.config['paths']['package_root'] = str(package_root)
            
            # Check for package data directories
            for data_type in ['Tests', 'Examples']:
                package_data_path = package_root / data_type
                if package_data_path.exists():
                    self.config['paths'][f'package_{data_type.lower()}'] = str(package_data_path)
                    print(f"Found package {data_type}: {package_data_path}")
            
        except Exception as e:
            print(f"Could not discover package paths: {e}")
    
    def get_package_root(self) -> Optional[Path]:
        """Get package installation root"""
        package_root = self.config.get('paths.package_root')
        if package_root:
            return Path(package_root)
        return None
    
    def get_test_data_path(self, module: str = "") -> Path:
        """Get path to test data"""
        # First check user directory
        user_path = self.dirs['Tests']
        if module:
            user_path = user_path / module
        
        if user_path.exists():
            return user_path
        
        # Fallback to package directory
        package_tests = self.config.get('paths.package_tests')
        if package_tests:
            package_path = Path(package_tests)
            if module:
                package_path = package_path / module
            if package_path.exists():
                return package_path
        
        # Return user path even if it doesn't exist yet
        return user_path
    
    def get_examples_path(self, module: str = "") -> Path:
        """Get path to examples"""
        # First check user directory
        user_path = self.dirs['Examples']
        if module:
            user_path = user_path / module
        
        if user_path.exists():
            return user_path
        
        # Fallback to package directory
        package_examples = self.config.get('paths.package_examples')
        if package_examples:
            package_path = Path(package_examples)
            if module:
                package_path = package_path / module
            if package_path.exists():
                return package_path
        
        # Return user path even if it doesn't exist yet
        return user_path
    
    def get(self, key: str, default: Any = None) -> Any:
        """Get configuration value"""
        keys = key.split('.')
        value = self.config
        
        for k in keys:
            if isinstance(value, dict) and k in value:
                value = value[k]
            else:
                return default
        
        return value
    
    def set(self, key: str, value: Any):
        """Set configuration value"""
        keys = key.split('.')
        config = self.config
        
        for k in keys[:-1]:
            if k not in config:
                config[k] = {}
            config = config[k]
        
        config[keys[-1]] = value
        self._save_config()
    
    def show(self):
        """Display configuration"""
        print("\n" + "="*60)
        print("OOCCuPY Configuration")
        print("="*60)
        
        import yaml
        print(yaml.dump(self.config, default_flow_style=False))
        
        print("Directory Paths:")
        print("="*60)
        for name, path in self.dirs.items():
            print(f"{name:10}: {path}")


# Global instance
_config_instance = None

def get_config(config_dir: Optional[Path] = None) -> OOCCuPYConfig:
    """Get or create global configuration"""
    global _config_instance
    if _config_instance is None:
        _config_instance = OOCCuPYConfig(config_dir)
    return _config_instance

def setup_environment():
    """Setup environment - can be called manually"""
    config = get_config()
    print(f"\nOOCCuPY Environment:")
    print(f"  Config: {config.config_dir}")
    print(f"  Tests: {config.get_test_data_path()}")
    print(f"  Examples: {config.get_examples_path()}")
    return config


# Simple utility for path discovery
def find_data_file(filename: str, module: str = "") -> Optional[Path]:
    """Find a data file searching user then package directories"""
    config = get_config()
    
    # Search locations in order
    search_locations = [
        config.get_test_data_path(module) / filename,
        config.get_examples_path(module) / filename,
    ]
    
    # Add package locations if available
    package_tests = config.get('paths.package_tests')
    package_examples = config.get('paths.package_examples')
    
    if package_tests:
        search_locations.append(Path(package_tests) / module / filename)
        search_locations.append(Path(package_tests) / filename)
    
    if package_examples:
        search_locations.append(Path(package_examples) / module / filename)
        search_locations.append(Path(package_examples) / filename)
    
    for location in search_locations:
        if location.exists():
            return location
    
    return None