#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Logging configuration module for OOCCuPY pDynamoWrapper library.
Provides centralized logging with support for:
- Console output
- File logging
- Separate DEBUG file for verbose messages
- Verbosity levels (CRITICAL, ERROR, WARNING, INFO, DEBUG)
"""

import logging
import os
import sys
from pathlib import Path


class DebugLogger:
    """
    Centralized logging system for OOCCuPY with console and file output.
    """
    
    _instance = None
    _loggers = {}
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(DebugLogger, cls).__new__(cls)
            cls._instance._initialized = False
        return cls._instance
    
    def __init__(self) -> None:
        if self._initialized:
            return
        self._initialized = True
        
        self.debug_file_path = None
        self.enable_debug_file = False
        self.verbosity_level = logging.INFO
        
    @classmethod
    def get_logger(cls, name, output_dir=None, enable_debug_file=False, verbosity="INFO"):
        """
        Get or create a logger for a specific module/class.
        
        Parameters:
        -----------
        name : str
            Logger name (typically __name__ or class name)
        output_dir : str, optional
            Directory to write log files to
        enable_debug_file : bool
            If True, creates a separate DEBUG file with verbose messages
        verbosity : str
            Log level: "CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"
        
        Returns:
        --------
        logging.Logger
            Configured logger instance
        """
        
        instance = cls()
        
        if name in cls._loggers:
            return cls._loggers[name]
        
        # Convert verbosity string to logging level
        level_map = {
            "CRITICAL": logging.CRITICAL,
            "ERROR": logging.ERROR,
            "WARNING": logging.WARNING,
            "INFO": logging.INFO,
            "DEBUG": logging.DEBUG,
        }
        log_level = level_map.get(verbosity.upper(), logging.INFO)
        instance.verbosity_level = log_level
        
        # Create logger
        logger = logging.getLogger(name)
        logger.setLevel(log_level)
        logger.propagate = False
        
        # Clear existing handlers
        logger.handlers.clear()
        
        # Console handler (always enabled)
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(log_level)
        console_formatter = logging.Formatter(
            '[%(name)s] [%(levelname)s] %(message)s'
        )
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
        
        # File handlers (if output_dir is specified)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            
            # Main log file
            log_file = os.path.join(output_dir, f"{name}.log")
            file_handler = logging.FileHandler(log_file, mode='w')
            file_handler.setLevel(log_level)
            file_formatter = logging.Formatter(
                '%(asctime)s [%(name)s] [%(levelname)s] %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
            file_handler.setFormatter(file_formatter)
            logger.addHandler(file_handler)
            
            # Debug file (if enabled) - captures DEBUG and above
            if enable_debug_file:
                debug_file = os.path.join(output_dir, f"{name}_DEBUG.log")
                debug_handler = logging.FileHandler(debug_file, mode='w')
                debug_handler.setLevel(logging.DEBUG)
                debug_formatter = logging.Formatter(
                    '%(asctime)s [%(name)s] [%(levelname)s] [%(funcName)s:%(lineno)d] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S'
                )
                debug_handler.setFormatter(debug_formatter)
                logger.addHandler(debug_handler)
                instance.debug_file_path = debug_file
                instance.enable_debug_file = True
        
        cls._loggers[name] = logger
        return logger


def get_logger(module_name, output_dir=None, enable_debug_file=False, verbosity="INFO"):
    """
    Convenience function to get a logger.
    
    Usage:
    ------
    logger = get_logger(__name__, output_dir="/path/to/logs", enable_debug_file=True)
    logger.info("Starting process")
    logger.debug("Detailed debug information")
    """
    return DebugLogger.get_logger(
        module_name,
        output_dir=output_dir,
        enable_debug_file=enable_debug_file,
        verbosity=verbosity
    )


class FunctionLogger:
    """
    Decorator to automatically log function entry and exit with parameters and return values.
    
    Usage:
    ------
    @FunctionLogger.log_function
    def my_function(param1, param2):
        return result
    """
    
    @staticmethod
    def log_function(func):
        """
        Decorator that logs function entry/exit.
        """
        def wrapper(*args, **kwargs):
            logger = get_logger(func.__module__)
            
            # Log entry
            args_repr = [repr(a) for a in args]
            kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]
            signature = ", ".join(args_repr + kwargs_repr)
            logger.debug(f"ENTERING: {func.__name__}({signature})")
            
            try:
                result = func(*args, **kwargs)
                logger.debug(f"EXITING: {func.__name__} (returned successfully)")
                return result
            except Exception as e:
                logger.error(f"EXCEPTION in {func.__name__}: {type(e).__name__}: {str(e)}")
                raise
        
        return wrapper


# =============================================================================
#  Quick logging templates for common tasks
# =============================================================================

def log_step_start(logger, step_name, **context):
    """
    Log the start of a processing step with context information.
    
    Usage:
    ------
    log_step_start(logger, "NEB Initialization", 
                   spring_constant=500.0, 
                   bins=12)
    """
    context_str = " | ".join(f"{k}={v}" for k, v in context.items())
    if context_str:
        logger.info(f"{'='*70}")
        logger.info(f">>> STARTING: {step_name}")
        logger.info(f"    Context: {context_str}")
        logger.info(f"{'='*70}")
    else:
        logger.info(f">>> STARTING: {step_name}")


def log_step_end(logger, step_name, **results):
    """
    Log the end of a processing step with results.
    
    Usage:
    ------
    log_step_end(logger, "NEB Initialization", 
                 frames_created=12,
                 initial_rms=1.5)
    """
    results_str = " | ".join(f"{k}={v}" for k, v in results.items())
    if results_str:
        logger.info(f"{'='*70}")
        logger.info(f"<<< COMPLETED: {step_name}")
        logger.info(f"    Results: {results_str}")
        logger.info(f"{'='*70}")
    else:
        logger.info(f"<<< COMPLETED: {step_name}")


def log_checkpoint(logger, checkpoint_name, status="OK", message=""):
    """
    Log a checkpoint for debugging state at specific points.
    
    Usage:
    ------
    log_checkpoint(logger, "Trajectory Loading", status="OK", 
                   message="Loaded 25 frames from source")
    """
    status_symbol = "✓" if status == "OK" else "✗" if status == "ERROR" else "⚠"
    msg = f"[{status_symbol}] CHECKPOINT: {checkpoint_name}"
    if message:
        msg += f" - {message}"
    
    if status == "ERROR":
        logger.error(msg)
    elif status == "WARNING":
        logger.warning(msg)
    else:
        logger.info(msg)
