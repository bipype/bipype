import os
import subprocess
import shutil 


def files_identical(filename_1, filename_2):
    """
    Compare two files with diff
    """
    return_code = subprocess.call(['diff', filename_1, filename_2])
    # change from unix convention (0 = True) to Python convention 
    return not return_code


def path_exists(path):
    """
    Check if given path exists.
    """
    return os.path.exists(path)


def module_exists(module_name):
    """
    Check if module of module_name exists.
    """
    import importlib

    try:
        importlib.import_module(module_name)
        return True
    except ImportError:
        return False


def remove_file(filename):
    return os.remove(filename)


def remove_entire_directory(path):
    """
    Delete an >>entire<< directory tree
    """
    return shutil.rmtree(path)


def create_directory(path):
    return os.mkdir(path)
