import os
import subprocess
import shutil 
from contextlib import contextmanager


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


def get_args_from_file(path):
    args = ''
    with open(path, 'r') as args_file:
        args = args_file.read()
    return args.split()


def remove_file(filename):
    return os.remove(filename)


def create_directory(path):
    """
    Creates directory of name given by 'path',
    but only if there is no (file or directory) with this name.
    """
    if not path_exists(path):
        return os.mkdir(path)
    else:
        print "Warning: path %s already exists" % path


def remove_entire_directory(path):
    """
    Delete an >>entire<< directory tree
    """
    return shutil.rmtree(path)


def pre_parse_args(args_string):
    """
    Simulates argparse-like arguments handling,
    but does _not_ require running from command line. 
    """
    args = args_string.split()
    pre_parsed_args = []
    for arg in args:
        if arg.startswith('@'):
            pre_parsed_args += get_args_from_file(arg[1:])
        else:
            pre_parsed_args += [arg]
    return pre_parsed_args


def copy_files_to_cwd(path, list_of_files):
    """
    Copies given files to cwd.
    """
    cwd = os.getcwd()

    for item in list_of_files:
        source = os.path.join(path, item)
        destination = os.path.join(cwd, item)
        shutil.copy(source, destination)


def remove_files_from_cwd(list_of_files):
    """
    Removes files from given list from current working directory
    """
    cwd = os.getcwd()
    
    for item in list_of_files:
        full_path = os.path.join(cwd, item)
        remove_file(full_path)


def get_files_list(path):
    """
    Gets all >files< from directory specified by 'path'
    """
    list_of_files = []
    for item in os.listdir(path):
        full_path = os.path.join(path, item)
        if os.path.isfile(full_path):
            list_of_files += [item]

    return list_of_files


@contextmanager
def copy_of_folder_in_cwd(path):
    """
    A wrapper designed to use with 'with' statement,
    to allow easy running bipype on files from other place,
    by creation a copy in current working directory,
    and then cleaning by removing these files.
    """
    list_of_files = get_files_list(path)
    
    copy_files_to_cwd(path, list_of_files)
    yield
    remove_files_from_cwd(list_of_files)
