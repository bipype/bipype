from common import *
"""
This is a proposed module for test creation.

It is only scaffold at the moment, without any input and output files.
It's essential to rename bipype to bipype.py to run these tests.

I highly recommend use of py.test to run these tests! 
"""


def test_self_integrity():
    """
    Tests if the most basic parts of bipype exists
    """
    assert module_exists('bipype')
    assert module_exists('refseq_bipype')


def test_paths_from_settings():
    """
    Tests existence of all files defined by non-private,
    prefixed with 'PATH' variables from settings_bipype module.
    """
    import settings_bipype

    namespace = settings_bipype.__dict__
    
    variables = { key: namespace[key]
                  for key in namespace
                  if key.startswith('PATH') }
    
    for var in variables.values():
        assert path_exists(var)


def test_prepare_taxonomy_stats():
    """
    Tests different scenarios of running prepare_taxonomy_stats function
    """
    
    import bipype

    # Scenario 1 - put there some description of scenario
    opts = bipype.parse_arguments('@tests/prepare_taxonomy_stats_1.in')
    bipype.prepare_taxonomy_stats(opts)
    assert files_identical('temp.out', 'tests/prepare_taxonomy_stats_1.out')
    os.remove('temp.out')

    # Scenario 2 - put there some description of scenario
    opts = bipype.parse_arguments('@tests/prepare_taxonomy_stats_2.in')
    bipype.prepare_taxonomy_stats(opts)
    assert files_identical('temp.out', 'tests/prepare_taxonomy_stats_2.out')
    os.remove('temp.out')
    

def test_sample():
    """
    Tests different scenarios of running sample function
    """
    
    import bipype

    # Scenario 1 - put there some description of scenario
    opts = bipype.parse_arguments('@tests/bipype_1.in')
    bipype.sample(opts)
    assert files_identical('temp.out', 'tests/bipype_1.out')
    os.remove('temp.out')

    # Scenario 2 - put there some description of scenario
    opts = bipype.parse_arguments('@tests/bipype_2.in')
    bipype.sample(opts)
    assert files_identical('temp.out', 'tests/bipype_2.out')
    remove_file('temp.out')
