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
    
    create_directory('tests/temp')
    create_directory('tests/temp/out')
    create_directory('tests/temp/other')
    create_directory('tests/temp/out/shotgun')
    create_directory('tests/temp/out/amplicons')
    create_directory('tests/temp/out/amplicons_onlyITS')
    create_directory('tests/temp/out/amplicons_only16S')
    create_directory('tests/temp/out/biodiversity')
    
    # Scenario 2 - Amplicons
    args_2 = preparse_args('@tests/amplicons.opts --out_dir tests/temp/out/amplicons')
    opts_2 = bipype.parse_arguments(args_2)

    with keeping_directory_clean('tests/data_do_opts/out/amplicons', move_to='tests/temp/other'):
        bipype.prepare_taxonomy_stats(opts_2)
        
    #Scenario 3 - Amplicons only ITS
    
    args_3 = preparse_args('@tests/amplicons_onlyITS.opts --out_dir tests/temp/out/amplicons_onlyITS')
    opts_3 = bipype.parse_arguments(args_3)

    with keeping_directory_clean('tests/data_do_opts/out/amplicons_onlyITS', move_to='tests/temp/other'):
        bipype.prepare_taxonomy_stats(opts_3)
        
    #Scenario 4 - Amplicons only 16S
    
    args_4 = preparse_args('@tests/amplicons_only16S.opts --out_dir tests/temp/out/amplicons_only16S')
    opts_4 = bipype.parse_arguments(args_4)

    with keeping_directory_clean('tests/data_do_opts/out/amplicons_only16S', move_to='tests/temp/other'):
        bipype.prepare_taxonomy_stats(opts_4)
        
    #Scenario 5 - Biodiversity
    
    args_5 = preparse_args('@tests/biodiversity.opts --out_dir tests/temp/out/biodiversity')
    opts_5 = bipype.parse_arguments(args_5)

    with keeping_directory_clean('tests/data_do_opts/out/biodiversity', move_to='tests/temp/other'):
        bipype.prepare_taxonomy_stats(opts_5)

    # create_directory('tests/temp')

    # How to create test?
    # 1. create input files for testing exculsively prepare_taxonomy_stats()
    #    for example: results of shotgun analysis.
    # 2. create opts files with proper commands to pass in "opts"
    # 3. run test as in the template Scenarios below and
    #    take a look on output files generated by prepare_taxonomy_stats().
    # 4. If an output file looks ok, save it in 'tests' directory,
    #    adding to filename information allowing to
    #    recognise later to whith test it belongs.
    #    For example:  if you have output ITS.krona, then save it as my_test_ITS.krona
    # 4. replace 'test_1.out' with your filename
    # 5. replace 'temp.out' with orginal filename
    # 6. repeat  steps 4-6 for all output files from prepare_taxonomy_stats()
    
    """
    # Scenario 1 - put there some description of scenario
    args = preparse_args('@tests/test_1.opts --out_dir tests/temp')
    opts = bipype.parse_arguments(args)
    with keeping_directory_clean('tests/input_or_some_subdirectory', move_to='some_waste_container'):
        bipype.prepare_taxonomy_stats(opts)
    assert files_identical('tests/temp/temp.out', 'tests/test_1.out')
    assert files_identical('tests/temp/temp_2.out', 'tests/test_1_b.out')
    
    # Scenario 2 - put there some description of scenario
    args = preparse_args('@tests/test_2.opts --out_dir tests/temp')
    opts = bipype.parse_arguments(args)
    with keeping_directory_clean('tests/input', move_to='some_bin'):
        bipype.prepare_taxonomy_stats(opts)
    assert files_identical('tests/temp/temp.out', 'tests/test_2.out')
    """
    
    # remove_entire_directory('tests/temp')
    

def omit_test_sample():
    """
    Tests different scenarios of running sample function
    """
    
    import bipype

    create_directory('tests/temp')
    create_directory('tests/temp/out')
    create_directory('tests/temp/other')
    create_directory('tests/temp/out/shotgun')
    create_directory('tests/temp/out/amplicons')
    create_directory('tests/temp/out/amplicons_onlyITS')
    create_directory('tests/temp/out/amplicons_only16S')
    create_directory('tests/temp/out/biodiversity')
    
    # Scenario 1 - Shotgun
    args_1 = preparse_args('@tests/shotgun.opts --out_dir tests/temp/out/shotgun')
    opts_1 = bipype.parse_arguments(args_1)

    with keeping_directory_clean('tests/data', move_to='tests/temp/other'):
        bipype.sample(opts_1)
        
    # Scenario 2 - Amplicons
    args_2 = preparse_args('@tests/amplicons.opts --out_dir tests/temp/out/amplicons')
    opts_2 = bipype.parse_arguments(args_2)

    with keeping_directory_clean('tests/data', move_to='tests/temp/other'):
        bipype.sample(opts_2)
        
    #Scenario 3 - Amplicons only ITS
    
    args_3 = preparse_args('@tests/amplicons_onlyITS.opts --out_dir tests/temp/out/amplicons_onlyITS')
    opts_3 = bipype.parse_arguments(args_3)

    with keeping_directory_clean('tests/data', move_to='tests/temp/other'):
        bipype.sample(opts_3)
        
    #Scenario 4 - Amplicons only 16S
    
    args_4 = preparse_args('@tests/amplicons_only16S.opts --out_dir tests/temp/out/amplicons_only16S')
    opts_4 = bipype.parse_arguments(args_4)

    with keeping_directory_clean('tests/data', move_to='tests/temp/other'):
        bipype.sample(opts_4)
        
    #Scenario 5 - Biodiversity
    
    args_5 = preparse_args('@tests/biodiversity.opts --out_dir tests/temp/out/biodiversity')
    opts_5 = bipype.parse_arguments(args_5)

    with keeping_directory_clean('tests/data', move_to='tests/temp/other'):
        bipype.sample(opts_5)

    # assert files_identical('tests/temp/temp.out', 'tests/shotgun.out')


    # at the moment we need all ouptut files to "peer review"
    # remove_entire_directory('tests/temp')
