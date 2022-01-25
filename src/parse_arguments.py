
import os
import sys
import argparse

import src.filesystem as fs
from src.printing import print_err
from src.input_modes import InputModes
from src.fatal_errors import FatalError
from src.arguments import OarsmanArguments
from src.dependencies import OarsmanDependencies


def parse_arguments():

    parser = argparse.ArgumentParser()

    # Pipeline input files

    parser.add_argument(
        '-1',
        '--reads-R1',
        help='path to a fastq file of short forward (R1) or unpaired reads.',
        required=False,
        action='append',
        nargs='+'
    )

    parser.add_argument(
        '-2',
        '--reads-R2',
        help='path to a fastq file of short reverse (R2) reads.',
        required=False,
        action='append',
        nargs='*'
    )

    parser.add_argument(
        '-l',
        '--reads-long',
        help='path to a fastq file of reverse long reads.',
        required=False,
        action='append',
        nargs='*'
    )

    parser.add_argument(
        '-r',
        '--reference-genome',
        help='path to a fasta file of sequence(s) of a reference genome',
        required=True
    )

    parser.add_argument(
        '-p',
        '--primers-fpath',
        help='path to a CSV file of primers names and sequences',
        required=True
    )

    # kromsatel arguments

    parser.add_argument(
        '--kromsatel-args',
        help='argument string for kromsatel',
        required=False
    )

    # Dependencies

    parser.add_argument(
        '--kromsatel',
        help='path to the kromsatel script ("kromsatel.py")',
        required=False
    )

    parser.add_argument(
        '--highlighter',
        help='path to the consensus-highlighter script ("consensus-highlighter.py")',
        required=False
    )

    parser.add_argument(
        '--bwa',
        help='path to bwa executable',
        required=False
    )

    # parser.add_argument(
    #     '--bowtie2',
    #     help='path to bowtie2 executable',
    #     required=False
    # )

    parser.add_argument(
        '--samtools',
        help='path to samtools executable',
        required=False
    )

    parser.add_argument(
        '--bcftools',
        help='path to bcftools executable',
        required=False
    )

    # Misc

    parser.add_argument(
        '-o',
        '--outdir',
        help='output directory',
        required=False
    )

    parser.add_argument(
        '-t',
        '--threads',
        help='number of CPU threads to use',
        required=False
    )

    argparse_args = parser.parse_args()


    oarsman_args, arguments_errors = _configure_oarsman_args(argparse_args)
    oarsman_dependencies, dependencies_errors = _configure_oarsman_dependencies(argparse_args)

    if len(arguments_errors) != 0 or len(dependencies_errors) != 0:
        print_err('\nErrors:')
        complete_error_list = arguments_errors + dependencies_errors
        for i, error in enumerate(complete_error_list):
            print_err(f'\nError #{i+1}: {error}')
        # end if
        raise FatalError
    # end if

    return oarsman_args, oarsman_dependencies
# end def


def _configure_oarsman_args(argparse_args):

    oarsman_args = OarsmanArguments()
    errors = list()

    oarsman_args = _configure_read_files(argparse_args, oarsman_args, errors)

    # Set primers fpath (mandatory)
    oarsman_args.primers_fpath = os.path.abspath(argparse_args.primers_fpath)
    if not os.path.isfile(oarsman_args.primers_fpath):
        errors.append(f'File `{oarsman_args.primers_fpath}` does not exist')
    # end if

    # Set genome fpath (mandatory)
    oarsman_args.reference_fpath = os.path.abspath(argparse_args.reference_genome)
    if not os.path.isfile(oarsman_args.reference_fpath):
        errors.append(f'File `{oarsman_args.reference_fpath}` does not exist')
    # end if

    # Set temporary directory (optional, has a default value)
    if not argparse_args.outdir is None:
        oarsman_args.outdir_path = os.path.abspath(argparse_args.outdir)
    # end if
    try:
        if not os.path.isdir(oarsman_args.outdir_path):
            os.makedirs(oarsman_args.outdir_path)
        # end if
    except OSError as err:
        errors.append('Cannot create temporary directory `{}`: {}' \
            .format(oarsman_args.outdir_path, err)
        )
    # end try

    # Set advanced kromsatel args
    if not argparse_args.kromsatel_args is None:
        oarsman_args.kromsatel_args = argparse_args.kromsatel_args
    # end if

    # Set number of CPU threads to use (optional, has a default value)
    if not argparse_args.threads is None:
        oarsman_args.n_threads = argparse_args.threads
        try:
            oarsman_args.n_threads = int(oarsman_args.n_threads)
            if oarsman_args.n_threads < 1:
                raise ValueError
        except ValueError:
            errors.append("""Invalid value passed with option `-t/--threads`: {}
        This value must be an integer > 0.""".format(oarsman_args.n_threads))
        else:
            if oarsman_args.n_threads > len(os.sched_getaffinity(0)):
                print(f'Your system has only {len(os.sched_getaffinity(0))} CPU threads available')
                print('Switching number of threads from {} to {}' \
                    .format(oarsman_args.n_threads), len(os.shed_getaffinity(0)))
                oarsman_args.n_threads = len(os.sched_getaffinity(0))
            # end if
        # end try
    # end if

    return oarsman_args, errors
# end def


def _configure_read_files(argparse_args, oarsman_args, errors):

    input_mode_str = _make_input_mode_str(argparse_args)
    _check_input_mode(input_mode_str)

    if input_mode_str[0] == 'F':
        oarsman_args.reads_R1_fpaths = \
            _configure_fpath_collections(argparse_args.reads_R1)
        _check_input_fpaths(oarsman_args.reads_R1_fpaths)
        oarsman_args.input_mode = InputModes.IlluminaSE

        if input_mode_str == 'FRl':
            oarsman_args.reads_R2_fpaths = \
                _configure_fpath_collections(argparse_args.reads_R2)
            _check_input_fpaths(oarsman_args.reads_R2_fpaths)
            _check_paired_reads(
                oarsman_args.reads_R1_fpaths,
                oarsman_args.reads_R2_fpaths,
            )
            oarsman_args.input_mode = InputModes.IlluminaPE
        # end if

    elif input_mode_str == 'frL':
        oarsman_args.reads_long_fpaths = \
            _configure_fpath_collections(argparse_args.reads_long)
        _check_input_fpaths(oarsman_args.reads_long_fpaths)
        oarsman_args.input_mode = InputModes.Nanopore
    # end if

    return oarsman_args
# end def


def _configure_fpath_collections(fpath_collections):
    return list(
        map(
            fs.abspath_collection,
            fpath_collections
        )
    )
# end def


def _check_input_fpaths(oarsman_collection):

    errors = list()
    for sample_read_fpaths in oarsman_collection:
        errors += _check_files_existance(sample_read_fpaths)
    # end for

    if len(errors) != 0:
        error_msg = '\nError: the following files do not exist:\n' \
            + '\n'.join(
                map(lambda s: '  `{}`'.format(s), errors)
            )
        raise FatalError(error_msg)
    # end if
# end def


def _check_files_existance(file_path_collection):

    errors = list()

    # Check existance
    non_extant_files = fs.check_files_exist(
        *file_path_collection
    )
    if len(non_extant_files) != 0:
        for fpath in non_extant_files:
            errors.append(f'File `{fpath}` does not exist')
        # end for
    # end if
    return errors
# end def


def _check_paired_reads(frw_fpath_collections, rvr_fpath_collections):

    if len(frw_fpath_collections) != len(rvr_fpath_collections):
        error_msg = '\nError: the number of samples having forward reads ({} samples)\n' \
            ' is not equal to the number of samples having reverse reads ({} samples).' \
            'Tip: you cannot mix paired-end and unpaired libraries during a single oarsman run.'
        raise FatalError(error_msg)
    # end if

    frw_rvr_zip = zip(
        frw_fpath_collections,
        rvr_fpath_collections
    )
    sample_number = 1
    for frw_fpaths, rvr_fpaths in frw_rvr_zip:
        if len(frw_fpaths) != len(rvr_fpaths):
            error_msg = 'Invalid number of input files of reads for the sample ' \
                'with ordinal number {}.\n' \
                'Forward-read files ({} files): {}.\n' \
                'Reverse-read files ({} files): {}.\n' \
                'The number of "forward" files must be equal' \
                ' to the number of "reverse" files.' \
                .format(
                    sample_number,
                    len(frw_fpaths), ', '.join(fs.basename_collection(frw_fpaths)),
                    len(rvr_fpaths), ', '.join(fs.basename_collection(rvr_fpaths))
                )
            raise FatalError(error_msg)
        # end if
        sample_number += 1
    # end for
# end def


def _check_input_mode(input_mode_str):

    no_data_mode = 'frl'
    if input_mode_str == no_data_mode:
        error_msg = '\nError: no input data passed to the program.'
        raise FatalError(error_msg)
    # end if

    allowed_modes = {'FRl', 'Frl', 'frL'}
    if not input_mode_str in allowed_modes:
        error_msg = '\nError: the program can work only in single mode:\n' \
            '  either in "short-single-end", "short-paired-end" or "long" mode.'
        raise FatalError(error_msg)
    # end if
# end def


def _make_input_mode_str(argparse_args):
    frw_char  = 'f' if argparse_args.reads_R1   is None else 'F'
    rvr_char  = 'r' if argparse_args.reads_R2   is None else 'R'
    long_char = 'l' if argparse_args.reads_long is None else 'L'

    input_mode_str = '{}{}{}'.format(frw_char, rvr_char, long_char)

    return input_mode_str
# end def


def _configure_oarsman_dependencies(argparse_args):
    oarsman_dependencies = OarsmanDependencies()
    errors = list()

    # We will call this function often, so get rid of dots
    abspath = os.path.abspath

    # Set kromsatel location (mandatory)
    if not argparse_args.kromsatel is None:
        oarsman_dependencies.kromsatel_fpath = abspath(argparse_args.kromsatel)
        if not os.path.isfile(oarsman_dependencies.kromsatel_fpath):
            errors.append(f'File `{oarsman_dependencies.kromsatel_fpath}` does not exist')
        # end if
    else:
        # This block will be run if a path to kromsatel.py executable is not specified in the command line
        # So, we will search for it in environment variables
        kromsatel_in_path = fs.util_is_in_path('kromsatel.py')
        if not kromsatel_in_path:
            errors.append('Cannot find kromsatel.py executable in the PATH environment variable')
        else:
            oarsman_dependencies.kromsatel_fpath = 'kromsatel.py'
        # end if
    # end if


    # Set consensus-highlighter script (mandatory)
    if not argparse_args.highlighter is None:
        oarsman_dependencies.highlighter_fpath = abspath(argparse_args.highlighter)
        if not os.path.isfile(highlighter_fpath):
            errors.append('File `{}` does not exist' \
                .format(oarsman_dependencies.highlighter_fpath))
        # end if
    else:
        # This block will be run if a path to kromsatel.py executable is not specified in the command line
        # So, we will search for it in environment variables
        highlighter_in_path = fs.util_is_in_path('consensus-highlighter.py')
        if not highlighter_in_path:
            errors.append(
                'Cannot find consensus-highlighter.py executable' \
                ' in the PATH environment variable'
            )
        else:
            oarsman_dependencies.highlighter_fpath = 'consensus-highlighter.py'
        # end if
    # end if

    # Set bwa executable path (optional, has a default value)
    oarsman_dependencies.bwa_fpath = argparse_args.bwa

    if not oarsman_dependencies.bwa_fpath is None:
        # This block will be run if a path to bwa executable is specified in the command line
        if not os.path.exists(oarsman_dependencies.bwa_fpath):
            errors.append(f'File `{oarsman_dependencies.bwa_fpath}` does not exist')
        elif not os.access(oarsman_dependencies.bwa_fpath, os.X_OK):
            errors.append(f"""bwa file `{oarsman_dependencies.bwa_fpath}` is not executable
    (please change permissions for it)""")
        # end if
    else:
        # This block will be run if a path to bwa executable is not specified in the command line
        # So, we will search for it in environment variables
        bwa_in_path = fs.util_is_in_path('bwa')
        if not bwa_in_path:
            errors.append('Cannot find bwa executable in the PATH environment variable')
        else:
            oarsman_dependencies.bwa_fpath = 'bwa'
        # end if
    # end if

    # # Set bowtie2 executable path (optional, has a default value)
    # oarsman_dependencies.bowtie2_fpath = argparse_args.bowtie2

    # if not oarsman_dependencies.bowtie2_fpath is None:
    #     # This block will be run if a path to bowtie2 executable is specified in the command line
    #     bowtie2_binaries = (
    #         oarsman_dependencies.bowtie2_fpath,
    #         oarsman_dependencies.bowtie2_fpath+'-build'
    #     )
    #     for util_fpath in bowtie2_binaries:
    #         if not os.path.exists(util_fpath):
    #             errors.append(f'File `{util_fpath}` does not exist')
    #         elif not os.access(util_fpath, os.X_OK):
    #             errors.append(f"""File `{util_fpath}` is not executable
    #     (please change permissions for it)""")
    #         # end if
    #     # end for
    # else:
    #     # This block will be run if a path to bowtie2 executable is not specified in the command line
    #     # So, we will search for it in environment variables
    #     for util_name in ('bowtie2', 'bowtie2-build'):
    #         util_in_path = fs.util_is_in_path(util_name)
    #         if not util_in_path:
    #             errors.append(f'Cannot find {util_name} executable in the PATH environment variable')
    #         else:
    #             oarsman_dependencies.bowtie2_fpath = 'bowtie2'
    #         # end if
    #     # end for
    # # end if


    # Set samtools executable path (optional, has a default value)
    oarsman_dependencies.samtools_fpath = argparse_args.samtools

    if not oarsman_dependencies.samtools_fpath is None:
        # This block will be run if a path to samtools executable is specified in the command line
        if not os.path.exists(oarsman_dependencies.samtools_fpath):
            errors.append(f'File `{oarsman_dependencies.samtools_fpath}` does not exist')
        elif not os.access(oarsman_dependencies.samtools_fpath, os.X_OK):
            errors.append(f"""samtools file `{oarsman_dependencies.samtools_fpath}` is not executable
    (please change permissions for it)""")
        # end if
    else:
        # This block will be run if a path to samtools executable is not specified in the command line
        # So, we will search for it in environment variables
        samtools_in_path = fs.util_is_in_path('samtools')
        if not samtools_in_path:
            errors.append('Cannot find samtools executable in the PATH environment variable')
        else:
            oarsman_dependencies.samtools_fpath = 'samtools'
        # end if
    # end if


    # Set bcftools executable path (optional, has a default value)
    oarsman_dependencies.bcftools_fpath = argparse_args.bcftools

    if not oarsman_dependencies.bcftools_fpath is None:
        # This block will be run if a path to samtools executable is specified in the command line
        if not os.path.exists(oarsman_dependencies.bcftools_fpath):
            errors.append(f'File `{oarsman_dependencies.bcftools_fpath}` does not exist')
        elif not os.access(oarsman_dependencies.bcftools_fpath, os.X_OK):
            errors.append(f"""bcftools file `{oarsman_dependencies.samtools_fpath}` is not executable
    (please change permissions for it)""")
        # end if
    else:
        # This block will be run if a path to samtools executable is not specified in the command line
        # So, we will search for it in environment variables
        bcftools_in_path = fs.util_is_in_path('bcftools')
        if not bcftools_in_path:
            errors.append('Cannot find bcftools executable in the PATH environment variable')
        else:
            oarsman_dependencies.bcftools_fpath = 'bcftools'
        # end if
    # end if


    # Set blastn executable path (optional, has a default value)
    # So, we will search for blastn in environment variables
    blastn_in_path = fs.util_is_in_path(oarsman_dependencies.blastn_fpath)
    if not blastn_in_path:
        errors.append('Cannot find blastn executable in the PATH environment variable')
    # end if

    # Set makeblastdb executable path (optional, has a default value)
    # So, we will search for makeblastdb in environment variables
    makeblastdb_in_path = fs.util_is_in_path(oarsman_dependencies.makeblastdb_fpath)
    if not makeblastdb_in_path:
        errors.append('Cannot find makeblastdb executable in the PATH environment variable')
    # end if

    return oarsman_dependencies, errors
# end def
