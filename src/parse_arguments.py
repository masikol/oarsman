
import os
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

    parser.add_argument(
        '--lofreq',
        help='path to lofreq executable',
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

    oarsman_args = _configure_read_files(argparse_args, oarsman_args)

    # Set primers fpath
    oarsman_args.primers_fpath = os.path.abspath(argparse_args.primers_fpath)
    if not os.path.isfile(oarsman_args.primers_fpath):
        errors.append(
            'File `{}` does not exist'.format(oarsman_args.primers_fpath)
        )
    # end if

    # Set reference fpath
    oarsman_args.reference_fpath = os.path.abspath(argparse_args.reference_genome)
    if not os.path.isfile(oarsman_args.reference_fpath):
        errors.append(
            'File `{}` does not exist'.format(oarsman_args.reference_fpath)
        )
    # end if

    # Set and create output directory
    if not argparse_args.outdir is None:
        oarsman_args.outdir_path = os.path.abspath(argparse_args.outdir)
    # end if
    oarsman_args.create_output_dir()

    # Set advanced kromsatel args
    if not argparse_args.kromsatel_args is None:
        oarsman_args.kromsatel_args = argparse_args.kromsatel_args
    # end if

    oarsman_args.n_threads = _get_n_threads(
        argparse_args,
        oarsman_args.n_threads,
        errors
    )

    return oarsman_args, errors
# end def


def _configure_read_files(argparse_args, oarsman_args):

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

    # # TODO : add long reads support
    # allowed_modes = {'FRl', 'Frl', 'frL'}
    allowed_modes = {'FRl', 'Frl'}
    if not input_mode_str in allowed_modes:
        error_msg = '\nError: the program can work only in one of the following modes:\n' \
            '  "short-single-end" or "short-paired-end" mode.'
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


def _get_n_threads(argparse_args, default_value, errors):

    n_threads_specified = not argparse_args.threads is None
    if not n_threads_specified:
        return default_value
    # end if

    try:
        n_threads = int(argparse_args.threads)
        if n_threads < 1:
            raise ValueError
    except ValueError:
        error_msg = 'Invalid value passed with option `-t/--threads`: {}\n' \
            'This value must be an integer > 0.'.format(argparse_args.threads)
        errors.append(error_msg)
        return None
    else:
        if n_threads > len(os.sched_getaffinity(0)):
            print(
                'Your system has only {} CPU threads available' \
                    .format(len(os.sched_getaffinity(0)))
            )
            print(
                'Switching number of threads from {} to {}' \
                    .format(
                        n_threads,
                        len(os.sched_getaffinity(0))
                    )
            )
            n_threads = len(os.sched_getaffinity(0))
        # end if
    # end try

    return n_threads
# end def


def _configure_oarsman_dependencies(argparse_args):

    oarsman_dependencies = OarsmanDependencies()
    errors = list()

    oarsman_dependencies.kromsatel_fpath = \
        _get_dependency(
            argparse_args.kromsatel,
            'kromsatel.py',
            errors,
            binary=False
        )

    oarsman_dependencies.highlighter_fpath = \
        _get_dependency(
            argparse_args.highlighter,
            'consensus-highlighter.py',
            errors,
            binary=False
        )

    oarsman_dependencies.bwa_fpath = \
        _get_dependency(
            argparse_args.bwa,
            'bwa',
            errors
        )

    oarsman_dependencies.samtools_fpath = \
        _get_dependency(
            argparse_args.samtools,
            'samtools',
            errors
        )

    oarsman_dependencies.bcftools_fpath = \
        _get_dependency(
            argparse_args.bcftools,
            'bcftools',
            errors
        )

    oarsman_dependencies.lofreq_fpath = \
        _get_dependency(
            argparse_args.lofreq,
            'lofreq',
            errors,
            search_in_PATH=False,
            required=False
        )

    _check_blastplus_in_path(oarsman_dependencies, errors)
    _check_var_call_dependencies(oarsman_dependencies, errors)

    return oarsman_dependencies, errors
# end def


def _get_dependency(argparse_dependency_arg,
                    executable_name,
                    errors,
                    binary=True,
                    search_in_PATH=True,
                    required=True):

    dependency_passed_with_cmd = not argparse_dependency_arg is None
    dependency_fpath = argparse_dependency_arg

    if dependency_passed_with_cmd:
        dependency_fpath = os.path.abspath(dependency_fpath)
        if not os.path.exists(dependency_fpath):
            errors.append('File `{}` does not exist'.format(dependency_fpath))
        elif binary and not _is_executable(dependency_fpath):
            errors.append('File `{}` is not executable ' \
                '(please change permissions for it)'.format(dependency_fpath))
        # end if
    elif search_in_PATH:
        dependency_fpath = _check_in_path(executable_name, errors, required)
    # end if

    return dependency_fpath
# end def


def _is_executable(file_path):
    return os.access(file_path, os.X_OK)
# end def


def _check_in_path(executable_name,
                   errors,
                   required=True):
    dependency_fpath = None
    dependency_in_path = fs.util_is_in_path(executable_name)

    if dependency_in_path:
        dependency_fpath = executable_name
    elif required:
        errors.append(
            'Cannot find {} executable in the PATH environment variable' \
                .format(executable_name)
        )
    # end if

    return dependency_fpath
# end def


# def _get_dependency_fpath_binary(argparse_dependency_arg,
#                                  executable_name,
#                                  errors,
#                                  required=True):

#     dependency_passed_with_cmd = not argparse_dependency_arg is None
#     dependency_fpath = argparse_dependency_arg

#     if dependency_passed_with_cmd:
#         dependency_fpath = os.path.abspath(dependency_fpath)
#         if not os.path.exists(dependency_fpath):
#             errors.append('File `{}` does not exist'.format(dependency_fpath))
#         elif not os.access(dependency_fpath, os.X_OK):
#             errors.append('File `{}` is not executable ' \
#                 '(please change permissions for it)'.format(dependency_fpath))
#         # end if
#     else:
#         dependency_in_path = fs.util_is_in_path(executable_name)
#         if dependency_in_path:
#             dependency_fpath = executable_name
#         elif required:
#             errors.append(
#                 'Cannot find {} executable in the PATH environment variable' \
#                     .format(executable_name)
#             )
#         # end if
#     # end if

#     return dependency_fpath
# # end def


def _get_dependency_fpath_script(argparse_dependency_arg,
                                 executable_name,
                                 errors):

    dependency_passed_with_cmd = not argparse_dependency_arg is None
    dependency_fpath = argparse_dependency_arg

    if dependency_passed_with_cmd:
        dependency_fpath = os.path.abspath(dependency_fpath)
        if not os.path.exists(dependency_fpath):
            errors.append('File `{}` does not exist'.format(dependency_fpath))
        # end if
    else:
        dependency_in_path = fs.util_is_in_path(executable_name)
        if dependency_in_path:
            dependency_fpath = executable_name
        else:
            errors.append(
                'Cannot find {} executable in the PATH environment variable' \
                    .format(executable_name)
            )
        # end if
    # end if

    return dependency_fpath
# end def


def _check_blastplus_in_path(oarsman_dependencies, errors):

    executables = (
        oarsman_dependencies.blastn_fpath,
        oarsman_dependencies.makeblastdb_fpath,
    )

    for executable in executables:
        found_in_path = fs.util_is_in_path(executable)
        if not found_in_path:
            errors.append(
                'Cannot find {} executable in the PATH environment variable' \
                    .format(executable)
            )
        # end if
    # end for
# end def


def _check_var_call_dependencies(oarsman_dependencies, errors):

    only_lofreq_available = oarsman_dependencies.bcftools_fpath is None \
                            and oarsman_dependencies.lofreq_fpath is None

    if only_lofreq_available:
        errors.append(
            'Error: lofreq is available, but bcftools is not. Please, install bcftools.'
        )
    # end if
# end if