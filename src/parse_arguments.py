
import os
import sys
import argparse

from src.oarsman_arguments import OarsmanArguments
from src.oarsman_dependencies import OarsmanDependencies

from src.filesystem import util_is_in_path


def parse_arguments():

    parser = argparse.ArgumentParser()

    # Pipeline input files

    parser.add_argument(
        '-p',
        '--primers-fpath',
        help='path to a CSV file of primers names and sequences',
        required=True
    )

    parser.add_argument(
        '-r',
        '--reference-genome',
        help='path to a fasta file of sequence(s) of a reference genome',
        required=True
    )

    # kromsatel arguments

    parser.add_argument(
        '--kromsatel-min-len',
        help='kromsatel minimum amplicon length',
        required=False
    )

    # Dependencies

    parser.add_argument(
        '--kromsatel-dir',
        help='path to kromsatel source directory',
        required=True
    )

    parser.add_argument(
        '--seqkit',
        help='path to seqkit executable',
        required=False
    )

    parser.add_argument(
        '--blastn',
        help='path to blastn executable',
        required=False
    )

    # Misc

    parser.add_argument(
        '--tmp-dir',
        help='path to a temporary directory for the run',
        required=False
    )

    command_line_args = parser.parse_args()


    oarsman_arguments, arguments_errors = _configure_oarsman_arguments(command_line_args)
    oarsman_dependencies, dependencies_errors = _configure_oarsman_dependencies(command_line_args)

    if len(arguments_errors) != 0 or len(dependencies_errors) != 0:
        print('\nErrors:')
        complete_error_list = arguments_errors + dependencies_errors
        for i, error in enumerate(complete_error_list):
            print(f'Error #{i+1}: {error}')
        # end if
    # end if

    return oarsman_arguments, oarsman_dependencies
# end def parse_arguments



def _configure_oarsman_arguments(command_line_args):

    oarsman_arguments = OarsmanArguments()
    errors = list()

    abspath = os.path.abspath

    # Set primers fpath (mandatory)
    oarsman_arguments.primers_fpath = abspath(command_line_args.primers_fpath)
    if not os.path.isfile(oarsman_arguments.primers_fpath):
        errors.append(f'File `{oarsman_arguments.primers_fpath}` does not exist')
    # end if

    # Set genome fpath (mandatory)
    oarsman_arguments.ref_genome_seq_fpath = abspath(command_line_args.reference_genome)
    if not os.path.isfile(oarsman_arguments.ref_genome_seq_fpath):
        errors.append(f'File `{oarsman_arguments.ref_genome_seq_fpath}` does not exist')
    # end if

    # Set temporary directory (optional, has a default value)
    if not command_line_args.tmp_dir is None:
        oarsman_arguments.tmp_dir_path = abspath(command_line_args.tmp_dir)
    # end if
    try:
        if not os.path.isdir(oarsman_arguments.tmp_dir_path):
            os.makedirs(oarsman_arguments.tmp_dir_path)
        # end if
    except OSError as err:
        errors.append(f'Cannot create temporary directory `{oarsman_arguments.tmp_dir_path}`: {err}')
    # end try

    # Set kromsatel min amplicon length (optional, has a default value)
    if not command_line_args.kromsatel_min_len is None:
        oarsman_arguments.min_amplicon_len = command_line_args.kromsatel_min_len
        try:
            oarsman_arguments.min_amplicon_len = int(oarsman_arguments.min_amplicon_len)
            if oarsman_arguments.min_amplicon_len < 1:
                raise ValueError
        except ValueError:
            errors.append(f"""Invalid value passed with option `--kromsatel-min-len`: {oarsman_arguments.min_amplicon_len}
        This value must be an integer > 0.""")
        # end try
    # end if

    return oarsman_arguments, errors
# end def _configure_oarsman_arguments


def _configure_oarsman_dependencies(command_line_args):
    oarsman_dependencies = OarsmanDependencies()
    errors = list()

    abspath = os.path.abspath

    # Set kromsatel source dir (mandatory)
    oarsman_dependencies.kromsatel_dirpath = abspath(command_line_args.kromsatel_dir)
    if not os.path.isdir(oarsman_dependencies.kromsatel_dirpath):
        errors.append(f'Directory `{oarsman_dependencies.kromsatel_dirpath}` does not exist')
    # end if

    # Set seqkit executable path (optional, has a default value)
    oarsman_dependencies.seqkit_fpath = command_line_args.seqkit

    if not oarsman_dependencies.seqkit_fpath is None:
        # This block will be run if a path to seqkit executable is specified in the command line
        if not os.path.exists(oarsman_dependencies.seqkit_fpath):
            errors.append('Cannot find seqkit executable in the PATH environment variable')
        elif not os.access(oarsman_dependencies.seqkit_fpath, os.X_OK):
            errors.append(f"""seqkit file `{oarsman_dependencies.seqkit_fpath}` is not executable
    (please change permissions for it)""")
        # end if
    else:
        # This block will be run if a path to seqkit executable is not specified in the command line
        # So, we will search for it in environment variables
        seqkit_in_path = util_is_in_path('seqkit')
        if not seqkit_in_path:
            errors.append('Cannot find seqkit executable in the PATH environment variable')
        else:
            oarsman_dependencies.seqkit_fpath = 'seqkit'
        # end if
    # end if

    # Set blastn executable path (optional, has a default value)
    oarsman_dependencies.blastn_fpath = 'blastn' # it must be in PATH

    # So, we will search for blastn in environment variables
    blastn_in_path = util_is_in_path(oarsman_dependencies.blastn_fpath)
    if not blastn_in_path:
        errors.append('Cannot find blastn executable in the PATH environment variable')
    # end if

    return oarsman_dependencies, errors
# end def _configure_oarsman_dependencies
