
import os
import sys
import argparse

from src.oarsman_arguments import OarsmanArguments
from src.oarsman_dependencies import OarsmanDependencies

from src.filesystem import util_is_in_path, check_files_exist


def parse_arguments():

    parser = argparse.ArgumentParser()

    # Pipeline input files

    parser.add_argument(
        '-1',
        '--reads-R1',
        help='path to a fastq file of forward (R1) reads or unpaired reads.',
        required=True,
        action='append',
        nargs='+'
    )

    parser.add_argument(
        '-2',
        '--reads-R2',
        help='path to a fastq file of reverse (R2) reads.',
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
        '--kromsatel-mi',
        help='minimum length of a minor amplicon for kromsatel',
        required=False
    )

    parser.add_argument(
        '--kromsatel-ma',
        help='minimum length of a minor amplicon for kromsatel',
        required=False
    )

    parser.add_argument(
        '--kromsatel-chunk',
        help='chunk size for kromsatel',
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

    parser.add_argument(
        '-t',
        '--threads',
        help='number of CPU threads to use',
        required=False
    )

    command_line_args = parser.parse_args()


    oarsman_args, arguments_errors = _configure_oarsman_args(command_line_args)
    oarsman_dependencies, dependencies_errors = _configure_oarsman_dependencies(command_line_args)

    if len(arguments_errors) != 0 or len(dependencies_errors) != 0:
        print('\nErrors:')
        complete_error_list = arguments_errors + dependencies_errors
        for i, error in enumerate(complete_error_list):
            print(f'\nError #{i+1}: {error}')
        # end if
        sys.exit(1)
    # end if

    return oarsman_args, oarsman_dependencies
# end def parse_arguments



def _configure_oarsman_args(command_line_args):

    oarsman_args = OarsmanArguments()
    errors = list()

    # We will call this function often, so get rid of dots for speed
    abspath = os.path.abspath

    # Set forward reads (R1) fpath (mandatory)
    # Iterate over samples -- each sample has its own collection of input files with reads
    # `i` is a sample index
    oarsman_args.reads_R1_fpaths = command_line_args.reads_R1
    for i, sample_read_fpaths in enumerate(oarsman_args.reads_R1_fpaths):
        # Make paths absolute
        oarsman_args.reads_R1_fpaths[i] = list(
            map(
                abspath,
                sample_read_fpaths
            )
        )
        # Check existance
        non_extant_files = check_files_exist(
            *oarsman_args.reads_R1_fpaths[i]
        )
        if len(non_extant_files) != 0:
            for fpath in non_extant_files:
                errors.append(f'File `{fpath}` does not exist')
            # end for
        # end if
    # end for

    # Set reverse reads (R2) fpath (optional, might be empty)
    if not command_line_args.reads_R2 is None:
        oarsman_args.reads_R2_fpaths = command_line_args.reads_R2
        # Iterate over samples -- each sample has its own collection of input files with reads
        # `i` is a sample index
        for i, sample_read_fpaths in enumerate(oarsman_args.reads_R2_fpaths):
            # Make paths absolute
            oarsman_args.reads_R2_fpaths[i] = list(
                map(
                    abspath,
                    sample_read_fpaths
                )
            )
            # Check existance
            non_extant_files = check_files_exist(
                *oarsman_args.reads_R2_fpaths[i]
            )
            if len(non_extant_files) != 0:
                for fpath in non_extant_files:
                    errors.append(f'File `{fpath}` does not exist')
                # end for
            # end if
        # end for

        # Check if number of forward reads is equal to number of reverse reads
        # (one cannot mix paired-end and unpaired libraries)
        if len(oarsman_args.reads_R1_fpaths) != len(oarsman_args.reads_R2_fpaths):

            error_msg = """The number of samples having forward reads ({} samples) is not equal to
    the number of samples having reverse reads ({} samples).
Tip: you cannot mix paired-end and unpaired libraries during a single oarsman run.""" \
                .format(
                    len(oarsman_args.reads_R1_fpaths),
                    len(oarsman_args.reads_R2_fpaths)
                )
            errors.append(error_msg)
        else:
            forw_revr_zip = zip(
                oarsman_args.reads_R1_fpaths,
                oarsman_args.reads_R2_fpaths
            )
            for i, (forw_fpaths, revr_fpaths) in enumerate(forw_revr_zip):
                if len(forw_fpaths) != len(revr_fpaths):
                    forw_basenames = tuple(map(os.path.basename, forw_fpaths))
                    revr_basenames = tuple(map(os.path.basename, revr_fpaths))

                    error_msg = """Invalid number of input files with reads for the sample
    with ordinal number {}.
The number of "forward" files must be equal to the number of "reverse" files
Forward-read files ({} files): {}.
Reverse-read files ({} files): {}.""".format(
                        i+1,
                        len(forw_fpaths), ', '.join(forw_basenames),
                        len(revr_fpaths), ', '.join(revr_basenames)
                    )

                    errors.append(error_msg)
                # end if
            # end for
        # end if
    # end if

    # Set primers fpath (mandatory)
    oarsman_args.primers_fpath = abspath(command_line_args.primers_fpath)
    if not os.path.isfile(oarsman_args.primers_fpath):
        errors.append(f'File `{oarsman_args.primers_fpath}` does not exist')
    # end if

    # Set genome fpath (mandatory)
    oarsman_args.ref_genome_seq_fpath = abspath(command_line_args.reference_genome)
    if not os.path.isfile(oarsman_args.ref_genome_seq_fpath):
        errors.append(f'File `{oarsman_args.ref_genome_seq_fpath}` does not exist')
    # end if

    # Set temporary directory (optional, has a default value)
    if not command_line_args.tmp_dir is None:
        oarsman_args.tmp_dir_path = abspath(command_line_args.tmp_dir)
    # end if
    try:
        if not os.path.isdir(oarsman_args.tmp_dir_path):
            os.makedirs(oarsman_args.tmp_dir_path)
        # end if
    except OSError as err:
        errors.append('Cannot create temporary directory `{}`: {}' \
            .format(oarsman_args.tmp_dir_path, err)
        )
    # end try

    # Set minimum length for a major amplicon for kromsatel (optional, has a default value)
    if not command_line_args.kromsatel_ma is None:
        oarsman_args.min_major_len = command_line_args.kromsatel_ma
        try:
            oarsman_args.min_major_len = int(oarsman_args.min_major_len)
            if oarsman_args.min_major_len < 1:
                raise ValueError
        except ValueError:
            errors.append("""Invalid value passed with option `--kromsatel-ma`: {}
        This value must be an integer > 0.""".format(oarsman_args.min_major_len))
        # end try
    # end if

    # Set minimum length for a minor amplicon for kromsatel (optional, has a default value)
    if not command_line_args.kromsatel_mi is None:
        oarsman_args.min_minor_len = command_line_args.kromsatel_mi
        try:
            oarsman_args.min_minor_len = int(oarsman_args.min_minor_len)
            if oarsman_args.min_minor_len < 1:
                raise ValueError
        except ValueError:
            errors.append("""Invalid value passed with option `--kromsatel-mi`: {}
        This value must be an integer > 0.""".format(oarsman_args.min_minor_len))
        # end try
    # end if

    # Set chunk size for kromsatel (optional, has a default value)
    if not command_line_args.kromsatel_chunk is None:
        oarsman_args.chunk_size = command_line_args.kromsatel_chunk
        try:
            oarsman_args.chunk_size = int(oarsman_args.chunk_size)
            if oarsman_args.chunk_size < 1:
                raise ValueError
        except ValueError:
            errors.append("""Invalid value passed with option `--kromsatel-chunk`: {}
        This value must be an integer > 0.""".format(oarsman_args.chunk_size))
        # end try
    # end if

    # Set number of CPU threads to use (optional, has a default value)
    if not command_line_args.threads is None:
        oarsman_args.n_threads = command_line_args.threads
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
                    .format(oarsman_args.n_threads), len(os.sched_getaffinity(0)))
                oarsman_args.n_threads = len(os.sched_getaffinity(0))
            # end if
        # end try
    # end if
    

    return oarsman_args, errors
# end def _configure_oarsman_args


def _configure_oarsman_dependencies(command_line_args):
    oarsman_dependencies = OarsmanDependencies()
    errors = list()

    # We will call this function often, so get rid of dots for speed
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
    # So, we will search for blastn in environment variables
    blastn_in_path = util_is_in_path(oarsman_dependencies.blastn_fpath)
    if not blastn_in_path:
        errors.append('Cannot find blastn executable in the PATH environment variable')
    # end if

    # Set makeblastdb executable path (optional, has a default value)
    # So, we will search for makeblastdb in environment variables
    makeblastdb_in_path = util_is_in_path(oarsman_dependencies.makeblastdb_fpath)
    if not makeblastdb_in_path:
        errors.append('Cannot find makeblastdb executable in the PATH environment variable')
    # end if

    return oarsman_dependencies, errors
# end def _configure_oarsman_dependencies
