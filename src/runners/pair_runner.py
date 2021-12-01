
import os
import sys
import subprocess as sp

from src.oarsman_arguments import PairArguments
from src.oarsman_dependencies import PairDependencies
from src.output_data import PairOutput

from src.filesystem import rm_fastq_extention, rm_file, check_files_exist


def _configure_pair_command(args: PairArguments, dependencies: PairDependencies):

    command = ' '.join(
        [
            dependencies.seqkit_fpath, 'pair',
            f'-1 {args.reads_R1_fpath}',
            f'-2 {args.reads_R2_fpath}',
            '-u'
        ]
    )

    return command
# end def _configure_pair_command


def run_pair(args, dependencies):

    command_str = _configure_pair_command(args, dependencies)

    pipe = sp.Popen(command_str, shell=True, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'The program seqkit pair returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if

    outdpath = os.path.dirname(args.reads_R1_fpath)

    paired_R1_fpath = os.path.join(
        outdpath,
        rm_fastq_extention(args.reads_R1_fpath) + '.paired.fastq.gz'
    )
    paired_R2_fpath = os.path.join(
        outdpath,
        rm_fastq_extention(args.reads_R2_fpath) + '.paired.fastq.gz'
    )
    unpaired_R1_fpath = os.path.join(
        outdpath,
        rm_fastq_extention(args.reads_R1_fpath) + '.unpaired.fastq.gz'
    )
    unpaired_R2_fpath = os.path.join(
        outdpath,
        rm_fastq_extention(args.reads_R2_fpath) + '.unpaired.fastq.gz'
    )

    # Check if all output files exist
    non_extant_fpaths = check_files_exist(
        paired_R1_fpath,
        paired_R2_fpath
    )
    if len(non_extant_fpaths) != 0:
        for i, fpath in enumerate(non_extant_fpaths):
            print(f'\nError #{i+1}: file `{fpath}` does not exist after `seqkit pair` program has ended.')
            print('This file is the output of the program `seqkit pair`, so it must exist.')
        # end for
        sys.exit(1)
    # end if

    # Add non-empty unpaired files
    unpaired_reads_fpaths = list()
    for unpaired_fpath in (unpaired_R1_fpath, unpaired_R2_fpath):
        if os.path.exists(unpaired_fpath) and os.path.getsize(unpaired_fpath) > 0:
            unpaired_reads_fpaths.append(unpaired_fpath)
        # end if
    # end for

    return PairOutput(
        paired_R1_fpath,
        paired_R2_fpath,
        unpaired_reads_fpaths
    )
# end def run_pair
