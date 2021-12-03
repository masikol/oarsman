
import os
import sys
import subprocess as sp

from src.oarsman_arguments import ReadMappingArguments
from src.oarsman_dependencies import BwaDependencies
from src.output_data import ReadMappingOutput


def _configure_bwa_index_command(
        args: ReadMappingArguments,
        dependencies: BwaDependencies
    ):

    command = ' '.join(
        [
            dependencies.bwa_fpath, 'index',
            f'-p {args.genome_index_base_fpath}',
            args.ref_genome_fpath
        ]
    )

    return command
# end def _configure_bwa_index_command


def _configure_bwa_command(
        args: ReadMappingArguments,
        dependencies: BwaDependencies,
        sam_outfpath: str
    ):

    if not args.reads_R1_fpath is None and not args.reads_R2_fpath is None:
        bwa_read_arguments = ' '.join([args.reads_R1_fpath, args.reads_R2_fpath])
    elif len(args.unpaired_reads_fpaths) == 1:
        bwa_read_arguments = args.unpaired_reads_fpaths[0]
    else:
        print('Error: invalied configuration of input files for bwa!')
        print('Please, contact the developer.')
        sys.exit(1)
    # end if

    command = ' '.join(
        [
            dependencies.bwa_fpath, 'mem',
            f'-t {args.n_threads}',
            f'-o {sam_outfpath}',
            f'{args.genome_index_base_fpath}',
            bwa_read_arguments,
        ]
    )

    return command
# end def _configure_bwa_command


def run_bwa(args, dependencies):

    if not os.path.isdir(args.outdir_path):
        try:
            os.makedirs(args.outdir_path)
        except OSError as err:
            print(f'\nError: cannot create directory `{args.outdir_path}`')
            print(str(err))
            sys.exit(1)
        # end try
    # end if

    sam_outfpath = os.path.join(
        args.outdir_path,
        args.sample_name + '.sam'
    )

    # Create index of reference fasta file
    command_str = _configure_bwa_index_command(args, dependencies)

    print('Building index of the reference genome...')
    pipe = sp.Popen(command_str, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'Program `bwa-build` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if

    # fasta_index_extention = '.fai'
    # ref_genome_index_fpath = args.ref_genome_fpath + fasta_index_extention

    # if not os.path.exists(ref_genome_index_fpath):
    #     print(f"""\nError: the index file `{ref_genome_index_fpath}`
    # of the genome sequence file `{args.ref_genome_fpath}` does not exist after indexing""")
    #     print('This file must exist, though. Exitting...')
    #     sys.exit(1)
    # # end if


    # Perform mapping
    command_str = _configure_bwa_command(args, dependencies, sam_outfpath)

    print('Mapping the reads...')
    # print(command_str)
    pipe = sp.Popen(command_str, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'Program `bwa` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if

    if not os.path.exists(sam_outfpath):
        print(f'\nError: SAM file `{sam_outfpath}` does not exist after bwa has mapped the reads')
        print('This file must exist, though. Exitting...')
        sys.exit(1)
    # end if

    return ReadMappingOutput(sam_outfpath)
# end def run_bwa