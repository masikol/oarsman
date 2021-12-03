
import os
import sys
import subprocess as sp

from src.oarsman_arguments import ReadMappingArguments
from src.oarsman_dependencies import Bowtie2Dependencies
from src.output_data import ReadMappingOutput


def _configure_bow1tie2_build_command(
        args: ReadMappingArguments,
        dependencies: Bowtie2Dependencies
    ):

    command = ' '.join(
        [
            dependencies.bowtie2_build_fpath,
            args.ref_genome_fpath,
            args.genome_index_base_fpath,
        ]
    )

    return command
# end def _configure_bow1tie2_build_command


def _configure_bow1tie2_command(
        args: ReadMappingArguments,
        dependencies: Bowtie2Dependencies,
        sam_outfpath: str
    ):

    forw_reads_option = '' if args.reads_R1_fpath is None else f'-1 {args.reads_R1_fpath}'
    revr_reads_option = '' if args.reads_R2_fpath is None else f'-2 {args.reads_R2_fpath}'

    unpaired_options = ' '.join(
        map(
            lambda s: f'-U {s}',
            args.unpaired_reads_fpaths
        )
    )

    command = ' '.join(
        [
            dependencies.bowtie2_fpath,
            f'-x {args.genome_index_base_fpath}',
            forw_reads_option,
            revr_reads_option,
            unpaired_options,
            f'--threads {args.n_threads}',
            f'-S {sam_outfpath}',
        ]
    )

    return command
# end def _configure_bow1tie2_build_command


def run_bowtie2(args, dependencies):

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
    command_str = _configure_bow1tie2_build_command(args, dependencies)

    print('Building index of the reference genome...')
    pipe = sp.Popen(command_str, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'Program `bowtie2-build` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if

    fasta_index_extention = '.fai'
    ref_genome_index_fpath = args.ref_genome_fpath + fasta_index_extention

    if not os.path.exists(ref_genome_index_fpath):
        print(f"""\nError: the index file `{ref_genome_index_fpath}`
    of the genome sequence file `{args.ref_genome_fpath}` does not exist after indexing""")
        print('This file must exist, though. Exitting...')
        sys.exit(1)
    # end if


    # Perform mapping
    command_str = _configure_bow1tie2_command(args, dependencies, sam_outfpath)

    print('Mapping the reads...')
    # print(command_str)
    pipe = sp.Popen(command_str, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'Program `bowtie2` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if

    if not os.path.exists(sam_outfpath):
        print(f'\nError: SAM file `{sam_outfpath}` does not exist after bowtie2 has mapped the reads')
        print('This file must exist, though. Exitting...')
        sys.exit(1)
    # end if

    return ReadMappingOutput(sam_outfpath)
# end def run_bowtie2