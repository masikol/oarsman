
import os
import sys
import subprocess as sp

from src.mapping import Mapping
from src.oarsman_arguments import AlnPreprocessArguments
from src.oarsman_dependencies import SamtoolsDependencies

from src.filesystem import rm_file


def _configure_sam2bam_command(
    args: AlnPreprocessArguments,
    dependencies: SamtoolsDependencies,
    sam_fpath,
    bam_fpath
):

    command = ' '.join(
        [
            dependencies.samtools_fpath, 'view',
            '-F 4 -b -O BAM',
            f'-@ {args.n_threads}',
            f'-T {args.ref_genome_fpath}',
            f'-o {bam_fpath}',
            sam_fpath,
        ]
    )

    return command
# def _configure_sam2bam_command


def _configure_sort_bam_command(
    args: AlnPreprocessArguments,
    dependencies: SamtoolsDependencies,
    bam_fpath,
    sorted_bam_fpath
):

    command = ' '.join(
        [
            dependencies.samtools_fpath, 'sort',
            '-O BAM',
            f'-@ {args.n_threads}',
            f'-o {sorted_bam_fpath}',
            bam_fpath,
        ]
    )

    return command
# def _configure_sort_bam_command


def _configure_index_bam_command(
    args: AlnPreprocessArguments,
    dependencies: SamtoolsDependencies,
    sorted_bam_fpath
):

    command = ' '.join(
        [
            dependencies.samtools_fpath, 'index',
            '-b',
            f'-@ {args.n_threads}',
            sorted_bam_fpath,
        ]
    )

    return command
# def _configure_index_bam_command


def run_aln_preprocess(args: AlnPreprocessArguments, dependencies: SamtoolsDependencies):


    outdir_path = os.path.dirname(args.raw_alignment_fpath)

    bam_outfpath = os.path.join(
        outdir_path,
       '{}_{}.bam'.format(args.sample_name, args.output_suffix)
    )
    sorted_bam_outfpath = os.path.join(
        outdir_path,
        '{}_{}.sorted.bam'.format(args.sample_name, args.output_suffix)
    )


    # Transform SAM to BAM
    command_str = _configure_sam2bam_command(
        args,
        dependencies,
        args.raw_alignment_fpath,
        bam_outfpath
    )

    print('Converting SAM to BAM (`{}` to `{}`)...'.format(
        args.raw_alignment_fpath, bam_outfpath
    ))
    pipe = sp.Popen(command_str, shell=True, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'Program `samtools view` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if
    rm_file(args.raw_alignment_fpath)

    if not os.path.exists(bam_outfpath):
        print(f'\nError: BAM file `{bam_outfpath}` does not exist after SAM->BAM conversion')
        print('This file must exist, though. Exitting...')
        sys.exit(1)
    # end if


    # Sort BAM file
    command_str = _configure_sort_bam_command(args, dependencies, bam_outfpath, sorted_bam_outfpath)

    print('Sorting BAM (`{}` -> `{}`)...'.format(
        bam_outfpath, sorted_bam_outfpath
    ))
    pipe = sp.Popen(command_str, shell=True, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'Program `samtools sort` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if
    rm_file(bam_outfpath)

    if not os.path.exists(sorted_bam_outfpath):
        print(f'\nError: sorted BAM file `{sorted_bam_outfpath}` does not exist after BAM sorting')
        print('This file must exist, though. Exitting...')
        sys.exit(1)
    # end if


    # Index BAM file
    command_str = _configure_index_bam_command(args, dependencies, sorted_bam_outfpath)

    print('Indexing sorted BAM (`{}`)...'.format(
        sorted_bam_outfpath
    ))
    pipe = sp.Popen(command_str, shell=True, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'Program `samtools index` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if

    bam_index_extention = '.bai'
    sorted_bam_index_fpath = sorted_bam_outfpath + bam_index_extention

    if not os.path.exists(sorted_bam_outfpath):
        print(f"""\nError: index file `{sorted_bam_index_fpath}`
    of the sorted BAM file `{sorted_bam_outfpath}` does not exist after indexing""")
        print('This file must exist, though. Exitting...')
        sys.exit(1)
    # end if

    return Mapping(sorted_bam_outfpath)
# end def run_aln_preprocess
