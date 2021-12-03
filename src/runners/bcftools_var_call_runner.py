
import os
import sys
import subprocess as sp

from src.oarsman_arguments import CallVariantsArguments
from src.oarsman_dependencies import BcfVarCallDependencies
from src.output_data import VariantCallingOutput


def _configure_variant_call_command(
    args: CallVariantsArguments,
    dependencies: BcfVarCallDependencies,
    variants_fpath: str
):

    max_coverage = 50000 # reads (for now)
    ploidy = 1 # haploid

    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'mpileup',
            '-Ob',
            f'--max-depth {max_coverage}', f'--max-idepth {max_coverage}',
            f'-f {args.reference_fpath}',
            f'--threads {args.n_threads}',
            args.alignment_fpath,
            '|',
            dependencies.bcftools_fpath, 'call',
            '-mv', f'--ploidy {ploidy}', '-Ou',
            f'--threads {args.n_threads}',
            f'-o {variants_fpath}'
        ]
    )

    return command
# end def _configure_bow1tie2_build_command


def _configure_norm_indels_command(
    args: CallVariantsArguments,
    dependencies: BcfVarCallDependencies,
    variants_fpath: str,
    norm_variants_fpath: str
):

    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'norm',
            '-Ob',
            f'--threads {args.n_threads}',
            f'-f {args.reference_fpath}',
            f'-o {norm_variants_fpath}',
            variants_fpath
        ]
    )

    return command
# end def _configure_norm_indels_command


def _configure_filter_command(
    args: CallVariantsArguments,
    dependencies: BcfVarCallDependencies,
    norm_variants_fpath: str,
    filt_variants_fpath: str
):

    indel_gap = 10
    spn_gap = 1

    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'filter',
            '-Ob',
            f'--threads {args.n_threads}',
            f"-e '%QUAL<{args.min_variant_qual}'",
            f'--IndelGap {indel_gap}',
            f'--SnpGap {spn_gap}',
            f'-o {filt_variants_fpath}',
            norm_variants_fpath
        ]
    )

    return command
# end def _configure_filter_command


def _configure_consensus_command(
    args: CallVariantsArguments,
    dependencies: BcfVarCallDependencies,
    filt_variants_fpath: str,
    consensus_fpath: str
):

    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'consensus',
            f'-f {args.reference_fpath}',
            f'--prefix {args.sample_name}_',
            f'-o {consensus_fpath}',
            filt_variants_fpath
        ]
    )

    return command
# end def _configure_consensus_command


def _configure_index_bcf_command(
    args: CallVariantsArguments,
    dependencies: BcfVarCallDependencies,
    file_path_to_index: str,
):

    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'index',
            f'--threads {args.n_threads}',
            file_path_to_index
        ]
    )

    return command
# end def _configure_bow1tie2_build_command


def run_bcftools_var_call(args, dependencies):

    for outdir_path in (args.var_calls_dirpath, args.consensus_dirpath):
        if not os.path.isdir(outdir_path):
            try:
                os.makedirs(outdir_path)
            except OSError as err:
                print(f'\nError: cannot create directory `{outdir_path}`')
                print(str(err))
                sys.exit(1)
            # end try
        # end if
    # end for

    baseline_variants_fpath = os.path.join(
        args.var_calls_dirpath,
        f'{args.sample_name}.bcf'
    )
    norm_variants_fpath = os.path.join(
        args.var_calls_dirpath,
        f'{args.sample_name}.norm.bcf'
    )
    filt_variants_fpath = os.path.join(
        args.var_calls_dirpath,
        f'{args.sample_name}.norm.filt.bcf'
    )
    bcf_index_extention = '.csi'

    consensus_fpath = os.path.join(
        args.consensus_dirpath,
        f'{args.sample_name}_consensus.fasta'
    )

    # Call variants
    command_str = _configure_variant_call_command(
        args,
        dependencies,
        baseline_variants_fpath
    )

    print('Calling variants')
    pipe = sp.Popen(command_str, shell=True, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'Program `bcftools` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if

    if not os.path.exists(baseline_variants_fpath):
        print(f"""\nError: the variant file `{baseline_variants_fpath}`
    does not exist after variant calling""")
        print('This file must exist, though. Exitting...')
        sys.exit(1)
    # end if


    # Index baseline variants file
    command_str = _configure_index_bcf_command(
        args,
        dependencies,
        baseline_variants_fpath
    )

    print(f'Indexing file `{baseline_variants_fpath}`')
    pipe = sp.Popen(command_str, shell=True, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'Program `bcftools` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if

    baseline_var_call_index = baseline_variants_fpath + bcf_index_extention

    if not os.path.exists(baseline_var_call_index):
        print(f"""\nError: the index file `{baseline_var_call_index}`
    of the variant file `{baseline_variants_fpath}` does not exist after indexing""")
        print('This file must exist, though. Exitting...')
        sys.exit(1)
    # end if


    # Normalize indels
    command_str = _configure_norm_indels_command(
        args,
        dependencies,
        baseline_variants_fpath,
        norm_variants_fpath
    )

    print('Normalizing indels')
    pipe = sp.Popen(command_str, shell=True, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'Program `bcftools` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if

    if not os.path.exists(norm_variants_fpath):
        print(f"""\nError: the variant file `{baseline_var_call_index}`
    does not exist after indel normalization""")
        print('This file must exist, though. Exitting...')
        sys.exit(1)
    # end if


    # Index normalized variants file
    command_str = _configure_index_bcf_command(
        args,
        dependencies,
        norm_variants_fpath
    )

    print(f'Indexing file `{norm_variants_fpath}`')
    pipe = sp.Popen(command_str, shell=True, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'Program `bcftools` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if

    norm_var_call_index = norm_variants_fpath + bcf_index_extention

    if not os.path.exists(norm_var_call_index):
        print(f"""\nError: the index file `{norm_var_call_index}`
    of the variant file `{norm_variants_fpath}` does not exist after indexing""")
        print('This file must exist, though. Exitting...')
        sys.exit(1)
    # end if


    # Filter variants
    command_str = _configure_filter_command(
        args,
        dependencies,
        norm_variants_fpath,
        filt_variants_fpath
    )

    print('Filtering variants')
    pipe = sp.Popen(command_str, shell=True, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'Program `bcftools` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        print(f'Command: "{command_str}"')
        sys.exit(1)
    # end if

    if not os.path.exists(filt_variants_fpath):
        print(f"""\nError: the variant file `{filt_variants_fpath}`
    does not exist after filtering.""")
        print('This file must exist, though. Exitting...')
        sys.exit(1)
    # end if


    # Index filtered variants file
    command_str = _configure_index_bcf_command(
        args,
        dependencies,
        filt_variants_fpath
    )

    print(f'Indexing file `{filt_variants_fpath}`')
    pipe = sp.Popen(command_str, shell=True, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'Program `bcftools` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if

    filt_var_call_index = filt_variants_fpath + bcf_index_extention

    if not os.path.exists(filt_var_call_index):
        print(f"""\nError: the index file `{filt_var_call_index}`
    of the variant file `{filt_variants_fpath}` does not exist after indexing""")
        print('This file must exist, though. Exitting...')
        sys.exit(1)
    # end if


    # Make a consensus
    command_str = _configure_consensus_command(
        args,
        dependencies,
        filt_variants_fpath,
        consensus_fpath
    )

    print(f'Making consensus...')
    pipe = sp.Popen(command_str, shell=True, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'Program `bcftools` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if

    if not os.path.exists(consensus_fpath):
        print(f"""\nError: the consensus file `{consensus_fpath}`
    does not exist after consensus generation""")
        print('This file must exist, though. Exitting...')
        sys.exit(1)
    # end if


    return VariantCallingOutput(
        consensus_fpath
    )
# end def run_bcftools_var_call