
import os

from src.runners.shell import launch_command

from src.arguments import CallVariantsArguments
from src.dependencies import BcfVarCallDependencies
from src.data_transfer_objects import VariantCall, \
                                      VariantCallIndex, \
                                      SequenceFile


def run_bcftools_var_call(args, dependencies):

    print('Calling variants...')
    baseline_var_call = _call_variants(args, dependencies)
    baseline_var_call.check_existance()

    baseline_var_call_index = \
        _index_var_call(baseline_var_call, args, dependencies)
    baseline_var_call_index.check_existance()

    print('Normalizing variants...')
    normalized_var_call = _normalize_indels(
        baseline_var_call,
        args,
        dependencies
    )
    normalized_var_call.check_existance()

    normalized_var_call_index = \
        _index_var_call(normalized_var_call, args, dependencies)
    normalized_var_call_index.check_existance()

    print('Filtering variants...')
    filtered_var_call = _filter_variants(normalized_var_call, args, dependencies)
    filtered_var_call.check_existance()

    filtered_var_call_index = \
        _index_var_call(filtered_var_call, args, dependencies)
    filtered_var_call_index.check_existance()

    print('Finally making consensus...')
    consensus_seq = _make_consensus_seq(
        filtered_var_call,
        args,
        dependencies
    )
    consensus_seq.check_existance()

    return consensus_seq
# end def


def _call_variants(var_call_args, dependencies):
    baseline_variants_fpath = _configure_baseline_var_call_fpath(var_call_args)
    command_str = _configure_variant_call_command(
        var_call_args,
        dependencies,
        baseline_variants_fpath
    )
    launch_command(command_str, 'bcftools mpileup|call')
    return VariantCall(baseline_variants_fpath)
# end def


def _configure_baseline_var_call_fpath(var_call_args):
    baseline_variants_fpath = os.path.join(
        var_call_args.outdir_path,
        f'{var_call_args.sample_name}.bcf'
    )
    return baseline_variants_fpath
# end def


def _configure_variant_call_command(var_call_args, dependencies, variants_fpath):

    max_coverage = 50000 # reads (for now)
    ploidy = 1 # haploid

    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'mpileup',
            '-Ob',
            f'--max-depth {max_coverage}', f'--max-idepth {max_coverage}',
            f'-f {var_call_args.reference_fpath}',
            f'--threads {var_call_args.n_threads}',
            var_call_args.alignment_fpath,
            '|',
            dependencies.bcftools_fpath, 'call',
            '-mv', f'--ploidy {ploidy}', '-Ou',
            f'--threads {var_call_args.n_threads}',
            f'-o {variants_fpath}'
        ]
    )

    return command
# end def


def _index_var_call(baseline_var_call, var_call_args, dependencies):
    command_str = _configure_index_bcf_command(
        baseline_var_call,
        var_call_args,
        dependencies
    )
    launch_command(command_str, 'bcftools index')

    return VariantCallIndex(baseline_var_call.var_call_fpath)
# end def


def _configure_index_bcf_command(var_call, var_call_args, dependencies):

    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'index',
            f'--threads {var_call_args.n_threads}',
            var_call.var_call_fpath
        ]
    )

    return command
# end def


def _normalize_indels(baseline_var_call, var_call_args, dependencies):
    normalized_var_call_fpath = \
        _configure_normalized_var_call_fpath(var_call_args)

    command_str = _configure_norm_indels_command(
        var_call_args,
        dependencies,
        baseline_var_call,
        normalized_var_call_fpath
    )
    launch_command(command_str, 'bcftools norm')

    return VariantCall(normalized_var_call_fpath)
# end def


def _configure_normalized_var_call_fpath(var_call_args):
    normalized_var_call_fpath = os.path.join(
        var_call_args.outdir_path,
        f'{var_call_args.sample_name}.norm.bcf'
    )
    return normalized_var_call_fpath
# end def


def _configure_norm_indels_command(var_call_args, dependencies, var_call, outfpath):

    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'norm',
            '-Ob',
            f'--threads {var_call_args.n_threads}',
            f'-f {var_call_args.reference_fpath}',
            f'-o {outfpath}',
            var_call.var_call_fpath
        ]
    )

    return command
# end def


def _filter_variants(normalized_var_call, var_call_args, dependencies):
    normalized_var_call_fpath = \
        _configure_filtered_var_call_fpath(var_call_args)

    command_str = _configure_filter_command(
        var_call_args,
        dependencies,
        normalized_var_call,
        normalized_var_call_fpath
    )
    launch_command(command_str, 'bcftools filter')

    return VariantCall(normalized_var_call_fpath)
# end def


def _configure_filtered_var_call_fpath(var_call_args):
    normalized_var_call_fpath = os.path.join(
        var_call_args.outdir_path,
        f'{var_call_args.sample_name}.filt.bcf'
    )
    return normalized_var_call_fpath
# end def


def _configure_filter_command(var_call_args,
                              dependencies,
                              normalized_var_call,
                              outfpath):

    indel_gap = 10
    spn_gap = 1

    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'filter',
            '-Ob',
            f'--threads {var_call_args.n_threads}',
            f"-e '%QUAL<{var_call_args.min_variant_qual}'",
            f'--IndelGap {indel_gap}',
            f'--SnpGap {spn_gap}',
            f'-o {outfpath}',
            normalized_var_call.var_call_fpath
        ]
    )

    return command
# end def


def _make_consensus_seq(filtered_var_call, var_call_args, dependencies):

    consensus_outfpath = _configure_consensus_outfpath(var_call_args)

    command_str = _configure_consensus_command(
        var_call_args,
        dependencies,
        filtered_var_call,
        consensus_outfpath
    )

    launch_command(command_str, 'bcftools consensus')

    return SequenceFile(consensus_outfpath)
# end def


def _configure_consensus_outfpath(var_call_args):
    consensus_outfpath = os.path.join(
        var_call_args.consensus_dirpath,
        f'{var_call_args.sample_name}_consensus.fasta'
    )
    return consensus_outfpath
# end def


def _configure_consensus_command(var_call_args,
                                 dependencies,
                                 filtered_var_call,
                                 consensus_outfpath):

    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'consensus',
            f'-f {var_call_args.reference_fpath}',
            f'--prefix {var_call_args.sample_name}_',
            f'-o {consensus_outfpath}',
            filtered_var_call.var_call_fpath
        ]
    )

    return command
# end def
