
import os

from src.shell import launch_command

from src.data_transfer_objects import VariantCall, \
                                      VariantCallIndex


def run_var_vall_postprocess(args, dependencies):

    raw_var_call_index = \
        _index_var_call(args.raw_var_call, args, dependencies)
    raw_var_call_index.check_existance()

    print('Normalizing variants...')
    normalized_var_call = _normalize_indels(
        args.raw_var_call,
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

    return filtered_var_call
# end def


def _index_var_call(var_call, postproc_args, dependencies):
    command_str = _configure_index_bcf_command(
        var_call,
        postproc_args,
        dependencies
    )
    launch_command(command_str, 'bcftools index')

    return VariantCallIndex(var_call.var_call_fpath)
# end def


def _configure_index_bcf_command(var_call, postproc_args, dependencies):

    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'index',
            f'--threads {postproc_args.n_threads}',
            var_call.var_call_fpath
        ]
    )

    return command
# end def


def _normalize_indels(raw_var_call, postproc_args, dependencies):
    normalized_var_call_fpath = \
        _configure_normalized_var_call_fpath(postproc_args)

    command_str = _configure_norm_indels_command(
        postproc_args,
        dependencies,
        raw_var_call,
        normalized_var_call_fpath
    )
    launch_command(command_str, 'bcftools norm')

    return VariantCall(normalized_var_call_fpath)
# end def


def _configure_normalized_var_call_fpath(postproc_args):
    normalized_var_call_fpath = os.path.join(
        postproc_args.outdir_path,
        f'{postproc_args.sample_name}.norm.bcf'
    )
    return normalized_var_call_fpath
# end def


def _configure_norm_indels_command(postproc_args, dependencies, var_call, outfpath):

    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'norm',
            '-Ob',
            f'--threads {postproc_args.n_threads}',
            f'-f {postproc_args.reference_fpath}',
            f'-o {outfpath}',
            var_call.var_call_fpath
        ]
    )

    return command
# end def


def _filter_variants(normalized_var_call, postproc_args, dependencies):
    normalized_var_call_fpath = \
        _configure_filtered_var_call_fpath(postproc_args)

    command_str = _configure_filter_command(
        postproc_args,
        dependencies,
        normalized_var_call,
        normalized_var_call_fpath
    )
    launch_command(command_str, 'bcftools filter')

    return VariantCall(normalized_var_call_fpath)
# end def


def _configure_filtered_var_call_fpath(postproc_args):
    normalized_var_call_fpath = os.path.join(
        postproc_args.outdir_path,
        f'{postproc_args.sample_name}.filt.bcf'
    )
    return normalized_var_call_fpath
# end def


def _configure_filter_command(postproc_args,
                              dependencies,
                              normalized_var_call,
                              outfpath):

    # indel_gap = 10
    # spn_gap = 1

    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'filter',
            '-Ob',
            f'--threads {postproc_args.n_threads}',
            f"-e '%QUAL<{postproc_args.min_variant_qual}'",
            # f'--IndelGap {indel_gap}',
            # f'--SnpGap {spn_gap}',
            f'-o {outfpath}',
            normalized_var_call.var_call_fpath
        ]
    )

    return command
# end def
