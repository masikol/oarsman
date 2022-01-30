
import os

from src.shell import launch_command

from src.arguments import CallVariantsArguments
from src.dependencies import BcfVarCallDependencies
from src.data_transfer_objects import VariantCall


def run_bcftools_var_call(args, dependencies):

    print('Calling variants...')
    var_call = _call_variants(args, dependencies)
    var_call.check_existance()

    return var_call
# end def


def _call_variants(var_call_args, dependencies):
    var_call_fpath = _configure_var_call_fpath(var_call_args)
    command_str = _configure_variant_call_command(
        var_call_args,
        dependencies,
        var_call_fpath
    )
    launch_command(command_str, 'bcftools mpileup|call')
    return VariantCall(var_call_fpath)
# end def


def _configure_var_call_fpath(var_call_args):
    var_call_fpath = os.path.join(
        var_call_args.outdir_path,
        f'{var_call_args.sample_name}.bcf'
    )
    return var_call_fpath
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
