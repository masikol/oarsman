
import os

from src.shell import launch_command

from src.data_transfer_objects import SequenceFile


def make_consensus(args, dependencies):

    consensus_seq = _run_consensus_making(
        args,
        dependencies
    )
    consensus_seq.check_existance()

    return consensus_seq
# end def


def _run_consensus_making(cons_args, dependencies):

    consensus_outfpath = _configure_consensus_outfpath(cons_args)

    command_str = _configure_consensus_command(
        cons_args,
        dependencies,
        consensus_outfpath
    )

    launch_command(command_str, 'bcftools consensus')

    return SequenceFile(consensus_outfpath)
# end def


def _configure_consensus_outfpath(cons_args):
    consensus_outfpath = os.path.join(
        cons_args.outdir_path,
        f'{cons_args.sample_name}_consensus.fasta'
    )
    return consensus_outfpath
# end def


def _configure_consensus_command(cons_args,
                                 dependencies,
                                 consensus_outfpath):

    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'consensus',
            f'-f {cons_args.reference_fpath}',
            f'--prefix {cons_args.sample_name}_',
            f'-o {consensus_outfpath}',
            cons_args.var_call.var_call_fpath
        ]
    )

    return command
# end def
