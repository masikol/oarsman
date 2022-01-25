
import os

import src.filesystem as fs
from src.input_modes import InputModes
from src.fatal_errors import FatalError
from src.runners.shell import launch_command

from src.arguments import KromsatelArguments
from src.dependencies import KromsatelDependencies
from src.data_transfer_objects import UnpairedKromsatelResult, \
                                      PairedKromsatelResult



def run_kromsatel(args, dependencies):

    sample_chunks = _get_sample_chunks(args)

    per_chunk_outputs = list()

    for sample_chunk in sample_chunks:
        curr_output = _launch_kromsatel_process(sample_chunk, args, dependencies)
        curr_output.check_existance()
        per_chunk_outputs.append(curr_output)
    # end for

    if len(per_chunk_outputs) > 1:
        output = _merge_outputs(per_chunk_outputs)
    else:
        output = next(iter(per_chunk_outputs))
    # end if

    return output
# end def run_kromnsatel


def _get_sample_chunks(kromsatel_args):

    input_mode = kromsatel_args.input_mode

    if input_mode == InputModes.IlluminaSE:
        sample_chunks = kromsatel_args.reads_R1_fpaths
    elif input_mode == InputModes.IlluminaPE:
        sample_chunks = zip(
            kromsatel_args.reads_R1_fpaths,
            kromsatel_args.reads_R2_fpaths
        )
    elif input_mode == InputModes.Nanopore:
        sample_chunks = kromsatel_args.reads_long_fpaths
    else:
        raise FatalError(
            '\nError: Internal error.\n' \
            'Please, contact the developer.\n' \
            'Error description: invalid mode in "_get_sampe_chunks".'
        )
    # end if

    return sample_chunks
# end def


def _launch_kromsatel_process(sample_chunk, kromsatel_args, dependencies):
    command_str = _configure_kromsatel_command(
        sample_chunk,
        kromsatel_args,
        dependencies
    )

    launch_command(command_str, 'kromsatel', print_stdout=True)

    return _make_output(sample_chunk, kromsatel_args)
# end def


def _configure_kromsatel_command(sample_chunk, kromsatel_args, dependencies):

    if os.path.isfile(dependencies.kromsatel_fpath):
        kromsatel_script_str = 'python3 {}'.format(dependencies.kromsatel_fpath)
    else:
        kromsatel_script_str = dependencies.kromsatel_fpath
    # end if

    reads_command_part = _configure_reads_command_part(sample_chunk, kromsatel_args)

    command = ' '.join(
        [
            kromsatel_script_str,
            reads_command_part,
            f'-p {kromsatel_args.primers_fpath}',
            f'-r {kromsatel_args.reference_fpath}',
            f'-t {kromsatel_args.n_threads}',
            f'-o {kromsatel_args.outdir_path}',
            kromsatel_args.advanced_args
        ]
    )

    return command
# end def


def _configure_reads_command_part(sample_chunk, kromsatel_args):

    reads_command_part = ''
    input_mode = kromsatel_args.input_mode

    if input_mode == InputModes.IlluminaSE:
        reads_command_part = '-1 {}'.format(sample_chunk)
    elif input_mode == InputModes.IlluminaPE:
        reads_command_part = '-1 {} -2 {}'.format(sample_chunk[0], sample_chunk[1])
    elif input_mode == InputModes.Nanopore:
        reads_command_part = '-l {}'.format(sample_chunk)
    else:
        raise FatalError(
            '\nError: Internal error.\n' \
            'Please, contact the developer.\n' \
            'Error description: invalid mode in "_get_sampe_chunks".'
        )
    # end if
    return reads_command_part
# end def



def _make_output(sample_chunk, kromsatel_args):

    input_mode = kromsatel_args.input_mode

    unpaired_input_modes = (InputModes.IlluminaSE, InputModes.Nanopore)

    if input_mode in unpaired_input_modes:
        outfpath = _configure_outfpath(
            sample_chunk,
            'cleaned',
            kromsatel_args.outdir_path
        )
        output = UnpairedKromsatelResult(outfpath)

    elif input_mode == InputModes.IlluminaPE:
        frw_fpath, rvr_fpath = sample_chunk[0], sample_chunk[1]
        frw_outfpath, rvr_outfpath = _configure_outfpath_pair(
            frw_fpath,
            rvr_fpath,
            'cleaned',
            kromsatel_args.outdir_path
        )
        frw_upr_outfpath, rvr_upr_outfpath = _configure_outfpath_pair(
            frw_fpath,
            rvr_fpath,
            'cleaned_unpaired',
            kromsatel_args.outdir_path
        )
        output = PairedKromsatelResult(frw_outfpath, rvr_outfpath, frw_upr_outfpath, rvr_upr_outfpath)
    # end if

    return output
# end def


def _configure_outfpath_pair(frw_fpath, rvr_fpath, suffix, outdir_path):
    frw_outfpath = _configure_outfpath(
        frw_fpath,
        suffix,
        outdir_path
    )
    rvr_outfpath = _configure_outfpath(
        rvr_fpath,
        suffix,
        outdir_path
    )
    return frw_outfpath, rvr_outfpath
# end def


def _configure_outfpath(input_fpath, suffix, outdir_path):
    basename_no_ext = fs.rm_fastq_extention(os.path.basename(input_fpath))
    return os.path.join(
        outdir_path,
        '{}_{}.fastq.gz'.format(basename_no_ext, suffix)
    )
# end def


def _merge_outputs(per_chunk_outputs):

    final_output = per_chunk_outputs[0]

    for partial_output in per_chunk_outputs[1:]:
        partial_output.append_to(final_output)
        partial_output.cleanup()
    # end for

    return final_output
# end def
