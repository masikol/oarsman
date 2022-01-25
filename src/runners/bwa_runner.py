
import os

import src.filesystem as fs
from src.input_modes import InputModes
from src.fatal_errors import FatalError
from src.runners.shell import launch_command

from src.data_transfer_objects import BwaSeqIndex, ReadMapping
from src.arguments import ReadMappingArguments
from src.dependencies import BwaDependencies


def run_bwa(args, dependencies):

    reference_index = _create_reference_index(args, dependencies)

    if args.input_mode == InputModes.IlluminaPE:
        mapping = _map_paired_reads(args, dependencies, reference_index)
    else:
        mapping = _map_unpaired_reads(args, dependencies, reference_index)
    # end if

    mapping.check_existance()

    return mapping
# end def


def _create_reference_index(mapping_args, dependencies):

    index_base_fpath = _make_index_base(mapping_args)

    print('Building index of the reference genome...')
    command_str = _configure_bwa_index_command(
        mapping_args,
        dependencies, 
        index_base_fpath
    )
    launch_command(command_str, 'bwa index')

    indexing_output = BwaSeqIndex(index_base_fpath)
    indexing_output.check_existance()

    return indexing_output
# end def


def _configure_bwa_index_command(mapping_args,
                                 dependencies,
                                 ref_index_base_fpath):

    command = ' '.join(
        [
            dependencies.bwa_fpath, 'index',
            f'-p {ref_index_base_fpath}',
            mapping_args.reference_fpath
        ]
    )

    return command
# end def


def _make_index_base(mapping_args):
    reference_fpath_basename = os.path.basename(mapping_args.reference_fpath)
    index_base_fpath = os.path.join(
        mapping_args.index_dir_path,
        fs.rm_fasta_extention(reference_fpath_basename) + '_index'
    )
    return index_base_fpath
# end def


def _configure_sam_outfpath(mapping_args):
    sam_outfpath = os.path.join(
        mapping_args.outdir_path,
        '{}_{}.sam'.format(mapping_args.sample_name, mapping_args.output_suffix)
    )
    return sam_outfpath
# end def


def _map_unpaired_reads(mapping_args, dependencies, reference_index):

    reads_fpath = mapping_args.reads.reads_fpath
    sam_outfpath = _configure_sam_outfpath(mapping_args)

    command_str = _configure_bwa_unpaired_command(
        reads_fpath,
        reference_index,
        mapping_args,
        dependencies,
        sam_outfpath
    )

    print('Mapping the reads...')
    launch_command(command_str, 'bwa mem')

    return ReadMapping(sam_outfpath)
# end def


def _map_paired_reads(mapping_args, dependencies, reference_index):

    sam_outfpath = _configure_sam_outfpath(mapping_args)

    frw_reads_fpath = mapping_args.reads.reads_frw_fpath
    rvr_reads_fpath = mapping_args.reads.reads_rvr_fpath

    print('Mapping the paired reads...')
    command_str = _configure_bwa_paired_command(
        frw_reads_fpath,
        rvr_reads_fpath,
        reference_index,
        mapping_args,
        dependencies,
        sam_outfpath
    )
    launch_command(command_str, 'bwa mem')
    print('done.')


    # TODO: map unpaired reads as well. Problem: unable to process merged SAM files further (SAM headers)
    # unpaired_reads_fpaths = (
    #     mapping_args.reads.reads_frw_upr_fpaths,
    #     mapping_args.reads.reads_rvr_upr_fpaths,
    # )

    # print('Mapping the unpaired reads...')
    # for reads_fpath in unpaired_reads_fpaths:
    #     command_str = _configure_bwa_unpaired_command(
    #         reads_fpath,
    #         reference_index,
    #         mapping_args,
    #         dependencies,
    #         sam_outfpath,
    #         append=True
    #     )
    #     launch_command(command_str)
    # # end def
    # print('done.')

    return ReadMapping(sam_outfpath)
# end def


def _configure_bwa_unpaired_command(reads_fpath,
                                    reference_index,
                                    args,
                                    dependencies,
                                    sam_outfpath,
                                    append=False):

    index_fpath = reference_index.index_base_fpath

    if append:
        output_cmd_part = '>> {}'.format(sam_outfpath)
    else:
        output_cmd_part = '-o {}'.format(sam_outfpath)
    # end if

    command = ' '.join(
        [
            dependencies.bwa_fpath, 'mem',
            f'-t {args.n_threads}',
            f'{index_fpath}',
            reads_fpath,
            output_cmd_part,
        ]
    )

    return command
# end def


def _configure_bwa_paired_command(frw_reads_fpath,
                                  rvr_reads_fpath,
                                  reference_index,
                                  args,
                                  dependencies,
                                  sam_outfpath):

    paired_reads_cmd_part = ' '.join(
        (frw_reads_fpath, rvr_reads_fpath)
    )

    index_fpath = reference_index.index_base_fpath

    command = ' '.join(
        [
            dependencies.bwa_fpath, 'mem',
            f'-t {args.n_threads}',
            f'-o {sam_outfpath}',
            f'{index_fpath}',
            paired_reads_cmd_part,
        ]
    )

    return command
# end def
