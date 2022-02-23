
import os
import sys

from src.shell import launch_command

from src.data_transfer_objects import ReadMapping, AlignmentIndex
from src.arguments import AlnPostprocessArguments
from src.dependencies import SamtoolsDependencies


def run_aln_postprocess(args, dependencies):

    bam_mapping = _convert_sam2bam(args, dependencies)
    bam_mapping.check_existance()
    args.raw_mapping.cleanup()

    sorted_bam_mapping = _sort_bam(bam_mapping, args, dependencies)
    sorted_bam_mapping.check_existance()
    bam_mapping.cleanup()

    bam_index = _index_bam(sorted_bam_mapping, args, dependencies)
    bam_index.check_existance()

    return sorted_bam_mapping
# end def


def _convert_sam2bam(postproc_args, dependencies):

    bam_outfpath = _configure_bam_outfpath(postproc_args)
    _print_converting_status_message(
        postproc_args.raw_mapping.alignment_fpath,
        bam_outfpath
    )

    command = _configure_sam2bam_command(
        postproc_args,
        dependencies,
        bam_outfpath
    )
    launch_command(command, 'samtools view')

    return ReadMapping(bam_outfpath)
# end def


def _configure_bam_outfpath(postproc_args):
    outdir_path = os.path.dirname(postproc_args.raw_mapping.alignment_fpath)
    bam_outfpath = os.path.join(
        outdir_path,
       '{}_{}.bam'.format(
            postproc_args.sample_name,
            postproc_args.output_suffix
        )
    )
    return bam_outfpath
# end def


def _print_converting_status_message(sam_fpath, bam_fpath):
    sam_basename = os.path.basename(sam_fpath)
    bam_basename = os.path.basename(bam_fpath)
    print(
        'Converting SAM to BAM (`{}` -> `{}`)...' \
            .format(sam_basename, bam_basename)
    )
# end def


def _configure_sam2bam_command(postproc_args,
                               dependencies,
                               bam_outfpath):

    input_sam_fpath = postproc_args.raw_mapping.alignment_fpath

    command = ' '.join(
        [
            dependencies.samtools_fpath, 'view',
            '-F 4 -b -O BAM',
            f'-@ {postproc_args.n_threads}',
            f'-T {postproc_args.reference_fpath}',
            f'-o {bam_outfpath}',
            input_sam_fpath,
        ]
    )

    return command
# end def


def _sort_bam(bam_mapping, postproc_args, dependencies):
    sorted_bam_outfpath = _configure_sorted_bam_outfpath(
        bam_mapping,
        postproc_args
    )
    print(
        'Sorting BAM (`{}`)...' \
            .format(os.path.basename(bam_mapping.alignment_fpath))
    )
    command = _configure_sort_bam_command(
        postproc_args,
        dependencies,
        bam_mapping,
        sorted_bam_outfpath
    )
    launch_command(command, 'samtools sort')

    return ReadMapping(sorted_bam_outfpath)
# end def


def _configure_sorted_bam_outfpath(bam_mapping, postproc_args):
    outdir_path = os.path.dirname(bam_mapping.alignment_fpath)
    sorted_bam_outfpath = os.path.join(
        outdir_path,
        '{}_{}.sorted.bam'.format(
            postproc_args.sample_name,
            postproc_args.output_suffix
        )
    )
    return sorted_bam_outfpath
# end def


def _configure_sort_bam_command(postproc_args,
                                dependencies,
                                bam_mapping,
                                sorted_bam_fpath):

    command = ' '.join(
        [
            dependencies.samtools_fpath, 'sort',
            '-O BAM',
            f'-@ {postproc_args.n_threads}',
            f'-o {sorted_bam_fpath}',
            bam_mapping.alignment_fpath,
        ]
    )

    return command
# end def


def _index_bam(sorted_bam_mapping, postproc_args, dependencies):
    print(
        'Indexing BAM (`{}`)...' \
            .format(os.path.basename(sorted_bam_mapping.alignment_fpath))
    )
    command = _configure_index_bam_command(
        postproc_args,
        dependencies,
        sorted_bam_mapping
    )
    launch_command(command, 'samtools index')

    return AlignmentIndex(sorted_bam_mapping.alignment_fpath)
# end def


def _configure_index_bam_command(postproc_args,
                                 dependencies,
                                 sorted_bam_mapping):
    command = ' '.join(
        [
            dependencies.samtools_fpath, 'index',
            '-b',
            f'-@ {postproc_args.n_threads}',
            sorted_bam_mapping.alignment_fpath,
        ]
    )

    return command
# end def
