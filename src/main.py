
import sys

from src.fatal_errors import FatalError
from src.printing import print_err

from src.reads import Reads
from src.parse_arguments import parse_arguments
from src.seq_sample import Sample

from src.runners.kromsatel_runner import run_kromsatel
from src.arguments import KromsatelArguments
from src.dependencies import KromsatelDependencies
from src.data_transfer_objects import KromsatelResult

from src.arguments import ReadMappingArguments

from src.runners.bwa_runner import run_bwa
from src.dependencies import BwaDependencies

from src.runners.aln_postprocess_runner import run_aln_postprocess
from src.arguments import AlnPostprocessArguments
from src.dependencies import SamtoolsDependencies

from src.runners.bcftools_var_call_runner import run_bcftools_var_call
from src.arguments import CallVariantsArguments
from src.dependencies import BcfVarCallDependencies
from src.data_transfer_objects import SequenceFile

from src.runners.highlighter_runner import run_highlighter
from src.data_transfer_objects import SequenceFile
from src.arguments import ConsensAnnotArguments
from src.dependencies import HighlighterDependencies



def main():

    # Parse command line arguments
    try:
        oarsman_args, oarsman_dependencies = parse_arguments()
    except FatalError as err:
        print_err(str(err))
        sys.exit(1)
    # end try

    reference_fpath = oarsman_args.reference_fpath

    # From this moment we will iterate over samples
    n_samples = oarsman_args.get_number_of_samples()
    for i_sample in range(n_samples):

        sample = Sample(oarsman_args.reads_R1_fpaths[i_sample][0], i_sample)

        print(
            '  \n|=== Processing sample #{}/{}: {} ===|' \
                .format(sample.ordinal_num+1, n_samples, sample.name)
        )

        cleaned_reads: Reads = _clean_reads(
            oarsman_args,
            oarsman_dependencies,
            sample
        )

        mapping = _map_reads(
            oarsman_args,
            oarsman_dependencies,
            cleaned_reads,
            reference_fpath,
            sample,
            'reference'
        )

        # Call variants
        print('Making consensus...')
        var_call_args = oarsman_args.get_var_call_args(
            sample.name,
            mapping.alignment_fpath
        )
        var_call_dependencies = \
            oarsman_dependencies.get_bcftools_dependencies()

        consensus_seq = run_bcftools_var_call(
            var_call_args,
            var_call_dependencies
        )

        print(f'The consensus is created: `{consensus_seq.file_path}`')

        print('Mapping the reads to the consensus')
        # Map the reads to the consensus
        mapping_to_consensus = _map_reads(
            oarsman_args,
            oarsman_dependencies,
            cleaned_reads,
            consensus_seq.file_path,
            sample,
            'consensus'
        )

        consensus_annotation: SequenceFile = _annotate_consensus(
            oarsman_args,
            oarsman_dependencies,
            sample,
            consensus_seq.file_path,
            mapping_to_consensus
        )

        print('Annotated consensus: `{}`'.format(consensus_annotation.file_path))
    # end for
# end def


def _clean_reads(oarsman_args, oarsman_dependencies, sample):

    # Run kromsatel
    kromsatel_arguments: KromsatelArguments = \
        oarsman_args.get_kromsatel_args(sample.ordinal_num)
    kromsatel_dependencies: KromsatelDependencies = \
        oarsman_dependencies.get_kromsatel_dependencies()

    kromsatel_output: KromsatelResult = run_kromsatel(
        kromsatel_arguments,
        kromsatel_dependencies
    )

    return kromsatel_output
# end def


def _map_reads(
    oarsman_args,
    oarsman_dependencies,
    reads,
    reference_fpath,
    sample,
    output_suffix
):

    read_map_args = oarsman_args.get_read_mapping_args(
        sample.name,
        reads,
        reference_fpath,
        output_suffix
    )

    read_map_dependencies = oarsman_dependencies.get_bwa_dependencies()

    read_mapping = run_bwa(
        read_map_args,
        read_map_dependencies
    )

    # Postprocess alignment data
    aln_postproc_args = oarsman_args.get_aln_postprocess_args(
        sample.name,
        read_mapping,
        reference_fpath,
        output_suffix
    )
    aln_postproc_dependencies = oarsman_dependencies.get_samtools_dependencies()

    aln_postprocess_output = run_aln_postprocess(
        aln_postproc_args,
        aln_postproc_dependencies
    )

    return aln_postprocess_output
# end def


def _annotate_consensus(oarsman_args, oarsman_dependencies, sample, seq_fpath, mapping):

    cons_annot_args = oarsman_args.get_consens_annot_args(
        sample.name,
        seq_fpath,
        mapping
    )

    cons_annot_dependencies = oarsman_dependencies.get_highlighter_dependencies()

    annotated_seq: SequenceFile = run_highlighter(
        cons_annot_args,
        cons_annot_dependencies
    )

    return annotated_seq
# end def
