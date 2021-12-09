
import sys

from src.mapping import Mapping
from src.reads import Reads
from src.parse_arguments import parse_arguments
from src.seq_sample import Sample

from src.runners.make_amplicons_runner import run_make_amplicons
from src.oarsman_arguments import MakeAmpliconsArguments
from src.oarsman_dependencies import MakeAmpliconsDependencies
from src.output_data import MakeAmpliconsOutput

from src.runners.make_db_runner import run_make_db
from src.oarsman_arguments import MakeDbArguments
from src.oarsman_dependencies import MakeDbDependencies
from src.output_data import MakeDbOutput

from src.runners.kromsatel_runner import run_kromsatel
from src.oarsman_arguments import KromsatelArguments
from src.oarsman_dependencies import KromsatelDependencies
from src.output_data import KromsatelOutput

from src.runners.pair_runner import run_pair
from src.oarsman_arguments import PairArguments
from src.oarsman_dependencies import PairDependencies
# from src.output_data import PairOutput


from src.oarsman_arguments import ReadMappingArguments
# from src.output_data import ReadMappingOutput

# from src.runners.bowtie2_runner import run_bowtie2
# from src.oarsman_dependencies import Bowtie2Dependencies

from src.runners.bwa_runner import run_bwa
from src.oarsman_dependencies import BwaDependencies


from src.runners.aln_preprocess_runner import run_aln_preprocess
from src.oarsman_arguments import AlnPreprocessArguments
from src.oarsman_dependencies import SamtoolsDependencies
# from src.output_data import AlnPreprocessOutput

from src.runners.bcftools_var_call_runner import run_bcftools_var_call
from src.oarsman_arguments import CallVariantsArguments
from src.oarsman_dependencies import BcfVarCallDependencies
from src.output_data import VariantCallingOutput

from src.runners.highlighter_runner import run_highlighter
from src.output_data import AnnotatedSeq
from src.oarsman_arguments import ConsensAnnotArguments
from src.oarsman_dependencies import HighlighterDependencies



def main():

    # Parse command line arguments
    oarsman_args, oarsman_dependencies = parse_arguments()

    print(oarsman_args.reads_R1_fpaths)
    print(oarsman_args.reads_R2_fpaths)
    # sys.exit(0)

    # Make amplicons (ma) for kromsatel
    ma_arguments: MakeAmpliconsArguments = oarsman_args.get_make_amplicons_args()
    ma_dependencies: MakeAmpliconsDependencies = oarsman_dependencies.get_make_amplicons_dependencies()

    ma_output_data: MakeAmpliconsOutput = run_make_amplicons(
        ma_arguments,
        ma_dependencies
    )

    print(ma_output_data.no_primers_fpath)
    print(ma_output_data.with_primers_fpath)


    # Make-db for kromsatel
    make_db_arguments: MakeDbArguments = oarsman_args.get_make_db_args(
        ma_output_data.no_primers_fpath
    )
    make_db_dependencies: MakeDbDependencies = oarsman_dependencies.get_make_db_dependencies()

    make_db_output_data: MakeDbOutput = run_make_db(
        make_db_arguments,
        make_db_dependencies
    )

    print(make_db_output_data.db_fpath)


    reference_seq_fpath = oarsman_args.ref_genome_seq_fpath

    # From this moment we will iterate over samples

    n_samples = len(oarsman_args.reads_R1_fpaths)
    for i_sample in range(n_samples):

        sample = Sample(oarsman_args.reads_R1_fpaths[i_sample][0], i_sample)

        # sample_name = fastq_fpath_to_sample_name(
        #     oarsman_args.reads_R1_fpaths[i_sample][0]
        # )

        print(
            '  \n|=== Processing sample #{}/{}: {} ===|\n' \
                .format(sample.ordinal_num+1, n_samples, sample.name)
        )

        cleaned_reads: Reads = _clean_reads(
            oarsman_args,
            oarsman_dependencies,
            sample,
            make_db_output_data
        )

        mapping: Mapping = _map_reads(
            oarsman_args,
            oarsman_dependencies,
            cleaned_reads,
            reference_seq_fpath,
            sample,
            'reference'
        )


        # Call variants
        var_call_args: CallVariantsArguments = oarsman_args.get_var_call_args(
            sample.name,
            mapping.file_path
        )
        var_call_dependencies: BcfVarCallDependencies = oarsman_dependencies.get_bcftools_dependencies()

        var_call_output: VariantCalling_Output = run_bcftools_var_call(
            var_call_args,
            var_call_dependencies
        )

        print(f'Consensus: `{var_call_output.consensus_fpath}`')

        print('Mapping the reads to the consensus')
        # Map the reads to the consensus
        mapping_to_consensus: Mapping = _map_reads(
            oarsman_args,
            oarsman_dependencies,
            cleaned_reads,
            var_call_output.consensus_fpath,
            sample,
            'consensus'
        )

        consensus_annotation: AnnotatedSeq = _annotate_consensus(
            oarsman_args,
            oarsman_dependencies,
            sample,
            var_call_output.consensus_fpath,
            mapping_to_consensus
        )

        print('Annotated consensus: `{}`'.format(consensus_annotation.file_path))
    # end for
# end def main


def _clean_reads(oarsman_args, oarsman_dependencies, sample, make_db_output_data):

    # Run kromsatel
    kromsatel_arguments: KromsatelArguments = oarsman_args.get_kromsatel_args(
        sample.ordinal_num,
        make_db_output_data.db_fpath
    )
    kromsatel_dependencies: KromsatelDependencies = oarsman_dependencies.get_kromsatel_dependencies()

    kromsatel_output: KromsatelOutput = run_kromsatel(
        kromsatel_arguments,
        kromsatel_dependencies
    )

    # Match up paired reads
    if not kromsatel_output.reads_R2_fpath is None:
        # After kromsatel paired reads are shuffled and now have to be matched up together
        pair_arguments: PairArguments = PairArguments(
            kromsatel_output.reads_R1_fpath,
            kromsatel_output.reads_R2_fpath
        )
        pair_dependencies: PairDependencies = oarsman_dependencies.get_pair_dependencies()
        pair_output: Reads = run_pair(pair_arguments, pair_dependencies)
    else:
        # If we have unpaired data, "forward" reads will be processed as unpaired from now
        pair_output: Reads = Reads(
            None,
            None,
            [kromsatel_output.reads_R1_fpath]
        )
    # end if

    return pair_output
# end def _clean_reads


def _map_reads(
    oarsman_args,
    oarsman_dependencies,
    reads,
    reference_seq_fpath,
    sample,
    output_suffix
):

    read_map_args = oarsman_args.get_read_mapping_args(
        sample.name,
        reads,
        reference_seq_fpath,
        output_suffix
    )

    # read_map_dependencies = oarsman_dependencies.get_bowtie2_dependencies()

    # read_mapping_output: ReadMappingOutput = run_bowtie2(
    #     read_map_args,
    #     read_map_dependencies
    # )

    read_map_dependencies = oarsman_dependencies.get_bwa_dependencies()

    read_mapping_output: Mapping = run_bwa(
        read_map_args,
        read_map_dependencies
    )


    # Preprocess alignment data
    aln_preproc_args = oarsman_args.get_aln_preprocess_args(
        sample.name,
        read_mapping_output.file_path,
        output_suffix
    )
    aln_preproc_dependencies = oarsman_dependencies.get_samtools_dependencies()

    aln_preproces_output: Mapping = run_aln_preprocess(
        aln_preproc_args,
        aln_preproc_dependencies
    )

    return aln_preproces_output
# end def _map_reads


def _annotate_consensus(oarsman_args, oarsman_dependencies, sample, seq_fpath, mapping):

    cons_annot_args = oarsman_args.get_consens_annot_args(
        sample.name,
        seq_fpath,
        mapping
    )

    cons_annot_dependencies = oarsman_dependencies.get_highlighter_dependencies()

    annotated_seq: AnnotatedSeq = run_highlighter(
        cons_annot_args,
        cons_annot_dependencies
    )

    return annotated_seq
# end def _annotate_consensus
