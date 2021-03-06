
from src.parse_arguments import parse_arguments
from src.seq_sample import Sample

from src.runners.kromsatel import run_kromsatel
from src.runners.bwa import run_bwa
from src.runners.aln_postprocess import run_aln_postprocess
from src.runners.bcftools_var_call import run_bcftools_var_call
from src.runners.lofreq_var_call import run_lofreq_var_call
from src.runners.var_call_postprocess import run_var_vall_postprocess
from src.runners.consensus_making import make_consensus
from src.runners.highlighter import run_highlighter


def run_pipeline():

    # Parse command line arguments
    oarsman_args, oarsman_dependencies = parse_arguments()
    oarsman_dependencies.check_versions()

    # From this moment we will iterate over samples
    n_samples = oarsman_args.get_number_of_samples()

    for i_sample in range(n_samples):

        sample = Sample(oarsman_args.reads_R1_fpaths[i_sample][0], i_sample)

        print(
            '  \n|=== Processing sample #{}/{}: {} ===|' \
                .format(sample.ordinal_num+1, n_samples, sample.name)
        )

        cleaned_reads = _clean_reads(
            oarsman_args,
            oarsman_dependencies,
            sample
        )

        print('Mapping the reads to the reference')
        mapping = _map_reads(
            oarsman_args,
            oarsman_dependencies,
            cleaned_reads,
            oarsman_args.reference_fpath,
            sample,
            'reference'
        )

        consensus_seq = _create_consensus_seq(
            mapping,
            oarsman_args,
            oarsman_dependencies,
            sample
        )

        print('Mapping the reads to the consensus')
        mapping_to_consensus = _map_reads(
            oarsman_args,
            oarsman_dependencies,
            cleaned_reads,
            consensus_seq.file_path,
            sample,
            'consensus'
        )

        consensus_annotation = _annotate_consensus(
            oarsman_args,
            oarsman_dependencies,
            sample,
            consensus_seq.file_path,
            mapping_to_consensus
        )

        print('Annotated consensus: `{}`'.format(consensus_annotation.file_path))
        print('-' * 40)
    # end for
# end def


def _clean_reads(oarsman_args, oarsman_dependencies, sample):

    kromsatel_arguments = \
        oarsman_args.get_kromsatel_args(sample.ordinal_num)
    kromsatel_dependencies = \
        oarsman_dependencies.get_kromsatel_dependencies()

    kromsatel_result = run_kromsatel(
        kromsatel_arguments,
        kromsatel_dependencies
    )
    print() # separate kromsatel stdout from further stdout

    return kromsatel_result
# end def


def _map_reads(oarsman_args,
               oarsman_dependencies,
               reads,
               reference_fpath,
               sample,
               output_suffix):
    raw_mapping = _run_mapping_program(
        oarsman_args,
        oarsman_dependencies,
        reads,
        reference_fpath,
        sample,
        output_suffix
    )

    postprocessed_mapping = _postprocess_mapping(
        raw_mapping,
        oarsman_args,
        oarsman_dependencies,
        reference_fpath,
        sample,
        output_suffix
    )

    return postprocessed_mapping
# end def


def _run_mapping_program(oarsman_args,
                         oarsman_dependencies,
                         reads,
                         reference_fpath,
                         sample,
                         output_suffix):
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

    return read_mapping
# end def


def _postprocess_mapping(raw_mapping,
                         oarsman_args,
                         oarsman_dependencies,
                         reference_fpath,
                         sample,
                         output_suffix):
    aln_postproc_args = oarsman_args.get_aln_postprocess_args(
        sample.name,
        raw_mapping,
        reference_fpath,
        output_suffix
    )
    aln_postproc_dependencies = oarsman_dependencies.get_samtools_dependencies()

    postprocessed_mapping = run_aln_postprocess(
        aln_postproc_args,
        aln_postproc_dependencies
    )

    return postprocessed_mapping
# end def


def _create_consensus_seq(mapping,
                          oarsman_args,
                          oarsman_dependencies,
                          sample):

    var_call = _call_variants(mapping, oarsman_args, oarsman_dependencies, sample)

    print('Making consensus...')

    make_consensus_args = oarsman_args.get_make_consensus_args(
        sample.name,
        var_call
    )
    make_consensus_dependencies = \
        oarsman_dependencies.get_bcftools_dependencies()

    consensus_seq = make_consensus(
        make_consensus_args,
        make_consensus_dependencies
    )

    print(f'The consensus is created: `{consensus_seq.file_path}`')

    return consensus_seq
# end def


def _call_variants(mapping,
                   oarsman_args,
                   oarsman_dependencies,
                   sample):

    if oarsman_dependencies.lofreq_fpath is None:
        var_call_args = oarsman_args.get_var_call_args(
            sample.name,
            mapping.alignment_fpath
        )
        var_call_dependencies = \
            oarsman_dependencies.get_bcftools_dependencies()

        raw_var_call = run_bcftools_var_call(
            var_call_args,
            var_call_dependencies
        )

        var_call_postproc_args = oarsman_args.get_var_call_postproc_args(
            sample.name,
            raw_var_call
        )
        var_call_postproc_dependencies = \
            oarsman_dependencies.get_bcftools_dependencies()

        var_call = run_var_vall_postprocess(
            var_call_postproc_args,
            var_call_postproc_dependencies
        )
    else:
        var_call_args = oarsman_args.get_var_call_args(
            sample.name,
            mapping.alignment_fpath
        )
        var_call_dependencies = \
            oarsman_dependencies.get_lofreq_dependencies()

        var_call = run_lofreq_var_call(
            var_call_args,
            var_call_dependencies
        )
    # end if

    return var_call
# end def


def _annotate_consensus(oarsman_args,
                        oarsman_dependencies,
                        sample,
                        seq_fpath,
                        mapping):
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
