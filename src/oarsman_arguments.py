
import os

from src.filesystem import rm_fasta_extention

class OarsmanArguments:

    def __init__(self):

        # Input data
        # Reads
        self.reads_R1_fpaths = list()
        self.reads_R2_fpaths = list()
        # Reference genome
        self.ref_genome_seq_fpath = None
        # Primers
        self.primers_fpath = None

        # Arguments for kromsatel and co
        self.min_major_len = 100 # bp
        self.min_minor_len = 25  # bp
        self.chunk_size = 1000   # reads

        # Variant calling
        self.min_variant_qual = 20 # in Phred scale

        # Misc
        self.outdir_path = os.path.join(os.getcwd(), 'oarsman_outdir')
        self.n_threads = 1 # thread
    # end def __init__

    def get_make_amplicons_args(self):

        return MakeAmpliconsArguments(
            self.primers_fpath,
            self.ref_genome_seq_fpath,
            os.path.join(self.outdir_path, 'reference_amplicons', 'reference_amplicons'),
            self.min_minor_len
        )
    # end def get_make_amplicons_args

    def get_make_db_args(self, input_seqs_fpath):

        return MakeDbArguments(
            input_seqs_fpath,
            os.path.join(self.outdir_path, 'db_for_kromsatel')
        )
    # end def get_make_db_args

    def get_kromsatel_args(self, i_sample, db_dir_path):

        unpaired_data = len(self.reads_R2_fpaths) == 0

        if unpaired_data:
            reads_R2_fpaths = list()
        else:
            reads_R2_fpaths = self.reads_R2_fpaths[i_sample]
        # end if

        return KromsatelArguments(
            self.reads_R1_fpaths[i_sample],
            reads_R2_fpaths,
            self.primers_fpath,
            db_dir_path,
            os.path.join(self.outdir_path, 'kromsatel_outdir'),
            self.min_major_len,
            self.min_minor_len,
            self.chunk_size,
            self.n_threads
        )
    # end def get_kromsatel_args

    def get_read_mapping_args(
        self,
        sample_name,
        reads_R1_fpath,
        reads_R2_fpath,
        unpaired_reads_fpaths
    ):

        ref_genome_basename = os.path.basename(self.ref_genome_seq_fpath)

        index_base_fpath = os.path.join(
            os.path.dirname(self.ref_genome_seq_fpath),
            rm_fasta_extention(ref_genome_basename) + 'index'
        )

        return ReadMappingArguments(
            sample_name,
            reads_R1_fpath,
            reads_R2_fpath,
            unpaired_reads_fpaths,
            self.ref_genome_seq_fpath,
            index_base_fpath,
            os.path.join(self.outdir_path, 'read_mappings'),
            self.n_threads
        )
    # end def get_read_mapping_args

    def get_var_call_args(self, sample_name, alignment_fpath):

        var_calls_dirpath = os.path.join(
            self.outdir_path,
            'variant_calls'
        )

        consensus_dirpath = os.path.join(
            self.outdir_path,
            'consensus'
        )

        return CallVariantsArguments(
            sample_name,
            alignment_fpath,
            self.ref_genome_seq_fpath,
            self.min_variant_qual,
            var_calls_dirpath,
            consensus_dirpath
        )

    # end def get_var_call_args
# end class OarsmanArguments


class MakeAmpliconsArguments:

    def __init__(
        self,
        primers_fpath,
        ref_genome_seq_fpath,
        prefix,
        min_minor_len
    ):

        self.primers_fpath = primers_fpath
        self.ref_genome_seq_fpath = ref_genome_seq_fpath
        self.prefix = prefix
        self.min_minor_len = min_minor_len
    # end def __init__
# end class MakeAmpliconsArguments


class MakeDbArguments:

    def __init__(
        self,
        amplicons_seqs_fpath,
        db_dir_path
    ):

        self.amplicons_seqs_fpath = amplicons_seqs_fpath
        self.db_dir_path = db_dir_path
    # end def __init__
# end class MakeDbArguments


class KromsatelArguments:

    def __init__(
        self,
        reads_R1_fpaths,
        reads_R2_fpaths,
        primers_fpath,
        db_dir_path,
        outdir,
        min_major_len,
        min_minor_len,
        chunk_size,
        n_threads
    ):

        self.reads_R1_fpaths = reads_R1_fpaths
        self.reads_R2_fpaths = reads_R2_fpaths
        self.primers_fpath = primers_fpath
        self.db_dir_path = db_dir_path
        self.outdir = outdir
        self.min_major_len = min_major_len
        self.min_minor_len = min_minor_len
        self.chunk_size = chunk_size
        self.n_threads = n_threads
    # end def __init__
# end class KromsatelArguments


class PairArguments:

    def __init__(
        self,
        reads_R1_fpath,
        reads_R2_fpath
    ):

        self.reads_R1_fpath = reads_R1_fpath
        self.reads_R2_fpath = reads_R2_fpath
    # end def __init__
# end class PairArguments


class ReadMappingArguments:

    def __init__(
        self,
        sample_name,
        reads_R1_fpath,
        reads_R2_fpath,
        unpaired_reads_fpaths,
        ref_genome_fpath,
        genome_index_base_fpath,
        outdir_path,
        n_threads
    ):

        self.sample_name = sample_name
        self.reads_R1_fpath = reads_R1_fpath
        self.reads_R2_fpath = reads_R2_fpath
        self.unpaired_reads_fpaths = unpaired_reads_fpaths
        self.ref_genome_fpath = ref_genome_fpath
        self.genome_index_base_fpath = genome_index_base_fpath
        self.outdir_path = outdir_path
        self.n_threads = n_threads
    # end def __init__
# end class ReadMappingArguments


class CallVariantsArguments:

    def __init__(
        self,
        sample_name,
        alignment_fpath,
        reference_fpath,
        min_variant_qual,
        var_calls_dirpath,
        consensus_dirpath
    ):

        self.sample_name = sample_name
        self.alignment_fpath = alignment_fpath
        self.reference_fpath = reference_fpath
        self.min_variant_qual = min_variant_qual
        self.var_calls_dirpath = var_calls_dirpath
        self.consensus_dirpath = consensus_dirpath
    # end def __init__
# end class CallVariantsArguments
