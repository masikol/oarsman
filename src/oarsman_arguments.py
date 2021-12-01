
import os

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

        # Misc
        self.tmp_dir_path = os.path.join(os.getcwd(), 'oarsman_tmpdir')
        self.n_threads = 1 # thread
    # end def __init__

    def get_make_amplicons_args(self):

        return MakeAmpliconsArguments(
            self.primers_fpath,
            self.ref_genome_seq_fpath,
            os.path.join(self.tmp_dir_path, 'reference_amplicons', 'reference_amplicons'),
            self.min_minor_len
        )
    # end def get_make_amplicons_args

    def get_make_db_args(self, input_seqs_fpath):

        return MakeDbArguments(
            input_seqs_fpath,
            os.path.join(self.tmp_dir_path, 'db_for_kromsatel')
        )
    # end def get_make_db_args

    def get_kromsatel_args(self, i_sample, db_dir_path):
        return KromsatelArguments(
            self.reads_R1_fpaths[i_sample],
            self.reads_R2_fpaths[i_sample],
            self.primers_fpath,
            db_dir_path,
            os.path.join(self.tmp_dir_path, 'kromsatel_outdir'),
            self.min_major_len,
            self.min_minor_len,
            self.chunk_size,
            self.n_threads
        )
    # end def get_kromsatel_args
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
