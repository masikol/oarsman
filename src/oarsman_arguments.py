
import os

class OarsmanArguments:

    def __init__(self):

        # Own arguments
        self.tmp_dir_path = os.path.join(os.getcwd(), 'oarsman_tmpdir')

        # Arguments for make-amplicons
        self.primers_fpath = None
        self.ref_genome_seq_fpath = None
        self.min_amplicon_len = 25
    # end def __init__

    def get_make_amplicons_args(self):

        return MakeAmpliconsArguments(
            self.primers_fpath,
            self.ref_genome_seq_fpath,
            os.path.join(self.tmp_dir_path, 'reference_amplicons', 'reference_amplicons'),
            self.min_amplicon_len
        )
    # end def get_make_amplicons_args

    def get_make_db_args(self, input_seqs_fpath):

        return MakeDbArguments(
            input_seqs_fpath,
            os.path.join(self.tmp_dir_path, 'db_for_kromsatel')
        )
    # end def get_make_db_args

# end class OarsmanArguments


class MakeAmpliconsArguments:

    def __init__(
        self,
        primers_fpath,
        ref_genome_seq_fpath,
        prefix='my-amplicons',
        min_amplicon_len=25
    ):

        self.primers_fpath = primers_fpath
        self.ref_genome_seq_fpath = ref_genome_seq_fpath
        self.prefix = prefix
        self.min_amplicon_len = min_amplicon_len
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
