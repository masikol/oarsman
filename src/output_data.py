
class MakeAmpliconsOutput:

    def __init__(self, no_primers_fpath, with_primers_fpath):
        self.no_primers_fpath = no_primers_fpath
        self.with_primers_fpath = with_primers_fpath
    # end __init__
# end class MakeAmpliconsOutput


class MakeDbOutput:

    def __init__(self, db_fpath):
        self.db_fpath = db_fpath
    # end __init__
# end class MakeAmpliconsOutput


class KromsatelOutput:

    def __init__(self, reads_R1_fpath, reads_R2_fpath):
        self.reads_R1_fpath = reads_R1_fpath
        self.reads_R2_fpath = reads_R2_fpath
    # end __init__
# end class KromsatelOutput


class PairOutput:

    def __init__(
        self,
        reads_R1_paired_fpath, reads_R2_paired_fpath,
        unpaired_reads_fpaths
    ):
        self.reads_R1_paired_fpath = reads_R1_paired_fpath
        self.reads_R2_paired_fpath = reads_R2_paired_fpath
        self.unpaired_reads_fpaths = unpaired_reads_fpaths
    # end __init__
# end class PairOutput


class ReadMappingOutput:
    def __init__(self, alignment_fpath):
        self.alignment_fpath = alignment_fpath
    # end def __init__
# end class ReadMappingOutput


class AlnPreprocessOutput:
    def __init__(self, preproc_aln_fpath):
        self.preproc_aln_fpath = preproc_aln_fpath
    # end def __init__
# end class AlnPreprocessOutput


class VariantCallingOutput:
    def __init__(self, consensus_fpath):
        self.consensus_fpath = consensus_fpath
    # end def __init__
# end class VariantCallingOutput