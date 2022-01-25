
import os
import src.filesystem as fs
from src.fatal_errors import FatalError


class DataTransferObject:
    def check_existance(self):
        raise NotImplementedError
    # end def
# end class


class KromsatelResult(DataTransferObject):
    def append_to(self, dest_output):
        raise NotImplementedError
    # end def
    def cleanup(self):
        raise NotImplementedError
    # end def
# end class


class UnpairedKromsatelResult(KromsatelResult):

    def __init__(self, reads_fpath):
        self.reads_fpath = reads_fpath
    # end def

    def append_to(dest_output):
        fs.append_file_to_file(self.reads_fpath, dest_output.reads_fpath)
    # end def

    def cleanup(self):
        fs.rm_file(self.reads_fpath)
    # end def

    def check_existance(self):
        if not os.path.exists(self.reads_fpath):
            error_msg = '\nError: file `{}` must exist after kromsatel work' \
                ' but it does not.\nPlease, contact the developer.' \
                .format(self.reads_fpath)
            raise FatalError(error_msg)
        # end if
    # end def
# end class


class PairedKromsatelResult(KromsatelResult):

    def __init__(self,
                 reads_frw_fpath,
                 reads_rvr_fpath,
                 reads_frw_upr_fpaths,
                 reads_rvr_upr_fpaths):
        self.reads_frw_fpath = reads_frw_fpath
        self.reads_rvr_fpath = reads_rvr_fpath
        self.reads_frw_upr_fpaths = reads_frw_upr_fpaths
        self.reads_rvr_upr_fpaths = reads_rvr_upr_fpaths
    # end def

    def append_to(self, dest_output):
        fs.append_file_to_file(self.reads_frw_fpath, dest_output.reads_frw_fpath)
        fs.append_file_to_file(self.reads_rvr_fpath, dest_output.reads_rvr_fpath)
        fs.append_file_to_file(self.reads_frw_upr_fpaths, dest_output.reads_frw_upr_fpaths)
        fs.append_file_to_file(self.reads_rvr_upr_fpaths, dest_output.reads_rvr_upr_fpaths)
    # end def

    def cleanup(self):
        fs.rm_file(self.reads_frw_fpath)
        fs.rm_file(self.reads_rvr_fpath)
        fs.rm_file(self.reads_frw_upr_fpaths)
        fs.rm_file(self.reads_rvr_upr_fpaths)
    # end def

    def check_existance(self):

        file_paths = (
            self.reads_frw_fpath, self.reads_rvr_fpath,
            self.reads_frw_upr_fpaths, self.reads_rvr_upr_fpaths,
        )

        for fpath in file_paths:
            if not os.path.exists(fpath):
                error_msg = '\nError: file `{}` must exist after kromsatel work' \
                    ' but it does not.\nPlease, contact the developer.' \
                    .format(fpath)
                raise FatalError(error_msg)
            # end if
        # end for
    # end def
# end class


class SeqIndex(DataTransferObject):

    def __init__(self, index_base_fpath):
        self.index_base_fpath = index_base_fpath
    # end def

    def check_existance(self):
        raise NotImplementedError
    # end def
# end class


class BwaSeqIndex(SeqIndex):

    bwa_index_extentions = ('.amb', '.ann','.bwt', '.pac', '.sa',)

    def check_existance(self):
        for extention in self.bwa_index_extentions:
            self._check_index_part(extention)
        # end if
    # end def

    def _check_index_part(self, index_extention):
        index_fpath = self.index_base_fpath + index_extention
        if not os.path.exists(index_fpath):
            error_msg = '\nError: a fasta index file `{}` must exist after indexing' \
                ' with bwa-build, but it does not.\nPlease, contact the developer.' \
                .format(index_fpath)
            raise FatalError(error_msg)
        # end if
    # end def
# end class


class ReadMapping(DataTransferObject):

    def __init__(self, alignment_fpath):
        self.alignment_fpath = alignment_fpath
    # end def __init__

    def check_existance(self):
        if not os.path.exists(self.alignment_fpath):
            error_msg = '\nError: file `{}` must exist after read mapping' \
                ' but it does not.\nPlease, contact the developer.' \
                .format(self.alignment_fpath)
            raise FatalError(error_msg)
        # end if
    # end def

    def cleanup(self):
        fs.rm_file(self.alignment_fpath)
    # end def
# end class


class AlignmentIndex(DataTransferObject):

    index_extention = '.bai'

    def __init__(self, alignment_fpath):
        self.alignment_fpath = alignment_fpath
    # end def

    def check_existance(self):
        index_fpath = self.alignment_fpath + self.index_extention
        if not os.path.exists(index_fpath):
            error_msg = '\nError: file `{}` must exist after BAM indexing' \
                ' but it does not.\nPlease, contact the developer.' \
                .format(index_fpath)
            raise FatalError(error_msg)
        # end if
    # end def
# end class


class VariantCall(DataTransferObject):

    def __init__(self, var_call_fpath):
        self.var_call_fpath = var_call_fpath
    # end def

    def check_existance(self):
        if not os.path.exists(self.var_call_fpath):
            error_msg = '\nError: file `{}` must exist after variant calling' \
                ' but it does not.\nPlease, contact the developer.' \
                .format(self.var_call_fpath)
            raise FatalError(error_msg)
        # end if
    # end def
# end class


class VariantCallIndex(DataTransferObject):

    index_extention = '.csi'

    def __init__(self, var_call_fpath):
        self.var_call_fpath = var_call_fpath
    # end def

    def check_existance(self):
        index_fpath = self.var_call_fpath + self.index_extention
        if not os.path.exists(index_fpath):
            error_msg = '\nError: file `{}` must exist after indexing' \
                ' but it does not.\nPlease, contact the developer.' \
                .format(index_fpath)
            raise FatalError(error_msg)
        # end if
    # end def
# end class


class SequenceFile(DataTransferObject):

    def __init__(self, file_path):
        self.file_path = file_path
    # end def

    def check_existance(self):
        if not os.path.exists(self.file_path):
            error_msg = '\nError: file `{}` must exist after annotation' \
                ' but it does not.\nPlease, contact the developer.' \
                .format(self.file_path)
            raise FatalError(error_msg)
        # end if
    # end def
# end class
