
import os

from src.filesystem import rm_fastq_extention


class Sample:
    def __init__(self, reads_fpath, ordinal_num):
        self.name = _fastq_fpath_to_sample_name(reads_fpath)
        self.ordinal_num = ordinal_num # 0-based
    # end def __init__
# end class Sample


def _fastq_fpath_to_sample_name(fpath):

    basename_no_extention = os.path.basename(
        rm_fastq_extention(fpath)
    )

    # Remove _R1_001 and/or _R2_001 substrings from `basename_no_extention`,
    #   and here we have the sample name
    direction_str = None
    if '_R1_' in basename_no_extention:
        direction_str = '_R1_'
    elif '_R2_' in basename_no_extention:
        direction_str = '_R2_'
    # end if

    if not direction_str is None and basename_no_extention.find(direction_str) != 0:
        sample_name = basename_no_extention.partition(direction_str)[0]
    # end if

    return sample_name
# end def _fastq_fpath_to_sample_name
