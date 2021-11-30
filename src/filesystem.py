
import os
import re


def _file_does_not_exist(fpath):
    return not os.path.isfile(fpath)
# end def _file_does_not_exist


def check_files_exist(*fpaths):

    non_extant_files = tuple(
        filter(
            _file_does_not_exist,
            fpaths
        )
    )

    return non_extant_files
# end def check_files_exist


def util_is_in_path(util_name):
    util_found = False

    for dir_from_path in os.environ['PATH'].split(os.pathsep):
        if util_name in os.listdir(dir_from_path):
            util_found = True
            break
        # end if
    # end for

    return util_found
# end def util_is_in_path


def rm_fasta_extention(fpath):

    fasta_extention_pattern = r'^.+(\.fa(sta)?)$'
    extention_match = re.match(fasta_extention_pattern, fpath)

    if not extention_match is None:
        fasta_extention = extention_match.group(1)
        new_fpath = os.path.basename(fpath).replace(fasta_extention, '')
    else:
        new_fpath = os.path.basename(fpath)
    # end if

    return os.path.join(
        os.path.dirname(new_fpath),
        new_fpath
    )
# end def rm_fasta_extention
