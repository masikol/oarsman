
import os
import re
import gzip
import shutil

from src.fatal_errors import FatalError


def create_directory(dir_path):
    if not os.path.isdir(dir_path):
        try:
            os.makedirs(dir_path)
        except OSError as err:
            error_msg = '\nError: cannot create directory `{}`:\n{}' \
                .format(dir_path, str(err))
            raise FatalError(error_msg)
        # end try
    # end if
# end def


def check_files_exist(*fpaths):

    non_extant_files = tuple(
        filter(
            _file_does_not_exist,
            fpaths
        )
    )

    return non_extant_files
# end def check_files_exist


def _file_does_not_exist(fpath):
    return not os.path.isfile(fpath)
# end def _file_does_not_exist


def abspath_collection(file_path_collection):
    return list(
        map(
            os.path.abspath,
            file_path_collection
        )
    )
# end def

def basename_collection(file_path_collection):
    return list(
        map(
            os.path.basename,
            file_path_collection
        )
    )
# end def


def util_is_in_path(util_name):
    util_found = False

    for dir_from_path in os.environ['PATH'].split(os.pathsep):
        if os.path.isdir(dir_from_path):
            if util_name in os.listdir(dir_from_path):
                util_found = True
                break
            # end if
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


def rm_fastq_extention(fpath):

    fastq_extention_pattern = r'^.+(\.f(ast)?q(\.gz)?)$'
    extention_match = re.match(fastq_extention_pattern, fpath)

    if not extention_match is None:
        fastq_extention = extention_match.group(1)
        new_fpath = os.path.basename(fpath).replace(fastq_extention, '')
    else:
        new_fpath = os.path.basename(fpath)
    # end if

    return os.path.join(
        os.path.dirname(new_fpath),
        new_fpath
    )
# end def rm_fasta_extention


def gzip_file(fpath):

    gzipped_fpath = fpath + '.gz'

    try:
        with open(fpath, 'rb') as infile, \
             gzip.open(gzipped_fpath, 'wb') as gzipped_file:
                shutil.copyfileobj(infile, gzipped_file)
        # with gzip
    except OSError as err:
        print(f'Error: cannot gzip file `{fpath}`')
        print(str(err))
    # end try

    try:
        os.unlink(fpath)
    except OSError as err:
        print(f'Error: cannot remove file `{fpath}` after gzipping')
        print(str(err))
        sys.exit(1)
    # end try

    return gzipped_fpath
# end def gzip_file


def rm_file(fpath):
    try:
        os.unlink(fpath)
    except OSError as err:
        error_msg = 'Error: cannot remove file `{}`:\n  {}' \
            .format(fpath, str(err))
        raise FatalError(error_msg)
    # end try
# end def


def rm_if_exists(fpath):
    if os.path.exists(fpath):
        rm_file(fpath)
    # end if
# end def


def append_file_to_file(source_fpath, dest_fpath, buffer_size_kbytes=512):

    kbyte_size = 1024 # bytes
    buffer_size = buffer_size_kbytes * kbyte_size

    with open(source_fpath, 'rb') as source_file, \
         open(dest_fpath, 'ab') as dest_file:

        chunk = b'non_empty_bytes'
        while chunk != b'':
            chunk = source_file.read(buffer_size)
            dest_file.write(chunk)
        # end while
    # end with
# end def
