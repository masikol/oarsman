
import os
import subprocess as sp

from src.oarsman_arguments import KromsatelArguments
from src.oarsman_dependencies import KromsatelDependencies
from src.output_data import KromsatelOutput

from src.filesystem import rm_fastq_extention, gzip_file, rm_file


def _configure_kromsatel_command(
        read_fpath: str,
        args: KromsatelArguments,
        dependencies: KromsatelDependencies
    ):

    command = ' '.join(
        [
            'python3', dependencies.kromsatel_fpath,
            read_fpath,
            f'-p {args.primers_fpath}',
            f'-d {args.db_dir_path}',
            f'--am {args.min_major_len}',
            f'--im {args.min_minor_len}',
            f'-c {args.chunk_size}',
            f'-t {args.n_threads}',
            f'-o {args.outdir}'
        ]
    )

    return command
# end def _configure_make_amplicons_command


def _make_outfpath(fq_fpath, outdir):
    # Function configures path to output file.
    #
    # :param fq_fpath: path to input fastq file;
    # :type fq_fpath: str;
    # :param outdir: path to output directory;
    # :type outdir: str;
    #
    # Returns path to output file (str).

    name_with_no_ext = rm_fastq_extention(fq_fpath)

    return os.path.join(
        outdir,
        '{}_cleaned.fastq'.format(name_with_no_ext)
    )
# end def _make_outfpath


def run_kromsatel(args, dependencies):

    cleaned_reads_fpaths = list()

    for read_fpath_collection in (args.reads_R1_fpaths, args.reads_R2_fpaths):
        gzipped_outfile_paths = list()

        for read_fpath in read_fpath_collection:

            command_str = _configure_kromsatel_command(read_fpath, args, dependencies)

            pipe = sp.Popen(command_str, shell=True, stderr=sp.PIPE)
            stdout_stderr = pipe.communicate()

            if pipe.returncode != 0:
                print('\nError!')
                print(f'Script `kromsatel.py` returned a non-zero exit code: {pipe.returncode}')
                print('Error message:')
                stderr_index = 1
                print(stdout_stderr[stderr_index].decode('utf-8'))
                sys.exit(1)
            # end if

            outfpath = _make_outfpath(read_fpath, args.outdir)

            print(f'Gzipping file `{outfpath}`')
            gzipped_outfile_paths.append(gzip_file(outfpath))
        # end for

        if len(gzipped_outfile_paths) > 0:

            merged_sample_fpath = gzipped_outfile_paths[0]

            if len(gzipped_outfile_paths) > 1:

                kbyte_size = 1024 # bytes
                buffer_size = 512 * kbyte_size

                print(f'Merging all sample output files into a single file: `{merged_sample_fpath}`')

                with open(merged_sample_fpath, 'ab') as merged_sample_file:
                    for partial_outfpath in gzipped_outfile_paths[1:]:
                        with open(partial_outfpath, 'rb') as partial_outfile:
                            chunk = b'non_empty_bytes'
                            while chunk != b'':
                                chunk = partial_outfile.read(buffer_size)
                                merged_sample_file.write(chunk)
                            # end while
                    # end for
                    rm_file(partial_outfpath)
                # end with
            # end if

            cleaned_reads_fpaths.append(merged_sample_fpath)
        # end if

        print()
    # end for

    forward_output_read_fpath = cleaned_reads_fpaths[0]
    reverse_output_read_fpath = cleaned_reads_fpaths[1]

    return KromsatelOutput(
        forward_output_read_fpath,
        reverse_output_read_fpath
    )
# end def run_kromnsatel
