
import os
import sys
import subprocess as sp

from src.oarsman_arguments import MakeDbArguments
from src.oarsman_dependencies import MakeDbDependencies
from src.output_data import MakeDbOutput

from src.filesystem import rm_fasta_extention


def _configure_make_db_command(
    args: MakeDbArguments,
    dependencies: MakeDbDependencies
):

    command = ' '.join(
        [
            'bash',
            dependencies.make_db_path,
            args.amplicons_seqs_fpath,
            args.db_dir_path
        ]
    )

    return command
# end def _configure_make_db_command


def run_make_db(args: MakeDbArguments, dependencies: MakeDbDependencies):

    command_str = _configure_make_db_command(args, dependencies)

    pipe = sp.Popen(command_str, shell=True, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Error!')
        print(f'Script `make-db.sh` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if

    # Print an empty line after makeblastdb output
    print()

    db_path = os.path.join(
        args.db_dir_path,
        os.path.basename(rm_fasta_extention(args.amplicons_seqs_fpath))
    )

    return MakeDbOutput(db_path)
# end def run_make_db
