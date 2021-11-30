
import sys
import subprocess as sp

from src.oarsman_arguments import MakeAmpliconsArguments
from src.oarsman_dependencies import MakeAmpliconsDependencies
from src.output_data import MakeAmpliconsOutput

from src.filesystem import check_files_exist



def _configure_make_amplicons_command(
    args: MakeAmpliconsArguments,
    dependencies: MakeAmpliconsDependencies
):

    command = ' '.join(
        [
            'bash',
            dependencies.make_amplicons_fpath,
            f'-p {args.primers_fpath}',
            f'-g {args.ref_genome_seq_fpath}',
            f'-o {args.prefix}',
            f'-s {dependencies.seqkit_fpath}',
            f'-i {args.min_amplicon_len}',
        ]
    )

    return command
# end def _configure_make_amplicons_command


def run_make_amplicons(args: MakeAmpliconsArguments, dependencies: MakeAmpliconsDependencies):
    command_str = _configure_make_amplicons_command(args, dependencies)
    
    pipe = sp.Popen(command_str, shell=True, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Error!')
        print(f'Script `make-amplicons.sh` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if

    no_primers_fpath = args.prefix + '.fasta'
    with_primers_fpath = args.prefix + '_with-primers.fasta'

    # Check if all output files exist
    non_extant_fpaths = check_files_exist(no_primers_fpath, with_primers_fpath)
    if len(non_extant_fpaths) != 0:
        for i, fpath in enumerate(non_extant_fpaths):
            print(f'\nError #{i+1}: file `{fpath}` does not exist after make-amplicons script has ended.')
            print('This file is the output of the script make-amplicons, so it must exist.')
        # end for
        sys.exit(1)
    # end if

    return MakeAmpliconsOutput(
        no_primers_fpath,
        with_primers_fpath
    )
# end def _make_amplicons
