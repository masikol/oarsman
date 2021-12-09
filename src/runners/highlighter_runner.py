
import os
import sys
import subprocess as sp


from src.output_data import AnnotatedSeq
from src.oarsman_arguments import ConsensAnnotArguments
from src.oarsman_dependencies import HighlighterDependencies


def _configure_highlighter_command(
        args: ConsensAnnotArguments,
        dependencies: HighlighterDependencies
    ):

    command = ' '.join(
        [
            'python3', dependencies.highlighter_fpath,
            f'-f {args.seq_fpath}',
            f'-b {args.mapping.file_path}',
            f'-o {args.outfpath}',
            '-c {}'.format(','.join(args.low_coverages))
        ]
    )

    return command
# end def _configure_bwa_index_command


def run_highlighter(args, dependencies):

    command_str = _configure_highlighter_command(args, dependencies)

    print('Annotating the consensus...')
    pipe = sp.Popen(command_str, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError!')
        print(f'Script `consensus-highlighter.py` returned a non-zero exit code: {pipe.returncode}')
        print('Error message:')
        stderr_index = 1
        print(stdout_stderr[stderr_index].decode('utf-8'))
        sys.exit(1)
    # end if

    if not os.path.exists(args.outfpath):
        print(f'\nError: GenBank file `{args.outfpath}` does not exist after annotation')
        print('This file must exist, though. Exitting...')
        sys.exit(1)
    # end if

    return AnnotatedSeq(args.outfpath)
# end def run_highlighter
