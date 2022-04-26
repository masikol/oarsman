
import os

from src.shell import launch_command

from src.data_transfer_objects import SequenceFile
from src.arguments import ConsensAnnotArguments
from src.dependencies import HighlighterDependencies


def run_highlighter(args, dependencies):

    print('Annotating the consensus...')
    command_str = _configure_highlighter_command(args, dependencies)
    launch_command(command_str, 'consensus-highlighter')

    annotated_seq = SequenceFile(args.outfpath)
    annotated_seq.check_existance()

    return annotated_seq
# end def


def _configure_highlighter_command(highlighter_args, dependencies):

    if os.path.isfile(dependencies.highlighter_fpath):
        highlighter_script_str = 'python3 {}'.format(dependencies.highlighter_fpath)
    else:
        highlighter_script_str = dependencies.highlighter_fpath
    # end if

    command = ' '.join(
        [
            highlighter_script_str,
            f'-f {highlighter_args.seq_fpath}',
            f'-b {highlighter_args.mapping.alignment_fpath}',
            f'-o {highlighter_args.outfpath}',
            highlighter_args.advanced_args
            # '-c {}'.format(','.join(highlighter_args.low_coverages))
        ]
    )

    return command
# end def
