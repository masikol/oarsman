
from src.parse_arguments import parse_arguments

from src.runners.make_amplicons_runner import run_make_amplicons
from src.output_data import MakeAmpliconsOutput


def main():

    # Parse command line arguments
    oarsman_arguments, oarsman_dependencies = parse_arguments()

    # Make amplicons for kromsatel
    make_amplicons_arguments = oarsman_arguments.get_make_amplicons_args()
    make_amplicons_dependencies = oarsman_dependencies.get_make_amplicons_dependencies()

    make_amplicons_output_data = run_make_amplicons(
        make_amplicons_arguments,
        make_amplicons_dependencies
    )


    print(make_amplicons_output_data.no_primers_fpath)
    print(make_amplicons_output_data.with_primers_fpath)
# end def main
