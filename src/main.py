
import sys

from src.parse_arguments import parse_arguments

from src.runners.make_amplicons_runner import run_make_amplicons
from src.oarsman_arguments import MakeAmpliconsArguments
from src.oarsman_dependencies import MakeAmpliconsDependencies
from src.output_data import MakeAmpliconsOutput

from src.runners.make_db_runner import run_make_db
from src.oarsman_arguments import MakeDbArguments
from src.oarsman_dependencies import MakeDbDependencies
from src.output_data import MakeDbOutput

from src.runners.kromsatel_runner import run_kromsatel
from src.oarsman_arguments import KromsatelArguments
from src.oarsman_dependencies import KromsatelDependencies
from src.output_data import KromsatelOutput

def main():

    # Parse command line arguments
    oarsman_args, oarsman_dependencies = parse_arguments()

    print(oarsman_args.reads_R1_fpaths)
    print(oarsman_args.reads_R2_fpaths)
    # sys.exit(0)

    # Make amplicons (ma) for kromsatel
    ma_arguments: MakeAmpliconsArguments = oarsman_args.get_make_amplicons_args()
    ma_dependencies: MakeAmpliconsDependencies = oarsman_dependencies.get_make_amplicons_dependencies()

    ma_output_data: MakeAmpliconsOutput = run_make_amplicons(
        ma_arguments,
        ma_dependencies
    )

    print(ma_output_data.no_primers_fpath)
    print(ma_output_data.with_primers_fpath)


    # Make-db for kromsatel
    make_db_arguments: MakeDbArguments = oarsman_args.get_make_db_args(
        ma_output_data.no_primers_fpath
    )
    make_db_dependencies: MakeDbDependencies = oarsman_dependencies.get_make_db_dependencies()

    make_db_output_data: MakeDbOutput = run_make_db(
        make_db_arguments,
        make_db_dependencies
    )

    print(make_db_output_data.db_fpath)

    # From this moment we will iterate over samples

    n_samples = len(oarsman_args.reads_R1_fpaths)
    for i_sample in range(n_samples):
        print(f'Doing sample {i_sample+1}...')

        # Run kromsatel
        kromsatel_arguments: KromsatelArguments = oarsman_args.get_kromsatel_args(
            i_sample,
            make_db_output_data.db_fpath
        )
        kromsatel_dependencies: KromsatelDependencies = oarsman_dependencies.get_kromsatel_dependencies()

        kromsatel_output: KromsatelOutput = run_kromsatel(
            kromsatel_arguments,
            kromsatel_dependencies
        )

        print(kromsatel_output.reads_R1_fpath)
        print(kromsatel_output.reads_R2_fpath)
    # end for
# end def main
