
from src.parse_arguments import parse_arguments

from src.runners.make_amplicons_runner import run_make_amplicons
from src.oarsman_arguments import MakeAmpliconsArguments
from src.oarsman_dependencies import MakeAmpliconsDependencies
from src.output_data import MakeAmpliconsOutput

from src.runners.make_db_runner import run_make_db
from src.oarsman_arguments import MakeDbArguments
from src.oarsman_dependencies import MakeDbDependencies
from src.output_data import MakeDbOutput

def main():

    # Parse command line arguments
    oarsman_arguments, oarsman_dependencies = parse_arguments()

    # Make amplicons (ma) for kromsatel
    ma_arguments: MakeAmpliconsArguments = oarsman_arguments.get_make_amplicons_args()
    ma_dependencies: MakeAmpliconsDependencies = oarsman_dependencies.get_make_amplicons_dependencies()

    ma_output_data: MakeAmpliconsOutput = run_make_amplicons(
        ma_arguments,
        ma_dependencies
    )

    print(ma_output_data.no_primers_fpath)
    print(ma_output_data.with_primers_fpath)


    # Make-db for kromsatel
    make_db_arguments: MakeDbArguments = oarsman_arguments.get_make_db_args(
        ma_output_data.no_primers_fpath
    )
    make_db_dependencies: MakeDbDependencies = oarsman_dependencies.get_make_db_dependencies()

    mak_db_output_data: MakeDbOutput = run_make_db(
        make_db_arguments,
        make_db_dependencies
    )

    print(mak_db_output_data.db_fpath)
# end def main
