
import src.pipeline
from src.fatal_errors import FatalError
from src.printing import print_err


def main():
    try:
        src.pipeline.run_pipeline()
    except FatalError as err:
        print_err(str(err))
        result_status = 1
    else:
        result_status = 0
    # end try

    return result_status
# end def
