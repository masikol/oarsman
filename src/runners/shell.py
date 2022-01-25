
import subprocess as sp

from src.fatal_errors import FatalError


def launch_command(command_str, program_name_str, print_stdout=False):

    stdout = None if print_stdout else sp.PIPE
    pipe = sp.Popen(command_str, shell=True, stdout=stdout, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        stderr = stdout_stderr[1].decode('utf-8')
        error_msg = '\nError: program {} returned a non-zero exit code: {}' \
            '\nError message: {}' \
            .format(program_name_str, pipe.returncode, stderr)
        raise FatalError(error_msg)
    # end if
# end def
