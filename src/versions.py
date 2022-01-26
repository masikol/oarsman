"""
This module is for checking versions of dependencies.
"""

import os
import re

from src.fatal_errors import FatalError
from src.runners.shell import launch_command_get_stdout


def check_kromsatel_version(kromsatel_fpath):

    min_major_version = 2

    kromsatel_version_cmd = _configure_pyscript_version_cmd(kromsatel_fpath)
    stdout_str = launch_command_get_stdout(
        kromsatel_version_cmd,
        'kromsatel.py'
    )

    try:
        version_str = _parse_ternary_version(stdout_str)
        major_version = _parse_major_version(version_str)
    except (ValueError, AttributeError):
        error_msg = '\nError: cannot parse the version of kromsatel' \
            ' or it has an invalid pattern.\n' \
            'Please, contact the develiper.'
        raise FatalError(error_msg)
    # end try

    if major_version < min_major_version:
        error_msg = '\nError: invalid version of kromsatel: {}\n' \
            'Minimum version is 2.0.a'.format(version_str)
        raise FatalError(error_msg)
    # end if

    print('Kromsatel version: {}. Ok.'.format(version_str))
# end def


def check_highlighter_version(highlighter_fpath):

    min_major_version = 2
    min_minor_version = 1

    highlighter_version_cmd = _configure_pyscript_version_cmd(highlighter_fpath)
    stdout_str = launch_command_get_stdout(
        highlighter_version_cmd,
        'consensus-highlighter.py'
    )

    try:
        version_str = _parse_ternary_version(stdout_str)
        major_version = _parse_major_version(version_str)
        minor_version = _parse_minor_version(version_str)
    except (ValueError, AttributeError):
        error_msg = '\nError: cannot parse the version of consensus-highlighter' \
            ' or it has an invalid pattern.\n' \
            'Please, contact the develiper.'
        raise FatalError(error_msg)
    # end try

    if major_version < min_major_version:
        error_msg = '\nError: invalid version of consensus-highlighter: {}\n' \
            'Minimum version is 2.1.a'.format(version_str)
        raise FatalError(error_msg)
    # end if

    if minor_version < min_minor_version:
        error_msg = '\nError: invalid version of consensus-highlighter: {}\n' \
            'Minimum version is 2.1.a'.format(version_str)
        raise FatalError(error_msg)
    # end if

    print('consensus-highlighter version: {}. Ok.'.format(version_str))
# end def


def check_samtools_version(samtools_fpath):

    min_major_version = 1
    min_minor_version = 12

    samtools_version_cmd = _configure_samtools_version_cmd(samtools_fpath)
    stdout_str = launch_command_get_stdout(
        samtools_version_cmd,
        'samtools --version'
    )

    try:
        version_str = _parse_samtools_version(stdout_str)
        major_version = _parse_major_version(version_str)
        minor_version = _parse_minor_version(version_str)
    except (ValueError, AttributeError):
        error_msg = '\nError: cannot parse the version of samtools' \
            ' or it has an invalid pattern.\n' \
            'Please, contact the develiper.'
        raise FatalError(error_msg)
    # end try

    if major_version < min_major_version:
        error_msg = '\nError: invalid version of samtools: {}\n' \
            'Minimum version is {}.{}' \
            .format(version_str, min_major_version, min_minor_version)
        raise FatalError(error_msg)
    # end if

    if minor_version < min_minor_version:
        error_msg = '\nError: invalid version of samtools: {}\n' \
            'Minimum version is {}.{}' \
            .format(version_str, min_major_version, min_minor_version)
        raise FatalError(error_msg)
    # end if

    print('Samtools version: {}. Ok.'.format(version_str))
# end def


def _configure_pyscript_version_cmd(script_fpath):

    if os.path.exists(script_fpath):
        script_cmd_part = 'python3 {}'.format(script_fpath)
    else:
        # The script executable is in PATH
        script_cmd_part = script_fpath
    # end if

    command = ' '.join(
        [
            script_cmd_part,
            '--version'
        ]
    )

    return command
# end def


def _configure_samtools_version_cmd(samtools_fpath):
    command = ' '.join(
        [
            samtools_fpath,
            '--version',
        ]
    )
    return command
# end def


def _parse_ternary_version(stdout_str):
    version_pattern = r'([0-9]+\.[0-9]+\.[a-z]+)'
    version_str = re.search(version_pattern, stdout_str).group(1)
    return version_str
# end def


def _parse_samtools_version(samtools_stdout):
    version_pattern = r'samtools ([0-9]+\.[0-9]+)'
    version_str = re.search(version_pattern, samtools_stdout).group(1)
    return version_str
# end def


def _parse_major_version(version_str):
    major_version_str = version_str.partition('.')[0]
    return int(major_version_str)
# end def


def _parse_minor_version(version_str):
    minor_version_pattern = r'[0-9]+\.([0-9]+).*'
    minor_version_str = re.search(minor_version_pattern, version_str).group(1)
    return int(minor_version_str)
# end def
