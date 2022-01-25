#!/usr/bin/env python3

__version__ = '0.2.a'
__last_update_date__ = '2022-01-25'
# __author__ = 'Maxim Sikolenko'

# |===== Check python interpreter version. =====|

import sys

if sys.version_info.major + sys.version_info.minor * 0.1 < (3.5 - 1e-9):
    print( '\nYour Python interpreter version is ' + '%d.%d' % (sys.version_info.major,
        sys.version_info.minor) )
    print('   Please, use Python 3.5 or newer.\a')
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    if sys.platform.startswith('win'):
        if sys.version_info.minor == 2:
            raw_input('Press ENTER to exit:')
        else:
            input('Press ENTER to exit:')
        # end if
    # end if
    sys.exit(1)
# end if


import os

print(f'\n  |=== {os.path.basename(__file__)} v{__version__} ===|\n')

from src.main import main

if __name__ == '__main__':
    main()
# end if
