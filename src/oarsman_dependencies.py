
import os


class OarsmanDependencies:

    def __init__(self):

        # Dependencies for kromsatel and co
        self.kromsatel_dirpath = None
        self.seqkit_fpath = None
    # end def __init__

    def get_make_amplicons_dependencies(self):

        return MakeAmpliconsDependencies(
            os.path.join(self.kromsatel_dirpath, 'db-scripts', 'make-amplicons.sh'),
            self.seqkit_fpath
        )

    # end def get_make_amplicons_args

# end class OarsmanArguments


class MakeAmpliconsDependencies:

    def __init__(
        self,
        make_amplicons_fpath= None,
        seqkit_fpath='seqkit'
    ):

        self.make_amplicons_fpath = make_amplicons_fpath
        self.seqkit_fpath = seqkit_fpath
    # end def __init__
# end class MakeAmpliconsDependencies
