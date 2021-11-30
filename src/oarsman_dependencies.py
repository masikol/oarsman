
import os


class OarsmanDependencies:

    def __init__(self):

        # Dependencies for kromsatel and co
        self.kromsatel_dirpath = None
        self.seqkit_fpath = None
        self.blastn_fpath = None
    # end def __init__

    def get_make_amplicons_dependencies(self):

        return MakeAmpliconsDependencies(
            os.path.join(self.kromsatel_dirpath, 'db-scripts', 'make-amplicons.sh'),
            self.seqkit_fpath
        )

    # end def get_make_amplicons_args

    def get_make_db_dependencies(self):

        return MakeDbDependencies(
            os.path.join(self.kromsatel_dirpath, 'db-scripts', 'make-db.sh'),
            self.blastn_fpath
        )

    # end def get_make_db_dependencies

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


class MakeDbDependencies:

    def __init__(
        self,
        make_db_path= None,
        blastn_fpath='blastn'
    ):

        self.make_db_path = make_db_path
        self.blastn_fpath = blastn_fpath
    # end def __init__
# end class MakeDbDependencies
