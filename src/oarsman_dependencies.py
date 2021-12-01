
import os


class OarsmanDependencies:

    def __init__(self):

        # Dependencies for kromsatel and co
        self.kromsatel_dirpath = None
        self.seqkit_fpath = None
        self.blastn_fpath = 'blastn' # it must be in PATH
        self.makeblastdb_fpath = 'makeblastdb' # it must be in PATH
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
            self.makeblastdb_fpath
        )
    # end def get_make_db_dependencies

    def get_kromsatel_dependencies(self):

        return KromsatelDependencies(
            os.path.join(self.kromsatel_dirpath, 'kromsatel.py'),
            self.blastn_fpath
        )
    # end def get_kromsatel_dependencies

    def get_pair_dependencies(self):

        return PairDependencies(
            self.seqkit_fpath
        )
    # end def get_pair_dependencies

# end class OarsmanArguments


class MakeAmpliconsDependencies:

    def __init__(
        self,
        make_amplicons_fpath,
        seqkit_fpath
    ):

        self.make_amplicons_fpath = make_amplicons_fpath
        self.seqkit_fpath = seqkit_fpath
    # end def __init__
# end class MakeAmpliconsDependencies


class MakeDbDependencies:

    def __init__(
        self,
        make_db_path,
        makeblastdb_fpath
    ):

        self.make_db_path = make_db_path
        self.makeblastdb_fpath = makeblastdb_fpath
    # end def __init__
# end class MakeDbDependencies


class KromsatelDependencies:

    def __init__(
        self,
        kromsatel_fpath,
        blastn_fpath
    ):

        self.kromsatel_fpath = kromsatel_fpath
        self.blastn_fpath = blastn_fpath
    # end def __init__
# end class KromsatelDependencies


class PairDependencies:

    def __init__(
        self,
        seqkit_fpath
    ):

        self.seqkit_fpath = seqkit_fpath
    # end def __init__
# end class PairDependencies
