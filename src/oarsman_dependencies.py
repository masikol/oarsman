
import os


class OarsmanDependencies:

    def __init__(self):

        # Dependencies for kromsatel and co
        self.kromsatel_dirpath = None
        self.seqkit_fpath = None
        self.blastn_fpath = 'blastn' # it must be in PATH
        self.makeblastdb_fpath = 'makeblastdb' # it must be in PATH

        # Dependencies for read mapping
        self.bowtie2_fpath = None
        self.bwa_fpath = None
        self.samtools_fpath = None

        # Dependencies for variant calling
        self.bcftools_fpath = None

        # Depencencies for consensus annnotation
        self.highlighter_fpath = None
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

    def get_bwa_dependencies(self):

        return BwaDependencies(
            self.bwa_fpath
        )
    # end def get_bwa_dependencies

    def get_bowtie2_dependencies(self):

        return Bowtie2Dependencies(
            self.bowtie2_fpath,
            self.bowtie2_fpath + '-build'
        )
    # end def get_bowtie2_dependencies

    def get_samtools_dependencies(self):
        return SamtoolsDependencies(
            self.samtools_fpath
        )
    # end def get_samtools_dependencies

    def get_bcftools_dependencies(self):
        return BcfVarCallDependencies(self.bcftools_fpath)
    # end def get_bcftools_dependencies

    def get_highlighter_dependencies(self):
        return HighlighterDependencies(self.highlighter_fpath, self.samtools_fpath)
    # end def get_highlighter_dependencies
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


class BwaDependencies:

    def __init__(
        self,
        bwa_fpath
    ):

        self.bwa_fpath = bwa_fpath
# end class BwaDependencies


class Bowtie2Dependencies:

    def __init__(
        self,
        bowtie2_fpath,
        bowtie2_build_fpath
    ):

        self.bowtie2_fpath = bowtie2_fpath
        self.bowtie2_build_fpath = bowtie2_build_fpath
# end class Bowtie2Dependencies


class SamtoolsDependencies:
    def __init__(self, samtools_fpath):
        self.samtools_fpath = samtools_fpath
    # end def __init__
# end class SamtoolsDependencies


class BcfVarCallDependencies:

    def __init__(self, bcftools_fpath):
        self.bcftools_fpath = bcftools_fpath
    # end def __init__
# end class BcfVarCallDependencies


class HighlighterDependencies:

    def __init__(self, highlighter_fpath, samtools_fpath):
        self.highlighter_fpath = highlighter_fpath
        self.samtools_fpath = samtools_fpath
    # end def __init__
# end class HighlighterDependencies
