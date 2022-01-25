
import os


class OarsmanDependencies:

    def __init__(self):

        # Dependencies for kromsatel and co
        self.kromsatel_fpath = None
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

    def get_kromsatel_dependencies(self):

        return KromsatelDependencies(
            self.kromsatel_fpath,
            self.blastn_fpath,
            self.makeblastdb_fpath
        )
    # end def get_kromsatel_dependencies

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


class KromsatelDependencies:

    def __init__(self,
                 kromsatel_fpath,
                 blastn_fpath,
                 makeblastdb_fpath):

        self.kromsatel_fpath = kromsatel_fpath
        self.blastn_fpath = blastn_fpath
        self.makeblastdb_fpath = makeblastdb_fpath
    # end def __init__
# end class


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
