
import src.versions


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
        self.lofreq_fpath = None

        # Depencencies for consensus annnotation
        self.highlighter_fpath = None
    # end def

    def check_versions(self):
        src.versions.check_kromsatel_version(self.kromsatel_fpath)
        src.versions.check_highlighter_version(self.highlighter_fpath)
        src.versions.check_samtools_version(self.samtools_fpath)
    # end def

    def get_kromsatel_dependencies(self):
        return KromsatelDependencies(
            self.kromsatel_fpath,
            self.blastn_fpath,
            self.makeblastdb_fpath
        )
    # end def

    def get_bwa_dependencies(self):
        return BwaDependencies(
            self.bwa_fpath,
            self.samtools_fpath
        )
    # end def

    def get_bowtie2_dependencies(self):
        return Bowtie2Dependencies(
            self.bowtie2_fpath,
            self.bowtie2_fpath + '-build'
        )
    # end def

    def get_samtools_dependencies(self):
        return SamtoolsDependencies(
            self.samtools_fpath
        )
    # end def

    def get_bcftools_dependencies(self):
        return BcfVarCallDependencies(self.bcftools_fpath)
    # end def

    def get_lofreq_dependencies(self):
        return LofreqVarCallDependencies(
            self.lofreq_fpath, self.samtools_fpath, self.bcftools_fpath
        )
    # end def

    def get_highlighter_dependencies(self):
        return HighlighterDependencies(self.highlighter_fpath, self.samtools_fpath)
    # end def
# end class


class KromsatelDependencies:

    def __init__(self,
                 kromsatel_fpath,
                 blastn_fpath,
                 makeblastdb_fpath):
        self.kromsatel_fpath = kromsatel_fpath
        self.blastn_fpath = blastn_fpath
        self.makeblastdb_fpath = makeblastdb_fpath
    # end def
# end class


class BwaDependencies:

    def __init__(self, bwa_fpath, samtools_fpath):
        self.bwa_fpath = bwa_fpath
        self.samtools_fpath = samtools_fpath
    # end def
# end class


class Bowtie2Dependencies:

    def __init__(self,
                 bowtie2_fpath,
                 bowtie2_build_fpath):
        self.bowtie2_fpath = bowtie2_fpath
        self.bowtie2_build_fpath = bowtie2_build_fpath
# end class


class SamtoolsDependencies:

    def __init__(self, samtools_fpath):
        self.samtools_fpath = samtools_fpath
    # end def
# end class


class BcfVarCallDependencies:

    def __init__(self, bcftools_fpath):
        self.bcftools_fpath = bcftools_fpath
    # end def
# end class


class LofreqVarCallDependencies:

    def __init__(self, lofreq_fpath, samtools_fpath, bcftools_fpath):
        self.lofreq_fpath = lofreq_fpath
        self.samtools_fpath = samtools_fpath
        self.bcftools_fpath = bcftools_fpath
    # end def
# end class


class HighlighterDependencies:

    def __init__(self, highlighter_fpath, samtools_fpath):
        self.highlighter_fpath = highlighter_fpath
        self.samtools_fpath = samtools_fpath
    # end def
# end class
