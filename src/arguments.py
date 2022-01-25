
import os
import sys

import src.filesystem as fs
from src.reads import Reads
from src.input_modes import InputModes
from src.fatal_errors import FatalError


class OarsmanArguments:

    def __init__(self):

        # Input data
        # Reads
        self.reads_R1_fpaths = None
        self.reads_R2_fpaths = None
        self.reads_long_fpaths = None
        self.input_mode = None

        # Reference genome
        self.reference_fpath = None
        # Primers
        self.primers_fpath = None

        # Arguments for kromsatel
        self.kromsatel_args = ''

        # Variant calling
        self.min_variant_qual = 20 # in Phred scale

        # Annotation
        low_coverages_int = (10,)
        self.low_coverages = tuple(
            map(
                str,
                low_coverages_int
            )
        )

        # Misc
        self.outdir_path = os.path.join(os.getcwd(), 'oarsman_outdir')
        self.n_threads = 1 # thread

        self._set_output_subdirs()
        self._create_output_subdirs()
    # end def


    def get_number_of_samples(self):
        if self.input_mode in (InputModes.IlluminaSE, InputModes.IlluminaPE):
            return len(self.reads_R1_fpaths)
        elif self.input_mode == InputModes.Nanopore:
            return len(self.reads_long_fpaths)
        else:
            raise FatalError(
                '\nError: Internal error.\n' \
                'Please, contact the developer.\n' \
                'Error description: invalid mode in "get_number_of_samples".'
            )
        # end if
    # end def

    def _set_output_subdirs(self):
        self.output_subdirs = {
            'kromsatel_outdir': os.path.join(self.outdir_path, 'kromsatel_outdir'),
            'fasta_indicies':   os.path.join(self.outdir_path, 'fasta_indicies'),
            'read_mappings':    os.path.join(self.outdir_path, 'read_mappings'),
            'variant_calls':    os.path.join(self.outdir_path, 'variant_calls'),
            'consensus':        os.path.join(self.outdir_path, 'consensus'),
        }
    # end def


    def _create_output_subdirs(self):
        for subdir_path in self.output_subdirs.values():
            fs.create_directory(subdir_path)
        # end for
    # end def

    def get_kromsatel_args(self, i_sample):
        if self.reads_R1_fpaths is None:
            reads_R1_fpaths = None
        else:
            reads_R1_fpaths = self.reads_R1_fpaths[i_sample]
        # end if
        if self.reads_R2_fpaths is None:
            reads_R2_fpaths = None
        else:
            reads_R2_fpaths = self.reads_R2_fpaths[i_sample]
        # end if
        if self.reads_long_fpaths is None:
            reads_long_fpaths = None
        else:
            reads_long_fpaths = self.reads_long_fpaths[i_sample]
        # end if
        return KromsatelArguments(
            self,
            reads_R1_fpaths,
            reads_R2_fpaths,
            reads_long_fpaths,
            self.output_subdirs['kromsatel_outdir']
        )
    # end def

    def get_read_mapping_args(
        self,
        sample_name,
        reads,
        reference_seq_fpath,
        output_suffix
    ):

        return ReadMappingArguments(
            self,
            sample_name,
            reads,
            reference_seq_fpath,
            output_suffix
        )
    # end def

    def get_aln_postprocess_args(
        self,
        sample_name,
        mapping,
        reference_fpath,
        output_suffix
    ):
        return AlnPostprocessArguments(
            sample_name,
            mapping,
            reference_fpath,
            self.n_threads,
            output_suffix
        )
    # end def

    def get_var_call_args(self, sample_name, alignment_fpath):

        return CallVariantsArguments(
            sample_name,
            alignment_fpath,
            self.reference_fpath,
            self.min_variant_qual,
            self.output_subdirs['variant_calls'],
            self.output_subdirs['consensus'],
            self.n_threads
        )
    # end def

    def get_consens_annot_args(self, sample_name, seq_fpath, mapping):

        outfpath = os.path.join(
            self.output_subdirs['consensus'],
            f'{sample_name}_annotated_consensus.gbk'
        )

        return ConsensAnnotArguments(
            sample_name,
            seq_fpath,
            mapping,
            self.low_coverages,
            outfpath
        )
    # end def
# end class


class KromsatelArguments:

    def __init__(self,
                 oarsman_args,
                 reads_R1_fpaths,
                 reads_R2_fpaths,
                 reads_long_fpaths,
                 outdir_path):
        self.reads_R1_fpaths = reads_R1_fpaths
        self.reads_R2_fpaths = reads_R2_fpaths
        self.reads_long_fpaths = reads_long_fpaths
        self.input_mode = oarsman_args.input_mode
        self.primers_fpath = oarsman_args.primers_fpath
        self.reference_fpath = oarsman_args.reference_fpath
        self.outdir_path = outdir_path
        self.n_threads = oarsman_args.n_threads
        self.advanced_args = oarsman_args.kromsatel_args
    # end def
# end class


class ReadMappingArguments:

    def __init__(
        self,
        oarsman_args,
        sample_name,
        reads,
        reference_seq_fpath,
        output_suffix
    ):

        self.sample_name = sample_name
        self.reads = reads
        self.input_mode = oarsman_args.input_mode
        self.reference_fpath = reference_seq_fpath
        self.outdir_path = oarsman_args.output_subdirs['read_mappings']
        self.index_dir_path = oarsman_args.output_subdirs['fasta_indicies']
        self.output_suffix = output_suffix
        self.n_threads = oarsman_args.n_threads
    # end def
# end class


class AlnPostprocessArguments:

    def __init__(
        self,
        sample_name,
        raw_mapping,
        reference_fpath,
        n_threads,
        output_suffix
    ):

        self.sample_name = sample_name
        self.raw_mapping = raw_mapping
        self.reference_fpath = reference_fpath
        self.n_threads = n_threads
        self.output_suffix = output_suffix
    # end def
# end class


class CallVariantsArguments:

    def __init__(
        self,
        sample_name,
        alignment_fpath,
        reference_fpath,
        min_variant_qual,
        outdir_path,
        consensus_dirpath,
        n_threads
    ):

        self.sample_name = sample_name
        self.alignment_fpath = alignment_fpath
        self.reference_fpath = reference_fpath
        self.min_variant_qual = min_variant_qual
        self.outdir_path = outdir_path
        self.consensus_dirpath = consensus_dirpath
        self.n_threads = n_threads
    # end def
# end class


class ConsensAnnotArguments:

    def __init__(self, sample_name, seq_fpath, mapping, low_coverages, outfpath):
        self.sample_name = sample_name
        self.seq_fpath = seq_fpath
        self.mapping = mapping
        self.low_coverages = low_coverages
        self.outfpath = outfpath
    # end def
# end class