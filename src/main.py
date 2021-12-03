
import sys

from src.parse_arguments import parse_arguments
from src.seq_sample import fastq_fpath_to_sample_name

from src.runners.make_amplicons_runner import run_make_amplicons
from src.oarsman_arguments import MakeAmpliconsArguments
from src.oarsman_dependencies import MakeAmpliconsDependencies
from src.output_data import MakeAmpliconsOutput

from src.runners.make_db_runner import run_make_db
from src.oarsman_arguments import MakeDbArguments
from src.oarsman_dependencies import MakeDbDependencies
from src.output_data import MakeDbOutput

from src.runners.kromsatel_runner import run_kromsatel
from src.oarsman_arguments import KromsatelArguments
from src.oarsman_dependencies import KromsatelDependencies
from src.output_data import KromsatelOutput

from src.runners.pair_runner import run_pair
from src.oarsman_arguments import PairArguments
from src.oarsman_dependencies import PairDependencies
from src.output_data import PairOutput


from src.oarsman_arguments import ReadMappingArguments
from src.output_data import ReadMappingOutput

from src.runners.bowtie2_runner import run_bowtie2
from src.oarsman_dependencies import Bowtie2Dependencies

from src.runners.bwa_runner import run_bwa
from src.oarsman_dependencies import BwaDependencies


from src.runners.aln_preprocess_runner import run_aln_preprocess
from src.oarsman_arguments import AlnPreprocessArguments
from src.oarsman_dependencies import SamtoolsDependencies
from src.output_data import AlnPreprocessOutput

from src.runners.bcftools_var_call_runner import run_bcftools_var_call
from src.oarsman_arguments import CallVariantsArguments
from src.oarsman_dependencies import BcfVarCallDependencies
from src.output_data import VariantCallingOutput



def main():

    # Parse command line arguments
    oarsman_args, oarsman_dependencies = parse_arguments()

    print(oarsman_args.reads_R1_fpaths)
    print(oarsman_args.reads_R2_fpaths)
    # sys.exit(0)

    # Make amplicons (ma) for kromsatel
    ma_arguments: MakeAmpliconsArguments = oarsman_args.get_make_amplicons_args()
    ma_dependencies: MakeAmpliconsDependencies = oarsman_dependencies.get_make_amplicons_dependencies()

    ma_output_data: MakeAmpliconsOutput = run_make_amplicons(
        ma_arguments,
        ma_dependencies
    )

    print(ma_output_data.no_primers_fpath)
    print(ma_output_data.with_primers_fpath)


    # Make-db for kromsatel
    make_db_arguments: MakeDbArguments = oarsman_args.get_make_db_args(
        ma_output_data.no_primers_fpath
    )
    make_db_dependencies: MakeDbDependencies = oarsman_dependencies.get_make_db_dependencies()

    make_db_output_data: MakeDbOutput = run_make_db(
        make_db_arguments,
        make_db_dependencies
    )

    print(make_db_output_data.db_fpath)


    # From this moment we will iterate over samples

    n_samples = len(oarsman_args.reads_R1_fpaths)
    for i_sample in range(n_samples):

        sample_name = fastq_fpath_to_sample_name(
            oarsman_args.reads_R1_fpaths[i_sample][0]
        )

        print(f'  \n|=== Processing sample #{i_sample+1}/{n_samples}: {sample_name} ===|\n')

        # Run kromsatel
        kromsatel_arguments: KromsatelArguments = oarsman_args.get_kromsatel_args(
            i_sample,
            make_db_output_data.db_fpath
        )
        kromsatel_dependencies: KromsatelDependencies = oarsman_dependencies.get_kromsatel_dependencies()

        kromsatel_output: KromsatelOutput = run_kromsatel(
            kromsatel_arguments,
            kromsatel_dependencies
        )

        # Match up paired reads
        if not kromsatel_output.reads_R2_fpath is None:
            # After kromsatel paired reads are shuffled and now have to be matched up together

            pair_arguments: PairArguments = PairArguments(
                kromsatel_output.reads_R1_fpath,
                kromsatel_output.reads_R2_fpath
            )
            pair_dependencies: PairDependencies = oarsman_dependencies.get_pair_dependencies()

            pair_output: PairOutput = run_pair(pair_arguments, pair_dependencies)

            print(
                pair_output.reads_R1_paired_fpath,
                pair_output.reads_R2_paired_fpath,
                pair_output.unpaired_reads_fpaths
            )
        else:
            # If we have unpaired data, "forward" reads will be processed as unpaired from now
            pair_output: PairOutput = PairOutput(
                None,
                None,
                [kromsatel_output.reads_R1_fpath]
            )
        # end if


        # Map the reads
        read_map_args = oarsman_args.get_read_mapping_args(
            sample_name,
            pair_output.reads_R1_paired_fpath,
            pair_output.reads_R2_paired_fpath,
            pair_output.unpaired_reads_fpaths
        )
        # read_map_dependencies = oarsman_dependencies.get_bowtie2_dependencies()

        # read_mapping_output: ReadMappingOutput = run_bowtie2(
        #     read_map_args,
        #     read_map_dependencies
        # )

        read_map_dependencies = oarsman_dependencies.get_bwa_dependencies()

        read_mapping_output: ReadMappingOutput = run_bwa(
            read_map_args,
            read_map_dependencies
        )


        # Preprocess alignment data
        aln_preproc_args = oarsman_args.get_aln_preprocess_args(
            sample_name,
            read_mapping_output.alignment_fpath
        )
        aln_preproc_dependencies = oarsman_dependencies.get_samtools_dependencies()

        aln_preproces_output = run_aln_preprocess(
            aln_preproc_args,
            aln_preproc_dependencies
        )


        # Call variants
        var_call_args: CallVariantsArguments = oarsman_args.get_var_call_args(
            sample_name,
            aln_preproces_output.preproc_aln_fpath
        )
        var_call_dependencies: BcfVarCallDependencies = oarsman_dependencies.get_bcftools_dependencies()

        var_call_output: VariantCalling_Output = run_bcftools_var_call(
            var_call_args,
            var_call_dependencies
        )

        print(f'Consensus: `{var_call_output.consensus_fpath}`')
    # end for
# end def main
