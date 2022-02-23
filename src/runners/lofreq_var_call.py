
import os
import re
import gzip
from functools import partial

from src.shell import launch_command
from src.filesystem import rm_if_exists
from src.fatal_errors import FatalError

from src.arguments import CallVariantsArguments
from src.dependencies import BcfVarCallDependencies
from src.data_transfer_objects import VariantCall, ReadMapping, VariantCallIndex


def run_lofreq_var_call(args, dependencies):

    print('Performing Viterbi realignment, which corrects mapping errors...')
    realigned_alignment = _viterbi_realignment(args, dependencies)
    realigned_alignment.check_existance()

    print('Inserting indel qualities...')
    indelqualed_alignment = _insert_indel_quals(args, dependencies, realigned_alignment)
    indelqualed_alignment.check_existance()
    realigned_alignment.cleanup()

    print('Sorting realigned BAM file...')
    indelqualed_sorted_alignment = _sort_bam(args, dependencies, indelqualed_alignment)
    indelqualed_sorted_alignment.check_existance()
    indelqualed_alignment.cleanup()

    print('Calling variants with lofreq...')
    var_call = _call_variants(args, dependencies, indelqualed_sorted_alignment)
    var_call.check_existance()
    indelqualed_sorted_alignment.cleanup()

    print('Filtering variants showing >50% frequency...')
    filtered_var_call = _filter_variants(args, dependencies, var_call)
    filtered_var_call.check_existance()

    _amend_malformatted_vcf_index(args, filtered_var_call)

    print('Converting VCF -> BCF')
    bcf_filtered_var_call = _vcf2bcf(args, dependencies, filtered_var_call)
    bcf_filtered_var_call.check_existance()
    filtered_var_call.cleanup()

    print('Indexing filtered BCF file...')
    filt_bcf_index = _index_var_call(bcf_filtered_var_call, dependencies)
    filt_bcf_index.check_existance()

    return bcf_filtered_var_call
# end def


def _viterbi_realignment(var_call_args, dependencies):
    realigned_bam_fpath = _configure_realigned_fpath(var_call_args)
    command_str = _configure_viterbi_command(
        var_call_args,
        dependencies,
        realigned_bam_fpath
    )

    rm_if_exists(realigned_bam_fpath)
    launch_command(command_str, 'lofreq viterbi')
    return ReadMapping(realigned_bam_fpath)
# end def


def _configure_realigned_fpath(var_call_args):
    var_call_fpath = os.path.join(
        var_call_args.outdir_path,
        f'{var_call_args.sample_name}.realigned.bam'
    )
    return var_call_fpath
# end def


def _configure_viterbi_command(var_call_args, dependencies, outfpath):
    command = ' '.join(
        [
            dependencies.lofreq_fpath, 'viterbi',
            '-f {}'.format(var_call_args.reference_fpath),
            '-o {}'.format(outfpath),
            var_call_args.alignment_fpath,
        ]
    )
    return command
# end def


def _insert_indel_quals(var_call_args, dependencies, realigned_alignment):
    alignment_indelqualed_fpath = _configure_alignment_indelqualed_fpath(var_call_args)
    command_str = _configure_indelqual_command(
        var_call_args,
        dependencies,
        realigned_alignment,
        alignment_indelqualed_fpath
    )

    rm_if_exists(alignment_indelqualed_fpath)
    launch_command(command_str, 'lofreq indelqual')
    return ReadMapping(alignment_indelqualed_fpath)
# end def


def _configure_alignment_indelqualed_fpath(var_call_args):
    var_call_fpath = os.path.join(
        var_call_args.outdir_path,
        f'{var_call_args.sample_name}.indelqualed.bam'
    )
    return var_call_fpath
# end def


def _configure_indelqual_command(var_call_args,
                                 dependencies,
                                 realigned_alignment,
                                 outfpath):
    command = ' '.join(
        [
            dependencies.lofreq_fpath, 'indelqual',
            '--dindel',
            '-f {}'.format(var_call_args.reference_fpath),
            '-o {}'.format(outfpath),
            realigned_alignment.alignment_fpath,
        ]
    )
    return command
# end def


def _sort_bam(var_call_args, dependencies, read_mapping):
    sorted_mapping_fpath = _configure_sorted_bam_fpath(read_mapping)
    command_str = _configure_sort_bam_command(
        var_call_args,
        dependencies,
        read_mapping,
        sorted_mapping_fpath
    )

    launch_command(command_str, 'samtools sort')
    return ReadMapping(sorted_mapping_fpath)
# end def


def _configure_sorted_bam_fpath(read_mapping):

    if read_mapping.alignment_fpath.endswith('.bam'):
        return read_mapping.alignment_fpath[:-4] + '.sorted.bam'
    else:
        raise FatalError('\nInternal error 3. Please, contact the developer.')
    # end if
# end def


def _configure_sort_bam_command(var_call_args,
                                dependencies,
                                read_mapping,
                                sorted_mapping_fpath):

    command_str = ' '.join(
        [
            dependencies.samtools_fpath, 'sort',
            '-O BAM',
            '-@ {}'.format(var_call_args.n_threads),
            '-o {}'.format(sorted_mapping_fpath),
            read_mapping.alignment_fpath
        ]
    )
    
    return command_str
# end def


def _call_variants(var_call_args, dependencies, indelqualed_alignment):
    var_call_fpath = _configure_var_call_fpath(var_call_args)
    command_str = _configure_variant_call_command(
        var_call_args,
        dependencies,
        indelqualed_alignment,
        var_call_fpath
    )

    rm_if_exists(var_call_fpath)
    launch_command(command_str, 'lofreq call')
    return VariantCall(var_call_fpath)
# end def


def _configure_var_call_fpath(var_call_args):
    var_call_fpath = os.path.join(
        var_call_args.outdir_path,
        f'{var_call_args.sample_name}.vcf'
    )
    return var_call_fpath
# end def


def _configure_variant_call_command(var_call_args,
                                    dependencies,
                                    indelqualed_alignment,
                                    var_call_fpath):

    min_base_qual = 20
    min_alt_base_qual = 20

    command = ' '.join(
        [
            dependencies.lofreq_fpath, 'call',
            '-f {}'.format(var_call_args.reference_fpath),
            '-o {}'.format(var_call_fpath),
            '-q {} -Q {}'.format(min_base_qual, min_alt_base_qual),
            '--call-indels', '--no-default-filter',
            indelqualed_alignment.alignment_fpath,
        ]
    )

    return command
# end def


def _filter_variants(var_call_args, dependencies, var_call):
    filtered_var_call_fpath = _configure_filt_var_call_fpath(var_call_args)
    command_str = _configure_filter_variants_command(
        var_call_args,
        dependencies,
        var_call,
        filtered_var_call_fpath
    )

    rm_if_exists(filtered_var_call_fpath)
    launch_command(command_str, 'lofreq call')
    return VariantCall(filtered_var_call_fpath)
# end def


def _configure_filt_var_call_fpath(var_call_args):
    filt_var_call_fpath = os.path.join(
        var_call_args.outdir_path,
        f'{var_call_args.sample_name}.filtered.vcf'
    )
    return filt_var_call_fpath
# end def


def _configure_filter_variants_command(var_call_args,
                                       dependencies,
                                       var_call,
                                       filtered_var_call_fpath):

    min_freq = var_call_args.freq_threshold + 1e-6

    command = ' '.join(
        [
            dependencies.lofreq_fpath, 'filter',
            '-i {}'.format(var_call.var_call_fpath),
            '-o {}'.format(filtered_var_call_fpath),
            '--af-min {}'.format(min_freq),
            '--no-defaults'
        ]
    )

    return command
# end def


def _amend_malformatted_vcf_index(var_call_args, var_call):
    seq_id, seq_length = _get_fasta_seqid_and_len(var_call_args.reference_fpath)
    _insert_contig_header_field(var_call, seq_id, seq_length)
# end def


def _get_fasta_seqid_and_len(fasta_fpath):

    if fasta_fpath.endswith('.gz'):
        open_func = partial(gzip.open, mode='rt')
    else:
        open_func = partial(open, mode='rt')
    # end if

    seq_length = 0

    with open_func(fasta_fpath) as fasta_infile:
        header_line = fasta_infile.readline().strip()
        seq_id = header_line[1:].partition(' ')[0]

        line = 'some_string'
        while not line.startswith('>') and not line == '':
            line = fasta_infile.readline().strip()
            seq_length += len(line)
        # end with
    # end with

    return seq_id, seq_length
# end def


def _insert_contig_header_field(var_call, seq_id, seq_length):

    contig_field = '##contig=<ID={},length={}>\n'.format(seq_id, seq_length)

    with open(var_call.var_call_fpath, 'rt') as vc_in_file:
        vc_lines = vc_in_file.readlines()
    # end with

    reference_pattern = r'##reference\='
    with open(var_call.var_call_fpath, 'wt') as vc_out_file:
        for line in vc_lines:
            vc_out_file.write(line)
            if not re.search(reference_pattern, line) is None:
                vc_out_file.write(contig_field)
            # end if
        # end with
    # end with
# end def


def _vcf2bcf(var_call_args, dependencies, vcf_var_call):
    bcf_var_call_fpath = _configure_bcf_fpath(vcf_var_call)
    command_str = _configure_vcf2bcf_command(
        var_call_args,
        dependencies,
        vcf_var_call,
        bcf_var_call_fpath
    )

    launch_command(command_str, 'bcftools view')
    return VariantCall(bcf_var_call_fpath)
# end def


def _configure_bcf_fpath(vcf_var_call):
    if vcf_var_call.var_call_fpath.endswith('.vcf'):
        return vcf_var_call.var_call_fpath[:-4] + '.bcf'
    else:
        raise FatalError('\nInternal error 2. Please, contact the developer.')
    # end if
# end def


def _configure_vcf2bcf_command(var_call_args,
                               dependencies,
                               vcf_var_call,
                               bcf_var_call_fpath):
    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'view',
            '-Ob',
            '-o {}'.format(bcf_var_call_fpath),
            vcf_var_call.var_call_fpath
        ]
    )
    return command
# end def


def _index_var_call(var_call, dependencies):
    command_str = _configure_index_var_call_command(
        var_call,
        dependencies
    )
    launch_command(command_str, 'bcftools index')

    return VariantCallIndex(var_call.var_call_fpath)
# end def


def _configure_index_var_call_command(var_call, dependencies):

    command = ' '.join(
        [
            dependencies.bcftools_fpath, 'index',
            var_call.var_call_fpath
        ]
    )

    return command
# end def
