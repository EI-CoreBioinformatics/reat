version 1.0

import "subworkflows/common/tasks.wdl"
import "subworkflows/common/structs.wdl"
import "subworkflows/align_short/wf_align_short.wdl" as aln_s
import "subworkflows/assembly_short/wf_assembly_short.wdl" as assm_s
import "subworkflows/portcullis/wf_portcullis.wdl" as portcullis_s
import "subworkflows/align_long/wf_align_long.wdl" as aln_l
import "subworkflows/assembly_long/wf_assembly_long.wdl" as assm_l
import "subworkflows/sanitise/wf_sanitise.wdl" as san

workflow wf_align {
    input {
        File reference_genome
        Array[PRSample]? paired_samples
        Array[LRSample]? LQ_long_read_samples
        Array[LRSample]? HQ_long_read_samples
        File? reference_annotation
        File? extra_junctions
        Map[String, Array[String]]? group_to_samples

        Float min_identity = 0.9
        Int? min_intron_len = 20
        Int? max_intron_len = 200000
        Int? max_intron_len_ends = 100000

        String portcullis_merge_operator = "max"
        String? portcullis_extra_parameters

        Boolean? skip_scallop = false

        String LQ_aligner = "minimap2"
        String HQ_aligner = "minimap2"
        String HQ_assembler = "merge"
        String LQ_assembler = "stringtie"

        # Long read aligners optional extra parameters
        String? HQ_aligner_extra_parameters
        String? LQ_aligner_extra_parameters

        # Short read aligners optional extra parameters
        String? PR_hisat_extra_parameters
        String? PR_star_extra_parameters

        # Assemblers optional extra parameters
        String? HQ_assembler_extra_parameters
        String? LQ_assembler_extra_parameters
        String? PR_stringtie_extra_parameters
        String? PR_scallop_extra_parameters

        Int? HQ_min_coverage
        Int? HQ_min_identity

        Int? LQ_min_coverage
        Int? LQ_min_identity

        RuntimeAttr? short_read_alignment_resources
        RuntimeAttr? short_read_alignment_sort_resources
        RuntimeAttr? short_read_stats_resources
        RuntimeAttr? short_read_merge_resources
        RuntimeAttr? short_read_stringtie_assembly_resources
        RuntimeAttr? short_read_scallop_assembly_resources
        RuntimeAttr? long_read_indexing_resources
        RuntimeAttr? long_read_alignment_resources
        RuntimeAttr? long_read_twopass_resources
        RuntimeAttr? long_read_twopass_merge_resources
        RuntimeAttr? long_read_assembly_resources
        RuntimeAttr? portcullis_resources
    }
    
    parameter_meta {
        reference_genome: "Reference genome to align against"
        paired_samples: "Paired short read samples, each item is defined by a biological replicate name with one or more technical"
        LQ_long_read_samples: "Low quality long read samples, each item is defined by a name, it's strand and one or more long read files."
        HQ_long_read_samples: "High quality long read samples, each item is defined by a name, it's strand and one or more long read files."
        reference_annotation: "Pre-existing annotation that will used during the alignment process and to cleanup the splicing sites."
    }

    call san.wf_sanitise {
        input:
        reference_genome = reference_genome,
        in_annotation = reference_annotation
    }

    if (defined(paired_samples)) {
        Array[PRSample] def_paired_samples = select_first([paired_samples])
        call aln_s.wf_align_short {
            input:
            samples = def_paired_samples,
            reference_genome = wf_sanitise.reference,
            reference_annotation = wf_sanitise.annotation,
            hisat_extra_parameters = PR_hisat_extra_parameters,
            star_extra_parameters = PR_star_extra_parameters,
            min_intron_len = select_first([min_intron_len, 20]),
            max_intron_len = select_first([max_intron_len, 200000]),
            alignment_resources = short_read_alignment_resources,
            sort_resources = short_read_alignment_sort_resources,
            merge_resources = short_read_merge_resources,
            stats_resources = short_read_stats_resources
        }

        call assm_s.wf_assembly_short {
            input:
            aligned_samples = wf_align_short.aligned_samples,
            reference_annotation = wf_sanitise.annotation,
            skip_scallop = skip_scallop,
            stringtie_extra_parameters = PR_stringtie_extra_parameters,
            scallop_extra_parameters = PR_scallop_extra_parameters,
            scallop_assembly_resources = short_read_scallop_assembly_resources,
            stringtie_assembly_resources = short_read_stringtie_assembly_resources
        }

        call portcullis_s.portcullis {
            input:
            reference = wf_sanitise.reference,
            annotation = wf_sanitise.annotation,
            merge_operator = portcullis_merge_operator,
            group_to_samples = group_to_samples,
            aligned_samples = wf_align_short.aligned_samples,
            portcullis_resources = portcullis_resources
        }
    }

    if (defined(LQ_long_read_samples)) {
        Array[LRSample] def_lq_long_sample = select_first([LQ_long_read_samples])
        call aln_l.wf_align_long as LQ_align {
            input:
            indexed_reference = wf_sanitise.indexed_reference,
            is_hq = false,
            aligner = LQ_aligner,
            long_samples = def_lq_long_sample,
            min_identity = min_identity,
            min_intron_len = select_first([min_intron_len, 20]),
            max_intron_len = select_first([max_intron_len, 200000]),
            max_intron_len_ends = select_first([max_intron_len_ends, 100000]),
            aligner_extra_parameters = select_first([LQ_aligner_extra_parameters, ""]),
            annotation = reference_annotation,
            extra_junctions = extra_junctions,
            portcullis_junctions = portcullis.pass_bed,
            indexing_resources = long_read_indexing_resources,
            alignment_resources = long_read_alignment_resources,
            twopass_resources = long_read_twopass_resources,
            twopass_merge_resources = long_read_twopass_merge_resources
        }

        # Check what is defined for lq-long and run that

        call assm_l.wf_assembly_long as LQ_assembly {
            input:
            reference_annotation = wf_sanitise.annotation,
            aligned_samples = LQ_align.bams,
            assembler = LQ_assembler,
            min_coverage = LQ_min_coverage,
            min_identity = LQ_min_identity,
            assembler_extra_parameters = LQ_assembler_extra_parameters,
            assembly_resources = long_read_assembly_resources,
            stats_output_prefix = "LQ_read_samples"
        }
    }

    if (defined(HQ_long_read_samples)) {
        Array[LRSample] def_hq_long_sample = select_first([HQ_long_read_samples])
        call aln_l.wf_align_long as HQ_align {
            input:
            indexed_reference = wf_sanitise.indexed_reference,
            is_hq = true,
            aligner = HQ_aligner,
            long_samples = def_hq_long_sample,
            min_identity = min_identity,
            min_intron_len = min_intron_len,
            max_intron_len = max_intron_len,
            max_intron_len_ends = max_intron_len_ends,
            aligner_extra_parameters = select_first([HQ_aligner_extra_parameters, ""]),
            annotation = reference_annotation,
            extra_junctions = extra_junctions,
            portcullis_junctions = portcullis.pass_bed,
            indexing_resources = long_read_indexing_resources,
            alignment_resources = long_read_alignment_resources,
            twopass_resources = long_read_twopass_resources,
            twopass_merge_resources = long_read_twopass_merge_resources
        }

        # Check what is defined for hq-long and run that

        call assm_l.wf_assembly_long as HQ_assembly {
            input:
            reference_annotation = wf_sanitise.annotation,
            aligned_samples = HQ_align.bams,
            assembler = HQ_assembler,
            min_coverage = HQ_min_coverage,
            min_identity = HQ_min_identity,
            assembler_extra_parameters = HQ_assembler_extra_parameters,
            assembly_resources = long_read_assembly_resources,
            stats_output_prefix = "HQ_read_samples"
        }
    }

    output {
        IndexedReference clean_reference_index = wf_sanitise.indexed_reference
        File? clean_annotation = wf_sanitise.annotation

        File? pass_filtered_tab = portcullis.pass_tab
        File? pass_filtered_bed = portcullis.pass_bed
        File? pass_filtered_gff3 = portcullis.pass_gff3

        File? fail_filtered_tab = portcullis.fail_tab
        File? fail_filtered_bed = portcullis.fail_bed
        File? fail_filtered_gff3 = portcullis.fail_gff3

        Array[AlignedSample]? SR_bams = wf_align_short.aligned_samples
        Array[AlignedSample]? LQ_bams = LQ_align.bams
        Array[AlignedSample]? HQ_bams = HQ_align.bams

        Array[AssembledSample]? SR_gff = wf_assembly_short.assemblies
        Array[AssembledSample]? LQ_gff = LQ_assembly.gff
        Array[AssembledSample]? HQ_gff = HQ_assembly.gff

#        Array[Array[Array[File]]]? plots = wf_align_short.plots
        Array[Array[File]]? stats = wf_align_short.stats
        Array[Array[File]]? actg_cycles_plots = wf_align_short.actg_cycles_plots
        Array[Array[File]]? coverage_plots = wf_align_short.coverage_plots
        Array[Array[File]]? gc_content_plots = wf_align_short.gc_content_plots
        Array[Array[File]]? gc_depth_plots = wf_align_short.gc_depth_plots
        Array[Array[File]]? htmls = wf_align_short.htmls
        Array[Array[File]]? indel_cycles_plots = wf_align_short.indel_cycles_plots
        Array[Array[File]]? indel_dist_plots = wf_align_short.indel_dist_plots
        Array[Array[File]]? insert_size_plots = wf_align_short.insert_size_plots
        Array[Array[File]]? quals_plots = wf_align_short.quals_plots
        Array[Array[File]]? quals2_plots = wf_align_short.quals2_plots
        Array[Array[File]]? quals3_plots = wf_align_short.quals3_plots
        Array[Array[File]]? quals_hm_plots = wf_align_short.quals_hm_plots

        Array[File]? SR_summary_stats = wf_align_short.summary_stats
        Array[File]? LQ_summary_stats = LQ_align.summary_stats
        Array[File]? HQ_summary_stats = HQ_align.summary_stats

        File? SR_summary_stats_table = wf_align_short.summary_stats_table
        File? LQ_summary_stats_table = LQ_align.summary_stats_table
        File? HQ_summary_stats_table = HQ_align.summary_stats_table

        Array[File]? SR_assembly_stats = wf_assembly_short.stats
        Array[File]? LQ_assembly_stats = LQ_assembly.stats
        Array[File]? HQ_assembly_stats = HQ_assembly.stats

        File? SR_stringtie_summary_stats = wf_assembly_short.stringtie_summary_stats
        File? SR_scallop_summary_stats = wf_assembly_short.scallop_summary_stats
        File? LQ_assembly_summary_stats = LQ_assembly.summary_stats
        File? HQ_assembly_summary_stats = HQ_assembly.summary_stats
    }
}