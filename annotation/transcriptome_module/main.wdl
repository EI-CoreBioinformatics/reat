version 1.0

import "mikado.wdl" as wfm
import "align.wdl" as waln
import "subworkflows/common/structs.wdl"
# import "subworkflows/exonerate/wf_exonerate.wdl" as exonerate
# import "subworkflows/repeat_masker/wf_repeat_masker.wdl" as repeatmasker
import "subworkflows/sanitise/wf_sanitise.wdl" as san

workflow ei_annotation {
    input {
        File reference_genome
        Array[PRSample]? paired_samples
        Array[LRSample]? LQ_long_read_samples
        Array[LRSample]? HQ_long_read_samples
        File? annotation
        File mikado_scoring_file
        File orf_calling_proteins
        File homology_proteins
    }

    parameter_meta {
        reference_genome: "Reference genome to align against"
        paired_samples: "Paired short read samples, each item is defined by a biological replicate name with one or more technical replicates. Technical replicates are defined by a name, R1, R2 and strand."
        LQ_long_read_samples: "Low quality long read samples, each item is defined by a name, it's strand and one or more long read files."
        HQ_long_read_samples: "High quality long read samples, each item is defined by a name, it's strand and one or more long read files."
        annotation: "Pre-existing annotation that will used during the alignment process and to cleanup the splicing sites."
        mikado_scoring_file: "Mikado scoring file"
    }

    call san.wf_sanitise {
        input:
        reference_genome = reference_genome,
        in_annotation = annotation
    }

    call waln.wf_align {
        input:
        reference_genome = reference_genome,
        paired_samples = paired_samples,
        LQ_long_read_samples = LQ_long_read_samples,
        HQ_long_read_samples = HQ_long_read_samples,
        reference_annotation = wf_sanitise.annotation
    }

    call wfm.wf_main_mikado {
        input:
        mikado_scoring_file = mikado_scoring_file,
        reference_genome = wf_align.clean_reference_index,
        SR_assemblies = wf_align.SR_gff,
        LQ_assemblies = wf_align.LQ_gff,
        HQ_assemblies = wf_align.HQ_gff,
        orf_calling_proteins = orf_calling_proteins,
        homology_proteins = homology_proteins
    }

#    call repeatmasker.wf_repeat_masker as RepeatMasker {
#        input:
#        reference_fasta = wf_sanitize.indexed_reference.fasta
#    }
#
#    if (defined(protein_related_species)) {
#        Array[LabeledFasta] def_protein_related_species = select_first([protein_related_species])
#        call exonerate.wf_exonerate as Exonerate {
#            input:
#            related_species_protein = def_protein_related_species,
#            masked_reference_genome = RepeatMasker.masked_genome
#        }
#    }

    output {
        File? clean_annotation = wf_align.clean_annotation
        IndexedReference clean_reference_index = wf_align.clean_reference_index

        Array[AlignedSample]? sr_bams = wf_align.sr_bams
        Array[AlignedSample]? lq_bams = wf_align.lq_bams
        Array[AlignedSample]? hq_bams = wf_align.hq_bams

        Array[AssembledSample]? sr_asms = wf_align.SR_gff
        Array[AssembledSample]? lq_asms = wf_align.LQ_gff
        Array[AssembledSample]? hq_asms = wf_align.HQ_gff

        File? sr_stringtie_stats = wf_align.SR_stringtie_stats
        File? sr_scallop_stats = wf_align.SR_scallop_stats
        File? lq_stats = wf_align.LQ_stats
        File? hq_stats = wf_align.HQ_stats

        File? portcullis_pass_tab = wf_align.pass_filtered_tab
        File? portcullis_pass_bed = wf_align.pass_filtered_bed
        File? portcullis_pass_gff3 = wf_align.pass_filtered_gff3

        File? portcullis_fail_tab = wf_align.fail_filtered_tab
        File? portcullis_fail_bed = wf_align.fail_filtered_bed
        File? portcullis_fail_gff3 = wf_align.fail_filtered_gff3

        # Array[Array[File]]? stats = wf_align.stats
        Array[Array[Array[File]]]? plots = wf_align.plots

        File? mikado_long_orfs = wf_main_mikado.long_orfs
        File? mikado_long_loci = wf_main_mikado.long_loci
        File? mikado_long_scores = wf_main_mikado.long_scores
        File? mikado_long_metrics = wf_main_mikado.long_metrics
        File? mikado_long_stats = wf_main_mikado.long_stats

        File? mikado_short_orfs = wf_main_mikado.short_orfs
        File? mikado_short_loci = wf_main_mikado.short_loci
        File? mikado_short_scores = wf_main_mikado.short_scores
        File? mikado_short_metrics = wf_main_mikado.short_metrics
        File? mikado_short_stats = wf_main_mikado.short_stats

        File? mikado_short_and_long_noLQ_orfs = wf_main_mikado.short_and_long_noLQ_orfs
        File? mikado_short_and_long_noLQ_loci = wf_main_mikado.short_and_long_noLQ_loci
        File? mikado_short_and_long_noLQ_scores = wf_main_mikado.short_and_long_noLQ_scores
        File? mikado_short_and_long_noLQ_metrics = wf_main_mikado.short_and_long_noLQ_metrics
        File? mikado_short_and_long_noLQ_stats = wf_main_mikado.short_and_long_noLQ_stats

        File? mikado_longHQ_orfs = wf_main_mikado.longHQ_orfs
        File? mikado_longHQ_loci = wf_main_mikado.longHQ_loci
        File? mikado_longHQ_scores = wf_main_mikado.longHQ_scores
        File? mikado_longHQ_metrics = wf_main_mikado.longHQ_metrics
        File? mikado_longHQ_stats = wf_main_mikado.longHQ_stats

        File? mikado_longLQ_orfs = wf_main_mikado.longLQ_orfs
        File? mikado_longLQ_loci = wf_main_mikado.longLQ_loci
        File? mikado_longLQ_scores = wf_main_mikado.longLQ_scores
        File? mikado_longLQ_metrics = wf_main_mikado.longLQ_metrics
        File? mikado_longLQ_stats = wf_main_mikado.longLQ_stats

#        IndexedReference masked_genome = RepeatMasker.masked_genome
#        Array[Array[File]]? maybe_exonerate_hits = Exonerate.exonerate_results
    }
}