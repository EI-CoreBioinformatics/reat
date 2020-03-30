version 1.0

import "mikado.wdl" as wfm
import "align.wdl" as waln
import "subworkflows/common/structs.wdl"
# import "subworkflows/exonerate/wf_exonerate.wdl" as exonerate
# import "subworkflows/repeat_masker/wf_repeat_masker.wdl" as repeatmasker
import "subworkflows/sanitise/wf_sanitise.wdl" as san
import "subworkflows/index/wf_index.wdl" as idx

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
        RuntimeAttr? short_read_alignment_resources
        RuntimeAttr? short_read_alignment_sort_resources
        RuntimeAttr? short_read_stats_resources
        RuntimeAttr? short_read_assembly_resources
        RuntimeAttr? long_read_indexing_resources
        RuntimeAttr? long_read_alignment_resources
        RuntimeAttr? long_read_assembly_resources
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

    call idx.wf_index {
        input:
        reference = wf_sanitise.reference
    }

    call waln.wf_align {
        input:
        reference_genome = reference_genome,
        paired_samples = paired_samples,
        LQ_long_read_samples = LQ_long_read_samples,
        HQ_long_read_samples = HQ_long_read_samples,
        reference_annotation = wf_sanitise.annotation,
        short_read_alignment_resources = short_read_alignment_resources,
        short_read_alignment_sort_resources = short_read_alignment_sort_resources,
        short_read_stats_resources = short_read_stats_resources,
        short_read_assembly_resources = short_read_assembly_resources,
        long_read_indexing_resources = long_read_indexing_resources,
        long_read_alignment_resources = long_read_alignment_resources,
        long_read_assembly_resources = long_read_assembly_resources
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
        File clean_reference = wf_align.clean_reference
        IndexedReference clean_reference_index = wf_align.clean_reference_index
        File? clean_annotation = wf_align.clean_annotation
        Array[File] gsnap_index = wf_align.gsnap_index
        Array[File] hisat_index = wf_align.hisat_index
        Array[File] star_index = wf_align.star_index

        Array[AssembledSample]? sr_asms = wf_align.SR_gff

        File? mikado_long_config = wf_main_mikado.long_config
        File? mikado_long_orfs = wf_main_mikado.long_orfs
        File? mikado_short_config = wf_main_mikado.short_config
        File? mikado_short_orfs = wf_main_mikado.short_orfs

#        IndexedReference masked_genome = RepeatMasker.masked_genome
#        Array[Array[File]]? maybe_exonerate_hits = Exonerate.exonerate_results
    }
}
