version 1.0

import "mikado.wdl" as wfm
import "align.wdl" as waln
import "workflows/common/structs.wdl"
# import "workflows/exonerate/wf_exonerate.wdl" as exonerate
# import "workflows/repeat_masker/wf_repeat_masker.wdl" as repeatmasker
import "workflows/sanitize/wf_sanitize.wdl" as san
import "workflows/index/wf_index.wdl" as idx

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
        reference_annotation: "Pre-existing annotation that will used during the alignment process and to cleanup the splicing sites."
        mikado_scoring_file: "Mikado scoring file"
    }

    call san.wf_sanitize {
        input:
        reference_genome = reference_genome,
        in_annotation = annotation
    }

    call idx.wf_index {
        input:
        reference = wf_sanitize.reference
    }

    call waln.wf_align {
        input:
        reference_genome = wf_sanitize.reference,
        paired_samples = paired_samples,
        LQ_long_read_samples = LQ_long_read_samples,
        HQ_long_read_samples = HQ_long_read_samples,
        reference_annotation = wf_sanitize.annotation
    }

    call wfm.wf_main_mikado {
        input:
        mikado_scoring_file = mikado_scoring_file,
        reference_genome = wf_sanitize.indexed_reference,
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
        File clean_reference = wf_sanitize.reference
        IndexedReference clean_reference_index = wf_sanitize.indexed_reference
        File? clean_annotation = wf_sanitize.annotation
        Array[File] gsnap_index = wf_index.gsnap_index
        Array[File] hisat_index = wf_index.hisat_index
        Array[File] star_index = wf_index.star_index

        Array[AssembledSample]? sr_asms = wf_align.SR_gff
        File mikado_long_config = wf_main_mikado.mikado_long_config
        File? mikado_long_orfs = wf_main_mikado.mikado_long_orfs

        File mikado_short_config = wf_main_mikado.mikado_short_config
        File? mikado_short_orfs = wf_main_mikado.mikado_short_orfs

#        IndexedReference masked_genome = RepeatMasker.masked_genome
#        Array[Array[File]]? maybe_exonerate_hits = Exonerate.exonerate_results
    }
}
