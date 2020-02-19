version 1.0

import "mikado.wdl" as wfm
import "align.wdl" as waln
import "workflows/common/structs.wdl"
import "workflows/exonerate/wf_exonerate.wdl" as exonerate
import "workflows/repeat_masker/wf_repeat_masker.wdl" as repeatmasker
import "workflows/sanitize/wf_sanitize.wdl" as san
import "workflows/index/wf_index.wdl" as idx

workflow ei_annotation {
    input {
        File reference_genome
        Array[PRSample]? paired_samples
        Array[LRSample]? LQ_long_read_samples
        Array[LRSample]? HQ_long_read_samples
        Array[LabeledFasta]? protein_related_species
        File? annotation
        File? mikado_scoring_file
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

    call wfm.wf_mikado {
        input:
        mikado_scoring_file = mikado_scoring_file,
        reference_genome = wf_sanitize.indexed_reference,
        SR_align = wf_align.sr_gff
    }

    call repeatmasker.wf_repeat_masker as RepeatMasker {
        input:
        reference_fasta = wf_sanitize.indexed_reference.fasta
    }

    if (defined(protein_related_species)) {
        Array[LabeledFasta] def_protein_related_species = select_first([protein_related_species])
        call exonerate.wf_exonerate as Exonerate {
            input:
            related_species_protein = def_protein_related_species,
            masked_reference_genome = RepeatMasker.masked_genome
        }
    }

    output {
        File clean_reference = wf_sanitize.reference
        IndexedReference clean_reference_index = wf_sanitize.indexed_reference
        File? clean_annotation = wf_sanitize.annotation
        Array[File] gsnap_index = wf_index.gsnap_index
        Array[File] hisat_index = wf_index.hisat_index
        Array[File] star_index = wf_index.star_index

        Array[AssembledSample]? sr_asms = wf_align.sr_gff
        File mikado_long_config = wf_mikado.mikado_long_config
        File? mikado_long_orfs = wf_mikado.mikado_long_orfs

        File mikado_short_config = wf_mikado.mikado_short_config
        File? mikado_short_orfs = wf_mikado.mikado_short_orfs

        IndexedReference masked_genome = RepeatMasker.masked_genome

        Array[Array[File]]? maybe_exonerate_hits = Exonerate.exonerate_results
    }
}
