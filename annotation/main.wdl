version 1.0
import "workflows/structs/structs.wdl"
import "workflows/mikado/wf_mikado.wdl" as mikado
import "workflows/exonerate/wf_exonerate.wdl" as exonerate
import "workflows/repeat_masker/wf_repeat_masker.wdl" as repeatmasker
import "workflows/portcullis/wf_portcullis.wdl" as portcullis_s
import "workflows/assembly_short/wf_assembly_short.wdl" as assm_s
import "workflows/align_short/wf_align_short.wdl" as aln_s
import "workflows/align_long/wf_align_long.wdl" as aln_l
import "workflows/sanitize/wf_sanitize.wdl" as san
import "workflows/index/wf_index.wdl" as idx

workflow ei_annotation {
    input {
        Array[PRSample] paired_samples
        Array[LRSample]? long_read_samples
        File reference_genome
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

    call aln_s.wf_align_short {
        input:
        samples = paired_samples,
        annotation = wf_sanitize.annotation,
        gsnap_index = wf_index.gsnap_index,
        hisat_index = wf_index.hisat_index,
        star_index = wf_index.star_index
    }

    # Sort and index the bams

    call assm_s.wf_assembly_short {
        input:
        aligned_samples = wf_align_short.indexed_aligned_samples,
        annotation = wf_sanitize.annotation
    }

    if (defined(long_read_samples)) {
        Array[LRSample] def_long_sample = select_first([long_read_samples])
        call aln_l.wf_align_long {
            input:
            reference = wf_sanitize.reference,
            # annotation = wf_sanitize.annotation,
            long_sample = def_long_sample,
            star_index = wf_index.star_index
        }
    }

    call portcullis_s.portcullis {
        input:
        reference = wf_sanitize.reference,
        annotation = wf_sanitize.annotation,
        aligned_samples = wf_align_short.indexed_aligned_samples
    }

    Array[AssembledSample]? long_assemblies_valid = wf_align_long.assemblies
    call mikado.wf_mikado as Mikado_long {
        input:
        scoring_file = mikado_scoring_file,
        indexed_reference =  wf_sanitize.indexed_reference,
        assemblies = wf_assembly_short.assemblies,
        long_assemblies = long_assemblies_valid,
        junctions = portcullis.bed
    }

    call mikado.wf_mikado as Mikado_short {
        input:
        scoring_file = mikado_scoring_file,
        indexed_reference =  wf_sanitize.indexed_reference,
        assemblies = wf_assembly_short.assemblies,
        junctions = portcullis.bed
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

        Array[IndexedAlignedSample] bams = wf_align_short.indexed_aligned_samples
        Array[File] stats = wf_align_short.stats
        Array[Array[File]] plots = wf_align_short.plots

        Array[AssembledSample] short_assemblies = wf_assembly_short.assemblies

        File filtered_tab = portcullis.tab
        File filtered_bed = portcullis.bed
        File filtered_gff3 = portcullis.gff3

        Array[AlignedSample?]? l_bams = wf_align_long.bams
        Array[AssembledSample?]? l_gff = wf_align_long.assemblies

        File mikado_long_config = Mikado_long.mikado_config
        File? mikado_long_orfs = Mikado_long.orfs

        File mikado_short_config = Mikado_short.mikado_config
        File? mikado_short_orfs = Mikado_short.orfs

        IndexedReference masked_genome = RepeatMasker.masked_genome

        Array[Array[File]]? maybe_exonerate_hits = Exonerate.exonerate_results
    }
}
