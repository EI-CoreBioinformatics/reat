version 1.0
import "workflows/mikado/wf_mikado.wdl" as mikado
import "workflows/portcullis/wf_portcullis.wdl" as portcullis_s
import "workflows/assembly_short/wf_assembly_short.wdl" as assm_s
import "workflows/align_short/wf_align_short.wdl" as aln_s
import "workflows/align_long/wf_align_long.wdl" as aln_l
import "workflows/sanitize/wf_sanitize.wdl" as san
import "workflows/index/wf_index.wdl" as idx

workflow ei_annotation {
    input {
        File short_R1
        File? short_R2
        File? long_R
        File reference_genome
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
        R1 = short_R1,
        R2 = short_R2,
        annotation = wf_sanitize.annotation,
        gsnap_index = wf_index.gsnap_index,
        hisat_index = wf_index.hisat_index,
        star_index = wf_index.star_index
    }

    # Sort and index the bams

    call assm_s.wf_assembly_short {
        input:
        bams = wf_align_short.indexed_bams,
        annotation = wf_sanitize.annotation
    }

    if (defined(long_R)) {
        call aln_l.wf_align_long {
            input:
            reference = wf_sanitize.reference,
            # annotation = wf_sanitize.annotation,
            LR = long_R,
            star_index = wf_index.star_index
        }
    }

    call portcullis_s.portcullis {
        input:
        reference = wf_sanitize.reference,
        annotation = wf_sanitize.annotation,
        bams = wf_align_short.indexed_bams
    }

    Array[File]? long_assemblies_valid = wf_align_long.gff
    call mikado.wf_mikado {
        input:
        scoring_file = mikado_scoring_file,
        reference_annotation =  wf_sanitize.annotation,
        assemblies = wf_assembly_short.assemblies,
        long_assemblies = long_assemblies_valid
    }

    output {
        File clean_reference = wf_sanitize.reference
        File clean_reference_index = wf_sanitize.index
        File? clean_annotation = wf_sanitize.annotation
        Array[File] gsnap_index = wf_index.gsnap_index
        Array[File] hisat_index = wf_index.hisat_index
        Array[File] star_index = wf_index.star_index

        Array[Pair[File,File]] bams = wf_align_short.indexed_bams
        Array[File] stats = wf_align_short.stats
        Array[Array[File]] plots = wf_align_short.plots
        Array[File] short_assemblies = wf_assembly_short.assemblies

        File filtered_tab = portcullis.tab
        File filtered_bed = portcullis.bed
        File filtered_gff3 = portcullis.gff3

        Array[File?]? l_bams = wf_align_long.bams
        Array[File?]? l_gff = wf_align_long.gff

        File mikado_config = wf_mikado.mikado_config
    }
}