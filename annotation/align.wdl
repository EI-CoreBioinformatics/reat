version 1.0
import "workflows/common/structs.wdl"
import "workflows/portcullis/wf_portcullis.wdl" as portcullis_s
import "workflows/assembly_short/wf_assembly_short.wdl" as assm_s
import "workflows/align_short/wf_align_short.wdl" as aln_s
import "workflows/align_long/wf_align_long.wdl" as aln_l
import "workflows/assembly_long/wf_assembly_long.wdl" as assm_l
import "workflows/sanitize/wf_sanitize.wdl" as san
import "workflows/index/wf_index.wdl" as idx

workflow wf_align {
    input {
        File reference_genome
        Array[PRSample]? paired_samples
        Array[LRSample]? LQ_long_read_samples
        Array[LRSample]? HQ_long_read_samples
        File? reference_annotation
    }

    call san.wf_sanitize {
        input:
        reference_genome = reference_genome,
        in_annotation = reference_annotation
    }

    call idx.wf_index {
        input:
        reference = wf_sanitize.reference
    }

    if (defined(paired_samples)) {
        Array[PRSample] def_paired_samples = select_first([paired_samples])
        call aln_s.wf_align_short {
            input:
            samples = def_paired_samples,
            annotation = wf_sanitize.annotation,
            gsnap_index = wf_index.gsnap_index,
            hisat_index = wf_index.hisat_index,
            star_index = wf_index.star_index
        }
        call assm_s.wf_assembly_short {
            input:
            aligned_samples = wf_align_short.indexed_aligned_samples,
            annotation = wf_sanitize.annotation
        }
        call portcullis_s.portcullis {
            input:
            reference = wf_sanitize.reference,
            annotation = wf_sanitize.annotation,
            aligned_samples = wf_align_short.indexed_aligned_samples
        }
    }
    # Sort and index the bams

    if (defined(LQ_long_read_samples)) {
        Array[LRSample] def_lq_long_sample = select_first([LQ_long_read_samples])
        call aln_l.wf_align_long as LQ_align {
            input:
            reference = wf_sanitize.reference,
            long_samples = def_lq_long_sample
        }

        # Check what is defined for lq-long and run that

        call assm_l.wf_assembly_long as LQ_assembly {
            input:
            reference_annotation = wf_sanitize.annotation,
            aligned_samples = LQ_align.bams
        }
    }

    if (defined(HQ_long_read_samples)) {
        Array[LRSample] def_hq_long_sample = select_first([HQ_long_read_samples])
        call aln_l.wf_align_long as HQ_align {
            input:
            reference = wf_sanitize.reference,
            long_samples = def_hq_long_sample
        }

        # Check what is defined for hq-long and run that
        
        call assm_l.wf_assembly_long as HQ_assembly {
            input:
            reference_annotation = wf_sanitize.annotation,
            aligned_samples = HQ_align.bams
        }
    }

    output {
        File clean_reference = wf_sanitize.reference
        IndexedReference clean_reference_index = wf_sanitize.indexed_reference
        File? clean_annotation = wf_sanitize.annotation
        Array[File] gsnap_index = wf_index.gsnap_index
        Array[File] hisat_index = wf_index.hisat_index
        Array[File] star_index = wf_index.star_index

        Array[IndexedAlignedSample]? bams = wf_align_short.indexed_aligned_samples
        Array[File]? stats = wf_align_short.stats
        Array[Array[File]]? plots = wf_align_short.plots

        Array[AssembledSample]? sr_gff = wf_assembly_short.assemblies

        File? filtered_tab = portcullis.tab
        File? filtered_bed = portcullis.bed
        File? filtered_gff3 = portcullis.gff3

        Array[AlignedSample]? lq_bams = LQ_align.bams
        Array[AlignedSample]? hq_bams = HQ_align.bams

        Array[File]? lq_gff = LQ_assembly.gff
        Array[File]? hq_gff = HQ_assembly.gff
    }
}