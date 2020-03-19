version 1.0
import "workflows/sanitize/wf_sanitize.wdl" as san
import "workflows/index/wf_index.wdl" as idx
import "workflows/common/structs.wdl"
import "workflows/align_short/wf_align_short.wdl" as aln_s
import "workflows/assembly_short/wf_assembly_short.wdl" as assm_s
import "workflows/portcullis/wf_portcullis.wdl" as portcullis_s
import "workflows/align_long/wf_align_long.wdl" as aln_l
import "workflows/assembly_long/wf_assembly_long.wdl" as assm_l

workflow wf_align {
    input {
        File reference_genome
        Array[PRSample]? paired_samples
        Array[LRSample]? LQ_long_read_samples
        Array[LRSample]? HQ_long_read_samples
        File? reference_annotation
        String HQ_assembler = "merge"
        String LQ_assembler = "stringtie"
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
        paired_samples: "Paired short read samples, each item is defined by a biological replicate name with one or more technical"
        LQ_long_read_samples: "Low quality long read samples, each item is defined by a name, it's strand and one or more long read files."
        HQ_long_read_samples: "High quality long read samples, each item is defined by a name, it's strand and one or more long read files."
        reference_annotation: "Pre-existing annotation that will used during the alignment process and to cleanup the splicing sites."
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
            reference_annotation = wf_sanitize.annotation,
            gsnap_index = wf_index.gsnap_index,
            hisat_index = wf_index.hisat_index,
            star_index = wf_index.star_index,
            alignment_resources = short_read_alignment_resources,
            sort_resources = short_read_alignment_sort_resources,
            stats_resources = short_read_stats_resources
        }

        call assm_s.wf_assembly_short {
            input:
            aligned_samples = wf_align_short.aligned_samples,
            reference_annotation = wf_sanitize.annotation,
            assembly_resources = short_read_assembly_resources
        }

        call portcullis_s.portcullis {
            input:
            reference = wf_sanitize.reference,
            annotation = wf_sanitize.annotation,
            aligned_samples = wf_align_short.aligned_samples
        }
    }

    if (defined(LQ_long_read_samples)) {
        Array[LRSample] def_lq_long_sample = select_first([LQ_long_read_samples])
        call aln_l.wf_align_long as LQ_align {
            input:
            reference = wf_sanitize.reference,
            is_hq = false,
            long_samples = def_lq_long_sample,
            bed_junctions = portcullis.bed,
            indexing_resources = long_read_indexing_resources,
            alignment_resources = long_read_alignment_resources
        }

        # Check what is defined for lq-long and run that

        call assm_l.wf_assembly_long as LQ_assembly {
            input:
            reference_annotation = wf_sanitize.annotation,
            aligned_samples = LQ_align.bams,
            assembler = LQ_assembler,
            assembly_resources = long_read_assembly_resources
        }
    }

    if (defined(HQ_long_read_samples)) {
        Array[LRSample] def_hq_long_sample = select_first([HQ_long_read_samples])
        call aln_l.wf_align_long as HQ_align {
            input:
            reference = wf_sanitize.reference,
            is_hq = true,
            long_samples = def_hq_long_sample,
            bed_junctions = portcullis.bed,
            indexing_resources = long_read_indexing_resources,
            alignment_resources = long_read_alignment_resources
        }

        # Check what is defined for hq-long and run that

        call assm_l.wf_assembly_long as HQ_assembly {
            input:
            reference_annotation = wf_sanitize.annotation,
            aligned_samples = HQ_align.bams,
            assembler = HQ_assembler,
            assembly_resources = long_read_assembly_resources
        }
    }

    output {
        File clean_reference = wf_sanitize.reference
        IndexedReference clean_reference_index = wf_sanitize.indexed_reference
        File? clean_annotation = wf_sanitize.annotation
        Array[File] gsnap_index = wf_index.gsnap_index
        Array[File] hisat_index = wf_index.hisat_index
        Array[File] star_index = wf_index.star_index

        Array[Array[File]]? stats = wf_align_short.stats
        Array[Array[Array[File]]]? plots = wf_align_short.plots


        File? filtered_tab = portcullis.tab
        File? filtered_bed = portcullis.bed
        File? filtered_gff3 = portcullis.gff3

        Array[AlignedSample]? sr_bams = wf_align_short.aligned_samples
        Array[AlignedSample]? lq_bams = LQ_align.bams
        Array[AlignedSample]? hq_bams = HQ_align.bams

        Array[AssembledSample]? SR_gff = wf_assembly_short.assemblies
        Array[AssembledSample]? LQ_gff = LQ_assembly.gff
        Array[AssembledSample]? HQ_gff = HQ_assembly.gff
    }
}