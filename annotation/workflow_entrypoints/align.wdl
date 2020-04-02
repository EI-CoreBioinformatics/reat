version 1.0

import "subworkflows/index/wf_index.wdl" as idx
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
        String LQ_aligner = "minimap2"
        String HQ_aligner = "gmap"
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

    RuntimeAttr default_runtime_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
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

    call idx.wf_index {
        input:
        reference = wf_sanitise.reference
    }

    if (defined(paired_samples)) {
        Array[PRSample] def_paired_samples = select_first([paired_samples])
        call aln_s.wf_align_short {
            input:
            samples = def_paired_samples,
            reference_annotation = wf_sanitise.annotation,
            gsnap_index = wf_index.gsnap_index,
            hisat_index = wf_index.hisat_index,
            star_index = wf_index.star_index,
            alignment_resources = select_first([short_read_alignment_resources, default_runtime_attr]),
            sort_resources = select_first([short_read_alignment_sort_resources, default_runtime_attr]),
            stats_resources = select_first([short_read_stats_resources, default_runtime_attr])
        }

        call assm_s.wf_assembly_short {
            input:
            aligned_samples = wf_align_short.aligned_samples,
            reference_annotation = wf_sanitise.annotation,
            assembly_resources = select_first([default_runtime_attr, short_read_assembly_resources])
        }

        call portcullis_s.portcullis {
            input:
            reference = wf_sanitise.reference,
            annotation = wf_sanitise.annotation,
            aligned_samples = wf_align_short.aligned_samples
        }
    }

    if (defined(LQ_long_read_samples)) {
        Array[LRSample] def_lq_long_sample = select_first([LQ_long_read_samples])
        call aln_l.wf_align_long as LQ_align {
            input:
            reference = wf_sanitise.reference,
            is_hq = false,
            aligner = LQ_aligner,
            long_samples = def_lq_long_sample,
            bed_junctions = portcullis.bed,
            indexing_resources = select_first([long_read_indexing_resources,default_runtime_attr]),
            alignment_resources = select_first([long_read_alignment_resources,default_runtime_attr])
        }

        # Check what is defined for lq-long and run that

        call assm_l.wf_assembly_long as LQ_assembly {
            input:
            reference_annotation = wf_sanitise.annotation,
            aligned_samples = LQ_align.bams,
            assembler = LQ_assembler,
            assembly_resources = select_first([long_read_assembly_resources, default_runtime_attr])
        }
    }

    if (defined(HQ_long_read_samples)) {
        Array[LRSample] def_hq_long_sample = select_first([HQ_long_read_samples])
        call aln_l.wf_align_long as HQ_align {
            input:
            reference = wf_sanitise.reference,
            is_hq = true,
            aligner = HQ_aligner,
            long_samples = def_hq_long_sample,
            bed_junctions = portcullis.bed,
            indexing_resources = select_first([long_read_indexing_resources, default_runtime_attr]),
            alignment_resources = select_first([long_read_alignment_resources, default_runtime_attr])
        }

        # Check what is defined for hq-long and run that

        call assm_l.wf_assembly_long as HQ_assembly {
            input:
            reference_annotation = wf_sanitise.annotation,
            aligned_samples = HQ_align.bams,
            assembler = HQ_assembler,
            assembly_resources = select_first([long_read_assembly_resources, default_runtime_attr])
        }
    }

    output {
        File clean_reference = wf_sanitise.reference
        IndexedReference clean_reference_index = wf_sanitise.indexed_reference
        File? clean_annotation = wf_sanitise.annotation
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