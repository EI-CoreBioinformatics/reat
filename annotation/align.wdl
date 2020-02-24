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
    }

    parameter_meta {
        reference_genome: "Reference genome to align against"
        paired_samples:
        "
        Paired short read samples, each item is defined by a biological replicate name with one or more technical
        replicates. Technical replicates are defined by a name, R1, R2 and strand.
        Example:
        [
            {
                \"biological_sample_name\": \"Ara\",
        	    \"technical_samples\":
        	    [
        	        {
        		        \"name\": \"Ara1\",
        		        \"strand\": \"fr-firststrand\",
        		        \"R1\": \"Ara1_R1.fastq.gz\",
        		        \"R2\": \"Ara1_R2.fastq.gz\"
        		    },
        		    {
        		        \"name\": \"Ara2\",
        		        \"strand\": \"fr-firststrand\",
        		        \"R1\": \"Ara2_R1.fastq.gz\",
        		        \"R2\": \"Ara2_R2.fastq.gz\"
        		    },
        	    ]
        	}
        ]
        "
        LQ_long_read_samples:
        "
        Low quality long read samples, each item is defined by a name, it's strand and one or more long read files.
        They optionally contain a is_ref flag which will carry through any transcripts present in the files all the way
        to the final outputs and an optional score flag which will apply a weight to the transcripts score in mikado.
        "
        HQ_long_read_samples:
        "
        High quality long read samples, each item is defined by a name, it's strand and one or more long read files.
        They optionally contain a is_ref flag which will carry through any transcripts present in the files all the way
        to the final outputs and an optional score flag which will apply a weight to the transcripts score in mikado.
        "
        reference_annotation:
        "
        Pre-existing annotation that will used during the alignment process and to cleanup the splicing sites.
        "
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
            aligned_samples = wf_align_short.aligned_samples,
            annotation = wf_sanitize.annotation
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

        Array[Array[File]]? stats = wf_align_short.stats
        Array[Array[Array[File]]]? plots = wf_align_short.plots


        File? filtered_tab = portcullis.tab
        File? filtered_bed = portcullis.bed
        File? filtered_gff3 = portcullis.gff3

        Array[AlignedSample]? sr_bams = wf_align_short.aligned_samples
        Array[AlignedSample]? lq_bams = LQ_align.bams
        Array[AlignedSample]? hq_bams = HQ_align.bams

        Array[AssembledSample]? SR_gff = wf_assembly_short.assemblies
        Array[File]? LQ_gff = LQ_assembly.gff
        Array[File]? HQ_gff = HQ_assembly.gff
    }
}