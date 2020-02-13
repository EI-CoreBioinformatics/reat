version 1.0

import "../common/structs.wdl"

workflow wf_align_long {
    input {
        File reference
        Array[LRSample] long_samples
        # File? junctions
        String aligner = "minimap2"
    }
    
    # Add aligner option
    if (aligner == "minimap2") {
        scatter (sample in long_samples) {
            call Minimap2Long {
                input:
                long_sample = sample,
                reference = reference
            }
        }
    }

    if (aligner == "gmap") {
        scatter (sample in long_samples) {
            call GMapIndex {
                input:
                reference = reference
            }
            call GMapLong {
                input:
                gmap_index = GMapIndex.gmap_index,
                reference = reference,
                sample = sample
            }
        }
    }

    Array[AlignedSample] def_alignments = select_first([GMapLong.aligned_sample, Minimap2Long.aligned_sample])

    output {
        Array[AlignedSample] bams = def_alignments
    }

}

task GMapIndex {
    input {
        File reference
    }

    output {
        Array[File] gmap_index = glob("gmapIndex/test_genome/*")
    }

    command <<<
        gmap_build --dir=gmapIndex --db=test_genome ~{reference}
    >>>
}

task GMapExonsIIT {
    input {
        File? annotation
    }

    output {
        File iit = "gmap_exons.iit"
    }

    command <<<
        gtf_genes ~{annotation} | iit_store -o gmap_exons.iit
    >>>
}

task GMapLong {
    input {
        File reference
        Array[File] gmap_index
        LRSample sample
        File? iit
        Int? min_trimmed_coverage
        Int? min_identity
        String? strand
    }

    output {
        AlignedSample aligned_sample = {"name": sample.name, "strand": sample.strand, "aligner": "gmap", "bam": "gmap."+sample.name+".sam"}
    }

    command <<<
        gzcat ~{sample.LR} | $(determine_gmap.py ~{reference}) --dir="$(dirname ~{gmap_index[0]})" --db=test_genome \
        --min-intronlength=20 --intronlength=2000 \
        ~{"-m " + iit} \
        ~{"--min-trimmed-coverage=" + min_trimmed_coverage} \
        ~{"--min-identity" + min_identity} \
        ~{"-z " + strand} \
        --format=samse \
        --nthreads=4 | samtools view -F 4 -F 0x900 -bS - | samtools sort -@ 4 --reference ~{reference} -T gmap.sort -o gmap.~{sample.name}.bam -
    >>>
}

task Minimap2Long {
    input {
        LRSample long_sample
        File reference
    }

    output {
        AlignedSample aligned_sample = {"name": long_sample.name, "strand": long_sample.strand, "aligner": "minimap2", "bam": "minimap2."+long_sample.name+".bam"}
    }

    command <<<
        # Replace long_sample.LR with samtools fastq if suffix is bam
        
        minimap2 \
        -x splice \
        --cs=long \
        -G 2000 \
        -u b \
        -t 4 \
        -a -L --MD \
        --eqx -2 \
        --secondary=no \
        ~{reference} ~{long_sample.LR} | samtools view -F 4 -F 0x900 -bS - | samtools sort -@ 4 --reference ~{reference} -T minimap2.sort -o minimap2.~{long_sample.name}.bam -
    >>>
}