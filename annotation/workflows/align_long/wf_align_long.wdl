version 1.0

import "../common/structs.wdl"

workflow wf_align_long {
    input {
        File reference
        Array[LRSample] long_samples
        # File? annotation
    }

    scatter (sample in long_samples) {
        call Minimap2Long {
            input:
            long_sample = sample,
            reference = reference
        }

        call Minimap2Long2Gff {
            input:
            aligned_sample = Minimap2Long.aligned_sample
        }
    }

    output {
        Array[AlignedSample] bams = Minimap2Long.aligned_sample
        Array[AssembledSample] assemblies = Minimap2Long2Gff.assembled_sample
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
        File? LR
        File? iit
        Int? min_trimmed_coverage
        Int? min_identity
        String? strand
    }

    output {
        File gff = "gmap.out.gff"
    }

    command <<<
        gzcat ~{LR} | $(determine_gmap.py ~{reference}) --dir="$(dirname ~{gmap_index[0]})" --db=test_genome \
        --min-intronlength=20 --intronlength=2000 \
        ~{"-m " + iit} \
        ~{"--min-trimmed-coverage=" + min_trimmed_coverage} \
        ~{"--min-identity" + min_identity} \
        ~{"-z " + strand} \
        --format=gff3_match_cdna \
        --nthreads=4 > gmap.out.gff
    >>>
}

task Minimap2Long {
    input {
        LRSample long_sample
        File reference
    }

    output {
        AlignedSample aligned_sample = {"name": long_sample.name, "strand": long_sample.strand, "aligner": "minimap2", "bam": "minimap2.out.bam"}
    }

    command <<<
        minimap2 \
        -x splice \
        --cs=long \
        -G 2000 \
        -u b \
        -t 4 \
        -a -L --MD \
        --eqx -2 \
        --secondary=no \
        ~{reference} ~{long_sample.LR} | samtools view -bS - | samtools sort -@ 4 --reference ~{reference} -T minimap2.sort -o minimap2.out.bam -
    >>>
}

task Minimap2Long2Gff {
    input {
        AlignedSample aligned_sample
    }

    output {
        AssembledSample assembled_sample = {"name": aligned_sample.name, "strand": aligned_sample.strand, "assembly": "output.bed12"}
    }

    command <<<
        paftools.js splice2bed -m <(samtools view -h ~{aligned_sample.bam}) | correct_bed12_mappings.py > output.bed12
    >>>
}
