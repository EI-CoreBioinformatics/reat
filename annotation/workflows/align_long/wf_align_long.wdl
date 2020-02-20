version 1.0

import "../common/structs.wdl"
import "../common/rt_struct.wdl"
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
            scatter (LR in sample.LR) {
                call Minimap2Long {
                    input:
                    LR = LR,
                    strand = sample.strand,
                    name = sample.name,
                    reference = reference
                }
            }
            AlignedSample mm2_aligned_sample = {"name": sample.name, "strand":sample.strand, "aligner":"minimap2", "bam": Minimap2Long.bam}
        }
    }

    if (aligner == "gmap") {
        scatter (sample in long_samples) {
            call GMapIndex {
                input:
                reference = reference
            }
            scatter (LR in sample.LR) {
                call GMapLong {
                    input:
                    LR = LR,
                    gmap_index = GMapIndex.gmap_index,
                    strand = sample.strand,
                    name = sample.name,
                    reference = reference
                }
            }
            AlignedSample gmap_aligned_sample = {"name": sample.name, "strand":sample.strand, "aligner":"gmap", "bam": GMapLong.bam}
        }
    }

    Array[AlignedSample] def_alignments = select_first([mm2_aligned_sample, gmap_aligned_sample])

    output {
        Array[AlignedSample] bams = def_alignments
    }

}

task GMapIndex {
    input {
        File reference
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        Array[File] gmap_index = glob("gmapIndex/test_genome/*")
    }


    command <<<
        set -euxo pipefail
        gmap_build --dir=gmapIndex --db=test_genome ~{reference}
    >>>
}

task GMapExonsIIT {
    input {
        File? annotation
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        File iit = "gmap_exons.iit"
    }

    command <<<
        set -euxo pipefail
        gtf_genes ~{annotation} | iit_store -o gmap_exons.iit
    >>>
}

task GMapLong {
    input {
        File reference
        Array[File] gmap_index
        File LR
        String name
        String strand
        File? iit
        Int? min_trimmed_coverage
        Int? min_identity
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 16

    output {
        File bam = "gmap."+name+".bam"
    }

    command <<<
        set -euxo pipefail
        filename=$(basename -- "~{LR}")
        extension="${filename##*.}"

        in_pipe="cat ~{LR}"
        if [ "$extension" == "bam" ]
        then
            in_pipe="samtools fastq ~{LR}"
        elif [ "$extension" == "gz" ]
        then
            in_pipe="gunzip -c ~{LR}"
        fi

        $in_pipe | $(determine_gmap.py ~{reference}) --dir="$(dirname ~{gmap_index[0]})" --db=test_genome \
        --min-intronlength=20 --intronlength=2000 \
        ~{"-m " + iit} \
        ~{"--min-trimmed-coverage=" + min_trimmed_coverage} \
        ~{"--min-identity" + min_identity} \
        ~{"-z " + strand} \
        --format=samse \
        --nthreads="~{cpus}" | samtools view -F 4 -F 0x900 -bS - | samtools sort -@ 4 --reference ~{reference} -T gmap.sort -o gmap.~{name}.bam -
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 16,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

}

task Minimap2Long {
    input {
        File LR
        String name
        String strand
        File reference
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 16

    output {
        File bam = "minimap2."+name+".bam"
    }

    command <<<
        set -euxo pipefail
        # Replace long_sample.LR with samtools fastq if suffix is bam
        filename=$(basename -- "~{LR}")
        extension="${filename##*.}"

        in_pipe="cat ~{LR}"
        if [ "$extension" == "bam" ]
        then
            in_pipe="samtools fastq ~{LR}"
        elif [ "$extension" == "gz" ]
        then
            in_pipe="gunzip -c ~{LR}"
        fi
        
        $in_pipe | \
        minimap2 \
        -x splice \
        --cs=long \
        -G 2000 \
        -u b \
        -t ~{cpus} \
        -a -L --MD \
        --eqx -2 \
        --secondary=no \
        ~{reference} - | samtools view -F 4 -F 0x900 -bS - | samtools sort -@ 4 --reference ~{reference} -T minimap2.sort -o minimap2.~{name}.bam -
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 16,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}