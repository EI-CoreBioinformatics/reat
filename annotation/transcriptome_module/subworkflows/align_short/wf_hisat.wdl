version 1.0

import "../common/structs.wdl"
import "../common/rt_struct.wdl"

workflow wf_Hisat{
    input {
        PRSample sample
        File? splice_sites
        Array[File] index
        RuntimeAttr? alignment_resources
    }

    scatter(PR in sample.read_pair) {
        call Hisat {
            input:
            sites = splice_sites,
            strand = sample.strand,
            name = sample.name,
            sample = PR,
            index = index,
            runtime_attr_override = alignment_resources
        }
    }

    output {
        AlignedSample aligned_sample = object {
            bam: Hisat.bam, 
            strand: sample.strand, 
            aligner: "hisat", 
            name: sample.name, 
            merge: sample.merge
            }
    }
}


task Hisat {
    input {
        Array[File] index
        File? sites
        ReadPair sample
        String strand
        String name
        String rp_name = name+"."+basename(sample.R1)
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    
    output {
        File bam = rp_name + ".hisat.bam"
    }

    command <<<
        set -euxo pipefail
        strandness=""
        case "~{strand}" in
            fr-firststrand)
            strandness="--rna-strandness=RF"
            ;;
            fr-secondstrand)
            strandness="--rna-strandness=FR"
            ;;
            f)
            strandness="--rna-strandness=F"
            ;;
            r)
            strandness="--rna-strandness=R"
            ;;
        esac
    hisat2 -p ~{cpus} -x ~{sub(index[0], "\\.\\d\\.ht2l?", "")} \
    ${strandness} \
    --min-intronlen=20 \
    --max-intronlen=2000 \
    ~{"--known-splicesite-infile " + sites} \
    -1 ~{sample.R1} ~{"-2 " + sample.R2} | samtools sort -@ 4 - > "~{rp_name}.hisat.bam"
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