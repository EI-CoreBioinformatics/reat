version 1.0

import "./structs.wdl"
import "./rt_struct.wdl"

task SanitiseAnnotation {
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

    output {
        File sanitised_annotation = "reference.san.gtf"
    }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    command <<<
        set -euxo pipefail
        filepath=~{annotation}
        if [ ${filepath##*.} == "gff" ]
        then
            mikado util convert -of gtf ~{annotation} "reference.san.gtf"
        else
            ln ~{annotation} "reference.san.gtf"
        fi
    >>>
}

task SplitSequences {
    input {
        File sequences_file
        String prefix = "out"
        Int num_out_files = 10
    }

    output {
        Array[File] seq_files = glob(prefix+"*")
    }

    command <<<
        set -euxo pipefail
        seqtk split -n ~{num_out_files} ~{prefix} ~{sequences_file}
    >>>
}

task MergeFiles {
    input {
        Array[File] files_to_merge
        String output_filename
    }

    output {
        File out = output_filename
    }

    command <<<
    set -euxo pipefail
    cat ~{sep=" " files_to_merge} > ~{output_filename}
    >>>
}

task IndexFasta {
    input {
        File reference_fasta
        RuntimeAttr? runtime_attr_override
    }
    
    output {
        IndexedReference indexed_fasta = {"fasta": basename(reference_fasta), "fai": basename(reference_fasta)+".fai"}
    }

    command <<<
        set -euxo pipefail
        ln -s ~{reference_fasta} .
        samtools faidx ~{basename(reference_fasta)}
    >>>

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
}

task SanitiseFasta {
    input {
        File reference
        RuntimeAttr? runtime_attr_override
    }
    
    output {
        File sanitised_reference = "reference.san.fasta"
    }

    command <<<
        set -euxo pipefail
        sanitize_sequence_db.py -o "reference.san.fasta" ~{reference}
    >>>
    
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
}

task GSnapIndex {
    input {
        File reference
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        Array[File] index = glob("gsnapIndex/ref/*")
    }

    command <<<
        set -euxo pipefail
    mkdir gsnapIndex
    gmap_build --dir="gsnapIndex" --db=ref ~{reference}
    >>>
}

task Hisat2Index {
    input {
        File reference
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        Array[File] index = glob("ref*")
    }

    command <<<
        set -euxo pipefail
        hisat2-build ~{reference} "ref"
    >>>
}

task StarIndex {
    input {
        File reference
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 4
    
    output {
        Array[File] index = glob('starIndex/*')
    }

    command <<<
        set -euxo pipefail
        mkdir StarIndex
        STAR --runThreadN ~{cpus} --runMode genomeGenerate --genomeDir starIndex \
        --genomeFastaFiles ~{reference}
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

task tophatIndex {
    input {
        File reference
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 4
    
    output {
        Array[File] index = glob('ref*')
    }

    command <<<
        set -euxo pipefail
        bowtie2-build --threads ~{cpus} ~{reference} "ref"
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

task SummaryStats {
    input {
        Array[File] stats
    }

    output {
        File summary = "summary.stats.tsv"
    }

    command <<<
    cat ~{sep=" " stats} > summary.stats.tsv
    >>>
}

task Stats {
    input {
        File gff
    }

    output {
        File stats = basename(gff)+".stats.tsv"
    }

    command <<<
    mikado util stats ~{gff} --tab-stats ~{basename(gff)}.stats.tsv
    >>>

}
