version 1.0

import "./structs.wdl"
import "./rt_struct.wdl"

task SanitiseAnnotation {
    input {
        File? annotation
        String output_directory = "annotation"
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File sanitised_annotation = output_directory + "/reference.san.gtf"
    }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    command <<<
        set -euxo pipefail
        filepath=~{annotation}
        mkdir ~{output_directory}
        if [ ${filepath##*.} == "gff" ]
        then
            mikado util convert -of gtf ~{annotation} ~{output_directory}/reference.san.gtf
        else
            ln ~{annotation} ~{output_directory}/reference.san.gtf
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
        String output_directory = "reference"
        RuntimeAttr? runtime_attr_override
    }
    
    output {
        IndexedReference indexed_fasta = {
            "fasta": output_directory+"/"+basename(reference_fasta), 
            "fai": output_directory+"/"+basename(reference_fasta)+".fai"
        }
    }

    command <<<
        set -euxo pipefail
        mkdir ~{output_directory}
        cd ~{output_directory}
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

task TophatIndex {
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

task TranscriptAssemblySummaryStats {
    input {
        Array[File] stats
        String output_prefix
    }

    output {
        File summary = output_prefix + ".summary.stats.tsv"
    }

    command <<<
    mikado_summary_stats ~{sep=" " stats} > ~{output_prefix}.summary.stats.tsv
    >>>
}

task TranscriptAssemblyStats {
    input {
        File gff
    }

    output {
        File stats = "assembly_stats/" + basename(gff) + ".stats.tsv"
    }

    command <<<
    mkdir assembly_stats
    cd assembly_stats
    mikado util stats ~{gff} ~{basename(gff)}.stats.tsv
    >>>

}
