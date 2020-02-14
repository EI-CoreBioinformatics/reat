version 1.0

import "./structs.wdl"
import "./rt_struct.wdl"

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
        seqtk split -n ~{num_out_files} ~{prefix} ~{sequences_file}
    >>>
}

task IndexFasta {
    input {
        File reference_fasta
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
        IndexedReference indexed_fasta = {"fasta": basename(reference_fasta), "fai": basename(reference_fasta)+".fai"}
    }

    command <<<
        ln -s ~{reference_fasta} .
        samtools faidx ~{basename(reference_fasta)}
    >>>
}

task sanitizeFasta {
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
        File sanitised_reference = "reference.san.fasta"
    }

    command <<<
        sanitize_sequence_db.py -o "reference.san.fasta" ~{reference}
    >>>
}
