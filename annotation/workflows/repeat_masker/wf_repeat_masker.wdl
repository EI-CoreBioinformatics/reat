version 1.0

import "../common/structs.wdl"
import "../common/tasks.wdl"

workflow wf_repeat_masker {
    input {
        File reference_fasta
        Boolean run_modeller = true
        Boolean retrieve_known = false
        Array[File]? safe_proteins
        String? clade
        String? specie
    }

    if (run_modeller) {
        call BuildModellerDB {
            input:
            reference_fasta = reference_fasta
        }

        call RepeatModeller {
            input:
            db = BuildModellerDB.db
        }
    }

    if (retrieve_known) {
        if (defined(clade) || defined(specie)) {
            call RetrieveLibraries {
                input:
                clade = clade,
                specie = specie
            }
        }
    }

    call CreateLibrary {
        input:
        repeat_modeler_consensi = RepeatModeller.consensi,
        retrieved_libraries = RetrieveLibraries.retrieved,
    }

    call RepeatMasker {
        input:
        rm_library = CreateLibrary.repeat_library,
        reference_fasta = reference_fasta
    }

    call tasks.IndexFasta {
        input: 
        reference_fasta = RepeatMasker.masked_genome
    }

    output {
        IndexedReference masked_genome = IndexFasta.indexed_fasta
    }
}

task BuildModellerDB {
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
        Array[File] db = glob("genome.*")
    }

    command <<<
    BuildDatabase -name genome -engine ncbi ~{reference_fasta}
    >>>
}

task RepeatModeller {
    input {
        Array[File] db
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
        Array[File] consensi = glob("RM*/consensi.fa.classified")
    }

    command <<<
    ln -s ~{sep=" " db} -t .
    RepeatModeler -engine ncbi -pa 4 -database genome
    >>>
}

task RetrieveLibraries {
    input {
        String? clade
        String? specie
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
        File retrieved = "retrieved.fa"
    }

    command <<<
    queryRepeatDatabase.pl ~{"-species " + specie} ~{"-clade " + clade} > "retrieved.fa"
    >>>
}

task CreateLibrary {
    input {
        Array[File]? repeat_modeler_consensi
        File? retrieved_libraries
        File? extra
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
        File repeat_library = "rm_library.fa"
    }

    command <<<
    cat ~{sep=" " repeat_modeler_consensi} ~{retrieved_libraries} ~{extra} > rm_library.fa
    >>>
}

task RepeatMasker {
    input {
        File rm_library
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
        File masked_genome = basename(reference_fasta) + ".masked"
        File masked_table = basename(reference_fasta) + ".out"
        File table = basename(reference_fasta) + ".tbl"
    }

    command <<<
    RepeatMasker -nolow -xsmall -dir . -gff -lib ~{rm_library} -pa 4 ~{reference_fasta}
    >>>
}
