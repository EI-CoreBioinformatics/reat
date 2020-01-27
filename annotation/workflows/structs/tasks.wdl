version 1.0

import "./structs.wdl"

task SplitSequences {
    input {
        File sequences_file
        String prefix = "out"
        Int num_out_files = 10
    }

    output {
        Array[File] seq_files = glob("out*")
    }

    command <<<
        seqtk split -n ~{num_out_files} ~{prefix} ~{sequences_file}
    >>>
}

task IndexFasta {
    input {
        File reference_fasta
    }

    output {
        IndexedReference indexed_fasta = {"fasta": basename(reference_fasta), "fai": basename(reference_fasta)+".fai"}
    }

    command <<<
        ln -s ~{reference_fasta} .
        samtools faidx ~{basename(reference_fasta)}
    >>>
}
