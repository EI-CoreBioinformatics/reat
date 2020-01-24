version 1.0

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
