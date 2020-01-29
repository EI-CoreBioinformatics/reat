version 1.0

import "../common/structs.wdl"

workflow wf_assembly_long {
    input {
        File reference
        Array[AlignedSample] bams
    }

    call stringtie_long {
        input:
        reference = reference,
        bams = bams
    }

    output {
        File gtf = stringtie_long.results
    }
}

task stringtie_long {
    input {
        File reference
        Array[AlignedSample] bams
    }

    output {
        File results = "result.gff"
    }

    command <<<
    stringtie ~{"-G " + reference} -L ~{sep=" " bams} -o "result.gff"
    >>>
}