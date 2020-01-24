version 1.0

import "../../structs/structs.wdl"
import "../../structs/tasks.wdl" as tasks

workflow wf_exonerate {
    input {
        IndexedReference masked_reference_genome
        IndexedReference clean_indexed_protein_db
    }

    call tasks.SplitSequences {
        input:
        sequences_file = clean_indexed_protein_db.fasta
    }

    output {
        String gff3 = "gff3"
    }
}
