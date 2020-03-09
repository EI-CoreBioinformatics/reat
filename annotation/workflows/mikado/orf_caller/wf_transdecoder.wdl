version 1.0

import "../../common/structs.wdl"
import "../align_protein/wf_protein_aligner.wdl" as prt_aln

workflow wf_transdecoder {
    input {
        File prepared_transcripts
        String program = "blastp"
        File orf_proteins
    }

    call TransdecoderLongOrf {
        input:
            prepared_transcripts = prepared_transcripts
    }

#    if (defined(db)) {
#        File def_db = select_first([db])
        call prt_aln.SanitiseProteinBlastDB { # Check if this is an input for both homology + orf_calling or just the same task
            input:
                db = orf_proteins
#                db = def_db
        }

        if (program == "blastp") {
            call prt_aln.BlastIndex as BlastIndex {
                input:
                target = SanitiseProteinBlastDB.clean_db
            }

            call prt_aln.BlastAlign {
                input:
                query = TransdecoderLongOrf.pep,
                outfmt = "6",
                output_filename = "mikado_blast_orfcalling.txt",
                index = BlastIndex.index
            }
        }
        
        if (program == "diamond") {
            call prt_aln.DiamondIndex as DiamondIndex {
                input:
                target = SanitiseProteinBlastDB.clean_db
            }

            call prt_aln.DiamondAlign {
                input:
                query = TransdecoderLongOrf.pep,
                index = DiamondIndex.index
            }
        }

        File maybe_blast_aligned_proteins = select_first([BlastAlign.out, DiamondAlign.out])
#    }

    call TransdecoderPredict {
        input:
        prepared_transcripts = prepared_transcripts,
        transdecoder_wd = TransdecoderLongOrf.transdecoder_wd,
        blast_aligned_proteins = maybe_blast_aligned_proteins
    }

    output {
        File long_orfs = TransdecoderLongOrf.orfs
        File clean_db = SanitiseProteinBlastDB.clean_db
        File final_orfs = TransdecoderPredict.bed
    }
}

task TransdecoderLongOrf {
    input {
        File prepared_transcripts
        Int minprot = 20 # Minimum protein length (TODO make parameter)
        String gencode = "Universal" # Genetic code (TODO Make parameter) 
        # More gencode info https://github.com/TransDecoder/TransDecoder/blob/master/PerlLib/Nuc_translator.pm
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
        File orfs = "transdecoder/longest_orfs.gff3"
        File pep = "transdecoder/longest_orfs.pep"
        Array[File] transdecoder_wd = glob("transdecoder/*")
    }

    command <<<
        set -euxo pipefail
        Transdecoder.LongOrfs -m "~{minprot}" -t "~{prepared_transcripts}" -G "~{gencode}" -O transdecoder
    >>>
}

task TransdecoderPredict {
    input {
        File prepared_transcripts
        Array[File] transdecoder_wd
        File? blast_aligned_proteins
        RuntimeAttr? runtime_attr_override
        String workdir = sub(transdecoder_wd[0], "\/(?!.*\/).*", "")
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
        File bed = "mikado_prepared.fasta.transdecoder.bed"
    }

    command <<<
        set -euxo pipefail
        Transdecoder.Predict -t ~{prepared_transcripts} -O ~{workdir} ~{"--retain_blastp_hits " + blast_aligned_proteins}
    >>>

}