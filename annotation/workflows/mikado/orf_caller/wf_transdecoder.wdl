version 1.0

import "../../common/structs.wdl"
import "../../common/tasks.wdl"
import "../align_protein/wf_protein_aligner.wdl" as prt_aln

workflow wf_transdecoder {
    input {
        File prepared_transcripts
        String program = "blast"
        File? orf_proteins
        Boolean refine_start_codons = true
    }

# Prepare the protein alignment database
    if (defined(orf_proteins)) {
        File def_db = select_first([orf_proteins])
        call prt_aln.SanitiseProteinBlastDB {
            input:
                db = def_db
        }
        if (program == "blast") {
            call prt_aln.BlastIndex as BlastIndex {
                input:
                target = SanitiseProteinBlastDB.clean_db
            }
        if (program == "diamond") {
            call prt_aln.DiamondIndex as DiamondIndex {
                input:
                target = SanitiseProteinBlastDB.clean_db
            }
    }


# Step 1 Subdivide the prepared_transcripts
    call tasks.SplitSequences {
        input:
        sequences_file = prepared_transcripts,
        prefix = "transdecoder_split"
    }

# Step 2 Generate the longest_orfs and nucleotide frequency table (just call LongOrf in a scatter)
    scatter (chunk in SplitSequences.seq_files) {
        call TransdecoderLongOrf {
            input:
            prepared_transcripts = chunk
        }
    }

# Step 3 Select top longest CDS entries from the scatter, filter similar proteins (can use TDC exclude_similar), select top longest for the final set (can use TDC get_top_longest)
# Step 4 Train the markov model
    call Select_TrainMM_Score {
        input:
        cds_files = TransdecoderLongOrf.cds
    }

# Step 7 preparation
    if(refine_start_codons) {
        call Train_PWM {
            input:
            prepared_transcripts = prepared_transcripts,
            top_cds = Select_TrainMM_Score.top_cds
        }
    }

# Step 6.1 Scatter again to score all cds entries using "hexamer_scores_file"
    scatter (chunk in processed_chunk) {
        call Calculate_Scores {
            input:
            cds = chunk,
            scores = Select_TrainMM_Score.hexamer_scores_file
        }
    }

# Step 6.2 (optional) Run blastp on the chunks to filter
    if (orf_proteins) {
        scatter (pep_chunk in pep_chunks) {
            if (defined(orf_proteins)) {
                if (program == "blast") {
                    call prt_aln.BlastAlign {
                        input:
                        index = BlastIndex.db,
                        query = chunk,
                        blast_type = "blastp",
                        output_filename = "mikado_blast_orfcalling.txt",
                    }
                }
                if (program == "diamond") {
                    call prt_aln.DiamondAlign {
                        input:
                        query = TransdecoderLongOrf.pep,
                        output_filename = "mikado_diamond_orfcalling.txt",
                        blast_type = "blastp",
                        index = DiamondIndex.index
                    }
                }
            }
        }
    }

# Step 6.3 In the same scatter generate "best_candidates.gff3" per chunk
        # call selectBestORFs {}
# Step 7 (optional) Refine start codons
        # call startCodonRefinement {}
    # }





    # call TransdecoderLongOrf {
    #     input:
    #         prepared_transcripts = prepared_transcripts
    # }

   if (defined(orf_proteins)) {
       File def_db = select_first([orf_proteins])
        call prt_aln.SanitiseProteinBlastDB { # Check if this is an input for both homology + orf_calling or just the same task
            input:
                db = def_db
        }

        if (program == "blast") {
            call prt_aln.BlastIndex as BlastIndex {
                input:
                target = SanitiseProteinBlastDB.clean_db
            }

            call prt_aln.BlastAlign {
                input:
                query = TransdecoderLongOrf.pep,
                blast_type = "blastp",
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
                output_filename = "mikado_diamond_orfcalling.txt",
                blast_type = "blastp",
                index = DiamondIndex.index
            }
        }

        File maybe_blast_aligned_proteins = select_first([BlastAlign.out, DiamondAlign.out])
   }

    call TransdecoderPredict {
        input:
        prepared_transcripts = prepared_transcripts,
        transdecoder_wd = TransdecoderLongOrf.transdecoder_wd,
        blast_aligned_proteins = maybe_blast_aligned_proteins
    }

    output {
        File long_orfs = TransdecoderLongOrf.orfs
        File? clean_db = SanitiseProteinBlastDB.clean_db
        File final_orfs = TransdecoderPredict.bed
    }
}

task TransdecoderLongOrf {
    input {
        File prepared_transcripts
        Int minprot = 20 # Minimum protein length (TODO make parameter)
        String gencode = "Universal" # Genetic code (TODO Make parameter)
# Options for gencode: 
    # Universal
    # Tetrahymena
    # Acetabularia
    # Ciliate
    # Dasycladacean
    # Hexamita
    # Candida
    # Euplotid
    # SR1_Gracilibacteria
    # Pachysolen_tannophilus
    # Mesodinium
    # Peritrich
    # Mitochondrial-Vertebrates
    # Mitochondrial-Yeast
    # Mitochondrial-Invertebrates
    # Mitochondrial-Protozoan
    # Mitochondrial-Echinoderm
    # Mitochondrial-Ascidian
    # Mitochondrial-Flatworm
    # Mitochondrial-Chlorophycean
    # Mitochondrial-Trematode
    # Mitochondrial-Scenedesmus_obliquus
    # Mitochondrial-Thraustochytrium
    # Mitochondrial-Pterobranchia
# For more gencode info https://github.com/TransDecoder/TransDecoder/blob/master/PerlLib/Nuc_translator.pm
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
        File cds = "transdecoder/longest_orfs.cds"
        File gff3 = "transdecoder/longest_orfs.gff3"
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

task Select_TrainMM_Score {
    input {
        Array[File] cds_files
        File base_freqs
        Int topORFs_train = 50
        Int red_num = topORFs_train * 10
    }

    command <<<
    set -euxo pipefail
    # TODO Write some code to get the longest from an array of files
    get_top_longest_from_multiple_files.py ~{sep=" " cds_files} ~{red_num} > "cds.longest_~{red_num}"
    exclude_similar_proteins.pl "cds.longest_~{red_num}" > "cds.longest_~{red_num}.nr"
    get_top_longest_fasta_entries.pl "cds.longest_~{red_num}.nr" ~{topORFs_train} > "top_~{topORFs_train}_longest"

    seq_n_baseprobs_to_loglikelihood_vals.pl "top_~{topORFs_train}_longest" "~{base_freqs}" > "hexamer.scores"
    >>>

    output {
        File top_cds = "top_"+topORFs_train+"_longest"
        File hexamer_scores_file = "hexamer.scores"
    }
}

task Train_PWM {
    input {
        File prepared_transcripts
        File top_cds
    }

    output {
        Array[File] refinement_model = glob("refinement*")
    }

    command <<<
    set -euxo pipefail
    train_start_PWM.pl --transcripts ~{prepared_transcripts} --selected_orfs ~{top_cds} --out_prefix refinement
    >>>
}

task Calculate_Scores {
    input {
        File cds
        File scores
    }

    output {
        File scored_cds = cds+".scores"
    }

    command <<<
    set -euxo pipefail
    score_CDS_likelihood_all_6_frames.pl "~{cds}" "~{scores}" > "~{cds}.scores"
    >>>
}