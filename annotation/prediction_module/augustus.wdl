version development

import "structs.wdl"

workflow wf_augustus {
	input {
		IndexedReference reference_genome
		String species
		Directory augustus_config
		File extrinsic_config
		File? intron_hints
		IndexedBAM? expressed_exon_hints
		Array[File]? homology_models
		File? repeats_gff
		Boolean train_utr
		File? gold_models
		File? silver_models
		File? bronze_models
		File? all_models
		File hints_source_and_priority
		String run_id
	}

	call ReadSourceAndPriority {
		input:
		hints_source_and_priority = hints_source_and_priority
	}

		if (defined(intron_hints)) {
			call PrepareIntronHints {
				input:
				intron_gff = select_first([intron_hints]),
				gold_priority = ReadSourceAndPriority.gold_intron_hints_priority,
				gold_source = ReadSourceAndPriority.gold_intron_hints_source,
				silver_priority = ReadSourceAndPriority.silver_intron_hints_priority,
				silver_source = ReadSourceAndPriority.silver_intron_hints_source
			}
		}

		if (defined(expressed_exon_hints) && ReadSourceAndPriority.alignment_hints_source != '#') {
			IndexedBAM def_exon = select_first([expressed_exon_hints])
			if (!defined(def_exon.index)) {
				call IndexBAM {
					input:
					bam = def_exon.bam
				}
				IndexedBAM def_calc_index = object {bam: def_exon.bam, index: IndexBAM.index}
			}

			IndexedBAM def_indexed_bam = select_first([def_calc_index, def_exon])

			call PrepareExonHints {
				input:
				expression_bam = def_indexed_bam,
				dUTP = true,
				source = ReadSourceAndPriority.alignment_hints_source,
				priority = ReadSourceAndPriority.alignment_hints_priority
			}
		}

		if (defined(homology_models) && ReadSourceAndPriority.protein_hints_source != '#') {
			call PrepareProteinHints {
				input:
				aligned_proteins_gffs = select_first([homology_models]),
				priority = ReadSourceAndPriority.protein_hints_priority,
				source = ReadSourceAndPriority.protein_hints_source
			}
		}

		if (defined(gold_models) && ReadSourceAndPriority.gold_model_hints_source != '#') {
			call PrepareTranscriptHints as gold {
				input:
				transcripts = [select_first([gold_models])],
				category = "gold",
				source = ReadSourceAndPriority.gold_model_hints_source,
				priority = ReadSourceAndPriority.gold_model_hints_priority
			}
		}

		if (defined(silver_models) && ReadSourceAndPriority.silver_model_hints_source != '#') {
			call PrepareTranscriptHints as silver {
				input:
				transcripts = [select_first([silver_models])],
				category = "silver",
				source = ReadSourceAndPriority.silver_model_hints_source,
				priority = ReadSourceAndPriority.silver_model_hints_priority
			}
		}

		if (defined(bronze_models) && ReadSourceAndPriority.bronze_model_hints_source != '#') {
			call PrepareTranscriptHints as bronze {
				input:
				transcripts = [select_first([bronze_models])],
				category = "bronze",
				source = ReadSourceAndPriority.bronze_model_hints_source,
				priority = ReadSourceAndPriority.bronze_model_hints_priority
			}
		}

		if (defined(all_models) && ReadSourceAndPriority.all_model_hints_source != '#') {
			call PrepareTranscriptHints as all {
				input:
				transcripts = [select_first([all_models])],
				category = "all",
				source = ReadSourceAndPriority.all_model_hints_source,
				priority = ReadSourceAndPriority.all_model_hints_priority
			}
		}

		call cat as augustus_hints {
			input:
			files = select_all([gold.result, silver.result, bronze.result, all.result,
							   PrepareIntronHints.gold_intron_hints, PrepareIntronHints.silver_intron_hints,
							   PrepareProteinHints.protein_hints,
							   PrepareExonHints.expression_gff,
							   repeats_gff
							   ]
					),
			out_filename = "hints_" + run_id + ".gff"
		}

		call Augustus {
			input:
			reference = reference_genome,
			extrinsic_config = extrinsic_config,
			hints = augustus_hints.out,
			with_utr = train_utr,
			species = species,
			config_path = augustus_config,
			id = run_id
		}


	output {
		File hints = augustus_hints.out
		File predictions = Augustus.predictions
	}
}

task cat {
	input {
		Array[File]+ files
		String out_filename
	}

	output {
		File out = out_filename
	}

	command <<<
		cat ~{sep=' ' files} > ~{out_filename}
	>>>
}

task PrepareIntronHints {
	input {
		File intron_gff
		Int gold_priority = 6
		String gold_source
		Int silver_priority = 4
		String silver_source
	}

	output {
		File? gold_intron_hints = 'gold_junctions_S'+gold_source+'P'+gold_priority+'.gff'
		File? silver_intron_hints = 'silver_junctions_S'+silver_source+'P'+silver_priority+'.gff'
	}

	command <<<
		if [[ '#' != ~{gold_source} ]]; then
			cat ~{intron_gff} | awk -F "\t" '$6==1 {print $0";pri=~{gold_priority};"}' | \
			sed 's/\tportcullis\t/\tPortcullis_pass_gold_S~{gold_source}P~{gold_priority}\t/' > gold_junctions_S~{gold_source}P~{gold_priority}.gff
		fi

		if [[ '#' != ~{silver_source} ]]; then
			cat ~{intron_gff} | awk -F "\t" '$6<1 {print $0";pri=~{silver_priority};"}' | \
			sed 's/\tportcullis\t/\tPortcullis_pass_silver_S~{silver_source}P~{silver_priority}\t/' > silver_junctions_S~{silver_source}P~{silver_priority}.gff
		fi
	>>>
}

task ReadSourceAndPriority {
	input {
		File hints_source_and_priority
	}
	# NOTE: Check the generate_augustus_hint_parameters script if the outputs or the workflow requires changes!
	output {
		String 	gold_model_hints_source = read_string("gold.s")
		Int 	gold_model_hints_priority = read_int("gold.p")
		String 	silver_model_hints_source = read_string("silver.s")
		Int 	silver_model_hints_priority = read_int("silver.p")
		String 	bronze_model_hints_source = read_string("bronze.s")
		Int 	bronze_model_hints_priority = read_int("bronze.p")
		String 	all_model_hints_source = read_string("all.s")
		Int 	all_model_hints_priority = read_int("all.p")
		String 	alignment_hints_source = read_string("alignment.s")
		Int 	alignment_hints_priority = read_int("alignment.p")
		String 	repeat_hints_source = read_string("repeat.s")
		Int 	repeat_hints_priority = read_int("repeat.p")
		String 	protein_hints_source = read_string("protein.s")
		Int 	protein_hints_priority = read_int("protein.p")
		String 	gold_intron_hints_source = read_string("gold_intron.s")
		Int 	gold_intron_hints_priority = read_int("gold_intron.p")
		String 	silver_intron_hints_source = read_string("silver_intron.s")
		Int 	silver_intron_hints_priority = read_int("silver_intron.p")
	}

	command <<<
		generate_augustus_hint_parameters ~{hints_source_and_priority}
	>>>
}

task IndexBAM {
	input {
		File bam
	}

	output {
		File index = basename(bam)+'.csi'
	}

	command <<<
		ln -s ~{bam}
		samtools index -@ 4 -c ~{basename(bam)}
	>>>
}

task PrepareExonHints {
	input {
		IndexedBAM expression_bam
		String source
		Int priority
		Boolean dUTP
	}

	String bam_name = basename(expression_bam.bam, ".bam")
	String jid = basename(expression_bam.bam)

	output {
		File expression_gff = bam_name + '.exonhints.S'+source+'P'+priority+'.augustus.gff'
	}

	# Uses 2 cpus

	command <<<
		ln -s ~{expression_bam.bam}
		ln -s ~{expression_bam.index}

		if [ ~{dUTP} ];
		then
			strand='1+-,1-+,2++,2--'
		else
			strand='1++,1--,2+-,2-+'
		fi
		touch ~{bam_name}_Forward.wig ~{bam_name}_Reverse.wig ~{bam_name}.wig
		samtools view -H ~{basename(expression_bam.bam)} | grep '^\@SQ'| cut -f2,3 | sed -e 's/SN://g' -e 's/LN://g' > lengths.txt
		bam2wig.py -i ~{basename(expression_bam.bam)} -s lengths.txt -o ~{bam_name} --strand=$strand"

		cat ~{bam_name}.Forward.wig | \
		wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=~{source} --type=exonpart --radius=4.5 --pri=~{priority} --strand='+' | \
		sed \"s/\\tw2h\\t/\\tw2h_~{jid}\\t/\" > ~{bam_name}.Forward.exonhints.S~{source}P~{priority}.augustus.gff" &

		cat ~{bam_name}.Reverse.wig | \
		awk '/^variableStep/ {gsub(/-/,\"\");print;} !/^variableStep/ {print;}' \ |
		wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=~{source} --type=exonpart --radius=4.5 --pri=~{priority} --strand='-' | \
		sed \"s/\\tw2h\\t/\\tw2h_~{jid}\\t/\" > ~{bam_name}.Reverse.exonhints.S~{source}P~{priority}.augustus.gff" &

		wait

		cat ~{bam_name}.Forward.exonhints.S~{source}P~{priority}.augustus.gff ~{bam_name}.Reverse.exonhints.S~{source}P~{priority}.augustus.gff > ~{bam_name}.exonhints.S~{source}P~{priority}.augustus.gff
	>>>
}

task PrepareProteinHints {
	input {
		Array[File] aligned_proteins_gffs
		String source = "P"
		Int priority = 4
	}

	output {
		File protein_hints = 'aligned_proteins.S'+source+'P'+priority+'.gff'
	}

	command <<<
		for i in ~{sep=' ' aligned_proteins_gffs}; do
		name=$(echo ${basename i} | sed 's/.gff//')
		gffread $i | gff_to_aug_hints -P ~{priority} -S ~{source} -s ${name}.protein -t CDS >> aligned_proteins.S~{source}P~{priority}.gff;
		done
	>>>

}

task PrepareTranscriptHints {
	input {
		Array[File] transcripts
		String category
		String source = "E"
		Int priority = 4
	}

	output {
		File result = category+'.transcripts.S'+source+'P'+priority+'.gff'
	}

	command <<<
		for i in ~{sep=' ' transcripts}; do
		name=$(echo ${basename i} | sed 's/.gff//')
		cat $i | gff_to_aug_hints -P ~{priority} -S ~{source} -s ${name}.transcripts -t exon >> ~{category}.transcripts.S~{source}P~{priority}.gff;
		done
	>>>
}

task Augustus {
	input {
		IndexedReference reference
		File extrinsic_config
		Directory config_path
		Boolean with_utr
		String species
		File? hints
		String id
	}

	output {
		File predictions = "augustus_"+id+".predictions.gff"
	}

	command <<<
		ln -s ~{config_path} config
		augustus --AUGUSTUS_CONFIG_PATH=~{config_path} \
		--UTR=~{if with_utr then "ON" else "OFF"} --stopCodonExcludedFromCDS=true --genemodel=partial \
		--alternatives-from-evidence=true ~{'--hintsfile=' + hints} --noInFrameStop=true \
		--allow_hinted_splicesites=atac --errfile=run~{id}.log --extrinsicCfgFile=~{extrinsic_config} \
		--species=~{species} ~{reference.fasta} > augustus_~{id}.predictions.gff
	>>>
}
