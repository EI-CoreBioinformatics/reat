version 1.0

workflow ei_prediction {
	input {
		File reference_genome
		Array[File]? transcriptome_models
		Array[File]? homology_models

		File protein_validation_database
	}

	# Preprocess gene models
	if (defined(transcriptome_models)) {
		scatter (models in select_first([transcriptome_models])) {
			call PreprocessFiles as PreprocessTranscriptomic{
				input:
				models = models
			}
		}

		Array[File] processed_transcriptome = PreprocessTranscriptomic.out
	}

	if (defined(homology_models)) {
		scatter (models in select_first([homology_models])) {
			call PreprocessFiles as PreprocessHomology {
				input:
				models = models
			}
		}

		Array[File] processed_homology = PreprocessHomology.out
	}

	Array[File] all_models = flatten(select_all([processed_transcriptome, processed_homology]))

	# Generate protein files for the input models
	call GenerateModelProteins {
		input:
		genome = reference_genome,
		models = all_models
	}

	# TODO: Maybe split proteins for alignments?

	call IndexProteinsDatabase {
		input:
		db = protein_validation_database
	}

	call AlignProteins {
		input:
		proteins = GenerateModelProteins.proteins,
		db = IndexProteinsDatabase.diamond_index
	}

	# Check models for 'full-length'
	call LengthChecker {
		input:
		genome = reference_genome,
		models = all_models,
		protein_models = GenerateModelProteins.proteins,
		hits = AlignProteins.hits
	}

	call SelfBlastFilter {
		input:
		genome = reference_genome,
		clustered_models = LengthChecker.clustered_models,
		classification = LengthChecker.nr_classification
	}

	# If we have enough models we train with UTR, otherwise we use the 'extended' training set and train without UTR
	if (SelfBlastFilter.num_utr_models >= 1200) {
		File maybe_utr = SelfBlastFilter.non_redundant_models_with_utr
		Int maybe_num_utr = SelfBlastFilter.num_utr_models
		Boolean maybe_train_utr = true
	}

	if (SelfBlastFilter.num_utr_models < 1200) {
		File maybe_no_utr = SelfBlastFilter.non_redundant_models_without_utr
		Int maybe_num_no_utr = SelfBlastFilter.num_noutr_models
		Boolean maybe_train_no_utr = false
	}

	File def_training_models = select_first([maybe_utr, maybe_no_utr])
	Int num_models = select_first([maybe_num_utr, maybe_num_no_utr])
	Boolean train_utr = select_first([maybe_train_utr, maybe_train_no_utr])

	# Generate CodingQuarry predictions
	if (num_models > 1500) {
		call CodingQuarry {
			input:
			genome = reference_genome,
			transcripts = def_training_models
		}
	}

	# Generate SNAP predictions
	if (num_models > 1300) {
		call SNAP {
			input:
			genome = reference_genome,
			transcripts = def_training_models
		}
	}

	# Generate GlimmerHMM predictions
	if (num_models > 2000) {
		call GlimmerHMM {
			input:
			genome = reference_genome,
			transcripts = def_training_models
		}
	}

	# Feed all to Augustus in the various configurations


	output {
		File gold_models = LengthChecker.gold
		File silver_models = LengthChecker.silver
		File  bronze_models = LengthChecker.bronze
		File training_models = def_training_models
		File? coding_quarry_predictions = CodingQuarry.predictions
		File? snap_predictions = SNAP.predictions
		File? glimmerhmm_predictions = GlimmerHMM.predictions
	}
}

task CodingQuarry {
	input {
		File genome
		File transcripts
	}

	output {
		File predictions = "codingquarry.predictions.gff"
	}

	Int num_cpus = 8

	command <<<

		CodingQuarry -p ~{num_cpus} -f ~{genome} -t ~{transcripts}
		mv out/PredictedPass.gff3 codingquarry.predictions.gff
	>>>
}

task PreprocessFiles {
	input {
		File models
	}

	output {
		File out = basename(models) + ".clustered.gtf"
	}

	command <<<
		gffread -F --cluster-only --keep-genes ~{models} > ~{basename(models)}.clustered.gtf
	>>>
}

task SNAP {
	input {
		File genome
		File transcripts
	}

	output {
		File predictions = "snap.predictions.gff"
	}

	command <<<
		gff_to_zff -a ~{transcripts} > ~{transcripts}.ann
		fathom ~{transcripts}.ann ~{genome} -categorize 1000
		fathom uni.ann uni.dna -export 1000 -plus
		mkdir params
		cd params
		forge ../export.ann ../export.dna
		cd ..
		hmm-assembler.pl ~{genome} params > ~{genome}.hmm
		snap ~{genome}.hmm ~{genome} > snap.prediction.zff
		zff_to_gff -z snap.prediction.zff > snap.predictions.gff
	>>>
}

task GlimmerHMM {
	input {
		File genome
		File transcripts
	}

	output {
		File predictions = "glimmer.predictions.gff"
	}

	command <<<
		gff_to_glimmer -a ~{transcripts} > ~{transcripts}.cds
		trainGlimmerHMM ~{genome} ~{transcripts}.cds -d glimmer_training
		glimmhmm.pl $(which glimmerhmm) ~{genome} glimmer_training -g | gffread -g ~{genome} -vE --keep-genes -P > glimmer.predictions.gff
	>>>
}

task GenerateModelProteins {
	input {
		File genome
		Array[File] models
	}

	output {
		File proteins = "proteins.faa"
	}

	command <<<
		cat ~{sep=" " models} | gffread --stream -g ~{genome} -y proteins.faa
	>>>
}

task IndexProteinsDatabase {
	input {
		File db
	}

	output {
		File diamond_index = basename(db) + ".dmnd"
	}

	command <<<
		diamond makedb -d ~{basename(db)}.dmnd --in ~{db}
	>>>
}

task AlignProteins {
	input {
		File db
		File proteins
	}

	output {
		File hits = "diamond.hits.tsv"
	}

	command <<<
		diamond blastp -d ~{db} -q ~{proteins} -f6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop > diamond.hits.tsv
	>>>
}

task LengthChecker {
	input{
		File genome
		Array[File] models
		File protein_models
		File hits
	}

	output {
		File classification = "transcripts.classification.tsv"
		File nr_classification = "transcript.locus.classification.tsv"
		File clustered_models = "all_models.clustered.gff"
		File gold = "gold.gff"
		File silver = "silver.gff"
		File bronze = "bronze.gff"
	}

	command <<<
		gffread -g ~{genome} -F --cluster-only --keep-genes -P ~{sep=" " models} > all_models.clustered.gff
		classify_transcripts -b ~{hits} -t all_models.clustered.gff
	>>>
}

task SelfBlastFilter{
	input {
		File genome
		File clustered_models
		File classification
	}

	output {
		File non_redundant_models_with_utr = "with_utr.gff"
		File non_redundant_models_without_utr = "without_utr.gff"
		Int num_utr_models = read_int("num_models_utr.int")
		Int num_noutr_models = read_int("num_models_noutr.int")
	}

	command <<<
		gffread -y proteins.faa -g ~{genome} ~{clustered_models}
		diamond makedb --db self -p 8 --in proteins.faa
		diamond blastp -p 8 -d self -q proteins.faa -f6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop > self.hits.tsv
		filter_self_hits -b self.hits.tsv -t ~{clustered_models} -c ~{classification}
		awk '$3 == "gene"' with_utr.gff|wc -l > num_models_utr.int
		awk '$3 == "gene"' without_utr.gff|wc -l > num_models_noutr.int
	>>>
}