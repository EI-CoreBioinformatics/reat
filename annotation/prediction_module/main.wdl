version development

workflow ei_prediction {
	input {
		File reference_genome
		String species
		Directory augustus_config_path
		Array[File]? transcriptome_models
		Array[File]? homology_models
		File? intron_hints
		Int flank = 200
		Int kfold = 8
		Boolean optimise_augustus = false
		Boolean force_train = false

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

	call SelectAugustusTestAndTrain {
		input:
		models = def_training_models,
		reference = reference_genome,
		flank = flank
	}

	if (num_models > 1000) {
	# Feed all to Augustus in the various configurations
	# Generate training input for augustus (training + test sets)
	# Transform tranining set to GeneBank format

		call Species {
			input:
			base_config = augustus_config_path,
			species = species
		}

		if (!Species.existing_species && !force_train) {
			call etraining as base_training {
				input:
				models = SelectAugustusTestAndTrain.train,
				species = species,
				config_path = Species.config_path,
				with_utr = train_utr
			}

			call Augustus as AugustusTest {
				input:
				reference = reference_genome,
				models = SelectAugustusTestAndTrain.test,
				species = species,
				config_path = base_training.improved_config_path
			}

			if (optimise_augustus) {
				call OptimiseAugustus {
					input:
					num_fold = kfold,
					with_utr = train_utr,
					species = species,
					config_path = base_training.improved_config_path,
					models = SelectAugustusTestAndTrain.train
				}

				call etraining as train_after_optimise {
					input:
					models = SelectAugustusTestAndTrain.train,
					species = species,
					config_path = OptimiseAugustus.optimised_config_path,
					with_utr = train_utr
				}
			}
		}

		Directory final_augustus_config = select_first([train_after_optimise.improved_config_path, base_training.improved_config_path, Species.config_path])

		call Augustus as AugustusGold {
			input:
			reference = reference_genome,
			models = LengthChecker.gold,
			species = species,
			config_path = final_augustus_config
		}
	}


	output {
		File gold_models = LengthChecker.gold
		File silver_models = LengthChecker.silver
		File bronze_models = LengthChecker.bronze
		File training_models = def_training_models
		File? coding_quarry_predictions = CodingQuarry.predictions
		File? snap_predictions = SNAP.predictions
		File? glimmerhmm_predictions = GlimmerHMM.predictions
		File? augustus_predictions = AugustusGold.predictions
		Directory? augustus_config = final_augustus_config
	}
}

task OptimiseAugustus {
	input {
		Int num_fold
		String species
		Directory config_path
		File models
		Boolean with_utr
	}

	output {
		Directory optimised_config_path = "config"
	}

	command <<<
		set -euxo pipefail
		ln -s ~{config_path} config
		optimize_augustus.pl --AUGUSTUS_CONFIG_PATH=config --species=~{species} ~{models} --cpus=~{num_fold} --kfold=~{num_fold} ~{if with_utr then "--UTR=on" else ""}
	>>>
}

task Augustus {
	input {
		File reference
		Directory config_path
		String species
		File models
	}

	output {
		File predictions = "augustus.predictions.txt"
	}

	command <<<
		ln -s ~{config_path} config
		head ~{models}
		augustus --AUGUSTUS_CONFIG_PATH=~{config_path} --species=~{species} ~{reference} > augustus.predictions.txt
	>>>
}

task Species {
	input {
		Directory base_config
		String species
	}

	output {
		Directory config_path = "config"
		Boolean existing_species = read_boolean("existed")
	}

	command <<<
		set -euxo pipefail
		cp -r ~{base_config} config
		if [ ! -d "config/~{species}" ];
		then
			new_species.pl --species=~{species} --AUGUSTUS_CONFIG_PATH=config
			echo "false" > existed
		else
			echo "true" > existed
		fi
	>>>
}

task etraining {
	input {
		File models
		String species
		Boolean with_utr
		Directory config_path
	}

	output {
		Directory improved_config_path = "config"
	}

	command <<<
		ln -s ~{config_path} config
		etraining ~{models} --AUGUSTUS_CONFIG_PATH=config --species=~{species} --stopCodonExcludedFromCDS=~{with_utr}
	>>>
}

task SelectAugustusTestAndTrain {
	input {
		File models
		File reference
		Int flank = 200
	}

	output {
		File test = "test.gb"
		File train = "train.gb"
	}

	command <<<
		set -euxo pipefail
		generate_augustus_test_and_train ~{models}
		gff2gbSmallDNA.pl test.gff ~{reference} ~{flank} test.gb
		gff2gbSmallDNA.pl train.gff ~{reference} ~{flank} train.gb
	>>>
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
		set -euxo pipefail
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
		set -euxo pipefail
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
		set -euxo pipefail
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
		set -euxo pipefail
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
		set -euxo pipefail
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
		set -euxo pipefail
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
		set -euxo pipefail
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
		set -euxo pipefail
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
		set -euxo pipefail
		gffread -y proteins.faa -g ~{genome} ~{clustered_models}

		diamond makedb --db self -p 8 --in proteins.faa
		diamond blastp -p 8 -d self -q proteins.faa -f6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop > self.hits.tsv
		filter_self_hits -b self.hits.tsv -t ~{clustered_models} -c ~{classification}
		awk '$3 == "gene"' with_utr.gff|wc -l > num_models_utr.int
		awk '$3 == "gene"' without_utr.gff|wc -l > num_models_noutr.int
	>>>
}