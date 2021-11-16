version development

import "structs.wdl"
import "augustus.wdl"

workflow ei_prediction {
	input {
		IndexedReference reference_genome
		File extrinsic_config
		String species
		Directory augustus_config_path
		Boolean do_codingquarry = false
		Directory? codingquarry_training
		Boolean do_glimmer = false
		Directory? glimmer_training
		Boolean do_snap = false
		File? snap_training
		Boolean do_augustus = true
		Array[File]? transcriptome_models  # Classify and divide into Gold, Silver and Bronze
		Array[File]? homology_models  # Take as is
		File? HQ_protein_alignments
		File? LQ_protein_alignments
		File? HQ_assembly
		File? LQ_assembly
		File? intron_hints  # Separate into gold == 1.0 score and silver (the rest)
		IndexedBAM? expressed_exon_hints  # Transform into gff passing by bigwig
		File? repeats_gff  # These are passed through to augustus
		File? extra_training_models  # These models are taken as-is directly as results from the training model selection
		Int flank = 200
		Int kfold = 8
		Int chunk_size = 3000000
		Int overlap_size = 100000
		Boolean optimise_augustus = false
		Boolean force_train = false
		Array[File]? augustus_runs # File with SOURCE PRIORITY pairs defining the augustus configurations
		File protein_validation_database
		File EVM_weights
		String? mikado_utr_files
		File? mikado_config
		File? mikado_scoring
	}

	call SoftMaskGenome {
		input:
		genome = reference_genome,
		repeats_gff = repeats_gff
	}

	if (defined(repeats_gff)) {
		call PreprocessRepeats {
			input:
			gff = select_first([repeats_gff]),
			source = "repeat"
		}
	}

	if (defined(HQ_protein_alignments)) {
		call ChangeSource as hq_protein {
			input:
			gff = select_first([HQ_protein_alignments]),
			source = "hq_protein"
		}
	}
	if (defined(LQ_protein_alignments)) {
		call ChangeSource as lq_protein {
			input:
			gff = select_first([LQ_protein_alignments]),
			source = "lq_protein"
		}
	}
	if (defined(HQ_assembly)) {
		call ChangeSource as hq_assembly {
			input:
			gff = select_first([HQ_assembly]),
			source = "hq_assembly"
		}
	}
	if (defined(LQ_assembly)) {
		call ChangeSource as lq_assembly {
			input:
			gff = select_first([LQ_assembly]),
			source = "lq_assembly"
		}
	}

	if (defined(expressed_exon_hints)) {
		IndexedBAM def_exon = select_first([expressed_exon_hints])
		if (!defined(def_exon.index)) {
			call IndexBAM {
				input:
				bam = def_exon.bam
			}
			IndexedBAM def_calc_index = object {bam: def_exon.bam, index: IndexBAM.index}
		}

		IndexedBAM def_indexed_bam = select_first([def_calc_index, def_exon])

		call Bam2Hints {
			input:
			expression_bam = def_indexed_bam,
			dUTP = true
		}
	}


	if (! defined(reference_genome.index)) {
		call IndexGenome {
			input:
			genome = reference_genome
		}

		IndexedReference new_reference_genome = object {fasta: SoftMaskGenome.soft_masked_genome, index: IndexGenome.index}
	}

	IndexedReference def_reference_genome = select_first(
											[
											new_reference_genome,
											object {fasta: SoftMaskGenome.soft_masked_genome, index: reference_genome.index}
											])

	if (defined(homology_models)) {
		scatter (models in select_first([homology_models])) {
			call PreprocessFiles as PreprocessHomology {
				input:
				models = models
			}
		}
		Array[File] processed_homology = PreprocessHomology.out
	}

	# Preprocess gene models
	if (defined(transcriptome_models)) {
		scatter (models in select_first([transcriptome_models])) {
			call PreprocessFiles as PreprocessTranscriptomic {
				input:
				models = models
			}
		}

		Array[File] processed_models = flatten(select_all([PreprocessTranscriptomic.out, processed_homology]))
		# Generate protein files for the input models
		call GenerateModelProteins {
			input:
			genome = def_reference_genome,
			models = processed_models
		}

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
			genome = def_reference_genome,
			models = processed_models,
			protein_models = GenerateModelProteins.proteins,
			hits = AlignProteins.hits
		}

		call SelfBlastFilter {
			input:
			genome = def_reference_genome,
			clustered_models = LengthChecker.clustered_models,
			classification = LengthChecker.nr_classification,
			extra_training_models = extra_training_models
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
		Boolean train_utr_ = select_first([maybe_train_utr, maybe_train_no_utr])
	}

	Boolean train_utr = select_first([train_utr_, false])

	# Generate CodingQuarry predictions
	# Considers lowercase as masked by default
	if (num_models > 1500 && do_codingquarry) {
		call CodingQuarry {
			input:
			genome = def_reference_genome,
			transcripts = select_first([def_training_models]),
			species = species,
			codingquarry_training = codingquarry_training
		}

		call CodingQuarry as CodingQuarryFresh {
			input:
			genome = def_reference_genome,
			transcripts = select_first([def_training_models]),
			species = species,
			codingquarry_training = codingquarry_training,
			fresh_prediction = true
		}
	}

	# Generate SNAP predictions
	# With the -lcmask option treats lowercase as masked nts
	if (num_models > 1300 && do_snap) {
		call SNAP {
			input:
			genome = def_reference_genome,
			transcripts = select_first([def_training_models]),
			pretrained_hmm = snap_training
		}
	}

	# Generate GlimmerHMM predictions
	# Does not consider any masking
	if (num_models > 2000 && do_glimmer) {
		call GlimmerHMM {
			input:
			genome = object {fasta: SoftMaskGenome.hard_masked_genome, index: def_reference_genome.index},
			transcripts = select_first([def_training_models]),
			training_directory = glimmer_training
		}
	}

	# Augustus
	if (num_models > 1000 && do_augustus) {
		IndexedReference augustus_genome = object { fasta: SoftMaskGenome.unmasked_genome, index: def_reference_genome.index }

		call GenerateGenomeChunks {
			input:
			reference = augustus_genome,
		}
		# Feed all to Augustus in the various configurations
		# Generate training input for augustus (training + test sets)
		# Transform tranining set to GeneBank format

		call Species {
			input:
			base_config = augustus_config_path,
			species = species
		}

		if (!Species.existing_species || force_train) {
			call SelectAugustusTestAndTrain {
				input:
				models = select_first([def_training_models]),
				reference = def_reference_genome,
				flank = flank
			}

			call etraining as base_training {
				input:
				models = SelectAugustusTestAndTrain.train,
				species = species,
				config_path = Species.config_path,
				with_utr = train_utr
			}

			call augustus.wf_augustus as AugustusAbinitio {
				input:
				single_seqs = GenerateGenomeChunks.single_seqs,
				many_seqs = GenerateGenomeChunks.many_seqs,
				extrinsic_config = extrinsic_config,
				train_utr = train_utr,
				species = species,
				augustus_config = base_training.improved_config_path,
				repeats_gff = PreprocessRepeats.processed_gff,
				chunk_size = chunk_size,
				overlap_size = overlap_size,
				run_id = "_ABINITIO"
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

		if (defined(augustus_runs)) {
			Array[Int] counter = range(length(select_first([augustus_runs])))
			scatter (augustus_run_and_index in zip(select_first([augustus_runs]), counter)) {
				call augustus.wf_augustus as Augustus {
					input:
					single_seqs = GenerateGenomeChunks.single_seqs,
					many_seqs = GenerateGenomeChunks.many_seqs,
					species = species,
					augustus_config = final_augustus_config,
					extrinsic_config = extrinsic_config,
					intron_hints = intron_hints,
					expressed_exon_hints = Bam2Hints.expression_gff,
					homology_models = homology_models,
					repeats_gff = PreprocessRepeats.processed_gff,
					train_utr = train_utr,
					gold_models = LengthChecker.gold,
					silver_models = LengthChecker.silver,
					bronze_models = LengthChecker.bronze,
					all_models = LengthChecker.clustered_models,
					hq_assembly_models = hq_assembly.processed_gff,
					lq_assembly_models = lq_assembly.processed_gff,
					hq_protein_alignment_models = hq_protein.processed_gff,
					lq_protein_alignment_models = lq_protein.processed_gff,
					hints_source_and_priority = augustus_run_and_index.left,
					chunk_size = chunk_size,
					overlap_size = overlap_size,
					run_id = augustus_run_and_index.right + 1
				}
			}
			Array[File] def_augustus_predictions = select_all(Augustus.predictions)
		}
	}


	call EVM {
		input:
		genome = def_reference_genome.fasta,
		augustus_abinitio = AugustusAbinitio.predictions,
		augustus_predictions = def_augustus_predictions,
		snap_predictions = SNAP.predictions,
		glimmer_predictions = GlimmerHMM.predictions,
		codingquarry_predictions = CodingQuarry.predictions,
		codingquarry_fresh_predictions = CodingQuarryFresh.predictions,
		hq_protein_alignments = hq_protein.processed_gff,
		lq_protein_alignments = lq_protein.processed_gff,
		hq_assembly = hq_assembly.processed_gff,
		lq_assembly = lq_assembly.processed_gff,
		weights = EVM_weights,
		segment_size = 500000,
		overlap_size = 50000
	}

	scatter(emv_part in read_lines(EVM.evm_commands)) {
		call ExecuteEVMCommand {
			input:
			command_to_run = emv_part
		}
	}

	call CombineEVM {
		input:
		genome = def_reference_genome.fasta,
		dummy = ExecuteEVMCommand.done,
		partitions = EVM.partitions_list
	}

	# Configurable list of inputs from (gold, silver, bronze, all, hq/lq assemblies)
	call DefineMikadoUTRs {
		input:
			files_selection = select_first([mikado_utr_files, "gold silver"]),
			gold = LengthChecker.gold,
			silver = LengthChecker.silver,
			bronze = LengthChecker.bronze,
			all = LengthChecker.clustered_models,
			hq_assembly = hq_assembly.processed_gff,
			lq_assembly = lq_assembly.processed_gff,
#			files = select_all([LengthChecker.gold, LengthChecker.silver, LengthChecker.bronze, LengthChecker.clustered_models,
#							   hq_assembly.processed_gff, lq_assembly.processed_gff, hq_protein.processed_gff, lq_protein.processed_gff]),
			out_filename = "models_with_utrs.gff3"
	}

	call Mikado {
		input:
			reference = def_reference_genome,
			extra_config = mikado_config,
			prediction = CombineEVM.predictions,
			utrs = DefineMikadoUTRs.out,
			output_prefix = "mikado"
	}

	call MikadoPick {
		input:
			config_file = Mikado.mikado_config,
			scoring_file = mikado_scoring,
			mikado_db = Mikado.mikado_db,
			transcripts = Mikado.prepared_gtf,
			output_prefix = "mikado"
	}

	call MikadoSummaryStats {
		input:
			stats = [MikadoPick.stats],
			output_prefix = "mikado"
	}

	output {
		File? gold_models = LengthChecker.gold
		File? silver_models = LengthChecker.silver
		File? bronze_models = LengthChecker.bronze
		File? training_models = def_training_models
		File? coding_quarry_predictions = CodingQuarry.predictions
		File? snap_predictions = SNAP.predictions
		File? glimmerhmm_predictions = GlimmerHMM.predictions
		Array[File]? augustus_predictions = def_augustus_predictions
		File? augustus_abinitio_predictions = AugustusAbinitio.predictions
		File evm_predictions = CombineEVM.predictions
		Directory? augustus_config = final_augustus_config
	}
}

task CombineEVM {
	input {
		File genome
		Array[String] dummy
		File partitions
	}

	output {
		File predictions = "evm.out.gff3"
	}

	command <<<
		$EVM_HOME/EvmUtils/recombine_EVM_partial_outputs.pl --partitions ~{partitions} --output_file_name evm.out
		cat $($EVM_HOME/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions ~{partitions} --output evm.out --genome ~{genome} | awk -F',' '{printf $2"/evm.out.gff3"} END{print ""}') > evm.out.gff3
	>>>
}

task ExecuteEVMCommand {
	input {
		String command_to_run
		RuntimeAttr? resources
	}

	Int cpus = 1
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 8,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([resources, default_attr])
    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

	runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }

	output {
		String done = "yes"
	}

	command <<<
		~{command_to_run}
	>>>
}

task SoftMaskGenome {
	input {
		IndexedReference genome
		File? repeats_gff
	}

	output {
		File soft_masked_genome = sub(basename(genome.fasta), "\\.(fasta|fa)$", ".softmasked.fa")
		File hard_masked_genome = sub(basename(genome.fasta), "\\.(fasta|fa)$", ".hardmasked.fa")
		File unmasked_genome = sub(basename(genome.fasta), "\\.(fasta|fa)$", ".unmasked.fa")
	}

	command <<<
cat ~{genome.fasta} | python3 -c "
import sys;
for line in sys.stdin:
    if line.startswith('>'):
        print(line, end='')
        continue
    print(line.upper(), end='')" > ~{sub(basename(genome.fasta), "\\.(fasta|fa)", ".unmasked.fa")}
rep_file=~{repeats_gff}
if [[ "~{repeats_gff}" == "" ]];
then
    touch repeats.gff
    rep_file="repeats.gff"
fi
bedtools maskfasta -soft -fi ~{genome.fasta} -bed <(gffread --bed $rep_file) -fo ~{sub(basename(genome.fasta), "\\.(fasta|fa)", ".softmasked.fa")}
bedtools maskfasta -mc 'N' -fi ~{genome.fasta} -bed <(gffread --bed $rep_file) -fo ~{sub(basename(genome.fasta), "\\.(fasta|fa)", ".hardmasked.fa")}
	>>>
}

task IndexBAM {
	input {
		File bam
	}

	output {
		File index = basename(bam)+'.bai'
	}

	command <<<
		ln -s ~{bam}
		samtools index -@ 4 ~{basename(bam)}
	>>>
}

task PreprocessRepeats {
	input {
		File gff
		String source
	}

	output {
		File processed_gff = "repeats.gff"
	}

	command <<<
	awk 'BEGIN{OFS="\t"} $3=="match" {print $1, "repmask", "nonexonpart", $4, $5, $6, $7, $8, "src=RM;pri=0"}' ~{gff} > repeats.gff
	>>>
}
task ChangeSource {
	input {
		File gff
		String source
	}

	output {
		File processed_gff = sub(basename(gff), "\\.(gff|gtf)" , "") + ".post.gff"
	}

	command <<<
		cat ~{gff} | gffread --stream -vE | awk -v 'OFS=\t' '
		$3=="exon" {print $1, "~{source}", "exon", $4, $5, $6, $7, $8, $9";src=generic_source;pri=0"}
		$3=="CDS" {print $1, "~{source}", "CDS", $4, $5, $6, $7, $8, $9";src=generic_source;pri=0"}
		($3 != "exon" && $3 != "CDS") {print $1, "~{source}", $3, $4, $5, $6, $7, $8, $9}
		' > ~{sub(basename(gff), "\\.(gff|gtf)" , "") + ".post.gff"}
	>>>
}

task Bam2Hints {
	input {
		IndexedBAM expression_bam
		Boolean? dUTP
	}

	String bam_name = basename(expression_bam.bam, ".bam")
	String jid = basename(expression_bam.bam)

	output {
		File expression_gff = bam_name + '.exonhints.gff'
	}

	# Uses 2 cpus

	command <<<
		ln -s ~{expression_bam.bam}
		ln -s ~{expression_bam.index}

		touch ~{bam_name}.unstranded.exonhints.augustus.gff \
		~{bam_name}.Forward.exonhints.augustus.gff \
		~{bam_name}.Reverse.exonhints.augustus.gff

		if [[ "~{dUTP}" != "" ]];
		then
			if [ ~{dUTP} ];
			then
				strand='1+-,1-+,2++,2--'
			else
				strand='1++,1--,2+-,2-+'
			fi
			touch ~{bam_name}.Forward.wig ~{bam_name}.Reverse.wig ~{bam_name}.wig
			samtools view -H ~{basename(expression_bam.bam)} | grep '^\@SQ'| cut -f2,3 | sed -e 's/SN://g' -e 's/LN://g' > lengths.txt
			bam2wig -i ~{basename(expression_bam.bam)} -s lengths.txt -o ~{bam_name} --strand="$strand"

			cat ~{bam_name}.Forward.wig | \
			wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=generic_source --type=exonpart --radius=4.5 --pri=0 --strand='+' | \
			sed 's/\\tw2h\\t/\\tw2h_~{jid}\\t/' > ~{bam_name}.Forward.exonhints.augustus.gff &

			cat ~{bam_name}.Reverse.wig | \
			awk '/^variableStep/ {gsub(/-/,\"\");print;} !/^variableStep/ {print;}' \ |
			wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=generic_source --type=exonpart --radius=4.5 --pri=0 --strand='-' | \
			sed 's/\\tw2h\\t/\\tw2h_~{jid}\\t/' > ~{bam_name}.Reverse.exonhints.augustus.gff &
			wait
		else
			bam2wig ~{basename(expression_bam.bam)} | \
			wig2hints.pl --width=10 --margin=10 --minthresh=2 \
			--minscore=4 --prune=0.1 --src=generic_source --type=exonpart \
			--radius=4.5 --pri=0 --strand="." > ~{bam_name}.unstranded.exonhints.augustus.gff
		fi

		cat ~{bam_name}.unstranded.exonhints.augustus.gff  \
		~{bam_name}.Forward.exonhints.augustus.gff \
		~{bam_name}.Reverse.exonhints.augustus.gff > ~{bam_name}.exonhints.gff
	>>>
}

task IndexGenome {
	input {
		IndexedReference genome
	}

	output {
		File index = basename(genome.fasta) + ".fai"
	}

	command <<<
		ln -s ~{genome.fasta}
		samtools faidx ~{basename(genome.fasta)}
	>>>
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
		if [ ! -d "~{base_config}/species/~{species}" ];
		then
			mkdir -p config/species
			cp -r ~{base_config}/species/generic config/species/generic
			new_species.pl --species=~{species} --AUGUSTUS_CONFIG_PATH=config
			sed '/^UTR/s/.*/UTR\ton/' config/species/~{species}/~{species}_parameters.cfg > config/species/~{species}/~{species}_parameters.cfg.utr
			mv config/species/~{species}/~{species}_parameters.cfg.utr config/species/~{species}/~{species}_parameters.cfg
			echo "false" > existed
		else
			cp -r ~{base_config}/species/~{species} config/species/~{species}
			echo "true" > existed
		fi
		cp -R ~{base_config}/{cgp,extrinsic,model,profile} config/
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
		etraining --AUGUSTUS_CONFIG_PATH=config --species=~{species} ~{models}
	>>>
}

task SelectAugustusTestAndTrain {
	input {
		File models
		IndexedReference reference
		Int flank = 200
	}

	output {
		File test = "test.gb"
		File train = "train.gb"
	}

	command <<<
		set -euxo pipefail
		generate_augustus_test_and_train ~{models}
		gff2gbSmallDNA.pl test.gff ~{reference.fasta} ~{flank} test.gb
		gff2gbSmallDNA.pl train.gff ~{reference.fasta} ~{flank} train.gb
	>>>
}

task CodingQuarry {
	input {
		IndexedReference genome
		Boolean fresh_prediction = false
		File transcripts
		Directory? codingquarry_training
		String? species
	}

	output {
		File predictions = "codingquarry.predictions.gff"
	}

	Int num_cpus = 8

	command <<<
		set -euxo pipefail
		if [ "~{codingquarry_training}" == "" ]; then
			CodingQuarry -p ~{num_cpus} -f ~{genome.fasta} -a ~{transcripts} ~{if fresh_prediction then "-n" else ""}
		else
			export QUARRY_PATH=~{codingquarry_training}
			CodingQuarry -p ~{num_cpus} -f ~{genome.fasta} -a ~{transcripts} -s ~{species} ~{if fresh_prediction then "-n" else ""}
		fi
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

task GenerateGenomeChunks {
	input {
		IndexedReference reference
	}

	output {
		Array[File]? single_seqs = glob("single_seqs_*.fa")
		Array[File]? many_seqs = glob("many_seqs_*.fa")
	}

	command <<<
		partition_genome --reference ~{reference.fasta}
	>>>
}

task SNAP {
	input {
		IndexedReference genome
		File transcripts
		File? pretrained_hmm
	}

	output {
		File predictions = "snap.predictions.gff"
		File raw = "snap.predictions.zff"
	}

	command <<<
		set -euxo pipefail
		if [ "~{pretrained_hmm}" == "" ]; then
			gff_to_zff -a ~{transcripts} > ~{transcripts}.ann
			fathom ~{transcripts}.ann ~{genome.fasta} -categorize 1000
			fathom uni.ann uni.dna -export 1000 -plus
			mkdir params
			cd params
			forge ../export.ann ../export.dna
			cd ..
			hmm-assembler.pl ~{genome.fasta} params > ~{genome.fasta}.hmm
			snap ~{genome.fasta}.hmm ~{genome.fasta} > snap.predictions.zff
		else
			snap ~{pretrained_hmm} ~{genome.fasta} > snap.predictions.zff
		fi
		zff_to_gff -z snap.predictions.zff > snap.predictions.gff
	>>>
}

task GlimmerHMM {
	input {
		IndexedReference genome
		File transcripts
		Directory? training_directory
	}

	output {
		File predictions = "glimmer.predictions.gff"
		File raw = "glimmer.raw.gtf"
	}

	command <<<
		set -euxo pipefail
		ln -s ~{genome.fasta}
		ln -s ~{genome.index}
		if [ "~{training_directory}" == "" ]; then
			gff_to_glimmer -a ~{transcripts} > ~{transcripts}.cds
			trainGlimmerHMM ~{basename(genome.fasta)} ~{transcripts}.cds -d glimmer_training
			glimmhmm.pl $(which glimmerhmm) ~{basename(genome.fasta)} glimmer_training -g | tee glimmer.raw.gtf | gffread -g ~{basename(genome.fasta)} -vE --keep-genes -P > glimmer.predictions.gff
		else
			glimmhmm.pl $(which glimmerhmm) ~{basename(genome.fasta)} ~{training_directory} -g | tee glimmer.raw.gtf | gffread -g ~{basename(genome.fasta)} -vE --keep-genes -P > glimmer.predictions.gff
		fi

	>>>
}

task GenerateModelProteins {
	input {
		IndexedReference genome
		Array[File] models
	}

	output {
		File proteins = "proteins.faa"
	}

	command <<<
		set -euxo pipefail
		ln -s ~{genome.fasta}
		ln -s ~{genome.index}
		cat ~{sep=" " models} | gffread --stream -g ~{basename(genome.fasta)} -y proteins.faa
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
	input {
		IndexedReference genome
		Array[File] models
		File protein_models
		File hits
		Float? min_pct_cds_fraction
		Int? max_tp_utr_complete
		Int? max_tp_utr
		Int? min_tp_utr
		Int? max_fp_utr_complete
		Int? max_fp_utr
		Int? min_fp_utr
		Int? query_start_hard_filter_distance
		Int? query_start_score
		Int? query_start_scoring_distance
		Int? query_end_hard_filter_distance
		Int? query_end_score
		Int? query_end_scoring_distance
		Int? target_start_hard_filter_distance
		Int? target_start_score
		Int? target_start_scoring_distance
		Int? target_end_hard_filter_distance
		Int? target_end_score
		Int? target_end_scoring_distance
		Int? min_query_coverage_hard_filter
		Int? min_query_coverage_score
		Int? min_query_coverage_scoring_percentage
		Int? min_target_coverage_hard_filter
		Int? min_target_coverage_score
		Int? min_target_coverage_scoring_percentage
		Int? max_single_gap_hard_filter
		Int? max_single_gap_score
		Int? max_single_gap_scoring_length
		RuntimeAttr? resources
	}

	Int cpus = 1
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 8,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([resources, default_attr])
    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

	runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
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
		ln -s ~{genome.fasta}
		ln -s ~{genome.index}
		gffread -g ~{basename(genome.fasta)} -F --cluster-only --keep-genes -P ~{sep=" " models} > all_models.clustered.gff
		classify_transcripts ~{"--min_pct_cds_fraction " + min_pct_cds_fraction} ~{"--max_tp_utr_complete " + max_tp_utr_complete} ~{"--max_tp_utr " + max_tp_utr} ~{"--min_tp_utr " + min_tp_utr} ~{"--max_fp_utr_complete " + max_fp_utr_complete} ~{"--max_fp_utr " + max_fp_utr} ~{"--min_fp_utr " + min_fp_utr} \
		-b ~{hits} ~{"--query_start_hard_filter_distance " + query_start_hard_filter_distance} ~{"--query_start_score " + query_start_score} ~{"--query_start_scoring_distance " + query_start_scoring_distance} ~{"--query_end_hard_filter_distance " + query_end_hard_filter_distance} ~{"--query_end_score " + query_end_score} ~{"--query_end_scoring_distance " + query_end_scoring_distance} \
		-t all_models.clustered.gff ~{"--target_start_hard_filter_distance " + target_start_hard_filter_distance} ~{"--target_start_score " + target_start_score} ~{"--target_start_scoring_distance " + target_start_scoring_distance} ~{"--target_end_hard_filter_distance " + target_end_hard_filter_distance} ~{"--target_end_score " + target_end_score} ~{"--target_end_scoring_distance " + target_end_scoring_distance} ~{"--min_query_coverage_hard_filter " + min_query_coverage_hard_filter} ~{"--min_query_coverage_score " + min_query_coverage_score} ~{"--min_query_coverage_scoring_percentage " + min_query_coverage_scoring_percentage} ~{"--min_target_coverage_hard_filter " + min_target_coverage_hard_filter} ~{"--min_target_coverage_score " + min_target_coverage_score} ~{"--min_target_coverage_scoring_percentage " + min_target_coverage_scoring_percentage} ~{"--max_single_gap_hard_filter " + max_single_gap_hard_filter} ~{"--max_single_gap_score " + max_single_gap_score} ~{"--max_single_gap_scoring_length " + max_single_gap_scoring_length}
	>>>
}

task SelfBlastFilter {
	input {
		IndexedReference genome
		File clustered_models
		File classification
		Int? top_n
		Int? identity
		Int? coverage
		File? extra_training_models
		RuntimeAttr? resources
	}

	Int cpus = 1
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 8,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([resources, default_attr])
    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

	runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }

	output {
		File non_redundant_models_with_utr = "with_utr.extra.gff"
		File non_redundant_models_without_utr = "without_utr.extra.gff"
		Int num_utr_models = read_int("num_models_utr.int")
		Int num_noutr_models = read_int("num_models_noutr.int")
	}

	command <<<
		set -euxo pipefail
		ln -s ~{genome.fasta}
		ln -s ~{genome.index}
		gffread -y proteins.faa -g ~{basename(genome.fasta)} ~{clustered_models}

		diamond makedb --db self -p 8 --in proteins.faa
		diamond blastp -p 8 -d self -q proteins.faa -f6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop > self.hits.tsv
		filter_self_hits -b self.hits.tsv -t ~{clustered_models} -c ~{classification} ~{"--top_n " + top_n} ~{"--max_identity " + identity} ~{"--max_coverage " + coverage}
		gffread --keep-genes --keep-exon-attrs -F with_utr.gff ~{extra_training_models} > with_utr.extra.gff
		gffread --keep-genes --keep-exon-attrs -F without_utr.gff ~{extra_training_models} > without_utr.extra.gff
		awk '$3 == "gene"' with_utr.extra.gff|wc -l > num_models_utr.int
		awk '$3 == "gene"' without_utr.extra.gff|wc -l > num_models_noutr.int
	>>>
}

task EVM {
	input {
		File genome
		Array[File]? augustus_predictions
		File? augustus_abinitio
		File? snap_predictions
		File? glimmer_predictions
		File? codingquarry_predictions
		File? codingquarry_fresh_predictions
		File? hq_protein_alignments
		File? lq_protein_alignments
		File? hq_assembly
		File? lq_assembly
		Int segment_size
		Int overlap_size
		File weights
		Array[File]? homology_models
		Array[File]? transcriptome_models
	}

	output {
		File evm_commands = "commands.list"
		File partitions_list = "partitions_list.out"
	}

	command <<<
		cat ~{if defined(hq_protein_alignments) then hq_protein_alignments else "/dev/null"} | sed 's/CDS/nucleotide_to_protein_match/g' | sed 's/ID=.*;Parent=/ID=/g' >> protein_alignments.gff
		cat ~{if defined(lq_protein_alignments) then hq_protein_alignments else "/dev/null"} | sed 's/CDS/nucleotide_to_protein_match/g' | sed 's/ID=.*;Parent=/ID=/g' >> protein_alignments.gff

		cat ~{if defined(hq_assembly) then hq_assembly else "/dev/null"} >> transcript_alignments.gff
		cat ~{if defined(lq_assembly) then lq_assembly else "/dev/null"} >> transcript_alignments.gff

		if [ "~{if defined(homology_models) then length(select_first([homology_models])) else 0}" != "0" ];
		then
			cat ~{sep=" " homology_models} | awk '$2="homology_models"' >> predictions.gff
		fi

		if [ "~{if defined(transcriptome_models) then length(select_first([transcriptome_models])) else 0}" != "0" ];
		then
			cat ~{sep=" " transcriptome_models} | awk '$2="transcriptome_models"' >> predictions.gff
		fi

		transcript_alignments_param=''
		protein_alignments_param=''
		if [ -s protein_alignments.gff ];
		then
			protein_alignments_param='--protein_alignments '$(realpath protein_alignments.gff)
		fi

		if [ -s transcript_alignments.gff ];
		then
			transcript_alignments_param='--transcript_alignments '$(realpath transcript_alignments.gff)
		fi

		if [ "~{snap_predictions}" != "" ]; then
			cat ~{snap_predictions} | gff_to_evm snap >> predictions.gff
		fi
		if [ "~{glimmer_predictions}" != "" ]; then
			cat ~{glimmer_predictions} | gff_to_evm glimmer >> predictions.gff
		fi

		if [ "~{codingquarry_predictions}" != "" ]; then
			cat ~{codingquarry_predictions} | gff_to_evm codingquarry >> predictions.gff
		fi
		if [ "~{codingquarry_fresh_predictions}" != "" ]; then
			cat ~{codingquarry_fresh_predictions} | gff_to_evm codingquarry >> predictions.gff
		fi

		if [ "~{augustus_abinitio}" != "" ]; then
			cat ~{augustus_abinitio} | gff_to_evm augustus | awk -v OFS="\t" '($2="AUGUSTUS_RUN_ABINITIO") && NF>8' >> predictions.gff
		fi

		augustus_run=0
		for i in ~{sep=" " augustus_predictions}; do
			((augustus_run++))
			cat $i | gff_to_evm augustus | awk -v OFS="\t" -v RUN=${augustus_run} '($2="AUGUSTUS_RUN"RUN) && NF>8' >> predictions.gff
		done

		predictions_fp=$(realpath predictions.gff)
		$EVM_HOME/EvmUtils/partition_EVM_inputs.pl --genome ~{genome} \
		--gene_predictions ${predictions_fp} ${protein_alignments_param} ${transcript_alignments_param} \
		--segmentSize ~{segment_size} \
		--overlapSize ~{overlap_size} \
		--partition_listing partitions_list.out

		$EVM_HOME/EvmUtils/write_EVM_commands.pl -S --debug --genome ~{genome} \
		--gene_predictions ${predictions_fp} ${protein_alignments_param} ${transcript_alignments_param} \
		--weights ~{weights} \
		--search_long_introns 5000 \
		--output_file_name evm.out \
		--partitions partitions_list.out > commands.list
	>>>
}

task DefineMikadoUTRs {
	input {
		String files_selection
		File? gold
		File? silver
		File? bronze
		File? all
		File? hq_assembly
		File? lq_assembly
		String out_filename
	}

	output {
		File out = out_filename
	}

	command <<<
		combine_utr_hint_files ~{"--gold " + gold} ~{"--silver " + silver} ~{"--bronze " + bronze} ~{"--all " + all} ~{"--hq_assembly " + hq_assembly} ~{"--lq_assembly " + lq_assembly} --selection ~{files_selection} > ~{out_filename}
	>>>
}


task Mikado {
	input {
		IndexedReference reference
		File? extra_config
		Int min_cdna_length = 100
		Int max_intron_length = 1000000
		File prediction
		File utrs
		File? junctions
		String output_prefix
		RuntimeAttr? resources
	}

	output {
		File mikado_config = output_prefix+"-mikado.yaml"
		File prepared_fasta = output_prefix+"-mikado/mikado_prepared.fasta"
		File prepared_gtf = output_prefix+"-mikado/mikado_prepared.gtf"
		File mikado_db = output_prefix+"-mikado/mikado.db"
	}

	Int cpus = 8
	RuntimeAttr default_attr = object {
								   cpu_cores: "~{cpus}",
								   mem_gb: 16,
								   max_retries: 1,
								   queue: ""
							   }

	RuntimeAttr runtime_attr = select_first([resources, default_attr])

	Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

	command <<<
		set -euxo pipefail
		export TMPDIR=/tmp
		# Create the lists file
		bname=$(basename ~{prediction})
		label=${bname%%.*}
		echo -e "~{prediction}\t${label}\tTrue\t0\tTrue" >> list.txt

		if [ "" != "~{utrs}" ]
		then
			label="UTRs"
			echo -e "~{utrs}\t${label}\tTrue\t0\tFalse" >> list.txt
		fi

		# mikado configure
		mikado configure --only-reference-update --scoring plant.yaml --copy-scoring plant.yaml \
		--max-intron-length ~{max_intron_length} --minimum-cdna-length ~{min_cdna_length} \
		--reference ~{reference.fasta} --list=list.txt original-mikado.yaml

		yaml-merge -s original-mikado.yaml ~{"-m " + extra_config} -o ~{output_prefix}-mikado.yaml

		# mikado prepare
		mikado prepare --procs=~{task_cpus} --json-conf ~{output_prefix}-mikado.yaml -od ~{output_prefix}-mikado

		gffread --nc -T -o mikado_prepared.nc.gtf -w mikado_prepared.nc.fasta -g ~{reference.fasta} ~{output_prefix}-mikado/mikado_prepared.gtf

		if [ -s mikado_prepared.nc.fasta ];
		then
			prodigal -g 1 -f gff -i mikado_prepared.nc.fasta -o mikado_prepared.cds.gff
		else
			touch mikado_prepared.cds.gff
		fi

		# mikado serialise
		mikado serialise --force ~{"--junctions " + junctions} \
		--transcripts ~{output_prefix}-mikado/mikado_prepared.fasta \
		--orfs mikado_prepared.cds.gff \
		--json-conf=~{output_prefix}-mikado.yaml --start-method=spawn -od ~{output_prefix}-mikado --procs=~{task_cpus}

	>>>

	runtime {
		cpu: task_cpus
		memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
		queue: select_first([runtime_attr.queue, default_attr.queue])
	}

}

task MikadoPick {
	input {
		File? config_file
		File? extra_config
		File? scoring_file
		File transcripts
		File mikado_db
		String output_prefix
		RuntimeAttr? resources
	}

	Int cpus = 8
	RuntimeAttr default_attr = object {
								   cpu_cores: "~{cpus}",
								   mem_gb: 16,
								   max_retries: 1,
								   queue: ""
							   }

	RuntimeAttr runtime_attr = select_first([resources, default_attr])

	Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

	output {
		File index_log  = output_prefix + "-index_loci.log"
		File loci_index = output_prefix + ".loci.gff3.midx"
		File loci       = output_prefix + ".loci.gff3"
		File scores     = output_prefix + ".loci.scores.tsv"
		File metrics    = output_prefix + ".loci.metrics.tsv"
		File stats      = output_prefix + ".loci.gff3.stats"
		File subloci    = output_prefix + ".subloci.gff3"
		File monoloci   = output_prefix + ".monoloci.gff3"
	}

	command <<<
		set -euxo pipefail
		export TMPDIR=/tmp
		yaml-merge -s ~{config_file} ~{"-m " + extra_config} -o pick_config.yaml
		mikado pick --source Mikado_~{output_prefix} --procs=~{task_cpus} ~{"--scoring-file" + scoring_file} \
		--start-method=spawn --json-conf=pick_config.yaml \
		--loci-out ~{output_prefix}.loci.gff3 -lv INFO ~{"-db " + mikado_db} \
		--subloci-out ~{output_prefix}.subloci.gff3 --monoloci-out ~{output_prefix}.monoloci.gff3 \
		~{transcripts}
		mikado compare -r ~{output_prefix}.loci.gff3 -l ~{output_prefix}-index_loci.log --index
		mikado util stats  ~{output_prefix}.loci.gff3 ~{output_prefix}.loci.gff3.stats
	>>>

	runtime {
		cpu: task_cpus
		memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
		queue: select_first([runtime_attr.queue, default_attr.queue])
	}

}

task MikadoSummaryStats {
	input {
		Array[File] stats
		String output_prefix
	}

	output {
		File summary = output_prefix + ".summary.stats.tsv"
	}

	command <<<
		mikado_summary_stats ~{sep=" " stats} > ~{output_prefix}.summary.stats.tsv
	>>>
}
