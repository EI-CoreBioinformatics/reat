version 1.0
import "common/rt_struct.wdl"

struct GenomeAnnotation {
    File genome
    File annotation_gff
}

struct GenomeProteins {
    File genome
    File annotation_gff
    File protein_sequences
}

workflow ei_homology {
    input {
        Array[GenomeAnnotation]? annotations
        Array[File]? protein_sequence_files
        File genome_to_annotate
        File mikado_scoring
        File mikado_config
        String spaln_extra_parameters
		Boolean do_spaln = false
		Boolean do_miniprot = false
        String species
        Int codon_table
        String output_prefix = "xspecies"
        RuntimeAttr? index_attr
        RuntimeAttr? score_attr
        RuntimeAttr? aln_attr
        RuntimeAttr? mikado_attr
        RuntimeAttr? task_attr
    }

    call IndexGenomeSpaln {
        input:
        genome = genome_to_annotate,
        runtime_attr_override = index_attr
    }
    call IndexGenomeMiniprot {
        input:
        genome = genome_to_annotate,
        runtime_attr_override = index_attr
    }

    if (defined(protein_sequence_files)) {
        scatter (protein_sequence_file in select_first([protein_sequence_files])) {

            call AlignProteinSequencesSpaln {
                input:
                genome_index = IndexGenomeSpaln.genome_index,
                proteins = protein_sequence_file,
                species = species,
                spaln_extra_parameters = spaln_extra_parameters,
                out_prefix = sub(basename(protein_sequence_file), "\.[^/.]+$", ""),
                runtime_attr_override = aln_attr
            }

            call AlignProteinSequencesMiniprot {
                input:
                genome_index = IndexGenomeMiniprot.genome_index,
                proteins = protein_sequence_file,
                out_prefix = sub(basename(protein_sequence_file), "\.[^/.]+$", ""),
                runtime_attr_override = aln_attr
            }

            call PrepareProteinAlignments as PrepareProteinAlignmentsSpaln {
                input:
                raw_alignments = AlignProteinSequencesSpaln.raw_alignments,
                aln_prefix = sub(basename(AlignProteinSequencesSpaln.raw_alignments), "\\.alignment\\.gff$", ""),
                genome_to_annotate = genome_to_annotate,
                runtime_attr_override = task_attr
            }

            call PrepareProteinAlignments as PrepareProteinAlignmentsMiniprot {
                input:
                raw_alignments = AlignProteinSequencesMiniprot.raw_alignments,
                aln_prefix = sub(basename(AlignProteinSequencesMiniprot.raw_alignments), "\\.alignment\\.gff$", ""),
                genome_to_annotate = genome_to_annotate,
                runtime_attr_override = task_attr
            }
        }
    }

    if (defined(annotations)) {
        scatter (annotation in select_first([annotations])) {
            call PrepareAnnotations {
                input:
                annotation = annotation,
                codon_table = codon_table,
                runtime_attr_override = task_attr
            }

            GenomeProteins cleaned_up = object { genome: annotation.genome, annotation_gff: PrepareAnnotations.cleaned_up_gff, protein_sequences: PrepareAnnotations.proteins}
            String out_prefix = sub(basename(cleaned_up.annotation_gff), "\.[^/.]+$", "")
            String ref_prefix = sub(basename(cleaned_up.genome), "\.[^/.]+$", "")

            call AlignProteinsSpaln {
                input:
                    genome_index = select_first([IndexGenomeSpaln.genome_index]),
                    genome_to_annotate = genome_to_annotate,
                    genome_proteins = cleaned_up,
                    ref_prefix = ref_prefix,
                    out_prefix = out_prefix,
                    species = species,
                    spaln_extra_parameters = spaln_extra_parameters,
                    runtime_attr_override = aln_attr
            }

            call AlignProteinsMiniprot {
                input:
                    genome_index = IndexGenomeMiniprot.genome_index,
                    genome_to_annotate = genome_to_annotate,
                    genome_proteins = cleaned_up,
                    ref_prefix = ref_prefix,
                    out_prefix = out_prefix,
                    runtime_attr_override = aln_attr
            }

            call PrepareAlignmentsSpaln {
                input:
                annotation_cdnas = PrepareAnnotations.cdnas,
                annotation_bed = PrepareAnnotations.bed,
                raw_alignments = AlignProteinsSpaln.raw_alignments,
                aln_prefix = sub(basename(AlignProteinsSpaln.raw_alignments), "\\.alignment\\.gff$", ""),
                genome_to_annotate = genome_to_annotate,
                runtime_attr_override = task_attr
            }

            call PrepareAlignmentsMiniprot {
                input:
                annotation_cdnas = PrepareAnnotations.cdnas,
                annotation_bed = PrepareAnnotations.bed,
                raw_alignments = AlignProteinsMiniprot.raw_alignments,
                aln_prefix = sub(basename(AlignProteinsMiniprot.raw_alignments), "\\.alignment\\.gff$", ""),
                genome_to_annotate = genome_to_annotate,
                runtime_attr_override = task_attr
            }


            String spaln_aln_prefix = sub(basename(PrepareAlignmentsSpaln.clean_alignments), "\.[^/.]+$", "")
            String miniprot_aln_prefix = sub(basename(PrepareAlignmentsMiniprot.clean_alignments), "\.[^/.]+$", "")

            call ScoreAlignmentsSpaln {
                input:
                cdnas = PrepareAlignmentsSpaln.cdnas,
                bed = PrepareAlignmentsSpaln.bed,
                groups = PrepareAlignmentsSpaln.groups,
                aln_prefix = spaln_aln_prefix,
                ref_prefix = ref_prefix,
                runtime_attr_override = score_attr
            }

            call ScoreAlignmentsMiniprot {
                input:
                cdnas = PrepareAlignmentsMiniprot.cdnas,
                bed = PrepareAlignmentsMiniprot.bed,
                groups = PrepareAlignmentsMiniprot.groups,
                aln_prefix = miniprot_aln_prefix,
                ref_prefix = ref_prefix,
                runtime_attr_override = score_attr
            }

            call CombineResultsSpaln {
                input:
                alignment_compare_detail = ScoreAlignmentsSpaln.alignment_compare_detail,
                alignment_gff = PrepareAlignmentsSpaln.clean_alignments,
                runtime_attr_override = task_attr
            }

            call CombineResultsMiniprot {
                input:
                alignment_compare_detail = ScoreAlignmentsMiniprot.alignment_compare_detail,
                alignment_gff = PrepareAlignmentsMiniprot.clean_alignments,
                runtime_attr_override = task_attr
            }
        }
        call ScoreSummary {
            input:
            alignments_compare = flatten(select_all([ScoreAlignmentsSpaln.alignment_compare_detail, ScoreAlignmentsMiniprot.alignment_compare_detail])),
            runtime_attr_override = task_attr
        }
    }

    Array[File] all_spaln_alignments = flatten(select_all([CombineResultsSpaln.augmented_alignments_gff, PrepareProteinAlignmentsSpaln.clean_alignments]))
    scatter (alignment in all_spaln_alignments) {
        call CombineXspecies as CombineXspeciesSpaln {
            input:
            alignment = alignment,
            alignments = all_spaln_alignments,
            runtime_attr_override = task_attr
        }
    }
    Array[File] all_miniprot_alignments = flatten(select_all([CombineResultsMiniprot.augmented_alignments_gff, PrepareProteinAlignmentsMiniprot.clean_alignments]))
    scatter (alignment in all_miniprot_alignments) {
        call CombineXspecies as CombineXspeciesMiniprot {
            input:
            alignment = alignment,
            alignments = all_miniprot_alignments,
            runtime_attr_override = task_attr
        }
    }

    # With the results so far and some optional inputs such as assemblies and junctions
    # consolidate them in a mikado run that should incorporate UTRs and add junctions to the scoring

    Array[File] all_combined_xspecies = flatten(select_all([CombineXspeciesSpaln.xspecies_scored_alignment, CombineXspeciesMiniprot.xspecies_scored_alignment]))

    call Mikado {
        input:
        reference = genome_to_annotate,
        extra_config = mikado_config,
        xspecies = all_combined_xspecies,
        output_prefix = output_prefix,
        runtime_attr_override = mikado_attr
    }

    call MikadoPick {
        input:
        config_file = Mikado.mikado_config,
        scoring_file = mikado_scoring,
        mikado_db = Mikado.mikado_db,
        transcripts = Mikado.prepared_gtf,
        output_prefix = output_prefix,
        runtime_attr_override = mikado_attr
    }

    call MikadoSummaryStats {
        input:
        stats = [MikadoPick.stats],
        output_prefix = output_prefix,
        runtime_attr_override = task_attr
    }

    output {

        # cleaned up annotations
        Array[File]? clean_annotations = PrepareAnnotations.cleaned_up_gff
        Array[File]? clean_annotations_proteins = PrepareAnnotations.proteins
        # Intermediates and Alignments
        Array[File]? spaln_intermediates = flatten(select_all([AlignProteinSequencesSpaln.raw_alignments, AlignProteinsSpaln.raw_alignments]))
        Array[File]? miniprot_intermediates = flatten(select_all([AlignProteinSequencesMiniprot.raw_alignments, AlignProteinsMiniprot.raw_alignments]))

        # add both spaln and miniprot results to mgc_evaluation and mgc_evaluation_detail
        Array[File]? mgc_evaluation_spaln = ScoreAlignmentsSpaln.alignment_compare
        Array[File]? mgc_evaluation_spaln_detail = ScoreAlignmentsSpaln.alignment_compare_detail
        Array[File]? mgc_evaluation_miniprot = ScoreAlignmentsMiniprot.alignment_compare
        Array[File]? mgc_evaluation_miniprot_detail = ScoreAlignmentsMiniprot.alignment_compare_detail

        # add both spaln and miniprot results to xspecies_scored_alignment
        Array[File]? xspecies_combined_alignments_spaln = CombineXspeciesSpaln.xspecies_scored_alignment
        Array[File]? xspecies_combined_alignments_miniprot = CombineXspeciesMiniprot.xspecies_scored_alignment

        # add both spaln and miniprot results to combined_results
        Array[File]? annotation_filter_stats = PrepareAnnotations.stats
        Array[File]? spaln_alignment_filter_stats = flatten(select_all([PrepareAlignmentsSpaln.stats, PrepareProteinAlignmentsSpaln.stats]))
        Array[File]? miniprot_alignment_filter_stats = flatten(select_all([PrepareAlignmentsMiniprot.stats, PrepareProteinAlignmentsMiniprot.stats]))
        File? mgc_score_summary = ScoreSummary.summary_table

        # add mikado_results
        File loci = MikadoPick.loci
        File scores = MikadoPick.scores
        File metrics = MikadoPick.metrics
        File stats = MikadoPick.stats
        File stats_summary = MikadoSummaryStats.summary
    }
}

task Mikado {
    input {
        File reference
        File extra_config
        Int min_cdna_length = 100
        Int max_intron_length = 1000000
        Array[File] xspecies
        File? utrs
        File? junctions
        String output_prefix
        RuntimeAttr? runtime_attr_override
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

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    command <<<
        set -euxo pipefail
        export TMPDIR=/tmp
        # Create the lists file
        count=1
        for i in ~{sep=" " xspecies}
        do
        bname=$(basename ${i})
        label=${bname%%.*}
        echo -e "${i}\t${count}_${label}\tTrue\t1\tTrue" >> list.txt
        count=$((count + 1))
        done

        apply_pad="--no-pad "
        if [ "" != "~{utrs}" ]
        then
            apply_pad=""
            label="UTRs"
            gffread -vE ~{utrs} > preprocessed_utrs.gff
            echo -e "$(pwd -P)/preprocessed_utrs.gff\t${label}\tTrue\t0\tFalse" >> list.txt
        fi

        # mikado configure
        mikado configure $apply_pad --only-reference-update --check-references --reference-update \
        --max-intron-length ~{max_intron_length} --minimum-cdna-length ~{min_cdna_length} \
        --reference ~{reference} --list=list.txt --copy-scoring plant.yaml original-mikado.yaml

        yaml-merge -s original-mikado.yaml ~{"-m " + extra_config} -o ~{output_prefix}-mikado.yaml

        # mikado prepare
        mikado prepare --procs=~{task_cpus} --json-conf ~{output_prefix}-mikado.yaml -od ~{output_prefix}-mikado

        # mikado serialise
        mikado serialise --force ~{"--junctions " + junctions} --transcripts ~{output_prefix}-mikado/mikado_prepared.fasta \
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
    input{
        File config_file
        File? extra_config
        File scoring_file
        File transcripts
        File mikado_db
        String output_prefix
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 16,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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
    mikado pick --source Mikado_~{output_prefix} --procs=~{task_cpus} --scoring-file ~{scoring_file} \
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

task ScoreSummary {
    input {
        Array[File] alignments_compare
        RuntimeAttr? runtime_attr_override
    }

    output {
        File summary_table = "ScoreAlignments/Junction_precision_binning.tsv"
    }

    Int cpus = 1
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 20,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    command <<<
        set -euxo pipefail
        mkdir ScoreAlignments
        cd ScoreAlignments
        create_alignment_summary.py ~{sep=' ' alignments_compare} > Junction_precision_binning.tsv
    >>>

    runtime {
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }
}

task CombineXspecies {
    input {
        File alignment
        Array[File] alignments
        RuntimeAttr? runtime_attr_override
    }

    String out_prefix = sub(basename(alignment), "\.[^/.]+$", "")

    output {
        File xspecies_scored_alignment = out_prefix + ".xspecies_scores.gff"
    }

    Int cpus = 1
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 16,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    command <<<
        set -euxo pipefail
        for o in ~{sep=' ' alignments};
        do
            if [ "~{alignment}" != "${o}" ]
            then
                o_name=$(basename $o)
                mikado compare --processes 1 -r ~{alignment} -p $o -o ~{out_prefix}_vs_${o_name%.*}
            fi
        done

        combine_alignment_with_compare --alignment ~{alignment} --compare ~{out_prefix}_vs_*.refmap > ~{out_prefix}.xspecies_scores.gff
    >>>

    runtime {
        continueOnReturnCode: [0]
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }
}

task IndexGenomeSpaln {
    input {
        File genome
        RuntimeAttr? runtime_attr_override
    }

    output {
        Array[File] genome_index = glob("genome_to_annotate.*")
    }

    Int cpus = 6
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 16,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    command <<<
        set -euxo pipefail
        ln -s ~{genome} genome_to_annotate.fasta
        spaln -W -KP -t~{task_cpus} genome_to_annotate.fasta
    >>>

    runtime {
        continueOnReturnCode: [0, -1]
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }
}

task IndexGenomeMiniprot {
    input {
        File genome
        RuntimeAttr? runtime_attr_override
    }

    output {
        File genome_index = "genome_to_annotate.mpi"
    }

    Int cpus = 6
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 50,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    command <<<
        set -euxo pipefail
        ln -s ~{genome} genome_to_annotate.fasta
        miniprot -t~{task_cpus} -d genome_to_annotate.mpi genome_to_annotate.fasta
    >>>

    runtime {
        continueOnReturnCode: [0, -1]
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }
}

task PrepareAnnotations {
    input {
        GenomeAnnotation annotation
        String out_prefix = sub(basename(annotation.annotation_gff),  "\.[^/.]+$", "")
        Int min_cds_len = 20 # nts
        Int max_intron_len = 200000
        Int min_exon_len = 10
        Int codon_table
        Array[String] filters = "all"
        RuntimeAttr? runtime_attr_override
    }

    output {
        File cleaned_up_gff = 'Annotations/' + out_prefix + ".clean.extra_attr.gff"
        File stats = 'Annotations/' + out_prefix + ".annotation.stats"
        File proteins = 'Annotations/' + out_prefix + ".proteins.fa"
        File cdnas = 'Annotations/' + out_prefix + ".cdna.fa"
        File bed = 'Annotations/' + out_prefix + ".bed"
    }

    Int cpus = 4
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 30,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    command <<<
        set -euxo pipefail
        mkdir Annotations/
        cd Annotations
        case ~{annotation.annotation_gff} in
            *.gz) gunzip <~{annotation.annotation_gff};;
            *)    cat ~{annotation.annotation_gff};;
        esac | gffread --keep-genes -E -o sorted_annotation.gff > ~{annotation.annotation_gff}.gffread.out 2> ~{annotation.annotation_gff}.gffread.err
        xspecies_cleanup --merge --filters ~{sep=" " filters} \
        --annotation sorted_annotation.gff -g ~{annotation.genome} \
        --max_intron ~{max_intron_len} --min_exon ~{min_exon_len} --min_protein ~{min_cds_len} \
        -x ~{out_prefix}.cdna.fa --bed ~{out_prefix}.bed \
        -o ~{out_prefix}.clean.extra_attr.gff > ~{out_prefix}.annotation.stats
        cat ~{out_prefix}.cdna.fa | seqkit translate -T ~{codon_table} -o ~{out_prefix}.proteins.fa
    >>>

    runtime {
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }
}

task AlignProteinsSpaln {
    input {
        Array[File] genome_index
        File genome_to_annotate
        GenomeProteins genome_proteins
        String species
        String spaln_extra_parameters
        Int min_spaln_intron_length = 20
        Int min_coverage = 80
        Int min_identity = 50
        Int max_per_query = 4
        String clip = "no_clip"
        Int term5c_len = 10
        Int term5i_len = 200000
        Int term3c_len = 10
        Int term3i_len = 200000
        String ref_prefix
        String out_prefix
        RuntimeAttr? runtime_attr_override
    }
    Int cpus = 8
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 32,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    output {
        File unfiltered_alignments = "spaln_atn_"+ref_prefix+".raw.alignment.gff"
        File raw_alignments = "spaln_atn_"+ref_prefix+".alignment.gff"
    }

    command <<<
        set -euxo pipefail
        ln ~{sub(genome_index[0], "\.[^/.]+$", "")}.* .
        ln -s ~{genome_proteins.protein_sequences} "./~{basename(genome_proteins.protein_sequences)}"
        spaln ~{spaln_extra_parameters} -t~{task_cpus} -KP -O12 -M~{max_per_query}.~{max_per_query} ~{"-T"+species} -dgenome_to_annotate -o ~{out_prefix} -yL~{min_spaln_intron_length} "~{basename(genome_proteins.protein_sequences)+" P"}"
        sortgrcd -O4 ~{out_prefix}.grd | tee ~{out_prefix}.s | spaln2gff -s "~{ref_prefix}.spaln" > spaln_atn_~{ref_prefix}.raw.alignment.gff
        sortgrcd -O4 ~{out_prefix}.grd | tee ~{out_prefix}.s | spaln2gff --min_coverage ~{min_coverage} --min_identity ~{min_identity} -s "~{ref_prefix}.spaln" | filter_GFF3_terminals.py --protein_alias protein_alias --term5c_len ~{term5c_len}  --term5i_len ~{term5i_len} --term3c_len ~{term3c_len} --term3i_len ~{term3i_len} --clip ~{clip} ~{basename(genome_proteins.protein_sequences)} - > spaln_atn_~{ref_prefix}.alignment.gff
    >>>

    runtime {
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }
}


task AlignProteinsMiniprot {
    input {
        File genome_index
        File genome_to_annotate
        GenomeProteins genome_proteins
        Int max_intron_length = 20000
        Int min_coverage = 80
        Int min_identity = 50
        Int max_per_query = 4
        String ref_prefix
        String out_prefix
        RuntimeAttr? runtime_attr_override
    }
    Int cpus = 8
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 32,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    output {
        File unfiltered_alignments = "miniprot_atn_"+ref_prefix+".raw.alignment.gff"
        File raw_alignments = "miniprot_atn_"+ref_prefix+".alignment.gff"
    }

    command <<<
        set -euxo pipefail
        ln -s ~{genome_index}
        ln -s ~{genome_proteins.protein_sequences} "./~{basename(genome_proteins.protein_sequences)}"
        miniprot -t~{task_cpus} -G~{max_intron_length} --outn=~{max_per_query} --gff ~{genome_index} ~{genome_proteins.protein_sequences} > miniprot_atn_~{ref_prefix}.raw.alignment.gff
        miniprot2GFF3.py --source "~{ref_prefix}.miniprot" --identity ~{min_identity} --coverage ~{min_coverage} ~{basename(genome_proteins.protein_sequences)} miniprot_atn_~{ref_prefix}.raw.alignment.gff > miniprot_atn_~{ref_prefix}.alignment.gff
    >>>

    runtime {
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }

}


task AlignProteinSequencesSpaln {
    input {
        Array[File] genome_index
        File proteins
        String species
        String spaln_extra_parameters
        Int min_spaln_intron_length = 20
        Int min_coverage = 80
        Int min_identity = 50
        Int max_per_query = 4
        String clip = "no_clip"
        Int term5c_len = 10
        Int term5i_len = 200000
        Int term3c_len = 10
        Int term3i_len = 200000
        String out_prefix
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 32,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    output {
        File unfiltered_alignments = "spaln_ptn_"+basename(proteins)+".raw.alignment.gff"
        File raw_alignments = "spaln_ptn_"+basename(proteins)+".alignment.gff"
    }

    command <<<
        set -euxo pipefail
        ln ~{sub(genome_index[0], "\.[^/.]+$", "")}.* .
        ln -s ~{proteins} .
        spaln ~{spaln_extra_parameters} -t~{task_cpus} -KP -O12 -M~{max_per_query}.~{max_per_query} ~{"-T"+species} -dgenome_to_annotate -o ~{out_prefix} -yL~{min_spaln_intron_length} ~{basename(proteins)}
        sortgrcd -O4 ~{out_prefix}.grd | tee ~{out_prefix}.s | spaln2gff -s "~{out_prefix}.spaln" > spaln_ptn_~{basename(proteins)}.raw.alignment.gff
        sortgrcd -O4 ~{out_prefix}.grd | tee ~{out_prefix}.s | spaln2gff --min_coverage ~{min_coverage} --min_identity ~{min_identity} -s "~{out_prefix}.spaln" | filter_GFF3_terminals.py --protein_alias protein_alias --term5c_len ~{term5c_len}  --term5i_len ~{term5i_len} --term3c_len ~{term3c_len} --term3i_len ~{term3i_len} --clip ~{clip} ~{basename(proteins)} - > spaln_ptn_~{basename(proteins)}.alignment.gff
    >>>


    runtime {
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }

}

task AlignProteinSequencesMiniprot {
    input {
        File genome_index
        File proteins
        Int max_intron_length = 20000
        Int min_coverage = 80
        Int min_identity = 50
        Int max_per_query = 4
        String out_prefix
        RuntimeAttr? runtime_attr_override
    }
    Int cpus = 8
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 32,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    output {
        File unfiltered_alignments = "miniprot_ptn_"+basename(proteins)+".raw.alignment.gff"
        File raw_alignments = "miniprot_ptn_"+basename(proteins)+".alignment.gff"
    }

    command <<<
        set -euxo pipefail
        ln -s ~{genome_index}
        ln -s ~{proteins} .
        miniprot -t~{task_cpus} -G~{max_intron_length} --outn=~{max_per_query} --gff ~{genome_index} ~{basename(proteins)} > miniprot_ptn_~{basename(proteins)}.raw.alignment.gff
        miniprot2GFF3.py --source "~{basename(proteins)}.miniprot" --identity ~{min_identity} --coverage ~{min_coverage} ~{basename(proteins)} miniprot_ptn_~{basename(proteins)}.raw.alignment.gff > miniprot_ptn_~{basename(proteins)}.alignment.gff
    >>>

    runtime {
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }

}

task PrepareProteinAlignments {
    input {
        String aln_prefix

        Int max_intron_len = 200000
        Int min_filter_exon_len = 10
        Int min_cds_len = 20
        Boolean show_intron_len = false
        Array[String] filters = "none"

        File raw_alignments
        File genome_to_annotate

        RuntimeAttr? runtime_attr_override
    }

    output {
        File clean_alignments = aln_prefix + ".alignment.stop_extended.extra_attr.gff"
        File stats = aln_prefix + ".alignment.stats"
    }

    Int cpus = 4
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 20,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])


    command <<<
        set -euxo pipefail

        xspecies_cleanup ~{if show_intron_len then "--show_intron_len" else ""} \
        --filters ~{sep=" " filters} --max_intron ~{max_intron_len} --min_exon ~{min_filter_exon_len} --min_protein ~{min_cds_len} \
        -g ~{genome_to_annotate} -A ~{raw_alignments} --bed ~{aln_prefix}.alignment.bed \
        -x ~{aln_prefix}.alignment.cdna.fa \
        -o ~{aln_prefix}.alignment.stop_extended.extra_attr.gff > ~{aln_prefix}.alignment.stats

    >>>

    runtime {
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }

}

task PrepareAlignmentsSpaln {
    input {
        String aln_prefix

        Int max_intron_len = 200000
        Int min_filter_exon_len = 10
        Int min_cds_len = 20
        Boolean show_intron_len = false
        Array[String] filters = "none"

        File raw_alignments
        File genome_to_annotate
        File annotation_cdnas
        File annotation_bed

        RuntimeAttr? runtime_attr_override
    }

    output {
        File clean_alignments = aln_prefix + ".alignment.stop_extended.extra_attr.gff"
        File stats = aln_prefix + ".alignment.stats"
        File cdnas = "all_cdnas.fa"
        File bed = "all_cdnas.bed"
        File groups = "groups.txt"
    }

    Int cpus = 4
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 20,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])


    command <<<
        set -euxo pipefail

        xspecies_cleanup ~{if show_intron_len then "--show_intron_len" else ""} \
        --filters ~{sep=" " filters} --max_intron ~{max_intron_len} --min_exon ~{min_filter_exon_len} --min_protein ~{min_cds_len} \
        -g ~{genome_to_annotate} -A ~{raw_alignments} --bed ~{aln_prefix}.alignment.bed \
        -x ~{aln_prefix}.alignment.cdna.fa \
        -o ~{aln_prefix}.alignment.stop_extended.extra_attr.gff > ~{aln_prefix}.alignment.stats

        create_mgc_groups_from_gff.py --alignments_gff ~{aln_prefix}.alignment.stop_extended.extra_attr.gff

        cat ~{aln_prefix}.alignment.cdna.fa ~{annotation_cdnas}  > all_cdnas.fa
        cat ~{aln_prefix}.alignment.bed ~{annotation_bed} > all_cdnas.bed
    >>>

    runtime {
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }

}

task PrepareAlignmentsMiniprot {
    input {
        String aln_prefix

        Int max_intron_len = 200000
        Int min_filter_exon_len = 10
        Int min_cds_len = 20
        Boolean show_intron_len = false
        Array[String] filters = "none"

        File raw_alignments
        File genome_to_annotate
        File annotation_cdnas
        File annotation_bed

        RuntimeAttr? runtime_attr_override
    }

    output {
        File clean_alignments = aln_prefix + ".alignment.stop_extended.extra_attr.gff"
        File stats = aln_prefix + ".alignment.stats"
        File cdnas = "all_cdnas.fa"
        File bed = "all_cdnas.bed"
        File groups = "groups.txt"
    }

    Int cpus = 4
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 20,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])


    command <<<
        set -euxo pipefail

        xspecies_cleanup ~{if show_intron_len then "--show_intron_len" else ""} \
        --filters ~{sep=" " filters} --max_intron ~{max_intron_len} --min_exon ~{min_filter_exon_len} --min_protein ~{min_cds_len} \
        -g ~{genome_to_annotate} -A ~{raw_alignments} --bed ~{aln_prefix}.alignment.bed \
        -x ~{aln_prefix}.alignment.cdna.fa \
        -o ~{aln_prefix}.alignment.stop_extended.extra_attr.gff > ~{aln_prefix}.alignment.stats

        create_mgc_groups_from_gff.py --alignments_gff ~{aln_prefix}.alignment.stop_extended.extra_attr.gff

        cat ~{aln_prefix}.alignment.cdna.fa ~{annotation_cdnas}  > all_cdnas.fa
        cat ~{aln_prefix}.alignment.bed ~{annotation_bed} > all_cdnas.bed
    >>>

    runtime {
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }

}

task ScoreAlignmentsSpaln {
    input {
        File cdnas
        File bed
        File groups
        String aln_prefix
        String ref_prefix
        RuntimeAttr? runtime_attr_override
    }


    output {
        File alignment_compare = 'ScoreAlignments/' + "comp_"+aln_prefix+"_"+ref_prefix+".tab"
        File alignment_compare_detail = 'ScoreAlignments/' + "comp_"+aln_prefix+"_"+ref_prefix+"_detail.tab"
    }

    Int cpus = 4
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 32,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    command <<<
        set -euxo pipefail
        mkdir ScoreAlignments/
        cd ScoreAlignments
        multi_genome_compare -t ~{task_cpus} --groups ~{groups} \
        --cdnas ~{cdnas} --bed12 ~{bed} \
        -o comp_~{aln_prefix}_~{ref_prefix}.tab -d comp_~{aln_prefix}_~{ref_prefix}_detail.tab
    >>>

    runtime {
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }

}

task ScoreAlignmentsMiniprot {
    input {
        File cdnas
        File bed
        File groups
        String aln_prefix
        String ref_prefix
        RuntimeAttr? runtime_attr_override
    }


    output {
        File alignment_compare = 'ScoreAlignments/' + "comp_"+aln_prefix+"_"+ref_prefix+".tab"
        File alignment_compare_detail = 'ScoreAlignments/' + "comp_"+aln_prefix+"_"+ref_prefix+"_detail.tab"
    }

    Int cpus = 4
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 32,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    command <<<
        set -euxo pipefail
        mkdir ScoreAlignments/
        cd ScoreAlignments
        multi_genome_compare -t ~{task_cpus} --groups ~{groups} \
        --cdnas ~{cdnas} --bed12 ~{bed} \
        -o comp_~{aln_prefix}_~{ref_prefix}.tab -d comp_~{aln_prefix}_~{ref_prefix}_detail.tab
    >>>

    runtime {
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }

}

task CombineResultsSpaln {
    input {
        File alignment_compare_detail
        File alignment_gff
        Int? exon_f1_filter
        Int? junction_f1_filter
        RuntimeAttr? runtime_attr_override
    }

    String aln_prefix = sub(basename(alignment_gff), "\.[^/.]+$", "")

    output {
        File augmented_alignments_gff = aln_prefix+".mgc.gff"
    }

    Int cpus = 4
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 32,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    command <<<
        set -euxo pipefail
        combine_alignments_with_mgc --detail ~{alignment_compare_detail} ~{'--exon_f1_filter ' + exon_f1_filter} \
        ~{'--junction_f1_filter ' + junction_f1_filter} --gff ~{alignment_gff} -o ~{aln_prefix}.mgc.gff
    >>>

    runtime {
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }
}

task CombineResultsMiniprot {
    input {
        File alignment_compare_detail
        File alignment_gff
        Int? exon_f1_filter
        Int? junction_f1_filter
        RuntimeAttr? runtime_attr_override
    }

    String aln_prefix = sub(basename(alignment_gff), "\.[^/.]+$", "")

    output {
        File augmented_alignments_gff = aln_prefix+".mgc.gff"
    }

    Int cpus = 4
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 32,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    command <<<
        set -euxo pipefail
        combine_alignments_with_mgc --detail ~{alignment_compare_detail} ~{'--exon_f1_filter ' + exon_f1_filter} \
        ~{'--junction_f1_filter ' + junction_f1_filter} --gff ~{alignment_gff} -o ~{aln_prefix}.mgc.gff
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
        RuntimeAttr? runtime_attr_override
    }

    output {
        File summary = output_prefix + ".summary.stats.tsv"
    }

    Int cpus = 1
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 20,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    command <<<
    mikado_summary_stats ~{sep=" " stats} > ~{output_prefix}.summary.stats.tsv
    >>>

    runtime {
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }
}
