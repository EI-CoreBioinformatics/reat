version 1.0

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
        Array[GenomeAnnotation] annotations
        File genome_to_annotate
    }

    call IndexGenome {
        input:
        genome = genome_to_annotate
    }

    scatter (annotation in annotations) {
        call PrepareAnnotations {
            input:
            annotation = annotation
        }

        GenomeProteins cleaned_up = object { genome: annotation.genome, annotation_gff: PrepareAnnotations.cleaned_up_gff, protein_sequences: PrepareAnnotations.proteins}

        call AlignProteins {
            input:
                genome_index = IndexGenome.genome_index,
                genome_to_annotate = genome_to_annotate,
                genome_proteins = cleaned_up
        }

        call ScoreAlignments {
            input:
            alignments = AlignProteins.alignments,
            genome_to_annotate = genome_to_annotate,
            genome_proteins = cleaned_up
        }
    }
    output {
        Array[File] alignments = AlignProteins.alignments
        Array[File] mgc_evaluation = ScoreAlignments.alignment_compare
        Array[File] mgc_evaluation_detail = ScoreAlignments.alignment_compare_detail
    }
}

task IndexGenome {
    input {
        File genome
    }

    output {
        Array[File] genome_index = glob("genome_to_annotate.*")
    }

    runtime {
        continueOnReturnCode: [0, -1]
    }

    command <<<
        set -euxo pipefail
        ln ~{genome} genome_to_annotate.fasta
        spaln -W -KP -t12 genome_to_annotate.fasta
    >>>
}

task PrepareAnnotations {
    input {
        GenomeAnnotation annotation
        String out_prefix = sub(basename(annotation.annotation_gff),  "\.[^/.]+$", "")
        Int min_cds_len = 20 # nts
    }

    output {
        File cleaned_up_gff = out_prefix + ".gff"
        File proteins = out_prefix + ".proteins.fa"
    }

    command <<<
        set -euxo pipefail
        xspecies_cleanup --filters all --annotation ~{annotation.annotation_gff} --genome ~{annotation.genome} --min_protein ~{min_cds_len} -y ~{out_prefix}.proteins.fa -o ~{out_prefix}.gff
    >>>
}

task AlignProteins {
    input {
        Array[File] genome_index
        File genome_to_annotate
        GenomeProteins genome_proteins
        String species = "Eudicoty" # TODO: Need a lookup table from the original file
        Int min_exon_len = 20
        Int min_coverage = 80
        Int min_identity = 50
    }

    Int num_cpus = 12
#    String out_prefix = sub(basename(genome_index[0]), "\.[^/.]+$", "")
    String out_prefix = sub(basename(genome_proteins.annotation_gff), "\.[^/.]+$", "")

    output {
        File alignments = out_prefix+".alignment.gff"
    }

    command <<<
        set -euxo pipefail
        ln ~{sub(genome_index[0], "\.[^/.]+$", "")}.* .
        spaln -t~{num_cpus} -KP -O0,12 -Q7 ~{"-T"+species} -dgenome_to_annotate -o ~{out_prefix} -yL~{min_exon_len} ~{genome_proteins.protein_sequences}
        sortgrcd -O4 ~{out_prefix}.grd | tee ~{out_prefix}.s | spaln2gff --min_coverage ~{min_coverage} --min_identity ~{min_identity} -s "spaln" > ~{out_prefix}.alignment.gff
    >>>


}

task ScoreAlignments {
    input {
        File alignments
        File genome_to_annotate
        GenomeProteins genome_proteins
        Int max_intron_len = 2000
    }

    String aln_prefix = sub(basename(alignments), "\.[^/.]+$", "")
    String ref_prefix = sub(basename(genome_proteins.genome), "\.[^/.]+$", "")

    output {
        File scored_alignments = "scored_alignments.gff"
        File alignment_compare = "comp_"+aln_prefix+"_"+ref_prefix+".tab"
        File alignment_compare_detail = "comp_"+aln_prefix+"_"+ref_prefix+"_detail.tab"
    }

    # TODO:
    # - Annotate models with the count of introns exceeding max_intron_len per model
    # - Annotate exons with splice signal (to account for non-canonical splicing)
    # - Compare models to original structure (mgc approach preferably on-line)

    # For now extract using `gffread -g call-AlignProteins/shard-0/inputs/1074434829/GCF_000184785.3_Aflo_1.1_genomic.fna
    # -x GCF_000184785.3_Aflo_1.1_transcripts.fa --bed -o GCF_000184785.3_Aflo_1.1.bed`
    # And similar for the alignments
    # Then call mgc with both beds and the groups created by `create_mgc_groups`
    # Finally, collate all the metrics and scores into a final value or set of values
    command <<<
        set -euxo pipefail
        xspecies_cleanup --filters none -g ~{genome_to_annotate} -A ~{alignments} -o ~{aln_prefix}.stop_extended.extra_attr.gff
        xspecies_cleanup --filters all -g ~{genome_proteins.genome} -a ~{genome_proteins.annotation_gff} -o ~{ref_prefix}.clean.extra_attr.gff

        gffread -g ~{genome_to_annotate} -T -w ~{aln_prefix}.cdna.fa -o ~{aln_prefix}.gtf ~{alignments}
        gffread -g ~{genome_proteins.genome} -T -w ~{ref_prefix}.cdna.fa -o ~{ref_prefix}.gtf ~{genome_proteins.annotation_gff}

        mikado util convert -if gtf -of bed12 ~{aln_prefix}.gtf ~{aln_prefix}.bed12
        mikado util convert -if gtf -of bed12 ~{ref_prefix}.gtf ~{ref_prefix}.bed12

        create_mgc_groups -f ~{aln_prefix}.cdna.fa

        cat ~{aln_prefix}.cdna.fa ~{ref_prefix}.cdna.fa > all_cdnas.fa
        cat ~{aln_prefix}.bed12 ~{ref_prefix}.bed12 > all_cdnas.bed

        multi_genome_compare.py -t 12 --groups groups.txt --cdnas all_cdnas.fa --bed12 all_cdnas.bed -o comp_~{aln_prefix}_~{ref_prefix}.tab -d comp_~{aln_prefix}_~{ref_prefix}_detail.tab
    >>>
}