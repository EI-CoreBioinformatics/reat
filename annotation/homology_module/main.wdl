version 1.0
import "./rt_struct.wdl"

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
        RuntimeAttr? index_attr
        RuntimeAttr? score_attr
        RuntimeAttr? aln_attr
    }

    call IndexGenome {
        input:
        genome = genome_to_annotate,
        runtime_attr_override = index_attr
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
                genome_proteins = cleaned_up,
                runtime_attr_override = aln_attr
        }

        call PrepareAlignments {
            input:
            alignments = AlignProteins.alignments,
            genome_to_annotate = genome_to_annotate,
            genome_proteins = cleaned_up
        }
        String aln_prefix = sub(basename(AlignProteins.alignments), "\.[^/.]+$", "")
        String ref_prefix = sub(basename(cleaned_up.genome), "\.[^/.]+$", "")

        call ScoreAlignments {
            input:
            cdnas = PrepareAlignments.cdnas,
            bed = PrepareAlignments.bed,
            groups = PrepareAlignments.groups,
            aln_prefix = aln_prefix,
            ref_prefix = ref_prefix,
            runtime_attr_override = score_attr
        }
    }
    output {
        Array[File] clean_annotations = PrepareAnnotations.cleaned_up_gff
        Array[File] alignments = AlignProteins.alignments
        Array[File] mgc_evaluation = ScoreAlignments.alignment_compare
        Array[File] mgc_evaluation_detail = ScoreAlignments.alignment_compare_detail
    }
}

task IndexGenome {
    input {
        File genome
        RuntimeAttr? runtime_attr_override
    }

    output {
        Array[File] genome_index = glob("genome_to_annotate.*")
    }

    Int cpus = 12
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 16,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = runtime_attr.cpu_cores

    command <<<
        set -euxo pipefail
        ln ~{genome} genome_to_annotate.fasta
        spaln -W -KP -t~{task_cpus} genome_to_annotate.fasta
    >>>

    runtime {
        continueOnReturnCode: [0, -1]
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

}

task PrepareAnnotations {
    input {
        GenomeAnnotation annotation
        String out_prefix = sub(basename(annotation.annotation_gff),  "\.[^/.]+$", "")
        Int min_cds_len = 20 # nts
        String filters = "all"
    }

    output {
        File cleaned_up_gff = out_prefix + ".clean.extra_attr.gff"
        File proteins = out_prefix + ".proteins.fa"
    }

    command <<<
        set -euxo pipefail
        xspecies_cleanup --merge --filters none --annotation ~{annotation.annotation_gff} --genome ~{annotation.genome} --min_protein ~{min_cds_len} -y ~{out_prefix}.proteins.fa -o ~{out_prefix}.clean.extra_attr.gff
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
        String filters = "none"
        RuntimeAttr? runtime_attr_override
    }
    Int cpus = 6
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 32,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = runtime_attr.cpu_cores
#    String out_prefix = sub(basename(genome_index[0]), "\.[^/.]+$", "")
    String out_prefix = sub(basename(genome_proteins.annotation_gff), "\.[^/.]+$", "")
    String ref_prefix = sub(basename(genome_proteins.genome), "\.[^/.]+$", "")

    output {
        File alignments = ref_prefix+".alignment.stop_extended.extra_attr.gff"
    }

    command <<<
        set -euxo pipefail
        ln ~{sub(genome_index[0], "\.[^/.]+$", "")}.* .
        spaln -t~{task_cpus} -KP -O0,12 -Q7 ~{"-T"+species} -dgenome_to_annotate -o ~{out_prefix} -yL~{min_exon_len} ~{genome_proteins.protein_sequences}
        sortgrcd -O4 ~{out_prefix}.grd | tee ~{out_prefix}.s | spaln2gff --min_coverage ~{min_coverage} --min_identity ~{min_identity} -s "spaln" > ~{ref_prefix}.alignment.gff

        xspecies_cleanup --filters none -g ~{genome_to_annotate} -A ~{ref_prefix}.alignment.gff -o ~{ref_prefix}.alignment.stop_extended.extra_attr.gff
    >>>

    runtime {
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

}

task PrepareAlignments {
    input {
        File alignments
        File genome_to_annotate
        GenomeProteins genome_proteins
        Int max_intron_len = 200000
    }

    String aln_prefix = sub(basename(alignments), "\.[^/.]+$", "")
    String ref_prefix = sub(basename(genome_proteins.genome), "\.[^/.]+$", "")

    output {
        File cdnas = "all_cdnas.fa"
        File bed = "all_cdnas.bed"
        File groups = "groups.txt"
    }

    command <<<
        set -euxo pipefail
        gffread -g ~{genome_to_annotate} -T -w ~{aln_prefix}.cdna.fa -o ~{aln_prefix}.gtf ~{alignments}
        gffread -g ~{genome_proteins.genome} -T -w ~{ref_prefix}.cdna.fa -o ~{ref_prefix}.gtf ~{genome_proteins.annotation_gff}

        mikado util convert -if gtf -of bed12 ~{aln_prefix}.gtf ~{aln_prefix}.bed12
        mikado util convert -if gtf -of bed12 ~{ref_prefix}.gtf ~{ref_prefix}.bed12

        create_mgc_groups -f ~{aln_prefix}.cdna.fa

        cat ~{aln_prefix}.cdna.fa ~{ref_prefix}.cdna.fa > all_cdnas.fa
        cat ~{aln_prefix}.bed12 ~{ref_prefix}.bed12 > all_cdnas.bed
    >>>

    runtime {
        cpu: 1
        memory: "8 GB"
        maxRetries: 1
    }

}

task ScoreAlignments {
    input {
        File cdnas
        File bed
        File groups
        String aln_prefix
        String ref_prefix
        RuntimeAttr? runtime_attr_override
    }


    output {
        File alignment_compare = "comp_"+aln_prefix+"_"+ref_prefix+".tab"
        File alignment_compare_detail = "comp_"+aln_prefix+"_"+ref_prefix+"_detail.tab"
    }

    Int cpus = 6
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 32,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int task_cpus = runtime_attr.cpu_cores

    # TODO:
    # - Add mgc results to alignments

    command <<<
        set -euxo pipefail
        multi_genome_compare.py -t ~{task_cpus} --groups ~{groups} --cdnas ~{cdnas} --bed12 ~{bed} -o comp_~{aln_prefix}_~{ref_prefix}.tab -d comp_~{aln_prefix}_~{ref_prefix}_detail.tab
    >>>

    runtime {
        cpu: task_cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

}