version 1.0

import "../common/structs.wdl"
import "./orf_caller/wf_transdecoder.wdl" as tdc
import "./orf_caller/wf_prodigal.wdl" as pdg
import "./homology/wf_homology.wdl" as hml

workflow wf_mikado {
    input {
        IndexedReference indexed_reference
        Array[AssembledSample]? SR_assemblies
        Array[AssembledSample]? LQ_assemblies
        Array[AssembledSample]? HQ_assemblies

        File? annotation
        Int annotation_score
        String mode
        Boolean check_reference

        File scoring_file
        File? orf_calling_proteins
        File? homology_proteins
        String output_prefix


        RuntimeAttr? orf_calling_resources
        RuntimeAttr? orf_protein_index_resources
        RuntimeAttr? orf_protein_alignment_resources
        RuntimeAttr? homology_alignment_resources
        RuntimeAttr? homology_index_resources
        RuntimeAttr? mikado_pick_resources
        RuntimeAttr? mikado_serialise_resources
        RuntimeAttr? mikado_prepare_resources
        File? prepare_extra_config
        File? serialise_extra_config
        File? pick_extra_config
        File? junctions
        Int prodigal_gencode = 1
        String transdecoder_genetic_code = "universal"
        String? orf_caller = "prodigal"
        String transdecoder_alignment_program = "blast"
        String homology_alignment_program = "diamond"
        Boolean mikado_do_homology_assessment
    }

    Array[Array[AssembledSample]] jsamples = select_all([SR_assemblies, LQ_assemblies, HQ_assemblies])
    call MikadoPrepare {
        input:
        reference_fasta = indexed_reference.fasta,
        blast_targets = homology_proteins,
        mode = mode,
        jsamples = jsamples,
        annotation = annotation,
        annotation_score = annotation_score,
        check_reference = check_reference,
        extra_config = prepare_extra_config,
        runtime_attr_override = mikado_prepare_resources,
        output_prefix = output_prefix
    }

    if (defined(annotation)) {
        call FilterPrepare {
            input:
                prepared_transcripts = MikadoPrepare.prepared_fasta
        }
    }

    # ORF Calling
    if (defined(orf_caller)) {
        String def_orf_caller = select_first([orf_caller])
        if (def_orf_caller == "prodigal") {
            call pdg.wf_prodigal {
                input:
                gencode = prodigal_gencode,
                prepared_transcripts = select_first([FilterPrepare.filtered_prepared_fasta, MikadoPrepare.prepared_fasta]),
                output_directory = output_prefix,
                prodigal_runtime_attr = orf_calling_resources
            }
        }

        if (def_orf_caller == "transdecoder") {
            call tdc.wf_transdecoder as Transdecoder {
                input:
                prepared_transcripts = select_first([FilterPrepare.filtered_prepared_fasta, MikadoPrepare.prepared_fasta]),
                orf_proteins = orf_calling_proteins,
                orf_calling_resources = orf_calling_resources,
                genetic_code = transdecoder_genetic_code,
                index_resources = orf_protein_index_resources,
                output_prefix = output_prefix,
                orf_alignment_resources = orf_protein_alignment_resources,
                program = transdecoder_alignment_program
            }
        }

        File maybe_orfs = select_first([wf_prodigal.orfs, Transdecoder.gff])
    }

    # Mikado Homology
    if (mikado_do_homology_assessment && defined(homology_proteins)) {
        call hml.wf_homology as Homology {
            input:
            homology_alignment_program = homology_alignment_program,
            reference = select_first([FilterPrepare.filtered_prepared_fasta, MikadoPrepare.prepared_fasta]),
            protein_db = homology_proteins,
            index_resources = homology_index_resources,
            protein_alignment_resources = homology_alignment_resources
        }
    }
    
    call MikadoSerialise {
        input:
        homology_alignments = 
        if (mikado_do_homology_assessment && defined(homology_proteins)) then 
            select_first([Homology.homology]) else 
            Homology.homology, # These lines ensure the Homology task completed correctly if it was requested
        clean_seqs_db = Homology.homology_clean_db,
        junctions = junctions,
        orfs = maybe_orfs,
        transcripts = MikadoPrepare.prepared_fasta,
        indexed_reference = indexed_reference,
        config = MikadoPrepare.mikado_config,
        runtime_attr_override = mikado_serialise_resources,
        extra_config = serialise_extra_config
    }

    call MikadoPick {
        input:
        config_file = MikadoPrepare.mikado_config,
        extra_config = pick_extra_config,
        scoring_file = scoring_file,
        mikado_db = MikadoSerialise.mikado_db,
        transcripts = MikadoPrepare.prepared_gtf,
        runtime_attr_override = mikado_pick_resources,
        output_prefix = output_prefix
    }

    output {
        File mikado_config = MikadoPrepare.mikado_config
        File? orfs = maybe_orfs
        Array[File]? homologies = Homology.homology
        Array[File] serialise_out = MikadoSerialise.out
        File loci = MikadoPick.loci
        File scores = MikadoPick.scores
        File metrics = MikadoPick.metrics
        File stats = MikadoPick.stats
    }
}

task WriteModelsFile {
    input {
        Array[String] models
        File? annotation
        Int annotation_score
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }

    output {
        File result = "models.txt"
    }
# awk 'BEGIN{OFS="\t"} {$1=$1} 1' > models.txt transforms inputs looking like:
# "Ara.hisat.stringtie.gtf	Ara.hisat.stringtie	True	0" "Ara.hisat.scallop.gtf	Ara.hisat.scallop	True	0"
# into a proper models file, this is a task because using write_lines(flatten([])) did not work properly in the HPC
    command <<<
        set -euxo pipefail
        for i in "~{sep="\" \"" models}"; do
        echo $i; done | awk 'BEGIN{OFS="\t"} {$1=$1} 1' > models.txt;

        if [ "~{annotation}" != "" ]
        then
            echo -e ~{annotation}'\t'reference'\t'True'\t'~{annotation_score}'\t'True >> models.txt
        fi
    >>>
}

### # \"" Fix the parser

task MikadoPick {
    input {
        File config_file
        File? extra_config
        File scoring_file
        File transcripts
        File mikado_db
        Int flank = 200
        String output_prefix
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 8,
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
    ~{"--flank " + flank} --start-method=spawn --json-conf=pick_config.yaml \
    --loci-out ~{output_prefix}.loci.gff3 -lv INFO ~{"-db " + mikado_db} \
    --subloci-out ~{output_prefix}.subloci.gff3 --monoloci-out ~{output_prefix}.monoloci.gff3 \
    ~{transcripts}
    mikado compare -r ~{output_prefix}.loci.gff3 -l ~{output_prefix}-index_loci.log --index
    mikado util stats  ~{output_prefix}.loci.gff3 ~{output_prefix}.loci.gff3.stats
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }
}

task MikadoSerialise {
    input {
        IndexedReference indexed_reference
        File config
        File? extra_config
        File transcripts
        Array[File]? homology_alignments
        File? clean_seqs_db
        File? junctions
        File? orfs
        RuntimeAttr? runtime_attr_override
    }
    Int cpus = 8
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 8,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    output {
        Array[File] out = glob("mikado_serialise/*")
        File mikado_db = "mikado_serialise/mikado.db"
    }

# The following line is a workaround for having no mechanism to output a "prefix" for an optional array expansion, i.e
# See the usages of xml_prefix
    String xml_prefix = if defined(homology_alignments) then "--xml=" else ""

    command <<<
        set -euxo pipefail
    fasta=~{indexed_reference.fasta}
    fai=~{indexed_reference.fai}
    
    ln -s ${fasta} .
    ln -s ${fai} .
    yaml-merge -s ~{config} ~{"-m " + extra_config} -o serialise_config.yaml
    export TMPDIR=/tmp
    mikado serialise ~{xml_prefix}~{sep="," homology_alignments} ~{"--blast_targets="+clean_seqs_db} ~{"--junctions="+junctions} ~{"--orfs="+orfs} \
    ~{"--transcripts=" + transcripts} --genome_fai=${fai} \
    --json-conf=serialise_config.yaml --force --start-method=spawn -od mikado_serialise --procs=~{task_cpus}
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }
}
    
task GenerateModelsList {
    input {
        AssembledSample assembly
        Int long_score_bias = 0
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1,
        queue: ""
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        String models = read_string(stdout())
    }

    command <<<
        set -euxo pipefail
        strand="False"
        if [ ~{assembly.strand} != "fr-unstranded" ]; then
            strand="True"
        fi
        echo -e "~{assembly.assembly}\t~{assembly.name}\t${strand}\t~{assembly.score}\t~{assembly.is_ref}\t~{assembly.exclude_redundant}"
    >>>
}

task MikadoPrepare {
    input {
        File reference_fasta
        String mode
        Boolean check_reference

        Array[Array[AssembledSample]] jsamples

        File? annotation
        Int annotation_score

        File? blast_targets
        File? extra_config
        String output_prefix
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 8,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    output {
        File mikado_config = output_prefix+"-mikado.yaml"
        File prepared_fasta = output_prefix+"-mikado_prepare/mikado_prepared.fasta"
        File prepared_gtf = output_prefix+"-mikado_prepare/mikado_prepared.gtf"
    }

    command <<<
        set -euxo pipefail
        export TMPDIR=/tmp

        # Choose configure mode

        case "~{mode}" in
                basic)
                    mode_parameter=""
                    ;;

                update)
                    mode_parameter="--reference-update"
                    ;;

                only_update)
                    mode_parameter="--only-reference-update"
                    ;;
        esac


        # Generate models from input SR, LQ and HQ assemblies + annotation/_score

        generate_model_file ~{"--annotation " + annotation} \
        ~{if defined(annotation) then "--annotation_score " + annotation_score else ""} \
        ~{write_json(jsamples)} > model.list

        mikado configure ${mode_parameter} ~{if(check_reference) then "--check-references" else ""} \
        ~{"-bt " + blast_targets} \
        --list=model.list \
        ~{"--reference=" + reference_fasta} \
        ~{output_prefix}-mikado.yaml

        # Merge special configuration file for this run here
        yaml-merge -s ~{output_prefix}-mikado.yaml ~{"-m " + extra_config} -o prepare_config.yaml
        mikado prepare --procs=~{task_cpus} --json-conf=prepare_config.yaml -od ~{output_prefix}-mikado_prepare --strip_cds
    >>>


    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }
}

task FilterPrepare {
    input {
        File prepared_transcripts
    }

    output {
        File filtered_prepared_fasta = "filtered_prepare.fasta"
    }

    command <<<
    bioawk -c'fastx' '{split($name, parts, "_"); if (parts[1] != "reference") print ">"$name"\n"$seq}' ~{prepared_transcripts} > filtered_prepare.fasta
    >>>
}