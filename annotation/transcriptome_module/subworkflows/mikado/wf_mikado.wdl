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
        File scoring_file
        File? orf_calling_proteins
        File? homology_proteins
        String output_prefix
        RuntimeAttr? orf_calling_resources
        RuntimeAttr? orf_protein_index_resources
        RuntimeAttr? orf_protein_alignment_resources
        RuntimeAttr? homology_alignment_resources
        RuntimeAttr? homology_index_resources
        File? extra_config
        File? junctions
        Int prodigal_gencode = 1
        String transdecoder_genetic_code = "universal"
        String? orf_caller = "prodigal"
        String transdecoder_alignment_program = "blast"
        String homology_alignment_program = "diamond"
        Boolean mikado_do_homology_assessment
    }

    if (defined(SR_assemblies)) {
        Array[AssembledSample] def_SR_assemblies = select_first([SR_assemblies])
        scatter (sr_assembly in def_SR_assemblies) {
            call GenerateModelsList as sr_models {
                input:
                assembly = sr_assembly
            }
        }
    }

    if (defined(LQ_assemblies)) {
        Array[AssembledSample] def_LQ_assemblies = select_first([LQ_assemblies])

        scatter (lr_assembly in def_LQ_assemblies) {
            call GenerateModelsList as LQ_models {
                input:
                assembly = lr_assembly,
                long_score_bias = 1
            }
        }
    }

    if (defined(HQ_assemblies)) {
        Array[AssembledSample] def_HQ_assemblies = select_first([HQ_assemblies])

        scatter (lr_assembly in def_HQ_assemblies) {
            call GenerateModelsList as HQ_models {
                input:
                assembly = lr_assembly,
                long_score_bias = 1
            }
        }
    }

    call WriteModelsFile {
        input:
        models = flatten(select_all([sr_models.models, LQ_models.models, HQ_models.models]))
    }

    call MikadoPrepare {
        input:
        reference_fasta = indexed_reference.fasta,
        blast_targets = homology_proteins,
        models = WriteModelsFile.result,
        scoring_file = scoring_file,
        extra_config = extra_config,
        output_prefix = output_prefix
    }

    # ORF Calling
    if (defined(orf_caller)) {
        String def_orf_caller = select_first([orf_caller])
        if (def_orf_caller == "prodigal") {
            call pdg.wf_prodigal {
                input:
                gencode = prodigal_gencode,
                prepared_transcripts = MikadoPrepare.prepared_fasta,
                output_directory = output_prefix,
                prodigal_runtime_attr = orf_calling_resources
            }
        }

        if (def_orf_caller == "GTCDS") {
            call GTCDS {
                input:
                prepared_transcripts = MikadoPrepare.prepared_fasta,
                gtf = MikadoPrepare.prepared_gtf,
                runtime_attr_override = orf_calling_resources
            }
        }

        if (def_orf_caller == "transdecoder") {
            call tdc.wf_transdecoder as Transdecoder {
                input:
                prepared_transcripts = MikadoPrepare.prepared_fasta,
                orf_proteins = orf_calling_proteins,
                orf_calling_resources = orf_calling_resources,
                genetic_code = transdecoder_genetic_code,
                index_resources = orf_protein_index_resources,
                output_prefix = output_prefix,
                orf_alignment_resources = orf_protein_alignment_resources,
                program = transdecoder_alignment_program
            }
        }

        File maybe_orfs = select_first([wf_prodigal.orfs, GTCDS.orfs, Transdecoder.gff])
    }

    # Mikado Homology
    if (mikado_do_homology_assessment && defined(homology_proteins)) {
        call hml.wf_homology as Homology {
            input:
            homology_alignment_program = homology_alignment_program,
            reference = MikadoPrepare.prepared_fasta,
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
        config = MikadoPrepare.mikado_config
    }

    call MikadoPick {
        input:
        config_file = MikadoPrepare.mikado_config,
        mikado_db = MikadoSerialise.mikado_db,
        transcripts = MikadoPrepare.prepared_gtf,
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
        File result = "models.txt"
    }
# awk 'BEGIN{OFS="\t"} {$1=$1} 1' > models.txt transforms inputs looking like:
# "Ara.hisat.stringtie.gtf	Ara.hisat.stringtie	True	0" "Ara.hisat.scallop.gtf	Ara.hisat.scallop	True	0"
# into a proper models file, this is a task because using write_lines(flatten([])) did not work properly in the HPC
    command <<<
        set -euxo pipefail
        for i in "~{sep="\" \"" models}"; do
        echo $i; done | awk 'BEGIN{OFS="\t"} {$1=$1} 1' > models.txt;
    >>>
}

### # \"" Fix the parser

task MikadoPick {
    input {
        File config_file
        File transcripts
        File mikado_db
#mode options = ("permissive", "stringent", "nosplit", "split", "lenient")
        String mode = "permissive"
        Int flank = 200
        String output_prefix
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8

    output {
        File index_log  = output_prefix + "-index_loci.log"
        File loci_index = output_prefix + "-" +mode+".loci.gff3.midx"
        File loci       = output_prefix + "-" + mode + ".loci.gff3"
        File scores     = output_prefix + "-" + mode + ".loci.scores.tsv"
        File metrics    = output_prefix + "-" + mode + ".loci.metrics.tsv"
        File stats      = output_prefix + "-" + mode + ".loci.gff3.stats"
        File subloci    = output_prefix + "-" + mode + ".subloci.gff3"
        File monoloci   = output_prefix + "-" + mode + ".monoloci.gff3"
    }

    command <<<
    set -euxo pipefail
    export TMPDIR=/tmp
    mikado pick ~{"--source Mikado_" + mode} ~{"--mode " + mode} --procs=~{cpus} \
    ~{"--flank " + flank} --start-method=spawn ~{"--json-conf=" + config_file} \
    --loci-out ~{output_prefix}-~{mode}.loci.gff3 -lv INFO ~{"-db " + mikado_db} \
    --subloci-out ~{output_prefix}-~{mode}.subloci.gff3 --monoloci-out ~{output_prefix}-~{mode}.monoloci.gff3 \
    ~{transcripts}
    mikado compare -r ~{output_prefix}-~{mode}.loci.gff3 -l ~{output_prefix}-index_loci.log --index
    mikado util stats  ~{output_prefix}-~{mode}.loci.gff3 ~{output_prefix}-~{mode}.loci.gff3.stats
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 16,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MikadoSerialise {
    input {
        IndexedReference indexed_reference
        File config
        File transcripts
        Array[File]? homology_alignments
        File? clean_seqs_db
        File? junctions
        File? orfs
        RuntimeAttr? runtime_attr_override
    }
    Int cpus = 8

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
    mikado serialise ~{xml_prefix}~{sep="," homology_alignments} ~{"--blast_targets="+clean_seqs_db} ~{"--junctions="+junctions} ~{"--orfs="+orfs} \
    ~{"--transcripts=" + transcripts} --genome_fai=${fai} \
    ~{"--json-conf=" + config} --force --start-method=spawn -od mikado_serialise --procs=~{cpus}
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 8,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task GTCDS {
    input {
        File prepared_transcripts
        File gtf
        RuntimeAttr? runtime_attr_override
    }

    output {
        File orfs = "mikado_prepared.gt_cds.trans.bed12"
        File gff3 = "mikado_prepared.gt_cds.gff3"
    }

    command <<<
        set -euxo pipefail
        awk '$3!~\"(CDS|UTR)\"' ~{gtf} \
        | mikado util convert -if gtf -of gff3 - \
        | gt gff3 -tidy -retainids -addids | gt cds -seqfile ~{prepared_transcripts} - matchdesc \
        | gff3_name_to_id.py - mikado_prepared.gt_cds.gff3 && \
        mikado util convert -t -of bed12 "mikado_prepared.gt_cds.gff3" "mikado_prepared.gt_cds.trans.bed12"
    >>>

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
        max_retries: 1
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
        echo -e "~{assembly.assembly}\t~{assembly.name}\t${strand}\t~{long_score_bias}"
    >>>
}

task MikadoPrepare {
    input {
        File models
        File reference_fasta
        File? blast_targets
        File? scoring_file
        File? extra_config
        String output_prefix
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8

    output {
        File mikado_config = output_prefix+"-mikado.yaml"
        File prepared_fasta = output_prefix+"-mikado_prepare/mikado_prepared.fasta"
        File prepared_gtf = output_prefix+"-mikado_prepare/mikado_prepared.gtf"
    }

    command <<<
        set -euxo pipefail
        mikado configure \
        ~{"--scoring=" + scoring_file} \
        ~{"-bt " + blast_targets} \
        --list=~{models} \
        ~{"--reference=" + reference_fasta} \
        ~{output_prefix}-mikado.yaml

        # Merge special configuration file for this run here
        yaml-merge ~{output_prefix}-mikado.yaml ~{extra_config}
        mikado prepare --procs=~{cpus} --json-conf=~{output_prefix}-mikado.yaml -od ~{output_prefix}-mikado_prepare --strip_cds
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 8,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
