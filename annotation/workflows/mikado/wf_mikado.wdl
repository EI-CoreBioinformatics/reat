version 1.0
workflow wf_mikado {
    input {
        Array[File] assemblies
        Array[File]? long_assemblies
        File? reference_annotation
        File? scoring_file
    }

    call GenerateModelsList {
        input:
        assemblies = assemblies,
        long_assemblies = select_first([long_assemblies])
    }

    call MikadoConfigure {
        input:
        reference_annotation = reference_annotation,
        models = GenerateModelsList.models,
        scoring_file = scoring_file
    }

    output {
        File mikado_config = MikadoConfigure.mikado_config
    }
}

task GenerateModelsList {
    input {
        Array[File] assemblies
        Array[File] long_assemblies
    }

    output {
        File models = "models_list.txt"
    }

    command <<<
        cat ~{write_lines(assemblies)} > models_list.txt
        cat ~{write_lines(long_assemblies)} >> models_list.txt
    >>>
}

task MikadoConfigure {
    input {
        File? reference_annotation
        File models
        File? scoring_file
    }

    output {
        File mikado_config = "mikado.yaml"
    }

    command <<<
        mikado configure \
        ~{"--scoring=" + scoring_file} \
        --list=~{models} \
        ~{"--reference=" + reference_annotation} \
        mikado.yaml
    >>>
}