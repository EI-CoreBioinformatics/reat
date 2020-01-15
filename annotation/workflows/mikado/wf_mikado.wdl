version 1.0
workflow wf_mikado {
    input {
        Array[File] assemblies
        Array[File?]? long_assemblies
        File? reference_annotation
        File? scoring_file
    }

    call GenerateModelsList {
        input:
        assemblies = assemblies,
        long_assemblies = long_assemblies
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
        Array[File?]? long_assemblies
    }

    output {
        File models = "models_list.txt"
    }

    command <<<
        rm models_list.txt
        for i in ~{sep=" " assemblies}; do
          echo -e "$i\tlabel\tstranded" >> models_list.txt
        done;
        for i in ~{sep=" " long_assemblies}; do
          echo -e "$i\tlabel\tstranded" >> models_list.txt
        done;
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