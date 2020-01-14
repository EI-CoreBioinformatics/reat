workflow wf_mikado {
    Array[File] assemblies
    Array[File?]? long_assemblies
    File? reference_annotation
    File? scoring_file

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
    Array[File] assemblies
    Array[File?]? long_assemblies

    command {
        rm models_list.txt
        for i in ${sep=" " assemblies}; do
          echo -e "$i\tlabel\tstranded" >> models_list.txt
        done;
        for i in ${sep=" " long_assemblies}; do
          echo -e "$i\tlabel\tstranded" >> models_list.txt
        done;
    }

    output {
        File models = "models_list.txt"
    }
}

task MikadoConfigure {
    File? reference_annotation
    File models
    File? scoring_file

    command {
        mikado configure \
        ${"--scoring=" + scoring_file} \
        --list=${models} \
        ${"--reference=" + reference_annotation} \
        mikado.yaml
    }

    output {
        File mikado_config = "mikado.yaml"
    }
}