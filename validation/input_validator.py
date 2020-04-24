import argparse
import json
from jsonschema import ValidationError, validators, Draft7Validator


def is_valid_name(validator, value, instance, schema):
    if not isinstance(instance, str):
        yield ValidationError("%r is not a string" % instance)

    if value and set(instance).intersection("][/?\\\'\" .*$)(}{"):
        yield ValidationError("%r is not alphanumeric" % (instance,))


reat_schema = {
    "$schema": "http://json-schema.org/draft-07/schema#",

    "definitions": {
        "Boolean": {
            "enum": ["true", "false"]
        },

        "strand": {
            "enum": ["fr-firststrand", "fr-unstranded", "fr-secondstrand"]
        },

        "LR_Aligner": {
            "enum": ["minimap2", "gmap"]
        },

        "LR_Sample": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "name": {
                        "type": "string",
                        "is_name": True,
                    },
                    "strand": {
                        "$ref": "#/definitions/strand"
                    },
                    "score": {
                        "type": "number"
                    },
                    "LR": {
                        "type": "array",
                        "items": {
                            "type": "string",
                        }
                    },
                },
                "required": ["name", "strand", "LR"]
            }
        },

        "LR_Assembler": {
            "enum": ["filter", "merge", "stringtie", "stringtie_collapse"]
        },
    },

    "type": "object",
    "properties": {
        "ei_annotation.reference_genome": {
            "type": "string"
        },
        "ei_annotation.paired_samples": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "name": {
                        "type": "string",
                        "is_name": True,
                    },
                    "strand": {
                        "$ref": "#/definitions/strand"
                    },
                    "score": {
                        "type": "number"
                    },
                    "read_pair": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "R1": {
                                    "type": "string"
                                },
                                "R2": {
                                    "type": "string"
                                }
                            },
                            "required": ["R1", "R2"]
                        }
                    },
                },
                "required": ["name", "strand", "read_pair"]
            }
        },
        "ei_annotation.HQ_long_read_samples": {
            "$ref": "#/definitions/LR_Sample"
        },
        "ei_annotation.LQ_long_read_samples": {
            "$ref": "#/definitions/LR_Sample"
        },
        "ei_annotation.annotation": {
            "type": "string"
        },
        "ei_annotation.homology_proteins": {
            "type": "string"
        },
        "ei_annotation.mikado_scoring_file": {
            "type": "string"
        },
        "ei_annotation.orf_calling_proteins": {
            "type": "string"
        },
        "ei_annotation.wf_align.HQ_aligner": {
            "$ref": "#/definitions/LR_Aligner"
        },
        "ei_annotation.wf_align.LQ_aligner": {
            "$ref": "#/definitions/LR_Aligner"
        },
        "ei_annotation.wf_align.HQ_assembler": {
            "$ref": "#/definitions/LR_Assembler"
        },
        "ei_annotation.wf_align.LQ_assembler": {
            "$ref": "#/definitions/LR_Assembler"
        },
        "ei_annotation.wf_main_mikado.annotation_bed": {
            "type": "string"
        },
        "ei_annotation.wf_main_mikado.orf_calling_program": {
            "enum": ["Prodigal", "Transdecoder", "None"]
        },
        "ei_annotation.wf_main_mikado.run_mikado_homology": {
            "$ref": "#/defitions/Boolean"
        },
        "ei_annotation.wf_main_mikado.separate_LQ": {
            "$ref": "#/defitions/Boolean"
        },
    },
    "required": ["ei_annotation.reference_genome"],
    "anyOf": [
        {"required": ["ei_annotation.paired_samples"], },
        {"required": ["ei_annotation.HQ_long_read_samples"], },
        {"required": ["ei_annotation.LQ_long_read_samples"], },
    ]
}


def validate_reat_inputs(inputs):
    all_validators = dict(Draft7Validator.VALIDATORS)
    all_validators["is_name"] = is_valid_name
    reat_validator = validators.create(meta_schema=reat_schema, validators=all_validators)
    reat_validator(reat_schema).validate(inputs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="REAT JSON inputs Validator")
    parser.add_argument("-i", "--input_file")

    args = parser.parse_args()

    try:
        reat_inputs = json.load(args.input_file)
    except json.JSONDecodeError:
        print("Invalid JSON file")

    validate_reat_inputs(reat_inputs)