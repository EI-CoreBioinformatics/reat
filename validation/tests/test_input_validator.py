from unittest import TestCase
from validation.input_validator import *


class InputValidatorTest(TestCase):
    def test_valid_input(self):
        valid_example = json.loads("""
            {
                "ei_annotation.paired_samples": [
                    {
                        "name": "Ara",
                        "strand": "fr-firststrand",
                        "read_pair": [
                            {
                                "R1": "inputs/reads/Ara1.1.fastq",
                                "R2": "inputs/reads/Ara1.2.fastq"
                            },
                            {
                                "R1": "inputs/reads/Ara2.1.fastq",
                                "R2": "inputs/reads/Ara2.2.fastq"
                            },
                            {
                                "R1": "inputs/reads/Ara3.1.fastq",
                                "R2": "inputs/reads/Ara3.2.fastq"
                            },
                            {
                                "R1": "inputs/reads/Ara5.1.fastq",
                                "R2": "inputs/reads/Ara5.2.fastq"
                            },
                            {
                                "R1": "inputs/reads/Ara6.1.fastq",
                                "R2": "inputs/reads/Ara6.2.fastq"
                            }
                        ]
                    }
                ],
                "ei_annotation.LQ_long_read_samples": [
                    {
                        "name": "A01_1",
                        "strand": "fr-firststrand",
                        "LR": ["inputs/reads/A01_1.fastq"]
                    },
                    {
                        "name": "A01_2",
                        "strand": "fr-firststrand",
                        "LR": ["inputs/reads/A01_2.fastq"]
                    },
                    {
                        "name": "B01",
                        "strand": "fr-firststrand",
                        "LR": ["inputs/reads/B01.fastq"]
                    },
                    {
                        "name": "C01",
                        "strand": "fr-firststrand",
                        "LR": ["inputs/reads/C01.fastq"]
                    }
                ],
                "ei_annotation.wf_align.LQ_assembler": "stringtie",
                "ei_annotation.HQ_long_read_samples": [
                    {
                        "name": "CCS",
                        "strand": "fr-firststrand",
                        "LR": ["inputs/reads/CCS.fastq"]
                    },
                    {
                        "name": "polished",
                        "strand": "fr-firststrand",
                        "LR": ["inputs/reads/polished.fastq"]
                    }
                ],
                "ei_annotation.wf_align.HQ_assembler": "merge",
                "ei_annotation.reference_genome": "inputs/reference/Myzus_persicae_O_v2.0.scaffolds.fa",
                "ei_annotation.mikado_scoring_file": "inputs/plant.yaml",
                "ei_annotation.homology_proteins": "inputs/proteins_db/cross_species_all.protein_nl70.fasta",
                "ei_annotation.orf_calling_proteins": "inputs/proteins_db/Rmai.fa",
                "ei_annotation.wf_main_mikado.orf_calling_program": "Transdecoder"
            }
            """)
        self.assertIsNone(validate_reat_inputs(valid_example))

    def test_invalid_input(self):
        invalid_example = {
            "name": "Test", "other": "invalid"
        }
        with self.assertRaises(ValidationError):
            validate_reat_inputs(invalid_example)
