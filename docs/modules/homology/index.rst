Homology Workflow
-----------------

.. highlight:: none
.. include:: ./homology_help.txt
  :literal:


When there is protein evidence available from related species, the homology workflow can be used to generate gene models based on this evidence. This is achieved by aligning the proteins provided through a set of related species annotations and evaluating these alignments to generate a score.

After the proteins from related species are aligned to the reference, these alignments are filtered using a set of selectable criteria such as: exon length, intron length, protein length, canonical splicing and presence of internal stop codons. The resulting sequences are then evaluated.

Protein alignments are evaluated in two ways, coherence of the alignment structure with respect to the original model's structure and consensus structure from the multiple species. These scores are then used by Mikado to group and filter models, generating a set of predicted models.

In the following example we observe two cases of the comparison between alignments and the original structure, the example on the top of the image shows the original gene structure matching that of the protein alignment structure. On the other hand, on the bottom of the image the original structure does not match the aligned protein, specifically the first and second exons have different lengths.

.. image:: /_static/homology/alignment_vs_original.png
  :alt: Comparison of original structure vs alignment result

When comparing alignments from multiple species we evaluate the consensus structure and provide a score corresponding to the most commonly observed structure. In the following example, we observe a clear majority of protein alignments supporting a structure with five exons, this leads to a **xspecies** score of 70% for the models marked correct.

.. image:: /_static/homology/cross_species_concordance.png
  :alt: Comparison alignment structure from multiple species

Configurable computational resources available::

  "ei_homology.CombineXspecies.runtime_attr_override": " {
                 cpu_cores -> Int
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_homology.aln_attr": " {
                 cpu_cores -> Int
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_homology.index_attr": " {
                 cpu_cores -> Int
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_homology.mikado_attr": " {
                 cpu_cores -> Int
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_homology.score_attr": " {
                 cpu_cores -> Int
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)"


.. image:: /_static/homology_workflow.png
  :alt: Homology workflow diagram
