Prediction Workflow
========================

The intention of the prediction workflow is to use a variety of transcript evidence, from short reads and long reads based gene assemblies, protein alignments, homology alignments and other evidence such as expression, introns and repeats to generate gene predictions ab initio and evidence based gene predictions.

.. highlight:: none
.. include:: ./prediction_help.txt
  :literal:


The prediction module takes as input a genome file along with a set of evidences for annotations over the genome, these can come from homology protein or transcript alignments, rna-seq gene models, repeat annotations, rna-seq alignments which can provide evidence to the presence/absence of exons. Also, the user should provide a set of proteins to validate against, these proteins are used to score input models, categorize them into Gold, Silver or Bronze and select the best models for training of the ab initio gene predictors.

Input models from homology or transcriptomic sources are aligned to a protein database of the user's choice and the results of these alignments are used to classify and score each input model into Bronze, Silver and Gold. Models from the Gold and Silver category are defined by:

- Having complete but not excessively long UTR's.
- Being fully covered by multiple proteins from the database.
- Having a long enough CDS, where the length is user defined.

For scoring the models, a score is calculated for the following properties: the distance between the start and end of the model and the target protein are compared to the start and end of the alignment; the coverage of the model and the target protein; and the length of the longest gap in the alignment. The score is defined by three parameters that are user controlled, a 'hard filter' after which the criteria is considered failed, a 'soft filter' from where alignments receive a score relative to the difference between the 'best' possible value and the current level, and finally the maximum possible score.

Once models have been scored, a user defined number of models at a user defined ratio between mono-exonic and multi-exonic are randomly selected to train ab initio predictors. These models are selected from the classified models ordered by 'category' (Gold, Silver, Bronze, others in this order) and score (highest to lowest).

Each of the ab initio predictors the user selected is then trained and used to generate predictions. In the case of Augustus, there is an initial ab initio prediction made with limited evidence, but further rounds of prediction with different weights for each evidence type can then be configured using a file containing a SOURCE and a SCORE value for each criteria (more details). All these predictions are then combined using Evidence Modeler with configurable weights for each type of prediction and evidence (more details). Finally, the EVM output is processed through Mikado using the Gold and Silver category models (which contain UTRs) to add UTRs where evidence supports it.


Configurable computational resources available::

  "ei_prediction.AlignProteins.resources": " {
                 cpu_cores -> Int
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_prediction.Augustus.resources": " {
                 cpu_cores -> Int
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_prediction.AugustusAbinitio.resources": " {
                 cpu_cores -> Int
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_prediction.ExecuteEVMCommand.resources": " {
                 cpu_cores -> Int
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_prediction.IndexProteinsDatabase.resources": " {
                 cpu_cores -> Int
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_prediction.LengthChecker.resources": " {
                 cpu_cores -> Int
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_prediction.Mikado.resources": " {
                 cpu_cores -> Int
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_prediction.MikadoPick.resources": " {
                 cpu_cores -> Int
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_prediction.SelfBlastFilter.resources": " {
                 cpu_cores -> Int
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)"

.. image:: /_static/prediction_workflow.png
  :alt: Prediction workflow diagram