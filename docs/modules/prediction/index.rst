Prediction Workflow
========================

The intention of the prediction workflow is to use a variety of transcript evidence, from short reads and long reads based gene assemblies, protein alignments, homology alignments and other evidence such as expression, introns and repeats to generate gene predictions ab initio and evidence based gene predictions.

.. highlight:: none
.. include:: ./prediction_help.txt
  :literal:


The prediction module takes as input a genome file along with a set of evidences for annotations over the genome (these should have gene->mrna->{exon,CDS} structure, where CDS is required for protein inputs), these can come from homology proteins or transcript alignments, rna-seq gene models, repeat annotations, rna-seq alignments which can provide evidence to the presence/absence of exons. Also, the user should provide a set of proteins to validate against, these proteins are used to score input models, categorize them into Gold, Silver or Bronze and select the best models for training of the ab initio gene predictors.

Multiple sets of input models from homology proteins or transcriptomic sources are aligned to a protein database of the user's choice and the results of these alignments are used to classify and score each input model into Bronze, Silver and Gold. Models from the Gold and Silver category are defined by:

- Having complete but not excessively long UTR's.
- Being fully covered by multiple proteins from the database.
- Having a long enough CDS, where the length is user defined.

For scoring the models, a score is calculated for the following properties: the distance between the start and end of the model and the target protein are compared to the start and end of the alignment; the coverage of the model and the target protein; and the length of the longest gap in the alignment. The score is defined by three parameters that are user controlled, a 'hard filter' after which the criteria is considered failed, a 'soft filter' from where alignments receive a score relative to the difference between the 'best' possible value and the current level, and finally the maximum possible score.

Once models have been scored, models with more than a coverage and identity user defined threshold (80% by default) are filtered. From the similarity filtered models, a user defined number of models at a user defined ratio between mono-exonic and multi-exonic are randomly selected to train ab initio predictors. These models are selected from the classified models ordered by 'category' (Gold, Silver, Bronze, others in this order) and score (highest to lowest).

Each of the ab initio predictors the user selected is then trained and used to generate predictions. In the case of Augustus, there is an initial ab initio prediction made with limited evidence, but further rounds of prediction with different weights for each evidence type can then be configured using a file containing a SOURCE and a SCORE value for each criteria (:ref:`see <augustus-runs>`). These parameters depend on the extrinsic information configuration file used by Augustus, for more information about REAT's default :ref:`see the following section <augustus-configuration>`. All these predictions are then combined using Evidence Modeler with configurable weights for each type of prediction and evidence (:ref:`see <evm-weights>`). Finally, the EVM output is processed through Mikado using the Gold and Silver category models (which contain UTRs) to add UTRs where evidence supports it.

.. _augustus-runs:

Configuring Augustus runs
---------------------------

When generating predictions using Augustus, we need to choose the weight parameters for each type of evidence, whilst at the same time possibly wanting to have multiple options of weight sets and priorities as to predict a comprehensive set of models that will maximise our chances of predicting correct structures. In REAT we can decide the number of Augustus predictions and the weights for each prediction using a configuration file per prediction. This file contains a pair of SOURCE and SCORE for each of the evidence types available, which are: gold models, silver models, bronze models, all models, gold introns, silver introns, protein models, coverage hints, repeat hints, high quality assemblies, low quality assemblies, high quality proteins, and low quality proteins. Each file provided to the :code:`--augustus_runs` parameter will trigger a run of Augustus using the specific combination of weights and priorities defined for each evidence type, resulting in as many predictions as files provided.

The default Augustus configuration file can be overriden to make available for the user different 'SOURCE's which can then be used for the :code:`--augustus_runs` files, the following is an example of a 'run' file::

	M 10
	F 9
	E 8
	E 7
	E 6
	E 4
	P 4
	W 3
	RM 1
	E 2
	E 2
	E 2
	E 2


.. _augustus-configuration:

Extrinsic information configuration file
-----------------------------------------

.. include:: ./extrinsic.ei_augustus_generic.cfg
  :literal:


.. _evm-weights:

Evidence Modeler default weights file
--------------------------------------
.. include:: ./evm_default_weights.wgt
  :literal:

Configurable computational resources available
-----------------------------------------------
 ::

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