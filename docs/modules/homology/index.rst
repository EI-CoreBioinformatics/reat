Homology workflow
-----------------

When there is protein evidence available from related species, the homology workflow can be used to generate gene models based on this evidence. This is achieved by aligning the proteins provided through a set of related species annotations and evaluating these alignments to generate a score.

Protein alignments are evaluated in two ways, coherence of the alignment structure with respect to the original model's structure and consensus structure from the multiple species. These scores are then used by Mikado to group and filter models, generating a set of predicted models.

.. highlight:: none
.. include:: ./homology_help.txt
  :literal:

.. image:: /_static/homology_workflow.png
  :alt: Homology workflow diagram
