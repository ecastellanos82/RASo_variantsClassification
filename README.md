# RASo_variantsClassification

Genetic diagnosis has been an expanding field since the beginning of the century, and especially in the last decade due to the development of new sequencing tools, also known as Next Generation Sequencing or NGS.

Currently, most genetic diagnostics are performed by sequencing multiple genes from multiple patients simultaneously. The implementation of the NGS in the genetic diagnosis has led to multiple improvements in cost-effectiveness, increased knowledge, etc., but also an increase in the number of diagnoses in patients that d not fullfil clinical criteria and the analysis of genes that are so far unknown. This has led to an increase in the detection of variants of uncertain significance, or VUS, where their relation to a specific pathology cannot be confirmed or ruled out. That is why the scientific community is developing clear and systematic guides to interpret the pathogenicity of all variants detected in genetic tests by using different parameters for each of the groups of pathologies analyzed.

RASopathies are a clear example of the need for these more specific guides. This group of diseases show genetic heterogeneity, clinical overlap between them, and almos all variants detected in most genes are gain of function missenses (GOF) and not truncating mutations, which hinders the interpretation of their pathogenicity. For this reason, we propose to develop a script to help us to semi-automatyze its classification.

Implementing these guides manually, however, is quite tedious and expensive.
Nowadays, the diagnosis of RASopathies is being performed in the Clinical Genomics Unit of the Germans Trias & Pujol Hospital, among other centers. This unit shares a relational database and analysis system called Pandora with the ICO Hereditary Cancer Genetic Diagnosis Unit. This system is used for sample management, runs, variants filtering and variant classification.

The objective of this scritp is to extracts data from a relational database (Pandora) to automatically score RASopathies criteria. These are based on: population frequency, domains where the variant is located and the result of predictors in silico. In addition, this script will allow to incorporate non-automatizable criteria (de novo, number of patients already described harboring the same variant, etc.).

For this reason, will adapt this script to be used from a GUI in Shiny to facilitate the use of this variant classification system for RASopathies throughout the staff of the Clinical Genomics Unit
