# Biopython-project:
# Project Title
In-Silico Functional Annotation of PET46: A Putative Thermostable PET-Degrading Enzyme from Candidatus Bathyarchaeota
# Objective
The objective of this project is to computationally characterize and functionally annotate an uncharacterized protein (PET46) using bioinformatics approaches.
The study aims to:

Analyze physicochemical properties of the protein

Identify homologous proteins using BLASTp

Predict protein family, domains, and catalytic residues

Assess structural features and subcellular localization

Evaluate the industrial and environmental potential of the protein in plastic biodegradation
# Input Data

Protein ID: A0ACD6BAI4

Protein Name: Alpha/Beta Hydrolase (PET46)

Source Organism: Candidatus Bathyarchaeota archaeon

Sequence Length: 262 amino acids

Data Format: FASTA

Data Source: NCBI Protein Database (MAG sequence)
# Methodology

1.Sequence Retrieval

Protein sequence was retrieved in FASTA format from NCBI.

2.Physicochemical Characterization

Molecular weight, isoelectric point, instability index, and amino acid composition were computed using Biopython’s ProtParam module.

3.Homology Search (BLASTp)

BLASTp was performed against the NCBI non-redundant (nr) database to identify homologous proteins and infer function.

4.Domain and Motif Analysis

InterProScan results were used to identify conserved domains and catalytic residues.

5.Structural Prediction

AlphaFold structural model was analyzed to identify fold type, lid domain, and confidence score.

6.Subcellular Localization Prediction

Localization was predicted to assess functional relevance in polymer degradation.

7.Functional Interpretation

Combined evidence from sequence similarity, domain architecture, and structural features was used for functional annotation.

# Results
1.Physicochemical Properties

Molecular Weight: 29.43 kDa

Isoelectric Point (pI): 5.66

Instability Index: 41.68 (moderately stable)

Dominant amino acids: Leucine, Glycine, Valine, Arginine, Aspartate

2.BLASTp Analysis

Top hit: pdb|8B4U|A

Identity: 262/262 (100%)

E-value: 0.0

Homologs identified as alpha/beta hydrolase family proteins

3.Domain Analysis

Identified domain: PF00756 (Esterase/Lipase family)

Superfamily: Alpha/Beta Hydrolase Fold

Catalytic triad: Ser–His–Asp

4.Structural Features

High-confidence AlphaFold model (pLDDT > 90)

Presence of a lid domain covering the active site

Typical fold associated with ester bond hydrolysis

5.Subcellular Localization

Predicted localization: Extracellular / Secreted

Probability: 0.89

Supports interaction with solid plastic substrates

The present study successfully performed an in-silico functional annotation of the protein PET46 from Candidatus Bathyarchaeota using a multi-tool bioinformatics pipeline. Physicochemical analysis revealed that PET46 is a moderately stable, acidic protein with a molecular weight of 29.43 kDa, consistent with enzymes belonging to the alpha/beta hydrolase superfamily.

Although the automated sequence quality control step indicated insufficient quality, manual verification and multi-database annotation confirmed the protein’s identity as a feruloyl esterase / alpha-beta hydrolase. Domain analysis identified the conserved PF00756 esterase/lipase domain, while AlphaFold structural prediction demonstrated a high-confidence model (pLDDT > 90) with a characteristic lid domain covering the active site, a feature commonly associated with ester bond hydrolysis.

Subcellular localization analysis predicted PET46 to be extracellular and secreted, supporting its functional role in interacting with polymeric substrates such as PET. The deep-sea origin of the source organism suggests thermal adaptation, further strengthening the enzyme’s suitability for high-temperature industrial processes.

Overall, the integrated computational evidence supports the conclusion that PET46 is a heat-adapted alpha/beta hydrolase with strong potential for industrial PET biodegradation at temperatures around 70 °C. Experimental validation is recommended to confirm enzymatic activity and assess its efficiency in real-world plastic recycling applications.

# Conclusion

This in-silico study confirms that PET46 is a member of the alpha/beta hydrolase superfamily, closely related to microbial esterases.
Structural and domain analyses suggest that PET46 functions as a thermostable feruloyl esterase capable of hydrolyzing ester bonds in PET polymers.
The extracellular nature and predicted thermal stability make PET46 a promising candidate for industrial plastic biodegradation and environmental remediation.
Further experimental validation is required to confirm PET-degrading activity.


# Tools and Software Used
1.Tool / Software	Purpose: Python 3.14	Programming environment

2.Biopython	Sequence handling and analysis:BLASTp	Homology and functional inference

3.InterProScan	Domain and motif identification

4.ProtParam (Biopython)	Physicochemical analysis

5.NCBI Protein Database	Sequence retrieval

Author: Namratha H

Submitting Individually




course work project using basics of Bio-python
