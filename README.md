# Human RNA-Seq
RNA-Seq Pipeline designed for human data.
Code for publication Franch-Arroyo et al. 2024

# Description
This repository consists of of two components: 1) the analysis pipeline utilizing bash, Python, and R scripts, generating some visualizations, and 2) a single python script creating the scatter plot of the data as it is shown in the publication. The authors of the first part are Timothy Sullivan and Eric Galvez, the author of the second part is Knut Finstermeier.

## Important:
This code has been developed to conduct analyses as outlined in the publication above and is currently in development. The code might contain bugs and is also not written for speed or for readability. Future releases might improve on those aspects.

## Disclaimer:

The software is intended for research purposes only and is provided "as is" without warranty of any kind. The developer(s) and provider(s) of this software make no representations or warranties, express or implied, regarding the use or performance of this software.

You assume full responsibility and risk for the use of this software. The developer(s) and provider(s) shall not be liable for any direct, indirect, incidental, special, exemplary, or consequential damages arising out of the use of this software.

This software is not intended to be used for critical, medical, or life-saving purposes. It should not be relied upon as the sole basis for decision-making.

By using this software, you acknowledge and agree to this disclaimer. If you do not agree to these terms, you should refrain from using the software.

Last updated: 2024-07-05

## Copyright:
Code has been created by Timothy Sullivan, Eric Galvez, and Knut Finstermeier, Max Planck Unit for the Science of Pathogens, Berlin, Germany, member of the Max Planck Society, https://www.mpg.de/en. The code is licensed under the MIT license (see LICENSE file).


# Setup
- clone repository  
```
git clone https://github.com/MPUSP/human_RNA_seq_pipeline.git
```
## Part 1
- Set paths to reference, annotation in wrapper.sh
  - Edit wrapper.sh to add paths to a STAR index, a GTF file, and your python executable (3 total edits in file)
- Sample fastq paths and metadata
  - Create TSV files for fastq paths and sample metadata (see example folder for format)
## Part 2
- Run setup_pt2.sh
- If needed, adjust the paths in the python script. By default, the script refers to the result tables from the original publication

# Running
## Part 1
### Wrapper
- wrapper.sh contains the commands for the full pipeline.  Once the paths at the top are set, running wrapper.sh is all that is needed.
```
./wrapper.sh
```
### Alignment
- Alignment is performed with STAR.  The parameters are copied from the Encode/GTEx pipelines, so that data can be most comparable with public datasets
### RNA-SeQC
- RNA-SeQC counts reads that overlap features in the GTF file.  It also computes alignment metrics for each sample.  A separate script joins the per-sample results into a single table, for read counts and read TPMs
### DESeq
- DESeq pre-processing calculates size factors and surrogate variables for the entire study, and saves the DESeq object to a file.  r-log transformed reads are also saved to a file for downstream visualization.
- DESeq post-processing performs the differential expression calculation for the specified comparison.  It also generates geneset analysis and pathway diagrams, and outputs the DESeq results table of p-values and log2 fold changes.
### Bokeh
- The bokeh script creates interactive MA and Volcano plots given a DESeq results table containing adjusted p-values and log2 fold changes.

## Part 2
- after adjusting the paths in the python script, run the script
```
./venv/bin/python create_plot.v2.py
```

## Contact:
Please use software-dev@mpusp.mpg.de for any questions or comments.