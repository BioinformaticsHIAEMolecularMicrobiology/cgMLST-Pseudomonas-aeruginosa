# cgMLST to *Pseudomonas aeruginosa*

The objective of this repository is to describe how we created a cgMLST for *P. aeruginosa*. This scheme was created with the ChewBBACA pipeline (link below). Input files used to create and validate the scheme as well as output files from each step are available. 

Additionally, the steps to run the scheme and/or download the *P. aeruginosa* selected *loci* that compose the cgMLST scheme are also available in this repository. (Steps 8, 9 and 10)

## chewBBACA

To download ChewBBACA access the link:

* https://github.com/B-UMMI/chewBBACA_tutorial

## Workflow used to create the scheme

* Step 1: Scheme Creation
* Step 2: Allele calling
* Step 3: Scheme Validation (Allele call)
* Step 4: Extracting the Matrix *loci*
* Step 5: Minimum Spanning Tree (MST)
* Step 6: Evaluation of the cgMLST scheme
* Step 7: Analysis of the proteins/genes of the wgMLST
* Step 8: How to identify the allelic profile of the genomes of interest using this cgMLST database for *P. aeruginosa*
* Step 9: Analyze the results
* Step 10: Allele profile view by minimum spanning tree (MST)
 
## Softwares and Downloads (Main dependencies)

* BLAST 2.5.0+ ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.5.0/ or above
* Prodigal 2.6.0 https://github.com/hyattpd/prodigal/releases/ or above

## Other dependencies (for schema evaluation only):

* ClustalW2 http://www.clustal.org/download/current/
* MAFFT https://mafft.cbrc.jp/alignment/software/linux.html
* DATAMASH https://www.gnu.org/software/datamash/
* MLST https://github.com/sanger-pathogens/mlst_check

## Step 1: Schema Creation

## Selection of complete genomes for schema creation

For the selection of *Pseudomonas aeruginosa* genome sequences we consulted RefSeq at https://www.ncbi.nlm.nih.gov/assembly. RefSeq (NCBI Reference Sequence database) corresponds to a comprehensive, non-redundant, well-annotated set of reference sequences. A set of 142 complete genomes sequences of *P. aeruginosa* were publicly available in RefSeq (https://www.ncbi.nlm.nih.gov/assembly) in September 2018. The list of all the complete genomes used to create the scheme can be found in the folder ```Complete_Genomes/Complete_Genomes.xlsx```.

Multilocus sequence type (MLST) for the 142 complete genomes was determined using (https://github.com/sanger-pathogens/mlst_check) and the MLST scheme for *P. aeruginosa* (www.pubmlst.org). The STs obtained for each of the 142 complete genomes using the sanger-pathogens/mlst_check can be found in the folder ```Complete_Genomes/Complete_Genomes.xlsx```. 

Among the 142 genomes, *Pseudomonas aeruginosa* PAO1 reference genome (GCF_000006765.1) was used by Prodigal algorithm as reference to recognize coding sequences (CDs). Prodigal generated the PAO1.trn file at this step. 

**The PAO1 genome was then removed from further analysis**.

## Step 1: Definition of CDs sequences  

## Command: 

```bash
# create schema
chewBBACA.py CreateSchema -i Complete_Genomes --cpu 15 -o schema_seed --ptf PAO1.trn
```

The above command uses 15 CPU and creates a preliminary scheme (wgMLST) in the schema_seed folder using the trained product folder PAO1.trn that was generated using the reference genome PAO1 (GCF_000006765.1) and 141 complete genome sequences. The wgMLST scheme generated contained 13588 *loci* based on the 141 complete genomes. At this point the scheme is defined as a set of candidate *loci*, each containing a single allele. 

Due to the size of the file schema_seed it was not possible to upload it on GitHub, but a link to access the file schema_seed is available at: (https://drive.google.com/open?id=1WpsmbMC0awZ7BH8lazoiv-t56A2xnMIO).

**Note:** Complete_Genomes: Folder containing the list 141 complete genomes that created the schema.

## Step 2: Allele calling

In this step the allele calling is performed using the resulting set of *loci* determined in step 1.

## Command: 

```bash
# run allelecall
chewBBACA.py AlleleCall -i Complete_Genomes -g schema_seed/ -o results_cg --cpu 15 --ptf PAO1.trn
```

The allele calling used the default BLAST Score Ratio (BSR) threshold of 0.6.

## Step 2.1: Paralog detection

In this step *loci* considered paralogous from result of the allelecall (see above) are removed

## Command: 

```bash
# run remove genes
chewBBACA.py RemoveGenes -i results_cg/results_20190921T183955/results_alleles.tsv -g results_cg/results_20190921T183955/RepeatedLoci.txt -o alleleCallMatrix_cg
```
In this step, 282 *loci* were identified as possible paralogs and were removed from further analysis. The list RepeatedLoci.txt can be found at: ```results_cg/results_20190921T183955/RepeatedLoci.txt```

The output file can be found at: ```analysis_cg/alleleCallMatrix_cg.tsv```

## Step 2.2: Genome Quality Control

In this step we define a *Threshold* for the scheme that limits the loss of *loci* targets defined in the previous steps per genome and excludes genomes considered to be of low quality due to significant *loci* absence. 

With this analysis we define the percentage of *loci* that will constitute the scheme based on how many targets we want to keep in this phase. For example,  **100%, 99.5%, 99% and 95%** of the *loci* may present in the set of high quality genomes. This is one of the main steps in defining the cgMLST scheme targets.

## Command:

```bash
# run test genome quality
chewBBACA.py TestGenomeQuality -i alleleCallMatrix_cg.tsv -n 13 -t 300 -s 5
```

In the *Threshold* 120 a set of 3168 *loci* were found to be present in all the analyzed complete genomes, while 4776 *loci* were present in at least 95%. The output file can be found in the folder: ```analysis_cg/GenomeQualityPlot.html```. The list with the genes present in 95% of the genomes at the chosen *Threshold* can be retrieved in the folder ```analysis_cg/Genes_95%.txt```. 

In this stage we chose the *loci* present in 100% (*p1.0*) of the complete genomes and the *Threshold* 120 that limited the loss of the *loci* in the genomes. In this *Threshold* (120) 11 complete genomes were removed due to loss of *loci* . The list of genomes removed at each *Threshold* can be retrieved in the folder ```analysis_cg/removedGenomes.txt```. From this list we created another (GenomeRemoved120thr.txt) with only the genomes removed at *Threshold* (120). The list of genomes removed at *Threshold 120* can be retrieved in the folder: ```analysis_cg/GenomeRemoved120thr.txt```

## Command:

```bash
# run ExtractCgMLST
chewBBACA.py ExtractCgMLST -i alleleCallMatrix_cg.tsv -o cgMLST_120 -p 1.0 -g GenomeRemoved120thr.txt
```

This script selects all * loci * present in the selected * Threshold *. The value * p * is associated with the percentage of * loci * we want to be present in the set of genomes, for example: * p1.0 * selects all * loci * present in the * Threshold * chosen in all genomes ie those present in 100% of genomes at that * Threshold *. Subsequently a cgMLST_120 folder receives the result of the allelic profile for each of the 3168 candidate * loci * (allelic profile matrix). The file in this folder (cgMLST.tsv) contains the allelic profile of 3168 selected * loci * and will be used to create the core gene list. In addition, another 3 output files are created at this point: *cgMLSTschema.txt; mdata_stats.tsv and Presence_Absence.tsv* and they can be found in the folder ```cgMLST_120/```.

## Step 2.3: Creating the core gene list

This command selects all target genes from the "cgMLST.tsv" spreadsheet.

```bash
# 10 list
head -n 1 cgMLST.tsv > Genes_100%_Core_120.txt
```
This step generated the file *Genes_100%_Core_120.txt* that can be retrieved in the folder ```results_cg/Genes_100%_Core_120.txt```. This list needs to be transposed so that each core gene name is reported in a single line:

## Command:

```bash
# transpose table
datamash -W transpose < Genes_100%_Core_120.txt > Genes_Core_Al.txt 
```

This step generated the file > Genes_Core_Al.txt

You can see the list file with 3168 *loci* at ```results_cg/Genes_Core_Al.txt``` and for the subsequent steps we added the full path to each locus fasta file.

This list ```results_cg/Genes_Core_Al.txt``` was then modified so that each name was preceeded by *schema_seed*:

## Command:

```bash
tail -n+1 Genes_Core_Al.txt | cut -f2 | perl -pe 's: :\n:g' | sort -Vu | awk '{print("schema_seed/"$1)}' > list_genes_core.txt
```
This modified list can be found: ```results_cg/list_genes_core.txt```.

## Step 3: Scheme Validation (Allele calling)

For the validation step we selected 2759 unfinished *P. aeruginosa* genomes that were publicly available in RefSeq (https://www.ncbi.nlm.nih.gov/assembly) in September 2018. The list of all the draft genomes used can be found in the folder ```Genomes_Validation/Genomes_Validation.xlsx```. 

Multilocus sequence type (MLST) was determined as described in step 1. The sequence type (STs) obtained for each of the 2759 drafts genomes using the sanger-pathogens/mlst_check can be found in the folder ```Genomes_Validation/Genomes_Validation.xlsx```.

Of the 2759 drafts genomes available, it was possible to assign STs to 2686 drafts genomes and genomes that could not be assigned STs were not included for the validation of the cgMLST scheme.  A second filter was added to remove draft genomes that had ≥200 contigs and 502 genomes were removed. 

In the end, 73 drafts genomes were removed due to the absence of MLST *loci* and 502 were removed because the available sequences consisted of ≥200 contigs. Thus, of the 2759 drafts genomes obtained from RefSeq, 2184 genomes were used for the validation of the cgMLST scheme. The list of all the 2184 drafts genomes used to validation the schema obtained from RefSeq can be found in the folder: ```Genomes_Validation/Genomes_Validation.xlsx```.


We this set of genomes (2184 draft genomes) we repeated the allele call step using only the 3164 candidate *loci*.

## Command:

```bash
chewBBACA.py AlleleCall -i Genomes_Validation -g list_genes_core.txt -o results_all --cpu 15 --ptf PAO1.trn
```

The folder **Genomes_Validation** contains the 2184 validation draft genomes used for validation of the scheme.

The folder with the output file can be found at: ```results_all/ results_20191126T121343/```. This folder contains 5 files: "logging_info.txt; RepeatedLoci.txt; results_alleles.tsv; results_contigsInfo.tsv and results_statistics.tsv".

The ```results_all/ results_20191126T121343/results_alleles.tsv``` file contains the allelic profile of the 2184 typed drafts genomes.

Due to the size of the **results_contigsInfo.tsv** file, it was not possible to upload it to GitHub in the ```results_all/ results_20191126T121343/``` folder, but a link to access the file is available at: (https://drive.google.com/open?id=11_sZqOXK8bkWFFsvcZo1oq7-PqCGif_G).


## Step 3.1: Concatenate the allelic profiles

The purpose of this step is to concatenate the matrix of the *loci* that defined the scheme and the matrix of the *loci* from the validation genomes. Thus, to concatenate the allelic profile matrix obtained from the creation of the scheme ```cgMLST_120/cgMLST.tsv``` with the matrix obtained for the validation genomes ```results_all/ results_20191126T121343/results_alleles.tsv```.  The following command was used:

## Command:

```bash
# create header
head -n 1 cgMLST_120/cgMLST.tsv > cgMLST_all.tsv
```

## Command:

```bash
# concatenate
grep -v ^FILE cgMLST_120/cgMLST.tsv results_all/ results_20191126T121343/results_alleles.tsv >> cgMLST_all.tsv
```
The cgMLST_all.tsv file can be found in the folder: ```analysis_all/cgMLST_all.tsv```. This file (cgMLST_all.tsv) contains the allelic profile of the 2314 genomes.

## Step 3.2: Evaluation of genome quality

After concatenation, we used the *TestGenomeQuality* to assess the impact of each validation genome on candidate *loci*  in order to exclude low quality validation genomes. In this step you may need to define a new *Threshold*, as well as a new value of the parameter *p*, because *loci* that remain after the filters are the ones that will constituted the final scheme.

## Command:

```bash
 chewBBACA.py TestGenomeQuality -i cgMLST_all.tsv -n 13 -t 300 -s 5
```
The folder with the output file can be found at: ```analysis_all/GenomeQualityPlot.html```. The list with the *loci* present in 95% of the genomes at the chosen *Threshold* can be retrieved in the folder ```analysis_all/Genes_95%.txt```. The list of genomes removed at each *Threshold* can be retrieved in the folder ```analysis_all/removedGenomes.txt```.

In order to exclude validation genomes that have left the scheme it is necessary to follow the steps described in **Step 2.2**

## Step 4: Extracting the Matrix loci

At this step we chose *loci* present in 99% (*p0.99*) of the validation genomes and the *Threshold* 200 to limit the loss of the *loci* in the genomes. In *Threshold* 200 a set of 2653 gene targets were found to be present in 99% the validation genomes. 

From the original "removedGenomes.txt" file  that can be retrieved in the ```analysis_all/removedGenomes.txt``` folder we created another (removedGenomes200thr.txt) file with only the genomes removed at *Threshold* (200). The list of genomes removed at *Threshold 200* can be retrieved in the folder: ```analysis_all/removedGenomes200.txt```

Using *Threshold* (200) only 5 draft genomes were removed due to absence of *loci* targets.

To transpose (put the names of the validation genomes one in each line) we used the datamash and created the file > removedGenomes200thr.txt. This file can be found in the folder: ```analysis_all/removedGenomes200.txt```

## Command:

```bash
# transpose
datamash -W transpose < removedGenomes200.txt > removedGenomes200thr.txt
``` 
The genomes that were excluded in the *Threshold 200* have been placed in **removedGenomes200thr.txt**

This file can be found in the folder: ```analysis_all/removedGenomes200thr.txt```

## Command:

```bash
chewBBACA.py ExtractCgMLST -i cgMLST_all.tsv -o cgMLST_200 -p0.99 -g removedGenomes200thr.txt 
```

This script selects *loci* and genomes that remained in the *Threshold* 200 and excludes the validation genomes and *loci* that were excluded with this *Threshold*.

The folder with the output file can be found at: ```cgMLST_200```. This folder contains four files "cgMLST.tsv; cgMLSTschema.txt; mdata_stats.tsv and Presence_Absence.tsv".

The cgMLST targets can be found at: ```cgMLST_200/cgMLSTschema.txt``` It contains the list of 2653 genes in the core genome defined as targets for this scheme. The list with the 2653 cgMLST gene targets with the path added for each locus fasta file to be searched in the schema_seed folder is in the folder:```cgMLST_200/gene_targets.txt```

## Step 4.1: Checking the percentage of target genes found in the evaluated genomes

We assumed that for a well-defined cgMLST scheme 95% of the cgMLST gene targets should be present in all tested isolates. Thus, an isolate is considered to be "typed" with our cgMLST scheme if at least 95% of the 2653 targets were present.

In order to check if this was true,  we re-analyzed all 2325 genomes used in the creation and validation stages, even those that were in these different stages. 

## Command: 

```bash
chewBBACA.py AlleleCall -i genomes -g gene_targets.txt -o results_2195_genomes --cpu 15 --ptf PAO1.trn
```
**Note**: with the ```results_2195_genomes/``` folder, analyze the file ```results/results_statistics.tsv.```containing the number of genes that were found with 100% identity in the analyzed genomes (EXC - alleles which have exact matches (100% DNA identity) with previously identified alleles). **With the proposed schema we consider that an isolate is "typed" when at least 95% of the targets are found in the genome of interest**.

## Step 5: Minimum Spanning Tree

For the visualization of results, minimum spanning trees were buitl. Based on the allelic profiles obtained by the cgMLST scheme for each of the 2309 genomes minimum spanning trees (MST) were constructed using the software GrapeTree (version 1.5.0) (https://github.com/achtman-lab/GrapeTree/releases) with parameters implemented in MSTree v2 ignoring missing values for the entire strain collection. The ```cgMLST_200/cgMLST.tsv ``` file contains the allelic profile of the 2309 genomes typed by cgMLST.

## Step 6: Evaluation of the schema cgMLST

To assess the variability of the *loci* targets of cgMLST as well as the quality of the *loci* we run this script and graphically visualize the data in a series of html files.

## Command:

```bash
chewBBACA.py SchemaEvaluator -i schema_seed/ -l rms/RmS.html -ta 11 --title "cgMLST custom r sales" --cpu 6
```

Due to the size of the file rms/RmS.html it was not possible to upload it on GitHub, but a link to access the file rms/RmS.html is available at: (https://drive.google.com/open?id=1xzruUPM6yYpm0YVR_EPlcwnlfF_Zzk9H).

## Step 7: Analyze the proteins in the genes of the wgMLST

To check which protein encodes each *loci* found in the wg/cgMLST a list of proteins corresponding to all 13588 *loci* identified in the wg / cgMLST targets is found in the **new_protids.tsv** file. The list of proteins encoded by the 2653 cgMLST target genes can be found as: **Target_genes_cgMLST.xlsx**

## Command:

```bash
chewBBACA.py UniprotFinder -i schema_seed/ -t proteinID_Genome.tsv --cpu 10
```

## Step 8: How to identify the allelic profile of your *P. aeruginosa* genomes of interest using this cgMLST scheme?

In the steps below we describe all the necessary steps to identify the allelic profile of the genomes of interest using this cgMLST scheme.

## Step 8.1: A master directory to run ChewBBACA

A master directory ```Analyze_genomes/``` containing the folders that are needed to run ChewBBACA with the 2653 cgMLST target genes was created as an example. This folder contains a directory called ```example_genomes/``` representing the folder of the genomes to be typed. The second folder is the ```example_schema_seed/``` representing the ```schema_seed/``` folder that must be downloaded. In addition to two “gene_targets.txt” files containing the 2653 target genes of cgMLST and the file "PA01.trn" which is the prodigal's trained file to recognize CDs.

## Step 8.2: Install ChewBBACA and all its dependencies;

To download ChewBBACA and all its dependencies see the topic: Softwares and Downloads (Main dependencies); Other dependencies (for scheme evaluation only):

## Step 8.3: Download Scheme

To have access to the cgMLST scheme for *P. aeruginosa* it is necessary to download the ```schema_seed/``` folder. To access this ```schema_seed/``` folder just enter the link (https://drive.google.com/open?id=1WpsmbMC0awZ7BH8lazoiv-t56A2xnMIO).

**Note:** The ```schema_seed/``` folder was created in Step 1 where we identified all CDs of the 141 complete genomes creating the scheme.

## Step 8.4: Download the list of target genes

This list already contains the full path for each locus fasta file for ChewBBACA to fetch the target genes in the folder ```schema_seed/```. 

The list of cgMLST target genes can be obtained from the ```Analyze_genomes/gene_targets.txt```folder.

## Step 8.5: Download the Prodigal trained file

It is recommended to upload the trained folder (PA01.trn) because it is the reference folder for Prodigal to recognize the coding sequences (CDs). Prodigal was trained with the reference genome of *Pseudomonas aeruginosa* (PAO1) to recognize CDs.

## Step 8.6: Genomes of interest

The genomes of interest must be in a specific folder, for example the ```genomes/``` folder and must be in the same directory that contains the ```schema_seed/``` folder.

After this part, the next step is to run the command:

## Command: 

```bash
chewBBACA.py AlleleCall -i genomes -g gene_targets.txt -o results --cpu 15 --ptf PAO1.trn
```
**Note**:The folder “genomes" represents the folder with the genomes to be typed.

**Note**:The list "gene_targets.txt" containing the 2653 cgMLST target genes with the path to the schema_seed folder”

This command will release the output in the folder: ```results/```

## Step 9: Analyze the results

The allelic profile of the typed genomes will be in the folder: ```results/results_alleles.tsv``` which is the output of the file released by the script above. Other outputs will be in the folder ```results/```, as an example: RepeatedLoci.txt; logging_info.txt; results_contigsInfo.tsv and results_statistics.tsv.

**Note**: With the ```results/``` folder, analyze the file ```results/results_statistics.tsv.```It contains the number of genes that were found with 100% identity in the analyzed genomes (EXC - alleles which have exact matches (100% DNA identity) with previously identified alleles). **With the proposed scheme we consider that an isolate is considered typed when at least 95% of the targets are found in its genome**.

## Step 10: Allele profile view by minimum spanning tree (MST)

To view the allelic profile data of the typed genomes you need to access the output of the script present in the folder ```results/results_alleles.tsv```. It is possible to do so in two ways: the first is to download the free software GrapeTree (version1.5.0) with parameters implemented in MSTree v2 ignoring missing values for the entire strain collection available in (https://github.com/achtman-lab/GrapeTree/releases) or through the software Phyloviz, online.
