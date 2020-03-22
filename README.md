# cgMLST-Pseudomonas-aeruginosa

## chewBBACA

To download ChewBBACA access the link:

* https://github.com/B-UMMI/chewBBACA_tutorial

## Workflow

* Step 1: Schema Creation
* Step 2: Allele calling
* Step 3: Schema Validation (Allele call)
* Step 4: Extracting the Matrix *loci*
* Step 5: Minimum Spanning Tree (MST)
* Step 6: Evaluation of the schema cgMLST
* Step 7: Analyze the proteins in the genes of the wgMLST
* Step 8: How to identify the allelic profile of the genomes of interest using this cgMLST database for *P. aeruginosa*
* Step 9: Analyze the results
* Step 10: Allele profile view by minimum spanning tree (mst)
 
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

For *Pseudomonas aeruginosa* select the option RefSeq from GenBank at https://www.ncbi.nlm.nih.gov/assembly. RefSeq corresponds to a comprehensive, non-redundant, well-annotated set of reference sequences. A set of 142 complete genomes sequences of *P. aeruginosa* were publicly available in RefSeq (https://www.ncbi.nlm.nih.gov/assembly) in September 2018. The list of all the complete genomes used to create the schema obtained from RefSeq can be found in the folder "Complete_Genomes" in the Complete_Genomes.xlsx format.

Multilocus sequence type (MLST) for the 142 complete genomes was determined using (https://github.com/sanger-pathogens/mlst_check) and the MLST schema for *P. aeruginosa* (www.pubmlst.org; downloaded September 2019). New sequence types (STs) were assigned a unique internal identifier (STs ≥4000). The sequence type (STs) obtained for each of the 142 complete genomes using the sanger-pathogens/mlst_check can be found in the folder "Complete_Genomes/Complete_Genomes.xlsx". 

Among the 142 genomes, *Pseudomonas aeruginosa* PAO1 reference genome (GCF_000006765.1) was included so that the Prodigal algorithm could use it as reference to recognize coding sequences (CDs). Prodigal generated the PAO1.trn file at this step. 

**The PAO1 genome was then removed from further analysis**.

## Step 1: Definition of CDs sequences 
