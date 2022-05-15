
Last updated May 14, 2022 by Steve Formel

Description:  This is the notebook for the processing of the sequences into a feature table for ecological analysis. I needed to download the raw sequences from NCBI before I run them through the standard DADA2 pipeline.

### Download from NCBI

Bio project ID = PRJNA597162

	wget 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&term=PRJNA597162' -O - | tee S5_16S_SraRunInfo.csv
	
Cut download addresses into new txt file

	cut -d "," -f 10 S5_16S_SraRunInfo.csv > S6_16S_SRR_download_paths.txt
	
Download SRA toolkit from ncbi

Note that we need an older version to function on Cypress (Tulane HPC)

	wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.0/sratoolkit.2.8.0-2-centos_linux64.tar.gz --no-check-certificate
	
	tar -xzvf sratoolkit.2.8.0-2-centos_linux64.tar.gz

Download all the SRR files

	wget -i S5_16S_SRR_download_paths.txt 
	
Convert to fastq

For help: https://ncbi.github.io/sra-tools/fastq-dump.html

	./sratoolkit.2.8.0-2-centos_linux64/bin/fastq-dump --outdir ./fastq --split-files ./SRR*

Wow, as far as I can tell this is useless.  When you convert the SRR, it gives each sequence a new ID that has no info related to the original sample name.  What a terrible system.  So I'm going to make add to the data spreadsheet each 64 SRA and the amazon download links for each one.  I'm doing this by clicking Bioproject SRA in ncbi and then manually going to the SRA data page and getting the links for each original file.

But I did realize that NCBI adds the SRR read hash in their metadata (what I downloaded up top) so that was a nice touch I didn't anticipate.

This makes the argument pretty clear as to why each sample should have it's own biosample number.

To get the original file links:

1. Made a file with a list of links to the SRR webpage

		like: https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10755313
	
	You can get the SRR strings from the file you made above.

2. Download all the source code

		wget -i SRR_links.txt 
		
3. Grep the code for the amazon links

		module load gnuparallel/20180322

		find . -type f | parallel -X grep -o '\"https://sra-pub-src-2.s3.amazonaws.com/.*.R1_001.fastq.gz.1\"' > R1_amazon_links.txt
		
4. finished in excel, put the list back into a text file to download them all with wget

		cat S5_16S_amazon_links.txt | parallel wget
		
Took about 5 minutes.

5. copy to box.com and google drive with rclone


		rclone copy ./ VBlab_box:VB_lab_Data/Sequences/S5/S5_16S
 
Made checksum file.  Not really as good as getting the actual file from the Gunsch lab, but they're not handing out data right now...

	md5sum *.fastqq.gz.1 > ../SF1_MD5_Cypress_check.txt

Note that the one on the end of the file can be removed so these will be read as fastq.gz files.
		
### DADA2 pipeline

Started idev node on Cypress:

		idev -t 4:00:00

Remove .1 from end of every fastq file

	for i in *.fastq.gz.1
		do
  	mv -- "$i" "${i/%.1/}"
	done		
#### Fast QC

	module load fastqc/0.11.7                    

	ls ./*.fastq.gz | xargs -P 2 -n 1 fastqc

Runtime: about 10 min

Moved all files (html and zip) into folder called fastQC

###### Concatenate all files to get a sense of overall quality:

	cat *R1_001.fastq.gz > all_R1.fastq.gz
	cat *R2_001.fastq.gz > all_R2.fastq.gz
run fastQC

	ls ./all_R*.fastq.gz | xargs -P 2 -n 1 fastqc 

Copied fastQC folder to 16S_seqs folder in S5 data folder.

###Trimming decision

These are really nice sequences.  They used the 819F and 1115Rmod primers according to what they put on ncbi.  Since this was 250 PE sequencing, that give us plenty to work with.  Based off the fastqc I'm going to trim to 230 and see how that looks.

#### Begin dada2 pipeline

https://benjjneb.github.io/dada2/tutorial.html

	idev -t 4:00:00
	
	module load R/3.6.1-intel gcc/6.3.0
	module swap intel-psxe/2015-update1 intel-psxe/2019-update1
	
	export LD_LIBRARY_PATH=/share/apps/intel_parallel_studio_xe/2019_update1/intelpython3/lib:$LD_LIBRARY_PATH
	
	R
	
	library(dada2)
	
tutorial 1.12

	path <- "/lustre/project/svanbael/steve/S5/16S/fastq"
	list.files(path)

get matched lists of the forward and reverse fastq files.

	# Forward and reverse fastq filenames have format: SAMPLENAME_L001_R1_001.fastq and SAMPLENAME_L001_R2_001.fastq
	fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq", full.names = TRUE))
	fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq", full.names = TRUE))
	
	# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
	sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

Made quality plots in fastqc (see above)


Assign filenames for filtered files

	# Place filtered files in filtered/ subdirectory

	filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
	filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

	names(filtFs) <- sample.names
	names(filtRs) <- sample.names

Filter by standard parameters

	out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,230),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

	# took about 10 min

Learn error rates

	errF <- learnErrors(filtFs, multithread=TRUE)
	errR <- learnErrors(filtRs, multithread=TRUE)

#Come back to this

Sanity check: do error rates look ok?

	plotErrors(errF, nominalQ=TRUE)

Sample Inference

	dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
	dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
	

Inspect object

	dadaFs[[1]]
	dadaRs[[1]]

Merge paired reads	
	mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
	
	# Inspect the merger data.frame from the first sample
	head(mergers[[1]])

Construct Seq Table

	seqtab <- makeSequenceTable(mergers)
	dim(seqtab)

	# Inspect distribution of sequence lengths
	table(nchar(getSequences(seqtab)))

Remove Chimeras

	seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
	dim(seqtab.nochim)

Track Reads through the pipeline

	getN <- function(x) sum(getUniques(x))
	track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
	
	# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

	colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
	rownames(track) <- sample.names
	head(track)

Download Silva 132 taxonomy

	wget https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1

Assign taxonomy

	taxa <- assignTaxonomy(seqtab.nochim, "/lustre/project/svanbael/steve/S5/16S/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

Inspect taxonomy

	taxa.print <- taxa # Removing sequence rownames for display only
	rownames(taxa.print) <- NULL
	head(taxa.print)

#####Whole thing as script

	#!/bin/bash
    
    #SBATCH --qos=normal
    #SBATCH --error S5_16S_dada2.error
    #SBATCH --output S5_16S_dada2.output
    #SBATCH --job-name S5_16S_dada2
    #SBATCH --time=23:00:00
    #SBATCH --nodes=1
    #SBATCH --mem=64000
    #SBATCH --cpus-per-task=20
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=sformel@tulane.edu


	module load R/3.6.1-intel gcc/6.3.0
	module swap intel-psxe/2015-update1 intel-psxe/2019-update1
	
	export LD_LIBRARY_PATH=/share/apps/intel_parallel_studio_xe/2019_update1/intelpython3/lib:$LD_LIBRARY_PATH
	
	cd /lustre/project/svanbael/steve/S5/16S/

	Rscript S5_16S_dada2.R
	
	library(dada2)
	
	#tutorial 1.12

	path <- "/lustre/project/svanbael/steve/S5/16S/fastq"
	list.files(path)

	#get matched lists of the forward and reverse fastq files.

	# Forward and reverse fastq filenames have format: SAMPLENAME_L001_R1_001.fastq and SAMPLENAME_L001_R2_001.fastq
	fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq", full.names = TRUE))
	fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq", full.names = TRUE))
	
	# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
	sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

	#Made quality plots in fastqc (see above)


	#Assign filenames for filtered files

	# Place filtered files in filtered/ subdirectory

	filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
	filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

	names(filtFs) <- sample.names
	names(filtRs) <- sample.names

	#Filter by standard parameters

	out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,230),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

	# took about 10 min

	#Learn error rates

	errF <- learnErrors(filtFs, multithread=TRUE)
	errR <- learnErrors(filtRs, multithread=TRUE)

	#Sanity check: do error rates look ok?

	#plotErrors(errF, nominalQ=TRUE) #Come back to this

	#Sample Inference

	dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
	dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

	#Inspect object

	dadaFs[[1]]
	dadaRs[[1]]

	#Merge paired reads	
	mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
	
	# Inspect the merger data.frame from the first sample
	head(mergers[[1]])

	#Construct Seq Table

	seqtab <- makeSequenceTable(mergers)
	dim(seqtab)

	# Inspect distribution of sequence lengths
	table(nchar(getSequences(seqtab)))

	#Remove Chimeras

	seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
	dim(seqtab.nochim)

	#Track Reads through the pipeline

	getN <- function(x) sum(getUniques(x))
	track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
	
	# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

	colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
	rownames(track) <- sample.names
	head(track)

	#Download Silva 132 taxonomy

	#wget https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1

	#Assign taxonomy

	taxa <- assignTaxonomy(seqtab.nochim, "/lustre/project/svanbael/steve/S5/16S/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

	#Inspect taxonomy

	taxa.print <- taxa # Removing sequence rownames for display only
	rownames(taxa.print) <- NULL
	head(taxa.print)

	#make phyloseq object

	library(phyloseq); packageVersion("phyloseq")
	library(Biostrings); packageVersion("Biostrings")

	#has no sample data
	ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))

	dna <- Biostrings::DNAStringSet(taxa_names(ps))
	names(dna) <- taxa_names(ps)
	ps <- merge_phyloseq(ps, dna)
	taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
	ps

	saveRDS(ps, "S5_16S_dada2_output_no_sample_data.rds")

#### Below saved as: S5\_16S_dada2.R

	library(dada2)
	
	#tutorial 1.12
	
	path <- "/lustre/project/svanbael/steve/S5/16S/fastq"
	list.files(path)
	
	#get matched lists of the forward and reverse fastq files.
	
	# Forward and reverse fastq filenames have format: SAMPLENAME_L001_R1_001.fastq and SAMPLENAME_L001_R2_001.fastq
	fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq", full.names = TRUE))
	fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq", full.names = TRUE))
	
	# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
	sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
	
	#Made quality plots in fastqc (see above)
	
	
	#Assign filenames for filtered files
	
	# Place filtered files in filtered/ subdirectory
	
	filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
	filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
	
	names(filtFs) <- sample.names
	names(filtRs) <- sample.names
	
	#Filter by standard parameters
	
	out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,230), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
	
	# took about 10 min
	
	#Learn error rates
	
	errF <- learnErrors(filtFs, multithread=TRUE)
	errR <- learnErrors(filtRs, multithread=TRUE)
	
	#save environment
	save.image("S5_16S_image.RData")
	
	#Sanity check: do error rates look ok?
	
	#plotErrors(errF, nominalQ=TRUE) #Come back to this
	
	#Sample Inference
	
	dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
	dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
	
	
	#save environment
	save.image("S5_16S_image.RData")
	
	
	#Inspect object
	
	dadaFs[[1]]
	dadaRs[[1]]
	
	#Merge paired reads
	mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
	
	# Inspect the merger data.frame from the first sample
	head(mergers[[1]])
	
	#Construct Seq Table
	
	seqtab <- makeSequenceTable(mergers)
	dim(seqtab)
	
	# Inspect distribution of sequence lengths
	table(nchar(getSequences(seqtab)))
	
	#Remove Chimeras
	
	seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
	dim(seqtab.nochim)
	
	#Track Reads through the pipeline
	
	getN <- function(x) sum(getUniques(x))
	track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
	
	# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
	
	colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
	rownames(track) <- sample.names
	head(track)
	
	#save environment
	save.image("S5_16S_image.RData")
	
	
	#Download Silva 132 taxonomy
	
	#wget https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1
	
	#Assign taxonomy
	
	taxa <- assignTaxonomy(seqtab.nochim, "/lustre/project/svanbael/steve/S5/16S/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
	
	#Inspect taxonomy
	
	taxa.print <- taxa # Removing sequence rownames for display only
	rownames(taxa.print) <- NULL
	head(taxa.print)
	
	#save environment
	save.image("S5_16S_image.RData")
	
	
	#make phyloseq object
	
	library(phyloseq); packageVersion("phyloseq")
	library(Biostrings); packageVersion("Biostrings")
	
	#has no sample data
	ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa))
	
	dna <- Biostrings::DNAStringSet(taxa_names(ps))
	names(dna) <- taxa_names(ps)
	ps <- merge_phyloseq(ps, dna)
	taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
	ps
	
	#save environment
	save.image("S5_16S_image.RData")

#### Need to change:  important that save.image path be explicit

The above worked up until the assign taxonomy step and then threw some errors.  So I loaded the RData file and ran as job.  It turns out that multithreading is buggy.  (you get a segmentation fault error).  So run without multithreading.

S5_16S_assign_tax.sh

	#!/bin/bash
    
    #SBATCH --qos=normal
    #SBATCH --error S5_16S_dada2.error
    #SBATCH --output S5_16S_dada2.output
    #SBATCH --job-name S5_16S_dada2
    #SBATCH --time=23:00:00
    #SBATCH --nodes=1
    #SBATCH --mem=64000
    #SBATCH --cpus-per-task=20
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=sformel@tulane.edu


	module load R/3.6.1-intel gcc/6.3.0
	module swap intel-psxe/2015-update1 intel-psxe/2019-update1
	
	export LD_LIBRARY_PATH=/share/apps/intel_parallel_studio_xe/2019_update1/intelpython3/lib:$LD_LIBRARY_PATH
	
	cd /lustre/project/svanbael/steve/S5/16S/

	Rscript S5_16S_dada2.R

Rscript S5_16S_assign_tax.R

	#Assign taxonomy
	library(dada2)

	load("S5_16S_image.RData")

	taxa <- assignTaxonomy(seqtab.nochim, "/lustre/project/svanbael/steve/S5/16S/silva_nr_v132_train_set.fa.gz")

	save.image("S5_16S_image.RData")


The above still didn't work on Cypress.  So I moved it to the macpro and finished taxonomy assignment. Saved on VB Lab drive as: "S5\_16S\_assigntax\_macpro.R"

```

#Assign taxonomy for S5 16S

library(dada2)
library(phyloseq)
library(Biostrings)

#Last updated Feb 1, 2020

load("/Users/vanbaellab/Google Drive/VBL_data/Collections, Projects/Collections/S5 collection/data/16S/dada2_env/S5_16S_image.RData")

#Assign taxonomy

taxa <- assignTaxonomy(seqtab.nochim, "/Users/vanbaellab/Google Drive/VBL_data/Collections, Projects/Collections/S5 collection/data/16S/dada2_env/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

#Inspect taxonomy

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#save environment
save.image("/Users/vanbaellab/Google Drive/VBL_data/Collections, Projects/Collections/S5 collection/data/16S/dada2_env/S5_16S_image_macpro.RData")


#make phyloseq object

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")

#has no sample data
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#save environment
save.image("/lustre/project/svanbael/steve/S5/16S/S5_16S_image_macpro.RData")
```

Make phyloseq object saved for future work.  Script saved on VB Lab google drive as "make_16S_phyloseq.R"

```
#Script for S5
#Purpose: Import 16S R data made by DADA2 in Cypress and the macpro and replace/clean sample data

#Last updated by Steve Formel
#Feb 3, 2020

#subtitle
subtitle = "Exploratory Analysis, 16S_diversity.R"

#need to get extraction date and sequencing info to test for batch effects

library(dada2)

#import env made on Cypress
#load("data/16S/dada2_env/S5_16S_image.RData")  #original without taxa
load("data/16S/dada2_env/S5_16S_image_macpro.RData")  #image from Cypress with taxa assigned on MacPro

#make phyloseq object

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(readxl)

#has no sample data
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#read in sample data and clean-----

samdf <- as.data.frame(read_excel("data/S5_final_data_4Mar2020.xlsx", 
                                  sheet = "R_16S", na = "NA"))
#reorder to match phyloseq object
samdf <- samdf[match(rownames(seqtab.nochim), samdf$duke_sampleID),]

#match?
any(rownames(seqtab.nochim)!=samdf$duke_sampleID)

#match rownames
rownames(samdf) <- rownames(seqtab.nochim)

#double check that rownames matches sampleID
any(rownames(samdf)!=samdf$duke_sampleID)


#cleaning and fixing----

#reorder factors
samdf$sampling_period <- factor(samdf$sampling_period, levels = c("NA", 
                                                                  "050216",
                                                                  "N16",
                                                                  "20317",
                                                                  "40117",
                                                                  "51817",
                                                                  "J17",
                                                                  "91817",
                                                                  "N17",
                                                                  "121317",
                                                                  "30318",
                                                                  "51418",
                                                                  "J18"))

samdf[,38:46] <- lapply(samdf[,38:46],as.numeric)

#rename and reorder levels
library(plyr)

samdf$special_char <- factor(samdf$special_char, levels = c("Pre-Experiment",
                                                            "Exp. Sample",
                                                            "Prev-Oiled Inoc.", 
                                                            "Not Prev-Oiled Inoc"))

samdf$oil_added <- revalue(samdf$oil_added, c("Y" = "Oil Added", "N" = "No Oil Added"))

samdf$orig_soil <- revalue(samdf$orig_soil, c("Y" = "Prev-Oiled Inoc.", "N" = "Not Prev-Oiled Inoc"))

samdf$Treatment <- interaction(samdf$oil_added, samdf$orig_soil,sep = " ; ")

#rename categories
library(tidyverse)

colnames(samdf) <- samdf %>% 
  dplyr::rename("Hopanes" = hopanes, "Total Chrys." = total_methy_chry, "Total Naph." = total_methy_naph,"Total Dibenz." = total_methy_dibenz, "Total Phenan." = total_methy_phen) %>% 
  colnames()

samdf$pH <- round(as.numeric(samdf$pH),digits = 2)
samdf$conductivity <- round(as.numeric(samdf$conductivity),digits = 2)

#make select columns factors and integers

cols <- c("sampleID", "tissue", "DNA_ext_date", "plantID", "plant_trt", "orig_soil", "oil_added", "table", "block", "season", "Treatment")
samdf[cols] <- lapply(samdf[cols], factor)

cols <- c("num_nodes", "num_stems_live", "num_stems_outside", "num_infl")
samdf[cols] <- lapply(samdf[cols], as.integer)


#make phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa), sample_data(samdf))

#rename
S5 <- ps

#Save phyloseq object as R object for future use.

saveRDS(S5, "S5_16S_pseq.rds")

```