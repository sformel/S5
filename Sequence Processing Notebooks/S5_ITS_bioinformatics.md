# S5 Bioinformatics Notebook

Last updated Jan 6, 2020 by Steve Formel

###Description

This notebook is to keep track of the bioinformatics pipeline and analysis done on the ITS sequences.

###Downloading Files

Duke only allows sftp from their server.  Which is a pain because the Cypress linux system doesn't allow for recursive downloading with sftp.  So I'm trying to work around it (I don't want to use filezilla because it will take along time to download/upload everything to cypress and box.com. Plus multiple copies increases the chances of corruption, and I would have to check everything by hand.

1. Started idev node on Cypress:

		idev -t 4:00:00

2. Navigated to the folder where I wanted to put the sequences

3. Logged into the duke sftp using the info they sent me

		sftp bael_6074@dnaseq2.igsp.duke.edu

	*Note that the password will stay invisible when you type it.*

4. Navigate to the folder containing the sequence files.
5. Transfer all files

		get *

Took about 5-10 minutes.

##### Check file integrity after transfer

	md5sum *.fastq.gz > ./S5_SF1_MD5_Cypress.checksum

Sort the Duke checksum file and compare the checksum columns

	sort -k2 VanBael_6074_19121902.checksum > DS_sorted.checksum 

See if column one from both files matches

	awk 'NR==FNR{a[$1]=$1; next} a[$1] {print "MATCH"}' DS_sorted.checksum S5_SF1_MD5_Cypress.checksum | wc -l

Should print out the number of fastq.gz files there are if all match.  This should be twice the number of samples you ran, but if you're not convinced, you can count the number of files by:

	find *fastq.gz | wc -l

##### Copy to Box.com

For some reason my token for transferring files between Cypress and box.com isn't stable.  Every time I go to do it, I have to regenerate the token.  See Steve's github or rclone website for how to do this.

	module load rclone/1.49.3

Find correct folder to place files.

	rclone lsd VBlab_box:

Move SF1 files

	rclone copy ./SF1_ITS VBlab_box:VB_lab_Data/Sequences/SF1/SF1_ITS -P

Took 43s (471Mb)

Check for corruption

	rclone check ./SF1_ITS VBlab_box:VB_lab_Data/Sequences/SF1/SF1_ITS -P --one-way

Good. Now S5.

	rclone copy ./S5_ITS VBlab_box:VB_lab_Data/Sequences/S5/S5_ITS -P

Took about 9 min (4.883 Gb)

Check for corruption

	rclone check ./S5_ITS VBlab_box:VB_lab_Data/Sequences/S5/S5_ITS -P --one-way

Good.  Took about 24s.

##### Copy to Google Drive

Move SF1 files

	rclone copy ./SF1_ITS VBlab_google:VB_lab/'Boggs443 Data & User folders'/Data/Sequences/SF1/SF1_ITS -P

Took 30s (471Mb)

Check for corruption

	rclone check ./SF1_ITS VBlab_google:VB_lab/'Boggs443 Data & User folders'/Data/Sequences/SF1/SF1_ITS -P --one-way

Good. Now S5.

	rclone copy ./S5_ITS VBlab_google:VB_lab/'Boggs443 Data & User folders'/Data/Sequences/S5/S5_ITS -P

Took about 4.5 min (4.883 Gb)

Check for corruption

	rclone check ./S5_ITS VBlab_google:VB_lab/'Boggs443 Data & User folders'/Data/Sequences/S5/S5_ITS -P --one-way

Good.  Took about 20s.


##### Upload S5 seqs to NCBI

SF1 is not required until publishing.

I followed Candice's excellent protocol for submitting and used her S6 file and the submissions by Duke to try and be consistent in my descriptions.  Lastly, before submission I sorted the Biosamples by:

* organism
* env_local_scale
* env_medium
* collection_date
* sample_name

This was because Sunshine thought that Candice might have had some problems because her files weren't sorted before hand.


In any case, everything went fine, it was a tedious process, but not that hard.

#### Note on processing

It makes sense to process all the S5 and SF1 (different project that was sequenced with S5) files at the same time so that any errors that occurred in sequencing will be apparent and accounted for.  After calling sequence variants and using the decontam package I'll break them into separate tables for the two projects.  So I moved all the files into one folder on Cypress called ITS_seqs.

#### Fast QC

	#!/bin/bash
    
    #SBATCH --qos=normal
    #SBATCH --error job_reports/S5_ITS_FQC.error
    #SBATCH --output job_reports/S5_ITS_FQC.output
    #SBATCH --job-name S5_ITS_FQC
    #SBATCH --time=23:00:00
    #SBATCH --nodes=1
    #SBATCH --mem=64000
    #SBATCH --cpus-per-task=20
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=sformel@tulane.edu

	module load fastqc/0.11.7                    

	ls ITS_seqs/*.fastq.gz | xargs -P 2 -n 1 fastqc

Runtime:17 min

Moved all files (html and zip) into folder called fastQC

###### Concatenate all files to get a sense of overall quality:

cat *.fastq.gz > all_seqs.fastq.gz

run fastQC in idev mode.

	fastqc all_seqs.fastq.gz 

Copied fastQC folder to ITS_seqs folder in S5 data folder.


### Can't get DADA2 on Cypress

No matter what I try, it has a mess of problems.  I'm going to switch to the mac pro in the lab.

```

#Cleaning and Clustering S5 and SF1 ITS sequences
#Last updated by Steve Formel on Jan 10, 2020

#Description:  Following the standard DADA2 pipeline (v1.8) (https://benjjneb.github.io/dada2/ITS_workflow.html),  done on the mac pro in the VB lab because I couldn't get DADA2 to install properly on Cypress.

#Package Versions----
library(dada2)
library(ShortRead)
library(Biostrings)

#dada2 - 1.12.1
#tidyverse - 1.2.1
#ShortRead - 1.42.0
#Biostrings - 2.52.0

#Import files----

path <- "~/Desktop/ITS_seqs_BI/ITS_seqs/fq/" 
list.files(path)

#generate matched lists of forward and reverse files
fnFs <- sort(list.files(path, pattern = "1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "2_001.fastq.gz", full.names = TRUE))

#identify primers
FWD <- "CTTGGTCATTTAGAGGAAGTAA"  #ITS1f
REV <- "GCTGCGTTCTTCATCGATGC"  #ITS2r

#verify primer presence and orientation in data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

#pre-filter seqs with Ns (ambiguous bases)
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#Primers are found in appropriate places

#Forward Complement Reverse RevComp
#FWD.ForwardReads    7067          0       0       0
#FWD.ReverseReads       0          0       0    1333
#REV.ForwardReads       0          0       0    2256
#REV.ReverseReads    4938          0       0       0

#Remove Primers-----
cutadapt <- "/Users/vanbaellab/miniconda3/bin/cutadapt"
system2(cutadapt, args = "--version") # Run shell commands from R


#Create filenames and parameters
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt - adjusted by Steve to run in parallel on 8 cores
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i], # input files
                             "-j", 8,
                             "--minimum-length", 1)) 
}


#thoughtful discussion: https://forum.qiime2.org/t/trim-trunc-length-for-its/5048/36
#The above took about 20 min to run.  Some warnings which I feel are safe to ignore (see above link).

#Sanity check
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

#came up with one primer in REV.Reverse reads, I'll ignore for now.

#Begin DADA2 pipeline-----

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "2_001.fastq.gz", full.names = TRUE))


##THIS DOESN"T WORK< IT DELETES LEAF NAMES##############################################################

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_L")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))


###ABOVE DOESN"T WORK###########


#trim extra stuf added by sequencing company
sample.names <- gsub(pattern = "(.*_)(.*)_S.*", replacement = "\\1\\2", x = sample.names)
head(sample.names)

#Inspect Read quality profiles
plotQualityProfile(cutFs, aggregate = TRUE) #triggers scale warning which can be ignored, took 10 min
plotQualityProfile(cutRs, aggregate = TRUE)

#pdf saved by hand.  alternate way below:
#plot.quals <- plotQualityProfile(fn)
#ggsave("qualplot.pdf", plot.quals, device="pdf")

#The above plots (named R1_qual_plot.pdf and R2...) don't look amazing, but I'm going to proceed, trimming them at about 275.  This will keep everything above 20.  I may have to adjust this later.

#Filter and Trim

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

#took about 10 min
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

#learn error rates - took 2 min
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

#visualize error rates
plotErrors(errF, nominalQ = TRUE)  #error about transformation: no big deal https://github.com/benjjneb/dada2/issues/742

#pipeline to start before I go home-----

Start <- Sys.time()

#derep reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#sample inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#construct seq table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#inspect
table(nchar(getSequences(seqtab.nochim)))

#track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                       getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

#assign taxonomy

#downloaded UNITE general fasta per DADA2 webpage (https://plutof.ut.ee/#/doi/10.15156/BIO/786353)
#DOI: 10.15156/BIO/786353, v8.0, release Nov. 18, 2018
#DADA2 says: https://github.com/benjjneb/dada2/issues/702

unite.ref <- "/Users/vanbaellab/Desktop/ITS_seqs_BI/UNITE/sh_general_release_dynamic_s_02.02.2019.fasta"  # CHANGE ME to location on your machine
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

#inspect taxonomic assignments
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#save tables
saveRDS(seqtab.nochim, "/Users/vanbaellab/Desktop/ITS_seqs_BI/ITS_SVtab_done.rds")
saveRDS(taxa, "/Users/vanbaellab/Desktop/ITS_seqs_BI/ITS_TAXtab_done.rds")

end <- Sys.time()

#Took only 1 hour and 45 min

### Load data into New session
st <- readRDS("ITS_SVtab_done.rds")
taxtab <- readRDS("ITS_TAXtab_done.rds")

#Bring into Phyloseq-----
library(dada2); packageVersion("dada2") #version 1.12.1
library(phyloseq); packageVersion("phyloseq") #version 1.28.0
library(Biostrings); packageVersion("Biostrings") #version 2.52.0
library(ggplot2); packageVersion("ggplot2") #version 3.2.1

sampleID <- attributes(st)[[2]][1] #the names are messed up for leaf samples because I did  a bad job replacing strings above.  But I think I can replace them since they should be in the same order, and I'll rerun the DADA2 pipeline after I fix it.

library(readxl)
split_experimental_data <- as.data.frame(read_excel("ITS_seqs_exp_data.xlsx"))

rownames(split_experimental_data) <- split_experimental_data$sampleID

#fix dates
split_experimental_data$DNA_ext_date <- as.Date(split_experimental_data$DNA_ext_date)
  
#replace wonky sample names with correct from list
attributes(st)[[2]][1] <- list(split_experimental_data$sampleID)

#make into phyloseq object
ps <- phyloseq(otu_table(st, taxa_are_rows=FALSE), 
               sample_data(split_experimental_data), 
               tax_table(taxtab))

#Rename seqs to "otus" for display purposes
taxa_names(ps) <- paste0("Seq", seq(ntaxa(ps)))

#reorder tissue factor
split_experimental_data$tissue <- factor(split_experimental_data$tissue, levels = c("Leaf", "Root", "Soil","Inoc", "NegCon"))

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 6868 taxa and 238 samples ]
#sample_data() Sample Data:       [ 238 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 6868 taxa by 7 taxonomic ranks ]

#quick inspection of seq count by experiment and extraction date

df.reads <- data.frame(sample_data(ps))
df.reads$total_reads <- as.vector(sample_sums(ps))

ggplot(data = df.reads, aes(x= DNA_ext_date, y = total_reads)) +
  geom_point(aes(color = tissue),
             size = 2, alpha = 0.75) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Total Reads") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  facet_wrap(~ experiment, scales = "free")

ggsave(filename = "ITS_seqs_batch_effects.png", width = 14, height = 7, units = "in")

#COntrol for seqs in negative controls------

#https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

library(decontam); packageVersion("decontam") #version 1.4.0
library(phyloseq); packageVersion("phyloseq") #version 1.28.0
library(ggplot2); packageVersion("ggplot2") #version 3.2.1

#inspect library sizes
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

#ID contaminants via frequency method----
contamdf.freq <- isContaminant(ps, method="frequency", conc="purified_concentration")
head(contamdf.freq)

#count possible contaminants
table(contamdf.freq$contaminant)

#35 possible contaminants

#which ones are contaminants?  
which(contamdf.freq$contaminant)

#plot to help me decide if they are real - note that I'm plotting them all together because ther are only 35.  More than that would probably crash your computer.
ps.contam <- prune_taxa(contamdf.freq$contaminant, ps)

#how many sequences make up each contaminant?  Only a few are slightly abundant (> 100 seqs) and only one above 1000 seqs (#259)
sort(taxa_sums(ps.contam))

plot_frequency(ps, taxa_names(ps.contam), conc="purified_concentration") + 
  xlab("DNA Concentration (ng/ul)")

#The above evidence isn't totally convincing, but I don't know that it's a terrible idea to get rid of any of them.

#ID contaminants via prevalence method-----
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Negative Control"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")

table(contamdf.prev$contaminant) #15 contaminants

#more aggresive p-value
contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant) #32

#how many of the SVs (i.e. OTUs) were shared between the two methods?
intersect(taxa_names(ps.contam), rownames(contamdf.prev05[contamdf.prev05$contaminant,])) 

#only 5 are in both methods, but they include the 4 most abundant.  So at least I can feel confident getting rid of those 4.

#Let's see hwo they divvy up in negative controls and real samples

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Negative Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") +
  geom_text(label = rownames(df.pa))

#The good news is that 259 is prevalent in 5 out of 6 negative controls, so it's almost certainly real. However, there are a bunch that only show up in one sample (and recall this is with the more stringent p value).  I think the thing to do is get rid of the 5 contaminants that showed up in both methods.  And maybe do some quick comparisons of how results change if you don't remove contaminants.

#except there is a formal combined method!

contamdf.comb <- isContaminant(ps, method="combined", neg="is.neg", conc="purified_concentration")
table(contamdf.comb$contaminant) #10 contaminants, although not the 5 intersecting ones from before.  But it still gets #259.

#which ones?
rownames(contamdf.comb[contamdf.comb$contaminant,])

#how abundant are they?
ps.contam <- prune_taxa(contamdf.comb$contaminant, ps)
sort(taxa_sums(ps.contam))

#Wow, Seq29 has 49,638 seqs.  Let's plot them.
plot_frequency(ps, taxa_names(ps.contam), conc="purified_concentration") + 
  xlab("DNA Concentration (ng/ul)")

ps.pa <- transform_sample_counts(ps.contam, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Negative Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contamdf.comb[contamdf.comb$contaminant,])

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") +
  geom_text(label = rownames(df.pa))

#It does look like Seq29 might be a real contaminant since it's in 3 negative controls.  Although if it wasn't for it's presence there, it wouldn't be that convincing.

#there is a batch method (two experiments on the same sequencing run), but I don't have enough negative controls to break out SF1 into it's own control.  At least this way I'm controlling for contamination that happened during sequencing.

#Remove contaminants from phyloseq object-----
ps.clean <- prune_taxa(!contamdf.comb$contaminant, ps)

#Split into S5 and SF1 samples and remove negative controls----

S5 <- prune_samples(sample_data(ps.clean)$experiment=="S5" & 
                      sample_data(ps.clean)$tissue!="NegCon", 
                    ps.clean)

SF1 <- prune_samples(sample_data(ps.clean)$experiment=="SF1" & 
                       sample_data(ps.clean)$tissue!="NegCon", 
                     ps.clean)

#save objects

saveRDS(object = S5, file = "S5_ITS_ps.rds")
saveRDS(object = SF1, file = "SF1_ITS_ps.rds")

```

