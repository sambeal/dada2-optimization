###########################
#dada2-optimization
###########################
#this directory contains:
#input folder of gDNA from 16 AWF (B2 fish was removed due to failure of 1 PCR rep) + SINE taxonomy fasta
#individual output folders for each set of parameter trials

#output-1 = all default
#output-2 = cutadapt -e 0.2,  --discard-untrimmed,  --minimum-length 1
#output-3 = cutadapt -e 0.1,  --discard-untrimmed,  --minimum-length 1

###########################

###########################

#load libaries
library(dada2)
library(ShortRead)
library(Biostrings)
library(R.utils)

#make sure in correct directory
#want to be in "/Users/samanthabeal/Documents/MSc/Bioinformatics/dada2-optimization"
#setwd("..")
getwd()

#copy input files to new output directory
#directory should contain all sequence files (unzipped, primers removed)
R.utils::copyDirectory("input", "output3")

path <- "output3"
list.files(path)

#housekeeping
#generate matched lists of the forward and reverse read files, as well as parsing out the sample name
# Forward fastq file names have format: SAMPLENAME_L001_R1_001.fastq
fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq", full.names = TRUE))

#identify primers
FWD <- "TAGCTCAGCTGGTAGAGCAC"
REV <- "TGCCATTTAGCAGACGCTTTT"

#check
#verify the presence and orientation of these primers in the data
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
REV.orients

#pull out ambiguous bases -- not sure I need to but it's in the tutoral
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE)

#count the number of times the primers appear in the forward and reverse read, 
#while considering all possible primer orientations.
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]))

#                 Forward Complement Reverse RevComp
#FWD.ForwardReads   24479          0       0       0
#REV.ForwardReads       0          0       0   22501

#cutadapt -- remove primers
#run 'whereis cutadapt' in terminal to get the path
cutadapt <- "/Library/Frameworks/Python.framework/Versions/3.10/bin/cutadapt"
#system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-e", 0.1, # -e 0.2 = maximum error rate of 20% allowed (0.1 = default)
                             "-o", fnFs.cut[i], # output files
                             fnFs.filtN[i],# input files
                             "--discard-untrimmed",
                             "--minimum-length", 1)) 
}


#check where output files went
fnFs.cut 
#"output3/cutadapt/190222-A1-1-42_S489_L001_R1_001.fastq" 
#move to dada2-optimization.sh script to count seqs BEFORE moving on to other steps


#check
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]))

#i have 1 forward primer left somewhere
#from https://github.com/benjjneb/dada2/issues/675: a small number of primers left is ok to keep going w the subsequent analysis

#samples are now ready to be analyzed with dada2
# Forward and reverse fastq filenames have the format SAMPLENAME_L001_R1_001.fastq:
cutFs <- sort(list.files(path.cut, pattern="_L001_R1_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
#this code makes the sample names "190222-A1-1-42" (ie. removes everything after the first _)

#plot read quality profiles
plotQualityProfile(cutFs[1:2])
#quality drops off around 100bp (expected amplicon length = 90)
#quality plots weren't working so enforced minlength of 1 (0s won't plot)

#will not be trimming to a fixed length

#Assigning the filenames for the output of the filtered reads to be stored as fastq.gz files.
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
#first tutorial had unzipped files, this tutorial is zipped.
#my files are unzipped so will need to run again to see if un/zipped makes a difference 

#Filter and trim
out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = 2, 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out)
View(out)
#file names are not shortened here but they are in the tutorial @_@

#learn error rate
errF <- learnErrors(filtFs, multithread = TRUE)
#output-1 = 101568798 total bases in 885772 reads from 23 samples will be used for learning the error rates.
#output-2 = 102869524 total bases in 1057895 reads from 37 samples will be used for learning the error rates.

#plot error rates
plotErrors(errF, nominalQ=TRUE)
#A2A, C2C,& T2T points are not smooth

#Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

#Sample Inference - denoising
#At this step, the core sample inference algorithm is applied to the dereplicated data.
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaFs[[1]]
dadaFs

#output-1 = 84 sequence variants were inferred from 18131 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
#output-2 = 82 sequence variants were inferred from 18455 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16


#Construct sequence table
#from https://github.com/benjjneb/dada2/issues/795: 
#"if using single-end reads, make sequence table out of the dadaFs object, i.e. makeSequenceTable(dadaFs)"
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
#output-1 = 48 1652 (48 samples, 1652 ASVs)
#output-2 = 48 821 (48 samples, 821 ASVs)

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#output-1 = Identified 58 bimeras out of 1652 input sequences.
#output-2 = Identified 58 bimeras out of 821 input sequences.

#Sequence length distribution
table(nchar(getSequences(seqtab.nochim)))
#output-1 = range from 51-151bp, majority are 151 (upstream won't keep anything shorter than 50bp (default parameters))
#output-2 = range from 51-148bp, majority are 148


#SEQTAB.NOCHIM IS THE DATA I WANT! #READS/SAMPLE/ASV
#ASVs are not named but instead are actual sequences

sum(seqtab.nochim)/sum(seqtab)
#output-1 = 0.997258 (~1% of sequence reads = chimeras)
#output-2 = 0.9962435 (~1% of sequence reads = chimeras)

#track reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))

colnames(track) <- c("cutadapt", "filtered", "denoisedF", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

#think this tracks the reads after primer removal, not input
#without specifying -discard-untrimmed, cutadapt will keep everything
#output-3 onward = have changed the column to primer removal

#save tracking info as .csv for later analysis
write.csv(track, "output3/output3_track.csv",
          row.names = FALSE)

#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "output3/SINE_taxonomy.fasta", multithread=TRUE)

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

#check
track
taxa.print
#the track df tells me how many reads/sample, but not how many ASVs those reads belong to
#the taxa.print df tells me how many ASVs I have, but not the depth or which samples they came from

#got to here before lab meeting


##### assess seqtab.nochim data-----
#make a backup
seqtab.nochim1 <- seqtab.nochim
seqtab.nochim1 <- as.data.frame(seqtab.nochim1)

#transpose so can see the read depth of each ASV/sample
t_seqtab.nochim <- as.data.frame(t(seqtab.nochim1))

#write out as .csv
write.csv(t_seqtab.nochim, "output3/output3_ASVs.csv",
          row.names = TRUE)


