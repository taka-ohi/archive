
#Getting ready
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")

setwd("~/16S/dada2/")  ## CHANGE ME to the directory containing the fastq files.
filez <- list.files()
file.rename(from=filez, to=sub(pattern=".fastq", replacement=".fastq.gz", filez))

fnFs <- sort(list.files(getwd(), pattern = ".fastq", full.names = TRUE))

#Identify primers
FWD <-"GTGCCAGCMGCCGCGGTAA" ## 515F
REV <-"GGACTACVSGGGTATCTAAT" ## 806R


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


get.sample.name <- function(fname) strsplit(basename(fname), ".fastq")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name ))
head(sample.names)
fnFs.filtN <- file.path(getwd(), "filtN", paste0(sample.names, ".fastq.gz"))

#to “pre-filter” the sequences just to remove those with  ambiguous bases (Ns)
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE)

# Identifying and counting the primers on one set of paired end FASTQ files(Sufficient, because all the files were created using the same library preparation)
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]) )

# cutadapt
cutadapt <- "~/cutadapt" 
system2(cutadapt, args = "--version")


path.cut <- file.path(getwd(), "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnFs)) # dummy file

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, "-n", 2,
                             "-o", fnFs.cut[i],  # output file                             
                             fnFs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = ".fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)


#Filter and trim
filtFs <- file.path(path.cut, "filtered", basename(cutFs))

out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = 2,
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread = TRUE)
p <- plotErrors(errF, nominalQ = TRUE)

ggsave(plot=p, file="errors.png")

#Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

#Sample Inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)

#Construct Sequence Table
#see https://github.com/benjjneb/dada2/issues/384
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
table(nchar(getSequences(seqtab.nochim)))

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace

colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")

rownames(track) <- sample.names
head(track)
write.table(track,file="track.txt", quote=F)

# import database
DATABASE <- "~/database/silva/silva_nr_v132_train_set.fa.gz"

#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, DATABASE, multithread = TRUE, tryRC = TRUE)
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.table(taxa, file="taxonomy.txt", quote=F)
write.table(seqtab.nochim, file="seqtabnochim.txt", quote=F)

samples.out<-rownames(seqtab.nochim)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#To output OTU table
otu_table.t<-t(ps@otu_table)
ps.t<-cbind(otu_table.t,ps@tax_table)
write.table(ps.t,  file="ASV_table.txt", quote=F, sep = "\t")


