title: "IC-BPS-Longitudinal dada2"
author: "Seth Reasoner"
date: "2024-06-21"

##We would like to credit the following tutorial for guiding this code:
##"https://ycl6.github.io/16S-Demo/index.html"

##Raw sequencing data is uploaded in NCBI BioProject PRJNA1141486

##Create Folders
IC_fastq = "IC_fastq"        # raw fastq files
IC_filtered = "IC_filtered"  # dada2 trimmed fastq files
IC_output = "IC_output"      # output files
IC_images = "IC_images"       # output images


if(!dir.exists(IC_filtered)) dir.create(IC_filtered)
if(!dir.exists(IC_output)) dir.create(IC_output)
if(!dir.exists(IC_images)) dir.create(IC_images)

##Create List of File Names
fns_IC = sort(list.files(IC_fastq, full.names = TRUE))
fnFwds_IC = fns_IC[grep("R1.fastq", fns_IC)] # forward reads
fnRevs_IC = fns_IC[grep("R2.fastq", fns_IC)] # reverse reads

##Create list of sample names
IC_sample_names = gsub("_R1.fastq", "", basename(fnFwds_IC))


##Visualize Read Quality
ii = 1:length(IC_sample_names)
pdf(paste0(IC_images, "/plotQualityProfile.pdf"), width = 8, height = 8, pointsize = 12)
for(i in ii) {
  message(paste0("[", i ,"/", length(IC_sample_names), "] ", IC_sample_names[i]))
  print(plotQualityProfile(fnFwds_IC[i]) + ggtitle("Fwd"))
  print(plotQualityProfile(fnRevs_IC[i]) + ggtitle("Rev"))
}
invisible(dev.off())


##Filter and Trim
filtFwds_IC = file.path(IC_filtered, basename(fnFwds_IC))
filtRevs_IC = file.path(IC_filtered, basename(fnRevs_IC))

out_IC = filterAndTrim(fnFwds_IC, filtFwds_IC, fnRevs_IC, filtRevs_IC,
                        truncLen = c(240,240), minLen = 200, maxN = 0, truncQ = 2, maxEE =c(2,2),
                        rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = FALSE)

out_IC = as.data.frame(out_IC)
rownames(out_IC) = IC_sample_names
head(out_IC, 10)

##Extraction blank "IC-ExtBlank_G6_P1" had no reads that met the filtering criteria, reducing the sample count 

##create list of file names of filtered files 
filt_fns_IC = sort(list.files(IC_filtered, full.names = TRUE))
filtFwds_IC_final = filt_fns_IC[grep("R1.fastq", filt_fns_IC)]
filtRevs_IC_final = filt_fns_IC[grep("R2.fastq", filt_fns_IC)]


##Learn the Error Rates
errFwd_IC = learnErrors(filtFwds_IC, multithread = FALSE)
errRev_IC = learnErrors(filtRevs_IC, multithread = FALSE)


##Sample Inference
dadaFwds_IC = dada(filtFwds_IC_final, err = errFwd_IC, pool = FALSE, multithread = FALSE)
dadaRevs_IC = dada(filtRevs_IC_final, err = errRev_IC, pool = FALSE, multithread = FALSE)

##Merge Paired Reads
merged_IC = mergePairs(dadaFwds_IC, filtFwds_IC_final, dadaRevs_IC, filtRevs_IC_final, verbose = TRUE)

##Construct the Sequence Table
seqtab_IC <- makeSequenceTable(merged_IC)

table(nchar(getSequences(seqtab_IC)))

##Remove Chimeras
seqtab_nochim_IC <- removeBimeraDenovo(seqtab_IC, method="consensus", multithread=2, verbose=TRUE)
##Percent of sequences that are non-chimeric 
sum(seqtab_nochim_IC)/sum(seqtab_IC)

rownames(seqtab_nochim_IC)

IC_sample_names = gsub("_R1.fastq", "", rownames(seqtab_nochim_IC))

rownames(seqtab_nochim_IC) <- gsub("_R1.fastq", "", rownames(seqtab_nochim_IC))


##Track Reads through the Pipeline
getN <- function(x) sum(getUniques(x))
track_IC <- cbind(out_IC, sapply(dadaFwds_IC, getN), sapply(dadaRevs_IC, getN), sapply(merged_IC, getN), rowSums(seqtab_nochim_IC))
colnames(track_IC) <- c("input", "filtered", "denoisedF", "denoisedR", "merged","nonchim")
rownames(track_IC) <- IC_sample_names
head(track_IC)


##Assign Taxonomy
seqs_IC = getSequences(seqtab_nochim_IC)
dbpath = "16S_DB/"
ref_db = paste0(dbpath, "silva_nr99_v138.1_train_set.fa.gz")

taxonomy_tab_IC <- assignTaxonomy(seqs_IC, refFasta = ref_db)

taxa <- addSpecies(seqs_IC, "silva_v138.2_assignSpecies.fa.gz")


IC_Metadata <- read.csv("https://raw.github.com/reaset41/IC-Urobiome-Longitudinal/refs/heads/main/IC_Metadata.csv")

##Create Phyloseq Object
IC_phyloseq_unfiltered <- phyloseq(otu_table(seqtab_nochim_IC, taxa_are_rows=FALSE), 
                                    tax_table(taxonomy_tab_IC), sample_data(IC_Metadata))


all_tax <- tax_table(IC_phyloseq_unfiltered)
all_ASV <- otu_table(IC_phyloseq_unfiltered)
Table <- cbind(all_tax,t(all_ASV))

write.csv(Table, "Table_PreDecontam.csv", row.names=TRUE)


## Save Phyloseq Object
save(IC_phyloseq_unfiltered, file="IC_phyloseq_unfiltered.RData")
