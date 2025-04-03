
##Bioinformatic Decontamination

library(tidyverse)
library(here)
library(phyloseq)
library(decontam)
library(vegan)
library(ggpubr)
library(microViz)

# using the Phyloseq object produced from DADA2 code 
IC_phyloseq_unfiltered

##Create subset with only mock communities
mock_phyloseq <-subset_samples(IC_phyloseq_unfiltered, SampleType=="Mock")

mock_phyloseq_norm <- transform_sample_counts(mock_phyloseq,function(x) 100* x/sum(x))

#mock community includes 8 bacteria
#sort the 8 most abundant ASVs in the mock phyloseq
sorted <- head(sort(colSums(otu_table(mock_phyloseq)), decreasing = TRUE), 8) 
mock_taxa <- cbind(as.data.frame(tax_table(mock_phyloseq)[names(sorted),]), Count = sorted)

##Create subset of data with only mock and negative controls  
mock_blanks <- subset_samples(IC_phyloseq_unfiltered, !SampleType=="Urine")

mock_blanks_norm <- transform_sample_counts(mock_blanks,function(x) 100* x/sum(x))

##extract out original ASV table
original_ASVs <- as.matrix(as.data.frame(otu_table(mock_blanks_norm)))

##identify the original proportion of contaminants
contaminants_original <- rowSums(original_ASVs[,!colnames(original_ASVs) %in% rownames(mock_taxa)])

##identify the original proportion of mock community ASVs
mock_original <- rowSums(original_ASVs[,colnames(original_ASVs) %in% rownames(mock_taxa)])

##Evaluate different thresholds of Decontam
contam_Mocks <- isNotContaminant(mock_blanks, neg = "is.neg", method = "prevalence", threshold = 0.5, normalize = TRUE, detailed = TRUE)

##create subset of data excluding mock community samples
Patient_Blanks <- subset_samples(IC_phyloseq_unfiltered, !SampleType =="Mock")

##apply Decontam
IC_contam_Patient <- isNotContaminant(Patient_Blanks, neg = "is.neg", method ="prevalence",threshold = 0.5, normalize = TRUE, detailed = TRUE)

##how many ASVs are not contaminants 
sum(IC_contam_Patient$not.contaminant)

##create physloseq object without contaminants
IC_patient_true_blanks <- prune_taxa(IC_contam_Patient$not.contaminant, Patient_Blanks)

##remove chloroplast, mitochondrial, and unspecified taxa
IC_patient_true_final <-subset_taxa(IC_patient_true_blanks, Kingdom == "Bacteria"& Family != "Mitochondria" & Genus != "Mitochondria" & Order!="Chloroplast")

IC_patient_true_final <- subset_samples(IC_patient_true_final, !is.neg =="TRUE")

final_tax <- tax_table(IC_patient_true_final)
final_ASV <- otu_table(IC_patient_true_final)
FinalASVs <- cbind(final_tax,t(final_ASV))

write.csv(FinalASVs,file="FinalASVs_taxonomy.csv")
## This is Table S2
