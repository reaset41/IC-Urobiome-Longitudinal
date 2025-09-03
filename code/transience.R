# Load required libraries
library(dplyr)
library(tidyr)

#ASV abundance data (rows: ASVs, columns: samples)
#subject column indicates which subject each sample belongs to
#file TableASVs_rel_abund_sorted_transience.csv available on GitHub
#file genera_ASV.csv available on GitHub
set.seed(123)  # For reproducibility

asv_data <- data.frame(t(TableASVs_rel_abund_sorted_transience))

#Convert to long format for easier manipulation
#if ASV relative abundance is >0, convert to Presence=1
asv_long <- asv_data %>%
  pivot_longer(
    cols = starts_with("ASV"), 
    names_to = "ASV", 
    values_to = "Presence"
  ) %>%
  mutate(Presence = ifelse(Presence > 0, 1, 0))


# Summarize persistence across three samples per subject
asv_transience <- asv_long %>%
  group_by(SubjectID, ASV) %>%
  summarise(
    Total_Present = sum(Presence),  # Number of samples where ASV is present
    .groups = "drop"
  ) %>%
  mutate(Transience = case_when(
    Total_Present == 3 ~ "Persistent",
    Total_Present == 2 ~ "Persistent",
    Total_Present == 0 ~ "Absent",
    TRUE ~ "Transient"
  ))

#calculate overall frequency of persistence for each ASV across subjects

asv_transience_summary <- asv_transience %>%
  group_by(ASV, Transience) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Transience, values_from = Count, values_fill = 0) %>%
  mutate(
    n_subjects = Transient + Persistent,
    frequency = Persistent / n_subjects
  )

#view the transience summary
print(asv_transience_summary)

#limit to ASVs found in more than 5 subjects
#label with genera of ASV and average relative abundance
#file genera_ASV.csv is available on GitHub
asv_transience_filtered <- asv_transience_summary %>%
  filter(n_subjects > 5) %>%
  inner_join(genera_ASV, by = "ASV")

#this is Supplemental Table 8