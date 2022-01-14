library(data.table)
library(tidyverse)
library(pseudoDrift)
rm(list = ls())

# Read in raw data -----------------------------------------------------
setwd("~/signal_drift/m_lcms/data/")
dm_og = fread("00_raw_peak_area_matrix.txt")
dm_meta = names(dm_og)[1:15]
dm = dm_og %>%
  pivot_longer(!c(all_of(dm_meta)), names_to = "compound", values_to = "area")
iss = c("IS_NEGATIVE","IS_POSITIVE")

# Limits of detection (LOD) by sub-batch -------------------------------
bbb = dm %>%
  filter(grepl("Blank",sample)) %>%
  group_by(batch, compound) %>%
  mutate(my_mad = mad(area, na.rm = TRUE),
         high_thresh = median(area, na.rm = TRUE) + 3*my_mad,
         low_thresh = median(area, na.rm = TRUE) - 3*my_mad) %>%
  filter(!(area<low_thresh | area>high_thresh)) %>%
  group_by(sub_batch, compound) %>%
  summarise(area_blank = mean(area),
            LOD = mean(area)*3) %>%
  mutate(LOD = ifelse(compound%in%all_of(iss), NA, LOD))
bbb_write = bbb %>%
  filter(!compound%in%iss)

## Setting values below LOD thresholds to NA
dm = dm %>%
  left_join(., bbb) %>%
  mutate(area_raw = area) %>%
  group_by(batch, compound) %>%
  mutate(LOD = ifelse(is.na(LOD), median(LOD, na.rm = TRUE), LOD),
         area = ifelse(area<=LOD, NA, area),
         area = ifelse(compound%in%all_of(iss), area_raw, area))

## Make df which will remain unchanged
dm_out = dm_long = dm

# Matrix effect calc for batch 4 ------------------------------------------

## Calculate matrix effect where internal standard is present
dm_write = dm_out %>%
  mutate(matrix_effect = ifelse(compound%in%all_of(iss), (area_blank-area_raw)/area_raw, NA))

## Identify samples with matrix effect >= 100. This would suggest ~90% of ions are lost in these samples
sme = dm_write %>%
  filter(matrix_effect>=100) %>%
  pull(name)

# Area wt normalization to give AUA ---------------------------------------
dm_write = dm_write %>%
  mutate(weight = ifelse(is.na(weight), mean(weight, na.rm = TRUE), weight),
         volume_solvent = ifelse(is.na(volume_solvent), 1e3, volume_solvent),
         volume_solvent = volume_solvent*1e-3,
         volume_solventL = volume_solvent*1e-3,
         tmp_vw = weight/volume_solvent,
         area_norm = (area/tmp_vw))

## Remove blanks and standards (internal and external) from table and set all normalized values to NA for samples with matrix effect >= 100 identified above
dm_write = dm_write %>%
  filter(!grepl(paste("Blank", "Standard", sep = "|"), sample)) %>%
  mutate(area_norm = ifelse(name%in%sme, NA, area_norm))

## Remove unnecessary columns
dout  = dm_write %>%
  select(-c(sample_type, plate, plate_decoder, weight,
            volume_solvent, area_blank,
            LOD, volume_solventL, tmp_vw, area, matrix_effect))

## Count the number of samples per compound and proportion missing data
samp_check = dout %>%
  group_by(compound, batch, sub_batch) %>%
  mutate(NsubB = n(),
         NsubB_missing = sum(is.na(area_norm))/n()) %>%
  group_by(compound, batch) %>%
  mutate(NB = n(),
         NB_missing = sum(is.na(area_norm))/n()) %>%
  group_by(compound) %>%
  mutate(N = n(),
         N_missing = sum(is.na(area_norm))/n()) %>%
  select(batch, sub_batch, compound, N, N_missing, NB, NB_missing, NsubB, NsubB_missing) %>%
  distinct()

## Add column to LOD with proportion missing per batch
lod_prop = dm_write %>%
  filter(!compound%in%all_of(iss)) %>%
  select(c(batch, sub_batch, compound, area_blank, LOD)) %>%
  left_join(., samp_check) %>%
  distinct() %>%
  mutate(N_missing = round(N_missing, 2),
         NB_missing = round(NB_missing, 2),
         NsubB_missing = round(NsubB_missing, 2))

## Make wide table with LOD thresholds by sub-batch
lod_prop_wide = lod_prop %>%
  mutate(LOD = round(LOD, 3)) %>%
  pivot_wider(!c(batch, area_blank, N, N_missing, NB, NB_missing, NsubB, NsubB_missing), names_from = sub_batch, values_from = LOD)

fwrite(lod_prop, "../out/lod_and_missing_data.txt", sep = "\t")
# fwrite(lod_prop_wide, "../out/lod_wide.txt", sep = "\t")


## Compounds with more than 25% missing data in any given batch, set all values for them to NA
all_comps = unique(dout$compound)
all_comps = all_comps[!all_comps%in%iss]
comps_outA = lod_prop %>%
  filter(NB_missing > 0.25) %>%
  pull(compound)
comps_withinA = all_comps[!all_comps%in%comps_outA]
dout = dout %>%
  mutate(area_norm = ifelse(compound%in%comps_withinA, area_norm, NA))

# Prep data for pw_outlier ------------------------------------------------
## Sample Metadata
m = dout
m_keep = c("injection_order", "name","sample","batch","compound", "rep", "rep_tech", "expr_rep","expr_flat","expr_pos")
s_meta = m %>%
  select(all_of(m_keep)) %>%
  distinct()

## Reformat how pw_outlier expects the data
d_dat = m %>%
  select(all_of(m_keep), area_norm) %>%
  group_by(batch, compound) %>%
  arrange(injection_order) %>%
  mutate(batch_index = 1:n(), .before = rep) %>%
  ungroup() %>%
  rename(experiment_index = injection_order,
         area = area_norm) %>%
  relocate(experiment_index, .before = batch_index)

# Pairwise outlier removal ------------------------------------------------
tmp = d_dat

pw_out = pw_outlier(
  df = tmp,
  samps_exclude = c("Pool"),
  n.cores = 18,
  return_plot = FALSE)

## pw_cleaned data
m1 = pw_out[[2]]
colnames(m1)[12] = "area_cleaned"

## data removed
drm = pw_out[[3]]

## write the files
fwrite(m1, file = "01_pw_outlier_norm_peak_area_matrix.txt", col.names = TRUE, sep = "\t", quote = FALSE, na = NA)
fwrite(drm, file = "../out/data_removed_00_raw_peak_area_matrix.txt", col.names = TRUE, sep = "\t", quote = FALSE, na = NA)
print(paste0("Done with data normalization and cleaning"))
