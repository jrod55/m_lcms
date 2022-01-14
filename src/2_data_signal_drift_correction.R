library(pseudoDrift)
library(data.table)
library(tidyverse)
rm(list = ls())

# Read in the pw_outlier cleaned and normalized data ----------------------
setwd("~/signal_drift/m_lcms/data/")
nc_dat = fread("01_pw_outlier_norm_peak_area_matrix.txt")

## rename area_cleaned column to area as pseudo_sdc expects
togB = nc_dat
colnames(togB)[12] = "area"

# Set params and run pseudo_sdc -------------------------------------------
togB_comps = togB %>%
  select(area, compound) %>%
  drop_na()
integrated_compounds = sort(unique(togB_comps$compound))

## Make lists where results will go
tog_corrected = list()
crit_table = list()
tog_corrected_plots = list()
source("../src/my_plot_function.R")
for (j in integrated_compounds) {

  ## Subset data
  tmp1 = togB %>%
    filter(compound == all_of(j)) %>%
    mutate(batch = paste0("B",batch)) %>%
    mutate(area = log2(area+1))

  ## Set training batch
  training_batch = "B4"

  ## QC frequency in training batch
  qc_freq = tmp1 %>%
    filter(batch==all_of(training_batch)) %>%
    group_by(compound, sample, batch) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(compound, batch) %>%
    mutate(Freq = sum(n)/n) %>%
    filter(sample=="Pool")

  ## Count number of pools
  n_pools = sum(qc_freq$n)

  if (n_pools<2) {
    print(paste0("Not enough pools to do signal drift correction with ", j, " data"))
  }else{

    ## Mean frequency of pools in training batch
    mnqc = round(mean(qc_freq$Freq))*2

    ## params
    ttt = tmp1 %>%
      filter(compound==all_of(j),
             !sample=="Pool") %>%
      group_by(batch) %>%
      summarise(samps = n_distinct(name))

    s_perBatch = min(ttt$samps)

    ## Test breaks in smallest batch
    t.b = round(s_perBatch/mnqc)

    ## Proportional in large training set
    ttt = ttt %>%
      mutate(bre = round(samps*all_of(t.b)/all_of(s_perBatch)))

    t.b.train = ttt %>%
      filter(batch == all_of(training_batch)) %>%
      pull(bre)

    t.b.train = seq(t.b.train, t.b.train+5, 1)

    ## Test window for median smoothing
    w.n.max = round(min(s_perBatch/t.b))/2
    w.n = 2*floor(w.n.max/2)+1
    if (w.n>w.n.max) {
      w.n = w.n-2
    }
    w.n = seq(1,w.n,2)

    ## Test index offset
    ti.max = round(s_perBatch*0.05)
    t.i = seq(0,ti.max,1)

    sdc_out = pseudo_sdc(df = tmp1,
                         n.cores = 18,
                         train.batch = "B4",
                         criteria = "MSE",
                         test.breaks = t.b.train,
                         test.window = w.n,
                         test.index = t.i,
                         qc.label = "Pool",
                         min.qc = min(t.b.train),
                         log_transform = FALSE,
                         quantile.increment = 0.10
                         )

    tog_corrected[[j]] = sdc_out[[3]]
    tog_corrected_plots[[j]] = cplt(sdc_out, dat_train = "B4", testing = "AUA")
    crit_table[[j]] = sdc_out[[4]] %>%
      filter(criteria == "MSE") %>%
      rename(MSE = value) %>%
      select(-c(criteria))

    ggsave(tog_corrected_plots[[j]], file=paste0("../plots/",j,"_AUA_sdc.png"), dpi = 300,width = 12, height = 7.5, units = "in",bg = "white")
  }
}

dat_corrected = bind_rows(tog_corrected)
dat_criteria = bind_rows(crit_table)

fwrite(dat_corrected, file = "02_signalDrift_corrected_peak_area_matrix.txt", col.names = TRUE, sep = "\t", quote = FALSE, na = NA)
fwrite(dat_criteria, file = "../out/criteria_table_02_signalDrift_corrected_peak_area_matrix.txt", col.names = TRUE, sep = "\t", quote = FALSE, na = NA)

print(paste0("Done with signal drift correction"))

# convert ../out/*.png -resize 1200x750 ../out/00_maize_phenolic_compounds.pdf




