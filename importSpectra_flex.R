library(AssaySupport)
library(tidyverse)
library(MALDIquant)

parentDir <- "Y:/1-DrArbeitPhD/1 - phd/32-Thomas Enzlein/Target/181106_MK_activity_assays/20181019_Substrates_37_45/"
namePrefix <- "MK20181019"
tol = 8

Ablist <- readxl::read_xlsx("Y:/1-DrArbeitPhD/1 - phd/32-Thomas Enzlein/Target/Ab_mz_diff_mutation.xlsx")

spectra <- load_spectra(Dir = parentDir)

# preprocess
spectra_align <- preprocess_spectra(spec = spectra,
                                    smooth_meth = "SavitzkyGolay",
                                    smooth_halfWindowSize  = 5,
                                    norm_meth = "TIC",
                                    filter_spectra = "Cal|D1|D2|CAL",
                                    baseline_meth = "TopHat",
                                    align = FALSE)

write_msd(spectra = spectra_align, parentDir = parentDir, namePrefix = namePrefix)

# pick peaks
cat("\n", AlzTools::timeNowHM(), "Processing peaks...\n")
peaks <- MALDIquant::detectPeaks(spectra_align, method = "SuperSmoother", SNR = 3)
peak_df <- generate_peakdf(spectra_align, pick_meth = "SuperSmoother", SNR = 3, binpeaks = FALSE) 

# generate result data.frame 
cat("\n", AlzTools::timeNowHM(), "Generating result data.frame...\n")
label_df <- peak_df %>%
  separate(ID, into = c("Substrate", "Type", "Condition"), sep = "_") %>%
  assign_species(., Ablist, tol = tol, mzcol = "mz")

res_df <- label_df %>%
  filter(!is.na(species)) %>%
  group_by(Substrate, Condition, Type) %>%
  mutate(total.Ab = sum(int[!species %in% c("AICD49_99", "AICD50_99")]),
         norm.tot = int/total.Ab*100)

# plot
cat("\n", AlzTools::timeNowHM(), "Plotting...\n")

# normalized Abeta barplot
if(FALSE) {
normAb_p <- res_df %>%
  filter(!species %in% c("AICD49_99", "AICD50_99", "Substrate", "Ab1_37"),
         !Condition %in% c("Substrate")) %>%
  ggplot( aes(x = species, y = norm.tot, fill = Condition)) +
  geom_col(col = "black", position = "dodge") +
  facet_wrap(~Substrate, nrow = 1) +
  #coord_flip() +
  AlzTools::xtheme() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Abeta species", y = "Rel. Intensity as % \nof total Abeta")
normAb_p
ggsave(plot = normAb_p,
       path = parentDir,
       filename = paste0(namePrefix, "_",basename(parentDir), "_normAb_plot.png"),
       dpi = 300, width = 20, height = 10, units = "cm")

for(subst in unique(res_df$Substrate)) {
normAb_p_single <- res_df %>%
  filter(!species %in% c("AICD49_99", "AICD50_99", "Substrate", "Ab1_37"),
         !Condition %in% c("Substrate"),
         Substrate == subst) %>%
  ggplot( aes(x = species, y = norm.tot, fill = Condition)) +
  geom_col(col = "black", position = "dodge") +
  #coord_flip() +
  AlzTools::xtheme() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Abeta species", y = "Rel. Intensity as % \nof total Abeta", 
       title = subst)
ggsave(plot = normAb_p_single,
       path = parentDir,
       filename = paste0(namePrefix, "_",basename(parentDir), "_normAb_plot_", subst, ".png"),
       dpi = 300, width = 16, height = 12, units = "cm")
}
}
# total Abeta plot
if(FALSE) {
label_df %>%
  filter(!is.na(species),
         !species %in% c()) %>%
  group_by(Substrate, Condition, Type) %>%
  mutate(total.Ab = sum(int[!species %in% c("AICD49_99", "AICD50_99")])) %>%
  group_by(Substrate, Condition) %>%
  dplyr::summarise(total.Ab = mean(total.Ab)) %>%
  ungroup() %>%
  mutate(norm.total.Ab = total.Ab/max(total.Ab)*100) %>%
  ggplot( aes(x = Condition, y = norm.total.Ab, fill = Substrate)) +
  geom_col(col = "black", position = "dodge") +
  #coord_flip() +
  AlzTools::xtheme() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Temperature [°C]", y = "Total Abeta as % of max(totalAb)")
}

# AICD line plot
if(FALSE) {
AICD_p <- label_df %>%
  filter(!is.na(species),
         species %in% c("AICD49_99", "AICD50_99", "Substrate"),
         !Condition == "X") %>%
  select(Type, Condition, Substrate, int, species) %>%
  spread(species, int) %>%
  group_by(Type) %>%
  mutate(norm.AICD49_99 = 100*AICD49_99/AICD49_99[Condition == "DMSO"],
         norm.AICD50_99 = 100*AICD50_99/AICD50_99[Condition == "DMSO"],
         norm.Substrate = 100*Substrate/Substrate[Condition == "DMSO"],
         Conc = c(0,1,10,100)) %>%
  select(-AICD49_99, -AICD50_99, - Substrate) %>%
  ungroup() %>%
  rename(AICD49_99 = norm.AICD49_99,
         AICD50_99 = norm.AICD50_99,
         Substrate = norm.Substrate) %>%
  gather(norm, value, -Type, -Condition, - Conc) %>%
  ggplot(aes(x = Conc, y = value, col = Type)) +
  facet_wrap(~norm) +
  geom_line() +
  AlzTools::xtheme() +
  labs(x = "Conc. Wagner Compound [µM]", y = "Rel. Intensity as % \nof CTRL")
AICD_p
ggsave(plot = AICD_p,
       path = parentDir,
       filename = paste0(namePrefix, "_", basename(parentDir),"_normAICD_plot.png"),
       dpi = 300, width = 20, height = 10, units = "cm")
}


# AICD-ratio line plot
if(FALSE) {
AICD_ratio_p <- label_df %>%
  filter(!is.na(species),
         species %in% c("AICD49_99", "AICD50_99", "Substrate"),
         !Condition == "InhX") %>%
  select(Type, Condition, Substrate, int, species) %>%
  spread(species, int) %>%
  group_by(Type) %>%
  mutate("AICD49_99/AICD50_99.ratio" = AICD49_99/AICD50_99,
         Conc = c(0,1,10,100)) %>%
  select(-AICD49_99, -AICD50_99, - Substrate) %>%
  ungroup() %>%
  gather(norm, value, -Type, -Condition, - Conc) %>%
  ggplot(aes(x = Conc, y = value, col = Type)) +
  geom_line() +
  AlzTools::xtheme() +
  labs(x = "Conc. Wagner Compound [µM]", y = "AICD49-99/AICD50-99-ratio")
AICD_ratio_p
ggsave(plot = AICD_ratio_p,
       path = parentDir,
       filename = paste0(namePrefix, "_","AICDratio_plot.png"),
       dpi = 300, width = 12, height = 12, units = "cm")
}
# spectra plots, overview and conditions
for(i in 1:length(spectra_align)){
  muts <- c(1, 4, 7, 10, 13, 16, 19, 22, 25)
  params <- c(3, 3, 3, 3, 3, 3, 3, 3, 3)
  # png(file = file.path(parentDir, paste0(namePrefix, "_", names(spectra_align[i]), "_overview.png")),
  #     width = 1200, height = 1000,
  #     units = "px",
  #     pointsize = 7,
  #     res = 250)
  # overviewplot(idx = i, overview_adjust = 1)
  # dev.off()
  if(i %in% muts) {
    j <- which(muts == i)

    ggsave(plot = ggoverlay_spectra(spec_ID = i, n_overlay = params[j], normalizing = NA, xlim = c(3500, 5000)),
           path = parentDir,
           filename = paste0(namePrefix, "_", names(spectra_align[i]), "_abetaregion.png"),
           dpi = 300, width = 18, height = 12, units = "cm")
    # ggsave(plot = ggoverlay_spectra(spec_ID = i, n_overlay = params[j], normalizing = "max"),
    #        path = parentDir,
    #        filename = paste0(namePrefix, "_", names(spectra_align[i]), "_overlay_maxNorm.png"),
    #        dpi = 300, width = 12, height = 12, units = "cm")
    # ggsave(plot = ggoverlay_spectra(spec_ID = i, n_overlay = params[j], normalizing = 4138, tol = 3),
    #        path = parentDir,
    #        filename = paste0(namePrefix, "_", names(spectra_align[i]), "_overlay_ContamintionNorm.png"),
    #        dpi = 300, width = 12, height = 12, units = "cm")
    ggsave(plot = ggoverview(spec_ID = i, n_overlay = params[j], abeta_xlim = c(3500, 5000)),
           path = parentDir,
           filename = paste0(namePrefix, "_", names(spectra_align[i]), "_overview.png"),
           dpi = 300, width = 12, height = 12, units = "cm", scale = 1.75)
    

    # substrate <- strsplit(x = names(spectra_align[i]), "_")[[1]][1]
    # gsec <- strsplit(x = names(spectra_align[i]), "_")[[1]][2]
    # 
    # png(file = file.path(parentDir, paste0(namePrefix, "_", substrate,"_", gsec, "_conditions.png")),
    #     width = 1200, height = 1000,
    #     units = "px",
    #     pointsize = 7,
    #     res = 250)
    # mutation_plot(idx = i, plot_n = params[j])
    # dev.off()
  }

}

write_peak_list(labels_df = res_df, parentDir = parentDir, namePrefix = namePrefix)

cat("\n", AlzTools::timeNowHM(), "Done!\n")


