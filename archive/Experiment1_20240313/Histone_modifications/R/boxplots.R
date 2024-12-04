## Head

## Load libraries 
library(tidyverse)
library(ggplot2)

## vectorize filenames
filenames <- c("H3_H4.csv", "H3K9ac_H3K27me3.csv", "H3K27ac.csv", "H3K36me3_H2AXp139.csv", 
               "H3K56ac_H1.csv", "H4K8ac_H3K9me3.csv", "H4K16ac.csv", "H4K20me3_H4K12ac.csv")
filenames <- paste0("Data/", filenames)

## loop over filenames to load data
files <- lapply(filenames, function(file){
  read.csv(file)
})

names(files) <- str_remove(str_extract(filenames, "H.+[\\dac]\\."), "\\.")

## function to tidy up data
tidy_histone_data <- function(data, histone){
  hist_string <- paste0("Integrated", histone)
  tmp <- select(df, "Dosage", "Time_Point", "IR_Pattern", "IntegratedDNA", all_of(hist_string))
  tmp$Histone <- histone
  tmp$grpMean <- mean(tmp[[hist_string]][tmp$Dosage == "UnIR"] / tmp$IntegratedDNA[tmp$Dosage == "UnIR"])
  tmp$standardized <- (tmp[[hist_string]]/tmp$IntegratedDNA) / tmp$grpMean
  tmp <- rename(tmp, Integrated = all_of(hist_string))
  
  ds <- unique(tmp$Dosage)
  n <- length(ds)
  
  df <- tmp[tmp$Time_Point == "UnIR",]
  df2 <- do.call("rbind", replicate(n, df, simplify = FALSE))
  ds_new <- lapply(1:n, function(i){
    rep(ds[i],dim(df)[1])
  })
  df2$Dosage <- do.call("c", ds_new)
  rbind(tmp, df2)
}

## extract and tidy data for each respective histone
dat_H3 <- tidy_histone_data(files$H3_H4, histone = "H3")
dat_H4 <- tidy_histone_data(files$H3_H4, histone = "H4")
dat_H3K9ac <- tidy_histone_data(files$H3K9ac_H3K27me3, histone = "H3K9ac")
dat_H3K27me3 <- tidy_histone_data(files$H3K9ac_H3K27me3, histone = "H3K27me3")
dat_H3K27ac <- tidy_histone_data(files$H3K27ac, histone = "H3K27ac")
dat_H3K36me3 <- tidy_histone_data(files$H3K36me3_H2AXp139, histone = "H3K36me3")
dat_H2AXp139 <- tidy_histone_data(files$H3K36me3_H2AXp139, histone = "H2AXp139")
dat_H3K56ac <- tidy_histone_data(files$H3K56ac_H1, histone = "H3K56ac")
dat_H1 <- tidy_histone_data(files$H3K56ac_H1, histone = "H1")
dat_H4K8ac <- tidy_histone_data(files$H4K8ac_H3K9me3, histone = "H4K8ac")
dat_H3K9me3 <- tidy_histone_data(files$H4K8ac_H3K9me3, histone = "H3K9me3")
dat_H4K16ac <- tidy_histone_data(files$H4K16ac, histone = "H4K16ac")
dat_H4K20me3 <- tidy_histone_data(files$H4K20me3_H4K12ac, histone = "H4K20me3")
dat_H4K12ac <- tidy_histone_data(files$H4K20me3_H4K12ac, histone = "H4K12ac")

## bind data into one data.frame
dat <- rbind(dat_H3, dat_H4, dat_H3K9ac, dat_H3K27me3, dat_H3K27ac, dat_H3K36me3,dat_H2AXp139, dat_H3K56ac,dat_H1, dat_H4K8ac,
      dat_H3K9me3, dat_H4K16ac, dat_H4K20me3,dat_H4K12ac)

## reorder dosage and time points
dat$Dosage <- factor(dat$Dosage, levels = c("UnIR", "UVA 8 kJ", "UVA 24 kJ", "3x UVA 8 kJ", "UVB 300 J", "UVB 900 J", "3x UVB 300 J"))
dat$Time_Point <- factor(dat$Time_Point,levels = c("UnIR","0.5 H", "2 H", "4 H", "8 H", "16 H", "24 H", "48 H"),
                         labels = c("UnIR","0.5h", "2h", "4h", "8h", "16h", "24h", "48h"))

## loop over plot list 
histones <- unlist(str_split(names(files), "_"))
plotlist <- lapply(histones, function(histone){
  p <- dat %>% 
    filter(Dosage != "UnIR" & Histone == histone) %>%
    ggplot(aes(x = Time_Point, y = standardized, fill = Time_Point)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~Dosage, scales = "free", nrow = 2) +
    geom_hline(yintercept = 1, lty = 2) +
    theme_classic() +
    labs(y = paste0(histone," / DNA (Standardized to UnIR)")) +
    theme(axis.text.x = element_text(size = 8, face = "bold", angle = 45, hjust = 1, vjust = 1),
          strip.background = element_blank(),
          legend.position = "none") +
    coord_cartesian(ylim = c(0, 5)) +
    scale_fill_brewer(palette = "Greys")
  
  png(paste0("Results/boxplot_", histone,".png"), units="in", width=10, height=6.2, res=600)
  print(p)
  dev.off()
})
names(plotlist) <- histones


