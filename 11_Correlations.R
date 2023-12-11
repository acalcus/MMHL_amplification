require(ggplot2)
require(tidyverse)
require(tidyr)
require(dplyr)
require(nlme)
require(lme4)
require(lsmeans)
require(MuMIn)
require(plyr)
require(readxl)
require(cowplot)
require(ppcor)
require(GGally)

# January 2023
# clear workspace
rm(list = ls())

# load data file
setwd("~/Dropbox/Dossier de l'équipe InMignonetteWeTrust/UCL/iCARE project/stats/")
d<-read_xls("11MMHL_Correlations.xls")
MM<-subset(d, group == "MM")
CA<-subset(d, group == "CA")

###############################
###### SPEECH PERCEPTION IN NOISE
###############################
# Correlation: SiQ by E
cor.test(MM$quiet, MM$E_F0_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(MM$quiet, MM$E_F0_ampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(CA$quiet, CA$E_F0_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(MM$quiet, MM$E_H2_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(MM$quiet, MM$E_H2_ampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(CA$quiet, CA$E_H2_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)

# Correlation: SIN by TFS
cor.test(MM$MR, MM$TFS_H4_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(MM$MR, MM$TFS_H4_ampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(CA$MR, CA$TFS_H4_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)

cor.test(MM$steady, MM$TFS_H4_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(MM$steady, MM$TFS_H4_ampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(CA$steady, CA$TFS_H4_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)

cor.test(MM$fluct, MM$TFS_H4_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(MM$fluct, MM$TFS_H4_ampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(CA$fluct, CA$TFS_H4_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)

cor.test(MM$MR, MM$TFS_H6_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(MM$MR, MM$TFS_H6_ampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(CA$MR, CA$TFS_H6_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)

cor.test(MM$steady, MM$TFS_H6_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(MM$steady, MM$TFS_H6_ampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(CA$steady, CA$TFS_H6_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)

cor.test(MM$fluct, MM$TFS_H6_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(MM$fluct, MM$TFS_H6_ampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(CA$fluct, CA$TFS_H6_unampl_Peak, use = "pairwise.complete.obs", na.action = na.exclude)

# Cortical: P1/MMN
cor.test(MM$steady, MM$P1_unamplified, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(MM$steady, MM$P1_amplified, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(CA$steady, CA$P1_unamplified, use = "pairwise.complete.obs", na.action = na.exclude)

cor.test(MM$steady, MM$MMN_unamplified, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(MM$steady, MM$MMN_amplified, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(CA$steady, CA$MMN_unamplified, use = "pairwise.complete.obs", na.action = na.exclude)

###############################
###### SPEECH DISCRIMINATION (QUIET)
###############################
# Correlation: SpeechDiscrim by MMN amplitude
cor.test(MM$SpeechDiscrim, MM$MMN_unamplified, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(MM$SpeechDiscrim, MM$MMN_amplified, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(CA$SpeechDiscrim, CA$MMN_unamplified, use = "pairwise.complete.obs", na.action = na.exclude)

###############################
###### SUBCORTICAL <-> CORTICAL
###############################
# Correlation: SIN by TFS
cor.test(MM$E_F0_unampl_Peak, MM$P1_unamplified, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(MM$E_F0_ampl_Peak, MM$P1_amplified, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(CA$E_F0_unampl_Peak, CA$P1_unamplified, use = "pairwise.complete.obs", na.action = na.exclude)

cor.test(MM$E_F0_unampl_Peak, MM$MMN_unamplified, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(MM$E_F0_ampl_Peak, MM$MMN_amplified, use = "pairwise.complete.obs", na.action = na.exclude)
cor.test(CA$E_F0_unampl_Peak, CA$MMN_unamplified, use = "pairwise.complete.obs", na.action = na.exclude)

###############################
###### SUBCORTICAL <-> SUBCORTICAL
###############################
CA_pairs <- CA[c(11, 13, 15, 17, 19, 21)]
MM_pairs_unamplified <- MM[c(11, 13, 15, 17, 19, 21)]
MM_pairs_amplified <- MM[c(12, 14, 16, 18, 20, 22)]
ggpairs(CA_pairs)
ggpairs(MM_pairs_unamplified)
ggpairs(MM_pairs_amplified)


