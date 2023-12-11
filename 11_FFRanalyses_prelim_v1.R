require(nlme)
require(pastecs)
require(lme4)
require(MuMIn)
require(ggplot2)
require(graphics)
require(readxl)
require(lmerTest)
require(tidyverse)
require(magrittr)
require(varhandle)
require(emmeans)
require(xlsx)
require(naniar)
require(lawstat)
require(car)
require(dplyr)
require(stringr)

# April 2023
# clear workspace
rm(list = ls())

# load data file
setwd("~/Dropbox/Dossier de l'équipe InMignonetteWeTrust/UCL/iCARE project/stats/")
d<-read_xlsx("11MMHL_FFRdata_Cz.xlsx") %>% 
  mutate(SNR=Peak/Baseline, SNRdB=20*log10(SNR), SNdiff=Peak-Baseline) %>% 
  mutate(across(where(is.character), factor))


# filter out outliers by calculating z scores (ZPeak & ZBaseline),
# putting the original data values into oPeak & oBaseline
scaleSR <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}
dZ <- d %>% 
  group_by(New_Group, Component, Harmonic) %>% # 
  mutate(ZPeak = scaleSR(Peak), ZBaseline = scaleSR(Baseline)) %>% 
  ungroup() %>% 
  mutate(oPeak=Peak, oBaseline=Baseline) %>% 
  mutate(Peak = ifelse(abs(ZPeak) >=3, NA, Peak))%>%
  mutate(Baseline = ifelse(abs(ZBaseline) >=3, NA, Baseline))


# Which harmonics contain significant responses ?
# Pivot longer to put Peak & Baseline in separate rows, but first filter out any rows with either Peak or Baseline missing. Then run paired-samples t-tests on all 36 combinations of Component & Harmonic with Holm correction
dZL <- dZ %>% 
  filter(!is.na(Peak) & !is.na(Baseline)) %>% 
  select(code, age, Group, Amplification, New_Group, Component, Harmonic, Peak, Baseline) %>% 
  pivot_longer(cols=c("Peak", "Baseline"), names_to="SpectralComp", values_to = "Ampl") %>% 
  filter(!is.na(Ampl))

tTR <- dZL %>%
  group_by(New_Group, Component, Harmonic) %>%
  t_test(Ampl ~ SpectralComp, paired = TRUE) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance() %>% 
  arrange(Component, Harmonic) %>% 
  select(!c(.y., group1, group2))

# Analyse baseline levels
# filter out Baseline = NA
dZB <- dZ %>% 
  filter(!is.na(Baseline)) 

mB.sat <- lmer(Baseline ~ 1 + (Group + Amplification) *Harmonic*Component + (1|code),
                   data = dZB, REML=FALSE, na.action=na.fail,
                   control = lmerControl(optCtrl = list(maxfun = 1e6)))
anova(mB.sat)

mB.01 <- lmer(Baseline ~ 1 + Harmonic*Component + (1|code),
                  data = dZB, REML=FALSE, na.action=na.fail,
                  control = lmerControl(optCtrl = list(maxfun = 1e6)))
anova(mB.01)

anova(mB.sat, mB.01)






###########################
####### Q3: Is the data normally distributed ?
###########################
NHU_E_F0 <- filter(dZ, New_Group == "NH-U" & Component == "E" & Harmonic == "F0")
NHU_E_H2 <- filter(dZ, New_Group == "NH-U" & Component == "E" & Harmonic == "H2")

ks.test(NHU_E_F0$Peak, "pnorm", mean = mean(NHU_E_F0$Peak), sd = sd(NHU_E_F0$Peak))
ks.test(NHU_E_H2$Peak, "pnorm", mean = mean(NHU_E_H2$Peak), sd = sd(NHU_E_H2$Peak))

NHU_TFS_H3 <- filter(dZ, New_Group == "NH-U" & Component == "TFS" & Harmonic == "H3")
NHU_TFS_H4 <- filter(dZ, New_Group == "NH-U" & Component == "TFS" & Harmonic == "H4")
NHU_TFS_H5 <- filter(dZ, New_Group == "NH-U" & Component == "TFS" & Harmonic == "H5")
NHU_TFS_H6 <- filter(dZ, New_Group == "NH-U" & Component == "TFS" & Harmonic == "H6")

ks.test(NHU_TFS_H3$Peak, "pnorm", mean = mean(NHU_TFS_H3$Peak, na.rm = TRUE), sd = sd(NHU_TFS_H3$Peak, na.rm = TRUE))
ks.test(NHU_TFS_H4$Peak, "pnorm", mean = mean(NHU_TFS_H4$Peak, na.rm = TRUE), sd = sd(NHU_TFS_H4$Peak, na.rm = TRUE))
ks.test(NHU_TFS_H5$Peak, "pnorm", mean = mean(NHU_TFS_H5$Peak, na.rm = TRUE), sd = sd(NHU_TFS_H5$Peak, na.rm = TRUE))
ks.test(NHU_TFS_H6$Peak, "pnorm", mean = mean(NHU_TFS_H6$Peak, na.rm = TRUE), sd = sd(NHU_TFS_H6$Peak, na.rm = TRUE))

HLU_E_F0 <- filter(dZ, New_Group == "HL-U" & Component == "E" & Harmonic == "F0")
HLU_E_H2 <- filter(dZ, New_Group == "HL-U" & Component == "E" & Harmonic == "H2")

ks.test(HLU_E_F0$Peak, "pnorm", mean = mean(HLU_E_F0$Peak), sd = sd(HLU_E_F0$Peak))
ks.test(HLU_E_H2$Peak, "pnorm", mean = mean(HLU_E_H2$Peak), sd = sd(HLU_E_H2$Peak))

HLU_TFS_H3 <- filter(dZ, New_Group == "HL-U" & Component == "TFS" & Harmonic == "H3")
HLU_TFS_H4 <- filter(dZ, New_Group == "HL-U" & Component == "TFS" & Harmonic == "H4")
HLU_TFS_H5 <- filter(dZ, New_Group == "HL-U" & Component == "TFS" & Harmonic == "H5")
HLU_TFS_H6 <- filter(dZ, New_Group == "HL-U" & Component == "TFS" & Harmonic == "H6")

ks.test(HLU_TFS_H3$Peak, "pnorm", mean = mean(HLU_TFS_H3$Peak, na.rm = TRUE), sd = sd(HLU_TFS_H3$Peak, na.rm = TRUE))
ks.test(HLU_TFS_H4$Peak, "pnorm", mean = mean(HLU_TFS_H4$Peak, na.rm = TRUE), sd = sd(HLU_TFS_H4$Peak, na.rm = TRUE))
ks.test(HLU_TFS_H5$Peak, "pnorm", mean = mean(HLU_TFS_H5$Peak, na.rm = TRUE), sd = sd(HLU_TFS_H5$Peak, na.rm = TRUE))
ks.test(HLU_TFS_H6$Peak, "pnorm", mean = mean(HLU_TFS_H6$Peak, na.rm = TRUE), sd = sd(HLU_TFS_H6$Peak, na.rm = TRUE))

HLA_E_F0 <- filter(dZ, New_Group == "HL-A" & Component == "E" & Harmonic == "F0")
HLA_E_H2 <- filter(dZ, New_Group == "HL-A" & Component == "E" & Harmonic == "H2")

ks.test(HLA_E_F0$Peak, "pnorm", mean = mean(HLA_E_F0$Peak), sd = sd(HLA_E_F0$Peak))
ks.test(HLA_E_H2$Peak, "pnorm", mean = mean(HLA_E_H2$Peak), sd = sd(HLA_E_H2$Peak))

HLA_TFS_H3 <- filter(dZ, New_Group == "HL-A" & Component == "TFS" & Harmonic == "H3")
HLA_TFS_H4 <- filter(dZ, New_Group == "HL-A" & Component == "TFS" & Harmonic == "H4")
HLA_TFS_H5 <- filter(dZ, New_Group == "HL-A" & Component == "TFS" & Harmonic == "H5")
HLA_TFS_H6 <- filter(dZ, New_Group == "HL-A" & Component == "TFS" & Harmonic == "H6")

ks.test(HLA_TFS_H3$Peak, "pnorm", mean = mean(HLA_TFS_H4$Peak, na.rm = TRUE), sd = sd(HLA_TFS_H3$Peak, na.rm = TRUE))
ks.test(HLA_TFS_H4$Peak, "pnorm", mean = mean(HLA_TFS_H4$Peak, na.rm = TRUE), sd = sd(HLA_TFS_H4$Peak, na.rm = TRUE))
ks.test(HLA_TFS_H5$Peak, "pnorm", mean = mean(HLA_TFS_H5$Peak, na.rm = TRUE), sd = sd(HLA_TFS_H5$Peak, na.rm = TRUE))
ks.test(HLA_TFS_H6$Peak, "pnorm", mean = mean(HLA_TFS_H6$Peak, na.rm = TRUE), sd = sd(HLA_TFS_H6$Peak, na.rm = TRUE))


#######################
write.csv(dZ, "11MMHL_FFRdata_tidy_v1.csv")



#########################
######## Levene's test
#########################
leveneTest(Peak ~ New_Group, data = subset(dZ, Component = "E" & Harmonic == "F0"), center = mean)
leveneTest(Peak ~ New_Group, data = subset(dZ, Component = "E" & Harmonic == "H2"), center = mean)
leveneTest(Peak ~ New_Group, data = subset(dZ, Component = "TFS"), center = mean)


