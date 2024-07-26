# Load packages
require(nlme)
require(pastecs)
require(lme4)
require(lsmeans)
require(MuMIn)
require(ggplot2)
require(car)
require(lmerTest)
require(tidyverse)
require(magrittr)
require(varhandle)
require(readxl)
require(lsr)
require(emmeans)
require(jtools)
require(performance)

# January 2023
# Updated July 2023 - P1 latency + N2 (amplitude & latency)
# Updated June 2024 - following reviewers' comments
# clear workspace
rm(list = ls())

# load data file
setwd("~/Dropbox/Dossier de l'équipe InMignonetteWeTrust/UCL/iCARE project/stats/")
d<-read_xlsx("11MMHL_CTX_long_Fz.xlsx", sheet = "Sheet3")
d$Group<-as.factor(d$Group)

############################################################
####### TIDY DATASET ? -> no need
# Get outliers
scaleSR <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}
dZ <- d %>% 
  group_by(New_Group) %>% # 
  mutate(ZMMN_amplitude = scaleSR(MMN_amplitude), ZP1_amplitude = scaleSR(P1_amplitude)) %>% 
  mutate(ZP1_latency = scaleSR(P1_latency), ZN2_amplitude = scaleSR(N2_amplitude)) %>% 
  mutate(ZN2_latency = scaleSR(N2_latency)) %>%
  ungroup() %>% 
  mutate(MMN_amplitude = ifelse(abs(ZMMN_amplitude) >=3, NA, MMN_amplitude))%>%
  mutate(P1_amplitude = ifelse(abs(ZP1_amplitude) >=3, NA, P1_amplitude)) %>%
  mutate(P1_latency = ifelse(abs(ZP1_latency) >=3, NA, P1_latency))%>%
  mutate(N2_amplitude = ifelse(abs(ZN2_amplitude) >=3, NA, N2_amplitude)) %>%
  mutate(N2_latency = ifelse(abs(ZN2_latency) >=3, NA, N2_latency))

###########################
####### Is the data normally distributed ? -> yes
###########################
THU <- filter(d, New_Group == "TH-U")
HLU <- filter(d, New_Group == "HL-U")
HLA <- filter(d, New_Group == "HL-A")

ks.test(THU$MMN_amplitude, "pnorm", mean = mean(THU$MMN_amplitude), sd = sd(THU$MMN_amplitude))
ks.test(THU$P1_amplitude, "pnorm", mean = mean(THU$P1_amplitude), sd = sd(THU$P1_amplitude))
ks.test(THU$P1_latency, "pnorm", mean = mean(THU$P1_latency), sd = sd(THU$P1_latency))
ks.test(THU$N2_amplitude, "pnorm", mean = mean(THU$N2_amplitude), sd = sd(THU$N2_amplitude))
ks.test(THU$N2_latency, "pnorm", mean = mean(THU$N2_latency), sd = sd(THU$N2_latency))

ks.test(HLU$MMN_amplitude, "pnorm", mean = mean(HLU$MMN_amplitude), sd = sd(HLU$MMN_amplitude))
ks.test(HLU$P1_amplitude, "pnorm", mean = mean(HLU$P1_amplitude), sd = sd(HLU$P1_amplitude))
ks.test(HLU$P1_latency, "pnorm", mean = mean(HLU$P1_latency), sd = sd(HLU$P1_latency))
ks.test(HLU$N2_amplitude, "pnorm", mean = mean(HLU$N2_amplitude), sd = sd(HLU$N2_amplitude))
ks.test(HLU$N2_latency, "pnorm", mean = mean(HLU$N2_latency), sd = sd(HLU$N2_latency))

ks.test(HLA$MMN_amplitude, "pnorm", mean = mean(HLA$MMN_amplitude), sd = sd(HLA$MMN_amplitude))
ks.test(HLA$P1_amplitude, "pnorm", mean = mean(HLA$P1_amplitude), sd = sd(HLA$P1_amplitude))
ks.test(HLA$P1_latency, "pnorm", mean = mean(HLA$P1_latency), sd = sd(HLA$P1_latency))
ks.test(HLA$N2_amplitude, "pnorm", mean = mean(HLA$N2_amplitude), sd = sd(HLA$N2_amplitude))
ks.test(HLA$N2_latency, "pnorm", mean = mean(HLA$N2_latency), sd = sd(HLA$N2_latency))

# Levene's test
leveneTest(MMN_amplitude ~ New_Group, data = d, center = mean)
leveneTest(MMN_latency ~ New_Group, data = d, center = mean)
leveneTest(P1_amplitude ~ New_Group, data = d, center = mean)
leveneTest(P1_latency ~ New_Group, data = d, center = mean)
leveneTest(N2_amplitude ~ New_Group, data = d, center = mean)
leveneTest(N2_latency ~ New_Group, data = d, center = mean)

############################################################
####### CTX ANALYSES
### Q1: NH-U vs HL-U (between subject comparison)
d_between <- subset(d, Amplification == "unamplified")

#--- P1
# Amplitude
SaturatedModel <- lm(P1_amplitude ~ Group*age, data = d_between)
m1 <- lm(P1_amplitude ~ Group+age, data = d_between)
m2 <- lm(P1_amplitude ~ Group, data = d_between)
anova(SaturatedModel, m1, m2)
AIC(SaturatedModel, m1, m2)
summary(m1)
performance(m1)

# Eta squared
m1 <- aov(P1_amplitude ~ Group*age, data = d_between)
etaSquared(m1, type = 2, anova = FALSE)

by(d_between$P1_amplitude, list(d_between$Group), stat.desc, basic = FALSE)

# Latency
SaturatedModel <- lm(P1_latency ~ Group*age, data = d_between)
m1 <- lm(P1_latency ~ Group+age, data = d_between)
m2 <- lm(P1_latency ~ Group, data = d_between)
anova(SaturatedModel, m1, m2)
AIC(SaturatedModel, m1, m2)
summary(m2)
performance(m2)

by(d_between$P1_latency, list(d_between$Group), stat.desc, basic = FALSE)

# Eta squared
m1 <- aov(P1_latency ~ Group*age, data = d_between)
etaSquared(m1, type = 2, anova = FALSE)

by(d_between$P1_latency, list(d_between$Group), stat.desc, basic = FALSE)


#--- N2
# Amplitude
SaturatedModel <- lm(N2_amplitude ~ Group*age, data = d_between)
m1 <- lm(N2_amplitude ~ Group+age, data = d_between)
m2 <- lm(N2_amplitude ~ Group, data = d_between)
anova(SaturatedModel, m1, m2)
AIC(SaturatedModel, m1, m2)
summary(m1)
performance(m1)

# Eta squared
m1 <- aov(N2_amplitude ~ Group*age, data = d_between)
etaSquared(m1, type = 2, anova = FALSE)

by(d_between$N2_amplitude, list(d_between$Group), stat.desc, basic = FALSE)

# Latency
SaturatedModel <- lm(N2_latency ~ Group*age, data = d_between)
m1 <- lm(N2_latency ~ Group+age, data = d_between)
m2 <- lm(N2_latency ~ Group, data = d_between)
anova(SaturatedModel, m1, m2)
AIC(SaturatedModel, m1, m2)
summary(m2)
performance(m2)

# Eta squared
m1 <- aov(N2_latency ~ Group*age, data = d_between)
etaSquared(m1, type = 2, anova = FALSE)

#--- MMN
# Amplitude
SaturatedModel <- lm(MMN_amplitude ~ Group*age, data = d_between)
m1 <- lm(MMN_amplitude ~ Group+age, data = d_between)
m2 <- lm(MMN_amplitude ~ Group, data = d_between)
anova(SaturatedModel, m1, m2)
AIC(SaturatedModel, m1, m2)
summary(m2)
performance(m2)

# Eta squared
m1 <- aov(MMN_amplitude ~ Group*age, data = d_between)
etaSquared(m1, type = 2, anova = FALSE)

by(d_between$MMN_amplitude, list(d_between$Group), stat.desc, basic = FALSE)

# Latency
SaturatedModel <- lm(MMN_latency ~ Group*age, data = d_between)
m1 <- lm(MMN_latency ~ Group+age, data = d_between)
m2 <- lm(MMN_latency ~ Group, data = d_between)
anova(SaturatedModel, m1, m2)
AIC(SaturatedModel, m1, m2)
summary(m2)
performance(m2)

# Eta squared
m1 <- aov(MMN_latency ~ Group*age, data = d_between)
etaSquared(m1, type = 2, anova = FALSE)


#######
### Q2: HL-U vs HL-A (within subject comparison)
d_within <- subset(d, Group == "HL")

#--- P1
# Amplitude
lm = lmer(P1_amplitude ~ Amplification*age*BEPTA + (1|code), data = d_within)
anova(lm)
# Model reduction
sm <- step(lm, reduce.random = F)
fm <- get_model(sm)
summary(fm)
performance(fm)

# Decompose the Amplification x BEPTA interaction
em1 <- emmeans(fm, pairwise ~ Amplification | BEPTA, adjust = "tukey")
print(em1)

# Eta squared
m1 <- aov(P1_amplitude ~ Amplification*age*BEPTA, data = d_within)
etaSquared(m1, type = 2, anova = FALSE)

by(d_within$P1_amplitude, list(d_within$Amplification), stat.desc, basic = FALSE)

# Latency
lm = lmer(P1_latency ~ Amplification*age*BEPTA + (1|code), data = d_within)
anova(lm)
# Model reduction
sm <- step(lm, reduce.random = F)
fm <- get_model(sm)
summary(fm)
performance(fm)

# Eta squared
m1 <- aov(P1_latency ~ Amplification+age+BEPTA, data = d_within)
etaSquared(m1, type = 2, anova = FALSE)

by(d_within$P1_latency, list(d_within$Amplification), stat.desc, basic = FALSE)

#--- N2
# Amplitude
lm = lmer(N2_amplitude ~ Amplification*age*BEPTA + (1|code), data = d_within)
anova(lm)
# Model reduction
sm <- step(lm, reduce.random = F)
fm <- get_model(sm)
summary(fm)
performance(fm)

# Eta squared
m1 <- aov(N2_amplitude ~ age, data = d_within)
etaSquared(m1, type = 2, anova = FALSE)

# Latency
lm = lmer(N2_latency ~ Amplification*age*BEPTA + (1|code), data = d_within)
anova(lm)
# Model reduction
sm <- step(lm, reduce.random = F)
fm <- get_model(sm)
summary(fm)
performance(fm)

# Decompose the Amplification x BEPTA interaction
em1 <- emmeans(fm, pairwise ~ Amplification | BEPTA, adjust = "tukey")
print(em1)

# Eta squared
m1 <- aov(N2_latency ~ Amplification*BEPTA, data = d_within)
etaSquared(m1, type = 2, anova = FALSE)

by(d_within$N2_latency, list(d_within$Amplification), stat.desc, basic = FALSE)

#--- MMN
# Amplitude
lm = lmer(MMN_amplitude ~ Amplification*age*BEPTA + (1|code), data = d_within)
anova(lm)
# Model reduction
sm <- step(lm, reduce.random = F)
fm <- get_model(sm)
summary(fm)
performance(fm)

# Eta squared
m1 <- aov(MMN_amplitude ~ age, data = d_within)
etaSquared(m1, type = 2, anova = FALSE)

# Linear regression on MMN Latency
lm = lmer(MMN_latency ~ Amplification*age*BEPTA + (1|code), data = d_within)
anova(lm)
# Model reduction
sm <- step(lm, reduce.random = F)
fm <- get_model(sm)
summary(fm)
performance(fm)

# Eta squared
m1 <- aov(MMN_latency ~ Amplification, data = d_within)
etaSquared(m1, type = 2, anova = FALSE)



#######
### Q3: NH-U vs HL-A (between subject comparison)
# MMN
d_CA <- subset(d, Group == "TH")
d_MM <- subset(d, Amplification == "amplified")
d_CAMM <- rbind(d_CA, d_MM)

#--- P1
# Amplitude
SaturatedModel <- lm(P1_amplitude ~ Group*age, data = d_CAMM)
m1 <- lm(P1_amplitude ~ Group+age, data = d_CAMM)
m2 <- lm(P1_amplitude ~ Group, data = d_CAMM)
anova(SaturatedModel, m1, m2)
AIC(SaturatedModel, m1, m2)
summary(m1)
performance(m1)

# Eta squared
m1 <- aov(P1_amplitude ~ Group*age, data = d_CAMM)
etaSquared(m1, type = 2, anova = FALSE)

by(d_CAMM$P1_amplitude, list(d_CAMM$Group), stat.desc, basic = FALSE)

# Latency
SaturatedModel <- lm(P1_latency ~ Group*age, data = d_CAMM)
m1 <- lm(P1_latency ~ Group+age, data = d_CAMM)
m2 <- lm(P1_latency ~ Group, data = d_CAMM)
anova(SaturatedModel, m1, m2)
AIC(SaturatedModel, m1, m2)
summary(m2)
performance(m2)

by(d_CAMM$P1_latency, list(d_CAMM$Group), stat.desc, basic = FALSE)

# Eta squared
m1 <- aov(P1_latency ~ Group*age, data = d_CAMM)
etaSquared(m1, type = 2, anova = FALSE)

#--- N2
# Amplitude
SaturatedModel <- lm(N2_amplitude ~ Group*age, data = d_CAMM)
m1 <- lm(N2_amplitude ~ Group+age, data = d_CAMM)
m2 <- lm(N2_amplitude ~ Group, data = d_CAMM)
anova(SaturatedModel, m1, m2)
AIC(SaturatedModel, m1, m2)
summary(m1)
performance(m1)

# Eta squared
m1 <- aov(N2_amplitude ~ Group*age, data = d_CAMM)
etaSquared(m1, type = 2, anova = FALSE)

by(d_CAMM$N2_amplitude, list(d_CAMM$Group), stat.desc, basic = FALSE)

# Latency
SaturatedModel <- lm(N2_latency ~ Group*age, data = d_CAMM)
m1 <- lm(N2_latency ~ Group+age, data = d_CAMM)
m2 <- lm(N2_latency ~ Group, data = d_CAMM)
anova(SaturatedModel, m1, m2)
AIC(SaturatedModel, m1, m2)
summary(m2)
performance(m2)

# Eta squared
m1 <- aov(N2_latency ~ Group*age, data = d_CAMM)
etaSquared(m1, type = 2, anova = FALSE)

#--- MMN
# Amplitude
SaturatedModel <- lm(MMN_amplitude ~ Group*age, data = d_CAMM)
m1 <- lm(MMN_amplitude ~ Group+age, data = d_CAMM)
m2 <- lm(MMN_amplitude ~ Group, data = d_CAMM)
anova(SaturatedModel, m1, m2)
AIC(SaturatedModel, m1, m2)
summary(m2)
performance(m2)

# Eta squared
m1 <- aov(MMN_amplitude ~ Group*age, data = d_CAMM)
etaSquared(m1, type = 2, anova = FALSE)

by(d_CAMM$MMN_amplitude, list(d_CAMM$Group), stat.desc, basic = FALSE)

# Latency
SaturatedModel <- lm(MMN_latency ~ Group*age, data = d_CAMM)
m1 <- lm(MMN_latency ~ Group+age, data = d_CAMM)
m2 <- lm(MMN_latency ~ Group, data = d_CAMM)
anova(SaturatedModel, m1, m2)
AIC(SaturatedModel, m1, m2)
summary(m2)
performance(m2)

# Eta squared
m1 <- aov(MMN_latency ~ Group*age, data = d_CAMM)
etaSquared(m1, type = 2, anova = FALSE)