require(nlme)
require(pastecs)
require(lme4)
require(lsmeans)
require(MuMIn)
require(ggplot2)
require(graphics)
require(readxl)
require(lmerTest)
require(tidyverse)
require(magrittr)
require(varhandle)
require(emmeans)
require(lsr)
require(interactions)

# February 2023
# clear workspace
rm(list = ls())

# load data file
setwd("~/Dropbox/Dossier de l'équipe InMignonetteWeTrust/UCL/iCARE project/stats/")
d<-read.csv("11MMHL_FFRdata_tidy_v1.csv", header = TRUE)

###########################
####### Q1: IS THERE A GROUP EFFECT ON THE UNAMPLIFIED DATA ?
###########################
E_between <- subset(d, Amplification == "unamplified" & Component == "E")
E_between <- subset(E_between, Harmonic == "F0" | Harmonic == "H2")
TFS_between <- subset(d, Amplification == "unamplified" & Component == "TFS")
TFS_between <- subset(TFS_between, Harmonic == "H3" | Harmonic == "H4" | Harmonic == "H5" | Harmonic == "H6")

# ENVELOPE (only F0 and H2, based on noise floor analysis - see 11_FFRanalyses_prelim.R)
m1 <-  lmer(Peak ~ Group*age*Harmonic + (1|code), data = E_between, na.action = na.exclude)
anova(m1)
# Model reduction
sm <- step(m1, reduce.random = F)
fm <- get_model(sm)
anova(fm)

# Eta squared
m1 <- aov(Peak ~ Harmonic*Group, data = E_between)
etaSquared(m1, type = 2, anova = FALSE)

by(E_between$Peak, list(E_between$Group), stat.desc, basic = FALSE)

# TFS (H3 to H6 - see 11_FFRanalyses_prelim.R)
m1 <-  lmer(Peak ~ Group*age*Harmonic + (1|code), data = TFS_between, na.action = na.exclude)
anova(m1)
# Model reduction
sm <- step(m1, reduce.random = F)
fm <- get_model(sm)
anova(fm)

# Eta squared
m1eta <- aov(Peak ~ Harmonic*Group*age, data = TFS_between)
etaSquared(m1eta, type = 2, anova = FALSE)

# Plot
p1 <- ggplot(data = TFS_between, aes(x = Group, y = Peak)) + 
  geom_point(position = position_dodge(0.5), alpha = 0.4) +
  geom_boxplot(lwd = 1) +
  facet_grid(~ Harmonic)
p1

by(TFS_between$Peak, list(TFS_between$Group), stat.desc, basic = FALSE)

###########################
####### Q2: IS THERE A BENEFIT OF AMPLIFICATION IN MM ?
###########################
# Screening: is there an age effect on BEPTA ?
sc <- subset(d, Group == "HL" & Amplification == "unamplified")
fit <- aov(BEPTA~age, data = sc)
summary(fit)

# Plot
p1 <- ggplot(data = sc, aes(x = age, y = BEPTA)) +
  geom_point(size = 3, alpha = 0.4) + 
  geom_smooth(method = "lm", alpha = 0.1) +
  coord_cartesian(xlim = c(8, 16), ylim = c(20, 65))
p1


E_within <- subset(d, Group == "HL" & Component == "E")
E_within <- subset(E_within, Harmonic == "F0" | Harmonic == "H2")
TFS_within <- subset(d, Group == "HL" & Component == "TFS")
TFS_within <- subset(TFS_within, Harmonic == "H3" | Harmonic == "H4" | Harmonic == "H5" | Harmonic == "H6")

# ENVELOPE (only F0 and H2, based on noise floor analysis - see 11_FFRanalyses_prelim.R)
m2 <-  lmer(Peak ~ Amplification*age*Harmonic*BEPTA + (1|code), data = E_within, na.action = na.exclude)
anova(m2)
# Model reduction
sm <- step(m2, reduce.random = F)
fm <- get_model(sm)
anova(fm)

# Eta squared
m2 <- aov(Peak ~ Amplification*age*Harmonic*BEPTA, data = E_within)
etaSquared(m2, type = 2, anova = FALSE)

# Decompose Amplification x age
Unampl <- aov(Peak ~ age, data = subset(E_within, Amplification == "unamplified"))
summary(Unampl)
Ampl <- aov(Peak ~ age, data = subset(E_within, Amplification == "amplified"))
summary(Ampl)

by(E_within$Peak, list(E_within$Amplification), stat.desc, basic = FALSE)

p1 <- ggplot(E_within, aes(x = age, y = Peak, colour = Amplification)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", alpha = 0.1)
p1

# TFS (H3 to H6, based on noise floor analysis - see 11_FFRanalyses_prelim.R)
m3 <-  lmer(Peak ~ Amplification*age*Harmonic*BEPTA + (1|code), data = TFS_within, na.action = na.exclude)
anova(m3)
# Model reduction
sm <- step(m3, reduce.random = F)
fm <- get_model(sm)
anova(fm)

# Eta squared
m3 <- aov(Peak ~ Amplification*age*Harmonic*BEPTA, data = TFS_within)
etaSquared(m3, type = 2, anova = FALSE)

# Decompose the Amplification x Harmonic
Unampl <- lmer(Peak ~ Harmonic + (1|code), data = subset(TFS_within, Amplification == "unamplified"))
anova(Unampl)
Ampl <- lmer(Peak ~ Harmonic + (1|code), data = subset(TFS_within, Amplification == "amplified"))
anova(Ampl)

em1 <- emmeans(fm, pairwise ~ Amplification|Harmonic, adjust = "bonferroni")
print(em1)


############################################################
####### Q3: does amplification normalize the performance ?
####### BETWEEN-SUBJECT COMPARISON
d_MM<-subset(d, New_Group == "HL-A")
d_CA<-subset(d, Group == "NH")
d_CAMM<-rbind(d_MM, d_CA)
d_CAMM_E <- subset(d_CAMM, Component == "E")

# ENVELOPE (only F0 to H2, based on noise floor analysis - see 11_FFRanalyses_prelim.R)
m3 <- lmer(Peak ~ New_Group*age*Harmonic + (1|code), data = d_CAMM_E, na.action = na.exclude)
anova(m3)
# Model reduction
sm <- step(m3, reduce.random = F)
fm <- get_model(sm)
anova(fm)

# Eta squared
m2 <- aov(Peak ~ New_Group*age*Harmonic, data = d_CAMM_E)
etaSquared(m2, type = 2, anova = FALSE)

by(d_CAMM_E$Peak, list(d_CAMM_E$Group), stat.desc, basic = FALSE)

# Plot
p1 <- ggplot(data = d_CAMM_E, aes(x = Group, y = Peak, colour = Group)) + 
  geom_boxplot(position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5), alpha = 0.4) +
  facet_grid(~Harmonic)
p1

# TFS (H3 to H6, based on noise floor analysis - see 11_FFRanalyses_prelim.R)
d_CAMM_TFS <- subset(d_CAMM, Component == "TFS")
d_CAMM_TFS <- subset(d_CAMM_TFS, Harmonic == "H3" | Harmonic == "H4" | Harmonic == "H5" |Harmonic == "H6" )

m3 <- lmer(Peak ~ New_Group*age*Harmonic + (1|code), data = d_CAMM_TFS, na.action = na.exclude)
anova(m3)
# Model reduction
sm <- step(m3, reduce.random = F)
fm <- get_model(sm)
anova(fm)

# Eta squared
m2 <- aov(Peak ~ New_Group*age*Harmonic, data = d_CAMM_TFS)
etaSquared(m2, type = 2, anova = FALSE)

# Contrasts
em1 <- emmeans(fm, pairwise ~ New_Group|Harmonic, adjust = "bonferroni")
print(em1)

# Plot
p1 <- ggplot(data = d_CAMM_TFS, aes(x = Group, y = Peak, colour = Group)) + 
  geom_boxplot(position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5), alpha = 0.4) +
  facet_grid(~Harmonic)
p1

