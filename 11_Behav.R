# Load packages
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
require(stringr)
require(lsr)
require(car)
require(MASS)
require(lme4)


# February 2023
# clear workspace
rm(list = ls())

# load data file
setwd("~/Dropbox/Dossier de l'équipe InMignonetteWeTrust/UCL/iCARE project/stats/")
d<-read_xls("11MMHL_Behav.xls")
d$Group<-as.factor(d$Group)

###########################
####### Is the data normally distributed ? -> yes
###########################
NH <- filter(d, Group == "NH")
HL <- filter(d, Group == "HL")

ks.test(NH$SpeechDiscrim, "pnorm", mean = mean(NH$SpeechDiscrim, na.rm = TRUE), sd = sd(NH$SpeechDiscrim, na.rm = TRUE))
ks.test(HL$SpeechDiscrim, "pnorm", mean = mean(HL$SpeechDiscrim, na.rm = TRUE), sd = sd(HL$SpeechDiscrim, na.rm = TRUE))

ks.test(NH$Quiet, "pnorm", mean = mean(NH$Quiet, na.rm = TRUE), sd = sd(NH$Quiet, na.rm = TRUE))
ks.test(HL$Quiet, "pnorm", mean = mean(HL$Quiet, na.rm = TRUE), sd = sd(HL$Quiet, na.rm = TRUE))

ks.test(NH$Steady, "pnorm", mean = mean(NH$Steady, na.rm = TRUE), sd = sd(NH$Steady, na.rm = TRUE))
ks.test(HL$Steady, "pnorm", mean = mean(HL$Steady, na.rm = TRUE), sd = sd(HL$Steady, na.rm = TRUE))

ks.test(NH$Fluctuating, "pnorm", mean = mean(NH$Fluctuating, na.rm = TRUE), sd = sd(NH$Fluctuating, na.rm = TRUE))
ks.test(HL$Fluctuating, "pnorm", mean = mean(HL$Fluctuating, na.rm = TRUE), sd = sd(HL$Fluctuating, na.rm = TRUE))


# Leven's test
leveneTest(MMN_amplitude ~ New_Group, data = d, center = mean)
leveneTest(P1_amplitude ~ New_Group, data = d, center = mean)

######### SPEECH DISCRIM
# subselect Speech discrimination task + age
SpeechPerc<-d[,c(1, 2, 4, 11)]
SpeechPerc$Group<-as.factor(SpeechPerc$Group)

# Model
m1 <- lm(SpeechDiscrim ~ Age*Group, data = SpeechPerc, na.action = na.exclude)
anova(m1)

# Eta squared
m1 <- aov(SpeechDiscrim ~ Age*Group, data = SpeechPerc)
etaSquared(m1, type = 2, anova = FALSE)

######### SIN
# subselect SIN tasks + age
SiQ <- d[,c(1, 2, 4, 12, 38)]
SIN <- d[,c(1, 2, 4, 13, 14, 38)]


# Model quiet - GLMM --- THIS IS WRONG, see Stuart's Rmd 
PQL <- glmmPQL(Quiet ~ Age*Group, ~1|Code, family = gaussian(link = "log"), data = SiQ, verbose = FALSE)
summary(PQL)

# Model SIN
SIN_long <- gather(SIN, key = "Condition", value = "SRT", Steady, Fluctuating)
m3 <- lmer(SRT ~ Age*Group*Condition + (1|Code), data = SIN_long)
anova(m3)
# Model reduction
sm <- step(m3, reduce.random = F)
fm <- get_model(sm)
anova(fm)
# Contrasts
em1 <- emmeans(fm, pairwise ~ Condition|Group, adjust = "bonferroni")
print(em1)
em1 <- emmeans(fm, pairwise ~ Group|Condition, adjust = "bonferroni")
print(em1)

by(SIN$Steady, list(SIN$Group), stat.desc, basic = FALSE)

# Eta squared
m4 <- aov(SRT ~ Group*Condition, data = SIN_long)
etaSquared(m4, type = 2, anova = FALSE)