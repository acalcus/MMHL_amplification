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
require(emmeans)
require(readxl)
require(lsr)

# January 2023
# Updated July 2023 - P1 latency + N2 (amplitude & latency)
# clear workspace
rm(list = ls())

# load data file
setwd("~/Dropbox/Dossier de l'équipe InMignonetteWeTrust/UCL/iCARE project/stats/")
d<-read_xlsx("11MMHL_CTX_long_Fz.xlsx")
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
NHU <- filter(d, New_Group == "NH-U")
HLU <- filter(d, New_Group == "HL-U")
HLA <- filter(d, New_Group == "HL-A")

ks.test(NHU$MMN_amplitude, "pnorm", mean = mean(NHU$MMN_amplitude), sd = sd(NHU$MMN_amplitude))
ks.test(NHU$P1_amplitude, "pnorm", mean = mean(NHU$P1_amplitude), sd = sd(NHU$P1_amplitude))
ks.test(NHU$P1_latency, "pnorm", mean = mean(NHU$P1_latency), sd = sd(NHU$P1_latency))
ks.test(NHU$N2_amplitude, "pnorm", mean = mean(NHU$N2_amplitude), sd = sd(NHU$N2_amplitude))
ks.test(NHU$N2_latency, "pnorm", mean = mean(NHU$N2_latency), sd = sd(NHU$N2_latency))

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

# Leven's test
leveneTest(MMN_amplitude ~ New_Group, data = d, center = mean)
leveneTest(P1_amplitude ~ New_Group, data = d, center = mean)
leveneTest(P1_latency ~ New_Group, data = d, center = mean)
leveneTest(N2_amplitude ~ New_Group, data = d, center = mean)
leveneTest(N2_latency ~ New_Group, data = d, center = mean)

############################################################
####### CTX ANALYSES
### Q1: NH-U vs HL-U (between subject comparison)
d_between <- subset(d, Amplification == "unamplified")

#--- MMN
m1 <- lm(MMN_amplitude ~ Group*age, data = d_between)
anova(m1)

by(d_between$MMN_amplitude, list(d_between$Group), stat.desc, basic = FALSE)

p1<-ggplot(data = d_between)+
  geom_point(aes(x = age, y = MMN_amplitude, colour = Group)) + 
  geom_smooth(aes(x = age, y = MMN_amplitude, colour = Group), method = "lm")
p1

# Eta squared
m1 <- aov(MMN_amplitude ~ Group*age, data = d_between)
etaSquared(m1, type = 2, anova = FALSE)

m1 <- lm(MMN_latency ~ Group*age, data = d_between)
anova(m1)

# Eta squared
m1 <- aov(MMN_latency ~ Group*age, data = d_between)
etaSquared(m1, type = 2, anova = FALSE)

#--- P1
lm = lm(P1_amplitude ~ Group*age, data = d_between)
anova(lm)

# Eta squared
m1 <- aov(P1_amplitude ~ Group*age, data = d_between)
etaSquared(m1, type = 2, anova = FALSE)

p1<-ggplot(data = d_between)+
  geom_point(aes(x = age, y = P1_amplitude, colour = Group)) + 
  geom_smooth(aes(x = age, y = P1_amplitude, colour = Group), method = "lm")
p1

by(d_between$P1_amplitude, list(d_between$Group), stat.desc, basic = FALSE)

lm = lm(P1_latency ~ Group*age, data = d_between)
anova(lm)

# Eta squared
m1 <- aov(P1_latency ~ Group*age, data = d_between)
etaSquared(m1, type = 2, anova = FALSE)

p1<-ggplot(data = d_between)+
  geom_point(aes(x = age, y = P1_latency, colour = Group)) + 
  geom_smooth(aes(x = age, y = P1_latency, colour = Group), method = "lm")
p1

by(d_between$P1_latency, list(d_between$Group), stat.desc, basic = FALSE)


#--- N2
lm = lm(N2_amplitude ~ Group*age, data = d_between)
anova(lm)

# Eta squared
m1 <- aov(N2_amplitude ~ Group*age, data = d_between)
etaSquared(m1, type = 2, anova = FALSE)

p1<-ggplot(data = d_between)+
  geom_point(aes(x = age, y = N2_amplitude, colour = Group)) + 
  geom_smooth(aes(x = age, y = N2_amplitude, colour = Group), method = "lm")
p1

by(d_between$N2_amplitude, list(d_between$Group), stat.desc, basic = FALSE)

lm = lm(N2_latency ~ Group*age, data = d_between)
anova(lm)

# Eta squared
m1 <- aov(N2_latency ~ Group*age, data = d_between)
etaSquared(m1, type = 2, anova = FALSE)

p1<-ggplot(data = d_between)+
  geom_point(aes(x = age, y = N2_latency, colour = Group)) + 
  geom_smooth(aes(x = age, y = N2_latency, colour = Group), method = "lm")
p1

by(d_between$N2_latency, list(d_between$Group), stat.desc, basic = FALSE)

#######
### Q2: HL-U vs HL-A (within subject comparison)
d_within <- subset(d, Group == "HL")
d_within <- subset(d_within, code != "HL_020")

# MMN
lm = lmer(MMN_amplitude ~ Amplification*age*BEPTA + (1|code), data = d_within)
anova(lm)
# Model reduction
sm <- step(lm, reduce.random = F)
fm <- get_model(sm)
anova(fm)

# Eta squared
m1 <- aov(MMN_amplitude ~ Amplification*age*BEPTA, data = d_within)
etaSquared(m1, type = 2, anova = FALSE)

MMN<-ggplot(data = d_within)+
  geom_point(aes(x = age, y = MMN_amplitude, colour = Amplification)) + 
  geom_smooth(aes(x = age, y = MMN_amplitude, colour = Amplification), method = "lm")
MMN

# Linear regression on MMN Latency
lm = lmer(MMN_latency ~ Amplification*age*BEPTA + (1|code), data = d_within)
anova(lm)
# Model reduction
sm <- step(lm, reduce.random = F)
fm <- get_model(sm)
anova(fm)

# Eta squared
m1 <- aov(MMN_latency ~ Amplification*age*BEPTA, data = d_within)
etaSquared(m1, type = 2, anova = FALSE)

#--- P1
lm = lmer(P1_amplitude ~ Amplification*age*BEPTA + (1|code), data = d_within)
anova(lm)
# Model reduction
sm <- step(lm, reduce.random = F)
fm <- get_model(sm)
anova(fm)

# Decompose the Amplification x BEPTA interaction
m1_Unampl <- lm(P1_amplitude~BEPTA, data = subset(d_within, New_Group == "HL-U"))
anova(m1_Unampl)
m1_Ampl <- lm(P1_amplitude~BEPTA, data = subset(d_within, New_Group == "HL-A"))
anova(m1_Ampl)

# Eta squared
m1 <- aov(P1_amplitude ~ Amplification*age*BEPTA, data = d_within)
etaSquared(m1, type = 2, anova = FALSE)

p1<-ggplot(data = d_within, aes(x = BEPTA, y = P1_amplitude, colour = Amplification))+
  geom_point() + 
  geom_smooth(method = "lm")
p1

by(d_within$P1_amplitude, list(d_within$Amplification), stat.desc, basic = FALSE)


lm = lmer(P1_latency ~ Amplification*age*BEPTA + (1|code), data = d_within)
anova(lm)
# Model reduction
sm <- step(lm, reduce.random = F)
fm <- get_model(sm)
anova(fm)

# Eta squared
m1 <- aov(P1_latency ~ Amplification*age*BEPTA, data = d_within)
etaSquared(m1, type = 2, anova = FALSE)

p1<-ggplot(data = d_within, aes(x = BEPTA, y = P1_latency, colour = Amplification))+
  geom_point() + 
  geom_smooth(method = "lm")
p1

p1<-ggplot(data = d_within, aes(x = age, y = P1_latency, colour = Amplification))+
  geom_point() + 
  geom_smooth(method = "lm")
p1

#--- N2
lm = lmer(N2_amplitude ~ Amplification*age*BEPTA + (1|code), data = d_within)
anova(lm)
# Model reduction
sm <- step(lm, reduce.random = F)
fm <- get_model(sm)
anova(fm)

p1<-ggplot(data = d_within, aes(x = age, y = N2_amplitude, colour = Amplification))+
  geom_point() + 
  geom_smooth(method = "lm")
p1

# Eta squared
m1 <- aov(N2_amplitude ~ Amplification*age*BEPTA, data = d_within)
etaSquared(m1, type = 2, anova = FALSE)

lm = lmer(N2_latency ~ Amplification*age*BEPTA + (1|code), data = d_within)
anova(lm)
# Model reduction
sm <- step(lm, reduce.random = F)
fm <- get_model(sm)
anova(fm)

# Eta squared
m1 <- aov(N2_latency ~ Amplification*age*BEPTA, data = d_within)
etaSquared(m1, type = 2, anova = FALSE)

#######
### Q3: NH-U vs HL-A (between subject comparison)
# MMN
d_CA <- subset(d, Group == "NH")
d_MM <- subset(d, Amplification == "amplified")
d_CAMM <- rbind(d_CA, d_MM)

# MMN
lm = lm(MMN_amplitude ~ Group*age, data = d_CAMM)
anova(lm)

# Eta squared
m1 <- aov(MMN_amplitude ~ Group*age, data = d_CAMM)
etaSquared(m1, type = 2, anova = FALSE)

# Decompose the Group x Age
Unampl <- lm(MMN_amplitude ~ age, data = filter(d_CAMM, Amplification == "unamplified"))
Ampl <- lm(MMN_amplitude ~ age, data = filter(d_CAMM, Amplification == "amplified"))
anova(Unampl)
anova(Ampl)

p1<-ggplot(data = d_CAMM)+
  geom_point(aes(x = age, y = MMN_amplitude, colour = Group)) + 
  geom_smooth(aes(x = age, y = MMN_amplitude, colour = Group), method = "lm")
p1

# MMN latency
lm = lm(MMN_latency ~ Group*age, data = d_CAMM)
anova(lm)

# Eta squared
m1 <- aov(MMN_latency ~ Group*age, data = d_CAMM)
etaSquared(m1, type = 2, anova = FALSE)

# Decompose the Group x Age
Unampl <- lm(MMN_latency ~ age, data = filter(d_CAMM, Amplification == "unamplified"))
Ampl <- lm(MMN_latency ~ age, data = filter(d_CAMM, Amplification == "amplified"))
anova(Unampl)
anova(Ampl)

p1<-ggplot(data = d_CAMM)+
  geom_point(aes(x = age, y = MMN_latency, colour = Group)) + 
  geom_smooth(aes(x = age, y = MMN_latency, colour = Group), method = "lm")
p1

#--- P1
lm = lm(P1_amplitude ~ Group*age, data = d_CAMM)
anova(lm)
by(d_CAMM$P1_amplitude, list(d_CAMM$Group), stat.desc, basic = FALSE)

# Eta squared
m1 <- aov(P1_amplitude ~ Group*age, data = d_CAMM)
etaSquared(m1, type = 2, anova = FALSE)

lm = lm(P1_latency ~ Group*age, data = d_CAMM)
anova(lm)
by(d_CAMM$P1_latency, list(d_CAMM$Group), stat.desc, basic = FALSE)

# Eta squared
m1 <- aov(P1_latency ~ Group*age, data = d_CAMM)
etaSquared(m1, type = 2, anova = FALSE)

p1<-ggplot(data = d_CAMM)+
  geom_point(aes(x = age, y = P1_latency, colour = Group)) + 
  geom_smooth(aes(x = age, y = P1_latency, colour = Group), method = "lm")
p1

#--- N2
lm = lm(N2_amplitude ~ Group*age, data = d_CAMM)
anova(lm)
by(d_CAMM$N2_amplitude, list(d_CAMM$Group), stat.desc, basic = FALSE)

# Eta squared
m1 <- aov(N2_amplitude ~ Group*age, data = d_CAMM)
etaSquared(m1, type = 2, anova = FALSE)

p1<-ggplot(data = d_CAMM)+
  geom_point(aes(x = age, y = N2_amplitude, colour = Group)) + 
  geom_smooth(aes(x = age, y = N2_amplitude, colour = Group), method = "lm")
p1

lm = lm(N2_latency ~ Group*age, data = d_CAMM)
anova(lm)
by(d_CAMM$N2_latency, list(d_CAMM$Group), stat.desc, basic = FALSE)

# Eta squared
m1 <- aov(N2_latency ~ Group*age, data = d_CAMM)
etaSquared(m1, type = 2, anova = FALSE)

