library(ggplot2)
library(tidyverse)
library(readxl)
library(lme4) 
library(emmeans) 
library(lmerTest)
library(data.table)
library(performance)
library(car)

###whole placenta analysis
lbs <- read.csv("LargeBloodSpacesQuant.csv", header = TRUE)

lbs$gr <- as.factor(paste(lbs$Population, lbs$Treatment))
lbs$uniqueSite <-as.factor(paste(lbs$gr, lbs$DamID, lbs$Slide.ID))
lbs$uniqueDam <-as.factor(paste(lbs$gr, lbs$DamID))

simpleData = lbs %>%
  group_by(Population, Treatment, DamID, Slide.ID, Placenta.ID) %>%
  summarise(WholePlacentaArea = median(WPArea_sqmm),
            LZArea = median(Lzarea_sqmm),
            TotalArea.FBS = median(FBStotal_sqmm),
            TotalPerim.FBS = sum(FBSperimeter_mm),
            TotalBloodSpaces = n())

simpleData$gr <- as.factor(paste(simpleData$Population, simpleData$Treatment))
simpleData$uniqueID <-as.factor(paste(simpleData$gr, simpleData$Slide.ID))
simpleData$uniqueSite <-as.factor(paste(simpleData$gr, simpleData$DamID, simpleData$Slide.ID))
simpleData$uniqueDam <-as.factor(paste(simpleData$gr, simpleData$DamID))

ggplot(simpleData, aes(fill=Population, y=LZArea, x=gr)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

RawAreaLZ <-lmer(LZArea ~ Population*Treatment + (1|uniqueSite),
                 data = simpleData)
check_model(RawAreaLZ)
anova(RawAreaLZ)
#pairs(emmeans(RawAreaLZ, ~Population*Treatment), adjust = "BH")

simpleData$LZRatio <- simpleData$LZArea/simpleData$WholePlacentaArea
ggplot(simpleData, aes(fill=Population, y=LZRatio, x=gr)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

rel <-lmer(LZArea ~ Population*Treatment + WholePlacentaArea + (1|uniqueSite),
             data = simpleData)
check_model(rel)
anova(rel)

####macroscopic vasculature
TotalFBS <-lmer(log(TotalArea.FBS) ~ LZArea + Treatment*Population + (1|uniqueSite),
                data = simpleData)
check_model(TotalFBS)
anova(TotalFBS)

IndivFBS <-lmer(log(Area.FBS) ~ Population*Treatment + (1|uniqueSite),
                data = lbs)
check_model(IndivFBS)
anova(IndivFBS)

TotPerimFBS <-lmer(log(TotalPerim.FBS) ~ log(TotalArea.FBS) + Treatment*Population + (1|uniqueSite),
                data = simpleData)
check_model(TotPerimFBS)
anova(TotPerimFBS)

FBSperim <-lmer(log(FBSperimeter_mm) ~ log(FBSarea_sqmm) + Population*Treatment + (1|uniqueSite),
                data = lbs)
check_model(FBSperim)
anova(FBSperim)

ggplot(lbs, aes(colour=Population, y=FBSperimeter_mm, x=FBSarea_sqmm)) + 
  geom_point(alpha = 0.3) +
  geom_point(alpha = 0.3, shape = 1,colour = "black") + 
  scale_colour_manual(values = c("#5B72D6", "#EDED69"))+
  scale_x_continuous(trans='log2') + 
  scale_y_continuous(trans='log2') +
  geom_smooth(method = "lm") + 
  facet_wrap(~Treatment)

##sterology
d <- read.csv("Microvasculature_StereologyQuant.csv", header = TRUE)

d$gr <- as.factor(paste(d$Population, d$Treatment))
d$uniqueID <- paste(d$DamID, d$SlideID)

ggplot(d, aes(fill=Population, y=Pct_MaternalBloodVolume, x=gr)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.4) +
  ylim(0,76) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

mat <-lmer(Pct_MaternalBloodVolume ~ Population*Treatment + (1|Counter/uniqueID),
           data = d)
check_model(mat)
anova(mat)

ggplot(d, aes(fill=Population, y=Pct_FetalBloodVolume, x=gr)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.4) +
  ylim(0,76) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

fet <-lmer(Pct_FetalBloodVolume ~ Population*Treatment +  (1|Counter/uniqueID),
           data = d)
check_model(fet)
anova(fet)
pairs(emmeans(fet, ~Population))

ggplot(d, aes(fill=Population, y=Pct_TissNuc, x=gr)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.4) +
  ylim(0,76) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tiss <-lmer(Pct_TissNuc ~ Population*Treatment +  (1|Counter/uniqueID),
           data = d)
check_model(tiss)
anova(tiss)
pairs(emmeans(tiss, ~Treatment))

###individual vascular space measurements

imagej <- read.csv("Microvasculature_IndivSpaces.csv", header = TRUE)

imagej$Population[which(imagej$Population=="Lowlander")] <- "1Lowlander"
imagej$Population[which(imagej$Population=="Highlander")] <- "2Highlander"
imagej$Treatment[which(imagej$Treatment=="Hypoxic")] <- "2Hypoxic"
imagej$Treatment[which(imagej$Treatment=="Normoxic")] <- "1Normoxic"

imagej$Blood_Space_Type <- factor(imagej$Blood_Space_Type, levels = c("Maternal", "Fetal"))

imagej$gr <- as.factor(paste(imagej$Population, imagej$Treatment))
imagej$uniqueID <- paste(imagej$DamID, imagej$SlideID)

ggplot(imagej, aes(fill=Population, y=Area_BS, x=gr)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, colour = "black") + 
  scale_fill_manual(values = c("#5B72D6", "#EDED69"))+
  scale_y_continuous(trans='log2') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~Blood_Space_Type)

ggplot(imagej[which(imagej$Blood_Space_Type=="Fetal"),], aes(colour=Population, y=Pm_BS, x=Area_BS)) + 
  geom_point(alpha = 0.3) +
  geom_point(alpha = 0.3, shape = 1,colour = "black") + 
  scale_colour_manual(values = c("#5B72D6", "#EDED69"))+
  scale_x_continuous(trans='log2') + 
  scale_y_continuous(trans='log2') +
  geom_smooth(method = "lm") 

mat <-lmer(log(Area_BS) ~ Population*Treatment +  (1|Counter/uniqueID),
           data = imagej[which(imagej$Blood_Space_Type=="Maternal"),])
check_model(mat)
anova(mat)

matPm <-lmer(log(Pm_BS) ~ Population*Treatment + (1|Counter/uniqueID),
             data = imagej[which(imagej$Blood_Space_Type=="Maternal"),])
check_model(matPm)
anova(matPm)

relmatPm <-lmer(log(Pm_BS) ~ Population*Treatment + log(Area_BS) + (1|Counter/uniqueID),
             data = imagej[which(imagej$Blood_Space_Type=="Maternal"),])
check_model(relmatPm)
anova(relmatPm)

fet <-lmer(log(Area_BS) ~ Population*Treatment +  (1|Counter/uniqueID),
           data = imagej[which(imagej$Blood_Space_Type=="Fetal"),])
check_model(fet)
anova(fet)

fetPm <-lmer(log(Pm_BS) ~ Population*Treatment + (1|Counter/uniqueID),
             data = imagej[which(imagej$Blood_Space_Type=="Fetal"),])
check_model(fetPm)
anova(fetPm)

relfetPm <-lmer(log(Pm_BS) ~ Population*Treatment + log(Area_BS) + (1|Counter/uniqueID),
                data = imagej[which(imagej$Blood_Space_Type=="Fetal"),])
check_model(relfetPm)
anova(relfetPm)