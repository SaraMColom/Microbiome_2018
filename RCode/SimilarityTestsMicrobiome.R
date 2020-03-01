
# Investigating microbiome 'similarity' of compeitors on plant fitness 

# Set directory
setwd("~/Google Drive File Stream/My Drive/Chapter3MicrobiomeMothur/DataSets/MothurOutput_Mock_Microbiome/Microbiome_Output/")

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)



# Read in data
Bray=read.csv("stability.opti_mcc.braycurtis.0.03.column.dist", header = FALSE,sep="\t")
head(Bray)
colnames(Bray)=c("Sample_ID","Species2","BrayDistance")
Bray$Position=gsub("P|H|A","",Bray$Sample_ID)
Bray$Position2=gsub("P|H|A","",Bray$Species2)


jest=read.csv("stability.opti_mcc.jest.0.03.column.dist", header = FALSE,sep="\t")
head(jest)
colnames(jest)=c("Sample_ID","Species2","jestDistance")
jest$Position=gsub("P|H|A","",jest$Sample_ID)
jest$Position2=gsub("P|H|A","",jest$Species2)


Info=read.csv("Sample_Info.csv")
Info$Position=gsub("P|H|A","",Info$Sample_ID)

Bray1=merge(Bray,Info)
head(Bray1)

# Subset by treatment
Alone=Bray1%>% filter(TRT=="Alone")
Comp=Bray1%>% filter(TRT!="Alone")

ggplot(Bray1,aes(TRT,BrayDistance,fill=Species))+
  geom_violin()+
  #geom_jitter(aes(alpha=0.5))+
  theme_classic()+
  facet_grid(~Block)

aggregate(BrayDistance~Species+TRT+Block,Bray1,mean)

hist(Bray1$BrayDistance)

modelBray=lm(BrayDistance~Species+TRT*Block,Bray1)
anova(modelBray)
summary(modelBray)

# Species effect
# Treatment effect
# Block effect
# Treatment by block effect

# Examine average Bray Curtis similarity index on average fitness

  # Subset individuals with matching position
  CompSub=Comp[Comp$Position==Comp$Position2,]
  
  # Read in leaf data (fitness proxy)
  LeafData<-read.csv("~/Google Drive File Stream/My Drive/Field_2018/Field_2018_Rcode/CharacterDisplacementRootTraits/CleanData/Leaf_Number_SizeProxy_Field2018_7232018.csv")[c("Position","ML","Leaf.Number")]
  
  LeafData$Sample_ID=paste(LeafData$Position,ifelse(grepl("Ihed",LeafData$ML),"H","P"),sep="")
  LeafData$Position=gsub("P|H|A","",LeafData$Position)
  
  
  LeafCompdat=droplevels(merge(CompSub,LeafData))
  
ggplot(LeafCompdat,aes(BrayDistance,Leaf.Number,color=Species))+
  geom_point()+
  facet_grid(~Block)
  
LeafCompdat$Leaf.Number=as.numeric(as.character(LeafCompdat$Leaf.Number))
modelLeafBray=glm(Leaf.Number~Block*BrayDistance,LeafCompdat,family = poisson)

anova(modelLeafBray)
summary(modelLeafBray)

# Read in seed number (fitness measurement)

Fitness=read.csv("~/Google Drive File Stream/My Drive/Field_2018/Field_2018_Rcode/CharacterDisplacementRootTraits/CleanData/FitPA4.csv")

# Calculate relative fitness
# First calculate mean seed number by species and treatment---note* we only have seed output of I. purpurea
MeanSeedNumber=Fitness%>%
  group_by(Species,Trt)%>%
  summarise("MeanSeedNumber"=mean(SeedNumber))

Ipurp.Fit=Fitness%>%
  filter(Species=="Ip")


FitnessPurp=merge(Ipurp.Fit,MeanSeedNumber)
FitnessPurp$RelativeFit=FitnessPurp$SeedNumber/FitnessPurp$MeanSeedNumber
FitnessPurp$Block=as.factor(FitnessPurp$Block)

# Average fitness by block
FitnesPurpBlock=aggregate(RelativeFit~Block+Combos,FitnessPurp,mean)


FitnessPurpComp=unique(FitnesPurpBlock[which(FitnesPurpBlock$Combos!=" "),c("Combos","RelativeFit","Block")])


FitCompBray=merge(FitnesPurpBlock,CompSub,by=c("Combos","Block"))
head(FitCompBray)

ggplot(FitCompBray,aes(BrayDistance,RelativeFit))+
  geom_point()+
  facet_grid(~Block)

modelFitBray=lm(RelativeFit~Block*BrayDistance,FitCompBray)

summary(modelFitBray) # No relationship

# Try polynomial
FitCompBray$BrayDistance2=FitCompBray$BrayDistance*FitCompBray$BrayDistance

modelFitBray2=lm(RelativeFit~Block+BrayDistance2+BrayDistance,FitCompBray)
summary(modelFitBray2) # No relationship


# Repeat for Jacaard

jest1=merge(jest,Info)
head(jest1)

# Subset by treatment
Alone=jest1%>% filter(TRT=="Alone")
Comp=jest1%>% filter(TRT!="Alone")

ggplot(jest1,aes(TRT,jestDistance,fill=Species))+
  geom_violin()+
  #geom_jitter(aes(alpha=0.5))+
  theme_classic()+
  facet_grid(~Block)

aggregate(jestDistance~Species+TRT+Block,jest1,mean)

hist(jest1$jestDistance)

modeljest=lm(jestDistance~Species+TRT*Block,jest1)
anova(modeljest)
summary(modeljest)

# Species effect
# Treatment effect
# Block effect
# Treatment by block effect

# Examine average jest Curtis similarity index on average fitness

# Subset individuals with matching position
CompSub=Comp[Comp$Position==Comp$Position2,]

# Read in leaf data (fitness proxy)
LeafData<-read.csv("~/Google Drive File Stream/My Drive/Field_2018/Field_2018_Rcode/CharacterDisplacementRootTraits/CleanData/Leaf_Number_SizeProxy_Field2018_7232018.csv")[c("Position","ML","Leaf.Number")]

LeafData$Sample_ID=paste(LeafData$Position,ifelse(grepl("Ihed",LeafData$ML),"H","P"),sep="")
LeafData$Position=gsub("P|H|A","",LeafData$Position)


LeafCompdat=droplevels(merge(CompSub,LeafData))

ggplot(LeafCompdat,aes(jestDistance,Leaf.Number,color=Species))+
  geom_point()+
  facet_grid(~Block)

LeafCompdat$Leaf.Number=as.numeric(as.character(LeafCompdat$Leaf.Number))
modelLeafjest=glm(Leaf.Number~Block*jestDistance,LeafCompdat,family = poisson)

anova(modelLeafjest)
summary(modelLeafjest)

# Read in seed number (fitness measurement)

Fitness=read.csv("~/Google Drive File Stream/My Drive/Field_2018/Field_2018_Rcode/CharacterDisplacementRootTraits/CleanData/FitPA4.csv")

# Calculate relative fitness
# First calculate mean seed number by species and treatment---note* we only have seed output of I. purpurea
MeanSeedNumber=Fitness%>%
  group_by(Species,Trt)%>%
  summarise("MeanSeedNumber"=mean(SeedNumber))

Ipurp.Fit=Fitness%>%
  filter(Species=="Ip")


FitnessPurp=merge(Ipurp.Fit,MeanSeedNumber)
FitnessPurp$RelativeFit=FitnessPurp$SeedNumber/FitnessPurp$MeanSeedNumber
FitnessPurp$Block=as.factor(FitnessPurp$Block)

# Average fitness by block
FitnesPurpBlock=aggregate(RelativeFit~Block+Combos,FitnessPurp,mean)


FitnessPurpComp=unique(FitnesPurpBlock[which(FitnesPurpBlock$Combos!=" "),c("Combos","RelativeFit","Block")])


FitCompjest=merge(FitnesPurpBlock,CompSub,by=c("Combos","Block"))
head(FitCompjest)

ggplot(FitCompjest,aes(jestDistance,RelativeFit))+
  geom_point()+
  facet_grid(~Block)

modelFitjest=lm(RelativeFit~Block*jestDistance,FitCompjest)

summary(modelFitjest) # No relationship

# Try polynomial
FitCompjest$jestDistance2=FitCompjest$jestDistance*FitCompjest$jestDistance

modelFitjest2=lm(RelativeFit~Block+jestDistance2+jestDistance,FitCompjest)
summary(modelFitjest2) # No relationship
