# Count samples
library(tidyr)
library(dplyr)

df<-read.csv("~/Google Drive File Stream/My Drive/Chapter3MicrobiomeMothur/DataCode/MetaDataMicrobeProj.csv")

head(df)



# Remove unknowns
df1<-data.frame(
  df%>%
    filter(TRT!="#N/A")
)


# Remove control
df1<-data.frame(
  df1%>%
    filter(TRT!="Control")
)

df1$Species<-ifelse(grepl("H",df1$Sample_ID),"Ihed","Ip")

df1%>%count(TRT,Species)

head(df1)

df$Combos=as.character(df$Combos)
df[which(df$TRT !="Inter"),]$Combos="none"
df$Combos=as.factor(df$Combos)

setwd("~/Google Drive File Stream/My Drive/Chapter3MicrobiomeMothur/DataSets/")
write.csv(df,"MetaDataMicrobiome.csv",row.names=F)
