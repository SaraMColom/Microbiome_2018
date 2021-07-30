RhizMicrobiome\_GenusLevel
================
Sara Colom
2/8/2020

  - [Sample sizes](#sample-sizes)
      - [Table sample number by species and
        treatment](#table-sample-number-by-species-and-treatment)
      - [Table number of maternal line per
        species](#table-number-of-maternal-line-per-species)
  - [Load Libraries](#load-libraries)
  - [Read in Data](#read-in-data)
      - [Sequence depth and pruning](#sequence-depth-and-pruning)
      - [Test for differently abundant
        OTUs](#test-for-differently-abundant-otus)
          - [Visualize the different groups at the family
            level](#visualize-the-different-groups-at-the-family-level)
      - [Testing if different at the phylum
        level](#testing-if-different-at-the-phylum-level)
      - [Family level](#family-level)
      - [Alpha Diversity](#alpha-diversity)
      - [Alpha Diversity](#alpha-diversity-1)
          - [Note: Measures of
            alpha-diversity](#note-measures-of-alpha-diversity)
      - [Test for differences](#test-for-differences)
  - [Alpha Diversity distribution](#alpha-diversity-distribution)
  - [Linear mixed models](#linear-mixed-models)
  - [ANOVA Test for treatment within I. purpurea (Table
    1)](#anova-test-for-treatment-within-i-purpurea-table-1)
  - [Community Composition](#community-composition)
      - [Beta Diversity](#beta-diversity)
  - [PERMANOVA (Table 2)](#permanova-table-2)
  - [Correlations with root traits](#correlations-with-root-traits)
      - [Prep root data](#prep-root-data)
      - [Within species (root traits and
        alphadiv)](#within-species-root-traits-and-alphadiv)
      - [Plotting significant linear
        associations](#plotting-significant-linear-associations)
      - [Linear mixed models (not
        reported)](#linear-mixed-models-not-reported)
      - [Selection on microbe
        variables](#selection-on-microbe-variables)
          - [Selection on richness](#selection-on-richness)
          - [Selection on Inverse
            Simpson](#selection-on-inverse-simpson)
      - [Selection on Simpson](#selection-on-simpson)
      - [Selection on Evenness](#selection-on-evenness)
  - [ANCOVA](#ancova)
  - [MANTEL (Table 4)](#mantel-table-4)

## Sample sizes

### Table sample number by species and treatment

| Species      | Treatment   | N  |
| ------------ | ----------- | -- |
| I. purpurea  | Alone       | 27 |
| I. purpurea  | Competition | 78 |
| I. hederacea | Competition | 78 |

### Table number of maternal line per species

| Species      | Number of ML |
| ------------ | ------------ |
| I. purpurea  | 10           |
| I. hederacea | 5            |

# Load Libraries

``` r
library(phyloseq)
library(ggplot2)
library(ape)
library(vegan)
library(plyr)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(pegas)
library(pgirmess)
library(multcomp)
library(multcompView)
library(ggpubr)
library(ggcorrplot)
library(RColorBrewer)
library(plotly)
library(ggthemes)
library(corrplot)
library(Hmisc)
library(emmeans)
library(lmerTest)


source("miSeq.R")

# Aesthetics
Tx<-theme(axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 20)) +
          theme(axis.text.x = element_text(vjust = 1, hjust=1, angle=0, size = 20),
          axis.title.x = element_text(angle=0, size = 20),
          plot.title=element_text(size = 25,hjust=0))

# Margins
Margin<-theme(
  panel.background = element_rect(fill = "white"),
  plot.margin = margin(2, 2, 2, 2, "cm"),
  plot.background = element_rect(
    fill = "white",
    colour = "black",
    size = 1
  )
)

# Aesthetics
Tx2<-theme(axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25)) +
          theme(axis.text.x = element_text(vjust = 1, hjust=1, size = 25),
          axis.title.x = element_text(size = 25),
          plot.title=element_text(size = 25,hjust=0))

GoldGrey=c("#F1CE63", "#79706E")
GreenBlue=c("#59A14F", "#4E79A7")
```

# Read in Data

``` r
### loading mothur output with FWDB+silva taxonomy and sample metadata. 
### Experiments run in 


sharedfile <- "../DataSets/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"
taxfile <- "../DataSets/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"

mothurdata <- import_mothur(mothur_shared_file = sharedfile, mothur_constaxonomy_file = taxfile)

sampledata <- read.csv('../DataSets/MetaDataTest.csv')

SAMPLE <- sampledata
row.names(SAMPLE)=SAMPLE$Sample_ID

SAMPLE <- subset(SAMPLE, SAMPLE$TRT == "Alone"|SAMPLE$TRT == "Inter")
SAMPLE <- sample_data(SAMPLE)

### create phyloseq object
physeq.all <- merge_phyloseq(mothurdata, SAMPLE) # Modified version worked


### We need to change the taxonomy names: when using the fwdb taxonomy we need to add different headers after removing the last column, which contains no information except for 1 Verrucomicrobia taxon

#   tax_table(physeq.all) <- tax_table(physeq.all)[,-7] # this removes the final column
colnames(tax_table(physeq.all))<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")


#   tax_table(physeq.all) <-cbind(tax_table(physeq.all),row.names(tax_table(physeq.all)))

### removing non-bacterial reads (was already done in mothur, but just to be safe after merging taxonomies)
#   physeq.all <- subset_taxa(physeq.all, Kingdom == "Bacteria")
```

``` r
## Aggregate at the genus level, then saving ea componant externally to reload and save time.

        # physeq_genus <- physeq.all %>%
         #tax_glom(taxrank = "Genus") 

        # otu=data.frame(otu_table(physeq_genus))
        
        # meta=data.frame(sample_data(physeq_genus))
        # tax=data.frame(tax_table(physeq_genus))
        
        # write.csv(x = otu,file = "otu_table.csv",row.names = TRUE)
        # write.csv(x = tax,file = "tax_table.csv")# 
        # write.csv(x = meta,file = "sample_table.csv",row.names = F)


otu <- read.csv("../DataSets/otu_table.csv", row.names = 1)

colnames(otu) <- gsub('X', "", colnames(otu))

otu <- otu_table(otu, taxa_are_rows = TRUE)


tax <- read.csv('../DataSets/tax_table.csv', row.names = 1)
taxRows <- row.names(tax)
taxCols <- colnames(tax)
tax <- tax_table(as.matrix(tax))

#row.names(tax)=taxRows
#colnames(tax)=taxCols


meta <- read.csv('../DataSets/MetaDataTest.csv')
namesKeep <- colnames(otu)

row.names(meta) <- meta$Sample_ID

meta <- meta[which(meta$Sample_ID %in% namesKeep), ]

meta <- sample_data(meta)
   
physeq_genus <- merge_phyloseq(otu, tax, meta)
```

``` r
### ADD THE PROTEOBACTERIA CLASSES TO THE PHYLA NAME FIELD IN PHYLOSEQ OBJECT TAXONOMY 

phy <- data.frame(tax_table(physeq_genus))
Phylum <- as.character(phy$Phylum)
Class <- as.character(phy$Class)
for  (i in 1:length(Phylum)){ 
  if (Phylum[i] == "Proteobacteria"){
    if (Class[i] == "unclassified"){
      Phylum[i] <- Phylum[i]       
    } else {
      Phylum[i] <- Class[i]
    }
  } 
}


Phylum <- as.character(phy$Phylum)
Class <- as.character(phy$Class)
phy$Phylum <- Phylum
t <- tax_table(as.matrix(phy))
```

## Sequence depth and pruning

``` r
physeq <- (physeq.all2)
physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq) # Remove taxa with no counts

#check number of reads in each sample, differences in count are in part due differet numbers of chlorophyl reads depending on time of experiment

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(physeq))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "#59A14F", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) +
  theme_classic()
```

![](README_files/figure-gfm/data%20scaling-1.png)<!-- -->

``` r
# Scales reads to smallest library size 
source("https://raw.githubusercontent.com/michberr/MicrobeMiseq/master/R/miseqR.R")
#physeq.scale <- scale_reads(physeq, min(sample_sums(physeq)))

##### Normalization #######

# Scales reads by 
# 1) taking proportions,
# 2) multiplying by a given library size of n
# 3) rounding down

physeq1 <- prune_samples(sample_sums(physeq) > 20000, physeq)

sample_sum_df1 <- data.frame(sum = sample_sums(physeq1))

# Histogram of sample read counts
ggplot(sample_sum_df1, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "#59A14F", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth after pruning") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) +
  theme_classic()
```

![](README_files/figure-gfm/data%20scaling-2.png)<!-- -->

``` r
n <- min(sample_sums(physeq1))
  physeq.scale <-
    transform_sample_counts(physeq1, function(x) {
      (n * x/sum(x))  # Transform to rel. abundance
    })
  
physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)


physeq1
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1098 taxa and 173 samples ]
    ## sample_data() Sample Data:       [ 173 samples by 9 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1098 taxa by 6 taxonomic ranks ]

``` r
physeq.bray <- phyloseq::distance(physeq = physeq.scale, method = "bray")

# # # # # # # # # # # # # # # # # # # # # 
# Subsample for within I. purpurea only
# # # # # # # # # # # # # # # # # # # # # 

physeq.Purp <- subset_samples(physeq.scale, Species == "Ip")
sampledf.Purp<- data.frame(sample_data(physeq.Purp))

# Calculate bray curtis for summer samples only
physeq.Purp.bray <- phyloseq::distance(physeq = physeq.Purp, method = "bray")
```

## Test for differently abundant OTUs

### Visualize the different groups at the family level

``` r
# Check for 'core' microbiome members as the phylum level which taxa are present
# at 1 % overall abundance and at least 75% of samples 

PhylumGlom <- tax_glom(physeq.Purp,taxrank = "Phylum")
FamGlom <- tax_glom(physeq.Purp,taxrank = "Family")
ClassGlom <- tax_glom(physeq.Purp,taxrank = "Class")


coreTaxa <- filter_taxa(PhylumGlom, function(x) sum(x > 1) > (0.75 * length(x)), TRUE)

coreTaxaFamily <- filter_taxa(FamGlom, function(x) sum(x > 1) > (0.75 * length(x)), TRUE)

## Plot at phylum level


physeq1_phylumAlone <- coreTaxa  %>%
  subset_samples(TRT == "Alone")%>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()      


physeq1_phylumComp <- coreTaxa  %>%
  subset_samples(TRT != "Alone")%>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()       

physeq1_phylum <- rbind(physeq1_phylumAlone,physeq1_phylumComp)

# Melt to long format
physeq1_phylum <- physeq1_phylum[order(physeq1_phylum$Phylum), ] 

#colnames(physeq1_phylum)[13]="Phylum"

# Sort data frame alphabetically by phylum
# Set colors for plotting
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange", "#DA5724", "#508578", "#CD9BCD",
   "#AD6F3B", "#673770", "#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"
)


#colnames(physeq1_phylum)[which(names(physeq1_phylum)%in%"Rank2")]="Phylum"

# Plot 
ggplot(physeq1_phylum, aes(x = TRT, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity",position='dodge') +
  #scale_fill_manual(values = phylum_colors) +
  theme(axis.title.x = element_blank()) + 
  ylab("Relative Abundance (Phyla > 75 %) \n") +
  ggtitle("Phylum Composition  \n Bacterial Communities by Sampling per Treatment") +
  theme_classic() +
  ylab("Relative Abundance")
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## Testing if different at the phylum level

``` r
PhlyaList <- as.list(unique(physeq1_phylum$Phylum))

KrusMicrTest <- function(Microbe){
  Result = kruskal.test(Abundance ~ TRT, data = physeq1_phylum[which(physeq1_phylum$Phylum == Microbe),]) 

ResDf <- data.frame("ChiSq" = Result$statistic, "Pval"=Result$p.value, "Phylum" = paste(Microbe))  
  return(ResDf)  
}

AllResults <- lapply(PhlyaList, KrusMicrTest)
KruskalRes <- do.call('rbind', AllResults)
```

## Family level

``` r
# Family level
physeq1_FamilyAlone <- coreTaxaFamily  %>%
  subset_samples(TRT == "Alone")%>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()                                      # Melt to long format

physeq1_FamilyComp <- coreTaxaFamily  %>%
  subset_samples(TRT != "Alone")%>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()       

physeq1_Family <- rbind(physeq1_FamilyAlone,physeq1_FamilyComp)


FamList <- as.list(unique(physeq1_Family$Family))

KrusMicrTest<-function(Microbe){
  Result <- kruskal.test(Abundance ~ TRT, data = physeq1_Family[which(physeq1_Family$Family == Microbe), ]) 

ResDf <- data.frame("ChiSq" = Result$statistic, "Pval"= Result$p.value, "Family" = paste(Microbe))  
  return(ResDf)  
}

AllResultsFamily <- lapply(FamList, KrusMicrTest)
KruskalResFamily <- do.call('rbind', AllResultsFamily)

KruskalResFamily %>% 
  filter(Pval < 0.055)
```

    ##                                  ChiSq        Pval          Family
    ## Kruskal-Wallis chi-squared54  5.192655 0.022682540   Subgroup_2_fa
    ## Kruskal-Wallis chi-squared69  3.904351 0.048161234     Gaiellaceae
    ## Kruskal-Wallis chi-squared139 3.752445 0.052730323  Lineage_IIa_fa
    ## Kruskal-Wallis chi-squared177 4.186783 0.040740338          SM2D12
    ## Kruskal-Wallis chi-squared180 4.509080 0.033715377        AKYG1722
    ## Kruskal-Wallis chi-squared181 3.904351 0.048161234 Microtrichaceae
    ## Kruskal-Wallis chi-squared207 8.682062 0.003213572           A0839

``` r
# Do false discovery rate for mulitple corrections

KruskalResFamily$AdjPval <- p.adjust(KruskalResFamily$Pval, method = "fdr", n = 213)
```

## Alpha Diversity

``` r
hist(sample_sums(physeq1))
```

![](README_files/figure-gfm/alpha%20diversity-1.png)<!-- -->

``` r
# Remove sample with less than 20K reads, this looks WAAAY off
# Rarify first
ps.rarefied <- rarefy_even_depth(physeq1, rngseed = 1, sample.size = min(sample_sums(physeq1)), replace = F)

plot_richness(ps.rarefied,x = "TRT",measures = c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic()
```

![](README_files/figure-gfm/alpha%20diversity-2.png)<!-- -->

## Alpha Diversity

#### Note: Measures of alpha-diversity

*_Inverse Simpson_* it is an indication of the richness in a community
with uniform evenness that would have the same level of diversity. So
while measures such as the Shannon index are somewhat abstract, the
inverse of the Simpson index has some biological interpretation. Other
advantages of the Simpson-based metrics are that they do not tend to be
as affected by sampling effort as the Shannon index.

  - *_Species richness_* is simply a count of species, and it does not
    take into account the abundances of the species or their relative
    abundance distributions.

  - *_Simpson’s Diversity Index_* is a measure of diversity which takes
    into account the number of species present, as well as the relative
    abundance of each species. As species richness and evenness
    increase, so diversity increases.

<!-- end list -->

``` r
### I commented out the subsampling to save time. The matrices produced were saved and then reoponed from their tab delmited formats into corresponding matrices.

set.seed(711)

r <- rarefy_even_depth(physeq1, sample.size = n, verbose = FALSE, replace = TRUE, rngseed = 711)
         
 ## Calculate richness
rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))

## Calculate Inverse Simpson

simp <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))

sim<- as.numeric(as.matrix(estimate_richness(r, measures = "Simpson")))

s <- rarefy_even_depth(physeq1, sample.size = n, verbose = FALSE, replace = TRUE, rngseed = 711)

shan <- as.numeric(as.matrix(estimate_richness(r, measures = "Shannon")))

InvSimp <- simp
```

``` r
# Create a new dataframe to hold the means and standard deviations of richness estimates

length(rich)
```

    ## [1] 173

``` r
Sample_ID <- sample_names(physeq1)
Block <- sample_data(physeq1)$Block
Species <-sample_data(physeq1)$Species
TRT <- sample_data(physeq1)$TRT
Combos <-sample_data(physeq1)$Combos
ML <- sample_data(physeq1)$ML

alpha <- data.frame(Sample_ID,ML, Block, TRT, Species, Combos, rich, InvSimp, sim, shan)

alpha$even <- alpha$shan/alpha$rich

# DO VIOLIN PLOT HERE!!!


p <- ggplot(alpha %>% filter(TRT == "Inter"), aes(x = Species, y = rich)) +
  geom_violin(trim = FALSE, aes(fill = Species), alpha = 0.3) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1, aes(color = Species, fill = Species)) +
  theme_classic() +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  ggtitle("Mean Species Richness")


q <- ggplot(alpha %>% filter(TRT == "Inter"), aes(x = Species, y = InvSimp)) +
  geom_violin(trim = FALSE, aes(fill = Species), alpha = 0.3) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1, aes(color = Species, fill = Species)) +
  theme_classic() +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  ggtitle("Mean Species Inverse Simpson") +
  ylab("")


t <- ggplot(alpha %>% filter(TRT == "Inter"), aes(x = Species, y = sim)) +
  geom_violin(trim = FALSE, aes(fill=Species), alpha = 0.3) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1, aes(color = Species, fill = Species)) +
  theme_classic() +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  ggtitle("Mean Simpson") +
  ylab("")

v <- ggplot(alpha %>% filter(TRT == "Inter"), aes(x = Species, y = even)) +
  geom_violin(trim = FALSE, aes(fill = Species), alpha = 0.3) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1, aes(color = Species, fill = Species)) +
  theme_classic() +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  ggtitle("Mean Evenness") +
  ylab("")



ggarrange(p, q, t, v, common.legend = T, ncol = 2, nrow = 2)
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
      ### Treatment

pT <- ggplot(alpha %>% filter(Species == "Ip"), aes(x = TRT, y = rich)) +
  geom_violin(trim = FALSE, aes(fill = TRT), alpha = 0.3) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1, aes(color = TRT,fill = TRT)) +
  theme_classic() +
  scale_color_manual(values = c("#B2DF8A", "#fee090")) +
  scale_fill_manual(values = c("#B2DF8A", "#fee090")) +
  ggtitle("Mean Treatment Richness")


qT <- ggplot(alpha %>% filter(Species == "Ip"), aes(x = TRT, y = InvSimp)) +
  geom_violin(trim = FALSE, aes(fill = TRT), alpha = 0.3) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize=1, aes(color = TRT, fill = TRT)) +
  theme_classic() +
  scale_color_manual(values = c("#B2DF8A", "#fee090")) +
  scale_fill_manual(values = c("#B2DF8A", "#fee090")) +
  ggtitle("Mean Treatment Inverse Simpson") +
  ylab(" ")

ggarrange(pT, qT, common.legend = T)
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

![](README_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

## Test for differences

# Alpha Diversity distribution

``` r
# first check for normality. To test for normalcy statistically, we can run the Shapiro-Wilk test of normality.

shapiro.test(alpha$rich) # Not normal
```

    FALSE 
    FALSE   Shapiro-Wilk normality test
    FALSE 
    FALSE data:  alpha$rich
    FALSE W = 0.98984, p-value = 0.2539

``` r
shapiro.test(alpha$even) # Normal
```

    FALSE 
    FALSE   Shapiro-Wilk normality test
    FALSE 
    FALSE data:  alpha$even
    FALSE W = 0.98648, p-value = 0.09413

``` r
shapiro.test(alpha$sim) # Normal
```

    FALSE 
    FALSE   Shapiro-Wilk normality test
    FALSE 
    FALSE data:  alpha$sim
    FALSE W = 0.90115, p-value = 2.331e-09

``` r
histogram(alpha$rich,xlab = "Richness") # Note richness though, not normal, looks close to normal.
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
histogram(alpha$InvSimp, xlab = "Inverse-Simpson")
```

![](README_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
histogram(alpha$sim, xlab = "Simpson")
```

![](README_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

# Linear mixed models

``` r
# Test for treatment within I. purpurea


RichLmm <- lmer(rich ~ TRT + Block + (1|Block:ML), alpha %>% filter(Species == "Ip"))

InvLmm <- lmer(InvSimp ~ TRT + Block + (1|Block:ML), alpha %>% filter(Species == "Ip"))

SimLmm <- lmer(sim ~ TRT + Block + (1|Block:ML), alpha %>% filter(Species == "Ip"))

EvenLmm <- lmer(even ~ TRT + Block + (1|Block:ML), alpha %>% filter(Species == "Ip"))
```

    ## boundary (singular) fit: see ?isSingular

``` r
anova(RichLmm)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##        Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
    ## TRT    1320.5  1320.5     1 87.801  1.4442 0.232690   
    ## Block 13858.8  4619.6     3 26.431  5.0524 0.006744 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
ranova(RichLmm)
```

    ## ANOVA-like table for random-effects: Single term deletions
    ## 
    ## Model:
    ## rich ~ TRT + Block + (1 | Block:ML)
    ##                npar  logLik    AIC    LRT Df Pr(>Chisq)
    ## <none>            7 -467.21 948.42                     
    ## (1 | Block:ML)    6 -467.22 946.44 0.0238  1     0.8774

``` r
anova(InvLmm)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##       Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)  
    ## TRT     0.27   0.273     1 82.613  0.0035 0.9533  
    ## Block 721.76 240.588     3 24.544  3.0419 0.0479 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
ranova(InvLmm)
```

    ## ANOVA-like table for random-effects: Single term deletions
    ## 
    ## Model:
    ## InvSimp ~ TRT + Block + (1 | Block:ML)
    ##                npar  logLik    AIC    LRT Df Pr(>Chisq)
    ## <none>            7 -354.94 723.88                     
    ## (1 | Block:ML)    6 -355.46 722.92 1.0486  1     0.3058

``` r
anova(SimLmm)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq    Mean Sq NumDF  DenDF F value Pr(>F)
    ## TRT   0.00002357 2.3571e-05     1 79.747  0.4399 0.5091
    ## Block 0.00033833 1.1278e-04     3 22.261  2.1049 0.1283

``` r
ranova(SimLmm)
```

    ## ANOVA-like table for random-effects: Single term deletions
    ## 
    ## Model:
    ## sim ~ TRT + Block + (1 | Block:ML)
    ##                npar logLik     AIC    LRT Df Pr(>Chisq)
    ## <none>            7 318.62 -623.24                     
    ## (1 | Block:ML)    6 317.94 -623.89 1.3513  1      0.245

``` r
anova(EvenLmm)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq    Mean Sq NumDF DenDF F value   Pr(>F)   
    ## TRT   2.2628e-07 2.2628e-07     1    95  1.1782 0.280475   
    ## Block 2.9697e-06 9.8990e-07     3    95  5.1542 0.002408 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
ranova(EvenLmm)
```

    ## ANOVA-like table for random-effects: Single term deletions
    ## 
    ## Model:
    ## even ~ TRT + Block + (1 | Block:ML)
    ##                npar logLik     AIC LRT Df Pr(>Chisq)
    ## <none>            7 591.89 -1169.8                  
    ## (1 | Block:ML)    6 591.89 -1171.8   0  1          1

# ANOVA Test for treatment within I. purpurea (Table 1)

``` r
alpha.purp <- alpha %>% filter(Species == "Ip")

aov.Richness <- lm(rich ~ TRT + Block, alpha.purp)
aov.simpsonInv <- lm(InvSimp ~ TRT + Block, alpha.purp)
aov.simpson <- lm(sim ~ TRT + Block, alpha.purp)
aov.evenness <- lm(even ~ TRT + Block, alpha.purp)
aov.shannon <- lm(shan ~ TRT + Block, alpha.purp)


#Call for the summary of that ANOVA, which will include P-values
anova(aov.Richness)
```

    FALSE Analysis of Variance Table
    FALSE 
    FALSE Response: rich
    FALSE           Df Sum Sq Mean Sq F value   Pr(>F)   
    FALSE TRT        1    962   961.9  1.0381 0.310860   
    FALSE Block      3  14694  4897.8  5.2858 0.002052 **
    FALSE Residuals 95  88028   926.6                    
    FALSE ---
    FALSE Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(aov.simpsonInv)
```

    FALSE Analysis of Variance Table
    FALSE 
    FALSE Response: InvSimp
    FALSE           Df Sum Sq Mean Sq F value  Pr(>F)  
    FALSE TRT        1    2.5    2.54  0.0288 0.86562  
    FALSE Block      3 1031.6  343.85  3.9020 0.01122 *
    FALSE Residuals 95 8371.5   88.12                  
    FALSE ---
    FALSE Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(aov.simpson)
```

    FALSE Analysis of Variance Table
    FALSE 
    FALSE Response: sim
    FALSE           Df    Sum Sq    Mean Sq F value  Pr(>F)  
    FALSE TRT        1 0.0000171 1.7100e-05  0.2785 0.59889  
    FALSE Block      3 0.0005091 1.6970e-04  2.7643 0.04617 *
    FALSE Residuals 95 0.0058321 6.1391e-05                  
    FALSE ---
    FALSE Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(aov.evenness)
```

    FALSE Analysis of Variance Table
    FALSE 
    FALSE Response: even
    FALSE           Df     Sum Sq    Mean Sq F value   Pr(>F)   
    FALSE TRT        1 1.9700e-07 1.9704e-07  1.0259 0.313685   
    FALSE Block      3 2.9697e-06 9.8990e-07  5.1542 0.002408 **
    FALSE Residuals 95 1.8245e-05 1.9206e-07                    
    FALSE ---
    FALSE Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# No significant differences
```

# Community Composition

## Beta Diversity

``` r
# Beta diversity pcoa Bray-Curtis DNA only 
physeq.pcoa <-
  ordinate(
    physeq = physeq.scale,
    method = "PCoA",
    distance = "bray"
  )

physeq.pcoa.vectors <- data.frame(physeq.pcoa$vectors[, 1:4])

physeq.pcoa.vectors$Duplicates <- row.names(physeq.pcoa.vectors)

SampData <- data.frame(sample_data(physeq))

colnames(SampData)[1] <- "Duplicates"

SampData <- subset(SampData, SampData$TRT == "Inter"|SampData$TRT == "Alone")

physeq.pcoa.df <- droplevels(merge(physeq.pcoa.vectors, SampData,by="Duplicates"))

bray_values <- physeq.pcoa$values
bray_rel_eigens <- bray_values$Relative_eig

bray_rel_eigen1 <- bray_rel_eigens[1]
bray_rel_eigen1_percent <- round(bray_rel_eigen1 * 100, digits = 1)

bray_rel_eigen2 <- bray_rel_eigens[2]
bray_rel_eigen2_percent <- round(bray_rel_eigen2 * 100, digits = 1)

bray_rel_eigen3 <- bray_rel_eigens[3]
bray_rel_eigen3_percent <- round(bray_rel_eigen3 * 100, digits = 1)

bray_rel_eigen4 <- bray_rel_eigens[4]
bray_rel_eigen4_percent <- round(bray_rel_eigen4 * 100, digits = 1)

bray_rel_eigen5 <- bray_rel_eigens[5]
bray_rel_eigen5_percent <- round(bray_rel_eigen5 * 100, digits = 1)

# Prep axis labels & title
bray_axis1 <- paste("PCoA 1:",bray_rel_eigen1_percent, "%")
bray_axis2 <- paste("PCoA 2:",bray_rel_eigen2_percent, "%")
bray_axis3 <- paste("PCoA 3:",bray_rel_eigen3_percent, "%")
bray_axis4 <- paste("PCoA 4:",bray_rel_eigen4_percent, "%")
PCoA_title <- paste("Bray-Curtis, ",ntaxa(physeq.scale), "OTUs")


pcoa_exp_trt <- ggplot(physeq.pcoa.df, aes(Axis.1, Axis.2, color = TRT, fill = TRT)) +
  xlab(bray_axis1) + 
  ylab(bray_axis2) +
  geom_point(alpha = 0.9) + 
  theme_classic() +
  scale_fill_manual(values = c("#00B050", "grey"), "Treatment", labels = c("Alone", "Competition")) +
  scale_color_manual(values = c("#00B050", "grey"), "Treatment", labels = c("Alone", "Competition")) +
  stat_ellipse(geom = "polygon", alpha = 1/6, aes(fill = TRT)) +
  Tx +
  theme(axis.text.x = element_text(angle = 45), axis.text = element_text(color="black", size = 15,vjust = 0.5,hjust = 1))
#  scale_fill_discrete(name = "Treatment", labels = c("Alone", "Competition")) +
 # scale_color_discrete(name = "Treatment", labels = c("Alone", "Competition"))



pcoa_exp_trtA <- ggplot(physeq.pcoa.df, aes(Axis.3, Axis.4, color = TRT,fill = TRT)) +
  xlab(bray_axis3) + 
  ylab(bray_axis4) + 
  geom_point(alpha = 0.9) + 
  theme_classic() +
  scale_fill_manual(values = GoldGrey, "Treatment", labels = c("Alone", "Competition")) +
  scale_color_manual(values = GoldGrey, "Treatment", labels = c("Alone", "Competition")) +
  stat_ellipse(geom = "polygon", alpha = 1/6, aes(fill = TRT)) +
  Tx


pcoa_exp_sp <- ggplot(physeq.pcoa.df, aes(Axis.1, Axis.2, color = Species,fill = Species)) +
  xlab("") + 
  ylab(bray_axis2) +
  geom_point(alpha = 0.9) + 
  theme_classic() +
 scale_fill_manual(values = GreenBlue, "Species", labels = c("I. hederacea", "I. purpurea")) +
  scale_color_manual(values = GreenBlue, "Species", labels = c("I. hederacea", "I. purpurea")) +
   stat_ellipse(geom = "polygon", alpha = 1/6, aes(fill = Species)) +
  Tx


pcoa_exp_spA <- ggplot(physeq.pcoa.df, aes(Axis.3, Axis.4, color = Species,fill = Species)) +
  xlab("") + 
  ylab(bray_axis4) +
  geom_point(alpha = 0.9) + 
  theme_classic() +
  scale_fill_manual(values = GreenBlue, "Species", labels = c("I. hederacea", "I. purpurea")) +
  scale_color_manual(values = GreenBlue, "Species", labels = c("I. hederacea", "I. purpurea")) +
  stat_ellipse(geom = "polygon", alpha = 1/6, aes(fill = Species)) +
  Tx

A <- ggarrange(pcoa_exp_trt,pcoa_exp_trtA,common.legend = T, Labels = c("C", "D"),font.label = list(size = 15, color = "black", face =
  "plain"),hjust = -8,vjust = 1)

B <- ggarrange(pcoa_exp_sp,pcoa_exp_spA,common.legend = T, labels = c("A", "B"),font.label = list(size = 15, color = "black", face =
  "plain"),hjust = -8,vjust = 1)

MainFig <- ggarrange(B, A, nrow = 2)

annotate_figure(MainFig,
                top = text_grob(PCoA_title, face = "plain", size = 25, lineheight = 2)
)
```

![](README_files/figure-gfm/PCoA%20Treatment%20and%20Species%20comparisoin%20Bray%20Curtis%20Estimate-1.png)<!-- -->

``` r
physeq.pcoa.df$TrtSp <- paste(physeq.pcoa.df$TRT, physeq.pcoa.df$Species, sep="")

p <- plot_ly(physeq.pcoa.df, x = ~Axis.1, y = ~Axis.2, z = ~Axis.3, color = ~TrtSp) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = bray_axis1),
                     yaxis = list(title = bray_axis2),
                     zaxis = list(title = bray_axis3)))


p1 <- plot_ly(physeq.pcoa.df, x = ~Axis.1, y = ~Axis.2, z = ~Axis.4, color = ~Species) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = bray_axis1),
                     yaxis = list(title = bray_axis2),
                     zaxis = list(title = bray_axis4)))

p3 <- plot_ly(physeq.pcoa.df, x = ~Axis.1, y = ~Axis.2, z = ~Axis.4, color = ~Species) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = bray_axis1),
                     yaxis = list(title = bray_axis3),
                     zaxis = list(title = bray_axis4)))
```

# PERMANOVA (Table 2)

``` r
## All experiment, add variable treatment, then species
sampledf <- data.frame(sample_data(physeq.scale))

physeq.purp = subset_samples(physeq.scale, Species == "Ip")

treatment <- sampledf %>% 
  filter(Species == "Ip") %>% 
  pull(TRT)

block <- sampledf %>% 
  filter(Species == "Ip") %>% 
  pull(Block)

adonis(otu_table(physeq.purp) %>% t ~ treatment + block, method = "bray")
```

    FALSE 
    FALSE Call:
    FALSE adonis(formula = otu_table(physeq.purp) %>% t ~ treatment + block,      method = "bray") 
    FALSE 
    FALSE Permutation: free
    FALSE Number of permutations: 999
    FALSE 
    FALSE Terms added sequentially (first to last)
    FALSE 
    FALSE           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
    FALSE treatment  1    0.0258 0.025835  0.8460 0.00796  0.551    
    FALSE block      3    0.3189 0.106312  3.4813 0.09826  0.001 ***
    FALSE Residuals 95    2.9011 0.030538         0.89378           
    FALSE Total     99    3.2459                  1.00000           
    FALSE ---
    FALSE Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Correlations with root traits

## Prep root data

``` r
RootData <- read.csv("../DataSets/RootTraits_PCs.csv")

RootAlphaObs <- merge(alpha, RootData[c("Sample_ID", "PC1", "PC2", "PC3", "PC4")])
```

## Within species (root traits and alphadiv)

``` r
library(multcomp)

#summary(glht(SimpInvPC1, mcp(rank="Tukey")))

#################################################################
##################### SUBSET for I.purpurea #####################
#################################################################


RootAlphaPurp <- droplevels(RootAlphaObs %>% filter(Species == "Ip"))
RootAlphaPurp$Comp <- sub(".*\\-", "", RootAlphaPurp$Combos)

### Linear regressions

SimpPC1 <- lm(sim ~ PC1 + Block + TRT, RootAlphaPurp) 
summary(SimpPC1)
```

    ## 
    ## Call:
    ## lm(formula = sim ~ PC1 + Block + TRT, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.027453 -0.003862  0.001500  0.005865  0.011300 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.977989   0.003869 252.755   <2e-16 ***
    ## PC1          0.001447   0.001258   1.150   0.2554    
    ## Block2      -0.005401   0.010369  -0.521   0.6047    
    ## Block3      -0.009643   0.010882  -0.886   0.3797    
    ## Block4      -0.002049   0.002635  -0.777   0.4405    
    ## TRTInter    -0.004216   0.002470  -1.707   0.0939 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.008187 on 51 degrees of freedom
    ## Multiple R-squared:  0.1337, Adjusted R-squared:  0.04877 
    ## F-statistic: 1.574 on 5 and 51 DF,  p-value: 0.1841

``` r
SimpInvPC1 <- lm(InvSimp ~ PC1 + Comp + TRT + Block, RootAlphaPurp) 
summary(SimpInvPC1)
```

    ## 
    ## Call:
    ## lm(formula = InvSimp ~ PC1 + Comp + TRT + Block, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -14.0281  -6.9591   0.3038   6.3289  18.8379 
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       42.7370     4.3928   9.729 6.18e-13 ***
    ## PC1                1.0167     1.4322   0.710   0.4812    
    ## CompPA 4.12 Ihed  -2.5142     3.6273  -0.693   0.4916    
    ## CompPA 4.15 Ihed  -0.2098     3.6612  -0.057   0.9545    
    ## CompPA 4.2 Ihed   -9.5479     4.1700  -2.290   0.0265 *  
    ## CompPA 4.3 Ihed   -3.6759     4.0327  -0.912   0.3666    
    ## TRTInter               NA         NA      NA       NA    
    ## Block2            -1.2232    11.8157  -0.104   0.9180    
    ## Block3            -5.0956    12.3933  -0.411   0.6828    
    ## Block4            -2.9551     2.9900  -0.988   0.3280    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 9.233 on 48 degrees of freedom
    ## Multiple R-squared:  0.2126, Adjusted R-squared:  0.08133 
    ## F-statistic:  1.62 on 8 and 48 DF,  p-value: 0.1441

``` r
RichPC1 <- lm(rich ~ PC1 + Block + TRT, RootAlphaPurp) 
summary(RichPC1)
```

    ## 
    ## Call:
    ## lm(formula = rich ~ PC1 + Block + TRT, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -68.694 -17.875   1.742  20.104  51.312 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  461.869     14.307  32.283   <2e-16 ***
    ## PC1            3.379      4.650   0.727   0.4707    
    ## Block2        -7.192     38.340  -0.188   0.8519    
    ## Block3        -2.043     40.236  -0.051   0.9597    
    ## Block4         1.186      9.744   0.122   0.9036    
    ## TRTInter     -17.545      9.133  -1.921   0.0603 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 30.27 on 51 degrees of freedom
    ## Multiple R-squared:  0.1664, Adjusted R-squared:  0.08472 
    ## F-statistic: 2.037 on 5 and 51 DF,  p-value: 0.08903

``` r
EvenPC1 <- lm(even ~ PC1 + Block + TRT , RootAlphaPurp) 
summary(EvenPC1)
```

    ## 
    ## Call:
    ## lm(formula = even ~ PC1 + Block + TRT, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##        Min         1Q     Median         3Q        Max 
    ## -8.015e-04 -3.058e-04  1.160e-06  2.792e-04  9.880e-04 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.033e-02  2.037e-04  50.696   <2e-16 ***
    ## PC1         -6.197e-06  6.620e-05  -0.094    0.926    
    ## Block2      -1.131e-04  5.458e-04  -0.207    0.837    
    ## Block3      -3.416e-04  5.728e-04  -0.596    0.554    
    ## Block4      -1.537e-04  1.387e-04  -1.108    0.273    
    ## TRTInter     1.765e-04  1.300e-04   1.358    0.180    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0004309 on 51 degrees of freedom
    ## Multiple R-squared:  0.1239, Adjusted R-squared:  0.03796 
    ## F-statistic: 1.442 on 5 and 51 DF,  p-value: 0.2253

``` r
SimpPC2 <- lm(sim ~ PC2 + Block + TRT, RootAlphaPurp) 
summary(SimpPC2)
```

    ## 
    ## Call:
    ## lm(formula = sim ~ PC2 + Block + TRT, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.027561 -0.003682  0.001345  0.005432  0.011845 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.9721778  0.0029179 333.173   <2e-16 ***
    ## PC2         -0.0009071  0.0006315  -1.436    0.157    
    ## Block2       0.0059736  0.0046160   1.294    0.201    
    ## Block3       0.0038107  0.0032096   1.187    0.241    
    ## Block4       0.0001549  0.0032208   0.048    0.962    
    ## TRTInter    -0.0028856  0.0024702  -1.168    0.248    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.00813 on 51 degrees of freedom
    ## Multiple R-squared:  0.1458, Adjusted R-squared:  0.06203 
    ## F-statistic: 1.741 on 5 and 51 DF,  p-value: 0.1421

``` r
SimpInvPC2 <- lm(InvSimp ~ PC2 + Comp + TRT + Block, RootAlphaPurp) 
summary(SimpInvPC2)
```

    ## 
    ## Call:
    ## lm(formula = InvSimp ~ PC2 + Comp + TRT + Block, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -13.8113  -6.6027  -0.0319   5.2411  20.5597 
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       38.8561     3.3502  11.598 1.59e-15 ***
    ## PC2               -0.5592     0.7312  -0.765   0.4481    
    ## CompPA 4.12 Ihed  -1.9226     3.6054  -0.533   0.5963    
    ## CompPA 4.15 Ihed   0.6735     3.6180   0.186   0.8531    
    ## CompPA 4.2 Ihed   -8.6322     4.2973  -2.009   0.0502 .  
    ## CompPA 4.3 Ihed   -2.4873     4.0792  -0.610   0.5449    
    ## TRTInter               NA         NA      NA       NA    
    ## Block2             6.6866     5.4195   1.234   0.2233    
    ## Block3             4.2155     3.7854   1.114   0.2710    
    ## Block4            -1.6399     3.7052  -0.443   0.6600    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 9.225 on 48 degrees of freedom
    ## Multiple R-squared:  0.2139, Adjusted R-squared:  0.08286 
    ## F-statistic: 1.632 on 8 and 48 DF,  p-value: 0.1405

``` r
RichPC2 <- lm(rich ~ PC2 + Block + TRT, RootAlphaPurp) 
summary(RichPC2)
```

    ## 
    ## Call:
    ## lm(formula = rich ~ PC2 + Block + TRT, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -65.527 -11.917   3.751  21.397  48.599 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  438.929     10.272  42.732  < 2e-16 ***
    ## PC2           -5.734      2.223  -2.579  0.01282 *  
    ## Block2        22.156     16.249   1.363  0.17872    
    ## Block3        35.137     11.298   3.110  0.00306 ** 
    ## Block4        17.473     11.338   1.541  0.12946    
    ## TRTInter     -11.287      8.696  -1.298  0.20012    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 28.62 on 51 degrees of freedom
    ## Multiple R-squared:  0.255,  Adjusted R-squared:  0.182 
    ## F-statistic: 3.491 on 5 and 51 DF,  p-value: 0.008644

``` r
car::Anova(mod = lm(rich ~ PC2 + Block + TRT*PC2, RootAlphaPurp), type = "III")
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: rich
    ##              Sum Sq Df   F value  Pr(>F)    
    ## (Intercept) 1331420  1 1609.8076 < 2e-16 ***
    ## PC2            3433  1    4.1509 0.04692 *  
    ## Block          8117  3    3.2714 0.02865 *  
    ## TRT             750  1    0.9064 0.34565    
    ## PC2:TRT         416  1    0.5030 0.48150    
    ## Residuals     41353 50                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
EvenPC2 <- lm(even ~ PC2 + Block + TRT , RootAlphaPurp) 
summary(EvenPC2) 
```

    ## 
    ## Call:
    ## lm(formula = even ~ PC2 + Block + TRT, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##        Min         1Q     Median         3Q        Max 
    ## -7.790e-04 -2.985e-04 -1.367e-05  2.307e-04  1.030e-03 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.053e-02  1.474e-04  71.425  < 2e-16 ***
    ## PC2          7.252e-05  3.190e-05   2.273  0.02725 *  
    ## Block2      -2.146e-04  2.332e-04  -0.920  0.36176    
    ## Block3      -5.085e-04  1.621e-04  -3.136  0.00284 ** 
    ## Block4      -3.746e-04  1.627e-04  -2.303  0.02541 *  
    ## TRTInter     1.110e-04  1.248e-04   0.890  0.37769    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0004107 on 51 degrees of freedom
    ## Multiple R-squared:  0.2043, Adjusted R-squared:  0.1263 
    ## F-statistic: 2.619 on 5 and 51 DF,  p-value: 0.03493

``` r
SimpPC3 <- lm(sim ~ PC3 + Block + TRT, RootAlphaPurp) 
summary(SimpPC3)
```

    ## 
    ## Call:
    ## lm(formula = sim ~ PC3 + Block + TRT, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.027765 -0.002862  0.001317  0.006049  0.012040 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.9737709  0.0024996 389.573   <2e-16 ***
    ## PC3          0.0009578  0.0007577   1.264    0.212    
    ## Block2       0.0062684  0.0046770   1.340    0.186    
    ## Block3       0.0027008  0.0030730   0.879    0.384    
    ## Block4      -0.0028575  0.0025835  -1.106    0.274    
    ## TRTInter    -0.0029872  0.0024794  -1.205    0.234    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.008166 on 51 degrees of freedom
    ## Multiple R-squared:  0.1382, Adjusted R-squared:  0.05373 
    ## F-statistic: 1.636 on 5 and 51 DF,  p-value: 0.1673

``` r
SimpInvPC3 <- lm(InvSimp ~ PC3 + Comp + TRT + Block, RootAlphaPurp) 
summary(SimpInvPC3)
```

    ## 
    ## Call:
    ## lm(formula = InvSimp ~ PC3 + Comp + TRT + Block, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -13.3897  -7.2801   0.0318   5.5277  19.0761 
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       39.8540     2.8711  13.881   <2e-16 ***
    ## PC3                0.5259     0.8919   0.590   0.5582    
    ## CompPA 4.12 Ihed  -1.5109     3.7682  -0.401   0.6902    
    ## CompPA 4.15 Ihed   0.3203     3.5946   0.089   0.9294    
    ## CompPA 4.2 Ihed   -8.9321     4.2637  -2.095   0.0415 *  
    ## CompPA 4.3 Ihed   -2.7888     4.0417  -0.690   0.4935    
    ## TRTInter               NA         NA      NA       NA    
    ## Block2             6.6726     5.4507   1.224   0.2269    
    ## Block3             3.6502     3.6522   0.999   0.3226    
    ## Block4            -3.4573     2.9385  -1.177   0.2452    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 9.247 on 48 degrees of freedom
    ## Multiple R-squared:   0.21,  Adjusted R-squared:  0.07836 
    ## F-statistic: 1.595 on 8 and 48 DF,  p-value: 0.1513

``` r
RichPC3 <- lm(rich ~ PC3 + Block + TRT, RootAlphaPurp) 
summary(RichPC3)
```

    ## 
    ## Call:
    ## lm(formula = rich ~ PC3 + Block + TRT, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -70.160 -16.679   2.978  20.134  52.500 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 452.2049     9.2684  48.790   <2e-16 ***
    ## PC3           1.9983     2.8096   0.711    0.480    
    ## Block2       19.8187    17.3422   1.143    0.258    
    ## Block3       26.7079    11.3946   2.344    0.023 *  
    ## Block4       -0.6489     9.5797  -0.068    0.946    
    ## TRTInter    -14.8465     9.1935  -1.615    0.113    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 30.28 on 51 degrees of freedom
    ## Multiple R-squared:  0.1661, Adjusted R-squared:  0.08432 
    ## F-statistic: 2.031 on 5 and 51 DF,  p-value: 0.08978

``` r
EvenPC3 <- lm(even ~ PC3 + Block + TRT , RootAlphaPurp) 
summary(EvenPC3)
```

    ## 
    ## Call:
    ## lm(formula = even ~ PC3 + Block + TRT, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##        Min         1Q     Median         3Q        Max 
    ## -8.037e-04 -3.090e-04 -1.690e-06  2.763e-04  9.947e-04 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.034e-02  1.319e-04  78.365   <2e-16 ***
    ## PC3          2.841e-06  3.999e-05   0.071   0.9436    
    ## Block2      -1.559e-04  2.468e-04  -0.632   0.5304    
    ## Block3      -3.921e-04  1.622e-04  -2.418   0.0192 *  
    ## Block4      -1.518e-04  1.364e-04  -1.113   0.2708    
    ## TRTInter     1.763e-04  1.309e-04   1.347   0.1839    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.000431 on 51 degrees of freedom
    ## Multiple R-squared:  0.1238, Adjusted R-squared:  0.03789 
    ## F-statistic: 1.441 on 5 and 51 DF,  p-value: 0.2256

``` r
SimpPC4 <- lm(sim ~ PC4 + Block + TRT, RootAlphaPurp) 
summary(SimpPC4)
```

    ## 
    ## Call:
    ## lm(formula = sim ~ PC4 + Block + TRT, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.020413 -0.004944  0.001874  0.005665  0.012687 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.9764140  0.0024003 406.792  < 2e-16 ***
    ## PC4          0.0024831  0.0009036   2.748  0.00827 ** 
    ## Block2       0.0012787  0.0046058   0.278  0.78242    
    ## Block3       0.0006558  0.0029679   0.221  0.82601    
    ## Block4      -0.0047031  0.0025559  -1.840  0.07158 .  
    ## TRTInter    -0.0028974  0.0023099  -1.254  0.21544    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.007739 on 51 degrees of freedom
    ## Multiple R-squared:  0.2259, Adjusted R-squared:   0.15 
    ## F-statistic: 2.976 on 5 and 51 DF,  p-value: 0.01969

``` r
car::Anova(mod = lm(sim ~ PC4 + Block + TRT*PC4, RootAlphaPurp), type = "III")
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: sim
    ##             Sum Sq Df    F value Pr(>F)    
    ## (Intercept) 9.1182  1 1.5163e+05 <2e-16 ***
    ## PC4         0.0000  1 1.9700e-01 0.6591    
    ## Block       0.0003  3 1.4874e+00 0.2293    
    ## TRT         0.0001  1 8.4520e-01 0.3623    
    ## PC4:TRT     0.0000  1 7.9930e-01 0.3756    
    ## Residuals   0.0030 50                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Repeat & remove outlier
SimpPC4_out <- lm(sim ~ PC4 + Block + TRT, RootAlphaPurp %>% filter(sim > 0.94)) 
summary(SimpPC4_out)
```

    ## 
    ## Call:
    ## lm(formula = sim ~ PC4 + Block + TRT, data = RootAlphaPurp %>% 
    ##     filter(sim > 0.94))
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.017811 -0.003811  0.001589  0.004450  0.011490 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.9747198  0.0022354 436.044   <2e-16 ***
    ## PC4          0.0008668  0.0009460   0.916    0.364    
    ## Block2       0.0035840  0.0042372   0.846    0.402    
    ## Block3       0.0016636  0.0027118   0.613    0.542    
    ## Block4      -0.0021819  0.0024349  -0.896    0.374    
    ## TRTInter    -0.0026408  0.0020996  -1.258    0.214    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.00703 on 50 degrees of freedom
    ## Multiple R-squared:  0.102,  Adjusted R-squared:  0.01223 
    ## F-statistic: 1.136 on 5 and 50 DF,  p-value: 0.3537

``` r
SimpInvPC4 <- lm(InvSimp ~ PC4 + Comp + TRT + Block, RootAlphaPurp) 
summary(SimpInvPC4)
```

    ## 
    ## Call:
    ## lm(formula = InvSimp ~ PC4 + Comp + TRT + Block, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -13.1982  -7.7521   0.3428   6.3656  21.7090 
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       41.6988     2.8055  14.863   <2e-16 ***
    ## PC4                1.9003     1.0889   1.745   0.0873 .  
    ## CompPA 4.12 Ihed  -2.0636     3.5041  -0.589   0.5587    
    ## CompPA 4.15 Ihed   0.2956     3.4983   0.084   0.9330    
    ## CompPA 4.2 Ihed   -8.0009     4.1457  -1.930   0.0595 .  
    ## CompPA 4.3 Ihed   -1.7723     3.9605  -0.448   0.6565    
    ## TRTInter               NA         NA      NA       NA    
    ## Block2             3.0685     5.5619   0.552   0.5837    
    ## Block3             2.1513     3.5765   0.602   0.5503    
    ## Block4            -4.8794     2.9836  -1.635   0.1085    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 9 on 48 degrees of freedom
    ## Multiple R-squared:  0.2518, Adjusted R-squared:  0.1271 
    ## F-statistic: 2.019 on 8 and 48 DF,  p-value: 0.06403

``` r
RichPC4 <- lm(rich ~ PC4 + Block + TRT, RootAlphaPurp) 
summary(RichPC4)
```

    ## 
    ## Call:
    ## lm(formula = rich ~ PC4 + Block + TRT, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -65.92 -15.72   0.98  20.30  53.09 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  456.445      9.346  48.836   <2e-16 ***
    ## PC4            3.503      3.518   0.996   0.3241    
    ## Block2        12.109     17.934   0.675   0.5026    
    ## Block3        23.597     11.557   2.042   0.0464 *  
    ## Block4        -3.106      9.953  -0.312   0.7563    
    ## TRTInter     -15.185      8.995  -1.688   0.0975 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 30.14 on 51 degrees of freedom
    ## Multiple R-squared:  0.1739, Adjusted R-squared:  0.09287 
    ## F-statistic: 2.147 on 5 and 51 DF,  p-value: 0.07469

``` r
EvenPC4 <- lm(even ~ PC4 + Block + TRT , RootAlphaPurp) 
summary(EvenPC4)
```

    ## 
    ## Call:
    ## lm(formula = even ~ PC4 + Block + TRT, data = RootAlphaPurp)
    ## 
    ## Residuals:
    ##        Min         1Q     Median         3Q        Max 
    ## -7.993e-04 -2.799e-04  2.209e-05  2.921e-04  9.660e-04 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.036e-02  1.333e-04  77.714   <2e-16 ***
    ## PC4          2.615e-05  5.018e-05   0.521   0.6045    
    ## Block2      -2.010e-04  2.558e-04  -0.786   0.4357    
    ## Block3      -4.111e-04  1.648e-04  -2.494   0.0159 *  
    ## Block4      -1.729e-04  1.420e-04  -1.218   0.2288    
    ## TRTInter     1.824e-04  1.283e-04   1.422   0.1611    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0004298 on 51 degrees of freedom
    ## Multiple R-squared:  0.1283, Adjusted R-squared:  0.04289 
    ## F-statistic: 1.502 on 5 and 51 DF,  p-value: 0.2057

## Plotting significant linear associations

``` r
# Richness and root architecture
P2.rich <- ggplot() +
  geom_point(data = RootAlphaPurp, aes(PC2,rich), alpha = 0.5, size = 5) +
  geom_smooth(data = RootAlphaPurp, method = "lm", aes(PC2,rich), fullrange = TRUE) +
  theme_classic() +
  ylab("Richness") +
  xlab("") +
  Tx2
#+ 
  #annotate("text", x = -1, y = 500, label = "paste(italic(R) ^ 2, \" = 0.15\")", parse = TRUE,hjust=0, size=5) + 
 # annotate("text", x = -1, y = 490, label = "paste(italic(B), \" = -4.65 +/-1.81\")", parse = TRUE,hjust=0, size=5) +
 # annotate("text", x = -1, y = 480, label = "paste(italic(P), \" = 0.01\")", parse = TRUE,hjust=0, size=5)

# Evenness and root architecture
P2.even <- ggplot() +
  geom_point(data = RootAlphaPurp, aes(PC2, even), alpha = 0.5, size = 5) +
  geom_smooth(data = RootAlphaPurp, method = "lm", aes(PC2, even), fullrange = TRUE) +
  theme_classic() +
  ylab("Evenness") +
  xlab("") +
  Tx2#+ 
  #annotate("text", x = -4.5, y = .0114, label = "paste(italic(R) ^ 2, \" = 0.12\")", parse = TRUE,hjust=0, size=5) + 
#  annotate("text", x = -4.5, y = .0112, label = "paste(italic(B), \" = -7.29E-5 +/-3.28E-5\")", parse = TRUE,hjust=0, size=5) +
 # annotate("text", x = -4.5, y =.011, label = "paste(italic(P), \" = 0.04\")", parse = TRUE,hjust=0, size=5)


# Root morphology on species diversity Simpson metric

P4.Sim <- ggplot() +
  geom_point(data = RootAlphaPurp, aes(PC4, sim), alpha = 0.5, size = 5) +
  geom_smooth(data = RootAlphaPurp, method = "lm", aes(PC4, sim), fullrange = TRUE) +
  theme_classic() +
  ylab("Simpson") +
  xlab("") +
  Tx2#+
  #annotate("text", x = -5, y = .988, label = "paste(italic(R) ^ 2, \" = 0.18\")", parse = TRUE,hjust=0, size=5) + 
  #annotate("text", x = -5, y = .985, label = "paste(italic(B), \" = 2.22E-4 +/-6.99E-4\")", parse = TRUE,hjust=0, size=5) +
 # annotate("text", x =-5, y = .982, label = "paste(italic(P), \" < 0.01\")", parse = TRUE,hjust=0, size=5)


P4.simIn <- ggplot() +
  geom_point(data = RootAlphaPurp, aes(PC4, InvSimp), alpha = 0.5, size = 5) +
  geom_smooth(data = RootAlphaPurp, method = "lm", aes(PC4, InvSimp), fullrange = TRUE) +
  theme_classic() +
  ylab("Inverse Simpson") +
  xlab("") +
  Tx2#+
 # annotate("text", x = -5, y = 60, label = "paste(italic(R) ^ 2, \" = 0.12\")", parse = TRUE,hjust=0, size=5) + 
 # annotate("text", x = -5, y = 57, label = "paste(italic(B), \" = 2.08 +/-1.05\")", parse = TRUE,hjust=0, size=5) +
 # annotate("text", x =-5, y = 54, label = "paste(italic(P), \" = 0.02\")", parse = TRUE,hjust=0, size=5)
```

## Linear mixed models (not reported)

``` r
### Simpson
SimpLMM <- lmer(sim ~ TRT + Block + (1|ML), alpha %>% filter(Species == "Ip")) 
```

    ## boundary (singular) fit: see ?isSingular

``` r
anova(SimpLMM)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq    Mean Sq NumDF DenDF F value  Pr(>F)  
    ## TRT   0.00003472 3.4721e-05     1    95  0.5656 0.45388  
    ## Block 0.00050911 1.6970e-04     3    95  2.7643 0.04617 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
ranova(SimpLMM)
```

    ## ANOVA-like table for random-effects: Single term deletions
    ## 
    ## Model:
    ## sim ~ TRT + Block + (1 | ML)
    ##          npar logLik     AIC         LRT Df Pr(>Chisq)
    ## <none>      7 317.94 -621.89                          
    ## (1 | ML)    6 317.94 -623.89 -4.5475e-13  1          1

``` r
### Inverse Simpson
SimpInvLMM <- lmer(InvSimp ~ TRT + Block + (1|ML), alpha %>% filter(Species == "Ip")) 
```

    ## boundary (singular) fit: see ?isSingular

``` r
anova(SimpInvLMM)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##        Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
    ## TRT      0.48    0.48     1    95  0.0055 0.94112  
    ## Block 1031.56  343.85     3    95  3.9020 0.01122 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
ranova(SimpInvLMM)
```

    ## ANOVA-like table for random-effects: Single term deletions
    ## 
    ## Model:
    ## InvSimp ~ TRT + Block + (1 | ML)
    ##          npar  logLik    AIC LRT Df Pr(>Chisq)
    ## <none>      7 -355.46 724.92                  
    ## (1 | ML)    6 -355.46 722.92   0  1          1

``` r
### Inverse Simpson
RichLMM <- lmer(rich ~ TRT + Block + (1|ML), alpha %>% filter(Species == "Ip")) 
anova(RichLMM)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##        Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
    ## TRT    1393.7  1393.7     1 91.101  1.5271 0.219729   
    ## Block 14143.2  4714.4     3 92.228  5.1656 0.002407 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
ranova(RichLMM)
```

    ## ANOVA-like table for random-effects: Single term deletions
    ## 
    ## Model:
    ## rich ~ TRT + Block + (1 | ML)
    ##          npar  logLik    AIC     LRT Df Pr(>Chisq)
    ## <none>      7 -467.16 948.31                      
    ## (1 | ML)    6 -467.22 946.44 0.12747  1     0.7211

``` r
LeafData <- read.csv("../DataSets/SizeFitData.csv")

LeafData$Sample_ID <- paste(LeafData$Position, ifelse(grepl("Ihed", LeafData$ML), "H", "P"), sep="")
```

``` r
Fitness = read.csv("../DataSets/FitPA4.csv")

library(dplyr)
# Calculate relative fitness
  # First calculate mean seed number by species and treatment---note* we only have seed output of I. purpurea
MeanSeedNumber <- aggregate(SeedNumber ~ Trt + Species, Fitness, mean)

colnames(MeanSeedNumber) <- c("Trt", "Species", "MeanSeedNumber")

Ipurp.Fit <- Fitness %>%
  filter(Species == "Ip")

Ipurp.Alpha <- alpha %>%
  filter(Species == "Ip")

FitnessPurp <- merge(Ipurp.Fit, MeanSeedNumber, by=c("Trt", "Species"))
FitnessPurp$RelativeFit <- FitnessPurp$SeedNumber/FitnessPurp$MeanSeedNumber
FitnessPurp$Block <- as.factor(FitnessPurp$Block)

FitnessPurp2 <- merge(FitnessPurp, LeafData)
str(FitnessPurp2)
```

    ## 'data.frame':    385 obs. of  24 variables:
    ##  $ Trt            : chr  "Alone" "Alone" "Alone" "Alone" ...
    ##  $ Species        : chr  "Ip" "Ip" "Ip" "Ip" ...
    ##  $ Position       : chr  "102" "114" "136" "138" ...
    ##  $ Block          : Factor w/ 4 levels "1","2","3","4": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Order          : int  102 114 136 138 151 153 172 203 204 213 ...
    ##  $ ML             : chr  "PA4.2Ip" "PA4.15Ip" "PA4.2Ip" "PA4.11Ip" ...
    ##  $ UniqId         : chr  "102Ip" "114Ip" "136Ip" "138Ip" ...
    ##  $ GerminationDate: chr  "43258" "43258" "43258" "43258" ...
    ##  $ Comment1       : chr  NA NA NA NA ...
    ##  $ Dead_plant     : chr  NA NA NA NA ...
    ##  $ DeathCause     : logi  NA NA NA NA NA NA ...
    ##  $ RootsHarvested : chr  "N" "N" "N" "N" ...
    ##  $ SeedsCounted   : chr  "Y" "Y" "Y" "Y" ...
    ##  $ Comp           : chr  NA NA NA NA ...
    ##  $ Combos         : chr  NA NA NA NA ...
    ##  $ Population     : chr  "PA4" "PA4" "PA4" "PA4" ...
    ##  $ SeedNumber     : int  96 75 93 56 383 321 256 202 79 172 ...
    ##  $ MeanSeedNumber : num  224 224 224 224 224 ...
    ##  $ RelativeFit    : num  0.428 0.334 0.415 0.25 1.708 ...
    ##  $ Block.1        : int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Leaf.Number    : int  6 5 6 7 9 10 11 7 7 8 ...
    ##  $ Comment        : logi  NA NA NA NA NA NA ...
    ##  $ X              : int  1 7 21 23 31 32 39 53 54 57 ...
    ##  $ Sample_ID      : chr  "102P" "114P" "136P" "138P" ...

``` r
FitnessPurp2$Leaf.Number <- as.numeric(as.character(FitnessPurp2$Leaf.Number))

SN1 <- lmer(SeedNumber~Trt + Block + Leaf.Number + Block:Trt + (1|ML), FitnessPurp2)
anova(SN1) # Treatment is significant effect on plant seed number; Significant Treatment by Block effect
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##              Sum Sq Mean Sq NumDF  DenDF  F value  Pr(>F)    
    ## Trt           83491   83491     1 370.05   3.9783 0.04682 *  
    ## Block        128423   42808     3 371.69   2.0398 0.10790    
    ## Leaf.Number 2425524 2425524     1 374.69 115.5752 < 2e-16 ***
    ## Trt:Block    164743   54914     3 371.54   2.6166 0.05082 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
ggplot(FitnessPurp, aes(Trt, RelativeFit, fill=Block)) +
  geom_boxplot() +
  scale_fill_brewer("Paired") +
  theme_classic() +
  ylab("Relative Fitness") +
  ggtitle("Relative Fitness by Treatment and Block")
```

![](README_files/figure-gfm/Read%20in%20fitness%20data-1.png)<!-- -->

``` r
# Remove size effects from fitness
StdFitness <- FitnessPurp2[c("Trt", "Species", "Block", "ML", "RelativeFit", "Leaf.Number")] # Subset fitness data for variables of interest

# Run one-way ANOVA to remove size effect--i.e., keep residuals
StdFitness$RelativeFitness <- residuals(lm(RelativeFit ~ Leaf.Number, FitnessPurp2)) 

# Compare residuals and non-standard values of fitness
plot(StdFitness$RelativeFit, StdFitness$RelativeFitness)
```

![](README_files/figure-gfm/Read%20in%20fitness%20data-2.png)<!-- -->

``` r
# Average fitness by block, maternal line and treatment--we use NON standardized fitness (ie size effects not removed)
FitAveraged = aggregate(RelativeFit ~ Block + Trt + ML, FitnessPurp, mean)
colnames(FitAveraged) <- c("Block", "Trt", "ML", "RelativeFitness")
FitAveraged$TRT <- FitAveraged$Trt
dim(FitAveraged)
```

    ## [1] 76  5

``` r
head(FitAveraged)
```

    ##   Block   Trt       ML RelativeFitness   TRT
    ## 1     1 Alone PA4.11Ip       0.4437172 Alone
    ## 2     2 Alone PA4.11Ip       1.0435158 Alone
    ## 3     3 Alone PA4.11Ip       1.0747321 Alone
    ## 4     4 Alone PA4.11Ip       1.4627059 Alone
    ## 5     1 Inter PA4.11Ip       0.6619616 Inter
    ## 6     2 Inter PA4.11Ip       0.6103389 Inter

``` r
# Average size
SizeAveraged <- aggregate(Leaf.Number ~ Block + Trt + ML + Species, LeafData,mean)
colnames(SizeAveraged) <- c("Block", "Trt", "ML", "Species", "Size")

SizeAveraged$TRT <- SizeAveraged$Trt
dim(SizeAveraged)
```

    ## [1] 76  6

``` r
head(SizeAveraged)
```

    ##   Block   Trt       ML Species      Size   TRT
    ## 1     1 Alone PA4.11Ip      Ip  9.000000 Alone
    ## 2     2 Alone PA4.11Ip      Ip 11.000000 Alone
    ## 3     3 Alone PA4.11Ip      Ip  8.000000 Alone
    ## 4     4 Alone PA4.11Ip      Ip  5.000000 Alone
    ## 5     1 Inter PA4.11Ip      Ip  7.545455 Inter
    ## 6     2 Inter PA4.11Ip      Ip  7.500000 Inter

``` r
SizePurp <- SizeAveraged %>% filter(Species == "Ip")
```

``` r
####  ####  ####  ####  ####  ####   ####
#  Examine selection on microbiome first
####  ####  ####  ####  ####  ####   ####
FitnessPurp$TRT <- FitnessPurp$Trt
FitnessPurp$Combos <- as.character(FitnessPurp$Combos)
FitnessPurp[which(FitnessPurp$TRT == "Alone"),]$Combos <- "none"
FitnessPurp$Combos <- as.factor(FitnessPurp$Combos)

BrayFit <- merge(physeq.pcoa.df,FitAveraged)

FitAlpha <- merge(FitAveraged, alpha)
dim(FitAlpha)
```

    ## [1] 97 13

``` r
ggplot(FitAlpha, aes(TRT, RelativeFitness)) +
  geom_boxplot() +
  scale_fill_brewer("Paired") +
  theme_classic() +
  ylab("Relative Fitness") +
  ggtitle("Relative Fitness by Treatment and Block") +
  facet_grid(~Block)
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# Combine with root data
library(dplyr)

RootAveraged <- aggregate(list(RootData[c("PC1", "PC2", "PC3", "PC4")]),by=list(RootData$Trt, RootData$ML),FUN=mean) 

colnames(RootAveraged) <- c("Trt", "ML", "PC1", "PC2", "PC3", "PC4")
head(RootAveraged)
```

    ##     Trt         ML         PC1        PC2        PC3         PC4
    ## 1 Alone   PA4.11Ip -0.69794646 -1.3274069  0.7756816 -0.75520608
    ## 2 Inter   PA4.11Ip  0.43519578  0.6954043 -0.1859846  0.10761866
    ## 3 Inter PA4.12Ihed  0.07075942 -0.9146971 -0.2814713 -0.01631047
    ## 4 Alone   PA4.13Ip -2.28254650  1.4831488  0.9299528 -0.64140759
    ## 5 Inter   PA4.13Ip  0.11317327 -0.1860671  0.1006711  0.43992946
    ## 6 Alone   PA4.14Ip  0.70071833 -0.4492977  0.9728631 -0.41237743

``` r
RootFitAlpha <- merge(FitAlpha, RootAveraged)
head(RootFitAlpha)
```

    ##         ML   Trt Block   TRT RelativeFitness Sample_ID Species
    ## 1 PA4.11Ip Alone     1 Alone       0.4437172      158P      Ip
    ## 2 PA4.11Ip Alone     1 Alone       0.4437172       84P      Ip
    ## 3 PA4.11Ip Alone     4 Alone       1.4627059      814P      Ip
    ## 4 PA4.11Ip Alone     3 Alone       1.0747321      674P      Ip
    ## 5 PA4.11Ip Inter     1 Inter       0.6619616      93AP      Ip
    ## 6 PA4.11Ip Inter     1 Inter       0.6619616       21P      Ip
    ##                    Combos rich  InvSimp       sim     shan        even
    ## 1                    none  434 31.43535 0.9681887 4.534697 0.010448610
    ## 2                    none  478 41.08728 0.9756616 4.753743 0.009945068
    ## 3                    none  466 48.22287 0.9792630 4.791226 0.010281601
    ## 4                    none  437 30.30030 0.9669970 4.491512 0.010278060
    ## 5 PA 4.11 Ip-PA 4.15 Ihed  427 35.29942 0.9716709 4.593859 0.010758451
    ## 6  PA 4.11 Ip-PA 4.3 Ihed  489 42.06259 0.9762259 4.752692 0.009719206
    ##          PC1        PC2        PC3        PC4
    ## 1 -0.6979465 -1.3274069  0.7756816 -0.7552061
    ## 2 -0.6979465 -1.3274069  0.7756816 -0.7552061
    ## 3 -0.6979465 -1.3274069  0.7756816 -0.7552061
    ## 4 -0.6979465 -1.3274069  0.7756816 -0.7552061
    ## 5  0.4351958  0.6954043 -0.1859846  0.1076187
    ## 6  0.4351958  0.6954043 -0.1859846  0.1076187

``` r
#RootFitLeafAlpha=merge(RootFitAlpha, LeafData[c("Sample_ID", "Leaf.Number")],by="Sample_ID")

#RootFitLeafAlpha$Leaf.Number=as.numeric(as.character(RootFitLeafAlpha$Leaf.Number))


# CombinE root,fitness/bray estimates
RootFitBray <- merge(BrayFit, RootAveraged)


# Plot boxplots of averaged root traits
PC2_Box <- ggplot(RootAveraged, aes(x = "", y = PC2)) +
  geom_boxplot() +
  #geom_jitter() +
  xlab("Root Architecture") +
  theme_classic() +
  Tx+ 
  coord_flip()

# Plot boxplots of averaged root traits
PC1_Box=ggplot(RootAveraged, aes(x = "", y = PC1)) +
  geom_boxplot() +
  #geom_jitter() +
  xlab("Root Topology") +
  theme_classic() +
  Tx+ 
  coord_flip()

# Plot boxplots of averaged root traits
PC4_Box=ggplot(RootAveraged, aes(x = "", y = PC4)) +
  geom_boxplot() +
  #geom_jitter() +
  xlab("Root Morphology") +
  theme_classic() +
  Tx+ 
  coord_flip()

#ggarrange(P2.rich,P2.even,P4.Sim,P4.simIn,nrow=2,ncol=2)
AB <- cowplot::plot_grid(P2.rich,P2.even, align = "hv", ncol = 2, Labels = c("A", "B"), label_size = 22 ,hjust = -11,vjust = 1, Label_x = -0.16)
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

    ## Warning in as_grob.default(plot): Cannot convert object of class character into
    ## a grob.

    ## Warning in as_grob.default(plot): Cannot convert object of class numeric into a
    ## grob.

    ## Warning: Graphs cannot be vertically aligned unless the axis parameter is set.
    ## Placing graphs unaligned.

    ## Warning: Graphs cannot be horizontally aligned unless the axis parameter is set.
    ## Placing graphs unaligned.

``` r
ab <- cowplot::plot_grid(P4.Sim,P4.simIn, align = "hv", ncol = 2, Labels =c("A", "B"), label_size =22 ,hjust = -11,vjust=1, Label_x = -0.17)
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

    ## Warning in as_grob.default(plot): Cannot convert object of class character into
    ## a grob.

    ## Warning in as_grob.default(plot): Cannot convert object of class numeric into a
    ## grob.

    ## Warning: Graphs cannot be vertically aligned unless the axis parameter is set.
    ## Placing graphs unaligned.

    ## Warning: Graphs cannot be horizontally aligned unless the axis parameter is set.
    ## Placing graphs unaligned.

``` r
P4.simIn +
  xlab("Root morphology (PC4)")
```

    ## `geom_smooth()` using formula 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
# Common x title

x.grob <- textGrob("Root morphology (PC4)", 
                   gp=gpar(col="black", fontsize = 25), rot=0)

gridExtra::grid.arrange(gridExtra::arrangeGrob(ab,bottom = x.grob, padding = unit(0.5,units = 'in'), nrow=1))
```

![](README_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->

``` r
# Common x title
x.grob1 <- textGrob("Root architecture (PC2)", 
                   gp = gpar(col="black", fontsize = 25), rot = 0)

gridExtra::grid.arrange(gridExtra::arrangeGrob(AB, bottom = x.grob1, padding = unit(0.5,units = 'in'), nrow=1))
```

![](README_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->

## Selection on microbe variables

``` r
# Quadratic microbe term
FitAlpha$even2 <- (FitAlpha$even*FitAlpha$even)
FitAlpha$rich2 <- (FitAlpha$rich*FitAlpha$rich)
FitAlpha$InSim2 <- (FitAlpha$InvSimp*FitAlpha$InvSimp)
FitAlpha$sim2 <- (FitAlpha$sim*FitAlpha$sim)
```

### Selection on richness

``` r
summary(lm(RelativeFitness ~ Block +rich+rich2, FitAlpha %>% filter(TRT == "Alone")))
```

    ## 
    ## Call:
    ## lm(formula = RelativeFitness ~ Block + rich + rich2, data = FitAlpha %>% 
    ##     filter(TRT == "Alone"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.80608 -0.16098  0.00677  0.17969  0.91777 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  6.384e+00  1.682e+01   0.379   0.7088  
    ## Block2       2.126e-01  2.649e-01   0.803   0.4326  
    ## Block3       4.705e-01  3.024e-01   1.556   0.1372  
    ## Block4       6.100e-01  2.469e-01   2.470   0.0237 *
    ## rich        -2.997e-02  7.405e-02  -0.405   0.6904  
    ## rich2        3.768e-05  8.153e-05   0.462   0.6495  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4318 on 18 degrees of freedom
    ## Multiple R-squared:  0.3867, Adjusted R-squared:  0.2164 
    ## F-statistic:  2.27 on 5 and 18 DF,  p-value: 0.09127

``` r
summary(lm(RelativeFitness ~ Block +rich+rich2, FitAlpha %>% filter(TRT != "Alone")))
```

    ## 
    ## Call:
    ## lm(formula = RelativeFitness ~ Block + rich + rich2, data = FitAlpha %>% 
    ##     filter(TRT != "Alone"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.38834 -0.11485 -0.01772  0.06251  0.67553 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -4.172e+00  3.263e+00  -1.278    0.205    
    ## Block2       7.177e-02  6.761e-02   1.062    0.292    
    ## Block3       7.769e-01  6.995e-02  11.107  < 2e-16 ***
    ## Block4       3.054e-01  7.026e-02   4.346 4.82e-05 ***
    ## rich         2.284e-02  1.496e-02   1.527    0.131    
    ## rich2       -2.663e-05  1.707e-05  -1.560    0.123    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2006 on 67 degrees of freedom
    ## Multiple R-squared:  0.7156, Adjusted R-squared:  0.6943 
    ## F-statistic: 33.71 on 5 and 67 DF,  p-value: < 2.2e-16

``` r
richness_sel <- lm(RelativeFitness ~ TRT*Block +rich, FitAlpha)
summary(richness_sel) # No evidence of indirect selection
```

    ## 
    ## Call:
    ## lm(formula = RelativeFitness ~ TRT * Block + rich, data = FitAlpha)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.71931 -0.13206 -0.04529  0.12517  0.97971 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      0.3412789  0.4489876   0.760 0.449222    
    ## TRTInter         0.1232141  0.1304234   0.945 0.347387    
    ## Block2           0.1579627  0.1580477   0.999 0.320311    
    ## Block3           0.6239940  0.1661699   3.755 0.000310 ***
    ## Block4           0.5893259  0.1513870   3.893 0.000192 ***
    ## rich             0.0005081  0.0009493   0.535 0.593839    
    ## TRTInter:Block2 -0.0729686  0.1823258  -0.400 0.689972    
    ## TRTInter:Block3  0.1227196  0.1872063   0.656 0.513836    
    ## TRTInter:Block4 -0.2617882  0.1777344  -1.473 0.144341    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2715 on 88 degrees of freedom
    ## Multiple R-squared:  0.5695, Adjusted R-squared:  0.5304 
    ## F-statistic: 14.55 on 8 and 88 DF,  p-value: 2.467e-13

``` r
anova(richness_sel)
```

    ## Analysis of Variance Table
    ## 
    ## Response: RelativeFitness
    ##           Df Sum Sq Mean Sq F value    Pr(>F)    
    ## TRT        1 0.1156 0.11565  1.5695    0.2136    
    ## Block      3 8.0999 2.69998 36.6409 1.842e-15 ***
    ## rich       1 0.0128 0.01276  0.1732    0.6783    
    ## TRT:Block  3 0.3508 0.11694  1.5870    0.1982    
    ## Residuals 88 6.4845 0.07369                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
ggplot(FitAlpha) +
  geom_point(aes(RelativeFitness, rich), alpha = 0.5, size = 5) +
  geom_smooth(method = "lm", aes(RelativeFitness, rich), fullrange = TRUE) +
  theme_classic() +
  ylab("Relative fitness") +
  xlab("Richness") +
  Tx2
```

    ## `geom_smooth()` using formula 'y ~ x'

![](README_files/figure-gfm/cor%20on%20richness-1.png)<!-- -->

### Selection on Inverse Simpson

``` r
summary(lm(RelativeFitness ~ Block +InvSimp+InSim2, FitAlpha %>% filter(TRT == "Alone")))
```

    ## 
    ## Call:
    ## lm(formula = RelativeFitness ~ Block + InvSimp + InSim2, data = FitAlpha %>% 
    ##     filter(TRT == "Alone"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.72488 -0.14886 -0.04733  0.14805  0.99197 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  0.8393042  2.4296123   0.345   0.7338  
    ## Block2       0.1485885  0.2694609   0.551   0.5881  
    ## Block3       0.6201396  0.2900918   2.138   0.0465 *
    ## Block4       0.5778902  0.2648341   2.182   0.0426 *
    ## InvSimp     -0.0160280  0.1260111  -0.127   0.9002  
    ## InSim2       0.0002334  0.0016387   0.142   0.8883  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4528 on 18 degrees of freedom
    ## Multiple R-squared:  0.3256, Adjusted R-squared:  0.1382 
    ## F-statistic: 1.738 on 5 and 18 DF,  p-value: 0.177

``` r
summary(lm(RelativeFitness ~ Block +InvSimp+InSim2, FitAlpha %>% filter(TRT != "Alone")))
```

    ## 
    ## Call:
    ## lm(formula = RelativeFitness ~ Block + InvSimp + InSim2, data = FitAlpha %>% 
    ##     filter(TRT != "Alone"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.39213 -0.11007 -0.02050  0.06454  0.71716 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  8.528e-01  2.654e-01   3.214  0.00201 ** 
    ## Block2       9.440e-02  6.875e-02   1.373  0.17427    
    ## Block3       7.751e-01  7.055e-02  10.986  < 2e-16 ***
    ## Block4       3.224e-01  7.012e-02   4.597 1.95e-05 ***
    ## InvSimp     -6.992e-03  1.337e-02  -0.523  0.60266    
    ## InSim2       6.198e-05  1.640e-04   0.378  0.70666    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2035 on 67 degrees of freedom
    ## Multiple R-squared:  0.7073, Adjusted R-squared:  0.6854 
    ## F-statistic: 32.38 on 5 and 67 DF,  p-value: < 2.2e-16

``` r
Invn_sel <- lm(RelativeFitness ~ TRT*Block + InvSimp + Block, FitAlpha)
summary(Invn_sel) # No evidence of indirect selection
```

    ## 
    ## Call:
    ## lm(formula = RelativeFitness ~ TRT * Block + InvSimp + Block, 
    ##     data = FitAlpha)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.70503 -0.14215 -0.02584  0.10658  0.99114 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      0.633031   0.163751   3.866 0.000212 ***
    ## TRTInter         0.108178   0.129252   0.837 0.404889    
    ## Block2           0.141277   0.157210   0.899 0.371291    
    ## Block3           0.640501   0.164571   3.892 0.000193 ***
    ## Block4           0.579859   0.151270   3.833 0.000237 ***
    ## InvSimp         -0.001500   0.003071  -0.488 0.626428    
    ## TRTInter:Block2 -0.051356   0.182252  -0.282 0.778769    
    ## TRTInter:Block3  0.133148   0.188112   0.708 0.480934    
    ## TRTInter:Block4 -0.256353   0.177467  -1.445 0.152146    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2715 on 88 degrees of freedom
    ## Multiple R-squared:  0.5693, Adjusted R-squared:  0.5301 
    ## F-statistic: 14.54 on 8 and 88 DF,  p-value: 2.523e-13

``` r
ggplot(FitAlpha) +
  geom_point(aes(RelativeFitness, InvSimp), alpha = 0.5, size = 5) +
  geom_smooth(method = "lm", aes(RelativeFitness, InvSimp), fullrange = TRUE) +
  theme_classic() +
  ylab("Relative fitness") +
  xlab("Inverse Simpson") +
  Tx2
```

    ## `geom_smooth()` using formula 'y ~ x'

![](README_files/figure-gfm/cor%20on%20Inverse%20Simpson-1.png)<!-- -->

## Selection on Simpson

``` r
summary(lm(RelativeFitness ~ Block + sim + sim2, FitAlpha %>% filter(TRT == "Alone")))
```

    ## 
    ## Call:
    ## lm(formula = RelativeFitness ~ Block + sim + sim2, data = FitAlpha %>% 
    ##     filter(TRT == "Alone"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.72245 -0.15050 -0.02453  0.14432  0.98295 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  185.7703  2769.9476   0.067   0.9473  
    ## Block2         0.1520     0.2688   0.565   0.5787  
    ## Block3         0.6314     0.2814   2.244   0.0377 *
    ## Block4         0.5850     0.2636   2.219   0.0395 *
    ## sim         -383.4087  5706.0154  -0.067   0.9472  
    ## sim2         198.4254  2938.5243   0.068   0.9469  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4531 on 18 degrees of freedom
    ## Multiple R-squared:  0.3247, Adjusted R-squared:  0.1372 
    ## F-statistic: 1.731 on 5 and 18 DF,  p-value: 0.1784

``` r
summary(lm(RelativeFitness ~ Block + sim + sim2, FitAlpha %>% filter(TRT != "Alone")))
```

    ## 
    ## Call:
    ## lm(formula = RelativeFitness ~ Block + sim + sim2, data = FitAlpha %>% 
    ##     filter(TRT != "Alone"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.39944 -0.12749 -0.02963  0.06538  0.70052 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -163.39346  171.28679  -0.954    0.344    
    ## Block2         0.08865    0.06850   1.294    0.200    
    ## Block3         0.78100    0.06940  11.253  < 2e-16 ***
    ## Block4         0.32009    0.06979   4.586 2.03e-05 ***
    ## sim          342.11116  355.00265   0.964    0.339    
    ## sim2        -178.28712  183.92949  -0.969    0.336    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2026 on 67 degrees of freedom
    ## Multiple R-squared:  0.7099, Adjusted R-squared:  0.6882 
    ## F-statistic: 32.78 on 5 and 67 DF,  p-value: < 2.2e-16

``` r
sim_sel <- lm(RelativeFitness ~ TRT*Block + sim, FitAlpha)
summary(sim_sel) # NS
```

    ## 
    ## Call:
    ## lm(formula = RelativeFitness ~ TRT * Block + sim, data = FitAlpha)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.70300 -0.13974 -0.02201  0.09864  0.99092 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      2.08771    3.55641   0.587 0.558688    
    ## TRTInter         0.10593    0.12992   0.815 0.417094    
    ## Block2           0.14195    0.15727   0.903 0.369198    
    ## Block3           0.63833    0.16450   3.880 0.000201 ***
    ## Block4           0.57940    0.15145   3.826 0.000243 ***
    ## sim             -1.55368    3.64892  -0.426 0.671299    
    ## TRTInter:Block2 -0.05088    0.18283  -0.278 0.781456    
    ## TRTInter:Block3  0.13200    0.18817   0.701 0.484844    
    ## TRTInter:Block4 -0.25461    0.17756  -1.434 0.155126    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2716 on 88 degrees of freedom
    ## Multiple R-squared:  0.569,  Adjusted R-squared:  0.5298 
    ## F-statistic: 14.52 on 8 and 88 DF,  p-value: 2.593e-13

``` r
ggplot(FitAlpha) +
  geom_point(aes(RelativeFitness, sim), alpha = 0.5, size = 5) +
  geom_smooth(method = "lm", aes(RelativeFitness, sim), fullrange = TRUE) +
  theme_classic() +
  ylab("Relative fitness") +
  xlab("Simpson") +
  Tx2
```

    ## `geom_smooth()` using formula 'y ~ x'

![](README_files/figure-gfm/cor%20on%20Simpson-1.png)<!-- -->

## Selection on Evenness

``` r
summary(lm(RelativeFitness ~ Block + even + even2 , FitAlpha %>% filter(TRT == "Alone")))
```

    ## 
    ## Call:
    ## lm(formula = RelativeFitness ~ Block + even + even2, data = FitAlpha %>% 
    ##     filter(TRT == "Alone"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.85960 -0.20018  0.00371  0.17956  0.91880 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  1.195e+01  3.982e+01   0.300   0.7675  
    ## Block2       2.294e-01  2.599e-01   0.883   0.3890  
    ## Block3       4.721e-01  2.847e-01   1.658   0.1147  
    ## Block4       5.963e-01  2.456e-01   2.428   0.0259 *
    ## even        -1.867e+03  7.813e+03  -0.239   0.8138  
    ## even2        7.376e+04  3.833e+05   0.192   0.8496  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4246 on 18 degrees of freedom
    ## Multiple R-squared:  0.4071, Adjusted R-squared:  0.2424 
    ## F-statistic: 2.472 on 5 and 18 DF,  p-value: 0.07149

``` r
summary(lm(RelativeFitness ~ Block + even + even2, FitAlpha %>% filter(TRT != "Alone"))) # Significant quadratic term for evenness on fitnes (super low slope tho)
```

    ## 
    ## Call:
    ## lm(formula = RelativeFitness ~ Block + even + even2, data = FitAlpha %>% 
    ##     filter(TRT != "Alone"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.41978 -0.11110 -0.01497  0.08058  0.67279 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -2.003e+01  8.585e+00  -2.334   0.0226 *  
    ## Block2       1.023e-01  6.608e-02   1.548   0.1264    
    ## Block3       7.702e-01  6.764e-02  11.388  < 2e-16 ***
    ## Block4       3.111e-01  6.810e-02   4.569 2.16e-05 ***
    ## even         3.961e+03  1.641e+03   2.414   0.0185 *  
    ## even2       -1.890e+05  7.838e+04  -2.411   0.0187 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1963 on 67 degrees of freedom
    ## Multiple R-squared:  0.7274, Adjusted R-squared:  0.7071 
    ## F-statistic: 35.76 on 5 and 67 DF,  p-value: < 2.2e-16

``` r
even_sel <- lm(RelativeFitness ~ TRT*Block + TRT*even2 + even, FitAlpha) # if we include block by treatment interaction there is marginal significance for a treatment by evenness squared effect
even_sel2 <- lm(RelativeFitness ~ Block + TRT*even2 + even, FitAlpha) # if we include block by treatment interaction there is marginal significance for a treatment by evenness squared effect

anova(even_sel, even_sel2) # Based on RSS--second model/simpler one is a bit better
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: RelativeFitness ~ TRT * Block + TRT * even2 + even
    ## Model 2: RelativeFitness ~ Block + TRT * even2 + even
    ##   Res.Df    RSS Df Sum of Sq      F  Pr(>F)  
    ## 1     86 5.9104                              
    ## 2     89 6.5628 -3   -0.6524 3.1643 0.02855 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(even_sel) # NS
```

    ## Analysis of Variance Table
    ## 
    ## Response: RelativeFitness
    ##           Df Sum Sq Mean Sq F value    Pr(>F)    
    ## TRT        1 0.1156 0.11565  1.6828  0.198026    
    ## Block      3 8.0999 2.69998 39.2863 4.363e-16 ***
    ## even2      1 0.0516 0.05161  0.7509  0.388607    
    ## even       1 0.0364 0.03638  0.5293  0.468869    
    ## TRT:Block  3 0.3686 0.12288  1.7880  0.155487    
    ## TRT:even2  1 0.4811 0.48105  6.9996  0.009692 ** 
    ## Residuals 86 5.9104 0.06873                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
ggplot(FitAlpha) +
  geom_point(aes(RelativeFitness, even), alpha = 0.5, size = 5) +
  stat_smooth(method = "lm",formula = y ~poly(x, 2), aes(RelativeFitness,even)) +
  theme_classic() +
  ylab("Relative fitness") +
  xlab("Evenness") +
  Tx2
```

![](README_files/figure-gfm/cor%20on%20Evenness-1.png)<!-- -->

``` r
# Remove block effects from fitness
RootFitAlpha$FitResidBlk <- resid(lm(RelativeFitness ~ Block, RootFitAlpha))

RootFitAlpha$Comp <- sub(".*\\-", "", RootFitAlpha$Combos)
```

# ANCOVA

``` r
#RootFitAlpha<-merge(SizePurp, RootFitAlpha)

library(interactions)
```

    ## Warning: package 'interactions' was built under R version 4.0.5

``` r
library(jtools)
```

    ## Warning: package 'jtools' was built under R version 4.0.5

    ## 
    ## Attaching package: 'jtools'

    ## The following object is masked from 'package:Hmisc':
    ## 
    ##     %nin%

``` r
# Scale the microbial variables
RootFitAlpha$richScaled<-scale(RootFitAlpha$rich)
RootFitAlpha$SimScaled<-scale(RootFitAlpha$sim)
RootFitAlpha$InvSimScaled<-scale(RootFitAlpha$InvSimp)
RootFitAlpha$EvenScaled<-scale(RootFitAlpha$even)


# ANCOVA MUltivariate Linear Regressions


ANCOVA<-(lm(RelativeFitness~TRT + Block+ TRT:Block + PC1*TRT+PC2*TRT+PC3*TRT+PC4*TRT+richScaled*TRT+InvSimScaled*TRT+EvenScaled*TRT+PC1*Block + PC2*Block + PC3*Block + PC4*Block+richScaled*Block+InvSimScaled*Block+EvenScaled*Block, RootFitAlpha)) # Full model reported

# Step wise backward regression
library(MASS)

step<-stepAIC(ANCOVA,direction = "backward",trace=FALSE)
step$anova
```

    ## Stepwise Model Path 
    ## Analysis of Deviance Table
    ## 
    ## Initial Model:
    ## RelativeFitness ~ TRT + Block + TRT:Block + PC1 * TRT + PC2 * 
    ##     TRT + PC3 * TRT + PC4 * TRT + richScaled * TRT + InvSimScaled * 
    ##     TRT + EvenScaled * TRT + PC1 * Block + PC2 * Block + PC3 * 
    ##     Block + PC4 * Block + richScaled * Block + InvSimScaled * 
    ##     Block + EvenScaled * Block
    ## 
    ## Final Model:
    ## RelativeFitness ~ TRT + Block + PC1 + PC2 + PC3 + PC4 + richScaled + 
    ##     InvSimScaled + EvenScaled + TRT:richScaled + Block:PC1 + 
    ##     Block:PC2 + Block:PC4 + Block:richScaled + Block:InvSimScaled + 
    ##     Block:EvenScaled
    ## 
    ## 
    ##                 Step Df    Deviance Resid. Df Resid. Dev       AIC
    ## 1                                          54   3.297351 -242.0144
    ## 2 - TRT:InvSimScaled  1 0.002711687        55   3.300062 -243.9347
    ## 3          - TRT:PC4  1 0.003105322        56   3.303168 -245.8434
    ## 4          - TRT:PC2  1 0.013989404        57   3.317157 -247.4335
    ## 5          - TRT:PC3  1 0.056177982        58   3.373335 -247.8045
    ## 6        - TRT:Block  3 0.154840232        61   3.528175 -249.4512
    ## 7   - TRT:EvenScaled  1 0.033234583        62   3.561410 -250.5418
    ## 8        - Block:PC3  3 0.180266978        65   3.741677 -251.7522
    ## 9          - TRT:PC1  1 0.070962567        66   3.812639 -251.9298

``` r
model<-lm(RelativeFitness ~ TRT + Block + PC1 + PC2 + PC3 + PC4 + richScaled + 
    InvSimScaled + EvenScaled + TRT:Block + TRT:PC1 + TRT:PC3 + 
    TRT:richScaled + TRT:EvenScaled + Block:PC1 + Block:PC2 + 
    Block:PC3 + Block:PC4 + Block:richScaled + Block:InvSimScaled + 
    Block:EvenScaled, RootFitAlpha
)


# anova(model) # Type I sum of squres--sequence order matters, interactions not accounted

car::Anova(model, type="III") # Report the type three sums of squares~
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: RelativeFitness
    ##                    Sum Sq Df F value    Pr(>F)    
    ## (Intercept)        0.1813  1  3.1151 0.0829261 .  
    ## TRT                0.0157  1  0.2692 0.6058680    
    ## Block              0.7129  3  4.0836 0.0107317 *  
    ## PC1                0.0342  1  0.5877 0.4464823    
    ## PC2                0.0082  1  0.1404 0.7092886    
    ## PC3                0.0008  1  0.0138 0.9067979    
    ## PC4                0.0791  1  1.3598 0.2484366    
    ## richScaled         0.0173  1  0.2967 0.5880598    
    ## InvSimScaled       0.0004  1  0.0075 0.9313327    
    ## EvenScaled         0.0002  1  0.0030 0.9568681    
    ## TRT:Block          0.2049  3  1.1738 0.3277565    
    ## TRT:PC1            0.1326  1  2.2793 0.1366337    
    ## TRT:PC3            0.0562  1  0.9653 0.3300002    
    ## TRT:richScaled     0.2928  1  5.0320 0.0287860 *  
    ## TRT:EvenScaled     0.0777  1  1.3358 0.2525962    
    ## Block:PC1          0.8474  3  4.8536 0.0044709 ** 
    ## Block:PC2          0.2576  3  1.4753 0.2308904    
    ## Block:PC3          0.2561  3  1.4670 0.2331593    
    ## Block:PC4          1.1515  3  6.5954 0.0006649 ***
    ## Block:richScaled   0.5966  3  3.4173 0.0232423 *  
    ## Block:InvSimScaled 0.5618  3  3.2181 0.0293533 *  
    ## Block:EvenScaled   0.6109  3  3.4991 0.0211239 *  
    ## Residuals          3.3172 57                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = RelativeFitness ~ TRT + Block + PC1 + PC2 + PC3 + 
    ##     PC4 + richScaled + InvSimScaled + EvenScaled + TRT:Block + 
    ##     TRT:PC1 + TRT:PC3 + TRT:richScaled + TRT:EvenScaled + Block:PC1 + 
    ##     Block:PC2 + Block:PC3 + Block:PC4 + Block:richScaled + Block:InvSimScaled + 
    ##     Block:EvenScaled, data = RootFitAlpha)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.49654 -0.09763  0.01365  0.09599  0.47450 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)          0.52615    0.29811   1.765  0.08293 . 
    ## TRTInter             0.15810    0.30471   0.519  0.60587   
    ## Block2               0.04877    0.32605   0.150  0.88161   
    ## Block3               0.91575    0.35930   2.549  0.01353 * 
    ## Block4               0.66711    0.32098   2.078  0.04219 * 
    ## PC1                 -0.07673    0.10010  -0.767  0.44648   
    ## PC2                  0.03186    0.08504   0.375  0.70929   
    ## PC3                  0.02631    0.22369   0.118  0.90680   
    ## PC4                  0.26781    0.22967   1.166  0.24844   
    ## richScaled           0.26580    0.48793   0.545  0.58806   
    ## InvSimScaled         0.02068    0.23890   0.087  0.93133   
    ## EvenScaled           0.02131    0.39233   0.054  0.95687   
    ## TRTInter:Block2      0.10631    0.34635   0.307  0.75999   
    ## TRTInter:Block3     -0.21020    0.37365  -0.563  0.57594   
    ## TRTInter:Block4     -0.35654    0.33461  -1.066  0.29112   
    ## TRTInter:PC1         0.14505    0.09607   1.510  0.13663   
    ## TRTInter:PC3         0.15937    0.16221   0.983  0.33000   
    ## TRTInter:richScaled -0.57273    0.25532  -2.243  0.02879 * 
    ## TRTInter:EvenScaled -0.29567    0.25582  -1.156  0.25260   
    ## Block2:PC1          -0.13807    0.13485  -1.024  0.31020   
    ## Block3:PC1           0.17929    0.12345   1.452  0.15191   
    ## Block4:PC1          -0.16788    0.12715  -1.320  0.19201   
    ## Block2:PC2          -0.17160    0.11705  -1.466  0.14814   
    ## Block3:PC2           0.03278    0.13890   0.236  0.81429   
    ## Block4:PC2          -0.16215    0.11344  -1.429  0.15836   
    ## Block2:PC3          -0.21900    0.26023  -0.842  0.40354   
    ## Block3:PC3          -0.34845    0.26467  -1.317  0.19327   
    ## Block4:PC3          -0.47937    0.24121  -1.987  0.05170 . 
    ## Block2:PC4          -0.04392    0.25558  -0.172  0.86416   
    ## Block3:PC4          -0.81588    0.26961  -3.026  0.00371 **
    ## Block4:PC4          -0.39335    0.25880  -1.520  0.13406   
    ## Block2:richScaled    0.20016    0.62161   0.322  0.74863   
    ## Block3:richScaled   -0.36966    0.58524  -0.632  0.53015   
    ## Block4:richScaled    1.34560    0.56207   2.394  0.01998 * 
    ## Block2:InvSimScaled -0.01345    0.31255  -0.043  0.96581   
    ## Block3:InvSimScaled  0.26252    0.30760   0.853  0.39698   
    ## Block4:InvSimScaled -0.70144    0.34142  -2.054  0.04452 * 
    ## Block2:EvenScaled    0.18169    0.44342   0.410  0.68352   
    ## Block3:EvenScaled   -0.23577    0.45699  -0.516  0.60790   
    ## Block4:EvenScaled    1.14761    0.43686   2.627  0.01105 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2412 on 57 degrees of freedom
    ## Multiple R-squared:  0.7798, Adjusted R-squared:  0.6291 
    ## F-statistic: 5.176 on 39 and 57 DF,  p-value: 1.293e-08

``` r
# summary(ANCOVA)




summary(lm(RelativeFitness~rich, RootFitAlpha %>% filter(TRT == "Alone"))) # Alone
```

    ## 
    ## Call:
    ## lm(formula = RelativeFitness ~ rich, data = RootFitAlpha %>% 
    ##     filter(TRT == "Alone"))
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.5824 -0.3509 -0.1523  0.2511  1.0696 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)
    ## (Intercept) -1.313978   1.400070  -0.939    0.358
    ## rich         0.004895   0.003069   1.595    0.125
    ## 
    ## Residual standard error: 0.4722 on 22 degrees of freedom
    ## Multiple R-squared:  0.1037, Adjusted R-squared:  0.06291 
    ## F-statistic: 2.544 on 1 and 22 DF,  p-value: 0.125

``` r
summary(lm(RelativeFitness~rich, RootFitAlpha %>% filter(TRT != "Alone"))) # Competition
```

    ## 
    ## Call:
    ## lm(formula = RelativeFitness ~ rich, data = RootFitAlpha %>% 
    ##     filter(TRT != "Alone"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.55517 -0.29823 -0.00738  0.21166  1.16446 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept) -0.322197   0.569532  -0.566   0.5734  
    ## rich         0.002953   0.001275   2.317   0.0234 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3523 on 71 degrees of freedom
    ## Multiple R-squared:  0.07029,    Adjusted R-squared:  0.0572 
    ## F-statistic: 5.368 on 1 and 71 DF,  p-value: 0.0234

``` r
# Plot
EvenPlot1<-ggplot() +
  geom_point(data=RootFitAlpha %>% filter(TRT == "Alone"), aes(EvenScaled,RelativeFitness), size = 3, alpha=0.5) +
#    geom_smooth(data=RootFitAlpha %>% filter(TRT == "Alone"), aes(even,RelativeFitness),method = "lm", fullrange = TRUE, se=FALSE) +
  geom_abline(slope=-0.47, intercept = 1.18,color="blue", size = 2) +
  theme_classic() +
  xlab("Evenness") +
  Tx+
  theme(axis.text.x = element_text(angle=45), axis.text = element_text(color="black", size=15,vjust = 0.5,hjust=1)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))


EvenPlot2<-ggplot() +
  geom_point(data=RootFitAlpha %>% filter(TRT=="Inter"), aes(EvenScaled,RelativeFitness), size = 3, alpha=0.5) +
#    geom_smooth(data=RootFitAlpha %>% filter(TRT=="Inter"), aes(EvenScaled,RelativeFitness),method = "lm", fullrange = TRUE, se=FALSE) +
    geom_abline(slope=-0.28, intercept = 1.03,color="blue", size = 2) +

  theme_classic() +
  xlab("Evenness") +
  Tx+
  theme(axis.text.x = element_text(angle=45), axis.text = element_text(color="black", size=15,vjust = 0.5,hjust=1)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))

ggplot() +
  geom_point(data=RootFitAlpha, aes(EvenScaled,RelativeFitness,color=TRT), size = 3, alpha=0.5) +
#    geom_smooth(data=RootFitAlpha %>% filter(TRT=="Inter"), aes(EvenScaled,RelativeFitness),method = "lm", fullrange = TRUE, se=FALSE) +
    geom_abline(slope=-0.47, intercept = 1.18,color="red", size = 2, Linetype="dashed") +

    geom_abline(slope=-0.28, intercept = 1.03,color="black", size = 2) +

  theme_classic() +
  xlab("Evenness") +
  Tx2+
  theme(axis.text.x = element_text(angle=45), axis.text = element_text(color="black", size=15,vjust = 0.5,hjust=1)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) + scale_color_manual(values = c("red", "black"), "Treatment", labels = c("Alone", "Competition"))
```

    ## Warning: Ignoring unknown parameters: Linetype

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
#ggarrange(P2.rich,P2.even,P4.Sim,P4.simIn,nrow=2,ncol=2)
cowplot::plot_grid(EvenPlot1,EvenPlot2, align = "hv", ncol = 2, Labels =c("A", "B"), label_size =22 ,hjust = -9,vjust=2, Label_x = -0.16)
```

    ## Warning in as_grob.default(plot): Cannot convert object of class character into
    ## a grob.

    ## Warning in as_grob.default(plot): Cannot convert object of class numeric into a
    ## grob.

    ## Warning: Graphs cannot be vertically aligned unless the axis parameter is set.
    ## Placing graphs unaligned.

    ## Warning: Graphs cannot be horizontally aligned unless the axis parameter is set.
    ## Placing graphs unaligned.

![](README_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
Block1<-ggplot() +
  geom_point(data=RootFitAlpha %>% filter(Block=="1"), aes(even,RelativeFitness), size = 3, alpha=0.5) +
    geom_smooth(data=RootFitAlpha %>% filter(Block=="1"), aes(even,RelativeFitness), se=FALSE,method = "lm", fullrange = TRUE) +
    theme_classic() +
  xlab("Eveness") +
  Tx+
  theme(axis.text.x = element_text(angle=45), axis.text = element_text(color="black", size=15,vjust = 0.5,hjust=1))
   
Block2<-ggplot() +
  geom_point(data=RootFitAlpha %>% filter(Block=="2"), aes(even,RelativeFitness), size = 3, alpha=0.5) +
    geom_smooth(data=RootFitAlpha %>% filter(Block=="2"), aes(even,RelativeFitness), se=FALSE,method = "lm", fullrange = TRUE) +
    theme_classic() +
  xlab("Eveness") +
  Tx+
  theme(axis.text.x = element_text(angle=45), axis.text = element_text(color="black", size=15,vjust = 0.5,hjust=1))


Block3<-ggplot() +
  geom_point(data=RootFitAlpha %>% filter(Block=="3"), aes(even,RelativeFitness), size = 3, alpha=0.5) +
    geom_smooth(data=RootFitAlpha %>% filter(Block=="3"), aes(even,RelativeFitness), se=FALSE,method = "lm", fullrange = TRUE) +
    theme_classic() +
  xlab("Eveness") +
  Tx+
  theme(axis.text.x = element_text(angle=45), axis.text = element_text(color="black", size=15,vjust = 0.5,hjust=1))


Block4<-ggplot() +
  geom_point(data=RootFitAlpha %>% filter(Block=="4"), aes(even,RelativeFitness), size = 3, alpha=0.5) +
    geom_smooth(data=RootFitAlpha %>% filter(Block=="4"), aes(even,RelativeFitness),method = "lm", fullrange = TRUE) +
    theme_classic() +
  xlab("Eveness") +
  Tx+
  theme(axis.text.x = element_text(angle=45), axis.text = element_text(color="black", size=15,vjust = 0.5,hjust=1))
```

# MANTEL (Table 4)

``` r
#  Here we use the Family mean values of root traits---this would indicate evidence for 'phenotypic selection' 

  # Load matrix Bray distances of OTU abundance
  OTU_table = t(otu_table(physeq1)) # Write out OTU table
  

  # Pull out Root traits of interest and save these to the Sampled data frame
  SampledFit = merge(FitAveraged, sampledf, by = c("TRT", "ML", "Block"))


# Pull out Root traits of interest and save these to the Sampled data frame
SampledRoots = merge(RootAlphaObs, sampledf)
SampledRootsIp = subset(SampledRoots, Species == "Ip")

OTU_table = OTU_table[row.names(OTU_table) %in% SampledRootsIp$Sample_ID,]

    OTU = OTU_table
    Roots = SampledRootsIp

  # Calculate root architecture distances with euclidean distance 

    PC2 = SampledRootsIp$PC2 # isolate PC2, i.e. root architecture
    PC2.dist = dist(PC2)

  # Calculate Bray distance matrix for OTU table

    Bray = vegdist(OTU, method = "bray")

    
# OTU Bray vs Root architecture

    OTU_pc2 = mantel(Bray, PC2.dist, method = "spearman", permutations = 9999, na.rm = TRUE)

    # Marginally significant--very low r value
    # Mantel statistic r: 0.06836  
    # Significance: 0.07
    
# Repeat for root topology, size and morphology

  PC1 = SampledRootsIp$PC1 # isolate 
  PC1.dist = dist(PC1)
  PC3 = SampledRootsIp$PC3 # isolate 
  PC3.dist=dist(PC3)
  PC4 = SampledRootsIp$PC1 # isolate 
  PC4.dist = dist(PC4)

# Examine overall distances in root system; i.e. use all PCs in distance calculation

    PCall = SampledRootsIp[grep("PC",names(SampledRootsIp))]
    PC.dist = dist(PCall)

    OTU_pc = mantel(Bray, PC.dist, method = "spearman", permutations = 9999, na.rm = TRUE)

###             examine mantel test using observed values of roots ###

OTU_pc1 = mantel(Bray, PC1.dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# OTU_pc2 = mantel(Bray,PC2.dist, method = "spearman", permutations = 9999, na.rm = TRUE)
OTU_pc3 = mantel(Bray, PC3.dist, method = "spearman", permutations = 9999, na.rm = TRUE)
OTU_pc4 = mantel(Bray, PC4.dist, method = "spearman", permutations = 9999, na.rm = TRUE)
```

PC1 Mantel

``` r
OTU_pc1
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = Bray, ydis = PC1.dist, method = "spearman", permutations = 9999,      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.04189 
    ##       Significance: 0.7636 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0729 0.0975 0.1179 0.1435 
    ## Permutation: free
    ## Number of permutations: 9999

PC2 Mantel

``` r
OTU_pc2
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = Bray, ydis = PC2.dist, method = "spearman", permutations = 9999,      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.06836 
    ##       Significance: 0.071 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0595 0.0777 0.0937 0.1109 
    ## Permutation: free
    ## Number of permutations: 9999

PC3 Mantel

``` r
OTU_pc3
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = Bray, ydis = PC3.dist, method = "spearman", permutations = 9999,      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.07133 
    ##       Significance: 0.1325 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0829 0.1064 0.1274 0.1491 
    ## Permutation: free
    ## Number of permutations: 9999

PC4 Mantel

``` r
OTU_pc4
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = Bray, ydis = PC4.dist, method = "spearman", permutations = 9999,      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.04189 
    ##       Significance: 0.7618 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0735 0.0954 0.1155 0.1411 
    ## Permutation: free
    ## Number of permutations: 9999
