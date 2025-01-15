library(decontam)
library(qiime2R)
library(readr)
library(phyloseq)
library(tidyverse)
library(biomformat)
library(dplyr)


countsTable <- qiime2R::read_qza('./Run6631-PTrimV1MixedV2-table-dada2.qza')

taxonomyTable <- qiime2R::read_qza('./Run6631-SILVA-V1MixedV2-Taxonomy.qza') 

countsMatrix <- countsTable$data
taxonomyData <- taxonomyTable$data
decontamMetadata <- readr::read_csv('DecontamMetadata.csv')


tax_tab <- taxonomyTable$data %>% 
  as.data.frame() %>%
  separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  select(-Confidence)

phyloseqDataset <- phyloseq::phyloseq(sample_data(column_to_rownames(decontamMetadata, "sample-id")), 
                                      otu_table(as.matrix(countsMatrix), taxa_are_rows = TRUE),
                                      tax_table(as.matrix(tax_tab)))

# Frequency filter
contamdf.freq <- isContaminant(subset_samples(phyloseqDataset, !is.na(`Total.16s.rRNA.copy.per.mL`)), method="frequency", conc="Total.16s.rRNA.copy.per.mL")

sample_data(phyloseqDataset)$is.neg <- sample_data(phyloseqDataset)$`SampleType` == "Control"
contamdf.prev <- isContaminant(phyloseqDataset, method="prevalence", neg="is.neg")

prevContaminants <- tax_table(phyloseqDataset)[contamdf.prev$contaminant == TRUE, ]
freqContaminants <- tax_table(phyloseqDataset)[contamdf.freq$contaminant == TRUE, ]

write.csv(prevContaminants, 'decontam-Filtered-Prevalence.csv')
write.csv(freqContaminants, 'decontam-Filtered-Frequence.csv')

## Plots 
library(ggplot2)

# Prevalence filter plot
ps.pa <- transform_sample_counts(phyloseqDataset, function(abund) 1*(abund > 0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$SampleType == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$SampleType == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
decontamPlotA <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") + theme_light()

ggsave('./DecontamPlot_Prevalence.png', decontamPlotB)


# Frequency filter plot
set.seed(100)
decontamPlotB <- plot_frequency(subset_samples(phyloseqDataset, !is.na(`Total.16s.rRNA.copy.per.mL`)), 
                                taxa_names(subset_samples(phyloseqDataset, !is.na(`Total.16s.rRNA.copy.per.mL`)))[sample(which(contamdf.freq$contaminant), 10)], 
                                conc="Total.16s.rRNA.copy.per.mL") +                   
xlab("Total 16s rRNA") + theme_light()

ggsave('./DecontamPlot_Frequency.png', decontamPlotB)

df <- as.data.frame(sample_data(phyloseqDataset)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(phyloseqDataset)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
decontamPlotA <- ggplot(data=df, aes(x=Index, y=LibrarySize, color=SampleType)) + geom_point()
ggsave('./DecontamPlotA.png', decontamPlotA)


taxaToFilter <- c(rownames(contamdf.prev[contamdf.prev$contaminant == TRUE, ]), rownames(contamdf.freq[contamdf.freq$contaminant ==TRUE, ]))

taxaSelection <- setdiff(taxa_names(phyloseqDataset), taxaToFilter)

phyloDFiltered <- prune_taxa(taxaSelection, phyloseqDataset)

otu_table(phyloDFiltered) %>%   as("matrix") %>%
  make_biom() %>%
  write_biom("./3361-ASV-Outputs-V1MixedV2-SILVA/table-decontam-filtered.biom")

