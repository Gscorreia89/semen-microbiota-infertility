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
df <- as.data.frame(sample_data(phyloseqDataset)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(phyloseqDataset)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
decontamPlotA <- ggplot(data=df, aes(x=Index, y=LibrarySize, color=SampleType)) + geom_point()

ggsave('./DecontamPlotA.png', decontamPlotA)

set.seed(100)
decontamPlotB <- plot_frequency(subset_samples(phyloseqDataset, !is.na(`Total.16s.rRNA.copy.per.mL`)), taxa_names(subset_samples(phyloseqDataset, !is.na(`Total.16s.rRNA.copy.per.mL`)))[sample(which(!contamdf.freq$contaminant), 81)], conc="Total.16s.rRNA.copy.per.mL") +                   
  xlab("Total 16s rRNA")

ggsave('./DecontamPlotB.png', decontamPlotB)

taxaToFilter <- c(rownames(contamdf.prev[contamdf.prev$contaminant == TRUE, ]), rownames(contamdf.freq[contamdf.freq$contaminant ==TRUE, ]))

taxaSelection <- setdiff(taxa_names(phyloseqDataset), taxaToFilter)

phyloDFiltered <- prune_taxa(taxaSelection, phyloseqDataset)

otu_table(phyloDFiltered) %>%   as("matrix") %>%
  make_biom() %>%
  write_biom("./3361-ASV-Outputs-V1MixedV2-SILVA/table-decontam-filtered.biom")
