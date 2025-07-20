# Purpose: Prepare and filter phyloseq object

# Load packages

library(phyloseq)
library(ggplot2)
library(decontam)
library(vegan)
library(dplyr)
library(stringr)
library(data.table)

# Load phyloseq object
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])}
ITS <-  loadRData("data/phyloseq_unitev10.Rdata") 
sample_data(ITS)

# 1.0 Add climate metadata -----

climatevars <- read.csv('data/climatevars.csv') # Load data
row.names(climatevars) <- climatevars$Sample 
sample_data <- data.frame(sample_data(ITS))

climate_sample <- merge(sample_data, climatevars, by = "row.names", all.x = TRUE) # Merge climate variables and sample data
col_names_indices <- data.frame(Index = seq_along(colnames(climate_sample)), Column = colnames(climate_sample))
print(col_names_indices)

climate_sample <- climate_sample[-c(1:3,6,8,9,14:20)] # Remove unwanted columns by index
row.names(climate_sample) <- climate_sample$Sample 
sample_data(ITS) <- sample_data(climate_sample) # Merge sample data back into the phyloseq object
sample_data(ITS) # Confirm successful integration


# 2.0 Remove contaminants -----


# Identify negative control samples 
  # 'Pool' samples are negative extraction controls
  # 'AGRF' samples are negative library prep & sequencing controls

df <- data.frame(sample_data(ITS))
df <- df %>% mutate(Sample_or_Control = case_when(str_detect(textbefore, "Pool") ~ 'Control', 
                                                  str_detect(textbefore, "AGRF") ~ 'Control',
                                                  is.na(textbefore) ~ NA_character_, TRUE ~ "Sample"))

df$LibrarySize <- sample_sums(ITS)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

# Visualise library size 
Control <- filter(df, Sample_or_Control == "Control")
Sample <- filter(df, Sample_or_Control == "Sample")

ggplot(data=Control, aes(x=Index, y=LibrarySize, color= Sample_or_Control))+ geom_point(data=Sample, aes(x=Index, y=LibrarySize, color= Sample_or_Control)) + geom_point(data=Control, aes(x=Index, y=LibrarySize, color= Sample_or_Control))

sample_data(ITS) <- sample_data(df)

# Identify contaminants by prevalence 
sample_data(ITS)$is.neg = sample_data(ITS)$Sample_or_Control == "Control"
contamdf.prev = isContaminant(ITS, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev$contaminant)

# Visualize prevalence of contaminants in dataset
ITS.pa = transform_sample_counts(ITS, function(abund) 1*(abund>0))
ITS.pa.neg = prune_samples(sample_data(ITS.pa)$Sample_or_Control == "Control", ITS.pa)
ITS.pa.pos = prune_samples(sample_data(ITS.pa)$Sample_or_Control == "Sample", ITS.pa)

df.pa = data.frame(pa.pos=taxa_sums(ITS.pa.pos), pa.neg=taxa_sums(ITS.pa.neg),
                   contaminant=contamdf.prev$contaminant)

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Extract the taxonomic classification of contaminants
list = row.names(contamdf.prev)[which(contamdf.prev$contaminant=="TRUE")]
contam = prune_taxa(list, ITS)

# Prune contaminant taxa from the phyloseq tax_table
ITS_decontam = prune_taxa(!contamdf.prev$contaminant, ITS)
ITS_decontam


# 3.0 Filter and transform data -----

## 3.1 Assign guild and filter -----

# Change names of tax table so that they correspond to FungalTraits
colnames(tax_table(ITS_decontam)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") 

# Function to assign traits based on taxonomy 
source("R_functions/assign_guild.R")
guild_table = assign_guild(object = ITS_decontam, database="FungalTraits")

# Subset samples collected at Bellbird  
ITS_decontam <- subset_samples(ITS_decontam, Location %in% "Bellbird")

# Filter OTUs
source("R_functions/filter_OTU_per_sample.R")
ITS_filtered = filter_OTU_per_sample(ITS_decontam, threshold = 0.01) # Remove OTUs with counts lower than 0.01% in relative abundance per sample
ITS_filtered = prune_taxa(taxa_sums(ITS_filtered) > 0, ITS_filtered) # Remove OTUs with zero counts
ITS_filtered

# Visualise read depth showing distribution of reads per sample
readcount = data.table(as(sample_data(ITS_filtered), "data.frame"),
                       total_reads = sample_sums(ITS_filtered), 
                       keep.rownames = TRUE)
setnames(readcount, "rn", "SampleID", skip_absent=TRUE)

ggplot(readcount, aes(total_reads)) + geom_histogram() + ggtitle("Sequencing Depth")

# Remove samples with read count lower than 1000
ITS_filtered<-prune_samples(sample_sums(ITS_filtered)>=1000, ITS_filtered)
sample_names(ITS_filtered)
ITS_filtered

# Remove control samples if still present
ITS_filtered = subset_samples(ITS_filtered, Sample_or_Control == "Sample") # Remove the remaining negative controls
ITS_filtered = subset_samples(ITS_filtered, textbefore != "MC") # Remove the positive control (mock community)

## 3.2 Transform data for further analysis ----

ITSrel = transform_sample_counts(ITS_filtered, function(x) x / sum(x) )
sample_sums(ITS_filtered)



# Relative abundance as count: 
ITSrel.count  = transform_sample_counts(ITS_filtered, function(x) 100000 * x/sum(x)) # Scale relative abundance to 100,000 to account for differences in sequence depth
otu_table(ITSrel.count) = ceiling(otu_table(ITSrel.count, "matrix")) # Transform to next integer 


# Presence/absence dataset
binary_table=decostand(otu_table(ITSrel),method="pa")
ITS.binary = phyloseq(otu_table(binary_table, taxa_are_rows=TRUE), 
                      sample_data(ITSrel), 
                      tax_table(tax_table(ITSrel)))
sample_sums(ITS.binary)



# Remove OTUs that are present in only 1 sample
ITS.binary2 = prune_taxa(taxa_sums(ITS.binary) > 1, ITS.binary)
list.ITS.binary2=taxa_names(ITS.binary2)
ITSrel2=prune_taxa(list.ITS.binary2, ITSrel)
ITSrel.count2=prune_taxa(list.ITS.binary2, ITSrel.count)

# Final datasets
ITSrel.count
ITSrel.count2
ITS.binary
ITS.binary2
guild_table

sample_data(ITSrel.count)

