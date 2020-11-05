# hecatomb-eukVir-readStatsTable.R

library("phyloseq")
library("data.table")
library("tidyverse")

#----- Description Files -----#

uniprot <- read.delim("./data/descriptions/uniprot.description.db", 
                      header = FALSE, col.names = c("target", "protein"),
                      stringsAsFactors = FALSE)

baltimore <- read.delim("./data/descriptions/2020_07_27_Viral_classification_table_ICTV2019.txt",
                        header = TRUE, stringsAsFactors = FALSE,
                        blank.lines.skip = TRUE,
                        col.names = c("Family", "Baltimore", "Baltimore_Group"))

#----- Eukaryotic Taxonomy Table -----#

# Load eukaryotic virus taxonomy table and remove any Bacteria sequences.
#   Replace any spaces and underscores which appear within a column with
#   hyphens "-" for easier parsing when summing counts across taxa of the
#   same name. Equivalent to Scott's viral_table.

eukTaxTable <- fread(file = "./data/hecatomb_out/viruses_tax_table.tsv", 
                     header = TRUE, stringsAsFactors = FALSE,
                     colClasses = "character") %>% 
  filter(Kingdom != "Bacteria") %>% 
  mutate(across(where(is.character), ~ str_replace_all(., "\\s|_", "-"))) %>% 
  as.data.frame()

#----- Full Count Table (Phage & Eukaryotic) -----#

# Contains all sequences - generated from all samples;
#   Remove sequence column and replace NAs with 0.
#   Equivalent to Scott's seqtable.

virusCountTable <- fread(file = "./data/hecatomb_out/seqtable.all", 
                         header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE,
                         colClasses = c("id" = "character")) %>% 
  as.data.frame() %>% 
  select(-sequence) %>% 
  mutate(across(where(is.numeric),~ replace(., is.na(.), 0)))

#----- Alignment Statistics -----#

# aa checked alignment:
aaAln <- fread(file = "./data/hecatomb_out/aa.aln.m8", header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, 
               colClasses = c("query" = "character")) %>% 
  rename(id = query) %>% 
  mutate(query_type = "aa") %>% 
  as.data.frame()

# nt checked alignment:
ntAln <- fread(file = "./data/hecatomb_out/nt.aln.m8", header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, 
               colClasses = c("query" = "character")) %>%  
  rename(id = query) %>% 
  mutate(query_type = "nt") %>% 
  as.data.frame()

# bind aa and alignment tables by row:
allAln <- rbind(aaAln, ntAln)

#----- Standardize Counts in virusCountTable to Per-Sample Library Size -----#

librarySize <- enframe(colSums(virusCountTable %>% select(-"id")),
                       name = "Sample", value = "library_size") %>% 
  as.data.frame()

#----- Generate Table of Individual Read Alignment Stats -----#

# Merge eukaryotic virus taxonomy table with alignment data and count table
mergedData <- Reduce(f = function(df1, df2) 
{base::merge(x = df1, y = df2, by = "id", all.x = TRUE)},
x = list(eukTaxTable, allAln, virusCountTable))

sampleData <- readRDS("./data/metadata/combinedMetadata.RDS")
sampleNames <- as.character(pull(sampleData, sample_names))

mergedPivot <- mergedData %>% 
  pivot_longer(cols = all_of(sampleNames), 
               names_to = "Sample", values_to = "Abundance") %>% 
  merge(librarySize, by = "Sample", all = TRUE) %>% 
  merge(baltimore, by = "Family", all = TRUE) %>% 
  merge(uniprot, by = "target", all = TRUE) %>% 
  filter(!is.na(Sample)) %>% 
  mutate(proportional_abundance = Abundance/library_size,
         scaled_abundance_min = (min(librarySize$library_size)*proportional_abundance),
         scaled_abundance_mean = (mean(librarySize$library_size)*proportional_abundance),
         scaled_abundance_median = (median(librarySize$library_size)*proportional_abundance)) %>% 
  mutate(alignment_length_adjusted = base::ifelse(query_type == "aa",
                                                  yes = alignment_length*3,
                                                  no = alignment_length*1)) %>% 
  merge(sampleData, by.x = "Sample", by.y = "sample_names", all = TRUE) 

mergedPivotFilt <- mergedPivot %>% 
  filter(Abundance > 0) %>% 
  droplevels # 30,206 x 50
