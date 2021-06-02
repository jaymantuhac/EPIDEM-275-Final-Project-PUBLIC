#Setup
library(tidyverse)
library(biomaRt)
library(readr)
library(gwascat)
library(stringr)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(kableExtra)

#Checking out TxDb package
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb

transcripts(txdb)

tx.by.gene <- transcriptsBy(txdb, "gene")
tx.by.gene

#Checking out org.Hw.eg.db package #Lactose intolerance SNP found in MCM6 gene
cols(org.Hs.eg.db)
select(org.Hs.eg.db, keys = "APOE", cols = c("ENTREZID", "SYMBOL", "GENENAME"),
       keytype = "SYMBOL")

#Load Reference Spreadsheets
pgp_ref_1 <- read_csv("PGPParticipantSurvey-20181010220019.csv")
pgp_ref_2 <- read_csv("PGPTrait&DiseaseSurvey2012_DigestiveSystem-20181010214607.csv")
pgp_ref_3 <- read_csv("PGPTrait&DiseaseSurvey2012_Endocrine,Metabolic,Nutritional,AndImmunity-20181010220044.csv")


#Load Subject Data
#Note: This section contains code for loading individual-level genome data into R. Due to the sensitive nature of this data, I have censored this section's code for display on this public repo.

#Convert "chrom" response datatype to factors
subject_1$chrom <- ordered(subject_1$chrom, levels = c(seq(1, 22), "X", "Y", "MT"))
subject_2$chrom <- ordered(subject_2$chrom, levels = c(seq(1, 22), "X", "Y", "MT"))
subject_3$chrom <- ordered(subject_3$chrom, levels = c(seq(1, 22), "X", "Y", "MT"))

#Update gwas reference data
updated_gwas_data <- as.data.frame(makeCurrentGwascat())

#Subject 1 ---------------------
#Subject 1 EDA
subject_1_EDA <- subject_1 %>% 
  ggplot(aes(chrom)) + 
  geom_bar() + 
  labs(x = "Chromosome", y = "Number of SNPs") + 
  ggtitle("Subject 1 SNP Distribution")

subject_1_joined <- inner_join(subject_1, updated_gwas_data, 
                               by = c('rsid' = "SNPS"))

subject_1_joined$risk_allele_clean <- str_sub(subject_1_joined$STRONGEST.SNP.RISK.ALLELE, -1)
subject_1_joined$my_allele_1 <- str_sub(subject_1_joined$genotype, 1, 1)
subject_1_joined$my_allele_2 <- str_sub(subject_1_joined$genotype, 2, 2)
subject_1_joined$have_risk_allele_count <- if_else(subject_1_joined$my_allele_1 == 
                                                     subject_1_joined$risk_allele_clean, 1, 0) +
  if_else(subject_1_joined$my_allele_2 == subject_1_joined$risk_allele_clean, 1, 0)

subject_1_risk_data <- dplyr::select(subject_1_joined, 
                                     rsid, 
                                     have_risk_allele_count,
                                     DISEASE.TRAIT, 
                                     risk_allele = risk_allele_clean,
                                     your_genotype = genotype,
                                     RISK.ALLELE.FREQUENCY,
                                     MAPPED_TRAIT,
                                     REPORTED.GENE.S.)
subject_1_risk_data

#Subject 1 Count of Disease Traits
subject_1_trait_count_mapped <- subject_1_risk_data %>% 
  group_by(MAPPED_TRAIT) %>% 
  summarise(risk_count = n())

#Subject 1 Count of Disease Traits for Risk Allele = 2
subject_1_trait_count_mapped2 <- subject_1_risk_data %>% 
  group_by(MAPPED_TRAIT, have_risk_allele_count) %>% 
  filter(have_risk_allele_count == 2) %>% 
  summarise(risk_count_2 = n()) %>% 
  subset(select = -c(have_risk_allele_count))

#Subject 1 Mutating Risk Counts of Allele = 1 and allele = 2
subject_1_trait_count_mapped3 <- 
  merge(subject_1_trait_count_mapped, subject_1_trait_count_mapped2,
        by = c("MAPPED_TRAIT", "MAPPED_TRAIT"), all = TRUE)

subject_1_trait_count_mapped3[is.na(subject_1_trait_count_mapped3)] <- 0

#Risk Count for Allele 1 Equation
risk_1 <- (subject_1_trait_count_mapped3$risk_count - subject_1_trait_count_mapped3$risk_count_2)

#Adding risk_1 column in mapped4 table
subject_1_trait_count_mapped4 <- dplyr::select(subject_1_trait_count_mapped3,
                                               MAPPED_TRAIT,
                                               risk_count,
                                               risk_count_2) %>% 
  mutate(risk_1)

#Adding Overall Risk Column and Equation
subject_1_overall_risk <- (subject_1_trait_count_mapped3$risk_count_2 / ((subject_1_trait_count_mapped3$risk_count_2 + subject_1_trait_count_mapped4$risk_1*2))) + 
  (subject_1_trait_count_mapped4$risk_1/((subject_1_trait_count_mapped3$risk_count_2+subject_1_trait_count_mapped4$risk_1*2)))

subject_1_overall_risk_count <- dplyr::select(subject_1_trait_count_mapped4,
                                              MAPPED_TRAIT,
                                              risk_count,
                                              risk_1,
                                              risk_count_2) %>% 
  mutate(subject_1_overall_risk * 100)

subject_1_overall_risk_count[is.na(subject_1_overall_risk_count)] <- 0

#Get Count of Disease Traits by DISEASE.TRAIT
subject_1_trait_count_disease <- subject_1_risk_data %>% 
  group_by(DISEASE.TRAIT) %>% 
  summarise(risk_count = n())

#Count of Disease Traits for Risk Allele = 2
subject_1_trait_count_disease2 <- subject_1_risk_data %>% 
  group_by(DISEASE.TRAIT, have_risk_allele_count) %>% 
  filter(have_risk_allele_count == 2) %>% 
  summarise(risk_count_2 = n()) %>% 
  subset(select = -c(have_risk_allele_count))

#Mutating Risk Counts for Allele = 1 and Allele = 2
subject_1_trait_count_disease3 <- 
  merge(subject_1_trait_count_disease, subject_1_trait_count_disease2,
        by = c("DISEASE.TRAIT", "DISEASE.TRAIT"), all = TRUE)
subject_1_trait_count_disease3[is.na(subject_1_trait_count_disease3)] <- 0 

#Risk Count for Allele 1 equation
risk_d_1 <- (subject_1_trait_count_disease3$risk_count - subject_1_trait_count_disease3$risk_count_2)

#Adding risk 1 column in 
subject_1_trait_count_disease4 <- dplyr::select(subject_1_trait_count_disease3,
                                                DISEASE.TRAIT,
                                                risk_count,
                                                risk_count_2) %>% 
  mutate(risk_d_1)

#Adding Overall Risk Column and Equation
overall_risk_d <- (subject_1_trait_count_disease4$risk_count_2 / ((subject_1_trait_count_disease4$risk_d_1 * 2) + subject_1_trait_count_disease4$risk_count_2)) + 
  (subject_1_trait_count_disease4$risk_d_1 / ((subject_1_trait_count_disease4$risk_d_1 * 2) + subject_1_trait_count_disease4$risk_count_2))

subject_1_overall_risk_count_d <- dplyr::select(subject_1_trait_count_disease4,
                                                DISEASE.TRAIT,
                                                risk_count,
                                                risk_1 = "risk_d_1",
                                                risk_count_2) %>%
  mutate(overall_risk_d*100) 

overall_risk_count_d[is.na(overall_risk_count_d)] <- 0

#Table of Disease Traits with over 100 risk alleles for MAPPED_TRAIT
subject_1_table1 <- subject_1_overall_risk_count %>% 
  filter(risk_count > 99) %>% 
  arrange(desc(risk_count)) %>% 
  subset(select = -c(risk_count, risk_1, risk_count_2))

names(subject_1_table1)[1] <- 'Mapped Trait'
names(subject_1_table1)[2] <- 'Overall Risk'

subject_1_table1 <- subject_1_table1 %>% 
  kable() %>% 
  kable_paper("hover", full_width = FALSE)

subject_1_table1

#Table of Disease Traits with over 100 risk alleles for DISEASE.TRAIT 
subject_1_table2 <- subject_1_overall_risk_count_d %>% 
  filter(risk_count > 99) %>% 
  arrange(desc(`overall_risk_d * 100`)) %>% 
  subset(select = -c(risk_count, risk_1, risk_count_2))

names(subject_1_table2)[1] <- 'Disease Trait'
names(subject_1_table2)[2] <- 'Overall Risk'

subject_1_table2 <- subject_1_table2 %>% 
  kable() %>% 
  kable_paper("hover", full_width = FALSE)

subject_1_table2

#Subject 2----------------
#Subject 2 EDA
subject_2_EDA <- subject_2 %>% 
  ggplot(aes(chrom)) + 
  geom_bar() + 
  labs(x = "Chromosome", y = "Number of SNPs") + 
  ggtitle("Subject 2 SNP Distribution")

#Subject 2 EDA
subject_2 %>% 
  ggplot(aes(chrom)) + 
  geom_bar()

#Subject 2 Risk Analysis
subject_2_joined <- inner_join(subject_2, updated_gwas_data, 
                               by = c('rsid' = "SNPS"))

subject_2_joined$risk_allele_clean <- str_sub(subject_2_joined$STRONGEST.SNP.RISK.ALLELE, -1)
subject_2_joined$my_allele_1 <- str_sub(subject_2_joined$genotype, 1, 1)
subject_2_joined$my_allele_2 <- str_sub(subject_2_joined$genotype, 2, 2)
subject_2_joined$have_risk_allele_count <- if_else(subject_2_joined$my_allele_1 == 
                                                     subject_2_joined$risk_allele_clean, 1, 0) +
  if_else(subject_2_joined$my_allele_2 == subject_2_joined$risk_allele_clean, 1, 0)

subject_2_risk_data <- dplyr::select(subject_2_joined, 
                                     rsid, 
                                     have_risk_allele_count,
                                     DISEASE.TRAIT, 
                                     risk_allele = risk_allele_clean,
                                     your_genotype = genotype,
                                     RISK.ALLELE.FREQUENCY,
                                     MAPPED_TRAIT,
                                     REPORTED.GENE.S.)
subject_2_risk_data

#Subject 2 Count of Disease Traits
subject_2_trait_count_mapped <- subject_2_risk_data %>% 
  group_by(MAPPED_TRAIT) %>% 
  summarise(risk_count = n())

#Subject 2 Count of Disease Traits for Risk Allele = 2
subject_2_trait_count_mapped2 <- subject_2_risk_data %>% 
  group_by(MAPPED_TRAIT, have_risk_allele_count) %>% 
  filter(have_risk_allele_count == 2) %>% 
  summarise(risk_count_2 = n()) %>% 
  subset(select = -c(have_risk_allele_count))

#Subject 2 Mutating Risk Counts of Allele = 1 and allele = 2
subject_2_trait_count_mapped3 <- 
  merge(subject_2_trait_count_mapped, subject_2_trait_count_mapped2,
        by = c("MAPPED_TRAIT", "MAPPED_TRAIT"), all = TRUE)

subject_2_trait_count_mapped3[is.na(subject_2_trait_count_mapped3)] <- 0

#Risk Count for Allele 1 Equation
subject_2_risk_1 <- (subject_2_trait_count_mapped3$risk_count - subject_2_trait_count_mapped3$risk_count_2)

#Adding risk_1 column in mapped4 table
subject_2_trait_count_mapped4 <- dplyr::select(subject_1_trait_count_mapped3,
                                               MAPPED_TRAIT,
                                               risk_count,
                                               risk_count_2) %>% 
  mutate(subject_2_risk_1)

#Adding Overall Risk Column and Equation
subject_2_overall_risk <- (subject_2_trait_count_mapped4$risk_count_2 / ((subject_2_trait_count_mapped4$subject_2_risk_1 * 2) + subject_2_trait_count_mapped4$risk_count_2)) + 
  (subject_2_trait_count_mapped4$subject_2_risk_1 / ((subject_2_trait_count_mapped4$subject_2_risk_1 * 2) + subject_2_trait_count_mapped4$risk_count_2))

subject_2_overall_risk_count <- dplyr::select(subject_2_trait_count_mapped4,
                                              MAPPED_TRAIT,
                                              risk_count,
                                              subject_2_risk_1,
                                              risk_count_2) %>% 
  mutate(subject_2_overall_risk * 100)

subject_2_overall_risk_count[is.na(subject_2_overall_risk_count)] <- 0

#Get Count of Disease Traits by DISEASE.TRAIT
subject_2_trait_count_disease <- subject_2_risk_data %>% 
  group_by(DISEASE.TRAIT) %>% 
  summarise(risk_count = n())

#Count of Disease Traits for Risk Allele = 2
subject_2_trait_count_disease2 <- subject_2_risk_data %>% 
  group_by(DISEASE.TRAIT, have_risk_allele_count) %>% 
  filter(have_risk_allele_count == 2) %>% 
  summarise(risk_count_2 = n()) %>% 
  subset(select = -c(have_risk_allele_count))

#Mutating Risk Counts for Allele = 1 and Allele = 2
subject_2_trait_count_disease3 <- 
  merge(subject_2_trait_count_disease, subject_2_trait_count_disease2,
        by = c("DISEASE.TRAIT", "DISEASE.TRAIT"), all = TRUE)

subject_2_trait_count_disease3[is.na(subject_2_trait_count_disease3)] <- 0 

#Risk Count for Allele 1 equation
subject_2_risk_d_1 <- (subject_2_trait_count_disease3$risk_count - subject_2_trait_count_disease3$risk_count_2)

#Adding risk 1 column in 
subject_2_trait_count_disease4 <- dplyr::select(subject_2_trait_count_disease3,
                                                DISEASE.TRAIT,
                                                risk_count,
                                                risk_count_2) %>% 
  mutate(subject_2_risk_d_1)

#Adding Overall Risk Column and Equation
subject_2_overall_risk_d <- (subject_2_trait_count_disease4$risk_count_2 / ((subject_2_trait_count_disease4$subject_2_risk_d_1 * 2) + subject_2_trait_count_disease4$risk_count_2)) + 
  (subject_2_trait_count_disease4$subject_2_risk_d_1 / ((subject_2_trait_count_disease4$subject_2_risk_d_1 * 2) + subject_2_trait_count_disease4$risk_count_2))

subject_2_overall_risk_count_d <- dplyr::select(subject_2_trait_count_disease4,
                                                DISEASE.TRAIT,
                                                risk_count,
                                                subject_2_risk_d_1,
                                                risk_count_2) %>%
  mutate(subject_2_overall_risk_d*100) 

subject_2_overall_risk_count_d[is.na(subject_2_overall_risk_count_d)] <- 0

#Table of Disease Traits with over 100 risk alleles for MAPPED_TRAIT
subject_2_table1 <- subject_2_overall_risk_count %>% 
  filter(risk_count > 99) %>% 
  arrange(desc(`subject_2_overall_risk * 100`)) %>% 
  subset(select = -c(risk_count, subject_2_risk_1, risk_count_2))

names(subject_2_table1)[1] <- 'Mapped Trait'
names(subject_2_table1)[2] <- 'Overall Risk'

subject_2_table1 <- subject_2_table1 %>% 
  kable() %>% 
  kable_paper("hover", full_width = FALSE)

subject_2_table1

#Table of Disease Traits with over 100 risk alleles for DISEASE.TRAIT 
subject_2_table2 <- subject_2_overall_risk_count_d %>% 
  filter(risk_count > 99) %>% 
  arrange(desc(`subject_2_overall_risk_d * 100`)) %>% 
  subset(select = -c(risk_count, subject_2_risk_d_1, risk_count_2))

names(subject_2_table2)[1] <- 'Disease Trait'
names(subject_2_table2)[2] <- 'Overall Risk'

subject_2_table2 <- subject_2_table2 %>% 
  kable() %>% 
  kable_paper("hover", full_width = FALSE)

subject_2_table2

# Subject 3 -----------------------
#Subject 3 EDA
subject_3_EDA <- subject_3 %>% 
  ggplot(aes(chrom)) + 
  geom_bar() + 
  labs(x = "Chromosome", y = "Number of SNPs") + 
  ggtitle("Subject 3 SNP Distribution")

#Subject 3 Risk Analysis
subject_3_joined <- inner_join(subject_3, updated_gwas_data, 
                               by = c('rsid' = "SNPS"))

subject_3_joined$risk_allele_clean <- str_sub(subject_3_joined$STRONGEST.SNP.RISK.ALLELE, -1)
subject_3_joined$my_allele_1 <- str_sub(subject_3_joined$genotype, 1, 1)
subject_3_joined$my_allele_2 <- str_sub(subject_3_joined$genotype, 2, 2)
subject_3_joined$have_risk_allele_count <- if_else(subject_3_joined$my_allele_1 == 
                                                     subject_3_joined$risk_allele_clean, 1, 0) +
  if_else(subject_3_joined$my_allele_2 == subject_3_joined$risk_allele_clean, 1, 0)

subject_3_risk_data <- dplyr::select(subject_3_joined, 
                                     rsid, 
                                     have_risk_allele_count,
                                     DISEASE.TRAIT, 
                                     risk_allele = risk_allele_clean,
                                     your_genotype = genotype,
                                     RISK.ALLELE.FREQUENCY,
                                     MAPPED_TRAIT,
                                     REPORTED.GENE.S.)
subject_3_risk_data

#Subject 3 Count of Disease Traits
subject_3_trait_count_mapped <- subject_3_risk_data %>% 
  group_by(MAPPED_TRAIT) %>% 
  summarise(risk_count = n())

#Subject 3 Count of Disease Traits for Risk Allele = 2
subject_3_trait_count_mapped2 <- subject_3_risk_data %>% 
  group_by(MAPPED_TRAIT, have_risk_allele_count) %>% 
  filter(have_risk_allele_count == 2) %>% 
  summarise(risk_count_2 = n()) %>% 
  subset(select = -c(have_risk_allele_count))

#Subject 3 Mutating Risk Counts of Allele = 1 and allele = 2
subject_3_trait_count_mapped3 <- 
  merge(subject_3_trait_count_mapped, subject_3_trait_count_mapped2,
        by = c("MAPPED_TRAIT", "MAPPED_TRAIT"), all = TRUE)

subject_3_trait_count_mapped3[is.na(subject_3_trait_count_mapped3)] <- 0

#Risk Count for Allele 1 Equation
subject_3_risk_1 <- (subject_3_trait_count_mapped3$risk_count - subject_3_trait_count_mapped3$risk_count_2)

#Adding risk_1 column in mapped4 table
subject_3_trait_count_mapped4 <- dplyr::select(subject_3_trait_count_mapped3,
                                               MAPPED_TRAIT,
                                               risk_count,
                                               risk_count_2) %>% 
  mutate(subject_3_risk_1)

#Adding Overall Risk Column and Equation
subject_3_overall_risk <- (subject_3_trait_count_mapped4$risk_count_2 / ((subject_3_trait_count_mapped4$subject_3_risk_1 * 2) + subject_3_trait_count_mapped4$risk_count_2)) + 
  (subject_3_trait_count_mapped4$subject_3_risk_1 / ((subject_3_trait_count_mapped4$subject_3_risk_1 * 2) + subject_3_trait_count_mapped4$risk_count_2))

subject_3_overall_risk_count <- dplyr::select(subject_3_trait_count_mapped4,
                                              MAPPED_TRAIT,
                                              risk_count,
                                              risk_count_2,
                                              subject_3_risk_1) %>% 
  mutate(subject_3_overall_risk * 100)

subject_3_overall_risk_count[is.na(subject_3_overall_risk_count)] <- 0

#Get Count of Disease Traits by DISEASE.TRAIT
subject_3_trait_count_disease <- subject_3_risk_data %>% 
  group_by(DISEASE.TRAIT) %>% 
  summarise(risk_count = n())

#Count of Disease Traits for Risk Allele = 2
subject_3_trait_count_disease2 <- subject_3_risk_data %>% 
  group_by(DISEASE.TRAIT, have_risk_allele_count) %>% 
  filter(have_risk_allele_count == 2) %>% 
  summarise(risk_count_2 = n()) %>% 
  subset(select = -c(have_risk_allele_count))

#Mutating Risk Counts for Allele = 1 and Allele = 2
subject_3_trait_count_disease3 <- 
  merge(subject_3_trait_count_disease, subject_3_trait_count_disease2,
        by = c("DISEASE.TRAIT", "DISEASE.TRAIT"), all = TRUE)

subject_3_trait_count_disease3[is.na(subject_3_trait_count_disease3)] <- 0 

#Risk Count for Allele 1 equation
subject_3_risk_d_1 <- (subject_3_trait_count_disease3$risk_count - subject_3_trait_count_disease3$risk_count_2)

#Adding risk 1 column in 
subject_3_trait_count_disease4 <- dplyr::select(subject_3_trait_count_disease3,
                                                DISEASE.TRAIT,
                                                risk_count,
                                                risk_count_2) %>% 
  mutate(subject_3_risk_d_1)

#Adding Overall Risk Column and Equation
subject_3_overall_risk_d <- (subject_3_trait_count_disease4$risk_count_2 / ((subject_3_trait_count_disease4$subject_3_risk_d_1 * 2) + subject_3_trait_count_disease4$risk_count_2)) + 
  (subject_3_trait_count_disease4$subject_3_risk_d_1 / ((subject_3_trait_count_disease4$subject_3_risk_d_1 * 2) + subject_3_trait_count_disease4$risk_count_2))

subject_3_overall_risk_count_d <- dplyr::select(subject_3_trait_count_disease4,
                                                DISEASE.TRAIT,
                                                risk_count,
                                                subject_3_risk_d_1,
                                                risk_count_2) %>%
  mutate(subject_3_overall_risk_d*100) 

subject_3_overall_risk_count_d[is.na(subject_3_overall_risk_count_d)] <- 0

#Table of Disease Traits with over 100 risk alleles for MAPPED_TRAIT
subject_3_table1 <- subject_3_overall_risk_count %>% 
  filter(risk_count > 99) %>% 
  arrange(desc(`subject_3_overall_risk * 100`)) %>% 
  subset(select = -c(risk_count, risk_count_2, subject_3_risk_1))

names(subject_3_table1)[1] <- 'Mapped Trait'
names(subject_3_table1)[2] <- 'Overall Risk'

subject_3_table1 <- subject_3_table1 %>% 
  kable() %>% 
  kable_paper("hover", full_width = FALSE)

subject_3_table1

#Table of Disease Traits with over 100 risk alleles for DISEASE.TRAIT 
subject_3_table2 <- subject_3_overall_risk_count_d %>% 
  filter(risk_count > 99) %>% 
  arrange(desc(`subject_3_overall_risk_d * 100`)) %>% 
  subset(select = -c(risk_count, subject_3_risk_d_1, risk_count_2))

names(subject_3_table2)[1] <- 'Disease Trait'
names(subject_3_table2)[2] <- 'Overall Risk'

subject_3_table2 <- subject_3_table2 %>% 
  kable() %>% 
  kable_paper("hover", full_width = FALSE)

subject_3_table2

# Final Output ----------

subject_1_EDA
subject_1_table1
subject_1_table2

subject_2_EDA
subject_2_table1
subject_2_table2

subject_3_EDA
subject_3_table1
subject_3_table2

