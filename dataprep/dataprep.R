# This script downloads and prepares the TCGA BRCA dataset and athe associated metadata 
library(conflicted)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(genefu)
library(tidyverse)
library(naniar)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
data("pam50")


d <-list()
d$clin <- clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
subtypes <- PanCancerAtlas_subtypes()

#----------------------------------------------------------------------------------------
intermediate_files <- c("./intermediate_files//brca_exp.RDS", "./intermediate_files//dataFilt.RDS", 
                        "./intermediate_files//dataSmNT.RDS", "./intermediate_files//dataSmTP.RDS",
                        "./intermediate_files/dataPrep.RDS","./intermediate_files/dataNorm.RDS", 
                        "./intermediate_files//dataSubt.RDS")


if(any(!file.exists( intermediate_files ))){
  query.exp <- GDCquery(project = "TCGA-BRCA", 
                        legacy = TRUE,
                        data.category = "Gene expression",
                        data.type = "Gene expression quantification",
                        platform = "Illumina HiSeq", 
                        file.type = "results",
                        experimental.strategy = "RNA-Seq",
                        sample.type = c("Primary solid Tumor","Solid Tissue Normal"))
  
  GDCdownload(query.exp)
  brca.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "brcaExp.rda")
  
  
  # get subtype information
  dataSubt <- TCGAquery_subtype(tumor = "BRCA")
  
  # Which samples are Primary Tumor
  dataSmTP <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"TP") 
  # which samples are solid tissue normal
  dataSmNT <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"NT")
  
  
  dataPrep <- TCGAanalyze_Preprocessing(object = brca.exp, cor.cut = 0.6)                      
  
  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                        geneInfo = geneInfo,
                                        method = "gcContent")                
  
  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                    method = "quantile", 
                                    qnt.cut =  0.25) 
  
  
  write_rds(x = brca.exp, path = "./intermediate_files/brca_exp.RDS")
  write_rds(x = dataSubt, path = "./intermediate_files/dataSubt.RDS")
  write_rds(x = dataSmTP, path = "./intermediate_files/dataSmTP.RDS")
  write_rds(x = dataSmNT, path = "./intermediate_files/dataSmNT.RDS")
  write_rds(x = dataPrep, path = "./intermediate_files/dataPrep.RDS")
  write_rds(x = dataNorm, path = "./intermediate_files/dataNorm.RDS")
  write_rds(x = dataFilt, path = "./intermediate_files/dataFilt.RDS")
  
}else{
  dataSmTP <- read_rds(path = "./intermediate_files/dataSmTP.RDS")
  dataSmNT <- read_rds(path ="./intermediate_files/dataSmNT.RDS")
  dataNorm <- read_rds(path = "./intermediate_files/dataNorm.RDS")
  dataSubt <- read_rds(path = "./intermediate_files/dataSubt.RDS")
  d$clin <- d$clin %>% rename(patient = submitter_id)
}


pam50_genes <- recode(rownames(pam50$centroids), CDCA1 = "NUF2", KNTC2 = "NDC80")
single_genes <- c("AR","ESR1","PGR","ERBB2", "MKI67","BRCA1", "BRCA2")
all_genes <- unique(c(pam50_genes, single_genes))

x <- t(dataNorm[all_genes,])


d$gx <- x %>% 
  data.frame(.,stringsAsFactors = F) %>% 
  mutate(id = rownames(.)) %>% 
  mutate(patient = substr(id,start = 1,12)) %>%
  mutate(tumor_or_normal =   case_when(id %in% dataSmTP ~ "Primary solid Tumor",
                                       id %in% dataSmNT ~ "Solid Tissue Normal")) %>%
  dplyr::select(c("id","patient","tumor_or_normal",all_of(all_genes)))


d$gx <- d$gx %>%
  filter(patient %in% intersect(d$gx$patient, d$clin$patient))

d$gx <- d$gx %>% split(x = .,f = .$tumor_or_normal)
names(d$gx) <- gsub(pattern = "\\.","_",x=make.names(names(d$gx)))

d$clin <- d$clin %>% dplyr::select(patient, gender, race, days_to_birth, days_to_diagnosis, age_at_diagnosis, days_to_death, days_to_last_follow_up ,vital_status, tumor_grade,morphology, prior_malignancy) 

d$clin <- d$clin %>%   filter(patient %in% intersect(d$gx$Primary_solid_Tumor$patient, d$clin$patient))  


d$clin <- d$clin %>% filter(patient %in% intersect(dataSubt$patient, d$gx$Primary_solid_Tumor$patient))
d$gx$Solid_Tissue_Normal <- d$gx$Solid_Tissue_Normal %>% filter(patient %in% intersect(dataSubt$patient, d$gx$Primary_solid_Tumor$patient))
d$gx$Primary_solid_Tumor <- d$gx$Primary_solid_Tumor %>% filter(patient %in% intersect(dataSubt$patient, d$gx$Primary_solid_Tumor$patient))
dataSubt <- dataSubt %>% filter(patient %in% intersect(dataSubt$patient, d$gx$Primary_solid_Tumor$patient))
d$molecular_subtype <- dataSubt

dir.exists("./output_data"){
  dir.create("./output_data")
}

write_rds(x = d,path = "./output_data/pacSci__TCGA_BRCA_data.RDS")
