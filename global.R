library(learnr)
library(glue)
library(conflicted)
library(tidyverse)
library(plotly)
library(ggsci)
library(patchwork)
library(mclust)
library(fontawesome)
library(survival)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")

# d <- d0 <- read_rds(path = "./pacSci__TCGA_BRCA_data2.RDS")

d <- d0 <- read_rds(path = "./dataprep/output_data/pacSci__TCGA_BRCA_data.RDS")

# d$pam50_genes <- (genefu::pam50)$centroids %>% rownames() %>% recode(., CDCA1 = "NUF2", KNTC2 = "NDC80")
# saveRDS(object = d,file = "./pacSci__TCGA_BRCA_data2.RDS")


d$gx$Primary_solid_Tumor %>% 
  right_join(d$molecular_subtype %>% dplyr::select(patient, BRCA_Subtype_PAM50), y = .) %>%
  mutate_if(.tbl = .,.predicate = is.numeric,.funs = function(x) log1p(x)) %>%
  mutate_if(.tbl = .,.predicate = is.numeric,.funs = function(x) scale(x,center = T,scale = T)[,1]) -> d$d_pam50


# Clean up: Remove column with too many NAs in subtype
d$molecular_subtype <- d$molecular_subtype %>% 
  mutate_if(.predicate = is.character, .funs = function(x) replace(x,x=="NA",NA)) %>%
  purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=50) 


# Add pre-calculated mclust results for computational efficiency
set.seed(12)

c(g5 = 5, g4 = 4, g3 = 3, g2 = 2) %>% 
  map(.x = ., .f = function(x) Mclust(data = d$gx$Primary_solid_Tumor %>% select(ESR1,PGR) %>% mutate_all(.funs = function(x) log1p(x)),G = x)) -> cls_mods

d$molecular_subtype <- d$gx$Primary_solid_Tumor %>% select(id, patient) %>% 
  cbind(mclust_esr1_pgr_g2 = glue("group {cls_mods$g2$classification}")) %>% 
  cbind(mclust_esr1_pgr_g3 = glue("group {cls_mods$g3$classification}")) %>% 
  cbind(mclust_esr1_pgr_g4 = glue("group {cls_mods$g4$classification}")) %>% 
  cbind(mclust_esr1_pgr_g5 = glue("group {cls_mods$g5$classification}")) %>% 
  left_join(d$molecular_subtype, y = .)


# Generate 11 random numbers
# runif11
rnorm11 <- rnorm(n = 11, mean = 75, sd = 10) %>% round

#' Given positive train and test sample sizes, this function returns valid patient ids
#' @param n_train a positive integer default is the entire dataset
#' @param n_test a positive integer default is 0
#' @param mutually_exclusive T/F determines if train and test should be mutually exclusive, default is T
get_train_test_patient_ids <- function(n_train = nrow(d$d_pam50), n_test = nrow(d$d_pam50) - n_train, mutually_exlusive = T){
  
  small_n_train_threshold <- 100
  N <-  nrow(d$d_pam50)
  
  # In case a vector or an matrix is provided, just get the first element
  n_train <- n_train[1]
  n_test <- n_test[1]
  
  if(!is.numeric(n_train) | !is.numeric(n_test))
    stop("n_train and n_test need to be positive integers")
  
  n_train <- as.integer(n_train)
  n_test <- as.integer(n_test)
  
  
  if(n_train < 0 | n_test < 0)
    stop("n_train and n_test need to be positive integers")
  
  if(mutually_exlusive){
    if((n_train + n_test) > N)
      stop("Invalid sample size")
  }else{
    if(n_train > N | n_test > N)
      stop("Invalid sample size")
  }
  
  
  if(n_train < small_n_train_threshold)
    warning("n_train is very small")
  
  
  patient_ids_train <- sample(x = d$d_pam50$patient, size = n_train, replace = F)
  
  if(n_test == 0)
    return(list(patient_ids_train = patient_ids_train, patient_ids_test = character(0)))
  
  if(mutually_exlusive){
    patient_ids_test <- sample(x = d$d_pam50 %>% filter(!patient %in% patient_ids_train) %>% pull(patient), size = n_test, replace = F)
  }else{
    patient_ids_test <- sample(x = d$d_pam50$patient %>% pull(patient), size = n_test, replace = F)
  }
  
  return(list(patient_ids_train = patient_ids_train, patient_ids_test = patient_ids_test))
  
}


#' Given patient ids, this estimates centroids for the population
#' It takes the d$d_pam50 from the global environment
#' @param patient_ids_train a vector of valid character ids e.g. TCGA-3C-AAAU. By default all ids in the dataset
#' @return trained_centroid a 50 x 5 data.frame 
train_centroids <-function(patient_ids_train = get_train_test_patient_ids()$patient_ids_train){
  
  if(!is.character(patient_ids_train))
    stop("Invalid patient_ids_train They need a be a vector of characters")
  
  training_dat <- d$d_pam50 %>% filter(patient %in% patient_ids_train)
  
  trained_centroid <- training_dat %>%
    group_by(BRCA_Subtype_PAM50) %>% 
    summarise_at(.tbl = ., .vars = vars(d$pam50_genes),.funs = median) %>%
    gather(gene,centroid,-c(BRCA_Subtype_PAM50)) %>% 
    spread(BRCA_Subtype_PAM50,centroid) %>%
    as.data.frame() %>%
    `rownames<-`(.$gene) %>% dplyr::select(-gene)
  
  return(trained_centroid)
}


#' This function returns a model (a function) given ids to train on
#' @param patient_ids_train a vector of valid character ids e.g. TCGA-3C-AAAU. By default all ids in the dataset
#' @return a function that takes a vector of valid charactter ids and predict their subtype
get_model <- function(patient_ids_train = get_train_test_patient_ids()$patient_ids_train){
  trained_centroid <- train_centroids(patient_ids_train = patient_ids_train)
  
  mat_pam50 <- d$d_pam50 %>% as.data.frame() %>% `rownames<-` (.$patient) %>% dplyr::select(rownames(trained_centroid)) %>% as.matrix
  
  
  model <- function(patient_ids_predict){
    cor(t(mat_pam50[patient_ids_predict, , drop = F]), as.matrix(trained_centroid), method = "spearman") %>% 
      as.data.frame()%>% mutate(patient = rownames(.)) %>% as_tibble() %>%
      gather(centroid, spearman, -patient) %>% 
      group_by(patient) %>%
      group_split() %>%
      map(~ .x %>% arrange(desc(spearman)) %>% mutate(predicted_subtype = .$centroid[1])) %>% 
      bind_rows() %>% spread(centroid,spearman) %>% 
      right_join(x = d$molecular_subtype %>% dplyr::select(patient, BRCA_Subtype_PAM50), y = .) %>%
      rename(subtype_pam50_actual = BRCA_Subtype_PAM50, subtype_pam50_predicted = predicted_subtype )
  }
  
  return(model)
  
}

#' This funciton evaluates and returns the quality of calls made
#' @param model_calls a data.frame generated by the trained model
#' @return a list with confusion matrix as well as tallied statistics like ppv and npv
evaluate_calls <- function(model_calls){
  confusion_mat <- model_calls %>% count(subtype_pam50_actual,subtype_pam50_predicted) %>% 
    spread(subtype_pam50_predicted, n) %>%
    mutate_if(.predicate = is.integer, function(x) replace_na(x,0))
  
  
  
  confusion_mat_tally <- confusion_mat %>%
    gather(subtype_pam50_predicted,n,-subtype_pam50_actual) %>% 
    mutate(subtype_pam50_actual = paste0(subtype_pam50_actual,"_actual")) %>%
    mutate(subtype_pam50_predicted = paste0(subtype_pam50_predicted,"_pred")) %>%
    mutate(value = T) %>% spread(subtype_pam50_actual, value, fill = F) %>%
    mutate(value = T) %>% spread(subtype_pam50_predicted, value, fill = F)
  
  
  confusion_mat_stat <- list(Basal = confusion_mat_tally %>% group_by(Basal_actual,Basal_pred) %>% summarise(n = sum(n)) %>% `colnames<-` (c("actual","predicted","n")),
                             Her2 = confusion_mat_tally %>% group_by(Her2_actual,Her2_pred) %>% summarise(n = sum(n)) %>% `colnames<-` (c("actual","predicted","n")),
                             LumA = confusion_mat_tally %>% group_by(LumA_actual,LumA_pred) %>% summarise(n = sum(n)) %>% `colnames<-` (c("actual","predicted","n")),
                             LumB = confusion_mat_tally %>% group_by(LumB_actual,LumB_pred) %>% summarise(n = sum(n)) %>% `colnames<-` (c("actual","predicted","n")),
                             Normal = confusion_mat_tally %>% group_by(Normal_actual,Normal_pred) %>% summarise(n = sum(n)) %>% `colnames<-` (c("actual","predicted","n"))) %>%
    bind_rows(.id = "subtype") %>%
    mutate(call_type = case_when(actual & predicted ~ "TP",
                                 !actual & !predicted ~ "TN",
                                 !actual & predicted ~ "FP",
                                 actual & !predicted ~ "FN")) %>%
    ungroup() %>%
    distinct(subtype,call_type,n) %>%
    spread(call_type,n)
       
  
  confusion_mat_stat %>%  `colnames<-` (tolower(colnames(.))) %>%
    mutate(ppv = tp/(tp+fp), npv = tn/(tn+fn))

  list(confusion_matrix = confusion_mat, confusion_mat_stat = confusion_mat_stat)
}


#' This function simply provides a simple interface for giving a number of grouping g between 2 to 5 
#' and getting the resultant classification plot
#' @param g Number of grouping: either 2, 3, 4 or 5
plot_g_cluster <- function(g){
  
  if(!is.numeric(g))
    stop(glue("The number of groups g is set to {g}. Please provide a number value between 2 and 5"))
  g <- as.integer(g)
  
  if(! g %in% c(2:5))
    stop(glue("The number of groups g is set to {g}. Please provide a number value between 2 and 5"))
  
  clust_col <- glue("mclust_esr1_pgr_g{g}")
  
  dg <- d$gx$Primary_solid_Tumor %>% 
    select(id, ESR1, PGR) %>% 
    left_join(x = ., y = d$molecular_subtype %>% select("id",all_of(clust_col)),by="id") 
  
  dg %>%
    ggplot(aes(x = log1p(ESR1), y = log1p(PGR), colour = mclust_esr1_pgr_g2)) +
    geom_point(aes_string(colour = clust_col)) + 
    ggsci::scale_color_jama() + 
    guides(colour=guide_legend(title="Algorithm classified groups")) + 
    theme(legend.position = "bottom")
}