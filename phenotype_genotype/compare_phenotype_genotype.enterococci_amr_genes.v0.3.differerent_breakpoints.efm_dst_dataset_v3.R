
# This R script is used to compare AST phenotypes with WGS-predicted enterococci_amr_genes phenotypes per antibiotic

###################################################################################################################################################
###                                                       LOADING DST METADATA                                                                  ###                     
###################################################################################################################################################

phen_metadata_file = "efm_dst_dataset_v3.dst_metadata.different_breakpoints.csv";
phen_table = read.delim(phen_metadata_file, sep = "\t", header = T);
dim(phen_table)
# [1] 4731  248

###################################################################################################################################################
###                                               DISCARDING SAMPLES BASED ON TYPE OF GENOMIC DATA                                              ###                     
###################################################################################################################################################

# For some antibiotics, assemblies cannot be used (e.g. to predict linezolid and tedizolid resistance)
tmp = which(grepl("GCA_", phen_table$sample_id))
criteria = c("chosen_R", "chosen_NonS", "EUCAST_CB_R", "EUCAST_CB_NonS", "CLSI_R","CLSI_NonS","EUCAST_ECOFF_R")
antibiotics_tmp = c("linezolid", "tedizolid")

for(c in 1:length(criteria))
{
  criterium = criteria[c]
  for(a in 1:length(antibiotics_tmp))
  {
    antibiotic = antibiotics_tmp[a]
    column_name = paste(antibiotic,"_RS_breakpoint_",criterium, sep = "")
    phen_table[tmp,column_name] = NA
  }
}

###################################################################################################################################################
###                                                         LOADING PREDICTIONS                                                                 ###                     
###################################################################################################################################################

predicted_phen_file = "efm_dst_dataset_v3.enterococci_amr_genes.v0.3.ariba.flag_vcf.pheno.csv";
predicted_phen_table = read.delim(predicted_phen_file, sep = "\t", header = T)
dim(predicted_phen_table)
# [1] 4730   29

predicted_geno_file = "efm_dst_dataset_v3.enterococci_amr_genes.v0.3.ariba.flag_vcf.geno.csv";
predicted_geno_table = read.delim(predicted_geno_file, sep = "\t", header = T)
dim(predicted_geno_table)
# [1] 4730   29


# Making sure samples have the same order in both tables
identical(as.vector(phen_table$sample_id), as.vector(predicted_phen_table$sample))
identical(as.vector(phen_table$sample_id), as.vector(predicted_geno_table$sample))
# [1] FALSE

xxx = match(as.vector(phen_table$sample_id), as.vector(predicted_phen_table$sample))
length(which(is.na(xxx)))
# [1] 0

predicted_phen_table = predicted_phen_table[xxx,]
predicted_geno_table = predicted_geno_table[xxx,]
identical(as.vector(phen_table$sample_id), as.vector(predicted_phen_table$sample))
# [1] TRUE
identical(as.vector(phen_table$sample_id), as.vector(predicted_geno_table$sample))
# [1] TRUE


###################################################################################################################################################
###                                               KEEPING COMMON ANTIBIOTIC COLUMNS BETWEEN TABLES                                              ###                     
###################################################################################################################################################

# selected_pAST_antibiotics contains the list of antibiotics for which there is enough R and S isolates, see efm_dst_dataset_v2.dst_metadata.summary_table.xlsx
# or antibiotic in both phenotypic metdata table and look-up table
selected_pAST_antibiotics = c("ampicillin", "penicillin", "fosfomycin", "vancomycin", "teicoplanin", "dalbavancin", "daptomycin", "ciprofloxacin", "levofloxacin","tetracycline", "doxycycline", 
                              "tigecycline", "streptomycin", "gentamicin", "kanamycin", "azithromycin", "erythromycin", "tylosin", "clindamycin", "streptogramin.A", "streptogramin.B", "quinupristin_dalfopristin", "virginiamycin",
                              "chloramphenicol", "linezolid", "tedizolid")

# excluded antibiotics: nitrofurantoin, cefotaxime, bacitracin, lincosamide, lincomycin, trimethoprim, sulfamethazine, salinomycin, avilamycin
# see enterococci_amr_genes.xlsx, tab 'not included', for reason of exclusion

# matching relevant antibiotics between tables
# efm_dst_dataset_v2.dst_metadata.csv   enterococci_amr_genes
#  ampicillin --> ampicillin
#  penicillin --> ampicillin
#  fosfomycin --> fosfomycin
#  vancomycin --> vancomycin
#  teicoplanin --> teicoplanin
#  dalbavancin --> dalbavancin
#  daptomycin --> daptomycin
#  ciprofloxacin --> fluoroquinolones  --> new column to add
#  levofloxacin --> fluoroquinolones  --> new column to add
#  tetracycline --> tetracycline
#  doxycycline --> tetracycline --> new column to add
#  tigecycline --> tigecycline
#  streptomycin --> streptomycin
#  gentamicin --> gentamicin
#  kanamycin --> kanamycin
#  azithromycin --> azithromycin
#  erythromycin --> erythromycin
#  tylosin --> tylosin
#  clindamycin --> clindamycin
#  quinupristin_dalfopristin --> quinupristin.dalfopristin --> new column to add
#  virginiamycin --> quinupristin_dalfopristin
#  chloramphenicol --> chloramphenicol
#  linezolid --> linezolid
#  tedizolid --> linezolid
#  nitrofurantoin --> not included

predicted_phen_table$penicillin = predicted_phen_table$ampicillin
predicted_phen_table$ciprofloxacin = predicted_phen_table$fluoroquinolones
predicted_phen_table$levofloxacin = predicted_phen_table$fluoroquinolones
predicted_phen_table$doxycycline = predicted_phen_table$tetracycline
predicted_phen_table$quinupristin_dalfopristin = predicted_phen_table$quinupristin.dalfopristin
predicted_phen_table$virginiamycin = predicted_phen_table$quinupristin_dalfopristin
predicted_phen_table$tedizolid = predicted_phen_table$linezolid


# Assigning quinupristin-dalfopristin resistance based on resistance to both streptogramin.A and streptogramin.B
res_pred1 = which(predicted_phen_table$quinupristin_dalfopristin == "resistance")
res_pred2 = which(predicted_phen_table$streptogramin.A == "resistance" & predicted_phen_table$streptogramin.B == "resistance")
quinupristin_dalfopristin = as.vector(predicted_phen_table$quinupristin_dalfopristin)
quinupristin_dalfopristin[c(res_pred1, res_pred2)] = "resistance"
predicted_phen_table$quinupristin_dalfopristin = quinupristin_dalfopristin
phen_table$streptogramin.A = NA
phen_table$streptogramin.B = NA
quinupristin_dalfopristin = paste(predicted_geno_table$streptogramin.A, predicted_geno_table$streptogramin.B, sep = ";")
predicted_geno_table$quinupristin_dalfopristin = quinupristin_dalfopristin


# ### matching sure ALL selected antibiotics are found in tables
# xxx = match(selected_pAST_antibiotics, colnames(predicted_phen_table))
# selected_pAST_antibiotics[which(is.na(xxx))]
# # character(0)
# # NOTE: newly assigned RS labels from EUCAST breakpoints are used (see script re-assign_SR_labels_from_breakpoints.R and file phenotype_AST_table.docx)
# colnames(phen_table) = gsub("_RS_eucast$", "", colnames(phen_table))
# xxx = match(selected_pAST_antibiotics, colnames(phen_table))
# selected_pAST_antibiotics[which(is.na(xxx))]
# # [1] "fosfomycin"      "levofloxacin"    "azithromycin"    "virginiamycin"   "chloramphenicol"
# # character(0)
# 
# ### keeping common antibiotic columns between the two tables
# colnames(predicted_phen_table)[1] = "sample_id"
# common_columns = match(colnames(phen_table), colnames(predicted_phen_table))
# phen_table = phen_table[,which(!is.na(common_columns))]
# predicted_phen_table = predicted_phen_table[,common_columns[which(!is.na(common_columns))]]
# identical(colnames(phen_table), colnames(predicted_phen_table))
# # [1] TRUE

dim(phen_table)
# [1] 4731  233
# [1] 4731  241
# [1] 4731  257
length(selected_pAST_antibiotics)
# [1] 26

###################################################################################################################################################
###                                                     EXTRACTING CONCONDANCE STATISTICS                                                      ###                     
###################################################################################################################################################

# Antibiotics with at least 50 susceptible and 50 resistant isolates across studies are presented in this table, antibiotics not meeting these criteria 
discarded_antibiotics = c("streptogramin.A", "streptogramin.B","penicillin","tylosin", "tedizolid", "dalbavancin", "cefotaxime", "fosfomycin", "bacitracin", "levofloxacin", "azithromycin", "lincomycin", "virginiamycin", "chloramphenicol", "sulfamethazine", "salinomycin", "avilamycin","lincosamide","nitrofurantoin","trimethoprim")
# NOTE: tigecycline is not discarded in efm_dst_dataset_v3 database
antibiotics_selected = selected_pAST_antibiotics[which(is.na(match(selected_pAST_antibiotics, discarded_antibiotics)))]
length(antibiotics_selected)
# [1] 15

summary_table = mat.or.vec(1,13)
colnames(summary_table) = c('antibiotic', 'breakpoint_used','num_R', 'num_I', 'num_S', 'num_NA', 'num_total', 'tp', 'fn', 'tn', 'fp', 'sen', 'spe')
criteria = c("chosen_R", "chosen_NonS", "EUCAST_CB_R", "EUCAST_CB_NonS", "CLSI_R","CLSI_NonS","EUCAST_ECOFF_R")

for(a in 1:length(antibiotics_selected))
{
  antibiotic = antibiotics_selected[a]
  gen_column = which(colnames(predicted_phen_table) == antibiotic)
  for(c in 1:length(criteria))
  {
    criterium = criteria[c]
    phen_column = paste(antibiotic,"_RS_breakpoint_",criterium, sep = "")
    num_R = length(which(phen_table[,phen_column] == 'R'))
    num_I = length(which(phen_table[,phen_column] == 'I'))
    num_S = length(which(phen_table[,phen_column] == 'S'))
    num_NA = length(which(is.na(phen_table[,phen_column])))
    num_total = num_R + num_I + num_S + num_NA
    to_print = c('antibiotic ', antibiotic, ' num_R ', num_R, ' num_I ', num_I, ' num_S ', num_S, ' num_NA ', num_NA, ' num_total ', num_total)
    print(paste(to_print, sep = ""))
    # sensitivity - ignoring I category
    tp = length(which(phen_table[,phen_column] == 'R' & predicted_phen_table[,gen_column] == 'resistance'))
    fn = length(which(phen_table[,phen_column] == 'R' & predicted_phen_table[,gen_column] == 'susceptible'))
    sen = (tp/(tp+fn))*100
    # specificity
    tn = length(which(phen_table[,phen_column] == 'S' & predicted_phen_table[,gen_column] == 'susceptible'))
    fp = length(which(phen_table[,phen_column] == 'S' & predicted_phen_table[,gen_column] == 'resistance'))
    spe = (tn/(tn+fp))*100
    print(paste(antibiotic, ' sen ', sen, ' spe ', spe, sep = ))
    newrow = c(antibiotic, criterium, num_R, num_I, num_S, num_NA, num_total, tp, fn, tn, fp, sen, spe)
    summary_table = rbind(summary_table, newrow)
  }
}

output_file = 'efm_dst_dataset_v3.enterococci_amr_genes.v0.3.ariba.flag_vcf.different_breakpoints.accuracy.csv'
write.table(summary_table, file = output_file, col.names = T, row.names = F, quote = F, sep = '\t')
