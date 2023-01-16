
# This R script is used to compare AST phenotypes with WGS-predicted ResFinder phenotypes per antibiotic
# Here, the latest ResFinder version (4.1.10) and ResFinder database were downloaded on 03/03/2022, AST breakpoints and dataset efm_dst_dataset_v3 are used.

###################################################################################################################################################
###                                                       LOADING DST METADATA                                                                  ###                     
###################################################################################################################################################

phen_metadata_file = "efm_dst_dataset_v3.dst_metadata.different_breakpoints.csv";
phen_table = read.delim(phen_metadata_file, sep = "\t", header = T);
dim(phen_table)
# [1] 4730  225
# [1] 4731  233


###################################################################################################################################################
###                                                   LOADING RESFINDER PREDICTIONS                                                             ###                     
###################################################################################################################################################

resfinder_phen_file = "efm_dst_dataset_v3.resfinder_phenotypes.csv";
resfinder_phen_table = read.delim(resfinder_phen_file, sep = "\t", header = T)
dim(resfinder_phen_table)
# [1] 4730   85

resfinder_efm_phen_file = "efm_dst_dataset_v3.resfinder_phenotypes.efm.csv";
resfinder_efm_phen_table = read.delim(resfinder_efm_phen_file, sep = "\t", header = T)
dim(resfinder_efm_phen_table)
# [1] 4730   13

resfinder_gen_file = "efm_dst_dataset_v3.resfinder_genotypes.csv";
resfinder_gen_table = read.delim(resfinder_gen_file, sep = "\t", header = T)
dim(resfinder_gen_table)
# [1] 4730   85

resfinder_efm_gen_file = "efm_dst_dataset_v3.resfinder_genotypes.efm.csv";
resfinder_efm_gen_table = read.delim(resfinder_efm_gen_file, sep = "\t", header = T)
dim(resfinder_efm_gen_table)
# [1] 4730   13

# Making sure samples have the same order in both tables
identical(as.vector(phen_table$sample_id), as.vector(resfinder_phen_table$sample))
identical(as.vector(phen_table$sample_id), as.vector(resfinder_efm_phen_table$sample))
# [1] FALSE

xxx = match(as.vector(phen_table$sample_id), as.vector(resfinder_phen_table$sample))
length(which(is.na(xxx)))
# [1] 0

resfinder_phen_table = resfinder_phen_table[xxx,]
resfinder_efm_phen_table = resfinder_efm_phen_table[xxx,]
resfinder_gen_table = resfinder_gen_table[xxx,]
resfinder_efm_gen_table = resfinder_efm_gen_table[xxx,]


identical(as.vector(phen_table$sample_id), as.vector(resfinder_phen_table$sample))
# [1] TRUE
identical(as.vector(phen_table$sample_id), as.vector(resfinder_efm_phen_table$sample))
# [1] TRUE
identical(as.vector(phen_table$sample_id), as.vector(resfinder_gen_table$sample))
# [1] TRUE
identical(as.vector(phen_table$sample_id), as.vector(resfinder_efm_gen_table$sample))
# [1] TRUE

###################################################################################################################################################
###                                               DISCARDING SAMPLES BASED ON TYPE OF GENOMIC DATA                                              ###                     
###################################################################################################################################################

# For some antibiotics, assemblies cannot be used (e.g. to predict linezolid and tedizolid resistance)
tmp = which(grepl("GCA_", phen_table$sample_id))
criteria = c("chosen_R", "chosen_NonS", "EUCAST_CB_R", "EUCAST_CB_NonS", "CLSI_R","CLSI_NonS")
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
###                             MAKING ANTIBIOTIC COMPARISONS POSSIBLE BETWEEN METADATA AND RESFINDER TABLES                                    ###                     
###################################################################################################################################################

# all antibiotics with phenotypic data
all_antibiotics = sort(gsub("_RIS", "", colnames(phen_table)[which(grepl("_RIS",colnames(phen_table))==TRUE)]))
length(all_antibiotics)
# [1] 33

# Antibiotics with at least 50 susceptible and 50 resistant isolates across studies are presented in this table, antibiotics not meeting these criteria 
discarded_antibiotics = c("streptogramin.A", "streptogramin.B","penicillin","tylosin", "tedizolid", "dalbavancin", "cefotaxime", "fosfomycin", "bacitracin", "levofloxacin", "azithromycin", "lincomycin", "virginiamycin", "chloramphenicol", "sulfamethazine", "salinomycin", "avilamycin","lincosamide","nitrofurantoin","trimethoprim","piperacillin_tazobactam","imipenem")
# NOTE: tigecycline is not discarded in efm_dst_dataset_v3 database
antibiotics_selected = all_antibiotics[which(is.na(match(all_antibiotics, discarded_antibiotics)))]
length(antibiotics_selected)
# [1] 15
selected_pAST_antibiotics = antibiotics_selected


# tedizolid/linezolid cross-resistance assumed
resfinder_phen_table$tedizolid = resfinder_phen_table$linezolid
resfinder_phen_table$quinupristin_dalfopristin = resfinder_phen_table$quinupristin.dalfopristin
xxx = match(selected_pAST_antibiotics, colnames(resfinder_phen_table))
selected_pAST_antibiotics[which(is.na(xxx))]
# [1] "daptomycin"
# NOTE: only daptomycin is not predicted by ResFinder
colnames(resfinder_phen_table)[1] = "sample_id"
resfinder_phen_table = resfinder_phen_table[,c(1,xxx[which(!is.na(xxx))])]
dim(resfinder_phen_table)
# [1] 4730   15

# same operations done for resfinder_efm_phen_table
resfinder_efm_phen_table$tedizolid = resfinder_efm_phen_table$linezolid
resfinder_efm_phen_table$quinupristin_dalfopristin = resfinder_efm_phen_table$quinupristin.dalfopristin
xxx = match(selected_pAST_antibiotics, colnames(resfinder_efm_phen_table))
selected_pAST_antibiotics[which(is.na(xxx))]
# [1] "clindamycin"  "daptomycin"   "doxycycline"  "kanamycin"    "streptomycin"
# NOTE: antibiotic above not predicted by ResFinder when choosing Efm as organism
colnames(resfinder_efm_phen_table)[1] = "sample_id"
resfinder_efm_phen_table = resfinder_efm_phen_table[,c(1,xxx[which(!is.na(xxx))])]
dim(resfinder_efm_phen_table)
# [1] 4397   11

# List of antibiotics matching both selected antibiotics from
resfinder_antibiotics_matched = colnames(resfinder_phen_table)[-1]
resfinder_efm_antibiotics_matched = colnames(resfinder_efm_phen_table)[-1]

# Making sure columns and row order (samples) match between metadata and ResFinder tables
# col = match(colnames(phen_table),colnames(resfinder_phen_table))
# phen_table_matched = phen_table[,which(!is.na(col))]
# resfinder_phen_table = resfinder_phen_table[,col[which(!is.na(col))]]
# identical(colnames(phen_table_matched),colnames(resfinder_phen_table))
# # [1] TRUE
# 
# # same operations done for resfinder_efm_phen_table
# col = match(colnames(phen_table),colnames(resfinder_efm_phen_table))
# phen_table_matched_efm = phen_table[,which(!is.na(col))]
# resfinder_efm_phen_table = resfinder_efm_phen_table[,col[which(!is.na(col))]]
# identical(colnames(phen_table_matched_efm),colnames(resfinder_efm_phen_table))
# # [1] TRUE


###################################################################################################################################################
###                                                     EXTRACTING CONCONDANCE STATISTICS                                                      ###                     
###################################################################################################################################################

summary_table = mat.or.vec(1,13)
colnames(summary_table) = c('antibiotic', 'breakpoint_used','num_R', 'num_I', 'num_S', 'num_NA', 'num_total', 'tp', 'fn', 'tn', 'fp', 'sen', 'spe')
criteria = c("chosen_R", "chosen_NonS", "EUCAST_CB_R", "EUCAST_CB_NonS", "CLSI_R","CLSI_NonS")

for(a in 1:length(resfinder_antibiotics_matched))
{
  antibiotic = resfinder_antibiotics_matched[a]
  gen_column = which(colnames(resfinder_phen_table) == antibiotic)
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
    tp = length(which(phen_table[,phen_column] == 'R' & resfinder_phen_table[,gen_column] == 'Resistant'))
    fn = length(which(phen_table[,phen_column] == 'R' & resfinder_phen_table[,gen_column] == 'No resistance'))
    sen = (tp/(tp+fn))*100
    # specificity
    tn = length(which(phen_table[,phen_column] == 'S' & resfinder_phen_table[,gen_column] == 'No resistance'))
    fp = length(which(phen_table[,phen_column] == 'S' & resfinder_phen_table[,gen_column] == 'Resistant'))
    spe = (tn/(tn+fp))*100
    print(paste(antibiotic, ' sen ', sen, ' spe ', spe, sep = ))
    newrow = c(antibiotic, criterium, num_R, num_I, num_S, num_NA, num_total, tp, fn, tn, fp, sen, spe)
    summary_table = rbind(summary_table, newrow)
  }
}
summary_table = summary_table[-1,]
output_file = 'efm_dst_dataset_v3.resfinder_v4.0.different_breakpoints.accuracy.csv'
write.table(summary_table, file = output_file, col.names = T, row.names = F, quote = F, sep = '\t')



summary_table = mat.or.vec(1,13)
colnames(summary_table) = c('antibiotic', 'breakpoint_used','num_R', 'num_I', 'num_S', 'num_NA', 'num_total', 'tp', 'fn', 'tn', 'fp', 'sen', 'spe')
criteria = c("chosen_R", "chosen_NonS", "EUCAST_CB_R", "EUCAST_CB_NonS", "CLSI_R","CLSI_NonS")

for(a in 1:length(resfinder_efm_antibiotics_matched))
{
  antibiotic = resfinder_efm_antibiotics_matched[a]
  gen_column = which(colnames(resfinder_efm_phen_table) == antibiotic)
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
    # sensitivity
    tp = length(which(phen_table[,phen_column] == 'R' & resfinder_efm_phen_table[,gen_column] == 'Resistant'))
    fn = length(which(phen_table[,phen_column] == 'R' & resfinder_efm_phen_table[,gen_column] == 'No resistance'))
    sen = (tp/(tp+fn))*100
    # specificity
    tn = length(which(phen_table[,phen_column] == 'S' & resfinder_efm_phen_table[,gen_column] == 'No resistance'))
    fp = length(which(phen_table[,phen_column] == 'S' & resfinder_efm_phen_table[,gen_column] == 'Resistant'))
    spe = (tn/(tn+fp))*100
    print(paste(antibiotic, ' sen ', sen, ' spe ', spe, sep = ))
    newrow = c(antibiotic, criterium, num_R, num_I, num_S, num_NA, num_total, tp, fn, tn, fp, sen, spe)
    summary_table = rbind(summary_table, newrow)
  }
}

summary_table = summary_table[-1,]
output_file = 'efm_dst_dataset_v3.resfinder_v4.0_efm.different_breakpoints.accuracy.csv'
write.table(summary_table, file = output_file, col.names = T, row.names = F, quote = F, sep = '\t')

