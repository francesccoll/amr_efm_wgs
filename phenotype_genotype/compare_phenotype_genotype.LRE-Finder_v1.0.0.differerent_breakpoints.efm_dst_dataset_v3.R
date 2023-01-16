# This R script is used to compare AST linezolid phenotypes with LRE-Finder predictions
# This is done for only for run accessions, not assemblies

###################################################################################################################################################
###                                                       LOADING DST METADATA                                                                  ###                     
###################################################################################################################################################

phen_metadata_file = "efm_dst_dataset_v3.dst_metadata.different_breakpoints.csv";
phen_table = read.delim(phen_metadata_file, sep = "\t", header = T);
dim(phen_table)
# [1] 4730  225
# [1] 4731  233

# Removing assemblies, keeping run_accessions only
tmp = which(grepl("GCA_", phen_table$sample_id))
phen_table = phen_table[-tmp,]
dim(phen_table)
# [1] 4186  225
# [1] 4187  233

###################################################################################################################################################
###                                                         LOADING PREDICTIONS                                                                 ###                     
###################################################################################################################################################

predicted_phen_file = "efm_dst_dataset_v3.run_accessions.lre_finder.txt";
predicted_phen_table = read.delim(predicted_phen_file, sep = "\t", header = T)
dim(predicted_phen_table)
# [1] 4186    3
# [1] 4186    3

# making sample ids order to match between files

xxx = match(phen_table$sample_id, predicted_phen_table$sample_id)
identical(phen_table$sample_id, predicted_phen_table$sample_id[xxx])
# [1] TRUE
predicted_phen_table = predicted_phen_table[xxx,]
identical(phen_table$sample_id, predicted_phen_table$sample_id)
# [1] TRUE

output_file = 'ordered.csv'
write.table(predicted_phen_table, file = output_file, col.names = T, row.names = F, quote = F, sep = '\t')



###################################################################################################################################################
###                                                     EXTRACTING CONCONDANCE STATISTICS                                                      ###                     
###################################################################################################################################################

summary_table = mat.or.vec(1,13)
colnames(summary_table) = c('antibiotic', 'breakpoint_used','num_R', 'num_I', 'num_S', 'num_NA', 'num_total', 'tp', 'fn', 'tn', 'fp', 'sen', 'spe')
criteria = c("chosen_R", "chosen_NonS", "EUCAST_CB_R", "EUCAST_CB_NonS", "CLSI_R","CLSI_NonS")

antibiotics_selected = c("linezolid")

for(a in 1:length(antibiotics_selected))
{
  antibiotic = antibiotics_selected[a]
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
    tp = length(which(phen_table[,phen_column] == 'R' & predicted_phen_table[,"predicted_phenotype"] == 'R'))
    fn = length(which(phen_table[,phen_column] == 'R' & predicted_phen_table[,"predicted_phenotype"] == 'S'))
    sen = (tp/(tp+fn))*100
    # specificity
    tn = length(which(phen_table[,phen_column] == 'S' & predicted_phen_table[,"predicted_phenotype"] == 'S'))
    fp = length(which(phen_table[,phen_column] == 'S' & predicted_phen_table[,"predicted_phenotype"] == 'R'))
    spe = (tn/(tn+fp))*100
    print(paste(antibiotic, ' sen ', sen, ' spe ', spe, sep = ))
    newrow = c(antibiotic, criterium, num_R, num_I, num_S, num_NA, num_total, tp, fn, tn, fp, sen, spe)
    summary_table = rbind(summary_table, newrow)
  }
}
summary_table = summary_table[-1,]
output_file = 'efm_dst_dataset_v3.lre_finder_v1.0.0.different_breakpoints.accuracy.csv'
write.table(summary_table, file = output_file, col.names = T, row.names = F, quote = F, sep = '\t')

