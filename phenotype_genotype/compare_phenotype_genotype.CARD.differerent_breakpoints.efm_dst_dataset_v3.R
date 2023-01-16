
# This R script is used to compare AST phenotypes with WGS-predicted RGI CARD phenotypes per antibiotic
# Here, the latest rgi version (v5.2.2) and CARD database v3.1.4 downloaded on 03/03/2022, AST breakpoints and dataset efm_dst_dataset_v3 are used.

###################################################################################################################################################
###                                                       LOADING DST METADATA                                                                  ###                     
###################################################################################################################################################

phen_metadata_file = "efm_dst_dataset_v3.dst_metadata.different_breakpoints.csv";
phen_table = read.delim(phen_metadata_file, sep = "\t", header = T);
dim(phen_table)
# [1] 4730  225
# [1] 4731  233

###################################################################################################################################################
###                                                    LOADING RGI-CARD PREDICTIONS                                                             ###                     
###################################################################################################################################################

rgicard_phen_file = "efm_dst_dataset_v3.rgi_card_phenotypes.csv";
rgicard_phen_table = read.delim(rgicard_phen_file, sep = "\t", header = T)
dim(rgicard_phen_table)
# [1] 4730   30

rgicard_gen_file = "efm_dst_dataset_v3.rgi_card_genotypes.csv";
rgicard_gen_table = read.delim(rgicard_gen_file, sep = "\t", header = T)
dim(rgicard_gen_table)
# [1] 4730   30

# Making sure samples have the same order in all three tables
identical(as.vector(phen_table$sample_id), as.vector(rgicard_phen_table$sample))
# [1] FALSE
identical(as.vector(phen_table$sample_id), as.vector(rgicard_gen_table$sample))
# [1] FALSE

xxx = match(as.vector(phen_table$sample_id), as.vector(rgicard_phen_table$sample))
length(which(is.na(xxx)))
# [1] 0

rgicard_phen_table = rgicard_phen_table[xxx,]
rgicard_gen_table = rgicard_gen_table[xxx,]
identical(as.vector(phen_table$sample_id), as.vector(rgicard_phen_table$sample))
# [1] TRUE
identical(as.vector(phen_table$sample_id), as.vector(rgicard_gen_table$sample))
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
###                             MAKING ANTIBIOTIC COMPARISONS POSSIBLE BETWEEN METADATA AND CARD-RGI TABLES                                    ###                     
###################################################################################################################################################

# all antibiotics with phenotypic data
all_antibiotics = sort(gsub("_RIS", "", colnames(phen_table)[which(grepl("_RIS",colnames(phen_table))==TRUE)]))
length(all_antibiotics)
# [1] 33
# [1] 35

# Antibiotics with at least 50 susceptible and 50 resistant isolates across studies are presented in this table, antibiotics not meeting these criteria 
discarded_antibiotics = c("streptogramin.A", "streptogramin.B","penicillin","tylosin", "tedizolid", "dalbavancin", "cefotaxime", "fosfomycin", "bacitracin", "levofloxacin", "azithromycin", "lincomycin", "virginiamycin", "chloramphenicol", "sulfamethazine", "salinomycin", "avilamycin","lincosamide","nitrofurantoin","trimethoprim","piperacillin_tazobactam","imipenem")
# NOTE: tigecycline is not discarded in efm_dst_dataset_v3 database
antibiotics_selected = all_antibiotics[which(is.na(match(all_antibiotics, discarded_antibiotics)))]
length(antibiotics_selected)
# [1] 15


#### Encoding for cross-resistance between CARD antibiotic classes and metadata antibiotics
# Antibiotic classes predicted by CARD-RGI
colnames(rgicard_phen_table)
# [1] "sample"                                     "acridine.dye"                               "aminocoumarin.antibiotic"                  
# [4] "aminoglycoside.antibiotic"                  "carbapenem"                                 "cephalosporin"                             
# [7] "cephamycin"                                 "diaminopyrimidine.antibiotic"               "disinfecting.agents.and.intercalating.dyes"
# [10] "elfamycin.antibiotic"                       "fluoroquinolone.antibiotic"                 "glycopeptide.antibiotic"                   
# [13] "glycylcycline"                              "isoniazid"                                  "lincosamide.antibiotic"                    
# [16] "macrolide.antibiotic"                       "monobactam"                                 "nitrofuran.antibiotic"                     
# [19] "nucleoside.antibiotic"                      "oxazolidinone.antibiotic"                   "penam"                                     
# [22] "penem"                                      "peptide.antibiotic"                         "phenicol.antibiotic"                       
# [25] "pleuromutilin.antibiotic"                   "rifamycin.antibiotic"                       "streptogramin.antibiotic"                  
# [28] "sulfonamide.antibiotic"                     "tetracycline.antibiotic"                    "triclosan"    

rgicard_phen_table$ampicillin = rgicard_phen_table$penam
rgicard_phen_table$penicillin = rgicard_phen_table$penam
rgicard_phen_table$vancomycin = rgicard_phen_table$glycopeptide.antibiotic
rgicard_phen_table$teicoplanin = rgicard_phen_table$glycopeptide.antibiotic
# rgicard_phen_table$dalbavancin = NA --> not added, not predicted by CARD
# rgicard_phen_table$daptomycin = NA --> not added, not predicted by CARD
rgicard_phen_table$ciprofloxacin = rgicard_phen_table$fluoroquinolone.antibiotic
rgicard_phen_table$tetracycline = rgicard_phen_table$tetracycline.antibiotic
rgicard_phen_table$doxycycline = rgicard_phen_table$tetracycline.antibiotic
rgicard_phen_table$tigecycline = rgicard_phen_table$glycylcycline
rgicard_phen_table$streptomycin = rgicard_phen_table$aminoglycoside.antibiotic 
rgicard_phen_table$gentamicin = rgicard_phen_table$aminoglycoside.antibiotic
rgicard_phen_table$kanamycin = rgicard_phen_table$aminoglycoside.antibiotic
rgicard_phen_table$erythromycin = rgicard_phen_table$macrolide.antibiotic
rgicard_phen_table$tylosin = rgicard_phen_table$macrolide.antibiotic 
rgicard_phen_table$clindamycin = rgicard_phen_table$lincosamide.antibiotic 
rgicard_phen_table$quinupristin_dalfopristin = rgicard_phen_table$streptogramin.antibiotic 
rgicard_phen_table$linezolid = rgicard_phen_table$oxazolidinone.antibiotic 
rgicard_phen_table$tedizolid = rgicard_phen_table$oxazolidinone.antibiotic


rgicard_gen_table$ampicillin = rgicard_gen_table$penam
rgicard_gen_table$penicillin = rgicard_gen_table$penam
rgicard_gen_table$vancomycin = rgicard_gen_table$glycopeptide.antibiotic
rgicard_gen_table$teicoplanin = rgicard_gen_table$glycopeptide.antibiotic
# rgicard_gen_table$dalbavancin = NA --> not added, not predicted by CARD
# rgicard_gen_table$daptomycin = NA --> not added, not predicted by CARD
rgicard_gen_table$ciprofloxacin = rgicard_gen_table$fluoroquinolone.antibiotic
rgicard_gen_table$tetracycline = rgicard_gen_table$tetracycline.antibiotic
rgicard_gen_table$doxycycline = rgicard_gen_table$tetracycline.antibiotic
rgicard_gen_table$tigecycline = rgicard_gen_table$glycylcycline
rgicard_gen_table$streptomycin = rgicard_gen_table$aminoglycoside.antibiotic 
rgicard_gen_table$gentamicin = rgicard_gen_table$aminoglycoside.antibiotic
rgicard_gen_table$kanamycin = rgicard_gen_table$aminoglycoside.antibiotic
rgicard_gen_table$erythromycin = rgicard_gen_table$macrolide.antibiotic
rgicard_gen_table$tylosin = rgicard_gen_table$macrolide.antibiotic 
rgicard_gen_table$clindamycin = rgicard_gen_table$lincosamide.antibiotic 
rgicard_gen_table$quinupristin_dalfopristin = rgicard_gen_table$streptogramin.antibiotic 
rgicard_gen_table$linezolid = rgicard_gen_table$oxazolidinone.antibiotic 
rgicard_gen_table$tedizolid = rgicard_gen_table$oxazolidinone.antibiotic

xxx = match(antibiotics_selected, colnames(rgicard_phen_table))
antibiotics_selected[which(is.na(xxx))]
# [1] "daptomycin" 
antibiotics_selected = antibiotics_selected[-which(is.na(xxx))]
length(antibiotics_selected)
# [1] 14

###################################################################################################################################################
###                                                     EXTRACTING CONCONDANCE STATISTICS                                                      ###                     
###################################################################################################################################################

summary_table = mat.or.vec(1,13)
colnames(summary_table) = c('antibiotic', 'breakpoint_used','num_R', 'num_I', 'num_S', 'num_NA', 'num_total', 'tp', 'fn', 'tn', 'fp', 'sen', 'spe')
criteria = c("chosen_R", "chosen_NonS", "EUCAST_CB_R", "EUCAST_CB_NonS", "CLSI_R","CLSI_NonS")

for(a in 1:length(antibiotics_selected))
{
  antibiotic = antibiotics_selected[a]
  gen_column = which(colnames(rgicard_phen_table) == antibiotic)
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
    tp = length(which(phen_table[,phen_column] == 'R' & rgicard_phen_table[,gen_column] == 'R'))
    fn = length(which(phen_table[,phen_column] == 'R' & rgicard_phen_table[,gen_column] == 'S'))
    sen = (tp/(tp+fn))*100
    # specificity
    tn = length(which(phen_table[,phen_column] == 'S' & rgicard_phen_table[,gen_column] == 'S'))
    fp = length(which(phen_table[,phen_column] == 'S' & rgicard_phen_table[,gen_column] == 'R'))
    spe = (tn/(tn+fp))*100
    print(paste(antibiotic, ' sen ', sen, ' spe ', spe, sep = ))
    newrow = c(antibiotic, criterium, num_R, num_I, num_S, num_NA, num_total, tp, fn, tn, fp, sen, spe)
    summary_table = rbind(summary_table, newrow)
  }
}
summary_table = summary_table[-1,]
output_file = 'efm_dst_dataset_v3.rgi_card.different_breakpoints.accuracy.csv'
write.table(summary_table, file = output_file, col.names = T, row.names = F, quote = F, sep = '\t')

