
# This R script is used to compare AST phenotypes with WGS-predicted AMRFinder phenotypes per antibiotic

###################################################################################################################################################
###                                                       LOADING DST METADATA                                                                  ###                     
###################################################################################################################################################

phen_metadata_file = "efm_dst_dataset_v3.dst_metadata.different_breakpoints.csv";
phen_table = read.delim(phen_metadata_file, sep = "\t", header = T);
dim(phen_table)
# [1] 4730  225



###################################################################################################################################################
###                                                   LOADING AMRFINDER PREDICTIONS                                                             ###                     
###################################################################################################################################################

amrfinder_phen_file = "efm_dst_dataset_v3.amrfinder_phenotypes.efm.csv";
amrfinder_phen_table = read.delim(amrfinder_phen_file, sep = "\t", header = T)
dim(amrfinder_phen_table)
# [1] 4730   35

amrfinder_gen_file = "efm_dst_dataset_v3.amrfinder_genotypes.efm.csv";
amrfinder_gen_table = read.delim(amrfinder_gen_file, sep = "\t", header = T)
dim(amrfinder_gen_table)
# [1] 4730   35

# Making sure samples have the same order in both tables
identical(as.vector(phen_table$sample_id), as.vector(amrfinder_phen_table$sample))
# [1] FALSE
identical(as.vector(phen_table$sample_id), as.vector(amrfinder_gen_table$sample))
# [1] FALSE

xxx = match(as.vector(phen_table$sample_id), as.vector(amrfinder_phen_table$sample))
length(which(is.na(xxx)))
# [1] 0

amrfinder_phen_table = amrfinder_phen_table[xxx,]
amrfinder_gen_table = amrfinder_gen_table[xxx,]
identical(as.vector(phen_table$sample_id), as.vector(amrfinder_phen_table$sample))
# [1] TRUE
identical(as.vector(phen_table$sample_id), as.vector(amrfinder_gen_table$sample))
# [1] TRUE

###################################################################################################################################################
###                             MAKING ANTIBIOTIC COMPARISONS POSSIBLE BETWEEN METADATA AND AMRFINDER TABLES                                    ###                     
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


### Matching AMRFinder and phenotype antibiotic names, and encoding for cross-resistance

colnames(amrfinder_phen_table)
# [1] "sample"              "AMIKACIN"            "AMINOGLYCOSIDE"      "APRAMYCIN"           "AVILAMYCIN"          "BACITRACIN"          "BETA.LACTAM"         "CARBAPENEM"          "CEPHALOSPORIN"      
# [10] "CHLORAMPHENICOL"     "DAPTOMYCIN"          "ERYTHROMYCIN"        "FLORFENICOL"         "FLUOROQUINOLONE"     "FOSFOMYCIN"          "GENTAMICIN"          "KANAMYCIN"           "LINCOSAMIDE"        
# [19] "LINEZOLID"           "MACROLIDE"           "OXAZOLIDINONE"       "PLEUROMUTILIN"       "QUATERNARY.AMMONIUM" "QUINOLONE"           "RIFAMYCIN"           "SPECTINOMYCIN"       "STREPTOGRAMIN"      
# [28] "STREPTOMYCIN"        "STREPTOTHRICIN"      "SULFONAMIDE"         "TETRACYCLINE"        "TIGECYCLINE"         "TOBRAMYCIN"          "TRIMETHOPRIM"        "VANCOMYCIN"   

selected_pAST_antibiotics
# [1] "ampicillin"                "ciprofloxacin"             "clindamycin"               "daptomycin"                "doxycycline"               "erythromycin"              "gentamicin"               
# [8] "kanamycin"                 "linezolid"                 "quinupristin_dalfopristin" "streptomycin"              "teicoplanin"               "tetracycline"              "tigecycline"              
# [15] "vancomycin"         

# Assuming cross-resistance between AMRFinder BETA.LACTAM and ampicillin
amrfinder_phen_table$ampicillin = amrfinder_phen_table$BETA.LACTAM
# Assuming cross-resistance between AMRFinder QUINOLONE and ciprofloxacin
# NOTE: FLUOROQUINOLONE predictions were all S/WT
amrfinder_phen_table$ciprofloxacin = amrfinder_phen_table$QUINOLONE
# Assuming cross-resistance between AMRFinder LINCOSAMIDE and clindamycin
amrfinder_phen_table$clindamycin = amrfinder_phen_table$LINCOSAMIDE
amrfinder_phen_table$daptomycin = amrfinder_phen_table$DAPTOMYCIN
# Assuming cross-resistance between doxycycline and TETRACYCLINE
amrfinder_phen_table$doxycycline = amrfinder_phen_table$TETRACYCLINE
amrfinder_phen_table$tetracycline = amrfinder_phen_table$TETRACYCLINE
# Assuming cross-resistance between AMRFinder MACROLIDE and erythromycin
# NOTE: ERYTHROMYCIN predictions were all S/WT
amrfinder_phen_table$erythromycin = amrfinder_phen_table$MACROLIDE
amrfinder_phen_table$gentamicin = amrfinder_phen_table$GENTAMICIN
amrfinder_phen_table$kanamycin = amrfinder_phen_table$KANAMYCIN
amrfinder_phen_table$linezolid = amrfinder_phen_table$LINEZOLID
amrfinder_phen_table$streptomycin = amrfinder_phen_table$STREPTOMYCIN
# Assuming cross-resistance between teicoplanin and VANCOMYCIN
amrfinder_phen_table$teicoplanin = amrfinder_phen_table$VANCOMYCIN
amrfinder_phen_table$tigecycline = amrfinder_phen_table$TIGECYCLINE
amrfinder_phen_table$vancomycin = amrfinder_phen_table$VANCOMYCIN
# Assuming cross-resistance between quinupristin_dalfopristin and STREPTOGRAMIN
amrfinder_phen_table$quinupristin_dalfopristin = amrfinder_phen_table$STREPTOGRAMIN

dim(amrfinder_phen_table)
# [1] 4730   50
length(selected_pAST_antibiotics)
# [1] 15


###################################################################################################################################################
###                                                     EXTRACTING CONCONDANCE STATISTICS                                                      ###                     
###################################################################################################################################################

length(antibiotics_selected)
# [1] 15

summary_table = mat.or.vec(1,13)
colnames(summary_table) = c('antibiotic', 'breakpoint_used','num_R', 'num_I', 'num_S', 'num_NA', 'num_total', 'tp', 'fn', 'tn', 'fp', 'sen', 'spe')
criteria = c("chosen_R", "chosen_NonS", "EUCAST_CB_R", "EUCAST_CB_NonS", "CLSI_R","CLSI_NonS")

for(a in 1:length(antibiotics_selected))
{
  antibiotic = antibiotics_selected[a]
  gen_column = which(colnames(amrfinder_phen_table) == antibiotic)
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
    tp = length(which(phen_table[,phen_column] == 'R' & amrfinder_phen_table[,gen_column] == 'R'))
    fn = length(which(phen_table[,phen_column] == 'R' & amrfinder_phen_table[,gen_column] == 'S'))
    sen = (tp/(tp+fn))*100
    # specificity
    tn = length(which(phen_table[,phen_column] == 'S' & amrfinder_phen_table[,gen_column] == 'S'))
    fp = length(which(phen_table[,phen_column] == 'S' & amrfinder_phen_table[,gen_column] == 'R'))
    spe = (tn/(tn+fp))*100
    print(paste(antibiotic, ' sen ', sen, ' spe ', spe, sep = ))
    newrow = c(antibiotic, criterium, num_R, num_I, num_S, num_NA, num_total, tp, fn, tn, fp, sen, spe)
    summary_table = rbind(summary_table, newrow)
  }
}
summary_table = summary_table[-1,]
output_file = 'efm_dst_dataset_v3.amrfinder_phenotypes.efm.different_breakpoints.accuracy.csv'
write.table(summary_table, file = output_file, col.names = T, row.names = F, quote = F, sep = '\t')
