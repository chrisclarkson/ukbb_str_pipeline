library(data.table)

#read and bind GT and clinical data
data=data.frame(fread('data_500k.tsv.gz',header=T),stringsAsFactors=F)
meta_data=data.frame(fread('HTT_project_all_features_meta_file.tsv.gz',sep='\t',header=T),stringsAsFactors=F)
data$VCF=basename(gsub('.vcf','',data$VCF))
data_gt36=data[((data$gene=='HTT') & (data$A2>=36)),]

#assemble ICD10 terms in from cohort with gt36 HTT repeats
icds=gsub(" ",'',gsub("\\]","",gsub("\\[","",gsub("'","",do.call('c',lapply(data_gt_36$p41202_Diagnoses_._main_ICD10,function(x){return(unlist(strsplit(x,',')))}))))))
icds=icds[order(icds[!is.na(icds)])]
counts=table(icds)
icds=names(counts)[counts>=2]
write.table(icds,'icds_gt36.tsv',row.names = F,col.names = F,quote = F)

#get columns that will date the ICD10 terms and prepare file so that ICD10 terms can be dated (using python script)
col_order=colnames(meta_data)[grep(pattern='41262',colnames(meta_data))]
col_order=cbind(col_order,unlist(lapply(col_order,function(x){return(strsplit(x,'._Array_')[[1]][2])})))
col_order=col_order[order(as.numeric(col_order[,2])),1]
out=cbind(meta_data$eid,gsub("'","",meta_data$p41202_Diagnoses_._main_ICD10),apply(meta_data[,col_order],2,as.character.Date))
out[,2]=gsub('\\[','',out[,2])
out[,2]=gsub('\\]','',out[,2])
out[,2]=gsub(' ','',out[,2])

colnames(out)[1:2]=c('eid','codes')
write.table(out,'ukbb_informal_dated_main_icd10_codes.tsv',sep='\t',row.names = F,quote=F)

#then use python parse_icd10_dates.py

