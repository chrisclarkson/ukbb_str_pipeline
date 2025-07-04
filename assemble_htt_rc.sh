for i in {10..60} 
do
path="/mnt/project/cc_analysis/output_dir/RC/"
file_v0="1-4_df_${i}.tsv"
file_v2="df_${i}.tsv"
full_path_v0="${path}${file_v0}"
full_path_v2="${path}${file_v2}"
Rscript clean_data_count_alleles.R $full_path_v0 cleaned_${file_v0}
Rscript format_RC_data.R $full_path_v2 fmtd_${file_v2}
Rscript format_RC_data_v0.R cleaned_${file_v0} fmtd_${file_v0}
Rscript phasing_args.R fmtd_${file_v2} fmtd_${file_v0} phased_${i}.tsv
done
rm fmt*
rm cleaned*
rm alleles_*
R
library(data.table);
data=do.call('rbind',lapply(list.files(pattern='phased_'),fread,sep='\t',header=T))
write.table(data,'rc_ukbb_htt_phased.tsv',sep='\t',row.names=F,col.names=T,quote=F)