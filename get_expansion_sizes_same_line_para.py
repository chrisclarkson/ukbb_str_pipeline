# open bash notebook: pip install PyVCF3; python get_expansion_sizes_same_line_para.py /mnt/project/cc_analysis/output_dir/
import vcf
import glob
import sys
import multiprocessing
try:
    directory=sys.argv[1]
except:
    directory='recalled_eh_bams'


def to_file(v,selected,a1,a2,gt,ci1,ci2,rt1,rt2,spanning_no1,spanning_no2,flanking_no1,flanking_no2,irr_no1,irr_no2,name):
    if selected=='na':
        gene='gene'
    else:
        gene=selected.INFO['VARID']
    with(open('expansion_sizes_omit_NAs_'+name+'_same_line.txt','a')) as out:
        out.write(v+'\t'+gene+'\t'+str(a1)+'\t'+str(a2)+'\t'+gt+'\t'+ci1+'\t'+ci2+'\t'+rt1+'\t'+rt2+'\t'+spanning_no1+'\t'+spanning_no2+'\t'+flanking_no1+'\t'+flanking_no2+'\t'+irr_no1+'\t'+irr_no2+'\n')

to_file('VCF','na','A1','A2','GT','CI1','CI2','RT1','RT2','spanning_no1','spanning_no2','flanking_no1','flanking_no2','irr_no1','irr_no2','header')
def parse(files,name):
    for v in files:
        try:
            vcf_reader = vcf.Reader(open(v, 'r'))
            for i in vcf_reader:
                selected,selected_info=i,i.INFO
                if selected.FILTER!=['LowDepth']:
                    hash_fields=dict(i.INFO)
                    hash_fields.update(dict(zip(i.samples[0].data._fields,i.samples[0].data)))
                    gt=hash_fields['GT']
                    gts=gt.split('/')
                    alleles=hash_fields['REPCN'].split('/')
                    if len(alleles)>1:
                        a1,a2=alleles[0],alleles[1]
                    else:
                        a1,a2=alleles[0],'NA'
                    cis=hash_fields['REPCI'].split('/')
                    if len(cis)>1:
                        ci1,ci2=cis[0],cis[1]
                    else:
                        ci1,ci2=cis[0],'NA'
                    rts=hash_fields['SO'].split('/')
                    if len(rts)>1:
                        rt1,rt2=rts[0],rts[1]
                    else:
                        rt1,rt2=rts[0],'NA'
                    spanning_reads=hash_fields['ADSP'].split('/')
                    if len(spanning_reads)>1:
                        spanning_no1,spanning_no2=spanning_reads[0],spanning_reads[1]
                    else:
                        spanning_no1,spanning_no2=spanning_reads[0],'NA'
                    flanking_reads=hash_fields['ADFL'].split('/')
                    if len(flanking_reads)>1:
                        flanking_no1,flanking_no2=flanking_reads[0],flanking_reads[1]
                    else:
                        flanking_no1,flanking_no2=flanking_reads[0],'NA'
                    irrs=hash_fields['ADIR'].split('/')
                    if len(irrs)>1:
                        irr_no1,irr_no2=irrs[0],irrs[1]
                    else:
                        irr_no1,irr_no2=irrs[0],'NA'
                    to_file(v,selected,a1,a2,gt,ci1,ci2,rt1,rt2,spanning_no1,spanning_no2,flanking_no1,flanking_no2,irr_no1,irr_no2,str(name))
        except:
            print('Something wrong with {0}'.format(v))

# for name in range(10,48):
#     files=glob.glob('{0}/*{1}.vcf'.format(directory,str(name)))
#     print('files loaded')
#     parse(files,str(name))



import concurrent.futures

def worker(name, directory):
    print('{0}{1}/*.vcf'.format(directory, str(name)))
    files = glob.glob('{0}{1}/*.vcf'.format(directory, str(name)))
    print(files[0:4])
    print('files loaded for {0}'.format(name))
    parse(files, str(name))
    print('finished')

def parallel_process(directory):
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        futures = [executor.submit(worker, name, directory) for name in range(10,61)]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as exc:
                print(f'Generated an exception: {exc}')


parallel_process(directory)


