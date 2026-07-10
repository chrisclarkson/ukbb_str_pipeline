# python get_coverage_para.py /mnt/project/cc_analysis/output_dir/EHv5/

import glob
import sys
import multiprocessing
import json
gene=sys.argv[1]
try:
    directory=sys.argv[2]
except:
    directory='recalled_eh_bams'


def to_file(v,gene,coverage,name):
    with(open('coverages_'+name+'.txt','a')) as out:
        out.write(v+'\t'+gene+'\t'+str(coverage)+'\n')

def parse(files,name,gene):
    for v in files:
        try:
            with(open(v,'r'))  as f:
                j=json.load(f)
                coverage=j['LocusResults'][gene]['Coverage']
                to_file(v,gene,coverage,str(name))
        except:
            print('Something wrong with {0}'.format(v))


import concurrent.futures

def worker(name, directory):
    print('{0}{1}/*.json'.format(directory, str(name)))
    files = glob.glob('{0}{1}/*.json'.format(directory, str(name)))
    print(files[0:4])
    print('files loaded for {0}'.format(name))
    parse(files, str(name),gene)

def parallel_process(directory):
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(worker, name, directory) for name in range(10,61)]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as exc:
                print(f'Generated an exception: {exc}')


parallel_process(directory)


