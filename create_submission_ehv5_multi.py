#chunker for ehv5_multi in ukbb
import sys
import math
 
input_file = sys.argv[1] # csv file containing sample EID, path to the cram file and sex of sample
batch_size = int(sys.argv[2])
project = "Matteo_Pietro" # change study name if required

def _parse_input(line, project):
    '''parse input_file into ID,readspath, and sex'''
    split_line = line.split(',')
    id = split_line[0]
    path = project + ":" + split_line[2].strip().replace(' ','\ ')
    sex = split_line[1].strip()
    return id, path, sex

fd = open(input_file)
lines = fd.readlines()

# sample_number=200000
input_number = 0
sample_number = len(lines)
number_of_batch = int(math.ceil(sample_number * 1.0 / batch_size))

for batch_number in range(number_of_batch):
    batch_crams = ''
    batch_sex = ''
    batch_output = ''
    folder = 'output_dir/EHv5/' # change output folder if required
    for member in range(batch_size):
        delim_line = lines[input_number]
        id, path, sex = _parse_input(delim_line, project)
        batch_crams += '-ireads="{}" '.format(path)
        batch_sex += '-isex={} '.format(sex)
        batch_output += '-ioutput_prefix={} '.format(id)
        input_number += 1
        if input_number == sample_number:
            break
    '''change reference and catalog if required'''
    print('dx run ehv5_multi_original {batch_crams}{batch_sex}{batch_output} \
-ireference=GRCh38_full_analysis_set_plus_decoy_hla.fa \
-ivariant_catalog=variant_catalog_final.json \
-ithreads=16 \
-ianalysis_mode=streaming \
-ioutput_folder={folder} -y --brief --priority low'.format(
          folder=folder,
          batch_crams=batch_crams,
          batch_sex=batch_sex,
          batch_output=batch_output))
          