import pandas as pd
import sys

workdir = '/mnt/scratch/DMP/DUDMP/TRANSGEN/transgen-mdx/ngs/9.scripts/nextflow/nextflow_test/'

samplesheet = sys.argv[1]
print(samplesheet)

df = pd.read_csv(samplesheet, header=0, sep=',', comment='#', encoding='utf8')

samples = df['samp_id'].tolist()

for line in range(len(df)):
	sample = df.iloc[line]['samp_id']
	data = {
		'sample_id': sample,
		'pool_id':  df.iloc[line]['proj_id'],
		'run_id':  df.iloc[line]['run_id'],
		'seq_type': df.iloc[line]['seq_type'],
		'target_bed': df.iloc[line]['target_bed'],
		'gatk_grp': df.iloc[line]['gatk_grp'],
		'vpanel': df.iloc[line]['vpanel'],
		'umi': df.iloc[line]['umi'],
		'tag': df.iloc[line]['tag'],
		'tumour_type': df.iloc[line]['tumour_type']
		}
	
	config = pd.DataFrame.from_dict(data, orient='index')
	print(config)
	config.to_csv('.'+sample+'.config', header=False, index=True, sep='=')
	
