"""One strict TSV and metadata contract for the GWAS teaching workflow."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
import pandas as pd

COLUMNS=['variant_id','chromosome','position','effect_allele','other_allele','beta',
         'standard_error','p_value','effect_allele_frequency','sample_size','imputation_info']

def sha256(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_inputs(folder):
    folder=Path(folder);path=folder/'summary-statistics.tsv'
    with path.open(newline='') as handle:header=next(csv.reader(handle,delimiter='\t'),[])
    if len(header)!=len(set(header)) or not set(COLUMNS)<=set(header):
        raise ValueError('Expected unique headers and all columns specified in the guide.')
    t=pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False)
    if t.empty or t[COLUMNS].apply(lambda x:x.str.strip().eq('')).any().any():raise ValueError('Empty table or blank required cells.')
    if t.variant_id.duplicated().any() or (t.variant_id!=t.variant_id.str.strip()).any():raise ValueError('Duplicate or padded variant IDs.')
    for c in COLUMNS[1:]:
        if c not in ['effect_allele','other_allele']:
            t[c]=pd.to_numeric(t[c],errors='raise')
            if not np.isfinite(t[c]).all():raise ValueError(f'{c}: nonfinite values.')
    for c in ['chromosome','position','sample_size']:
        if ((t[c]%1)!=0).any():raise ValueError(f'{c}: expected integers.')
    if not t.chromosome.between(1,22).all() or not (t.position>0).all():raise ValueError('Only autosomes 1–22 and positive positions are supported.')
    if not (t.sample_size>0).all() or not (t.standard_error>0).all():raise ValueError('N and standard errors must be positive.')
    if not ((t.p_value>0)&(t.p_value<=1)).all():raise ValueError('p_value must satisfy 0 < p <= 1. Recover underflowed p-values upstream; do not replace zero silently.')
    for c in ['effect_allele_frequency','imputation_info']:
        if not t[c].between(0,1).all():raise ValueError(f'{c}: expected values in [0,1].')
    if not t.effect_allele.isin(list('ACGT')).all() or not t.other_allele.isin(list('ACGT')).all() or (t.effect_allele==t.other_allele).any():
        raise ValueError('This starter supports distinct uppercase single-base A/C/G/T allele pairs only.')
    if t.duplicated(['chromosome','position','effect_allele','other_allele']).any():raise ValueError('Duplicate coordinate/allele records.')
    if t.duplicated(['chromosome','position']).any():raise ValueError('This starter expects one biallelic record per position; reconcile multiallelic records upstream.')
    metadata=json.loads((folder/'study-metadata.json').read_text())
    for c in ['study_id','trait','population','genome_build','effect_type','effect_units','summary_scope','source','synthetic']:
        if c not in metadata or metadata[c] is None or str(metadata[c]).strip()=='':raise ValueError(f'Missing study metadata: {c}')
    if not isinstance(metadata['synthetic'],bool):raise ValueError('synthetic must be a JSON boolean.')
    if metadata['effect_type']!='beta':raise ValueError('Starter accepts quantitative-trait beta only. Do not pass odds ratios or log-odds without extending and documenting the model.')
    if metadata['genome_build'] not in ['GRCh37','GRCh38','SYNTHETIC']:raise ValueError('Use GRCh37, GRCh38, or SYNTHETIC genome_build.')
    if (metadata['genome_build']=='SYNTHETIC') != metadata['synthetic']:raise ValueError('Synthetic flag and genome build disagree.')
    if metadata['summary_scope'] not in ['genome_wide','selected_variants']:raise ValueError('summary_scope must be genome_wide or selected_variants.')
    return t,metadata
