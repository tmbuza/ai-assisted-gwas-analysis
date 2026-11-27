"""Validate prepared quantitative-trait GWAS summaries without fabricating fields."""
from pathlib import Path
from decimal import Decimal,InvalidOperation
import csv,json,hashlib,gzip
import numpy as np
import pandas as pd

COLUMNS=['variant_id','chromosome','position','effect_allele','other_allele','beta',
         'standard_error','p_value','effect_allele_frequency','sample_size','imputation_info']

def sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()

def input_path(folder):
    folder=Path(folder)
    if (folder/'summary-statistics.tsv').exists():return folder/'summary-statistics.tsv'
    return folder/'summary-statistics.tsv.gz'

def open_text(path):return gzip.open(path,'rt',newline='') if str(path).endswith('.gz') else Path(path).open(newline='')

def negative_log10_p(values):
    # Preserve source strings; Decimal handles probabilities below binary64 range.
    numeric=pd.to_numeric(values,errors='coerce').to_numpy(dtype=float)
    if np.any(~np.isfinite(numeric)) or np.any(numeric>1) or np.any(numeric<0):raise ValueError('Invalid source p-values.')
    with np.errstate(divide='ignore'):result=-np.log10(numeric)
    for i in np.flatnonzero(numeric<1e-250):
        try:
            value=Decimal(str(values.iloc[i]))
            if not value.is_finite() or not 0<value<=1:raise ValueError('p_value must be positive; literal zero cannot be recovered automatically.')
            result[i]=float(-value.log10())
        except InvalidOperation as e:raise ValueError('Invalid source p-value text.') from e
    if not np.isfinite(result).all():raise ValueError('Nonfinite association strength.')
    return result

def load_inputs(folder):
    folder=Path(folder);path=input_path(folder)
    metadata=json.loads((folder/'study-metadata.json').read_text())
    for c in ['study_id','trait','population','genome_build','effect_type','effect_units','summary_scope','source','synthetic']:
        if c not in metadata or metadata[c] is None or str(metadata[c]).strip()=='':raise ValueError(f'Missing study metadata: {c}')
    if not isinstance(metadata['synthetic'],bool):raise ValueError('synthetic must be a JSON boolean.')
    if metadata['effect_type']!='beta':raise ValueError('Quantitative-trait beta only; do not substitute odds ratios or log-odds.')
    if metadata['genome_build'] not in ['GRCh37','GRCh38','SYNTHETIC']:raise ValueError('Declare a supported genome build.')
    if (metadata['genome_build']=='SYNTHETIC')!=metadata['synthetic']:raise ValueError('Synthetic flag and build disagree.')
    if metadata['summary_scope'] not in ['genome_wide','selected_variants']:raise ValueError('Invalid summary_scope.')
    info_available=metadata.get('info_available',True)
    if not isinstance(info_available,bool):raise ValueError('info_available must be a boolean.')
    required=COLUMNS if info_available else COLUMNS[:-1]
    with open_text(path) as handle:header=next(csv.reader(handle,delimiter='\t'),[])
    if len(header)!=len(set(header)) or not set(required)<=set(header):raise ValueError('Missing or duplicate headers.')
    t=pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False)
    if t.empty or t[required].apply(lambda s:s.str.strip().eq('')).any().any():raise ValueError('Empty table or blank required fields.')
    if t.variant_id.duplicated().any() or (t.variant_id!=t.variant_id.str.strip()).any():raise ValueError('Duplicate or padded variant IDs.')
    t['source_minus_log10_p']=negative_log10_p(t.p_value)
    for c in ['chromosome','position','beta','standard_error','effect_allele_frequency','sample_size']:
        t[c]=pd.to_numeric(t[c],errors='raise')
        if not np.isfinite(t[c]).all():raise ValueError(f'{c}: nonfinite values.')
    for c in ['chromosome','position']:
        if ((t[c]%1)!=0).any():raise ValueError(f'{c}: expected integers.')
        t[c]=t[c].astype('int64')
    if metadata.get('sample_size_type','integer')=='integer' and ((t.sample_size%1)!=0).any():raise ValueError('N must be integer unless source_reported is declared.')
    if not t.chromosome.between(1,22).all() or not (t.position>0).all():raise ValueError('Autosomes 1–22 and positive positions only.')
    if not (t.sample_size>0).all() or not (t.standard_error>0).all():raise ValueError('N and SE must be positive.')
    if not t.effect_allele_frequency.between(0,1).all():raise ValueError('Frequency outside [0,1].')
    if info_available:
        t['imputation_info']=pd.to_numeric(t.imputation_info,errors='raise')
        if not np.isfinite(t.imputation_info).all() or not t.imputation_info.between(0,1).all():raise ValueError('Invalid INFO.')
    else:
        if 'imputation_info' in t and not t.imputation_info.isin(['','NA']).all():raise ValueError('INFO declared unavailable but numerical scores supplied.')
        t['imputation_info']=np.nan
    if not t.effect_allele.isin(list('ACGT')).all() or not t.other_allele.isin(list('ACGT')).all() or (t.effect_allele==t.other_allele).any():raise ValueError('Distinct uppercase A/C/G/T SNVs only.')
    if t.duplicated(['chromosome','position']).any():raise ValueError('Reconcile multiple rows per coordinate upstream.')
    return t,metadata
