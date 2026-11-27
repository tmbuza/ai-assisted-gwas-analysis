"""Reproduce prepared full GLGC 2013 HDL summaries from a checksum-pinned source."""
from pathlib import Path
import argparse,json,urllib.request
import pandas as pd
import numpy as np
from common import sha256,negative_log10_p
ROOT=Path(__file__).resolve().parents[2]
URL='https://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_HDL.txt.gz'
SHA='7d03f9691db53ddd78d8545f85ebd11bf8a4d6686e334e44509b1ec447355278'
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',type=Path);a=ap.parse_args()
    raw=a.source or ROOT/'source/glgc2013/jointGwasMc_HDL.txt.gz';raw.parent.mkdir(parents=True,exist_ok=True)
    if not raw.exists():
        partial=raw.with_suffix('.partial');print('Downloading verified 59 MB source release.',flush=True)
        with urllib.request.urlopen(URL,timeout=120) as response,partial.open('wb') as out:
            for block in iter(lambda:response.read(1024*1024),b''):out.write(block)
        partial.replace(raw)
    if sha256(raw)!=SHA:raise ValueError('Source checksum differs; investigate before accepting a changed source.')
    parts=[];excluded=[];total=0
    for c in pd.read_csv(raw,sep='\t',dtype=str,keep_default_na=False,chunksize=200000):
        total+=len(c);loc=c.SNP_hg19.str.extract(r'^chr([0-9]+):([0-9]+)$')
        chrom=pd.to_numeric(loc[0],errors='coerce');pos=pd.to_numeric(loc[1],errors='coerce')
        freq=pd.to_numeric(c['Freq.A1.1000G.EUR'],errors='coerce')
        invalid=~chrom.between(1,22)|pos.isna()|~freq.between(0,1)
        if invalid.any():
            ex=c.loc[invalid,['SNP_hg19','rsid']].copy();ex['reason']=np.where(freq.loc[invalid].isna(),'reference_frequency_unavailable','unsupported_coordinate_or_frequency');excluded.append(ex)
        c=c.loc[~invalid].copy();ea=c.A1.str.upper();oa=c.A2.str.upper();ids=c.rsid.copy();missing=ids.isin(['NA','','.'])
        ids.loc[missing]=c.loc[missing,'SNP_hg19']+':'+ea.loc[missing]+':'+oa.loc[missing]
        out=pd.DataFrame(dict(variant_id=ids,chromosome=chrom.loc[c.index].astype(int),position=pos.loc[c.index].astype(int),
            effect_allele=ea,other_allele=oa,beta=c.beta,standard_error=c.se,p_value=c['P-value'],
            effect_allele_frequency=freq.loc[c.index],sample_size=c.N))
        negative_log10_p(out.p_value.reset_index(drop=True));parts.append(out)
    prepared=pd.concat(parts,ignore_index=True)
    if prepared.variant_id.duplicated().any() or prepared.duplicated(['chromosome','position']).any():raise ValueError('Unresolved duplicate IDs or coordinates.')
    prepared=prepared.sort_values(['chromosome','position']).reset_index(drop=True)
    dest=ROOT/'data/glgc-hdl';dest.mkdir(parents=True,exist_ok=True)
    prepared.to_csv(dest/'summary-statistics.tsv.gz',sep='\t',index=False,compression={'method':'gzip','compresslevel':6})
    meta=dict(study_id='GLGC_2013_HDL_joint',trait='HDL cholesterol',population='European-ancestry discovery meta-analysis',
        genome_build='GRCh37',effect_type='beta',effect_units='inverse-normal-transformed HDL units per A1 allele',summary_scope='genome_wide',
        source=URL,synthetic=False,info_available=False,sample_size_type='source_reported',frequency_basis='1000 Genomes EUR reference',
        source_p_method='Stouffer meta-analysis with genomic control; not beta/SE Wald p-values',
        effect_coding='A1 is the trait-increasing allele in this release; nonnegative betas reflect source coding.',
        source_sha256=SHA,publication_doi='10.1038/ng.2797',original_source_rows=total,
        preparation='Full autosomal release with available reference frequency. No p-value selection or chromosome subsampling. INFO absent, filter skipped. Source N decimals and p strings preserved.')
    (dest/'study-metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    prov=ROOT/'data/provenance/glgc-hdl';prov.mkdir(parents=True,exist_ok=True)
    exclusions=pd.concat(excluded,ignore_index=True) if excluded else pd.DataFrame(columns=['SNP_hg19','rsid','reason'])
    exclusions.to_csv(prov/'preparation-exclusions.tsv.gz',sep='\t',index=False,compression='gzip')
    audit=dict(source_url=URL,source_sha256=SHA,source_rows=total,prepared_rows=len(prepared),excluded_rows=len(exclusions),
        exclusion_counts=exclusions.reason.value_counts().to_dict(),derived_coordinate_ids=int(prepared.variant_id.str.startswith('chr').sum()),
        prepared_sha256=sha256(dest/'summary-statistics.tsv.gz'),source_p_strings_preserved=True,info_scores_fabricated=False,reported_N_rounded=False)
    (prov/'preparation-audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit),flush=True)
if __name__=='__main__':main()
