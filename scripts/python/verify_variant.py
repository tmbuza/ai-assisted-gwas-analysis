"""Independently recalculate one variant using Python's standard library."""
import argparse,csv,json,math
from pathlib import Path
from statistics import NormalDist


def get_row(path,variant):
    with path.open() as handle:
        matches=[r for r in csv.DictReader(handle,delimiter='\t') if r['variant_id']==variant]
    if len(matches)!=1:raise ValueError(f'Expected exactly one matching row in {path}; found {len(matches)}.')
    return matches[0]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('data_dir',type=Path);p.add_argument('output_dir',type=Path);p.add_argument('variant_id');a=p.parse_args()
    try:
        raw=get_row(a.data_dir/'summary-statistics.tsv',a.variant_id)
        saved=get_row(a.output_dir/'tables/variant-results.tsv',a.variant_id)
        beta=float(raw['beta']);se=float(raw['standard_error']);z=beta/se
        ci=NormalDist().inv_cdf(.975)*se
        checks={'z_from_beta_se':z,'wald_p_from_beta_se':math.erfc(abs(z)/math.sqrt(2)),
                'minus_log10_p':-math.log10(float(raw['p_value'])),'ci_lower':beta-ci,'ci_upper':beta+ci,
                'maf':min(float(raw['effect_allele_frequency']),1-float(raw['effect_allele_frequency']))}
        for key,value in checks.items():
            recorded=float(saved[key]);agrees=math.isclose(value,recorded,rel_tol=1e-9,abs_tol=1e-300)
            print(f'{key}: recalculated={value:.12g}; saved={recorded:.12g}; agrees={agrees}')
            if not agrees:raise ValueError(f'{key} disagrees.')
        print('Variant verification passed. Source p-values can differ from Wald p-values for other study models.')
    except (OSError,ValueError,KeyError) as error:p.exit(1,f'Verification failed: {error}\n')

if __name__=='__main__':main()
