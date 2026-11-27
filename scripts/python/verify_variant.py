"""Independently check one variant, including extreme source p-values, with stdlib."""
import argparse,csv,math
from decimal import Decimal
from pathlib import Path
from statistics import NormalDist
from common import open_text,input_path

def get_row(path,variant):
    found=None
    with open_text(path) as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            if row['variant_id']==variant:
                if found is not None:raise ValueError('Duplicate verification identifier.')
                found=row
    if found is None:raise ValueError(f'Variant {variant} not found in {path}.')
    return found

def normal_log10_p(z):
    x=abs(z)
    if x<10:return -math.log10(math.erfc(x/math.sqrt(2)))
    # Mills-ratio asymptotic expansion avoids normal-tail underflow.
    term=1.;series=1.
    for k in range(1,21):term*=-(2*k-1)/(x*x);series+=term
    logp=math.log(2)-x*x/2-math.log(x)-.5*math.log(2*math.pi)+math.log(series)
    return -logp/math.log(10)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('data_dir',type=Path);parser.add_argument('output_dir',type=Path);parser.add_argument('variant_id');a=parser.parse_args()
    try:
        raw=get_row(input_path(a.data_dir),a.variant_id)
        path=a.output_dir/'tables/variant-results.tsv'
        if not path.exists():path=path.with_suffix('.tsv.gz')
        saved=get_row(path,a.variant_id)
        beta=float(raw['beta']);se=float(raw['standard_error']);z=beta/se;interval=NormalDist().inv_cdf(.975)*se
        checks={'z_from_beta_se':z,'minus_log10_p':float(-Decimal(raw['p_value']).log10()),
                'minus_log10_wald_p':normal_log10_p(z),'ci_lower':beta-interval,'ci_upper':beta+interval,
                'maf':min(float(raw['effect_allele_frequency']),1-float(raw['effect_allele_frequency']))}
        for name,manual in checks.items():
            recorded=float(saved[name]);agrees=math.isclose(manual,recorded,rel_tol=1e-9,abs_tol=1e-12)
            print(f'{name}: recalculated={manual:.12g}; saved={recorded:.12g}; agrees={agrees}',flush=True)
            if not agrees:raise ValueError(f'{name} disagrees.')
        print('Source p-value text:',raw['p_value']);print('Variant verification passed. Source meta-analysis p-values need not equal derived Wald values.')
    except (OSError,ValueError,KeyError) as e:parser.exit(1,f'Verification failed: {e}\n')
if __name__=='__main__':main()
