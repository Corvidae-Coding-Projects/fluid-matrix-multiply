"""Independent checks and compact, explicitly sampled CFD data for the viewer."""
from pathlib import Path
import json
import numpy as np

ROOT = Path(__file__).resolve().parent
R = ROOT / 'results'

def meta(name):
    return json.loads((R / f'{name}.json').read_text())

def field(name):
    return np.fromfile(R / f'{name}-field.bin', dtype='<f4').reshape(-1, 7)

def analyze():
    names = ['e0-s1', 'e1-s1', 'e2-s1', 'e3-s1', 'mix-s1', 'signed-s1', 'mix-s2']
    cases = {n: meta(n) for n in names}
    matrix = np.array([np.array(cases[f'e{i}-s1']['output_mass_equivalent_nLs']) / cases[f'e{i}-s1']['pressure_Pa'] for i in range(4)]).T
    validation = {'matrix_nLs_per_Pa': matrix.tolist(), 'checks': {}}
    for name in ['mix-s1', 'signed-s1']:
        d = cases[name]
        predicted = matrix @ np.array(d['inputs']) * d['pressure_Pa']
        actual = np.array(d['output_mass_equivalent_nLs'])
        validation['checks'][name] = {
            'predicted_nLs': predicted.tolist(), 'computed_nLs': actual.tolist(),
            'relative_L2_superposition_error': float(np.linalg.norm(actual-predicted)/np.linalg.norm(actual)),
            'max_output_relative_error': float(np.max(np.abs((actual-predicted)/actual)))
        }
    q1 = np.array(cases['mix-s1']['output_mass_equivalent_nLs'])
    q2 = np.array(cases['mix-s2']['output_mass_equivalent_nLs'])
    validation['grid_refinement'] = {
        'coarse_dx_um': 10, 'fine_dx_um': 5,
        'relative_L2_change': float(np.linalg.norm(q2-q1)/np.linalg.norm(q2)),
        'per_output_relative_change': ((q2-q1)/q2).tolist(),
        'note': 'Two grids measure sensitivity, not a proven continuum error bound.'
    }
    ducts = []
    coeff = (1-192/np.pi**5*sum(np.tanh(n*np.pi/2)/n**5 for n in range(1,1000,2)))/12
    for scale in [1,2]:
        name = f'duct-s{scale}'; d = meta(name); f = field(name)
        xs = np.unique(f[:,0]); p = np.array([f[f[:,0]==x,3].mean(dtype=np.float64) for x in xs])
        interior = (xs>150)&(xs<490)
        gradient = -np.polyfit(xs[interior]*1e-6, p[interior], 1)[0]
        exact = coeff*(80e-6)**4/.001*gradient*1e12
        # Check volumetric, not mass-equivalent, axial flux in the center section.
        center = f[f[:,0]==xs[len(xs)//2]]
        volume = center[:,4].sum(dtype=np.float64)*1e-3*d['dx_m']**2*1e12
        ducts.append({'scale':scale, 'analytical_nLs':float(exact), 'computed_volumetric_nLs':float(volume),
                      'relative_error':float(abs(volume/exact-1)), 'converged':d['converged']})
    validation['square_duct'] = ducts
    zero = meta('zero-s1')
    validation['zero_input_max_output_nLs'] = max(map(abs,zero['output_mass_equivalent_nLs']))
    validation['all_device_runs_converged'] = all(d['converged'] for d in cases.values())
    validation['maximum_device_mass_imbalance'] = max(d['mass_balance_relative'] for d in cases.values())
    validation['maximum_device_mach'] = max(d['max_mach'] for d in cases.values())
    validation['numerical_tolerances'] = {'device_velocity_relative_change':1e-6, 'duct_velocity_relative_change':1e-8,
                                         'consecutive_checks_required':3, 'physical_check_interval_s':.002}
    (R/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    # Regular subsampling retains real volume-cell values. It does not construct
    # synthetic pipe profiles, interpolate between CFD cases, or animate flow.
    packed = {'matrix':np.round(matrix,8).tolist(), 'validation':validation, 'cases':{},
              'widths_um':(np.loadtxt(ROOT/'widths.txt')*10).tolist()}
    labels = {'mix-s1':'Mixed inputs · 10 µm grid','signed-s1':'Signed inputs · 10 µm grid',
              'mix-s2':'Mixed inputs · 5 µm grid', **{f'e{i}-s1':f'Input {i+1} only · 10 µm grid' for i in range(4)}}
    # Data budgets: ~4,000 volume samples per case, including all 3 velocity components.
    for name in names:
        d=cases[name]; f=field(name); scale=d['scale']
        ijk=np.floor(f[:,:3]/(10/scale)).astype(int)
        stride=3*scale
        keep=np.all(ijk%stride == stride//2,axis=1)
        a=f[keep].copy()
        a[:,:3]=np.round(a[:,:3],1);a[:,3:]=np.round(a[:,3:],5)
        packed['cases'][name]={'label':labels[name], 'meta':d,'samples':[[round(float(v),1 if k<3 else 4) for k,v in enumerate(row)] for row in a]}
    (R/'viewer-data.json').write_text(json.dumps(packed,separators=(',',':')))
    print(json.dumps(validation,indent=2))
    # These gates test real risks; mesh sensitivity is reported, not passed off as accuracy.
    assert validation['all_device_runs_converged'], 'An unconverged CFD run remains'
    assert validation['maximum_device_mass_imbalance'] < 1e-4, 'Device mass conservation failure'
    assert validation['maximum_device_mach'] < .03, 'Low-Mach assumption failure'
    assert ducts[1]['relative_error'] < .02 and ducts[1]['relative_error'] < ducts[0]['relative_error'], 'Duct benchmark failure'
    assert validation['zero_input_max_output_nLs'] < 1e-8, 'Spurious zero-input flow'
    assert validation['checks']['mix-s1']['relative_L2_superposition_error'] < .01, 'Multiplier is too nonlinear in this operating range'
    assert validation['checks']['signed-s1']['relative_L2_superposition_error'] < .01, 'Signed-pressure superposition failure'

if __name__ == '__main__':
    analyze()
