import numpy as np
import os 
import sys
import numba 
import matplotlib.pyplot as plt
import iminuit 
import scipy
import pandas

from template_fit import *


template = LXeTemplate.load('lxe.template')
# template = LXeTemplate.load('meg.template')


df = []

for e1 in np.linspace(50,5000,20):
    print(e1)
    for e2 in np.linspace(50,5000,40):
        print("   ->", e2)
        for t1 in [100]:
            for t2 in np.linspace(50,150,50):
                for noisefloor in [1,5,10]:
                    times = [t1,t2]
                    energies = [e1,e2]
                    # print(times,energies)
                    data = make_fake_data_from_template(template.template,times, energies , noisefloor=noisefloor)
                    fit = TemplateFit(data, template.template, verbose=False, chi2limit=100, minimum_energy=25)
                    fit.do_fit()
                    dfi = {
                        't1':t1,
                        't2':t2,
                        'e1':e1,
                        'e2':e2,
                        'noise':noisefloor,
                        'nfit':fit.npulses,
                        'params':fit.current_guess,
                        'energies':[fit.current_guess[i*fit.npar] for i in range(int(len(fit.current_guess)/fit.npar))],
                        'times':[fit.current_guess[i*fit.npar+2] for i in range(int(len(fit.current_guess)/fit.npar))]
                    }
                    for i in range(int(len(fit.current_guess)/fit.npar)):
                        # print(i)
                        dfi[f'fit_t{i}'] = -1*fit.current_guess[i*fit.npar+2]
                        dfi[f'fit_e{i}'] = fit.current_guess[i*fit.npar+0]
                    df.append(dfi)

                    if(len(dfi['times']) > 2):
                        fit.plot()
                        plt.title(f"Times: {t1} / {t2} | Energies: {e1} / {e2}")
                        plt.tight_layout()
                        plt.savefig(f"./images/abnormal_fit/abnormal_fit_{t1}_{t2}_{e1}_{e2}_{noisefloor}.png", bbox_inches='tight', facecolor='w')
                        plt.close()
                        fit.save(f"./images/abnormal_fit/abnormal_fit_{t1}_{t2}_{e1}_{e2}_{noisefloor}.pickle")
                        # break
            #         break 
            #     break 
            # break

df = pandas.DataFrame(df)
df.head()


df['nfound'] = df['times'].map(len)
df.head()


df.to_csv('pulse_seperation_standalone.csv', sep='|')

