import tenpy
from tenpy.networks.mps import MPS
from tenpy.models.tf_ising import TFIChain
from tenpy.algorithms import tebd
import numpy as np


def run_tebd_tfi(L,mindim,g, tf,dt,verbose = False,order=2):
    model_params = dict(L=L, J=1., g=g, bc_MPS='finite', conserve=None, verbose=verbose)
    M = TFIChain(model_params)
    product_state = ["up"] * M.lat.N_sites
    psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)
    tebd_params = {
        'order': order,
        'dt': dt,
        'N_steps': int(tf/dt),
        'trunc_params': {
            'trunc_cut': 1e-14,
            'svd_min': 1e-24,
            'chi_max': mindim
            # 'chi_min': mindim
        },
        'verbose': verbose,
    }
    eng = tebd.Engine(psi,M,tebd_params)
    eng.run()
    print(np.amax(psi.chi))
    return psi

