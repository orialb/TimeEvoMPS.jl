import tenpy
from tenpy.networks.mps import MPS
from tenpy.models.tf_ising import TFIChain
from tenpy.algorithms import tdvp
import numpy as np

def run_tenpy_tdvp(L,maxdim,tf,verbose=False):
    model_params = dict(L=L, J=1., g=1, bc_MPS='finite', conserve=None, verbose=verbose)
    M = TFIChain(model_params)
    product_state = ["up"] * M.lat.N_sites
    psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

    tdvp_params = {
        'start_time': 0,
        'dt': 0.1,
        'trunc_params': {
            'chi_max': maxdim,
            'svd_min': 0,
            'trunc_cut': None
        }
    }
    tdvp_engine = tdvp.Engine(psi=psi, model=M, TDVP_params=tdvp_params)
    tdvp_engine.run_two_sites(N_steps=int(tf/0.1))

    return np.amax(psi.chi)

if __name__ == "__main__":
    run_tenpy_tdvp()
