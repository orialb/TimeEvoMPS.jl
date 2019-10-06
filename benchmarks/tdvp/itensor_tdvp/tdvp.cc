// this code is based on the tdvp sample in the "tdvp" barnch of itensor
// (modified to run on transverse field Ising model)
#include "itensor/all.h"

using namespace itensor;

int
main()
    {
      int N = 20; //number of sites
      Real tstep = 0.1; //time step (smaller is generally more accurate)
      Real ttotal = 20.; //total time to evolve
      Real cutoff = 1E-14; //truncation error cutoff when restoring MPS form
      Real mindim = 100;
      Real maxdim = 100;

      //Define a site set object "sites" which lets us
      //easily obtain Site indices defining our Hilbert space
      //and S=1/2 single-site operators
      auto sites = SpinHalf(N, {"ConserveQNs=",false});

      auto state = InitState(sites);
      for(auto j : range1(N))
        {
          state.set(j,"Up");
        }
      auto psi = MPS(state);


      //Save initial state;
      auto psi0 = psi;

      auto ampo = AutoMPO(sites);

      // chain
      for(int j = 1; j < N; ++j)
        {
          ampo += -1,"Sz",j,"Sz",j+1;
          ampo +=  -0.5,"Sx",j;
        }
      ampo += -0.5,"Sx",N;

      auto H = toMPO(ampo);


      const int nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));
      if(fabs(nt*tstep-ttotal) > 1E-9)
        {
          Error("Timestep not commensurate with total time");
        }

      auto sweeps = Sweeps(nt);
      sweeps.maxdim() = maxdim;
      sweeps.cutoff() = cutoff;
      sweeps.niter() = 30;
      println(sweeps);

      auto energy = tdvp(psi0,H,tstep,sweeps,{"DoNormalize",true,"UseLanczos",true,"Ideg",6,
                                              "IsRealTevol",true,"Quiet",true,"NumCenter",2});

      printfln("Maximum MPS bond dimension after time evolution is %d",maxLinkDim(psi0));

      return 0;

    }



