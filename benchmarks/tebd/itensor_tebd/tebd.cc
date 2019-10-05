// this code is taken from http://itensor.org/docs.cgi?vers=cppv3&page=formulas/tevol_trotter
// (modified to run on transverse field Ising model)
#include "itensor/all.h"

using namespace itensor;

int 
main()
    {
      int N = 20; //number of sites
      Real tstep = 0.1; //time step (smaller is generally more accurate)
      Real ttotal = 5.; //total time to evolve
      Real cutoff = 1E-14; //truncation error cutoff when restoring MPS form
      Real mindim = 100;
      Real maxdim = 100;

      //Define a site set object "sites" which lets us
      //easily obtain Site indices defining our Hilbert space
      //and S=1/2 single-site operators
      auto sites = SpinHalf(N, {"ConserveQNs=",false});

      //Make initial MPS psi to be in the Neel state
      auto state = InitState(sites);
      for(auto j : range1(N))
        {
          state.set(j,"Up");
        }
      auto psi = MPS(state);

      //Create a std::vector (dynamically sizeable array)
      //to hold the Trotter gates
      auto gates = vector<BondGate>();

      //Create the gates exp(-i*tstep/2*hterm)
      //and add them to gates
      for(int b = 1; b <= N-1; ++b)
        {
          auto f1 = b==1 ? 1. : 0.5;
          auto f2 = b==N-1 ? 1. : 0.5;
          auto hterm = f1*0.5*op(sites,"Sx",b)*op(sites,"Id",b+1);
          hterm += f2*0.5*op(sites,"Sx",b+1)*op(sites,"Id",b);
          hterm += op(sites,"Sz",b)*op(sites,"Sz",b+1);

          auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
          gates.push_back(g);
        }
      //Create the gates exp(-i*tstep/2*hterm) in reverse
      //order (to get a second order Trotter breakup which
      //does a time step of "tstep") and add them to gates
      for(int b = N-1; b >= 1; --b)
        {
          auto f1 = b==1 ? 1. : 0.5;
          auto f2 = b==N-1 ? 1. : 0.5;
          auto hterm = f1*0.5*op(sites,"Sx",b)*op(sites,"Id",b+1);
          hterm += f2*0.5*op(sites,"Sx",b+1)*op(sites,"Id",b);

          hterm += op(sites,"Sz",b)*op(sites,"Sz",b+1);

          auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
          gates.push_back(g);
        }

      //Save initial state;
      auto psi0 = psi;

      //Time evolve, overwriting psi when done
      gateTEvol(gates,ttotal,tstep,psi,{"Cutoff=",cutoff,"Verbose=",true,"MaxDim",maxdim,"MinDim",mindim});

      printfln("Maximum MPS bond dimension after time evolution is %d",maxLinkDim(psi));

      //Print overlap of final state with initial state
      //(Will be complex so using innerC which can return complex);
      Print(innerC(psi,psi0));

      return 0;

    }



