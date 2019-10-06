// this code is based on the tdvp sample in the "tdvp" barnch of itensor
// (modified to run on transverse field Ising model)
#include "itensor/all.h"
#include "mylocalop.h"

int
main()
    {
      Real m = 100;

      auto i = Index(m);
      auto j = Index(m);
      auto A = randomITensor(prime(i),i);
      auto B = randomITensor(prime(j),j);
      A /= norm(A);
      B /= norm(B);
      auto M = MyLocalOp(A,B);
      auto v = randomITensor(i,j);
      auto v1 = v;


      for(int n=1; n<=1000; ++n)
        {
          auto v = randomITensor(i,j);
          applyExp(M,v,0.1,{"MaxIter", 30,  "ErrGoal",1E-12});
        }

      return 0;

    }



