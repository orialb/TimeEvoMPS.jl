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

      // exponentiate multiple times to make sure that the total runtime is completely
      // dominated by applyExp (I run `time ./exp` and divide by 1000)
      for(int n=1; n<=1000; ++n)
        {
          // applyExp updates the vector in-place so we need to
          // initialize it every time.
          auto v1 = v;
          applyExp(M,v1,0.1,{"MaxIter", 30,  "ErrGoal",1E-12});
        }

      return 0;

    }



