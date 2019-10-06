#include "itensor/all.h"

using namespace itensor;

// A quick hack for testing purposes
// this represents an operator which acts on a two-legged tensor phi, as
//  -- A -- phi -- B --
// I implemented it to be able to easily build a BigMatrixT object to pass to
// the Krylov exponentiation method.
class MyLocalOp
{
  ITensor const* A_;
  ITensor const* B_;

 public:

	MyLocalOp(ITensor const& A, ITensor const& B );


	void
    update(ITensor const& A, ITensor const& B);

  ITensor const&
    A() const
  {
    return *A_;
  }

  ITensor const&
    B() const
  {
    return *B_;
  }


  unsigned long
    size() const;

  void
    product(ITensor const& phi, ITensor & phip) const;
};

inline MyLocalOp::
MyLocalOp(const ITensor& A, const ITensor& B)
{
  update(A,B);
}

void inline MyLocalOp::
update(const ITensor& A, const ITensor& B)
{
  A_ = &A;
  B_ = &B;
}

void inline MyLocalOp::
product(ITensor const& phi, ITensor & phip) const
{
  phip = phi;
  phip *= A();
  phip *= B();
  phip.noPrime();
}

unsigned long inline MyLocalOp::
size() const
{
  return dim(A_ -> inds()(1))*dim(B_ -> inds()(1));
}

