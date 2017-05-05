#include "STKpp.h"
using namespace STK;
int main(int argc, char** argv)
{
  CArray<Real, 5, 8> A;
  CArrayPoint<Real, 8> mu, std, aux;
  Law::Normal law(1,2);
  A.rand(law);
  stk_cout << _T("mean(A) =") << (mu=Stat::mean(A));
  stk_cout << _T("variance(A) =") << Stat::varianceWithFixedMean(A, mu, false) << _T("\n\n");
  Stat::standardize(A, mu, std);
  stk_cout << _T("mean(A) =") << (aux=Stat::mean(A));
  stk_cout << _T("variance(A) =") << Stat::varianceWithFixedMean(A, aux, false) << _T("\n\n");
  Stat::unstandardize(A, mu, std);
  stk_cout << _T("mean(A) =") << (aux=Stat::mean(A));
  stk_cout << _T("variance(A) =") << Stat::varianceWithFixedMean(A, aux, false) << _T("\n\n");

  return 0;
}
