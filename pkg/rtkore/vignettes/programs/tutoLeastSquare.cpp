#include <STKpp.h>
using namespace STK;
/** @ingroup tutorial */
int main(int argc, char** argv)
{
  CArrayXX y(1000,3), x(1000,5), beta(5,3);
  Law::Normal l(0,1);
  x.rand(l);
  beta << 0, 1, 2,
          2, 3, 4,
          2, 1, 6,
          0, 3,-1,
          3,-1, 1;
  y = x * beta + CArrayXX(1000, 3).rand(l);
  stk_cout << _T("STK++ MultiLeastSquare:\n");
  stk_cout << _T("-----------------------\n");
  STK::MultiLeastSquare<CArrayXX, CArrayXX> ols(y,x); ols.run();
  stk_cout << _T("beta matrix:\n");
  stk_cout << ols.x();
  stk_cout << _T("lapack MultiLeastSquare:\n");
  stk_cout << _T("-----------------------\n");
  STK::lapack::MultiLeastSquare<CArrayXX, CArrayXX> ls(y,x); ls.run();
  stk_cout << _T("beta matrix:\n");
  stk_cout << ls.x();
  return 0;
}
