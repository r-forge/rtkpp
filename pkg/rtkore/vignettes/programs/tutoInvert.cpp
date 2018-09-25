#include <STKpp.h>
using namespace STK;
/** @ingroup tutorial */
int main(int argc, char** argv)
{
  CArray<Real, 5, 5> a(5,5);
  CArraySquare<Real, 5> b;
  a << 0, 1, 2, 3, 4,
       2, 3, 4, 5, 6,
       2, 1, 6, 7, 8,
       0, 3,-1, 2, 3,
       3,-1, 1, 1, 2;
  stk_cout << _T("Inverse general matrix:\n");
  stk_cout << _T("-----------------------\n");
  stk_cout << a*invert(a);
  stk_cout << _T("\nInverse upper-symmetric matrix:\n");
  stk_cout << _T("-------------------------------\n");
  stk_cout << a.upperSymmetrize();
  stk_cout << a.upperSymmetrize()*invert(a.upperSymmetrize());
  stk_cout << _T("\nInverse lower-symmetric matrix:\n");
  stk_cout << _T("-------------------------------\n");
  stk_cout << a.lowerSymmetrize();
  stk_cout << a.lowerSymmetrize()*invert(a.lowerSymmetrize());
  return 0;
}
