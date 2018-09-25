#include <STKpp.h>
using namespace STK;
/** @ingroup tutorial */
int main(int argc, char** argv)
{
  Array2DLowerTriangular<Real> a(5,5);
  a << 1,
       1, 2,
       3, 4, 3,
       4, 5, 6, 6,
       7, 8, 2, 3, 2;
  stk_cout << _T("\nInverse lower-triangular matrix:\n");
  stk_cout << _T("----------------------------------\n");
  stk_cout << a*invert(a);
  stk_cout << _T("\nInverse upper-triangular matrix:\n");
  stk_cout << _T("----------------------------------\n");
  stk_cout << a.transpose()*invert(a.transpose());
  return 0;
}
