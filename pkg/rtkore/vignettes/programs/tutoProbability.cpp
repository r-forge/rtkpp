#include "STKpp.h"
using namespace STK;
/** @ingroup tutorial */
int main(int argc, char** argv)
{
  Law::Normal l(2, 1);
  stk_cout << "l.cdf(3.96)= "   << l.cdf(3.96)  << _T("\n");
  stk_cout << "l.icdf(0.975)= " << l.icdf(0.975)  << _T("\n");
  stk_cout << "l.pdf(2)= "      << l.pdf(2)  << _T("\n");
  stk_cout << "l.lpdf(2)= "     << l.lpdf(2)  << _T("\n");
  stk_cout << "l.rand()= "      << l.rand()  << _T("\n");
  CArray33 a;
  a.rand(l);
  stk_cout << "a=\n"     << a  << _T("\n");
}
