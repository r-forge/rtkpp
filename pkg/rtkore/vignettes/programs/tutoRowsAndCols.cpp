#include "STKpp.h"
using namespace STK;
int main(int argc, char *argv[])
{
  ArrayXX A(4, 2);
  A << 1, 2,   1, 3,   2, 3,   2, 4;
  stk_cout << _T("A =\n") << A;
  // Adding a column with 1
  A.pushFrontCols(Const::VectorX(4));
  // Insert x^2 and x^3
  A.insertCols(2, 1);
  A.col(2) = A.col(1).square();
  A.pushBackRows(1);
  A.row(4) = A.row(0).square();
  stk_cout << _T("A =\n") << A;
}
