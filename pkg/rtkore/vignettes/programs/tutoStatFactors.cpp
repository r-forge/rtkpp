#include "STKpp.h"
using namespace STK;
int main(int argc, char** argv)
{
  CArray<char, 13, 3> datai;
  datai << 'b', 'a', 'a',   'c', 'b', 'a',   'b', 'a', 'a',   'c', 'a', 'a',   'd', 'a', 'a',
           'b', 'b', 'a',   'd', 'c', 'b',   'c', 'c', 'b',   'b', 'c', 'b',   'd', 'b', 'b',
           'b', 'a', 'b',   'd', 'a', 'b',   'c', 'c', 'b';

  Stat::Factor<CArrayVector<char, 13> > f1d(datai.col(1));
  f1d.run();
  stk_cout << _T("nbLevels= ")  << f1d.nbLevels() << _T("\n");
  stk_cout << _T("asInteger= ") << f1d.asInteger().transpose() << _T("\n");
  stk_cout <<_T("Levels: ")     << f1d.levels().transpose();
  stk_cout <<_T("Levels counts: ") << f1d.counts().transpose();

  Stat::MultiFactor<CArray<char, 13, 3> > f2d(datai);
  f2d.run();
  stk_cout << _T("nbLevels= ")   << f2d.nbLevels() << _T("\n");
  stk_cout << _T("asInteger=\n") << f2d.asInteger() << _T("\n");
  for (int i=f2d.levels().begin(); i<f2d.levels().end(); ++i)
  { stk_cout <<_T("Levels variable ")<< i << _T(": ")<< f2d.levels()[i].transpose();}
  stk_cout << _T("\n");
  for (int i=f2d.levels().begin(); i<f2d.levels().end(); ++i)
  { stk_cout <<_T("Levels counts ")<< i << _T(": ")<< f2d.counts()[i].transpose();}
  stk_cout << _T("\n");

 return 0;
}
