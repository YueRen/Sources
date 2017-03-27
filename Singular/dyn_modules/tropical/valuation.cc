#include <Singular/libsingular.h>
#include <gfanlib/gfanlib.h>

gfan::Rational puiseuxValuation(number c, const coeffs cf)
{
  assume(getCoeffType(cf) == n_transExt);

  number n=n_GetNumerator(c,cf);
  number d=n_GetDenom(c,cf);
  gfan::Rational v(n_ParDeg(n,cf)-n_ParDeg(d,cf));

  n_Delete(&n,cf);
  n_Delete(&d,cf);
  return v;
}
