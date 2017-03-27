#include <vector>
#include <set>

#include <Singular/libsingular.h>
#include <coeffs/bigintmat.h>
#include <coeffs/longrat.h>
#include <coeffs/numbers.h>
#include <gfanlib/gfanlib.h>

std::set<std::vector<int> > convertPermutations(lists generators0)
{
  std::set<std::vector<int> > generators1;
  for (int i=0; i<=lSize(generators0); i++)
  {
    if (generators0->m[i].Typ() != INTVEC_CMD)
    {
      WerrorS("generatePermutationGroup: illegal entries in list");
      return std::set<std::vector<int> >();
    }
    intvec* generator0 = (intvec*) generators0->m[i].Data();

    std::vector<int> generator1(generator0->rows());
    for (int j=0; j<generator0->rows(); j++)
      generator1[j] = (*generator0)[j]-1;

    generators1.insert(generator1);
  }
  return generators1;
}


lists convertPermutations(std::set<std::vector<int> > permutationGroup1)
{
  lists permutationGroup0 = (lists) omAllocBin(slists_bin);
  permutationGroup0->Init(permutationGroup1.size()); int i=0;
  for (std::set<std::vector<int> >::iterator it = permutationGroup1.begin(); it!=permutationGroup1.end(); ++it, ++i)
  {
    intvec* permutation0 = new intvec(it->size());
    for (std::size_t j=0; j<it->size(); j++)
      (*permutation0)[j] = (*it)[j]+1;

    permutationGroup0->m[i].rtyp = INTVEC_CMD;
    permutationGroup0->m[i].data = (void*) permutation0;
  }
  return permutationGroup0;
}


int rationalToInt(const gfan::Rational &r)
{
  mpq_t i0;
  mpq_init(i0);
  r.setGmp(i0);

  assume(mpz_cmp_si(mpq_denref(i0),1)==0);

  int i = mpz_get_si(mpq_numref(i0));

  mpq_canonicalize(i0);
  mpq_clear(i0);

  return i;
}


gfan::Rational numberToRational(number c, coeffs cf)
{
  number n=n_GetNumerator(c,cf);
  number d=n_GetDenom(c,cf);

  gfan::Rational a = n_Int(n,cf);
  gfan::Rational b = n_Int(d,cf);
  n_Delete(&n,cf);
  n_Delete(&d,cf);

  return (a/b);
}


gfan::ZVector intvecStarToZVector(intvec* v)
{
  gfan::ZVector w(v->length());
  for (int i=0; i<v->length(); i++)
  {
    w[i] = gfan::Integer((*v)[i]);
  }
  return w;
}


gfan::QMatrix matrixToQMatrix(matrix M, ring r)
{
  gfan::QMatrix N(M->rows(),M->cols());
  for (int i=0; i<M->rows(); i++)
  {
    for (int j=0; j<M->cols(); j++)
    {
      poly cpoly = MATELEM(M,i+1,j+1);
      if (cpoly==NULL)
      {
        N[i][j] = gfan::Rational((long) 0);
      }
      else
      {
        number cnumber = p_GetCoeff(cpoly,r);
        N[i][j] = numberToRational(cnumber,r->cf);
      }
    }
  }
  return N;
}
