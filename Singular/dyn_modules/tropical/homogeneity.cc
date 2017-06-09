#include <gfanlib/gfanlib.h>
#include <Singular/libsingular.h>

#include <Singular/dyn_modules/gfanlib/callgfanlib_conversion.h>
#include <Singular/dyn_modules/gfanlib/initial.h>

#include <singularWrappers.h>

bool isWeightedHomogeneous(const ideal I, const ring r, const gfan::ZVector w)
{
  for (int i=0; i<IDELEMS(I); i++)
  {
    poly g = I->m[i];
    long d = wDeg(g,r,w);
    pIter(g);
    for (; g!=NULL; pIter(g))
    {
      if (wDeg(g,r,w)!=d)
        return false;
    }
  }

  return true;
}


static gfan::ZCone homogeneitySpace(ideal I, ring r)
{
  int n = rVar(r);
  poly g;
  int* leadexpv = (int*) omAlloc((n+1)*sizeof(int));
  int* tailexpv = (int*) omAlloc((n+1)*sizeof(int));
  gfan::ZVector leadexpw = gfan::ZVector(n);
  gfan::ZVector tailexpw = gfan::ZVector(n);
  gfan::ZMatrix equations = gfan::ZMatrix(0,n);
  for (int i=0; i<IDELEMS(I); i++)
  {
    g = (poly) I->m[i];
    if (g)
    {
      p_GetExpV(g,leadexpv,r);
      leadexpw = intStar2ZVector(n,leadexpv);
      pIter(g);
      while (g)
      {
        p_GetExpV(g,tailexpv,r);
        tailexpw = intStar2ZVector(n,tailexpv);
        equations.appendRow(leadexpw-tailexpw);
        pIter(g);
      }
    }
  }
  omFreeSize(leadexpv,(n+1)*sizeof(int));
  omFreeSize(tailexpv,(n+1)*sizeof(int));
  return gfan::ZCone(gfan::ZMatrix(0, equations.getWidth()),equations);
}


gfan::ZVector homogeneityVector(const ideal I, const ring r)
{
  // 1. check if the current ordering is a weighted ordering.
  //    if so, check if the ideal is homogeneous with respect to the weight vector
  if (r->order[0]==ringorder_a || r->order[0]==ringorder_wp)
  {
    if (r->block0[0]==0 && r->block1[0]==rVar(r))
    {
      gfan::ZVector w = wvhdlEntryToZVector(rVar(r),r->wvhdl[0]);
      if (w.isPositive() && isWeightedHomogeneous(I,r,w))
        return w;
    }
  }

  // 2. check if the ideal is homogeneous with respect to the (1,...,1) vector
  gfan::ZVector w = gfan::ZVector::allOnes(rVar(r));
  if (isWeightedHomogeneous(I,r,w))
    return w;

  // 3. compute the homogeneity space of the ideal
  //    and check if it contains a strictly positive vector
  ideal stdI = tropical_kStd_wrapper(I,r);
  gfan::ZCone C0 = homogeneitySpace(stdI,r);
  gfan::ZCone C0pos = gfan::intersection(C0,gfan::ZCone::positiveOrthant(rVar(r)));
  w = C0pos.getRelativeInteriorPoint();

  id_Delete(&stdI,r);
  if (w.isPositive())
    return w;

  return gfan::ZVector();
}
