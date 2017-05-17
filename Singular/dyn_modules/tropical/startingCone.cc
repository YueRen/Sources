#include <gfanlib/gfanlib.h>
#include <Singular/libsingular.h>

#include <kernel/combinatorics/stairc.h> // for scDimInt
#include <Singular/dyn_modules/gfanlib/callgfanlib_conversion.h>
#include <Singular/dyn_modules/gfanlib/initial.h>

#include <groebnerCone.h>
#include <singularWrappers.h>
#include <symmetry.h>

ring createRingWeightedOrdering(const ring r, gfan::ZVector w)
{
  bool ok;

  ring s = rCopy0(r,FALSE,FALSE);
  int n = rVar(s);
  s->order = (rRingOrder_t*) omAlloc0(5*sizeof(rRingOrder_t));
  s->block0 = (int*) omAlloc0(5*sizeof(int));
  s->block1 = (int*) omAlloc0(5*sizeof(int));
  s->wvhdl = (int**) omAlloc0(5*sizeof(int**));


  s->order[0] = ringorder_a;
  s->block0[0] = 1;
  s->block1[0] = n;
  s->wvhdl[0] = ZVectorToIntStar(gfan::ZVector::allOnes(n),ok);
  s->order[1] = ringorder_a;
  s->block0[1] = 1;
  s->block1[1] = n;
  s->wvhdl[1] = ZVectorToIntStar(w,ok);
  s->order[2] = ringorder_lp;
  s->block0[2] = 1;
  s->block1[2] = n;
  s->order[3] = ringorder_C;
  rComplete(s);
  rTest(s);

  return s;
}


int dim(ideal I, ring r)
{
  ring origin = currRing;
  if (origin != r)
    rChangeCurrRing(r);
  int d = scDimInt(I,currRing->qideal);
  if (origin != r)
    rChangeCurrRing(origin);
  return d;
}


tropical::groebnerCone tropicalStartingCone(const ideal I, const ring r, const std::set<std::vector<int> > &symmetryGroup)
{
  bool done = false;
  tropical::groebnerCone startingCone;
  int d = -1;

  do
  {
    gfan::QVector tropicalStartingPoint0 = callTropicalPointNewton(I,r);
    // WARNING: homogeneous case only
    gfan::ZVector tropicalStartingPoint = QToZVectorPrimitive(tropicalStartingPoint0);

    ring polynomialRing = createRingWeightedOrdering(r,tropicalStartingPoint);
    nMapFunc identity = n_SetMap(r->cf,polynomialRing->cf);
    int k = IDELEMS(I);
    ideal Is = idInit(k);
    for (int i=0; i<k; i++)
      Is->m[i] = p_PermPoly(I->m[i],NULL,r,polynomialRing,identity,NULL,0);
    ideal polynomialIdeal = tropical_kStd_wrapper(Is,polynomialRing);
    if (d<0)
    {
      d = dim(polynomialIdeal,polynomialRing);
    }
    id_Delete(&Is,polynomialRing);

    startingCone = tropical::groebnerCone(polynomialIdeal,polynomialRing,tropicalStartingPoint,symmetryGroup);
    id_Delete(&polynomialIdeal,polynomialRing);
    rDelete(polynomialRing);
  }
  while (startingCone.getPolyhedralCone().dimension() < d);

  return startingCone;
}
