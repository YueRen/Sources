#include <Singular/libsingular.h>
#include <gfanlib/gfanlib.h>
#include <Singular/dyn_modules/gfanlib/callgfanlib_conversion.h>
#include <Singular/dyn_modules/gfanlib/initial.h>
#include <Singular/dyn_modules/gfanlib/bbcone.h>

#include <groebnerCone.h>
#include <symmetry.h>

#include <set>
#include <vector>

namespace tropical {

  groebnerCone::groebnerCone():
    polynomialRing(NULL),
    polynomialIdeal(NULL),
    polyhedralCone(gfan::ZCone(0)),
    uniquePoint(gfan::ZVector(0))
  {
  }

  groebnerCone::groebnerCone(const groebnerCone &sigma):
    polynomialRing(rCopy(sigma.getPolynomialRing())),
    polynomialIdeal(id_Copy(sigma.getPolynomialIdeal(),sigma.getPolynomialRing())),
    polyhedralCone(gfan::ZCone(sigma.getPolyhedralCone())),
    uniquePoint(gfan::ZVector(sigma.getUniquePoint()))
  {
  }


  groebnerCone::groebnerCone(ideal I, ring r, gfan::ZVector interiorPoint, std::set<std::vector<int> > symmetryGroup):
    polynomialIdeal(id_Copy(I,r)),
    polynomialRing(rCopy(r))
  {
    int n = rVar(polynomialRing);
    gfan::ZMatrix inequalities = gfan::ZMatrix(0,n);
    gfan::ZMatrix equations = gfan::ZMatrix(0,n);
    int* expv = (int*) omAlloc((n+1)*sizeof(int));
    for (int i=0; i<IDELEMS(polynomialIdeal); i++)
    {
      poly g = polynomialIdeal->m[i];
      if (g)
      {
        p_GetExpV(g,expv,polynomialRing);
        gfan::ZVector leadexpv = intStar2ZVector(n,expv);
        long d = wDeg(g,polynomialRing,interiorPoint);
        for (pIter(g); g; pIter(g))
        {
          p_GetExpV(g,expv,polynomialRing);
          gfan::ZVector tailexpv = intStar2ZVector(n,expv);
          if (wDeg(g,polynomialRing,interiorPoint)==d)
            equations.appendRow(leadexpv-tailexpv);
          else
          {
            assume(wDeg(g,polynomialRing,interiorPoint)<d);
            inequalities.appendRow(leadexpv-tailexpv);
          }
        }
      }
    }
    omFreeSize(expv,(n+1)*sizeof(int));

    polyhedralCone = gfan::ZCone(inequalities,equations);
    polyhedralCone.canonicalize();
    uniquePoint = minimalRepresentative(polyhedralCone.getUniquePoint(),symmetryGroup);
  }

  groebnerCone::groebnerCone(ideal I, ideal inI, ring r, std::set<std::vector<int> > symmetryGroup):
    polynomialIdeal(id_Copy(I,r)),
    polynomialRing(rCopy(r))
  {
    assume(IDELEMS(I)==IDELEMS(inI));

    int n = rVar(polynomialRing);
    int* expv = (int*) omAlloc((n+1)*sizeof(int));
    gfan::ZMatrix inequalities = gfan::ZMatrix(0,n);
    gfan::ZMatrix equations = gfan::ZMatrix(0,n);
    for (int i=0; i<IDELEMS(I); i++)
    {
      poly g = I->m[i];
      if (g)
      {
        p_GetExpV(g,expv,polynomialRing);
        gfan::ZVector leadexpv = intStar2ZVector(n,expv);
        for (pIter(g); g; pIter(g))
        {
          p_GetExpV(g,expv,polynomialRing);
          gfan::ZVector tailexpv = intStar2ZVector(n,expv);
          inequalities.appendRow(leadexpv-tailexpv);
        }
      }
      g = inI->m[i];
      if (g)
      {
        p_GetExpV(g,expv,polynomialRing);
        gfan::ZVector leadexpv = intStar2ZVector(n,expv);
        for (pIter(g); g; pIter(g))
        {
          p_GetExpV(g,expv,polynomialRing);
          gfan::ZVector tailexpv = intStar2ZVector(n,expv);
          equations.appendRow(leadexpv-tailexpv);
        }
      }
    }
    omFreeSize(expv,(n+1)*sizeof(int));

    polyhedralCone = gfan::ZCone(inequalities,equations);
    polyhedralCone.canonicalize();
    uniquePoint = minimalRepresentative(polyhedralCone.getUniquePoint(),symmetryGroup);
  }


  groebnerCone::~groebnerCone()
  {
    id_Delete(&polynomialIdeal,polynomialRing);
    rDelete(polynomialRing);
  }


  groebnerCone& groebnerCone::operator=(const groebnerCone& sigma)
  {
    polynomialRing = rCopy(sigma.getPolynomialRing());
    polynomialIdeal = id_Copy(sigma.getPolynomialIdeal(),sigma.getPolynomialRing());
    polyhedralCone = sigma.getPolyhedralCone();
    uniquePoint = sigma.getUniquePoint();
    return *this;
  }


  bool groebnerCone::operator<(const groebnerCone& sigma) const
  {
    return uniquePoint<sigma.getUniquePoint();
  }

  gfan::ZFan* groebnerConesToZFanStar(std::set<groebnerCone>& groebnerCones)
  {
    std::set<groebnerCone>::iterator groebnerCone = groebnerCones.begin();
    gfan::ZCone polyhedralCone = groebnerCone->getPolyhedralCone();
    gfan::ZFan* groebnerSubfan = new gfan::ZFan(polyhedralCone.ambientDimension());
    groebnerSubfan->insert(polyhedralCone);

    for (; groebnerCone!=groebnerCones.end(); ++groebnerCone)
      groebnerSubfan->insert(groebnerCone->getPolyhedralCone());

    return groebnerSubfan;
  }


  lists groebnerConesToListOfZCones(std::set<groebnerCone>& groebnerCones)
  {
    lists L = (lists)omAllocBin(slists_bin);
    L->Init(groebnerCones.size());

    int i=0;
    for (std::set<groebnerCone>::iterator groebnerCone = groebnerCones.begin(); groebnerCone!=groebnerCones.end(); ++groebnerCone)
    {
      L->m[i].rtyp = coneID;
      L->m[i].data = (void*) new gfan::ZCone(groebnerCone->getPolyhedralCone());
      i++;
    }

    return L;
  }
}
