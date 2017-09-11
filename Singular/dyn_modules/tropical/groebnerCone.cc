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
    polynomialRing(NULL),
    polynomialIdeal(NULL),
    polyhedralCone(gfan::ZCone(sigma.getPolyhedralCone())),
    uniquePoint(gfan::ZVector(sigma.getUniquePoint()))
  {
    if (sigma.getPolynomialRing()!=NULL)
      polynomialRing = rCopy(sigma.getPolynomialRing());
    if (sigma.getPolynomialIdeal()!=NULL)
      polynomialIdeal = id_Copy(sigma.getPolynomialIdeal(),sigma.getPolynomialRing());
  }

	groebnerCone::groebnerCone(ideal I, ring r):
		polynomialIdeal(id_Copy(I,r)),
    polynomialRing(rCopy(r))
 	{
		int n = rVar(polynomialRing);
    gfan::ZMatrix inequalities = gfan::ZMatrix(0,n);
    int* expv = (int*) omAlloc((n+1)*sizeof(int));
		std::cout << "constructing inequalities of Groebner cone";
    for (int i=0; i<IDELEMS(polynomialIdeal); i++)
    {
      poly g = polynomialIdeal->m[i];
      if (g)
      {
        p_GetExpV(g,expv,polynomialRing);
        gfan::ZVector leadexpv = intStar2ZVector(n,expv);

				// 1. Compute the vertices of the Newton polytope
				gfan::ZVector leadexpvLifted = leadexpv;
				leadexpvLifted.push_back((long)1);
				gfan::ZMatrix ExpvLifted(0,n+1);
				ExpvLifted.appendRow(leadexpvLifted);
        for (pIter(g); g; pIter(g))
        {
          p_GetExpV(g,expv,polynomialRing);
          gfan::ZVector tailexpvLifted = intStar2ZVector(n,expv);
					tailexpvLifted.push_back((long)1);
					ExpvLifted.appendRow(tailexpvLifted);
        }
				gfan::ZCone newtonPolytope = gfan::ZCone(ExpvLifted,gfan::ZMatrix(0,n+1));
				gfan::ZMatrix ExpvExtremal = newtonPolytope.getFacets();
				ExpvExtremal = ExpvExtremal.submatrix(0,0,ExpvExtremal.getHeight(),n);

				// 2. Read off the inequalities of the Groebner cone
				for (int i=0; i<ExpvExtremal.getHeight(); i++)
				{
					gfan::ZVector tailexpv = ExpvExtremal[i].toVector();
					inequalities.appendRow(leadexpv-tailexpv);		
				}
				std::cout << "." << std::flush;
			}
    }
    omFreeSize(expv,(n+1)*sizeof(int));

		std::cout << std::endl << "constructing facets of Groebner cone" << std::endl;
    polyhedralCone = gfan::ZCone(inequalities,gfan::ZMatrix(0,inequalities.getWidth()));
    polyhedralCone.canonicalize();
    uniquePoint = polyhedralCone.getRelativeInteriorPoint();
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
		if (symmetryGroup.empty())
			uniquePoint = polyhedralCone.getRelativeInteriorPoint();
		else
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
    if (polynomialIdeal!=NULL)
      id_Delete(&polynomialIdeal,polynomialRing);
    if (polynomialRing!=NULL)
    {
      rDelete(polynomialRing);
      polynomialRing = NULL;
    }
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


  void groebnerCone::deletePolynomialIdealAndRing()
  {
    if (polynomialIdeal!=NULL)
      id_Delete(&polynomialIdeal,polynomialRing);
    if (polynomialRing!=NULL)
    {
      rDelete(polynomialRing);
      polynomialRing = NULL;
    }
  }

  void groebnerCone::deletePolyhedralCone()
  {
    polyhedralCone = gfan::ZCone();
  }

	gfan::ZMatrix matrixMultiplication(const gfan::ZMatrix &M, const gfan::ZMatrix &N)
	{
		assert(M.getWidth()==N.getHeight());
		gfan::ZMatrix MN(M.getHeight(),N.getWidth());
		for (int i=0; i<M.getHeight(); i++)
		{
			for (int j=0; j<N.getWidth(); j++)
			{
				gfan::Integer MNij((long)0);
				for (int k=0; k<M.getWidth(); k++)
					MNij += M[i][k] * N[k][j];
				MN[i][j] = MNij;
			}
		}
		return MN;
	}

	std::pair<gfan::ZVector,gfan::ZVector> groebnerCone::facetPointingTo(const gfan::ZMatrix &P) const
	{
		gfan::ZMatrix outerFacetNormals = -(this->polyhedralCone.getFacets());

		gfan::ZMatrix M = matrixMultiplication(outerFacetNormals,P.transposed());
		gfan::ZVector maximalDotProduct(M.getWidth());
		int maximalIndex = -1;
		for (int i=0; i<outerFacetNormals.getHeight(); i++)
		{
			gfan::ZVector dotProduct = M[i].toVector();
			if (maximalDotProduct<dotProduct)
			{
				maximalDotProduct = dotProduct;
				maximalIndex = i;
			}
		}

		if (maximalIndex<0)
		{
			return std::make_pair(gfan::ZVector(),gfan::ZVector());
		}
		
		gfan::ZVector outerFacetNormal = outerFacetNormals[maximalIndex].toVector();
		int n = polyhedralCone.ambientDimension();
		gfan::ZMatrix facetEquation(0,n);
		facetEquation.appendRow(outerFacetNormal);
		gfan::ZMatrix identity(n,n);
		for (int i=0; i<n; i++)
			identity[i][i] = 1;
		gfan::ZCone supportingPositiveHyperplane = gfan::ZCone(identity,facetEquation);
		gfan::ZCone facet = gfan::intersection(polyhedralCone,supportingPositiveHyperplane);
		gfan::ZVector interiorFacetPoint = facet.getRelativeInteriorPoint();
		
		return std::make_pair(interiorFacetPoint,outerFacetNormal);
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
