#include <Singular/libsingular.h>
#include <gfanlib/gfanlib.h>

#include <kernel/GBEngine/kstd1.h>
#include <Singular/dyn_modules/gfanlib/callgfanlib_conversion.h>
#include <Singular/dyn_modules/gfanlib/bbcone.h>
#include <Singular/dyn_modules/gfanlib/initial.h>

#include <set>
#include <iostream>

#include <groebnerCone.h>
#include <singularWrappers.h>
#include <symmetry.h>


ring createNeighbouringRing(ring homeRing, gfan::ZVector interiorFacetPoint, gfan::ZVector outerFacetNormal)
{
  bool ok;
  ring neighbouringRing = rCopy0(homeRing,TRUE,FALSE);
  int n = rVar(neighbouringRing);
  neighbouringRing->order = (rRingOrder_t*) omAlloc0(6*sizeof(rRingOrder_t));
  neighbouringRing->block0 = (int*) omAlloc0(6*sizeof(int));
  neighbouringRing->block1 = (int*) omAlloc0(6*sizeof(int));
  neighbouringRing->wvhdl = (int**) omAlloc0(6*sizeof(int**));
  neighbouringRing->order[0] = ringorder_a;
  neighbouringRing->block0[0] = 1;
  neighbouringRing->block1[0] = n;
  neighbouringRing->wvhdl[0] = ZVectorToIntStar(gfan::ZVector::allOnes(n),ok);
  neighbouringRing->order[1] = ringorder_a;
  neighbouringRing->block0[1] = 1;
  neighbouringRing->block1[1] = n;
  neighbouringRing->wvhdl[1] = ZVectorToIntStar(interiorFacetPoint,ok);
  neighbouringRing->order[2] = ringorder_a;
  neighbouringRing->block0[2] = 1;
  neighbouringRing->block1[2] = n;
  neighbouringRing->wvhdl[2] = ZVectorToIntStar(outerFacetNormal,ok);
  neighbouringRing->order[3] = ringorder_lp;
  neighbouringRing->block0[3] = 1;
  neighbouringRing->block1[3] = n;
  neighbouringRing->order[4] = ringorder_C;
  rComplete(neighbouringRing);
  rTest(neighbouringRing);
  return neighbouringRing;
}


std::set<tropical::groebnerCone> tropicalTraversal(const tropical::groebnerCone &startingCone, const std::set<std::vector<int> > &symmetryGroup)
{
  std::set<tropical::groebnerCone> workingList;
  workingList.insert(startingCone);
  std::set<tropical::groebnerCone> finishedList;
  std::set<gfan::ZVector> knownTropicalLinks;

  while (!workingList.empty())
  {
    tropical::groebnerCone maximalGroebnerCone = *(workingList.begin());
    finishedList.insert(maximalGroebnerCone);
    workingList.erase(maximalGroebnerCone);

    gfan::ZCone maximalPolyhedralCone = maximalGroebnerCone.getPolyhedralCone();
    gfan::ZMatrix interiorFacetPoints = interiorPointsOfFacets(maximalPolyhedralCone);

    ideal polynomialIdeal = maximalGroebnerCone.getPolynomialIdeal();
    ring polynomialRing = maximalGroebnerCone.getPolynomialRing();

    for (int i=0; i<interiorFacetPoints.getHeight(); i++)
    {
      gfan::ZVector interiorFacetPoint = interiorFacetPoints[i].toVector();
      std::pair<std::set<gfan::ZVector>::iterator,bool> todo = knownTropicalLinks.insert(minimalRepresentative(interiorFacetPoint,symmetryGroup));
      if (todo.second)
      {
        ideal initialIdeal = initial(polynomialIdeal,polynomialRing,interiorFacetPoint);
        gfan::ZMatrix tropicalLink = callTropicalLinkNewton(initialIdeal,polynomialRing);
        for (int j=0; j<tropicalLink.getHeight(); j++)
        {
          /* create a ring with weighted ordering  */
          ring s = createNeighbouringRing(polynomialRing,interiorFacetPoint,tropicalLink[j].toVector());
          nMapFunc identity1 = n_SetMap(polynomialRing->cf,s->cf);

          int k = IDELEMS(initialIdeal);
          ideal inIs = idInit(k);
          for (int l=0; l<k; l++)
            inIs->m[l] = p_PermPoly(initialIdeal->m[l],NULL,polynomialRing,s,identity1,NULL,0);
          ideal inIsGB = tropical_kStd_wrapper(inIs,s);
          ideal inIsGBNF = tropical_kNF_wrapper(inIsGB,s,polynomialIdeal,polynomialRing);
          id_Delete(&inIs,s);

          k = IDELEMS(inIsGB);
          ideal IsGBnonred = idInit(k);
          for (int l=0; l<k; l++)
          {
            IsGBnonred->m[l] = p_Add_q(inIsGB->m[l],p_Neg(inIsGBNF->m[l],s),s);
            inIsGB->m[l] = NULL;
            inIsGBNF->m[l] = NULL;
          }
          ideal IsGB = tropical_kStd_wrapper(IsGBnonred,s);
          id_Delete(&IsGBnonred,s);
          id_Delete(&inIsGB,s);
          id_Delete(&inIsGBNF,s);

          gfan::ZMatrix outerFacetNormal(0,tropicalLink[j].toVector().size());
          outerFacetNormal.appendRow(tropicalLink[j].toVector());
          inIsGB = initial(IsGB,s,interiorFacetPoint,outerFacetNormal);


          tropical::groebnerCone neighbouringGroebnerCone(IsGB,inIsGB,s,symmetryGroup);
          id_Delete(&IsGB,s);
          id_Delete(&inIsGB,s);

          if (finishedList.count(neighbouringGroebnerCone)==0)
            workingList.insert(neighbouringGroebnerCone);
        }
      }
    }

    std::cout << "workingList.size(): " << workingList.size()
              << "   finishedList.size(): " << finishedList.size() << std::endl;
  }

  return finishedList;
}
