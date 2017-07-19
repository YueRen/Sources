#include <Singular/libsingular.h>
#include <gfanlib/gfanlib.h>

#include <kernel/GBEngine/kstd1.h>
#include <Singular/dyn_modules/gfanlib/callgfanlib_conversion.h>
#include <Singular/dyn_modules/gfanlib/bbcone.h>
#include <Singular/dyn_modules/gfanlib/initial.h>

#include <set>
#include <groebnerCone.h>
#include <singularWrappers.h>
#include <symmetry.h>

#define TRAVERSAL_TIMINGS_ON 1

#include <iostream>

#if TRAVERSAL_TIMINGS_ON
#include <iostream>
#include <ctime>
#endif


ring createNeighbouringRing(ring homeRing, const gfan::ZVector &interiorFacetPoint, const gfan::ZVector &outerFacetNormal, const gfan::ZVector &w)
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
  neighbouringRing->wvhdl[0] = ZVectorToIntStar(w,ok);
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


std::set<tropical::groebnerCone> tropicalTraversal(const tropical::groebnerCone &startingCone, const gfan::ZVector &w,
                                                   const std::set<std::vector<int> > &symmetryGroup)
{
  std::set<tropical::groebnerCone> workingList;
  workingList.insert(startingCone);
  std::set<tropical::groebnerCone> finishedList;
  std::set<gfan::ZVector> knownTropicalLinks;

  while (!workingList.empty())
  {
    tropical::groebnerCone maximalGroebnerCone = *(workingList.begin());

    gfan::ZCone maximalPolyhedralCone = maximalGroebnerCone.getPolyhedralCone();
    gfan::ZMatrix interiorFacetPoints = interiorPointsOfFacets(maximalPolyhedralCone);

    ideal polynomialIdeal = maximalGroebnerCone.getPolynomialIdeal();
    ring polynomialRing = maximalGroebnerCone.getPolynomialRing();

#if TRAVERSAL_TIMINGS_ON
    std::clock_t tropicalLinkTime = 0;
    std::clock_t initialStdTime = 0;
    std::clock_t initialRedTime = 0;
    std::clock_t entireIterationTime = 0;
    std::clock_t titerationstart = std::clock();
    std::clock_t ms = 0;
#endif

    for (int i=0; i<interiorFacetPoints.getHeight(); i++)
    {
      gfan::ZVector interiorFacetPoint = interiorFacetPoints[i].toVector();
      std::pair<std::set<gfan::ZVector>::iterator,bool> todo = knownTropicalLinks.insert(minimalRepresentative(interiorFacetPoint,symmetryGroup));
      if (todo.second)
      {
        ideal initialIdeal = initial(polynomialIdeal,polynomialRing,interiorFacetPoint);

#if TRAVERSAL_TIMINGS_ON
        std::clock_t tlinkstart = std::clock();
#endif

        gfan::ZMatrix tropicalLink = callTropicalLinkNewton(initialIdeal,polynomialRing);

#if TRAVERSAL_TIMINGS_ON
        std::clock_t tlinkend = std::clock();
        ms = 1000.0 * (tlinkend - tlinkstart)/ CLOCKS_PER_SEC;
        tropicalLinkTime += ms;
#endif

        for (int j=0; j<tropicalLink.getHeight(); j++)
        {
#if TRAVERSAL_TIMINGS_ON
          std::clock_t tloopstart = std::clock();
#endif
          /* create a ring with weighted ordering  */
          ring s = createNeighbouringRing(polynomialRing,interiorFacetPoint,tropicalLink[j].toVector(),w);
          nMapFunc identity1 = n_SetMap(polynomialRing->cf,s->cf);

          int k = IDELEMS(initialIdeal);
          ideal inIs = idInit(k);
          for (int l=0; l<k; l++)
            inIs->m[l] = p_PermPoly(initialIdeal->m[l],NULL,polynomialRing,s,identity1,NULL,0);
#if TRAVERSAL_TIMINGS_ON
          std::clock_t tstdstart = std::clock();
#endif
          ideal inIsGB = tropical_kStd_wrapper(inIs,s);
#if TRAVERSAL_TIMINGS_ON
          std::clock_t tstdend = std::clock();
          ms = 1000.0 * (tstdend - tstdstart)/ CLOCKS_PER_SEC;
          initialStdTime += ms;
#endif
          ideal inIsGBNF = tropical_kNF_wrapper(inIsGB,s,polynomialIdeal,polynomialRing); // todo: time this + for loop below
          id_Delete(&inIs,s);

          k = IDELEMS(inIsGB);
          ideal IsGBnonred = idInit(k);
          for (int l=0; l<k; l++)
          {
            IsGBnonred->m[l] = p_Add_q(inIsGB->m[l],p_Neg(inIsGBNF->m[l],s),s);
            inIsGB->m[l] = NULL;
            inIsGBNF->m[l] = NULL;
          }
#if TRAVERSAL_TIMINGS_ON
          std::clock_t tredstart = std::clock();
#endif
          ideal IsGB = tropical_kInterRed_wrapper(IsGBnonred,s);
#if TRAVERSAL_TIMINGS_ON
          std::clock_t tredend = std::clock();
          ms = 1000.0 * (tredend - tredstart)/ CLOCKS_PER_SEC;
          initialRedTime += ms;
#endif
          id_Delete(&IsGBnonred,s);
          id_Delete(&inIsGB,s);
          id_Delete(&inIsGBNF,s);

          gfan::ZMatrix outerFacetNormal(0,tropicalLink.getWidth());
          outerFacetNormal.appendRow(tropicalLink[j].toVector());
          inIsGB = initial(IsGB,s,interiorFacetPoint,outerFacetNormal);

          tropical::groebnerCone neighbouringGroebnerCone(IsGB,inIsGB,s,symmetryGroup); // todo: time this
          id_Delete(&IsGB,s);
          id_Delete(&inIsGB,s);
          rDelete(s);

          if (finishedList.count(neighbouringGroebnerCone)==0) // todo: time this
            workingList.insert(neighbouringGroebnerCone);
#if TRAVERSAL_TIMINGS_ON
          std::clock_t tloopend = std::clock();
          ms = 1000.0 * (tloopend - tloopstart)/ CLOCKS_PER_SEC;
          std::cout << "one iteration: " << ms << " milliseconds." << std::endl;
#endif
        }
        id_Delete(&initialIdeal,polynomialRing);
      }
    }

    maximalGroebnerCone.deletePolynomialIdealAndRing();
    maximalGroebnerCone.deletePolyhedralCone();
    finishedList.insert(maximalGroebnerCone);
    workingList.erase(maximalGroebnerCone);

#if TRAVERSAL_TIMINGS_ON
    std::clock_t titerationend = std::clock();
    entireIterationTime = 1000.0 *(titerationend - titerationstart)/ CLOCKS_PER_SEC;

    std::cout << "tropicalLinkTime: " << tropicalLinkTime << " ms" << std::endl;
    std::cout << "initialStdTime: " << initialStdTime << " ms" <<  std::endl;
    std::cout << "initialRedTime: " << initialRedTime << " ms" <<  std::endl;
    std::cout << "entireIterationTime: " << entireIterationTime << " ms" << std::endl;
#endif

    if (printlevel > 0)
      Print("cones finished: %lu   cones in working list: %lu\n",
      (unsigned long)finishedList.size(), (unsigned long)workingList.size());
  }

  return finishedList;
}
