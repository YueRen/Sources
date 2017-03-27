#include <kernel/mod2.h>

#include <Singular/libsingular.h>
#include <gfanlib/gfanlib.h>

#include <Singular/dyn_modules/gfanlib/callgfanlib_conversion.h>
#include <Singular/dyn_modules/gfanlib/bbcone.h>
#include <Singular/dyn_modules/gfanlib/bbfan.h>

#include <conversion.h>
#include <singularWrappers.h>
#include <symmetry.h>
#include <valuation.h>
#include <groebnerCone.h>
#include <startingCone.h>
#include <traversal.h>


BOOLEAN permutationGroup(leftv res, leftv args)
{
  leftv u = args;
  if ((u != NULL) && (u->Typ() == LIST_CMD) && (u->next == NULL))
  {
    lists generators0 = (lists) u->Data();

    std::set<std::vector<int> > generators1 = convertPermutations(generators0);
    std::set<std::vector<int> > permutationGroup1 = permutationGroup(generators1);
    lists permutationGroup0 = convertPermutations(permutationGroup1);

    res->rtyp = LIST_CMD;
    res->data = (void*) permutationGroup0;
    return FALSE;
  }
  WerrorS("permutationGroup: unexpected parameters");
  return TRUE;
}


BOOLEAN minimalRepresentative(leftv res, leftv args)
{
  leftv u = args;
  if ((u != NULL) && (u->Typ() == BIGINTMAT_CMD))
  {
    leftv v = u->next;
    if ((v != NULL) && (v->Typ() == LIST_CMD))
    {
      bigintmat* w0 = (bigintmat*) u->Data();
      lists generators0 = (lists) v->Data();

      gfan::ZVector* w = bigintmatToZVector(w0);

      std::set<std::vector<int> > generators1 = convertPermutations(generators0);
      gfan::ZVector wmin = minimalRepresentative(*w,generators1);
      delete w;

      res->rtyp = BIGINTMAT_CMD;
      res->data = (void*) zVectorToBigintmat(wmin);
      return FALSE;
    }
  }
  WerrorS("minimalRepresentative: unexpected parameters");
  return TRUE;
}


BOOLEAN valuation(leftv res, leftv args)
{
  leftv u = args;
  if ((u != NULL) && (u->Typ() == NUMBER_CMD))
  {
    number c = (number) u->Data();
    gfan::Rational v = puiseuxValuation(c,currRing->cf);

    res->rtyp = INT_CMD;
    res->data = (void*) rationalToInt(v);
    return FALSE;
  }
  WerrorS("valuation: unexpected parameters");
  return TRUE;
}


BOOLEAN tropicalStartingConeNewton(leftv res, leftv args)
{
  leftv u = args;
  if ((u != NULL) && (u->Typ() == IDEAL_CMD))
  {
    ideal I = (ideal) u->Data();
    tropical::groebnerCone sigma = tropicalStartingCone(I,currRing);

    res->rtyp = coneID;
    res->data = (void*) new gfan::ZCone(sigma.getPolyhedralCone());
    return FALSE;
  }
  WerrorS("tropicalStartingCone: unexpected parameters");
  return TRUE;
}


BOOLEAN tropicalDebug(leftv res, leftv args)
{
  leftv u = args;
  ideal I = (ideal) u->Data();

  std::cerr << "computing starting cone..." << std::endl;
  tropical::groebnerCone sigma = tropicalStartingCone(I,currRing,std::set<std::vector<int> >());
  std::cerr << "starting traversal..." << std::endl;
  std::set<tropical::groebnerCone> TropI = tropicalTraversal(sigma);

  res->rtyp = fanID;
  res->data = (void*) groebnerConesToZFanStar(TropI);
  return FALSE;
}

extern "C" int SI_MOD_INIT(tropical)(SModulFunctions* p)
{
  p->iiAddCproc("","permutationGroup",FALSE,permutationGroup);
  p->iiAddCproc("","minimalRepresentative",FALSE,minimalRepresentative);
  p->iiAddCproc("","valuation",FALSE,valuation);
  p->iiAddCproc("","tropicalStartingConeNewton",FALSE,tropicalStartingConeNewton);
  p->iiAddCproc("","tropicalDebug",FALSE,tropicalDebug);
  return MAX_TOK;
}
