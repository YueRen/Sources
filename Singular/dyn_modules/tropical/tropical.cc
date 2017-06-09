#include <kernel/mod2.h>

#include <Singular/libsingular.h>
#include <gfanlib/gfanlib.h>

#include <Singular/dyn_modules/gfanlib/callgfanlib_conversion.h>
#include <Singular/dyn_modules/gfanlib/bbcone.h>
#include <Singular/dyn_modules/gfanlib/bbfan.h>

#include <conversion.h>
#include <groebnerCone.h>
#include <homogeneity.h>
#include <singularWrappers.h>
#include <startingCone.h>
#include <symmetry.h>
#include <traversal.h>
#include <valuation.h>


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
    gfan::ZVector w = homogeneityVector(I,currRing);

    tropical::groebnerCone sigma = tropicalStartingCone(I,currRing,w);

    res->rtyp = coneID;
    res->data = (void*) new gfan::ZCone(sigma.getPolyhedralCone());
    return FALSE;
  }
  WerrorS("tropicalStartingCone: unexpected parameters");
  return TRUE;
}


BOOLEAN tropicalVarietyNew(leftv res, leftv args)
{
  leftv u = args;
  if ((u != NULL) && (u->Typ() == IDEAL_CMD))
  {
    ideal I = (ideal) u->Data();
    std::set<std::vector<int> > symmetryGroup1 = std::set<std::vector<int> >();

    leftv v = u->next;
    if ((v != NULL) && (v->Typ() == LIST_CMD))
    {
      lists symmetryGroup0 = (lists) v->Data();
      symmetryGroup1 = convertPermutations(symmetryGroup0);
    }

    gfan::ZVector w = homogeneityVector(I,currRing);

    tropical::groebnerCone sigma = tropicalStartingCone(I,currRing,w,symmetryGroup1);
    std::set<tropical::groebnerCone> TropI = tropicalTraversal(sigma,w,symmetryGroup1);

    res->rtyp = LIST_CMD;
    res->data = (void*) groebnerConesToListOfZCones(TropI);
    return FALSE;
  }
  WerrorS("tropicalVarietyNew: unexpected parameters");
  return TRUE;
}


std::vector<std::set<gfan::ZVector> > uniquePointsOfFacesByDimension(lists maximalCones)
{
  gfan::ZCone* maximalCone0 = (gfan::ZCone*) maximalCones->m[0].Data();
  std::vector<std::set<gfan::ZVector> > uniquePointsGroupedByDimension(maximalCone0->ambientDimension());

  for (int i=0; i<=lSize(maximalCones); i++)
  {
    gfan::ZCone* maximalCone = (gfan::ZCone*) maximalCones->m[i].Data();

    gfan::ZFan maximalConeFan = gfan::ZFan(maximalCone->ambientDimension());
    maximalConeFan.insert(*maximalCone);
    gfan::ZMatrix maximalConeRays = rays(&maximalConeFan);
    int ld = maximalConeFan.getLinealityDimension();

    for (int d=0; d<=maximalConeFan.getDimension()-ld; d++)
    {
      for (int j=0; j<maximalConeFan.numberOfConesOfDimension(d,0,0); j++)
      {
        gfan::IntVector rayIndices = maximalConeFan.getConeIndices(d,j,0,0);

        gfan::ZVector uniquePoint = gfan::ZVector(maximalCone->ambientDimension());
        for (int l=0; l<rayIndices.size(); l++)
        {
          uniquePoint = uniquePoint + maximalConeRays[rayIndices[l]].toVector();
        }
        uniquePointsGroupedByDimension[d+ld].insert(uniquePoint);
      }
    }
  }

  return uniquePointsGroupedByDimension;
}


BOOLEAN fVectorModuloSymmetry(leftv res, leftv args)
{
  leftv u = args;
  if ((u != NULL) && (u->Typ() == LIST_CMD))
  {
    leftv v = u->next;
    if ((u != NULL) && (u->Typ() == LIST_CMD))
    {
      lists maximalCones = (lists) u->Data();
      lists symmetryGroup0 = (lists) v->Data();
      std::set<std::vector<int> > symmetryGroup = convertPermutations(symmetryGroup0);

      gfan::ZCone* maximalCone = (gfan::ZCone*) maximalCones->m[0].Data();
      int d = maximalCone->dimension();
      int ld = maximalCone->dimensionOfLinealitySpace();

      bigintmat* fVector = new bigintmat(1,d+1,coeffs_BIGINT);
      number temp0 = n_Init(0,coeffs_BIGINT);
      for (int i=0; i<ld; i++)
        fVector->set(1,i+1,temp0);
      n_Delete(&temp0,coeffs_BIGINT);
      number temp1 = n_Init(1,coeffs_BIGINT);
      fVector->set(1,ld+1,temp1);
      n_Delete(&temp1,coeffs_BIGINT);

      std::vector<std::set<gfan::ZVector> > uniquePointsGroupedByDimension = uniquePointsOfFacesByDimension(maximalCones);
      for (int i=ld+1; i<=d; i++)
      {
        std::set<gfan::ZVector> uniquePoints = uniquePointsGroupedByDimension[i];
        std::set<gfan::ZVector> uniquePointsModuloSymmetry;
        for (std::set<gfan::ZVector>::iterator uniquePoint = uniquePoints.begin(); uniquePoint!=uniquePoints.end(); ++uniquePoint)
        {
          gfan::ZVector uniquePointMinimalRepresentative = minimalRepresentative(*uniquePoint,symmetryGroup);
          if (uniquePointsModuloSymmetry.count(uniquePointMinimalRepresentative)==0)
            uniquePointsModuloSymmetry.insert(uniquePointMinimalRepresentative);
        }

        number temp = integerToNumber(uniquePointsModuloSymmetry.size());
        fVector->set(1,i+1,temp);
        n_Delete(&temp,coeffs_BIGINT);
      }

      res->rtyp = BIGINTMAT_CMD;
      res->data = (void*) fVector;
      return FALSE;
    }
  }
  WerrorS("fVectorModuloSymmetry: unexpected parameters");
  return TRUE;
}


BOOLEAN fVectorIgnoringSymmetry(leftv res, leftv args)
{
  leftv u = args;
  if ((u != NULL) && (u->Typ() == LIST_CMD))
  {
    leftv v = u->next;
    if ((u != NULL) && (u->Typ() == LIST_CMD))
    {
      lists maximalCones = (lists) u->Data();
      lists symmetryGroup0 = (lists) v->Data();
      std::set<std::vector<int> > symmetryGroup = convertPermutations(symmetryGroup0);

      gfan::ZCone* maximalCone = (gfan::ZCone*) maximalCones->m[0].Data();
      int d = maximalCone->dimension();
      int ld = maximalCone->dimensionOfLinealitySpace();

      bigintmat* fVector = new bigintmat(1,d+1,coeffs_BIGINT);
      number temp0 = n_Init(0,coeffs_BIGINT);
      for (int i=0; i<ld; i++)
        fVector->set(1,i+1,temp0);
      n_Delete(&temp0,coeffs_BIGINT);
      number temp1 = n_Init(1,coeffs_BIGINT);
      fVector->set(1,ld+1,temp1);
      n_Delete(&temp1,coeffs_BIGINT);

      std::vector<std::set<gfan::ZVector> > uniquePointsGroupedByDimension = uniquePointsOfFacesByDimension(maximalCones);
      for (int i=ld+1; i<=d; i++)
      {
        std::set<gfan::ZVector> uniquePoints = uniquePointsGroupedByDimension[i];
        std::set<gfan::ZVector> uniquePointsIgnoringSymmetry;
        for (std::set<gfan::ZVector>::iterator uniquePoint = uniquePoints.begin(); uniquePoint!=uniquePoints.end(); ++uniquePoint)
        {
          std::set<gfan::ZVector> uniquePointOrbit = orbit(*uniquePoint,symmetryGroup);
          uniquePointsIgnoringSymmetry.insert(uniquePointOrbit.begin(),uniquePointOrbit.end());
        }

        number temp = integerToNumber(uniquePointsIgnoringSymmetry.size());
        fVector->set(1,i+1,temp);
        n_Delete(&temp,coeffs_BIGINT);
      }

      res->rtyp = BIGINTMAT_CMD;
      res->data = (void*) fVector;
      return FALSE;
    }
  }
  WerrorS("fVectorIgnoringSymmetry: unexpected parameters");
  return TRUE;
}


BOOLEAN tropicalDebug(leftv res, leftv args)
{
  res->rtyp = NONE;
  res->data = NULL;
  return FALSE;
}

extern "C" int SI_MOD_INIT(tropical)(SModulFunctions* p)
{
  p->iiAddCproc("","permutationGroup",FALSE,permutationGroup);
  p->iiAddCproc("","minimalRepresentative",FALSE,minimalRepresentative);
  p->iiAddCproc("","valuation",FALSE,valuation);
  p->iiAddCproc("","tropicalStartingConeNewton",FALSE,tropicalStartingConeNewton);
  p->iiAddCproc("","tropicalVarietyNew",FALSE,tropicalVarietyNew);
  p->iiAddCproc("","fVectorModuloSymmetry",FALSE,fVectorModuloSymmetry);
  p->iiAddCproc("","fVectorIgnoringSymmetry",FALSE,fVectorIgnoringSymmetry);
  p->iiAddCproc("","tropicalDebug",FALSE,tropicalDebug);
  return MAX_TOK;
}
