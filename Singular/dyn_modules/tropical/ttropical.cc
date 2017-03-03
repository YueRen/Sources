#include <kernel/mod2.h>

#include <Singular/ipid.h>
#include <Singular/mod_lib.h>

#include <Singular/dyn_modules/gfanlib/callgfanlib_conversion.h>

#include <conversion.h>
#include <symmetries.h>




BOOLEAN generatePermutationGroup(leftv res, leftv args)
{
  leftv u = args;
  if ((u != NULL) && (u->Typ() == LIST_CMD) && (u->next == NULL))
  {
    lists generators0 = (lists) u->Data();

    std::set<std::vector<int> > generators1 = convertPermutations(generators0);
    std::set<std::vector<int> > permutationGroup1 = generatePermutationGroup(generators1);
    lists permutationGroup0 = convertPermutations(permutationGroup1);

    res->rtyp = LIST_CMD;
    res->data = (void*) permutationGroup0;
    return FALSE;
  }
  WerrorS("generatePermutationGroup: unexpected parameters");
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



extern "C" int SI_MOD_INIT(tropical)(SModulFunctions* p)
{
  p->iiAddCproc("","generatePermutationGroup",FALSE,generatePermutationGroup);
  p->iiAddCproc("","minimalRepresentative",FALSE,minimalRepresentative);
  return MAX_TOK;
}
