#include <sstream>

#include <Singular/blackbox.h>
#include <Singular/ipid.h>
#include <Singular/ipshell.h>

#include <gfanlib/gfanlib.h>
#include <bbmarkedConeCollection.h>
#include <bbcone.h>
#include <callgfanlib_conversion.h>

int markedConeCollection_CMD;


markedCone::markedCone(gfan::ZCone& zc)
{
  zc.canonicalize();
  polyhedralCone = zc;
  interiorPoint = zc.getRelativeInteriorPoint();
}


bool markedCone::contains(const gfan::ZVector& p) const
{
  return polyhedralCone.contains(p);
}


void* markedConeCollection_Init(blackbox* /*b*/)
{
  return (void*) new markedConeCollection();
}


void* markedConeCollection_Copy(blackbox* /*b*/, void *d)
{
  markedConeCollection* Sigma = (markedConeCollection*)d;
  markedConeCollection* SigmaCopy = new markedConeCollection(*Sigma);
  return SigmaCopy;
}


void markedConeCollection_destroy(blackbox* /*b*/, void *d)
{
  if (d!=NULL)
  {
    markedConeCollection* Sigma = (markedConeCollection*) d;
    delete Sigma;
  }
}


BOOLEAN markedConeCollection_Assign(leftv l, leftv r)
{
  markedConeCollection* Sigma;
  if (r==NULL)
  {
    if (l->Data()!=NULL)
    {
      markedConeCollection* Delta = (markedConeCollection*) l->Data();
      delete Delta;
    }
    Sigma = new markedConeCollection();
  }
  else if (r->Typ()==l->Typ())
  {
    if (l->Data()!=NULL)
    {
      markedConeCollection* Delta = (markedConeCollection*) l->Data();
      delete Delta;
    }
    Sigma = (markedConeCollection*) r->CopyD();
  }
  else
  {
    Werror("assign Type(%d) = Type(%d) not implemented",l->Typ(),r->Typ());
    return TRUE;
  }

  if (l->rtyp==IDHDL)
    IDDATA((idhdl)l->data) = (char*) Sigma;
  else
    l->data = (void*) Sigma;

  return FALSE;
}


char* markedConeCollection_String(blackbox* /*b*/, void *d)
{
  if (d==NULL)
    return omStrDup("invalid object");
  else
  {
    markedConeCollection* Sigma = (markedConeCollection*) d;
    std::stringstream printout;
    printout << "set of "
             << Sigma->size()
             << " marked cones";
    std::string printoutString = printout.str();
    return omStrDup(printoutString.c_str());
  }
}


bool markedConeCollection::isEmpty() const
{
  return setOfMarkedCones.empty();
}


BOOLEAN mccIsEmpty(leftv res, leftv args)
{
  leftv u = args;
  if ((u != NULL) && (u->Typ() == markedConeCollection_CMD) && (u->next == NULL))
  {
    markedConeCollection* Sigma = (markedConeCollection*)u->Data();
    res->rtyp = INT_CMD;
    res->data = (void*) (long) Sigma->isEmpty();
    return FALSE;
  }
  WerrorS("mccIsEmpty: unexpected parameters");
  return TRUE;
}


int markedConeCollection::size() const
{
  return (int) setOfMarkedCones.size();
}


BOOLEAN mccSize(leftv res, leftv args)
{
  leftv u = args;
  if ((u != NULL) && (u->Typ() == markedConeCollection_CMD) && (u->next == NULL))
  {
    markedConeCollection* Sigma = (markedConeCollection*)u->Data();
    res->rtyp = INT_CMD;
    res->data = (void*) (long) Sigma->size();
    return FALSE;
  }
  WerrorS("mccSize: unexpected parameters");
  return TRUE;
}


void markedConeCollection::insert(gfan::ZCone& zc)
{
  setOfMarkedCones.insert(markedCone(zc));
}


BOOLEAN mccInsertCone(leftv res, leftv args)
{
  leftv u = args;
  if ((u != NULL) && (u->Typ() == markedConeCollection_CMD))
  {
    leftv v = u->next;
    if ((v != NULL) && (v->Typ() == coneID) && (v->next == NULL))
    {
      markedConeCollection* Sigma = (markedConeCollection*)u->Data();
      gfan::ZCone* zc = (gfan::ZCone*)v->Data();
      zc->canonicalize();
      Sigma->insert(*zc);

      res->rtyp = NONE;
      res->data = NULL;
      return FALSE;
    }
  }
  WerrorS("mccInsertCone: unexpected parameters");
  return TRUE;
}


void markedConeCollection::remove(gfan::ZCone& zc)
{
  setOfMarkedCones.erase(markedCone(zc));
}


BOOLEAN mccRemoveCone(leftv res, leftv args)
{
  leftv u = args;
  if ((u != NULL) && (u->Typ() == markedConeCollection_CMD))
  {
    leftv v = u->next;
    if ((v != NULL) && (v->Typ() == coneID) && (v->next == NULL))
    {
      markedConeCollection* Sigma = (markedConeCollection*)u->Data();
      gfan::ZCone* zc = (gfan::ZCone*)v->Data();
      zc->canonicalize();
      Sigma->remove(*zc);

      res->rtyp = NONE;
      res->data = NULL;
      return FALSE;
    }
  }
  WerrorS("mccInsertCone: unexpected parameters");
  return TRUE;
}


lists markedConeCollection::getListOfCones() const
{
  lists listOfCones = (lists)omAllocBin(slists_bin);
  listOfCones->Init(setOfMarkedCones.size());
  int i=0;
  for (std::set<markedCone,markedCone_compare>::iterator it = setOfMarkedCones.begin(); it!=setOfMarkedCones.end(); it++)
  {
    listOfCones->m[i].rtyp = coneID;
    listOfCones->m[i].data = (void*) new gfan::ZCone(it->getPolyhedralCone());
  }
  return listOfCones;
}


BOOLEAN mccGetListOfCones(leftv res, leftv args)
{
  leftv u = args;
  if ((u != NULL) && (u->Typ() == markedConeCollection_CMD) && (u->next == NULL))
  {
    markedConeCollection* Sigma = (markedConeCollection*)u->Data();
    res->rtyp = LIST_CMD;
    res->data = (void*) Sigma->getListOfCones();
    return FALSE;
  }
  WerrorS("mccGetListOfCones: unexpected parameters");
  return TRUE;
}


bool markedConeCollection::hasConeContaining(const gfan::ZVector& p) const
{
  for (std::set<markedCone,markedCone_compare>::iterator it = setOfMarkedCones.begin(); it!=setOfMarkedCones.end(); it++)
  {
    if (it->contains(p))
      return true;
  }
  return false;
}


BOOLEAN mccHasConeContaining(leftv res, leftv args)
{
  leftv u = args;
  if ((u != NULL) && (u->Typ() == markedConeCollection_CMD))
  {
    leftv v = u->next;
    if ((v != NULL) && (v->Typ() == BIGINTMAT_CMD) && (v->next == NULL))
    {
      markedConeCollection* Sigma = (markedConeCollection*)u->Data();
      bigintmat* pBim = (bigintmat*) v->Data();
      gfan::ZVector* p = bigintmatToZVector(pBim);
      res->rtyp = INT_CMD;
      res->data = (void*) (long) Sigma->hasConeContaining(*p);
      delete p;
      return FALSE;
    }
  }
  WerrorS("mccHasConeContaining: unexpected parameters");
  return TRUE;
}


markedCone markedConeCollection::getConeAndDelete()
{
  std::set<markedCone,markedCone_compare>::iterator it = setOfMarkedCones.begin();
  markedCone sigma = *it;
  setOfMarkedCones.erase(it);
  return sigma;
}


BOOLEAN mccGetConeAndDelete(leftv res, leftv args)
{
  leftv u = args;
  if ((u != NULL) && (u->Typ() == markedConeCollection_CMD) && (u->next == NULL))
  {
    markedConeCollection* Sigma = (markedConeCollection*)u->Data();
    markedCone sigma = Sigma->getConeAndDelete();
    res->rtyp = coneID;
    res->data = (void*) new gfan::ZCone(sigma.getPolyhedralCone());
    return FALSE;
  }
  WerrorS("mccGetConeAndDelete: unexpected parameters");
  return TRUE;
}


void bbmarkedConeCollection_setup(SModulFunctions* p)
{
  blackbox *b=(blackbox*)omAlloc0(sizeof(blackbox));
  // all undefined entries will be set to default in setBlackboxStuff
  // the default Print is quite usefule,
  // all other are simply error messages
  b->blackbox_Init=markedConeCollection_Init;
  b->blackbox_Copy=markedConeCollection_Copy;
  b->blackbox_destroy=markedConeCollection_destroy;
  b->blackbox_Assign=markedConeCollection_Assign;
  b->blackbox_String=markedConeCollection_String;
  //b->blackbox_Print=blackbox_default_Print;
  p->iiAddCproc("","mccIsEmpty",FALSE,mccIsEmpty);
  p->iiAddCproc("","mccSize",FALSE,mccSize);
  p->iiAddCproc("","mccInsertCone",FALSE,mccInsertCone);
  p->iiAddCproc("","mccRemoveCone",FALSE,mccRemoveCone);
  p->iiAddCproc("","mccGetListOfCones",FALSE,mccGetListOfCones);
  p->iiAddCproc("","mccHasConeContaining",FALSE,mccHasConeContaining);
  p->iiAddCproc("","mccGetConeAndDelete",FALSE,mccGetConeAndDelete);
  markedConeCollection_CMD = setBlackboxStuff(b,"markedConeCollection");
  //Print("created type %d (fan)\n",fanID);
}
