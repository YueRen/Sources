#include <bbcone.h>
#include <tropicalPrevarieties.h>
#include <callgfanlib_conversion.h>
#include <Singular/ipshell.h> // for printlevels
#include <Singular/lists.h>

markedCone::markedCone(gfan::ZCone &zc):
  polyhedralCone(),
  interiorPoint()
{
  zc.canonicalize();
  polyhedralCone = zc;
  interiorPoint = zc.getRelativeInteriorPoint();
}


markedCone::markedCone(const gfan::ZCone &zc, const gfan::ZVector &p):
  polyhedralCone(zc),
  interiorPoint(p)
{
}


tropicalPrevariety::tropicalPrevariety():
  maximalCones(),
  ambientDimension(0),
  dimension(0),
  expectedDimension(0)
{
}


tropicalPrevariety::tropicalPrevariety(int ad, int d, int ed):
  maximalCones(),
  ambientDimension(ad),
  dimension(d),
  expectedDimension(ed)
{
}


tropicalPrevariety::tropicalPrevariety(gfan::ZCone &fundamentalDomain)
{
  maximalCones.insert(markedCone(fundamentalDomain));
  ambientDimension = fundamentalDomain.ambientDimension();
  dimension = fundamentalDomain.dimension();
  expectedDimension = fundamentalDomain.ambientDimension();
}


tropicalPrevariety::tropicalPrevariety(const poly g, const ring r, const gfan::ZCone* fundamentalDomain=NULL)
{
  int n = rVar(r);
  ambientDimension = n;
  dimension = n-1;
  expectedDimension = n-1;

  if (g!=NULL && g->next!=NULL)
  {
    // 1. compute all vertices of the newton polytope
    gfan::ZMatrix homogenizedExponents = gfan::ZMatrix(pLength(g),n+1);
    int l = 0;
    for (poly s=g; s!=NULL; pIter(s))
    {
      homogenizedExponents[l][0] = 1;
      for (int i=1; i<=n; i++)
        homogenizedExponents[l][i] = p_GetExp(s,i,r);
      l++;
    }

    gfan::ZCone homogenizedNewtonPolytope = gfan::ZCone(homogenizedExponents,gfan::ZMatrix(0,n+1));
    gfan::ZMatrix homogenizedVertices = homogenizedNewtonPolytope.getFacets();
    l = homogenizedVertices.getHeight();
    gfan::ZMatrix vertices = gfan::ZMatrix(l,n);
    for (int i=0; i<l; i++)
      for (int j=0; j<n; j++)
        vertices[i][j] = homogenizedVertices[i][j+1];


    for (int i=0; i<l; i++)
    {
      // 2. for each vertex, compute the corresponding maximal Groebner cone and its facets
      gfan::ZMatrix inequalities = gfan::ZMatrix(0,n);
      for (int j=0; j<l; j++)
        if (j!=i) inequalities.appendRow(vertices[i]-vertices[j]);
      gfan::ZCone maximalGroebnerCone = gfan::ZCone(inequalities,gfan::ZMatrix(0,n),gfan::PCP_impliedEquationsKnown);
      if (fundamentalDomain!=NULL)
      {
        maximalGroebnerCone = gfan::intersection(maximalGroebnerCone,*fundamentalDomain);
        if (maximalGroebnerCone.dimension()<n) continue;
      }

      gfan::ZMatrix maximalGroebnerConeFacets = maximalGroebnerCone.getFacets();

      // 3. for each maximal Groebner cone, add every facet to
      //    the list of maximal cones in the tropical prevariety
      int k = maximalGroebnerConeFacets.getHeight();
      for (int j=0; j<k; j++)
      {
        gfan::ZMatrix tropicalEquation = gfan::ZMatrix(0,n);
        tropicalEquation.appendRow(maximalGroebnerConeFacets[j]);
        gfan::ZCone tropicalCone = gfan::ZCone(maximalGroebnerConeFacets,tropicalEquation,gfan::PCP_impliedEquationsKnown);
        tropicalCone.canonicalize();
        maximalCones.insert(markedCone(tropicalCone));
      }
    }
  }
}


void tropicalPrevariety::deleteNonMaximalCones()
{
  // 1. iterate over all cones of lower dimension
  for (markedCones::iterator sigma = maximalCones.begin(); sigma!=maximalCones.end(); )
  {
    int d = sigma->getPolyhedralCone().dimension();
    if (d==dimension) break;
    bool sigmaToBeDeleted = false;

    // 2. iterate over all cones of higher dimension
    for (markedCones::reverse_iterator delta = maximalCones.rbegin(); delta!=maximalCones.rend(); delta++)
    {
      if (d>=delta->getPolyhedralCone().dimension())
        break;

      // 3. if delta contains sigma, mark sigma as to be deleted
      if (delta->getPolyhedralCone().contains(sigma->getInteriorPoint()))
      {
        sigmaToBeDeleted = true;
        break;
      }
    }

    // 4. delete sigma if necessary
    if (sigmaToBeDeleted)
      maximalCones.erase(sigma++);
    else
      sigma++;
  }
}


tropicalPrevariety intersectTropicalPrevarieties(const tropicalPrevariety &Sigma, const tropicalPrevariety &Delta)
{
  tropicalPrevariety SigmaCapDelta;
  SigmaCapDelta.ambientDimension = Sigma.ambientDimension;
  int d = Sigma.expectedDimension + Delta.expectedDimension - SigmaCapDelta.ambientDimension;
  if (d>0)
    SigmaCapDelta.expectedDimension = d;
  else
    SigmaCapDelta.expectedDimension = 0;

  for (markedCones::iterator sigma = Sigma.maximalCones.begin(); sigma!=Sigma.maximalCones.end(); sigma++)
  {
    for (markedCones::iterator delta = Delta.maximalCones.begin(); delta!=Delta.maximalCones.end(); delta++)
    {
      gfan::ZCone zc = gfan::intersection(sigma->getPolyhedralCone(), delta->getPolyhedralCone());
      if (SigmaCapDelta.expectedDimension<=zc.dimension())
        SigmaCapDelta.maximalCones.insert(markedCone(zc));
    }
  }

  markedCones::reverse_iterator tau = SigmaCapDelta.maximalCones.rbegin();
  SigmaCapDelta.dimension = tau->getPolyhedralCone().dimension();

  return SigmaCapDelta;
}


lists tropicalPrevarietyToListOfMarkedCones(const tropicalPrevariety &Sigma)
{
  lists listOfMarkedCones = (lists) omAllocBin(slists_bin);
  listOfMarkedCones->Init(Sigma.maximalCones.size()+1);

  lists listOfDimensionalInformation = (lists) omAllocBin(slists_bin);
  listOfDimensionalInformation->Init(3);
  listOfDimensionalInformation->m[0].rtyp = INT_CMD;
  listOfDimensionalInformation->m[0].data = (void*) (long) Sigma.ambientDimension;
  listOfDimensionalInformation->m[1].rtyp = INT_CMD;
  listOfDimensionalInformation->m[1].data = (void*) (long) Sigma.dimension;
  listOfDimensionalInformation->m[2].rtyp = INT_CMD;
  listOfDimensionalInformation->m[2].data = (void*) (long) Sigma.expectedDimension;
  listOfMarkedCones->m[0].rtyp = LIST_CMD;
  listOfMarkedCones->m[0].data = (void*) listOfDimensionalInformation;

  int i = 1;
  for (markedCones::iterator sigma = Sigma.maximalCones.begin(); sigma!=Sigma.maximalCones.end(); sigma++)
  {
    lists markedCone = (lists) omAllocBin(slists_bin);
    markedCone->Init(2);
    markedCone->m[0].rtyp = coneID;
    markedCone->m[0].data = (void*) new gfan::ZCone(sigma->getPolyhedralCone());
    markedCone->m[1].rtyp = BIGINTMAT_CMD;
    markedCone->m[1].data = (void*) zVectorToBigintmat(sigma->getInteriorPoint());

    listOfMarkedCones->m[i].rtyp = LIST_CMD;
    listOfMarkedCones->m[i].data = (void*) markedCone;
    i++;
  }

  return listOfMarkedCones;
}


tropicalPrevariety listOfMarkedConesToTropicalPrevariety(const lists listOfMarkedCones)
{
  lists listOfDimensionalInformation = (lists) listOfMarkedCones->m[0].Data();
  int ad = (int)(long) listOfDimensionalInformation->m[0].Data();
  int d = (int)(long) listOfDimensionalInformation->m[1].Data();
  int ed = (int)(long) listOfDimensionalInformation->m[2].Data();
  tropicalPrevariety Sigma(ad,d,ed);

  for (int i=1; i<lSize(listOfMarkedCones); i++)
  {
    lists listOfMarkedCone = (lists) listOfMarkedCones->m[i].Data();
    gfan::ZCone* zc = (gfan::ZCone*) listOfMarkedCone->m[0].Data();
    bigintmat* p0 = (bigintmat*) listOfMarkedCone->m[1].Data();
    gfan::ZVector* p = bigintmatToZVector(p0);
    Sigma.maximalCones.insert(markedCone(*zc,*p));
    delete p;
  }
  return Sigma;
}


BOOLEAN gfanlib_tropicalPrevariety(leftv res, leftv args)
{
  leftv u = args;
  if ((u!=NULL) && (u->Typ()==POLY_CMD))
  {
    leftv v = u->next;

    poly g = (poly) u->Data();
    gfan::ZCone* fundamentalDomain=NULL;
    if ((v!=NULL) && (v->Typ()==coneID))
      fundamentalDomain = (gfan::ZCone*) v->Data();
    tropicalPrevariety Tg(g,currRing,fundamentalDomain);

    res->rtyp = LIST_CMD;
    res->data = (void*) tropicalPrevarietyToListOfMarkedCones(Tg);
    return FALSE;
  }
  // if ((u!=NULL) && (u->Typ()==IDEAL_CMD))
  // {
  //   leftv v = u->next;
  //   gfan::ZCone fundamentalDomain;
  //   if ((v!=NULL) && (v->Typ()==coneID))
  //     fundamentalDomain = *((gfan::ZCone*) v->Data());
  //   else if (v==NULL)
  //     fundamentalDomain = gfan::ZCone(rVar(currRing));   // = full space of specified dimension
  //   else
  //   {
  //     WerrorS("tropicalPrevariety: unexpected parameters");
  //     return TRUE;
  //   }
  //   ideal I = (ideal) u->Data();
  //   tropicalPrevariety Sigma(fundamentalDomain);
  //   for (int i=0; i<idSize(I); i++)
  //   {
  //     if (printlevel > 0)
  //       Print("tropicalPrevariety: parsing generator %d of %d", i+1, idSize(I)+1);
  //     Sigma = intersectTropicalPrevarieties(Sigma,tropicalPrevariety(I->m[i],currRing));
  //     Sigma.deleteNonMaximalCones();
  //     if (printlevel > 0)
  //       Print("  number of currently maximal cones: %lu\n",(unsigned long)Sigma.maximalCones.size());
  //   }
  //   res->rtyp = LIST_CMD;
  //   res->data = (void*) tropicalPrevarietyToListOfCones(Sigma);
  //   return FALSE;
  // }
  WerrorS("tropicalPrevariety: unexpected parameters");
  return TRUE;
}


BOOLEAN gfanlib_intersectTropicalPrevarieties(leftv res, leftv args)
{
  leftv u = args;
  if ((u!=NULL) && (u->Typ()==LIST_CMD))
  {
    leftv v = u->next;
    if ((u!=NULL) && (u->Typ()==LIST_CMD))
    {
      lists Sigma0 = (lists) u->Data();
      lists Delta0 = (lists) v->Data();
      tropicalPrevariety Sigma = listOfMarkedConesToTropicalPrevariety(Sigma0);
      tropicalPrevariety Delta = listOfMarkedConesToTropicalPrevariety(Delta0);

      tropicalPrevariety SigmaCapDelta = intersectTropicalPrevarieties(Sigma,Delta);
      SigmaCapDelta.deleteNonMaximalCones();

      res->rtyp = LIST_CMD;
      res->data = (void*) tropicalPrevarietyToListOfMarkedCones(SigmaCapDelta);
      return FALSE;
    }
  }
  WerrorS("intersectTropicalPrevarieties: unexpected parameters");
  return TRUE;
}
