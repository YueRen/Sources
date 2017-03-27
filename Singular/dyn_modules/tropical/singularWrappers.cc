#include <Singular/libsingular.h>
#include <gfanlib/gfanlib.h>

#include <Singular/dyn_modules/gfanlib/callgfanlib_conversion.h>
#include <conversion.h>

gfan::QVector callTropicalPointNewton(ideal I, ring r)
{
  ring origin = currRing;
  if (origin != r)
  {
    std::cerr << "changing ring!!!" << std::endl;
    rChangeCurrRing(r);
  }

  idhdl h = ggetid("tropicalPointNewton");
  if ((h==NULL) || (h->typ!=PROC_CMD))
  {
    WerrorS("procedure tropicalPointNewton not found");
    return gfan::QVector();
  }
  sleftv args;
  args.Init();
  args.rtyp=IDEAL_CMD;
  args.data=(void*) id_Copy(I,r);
  args.next=NULL;

  BOOLEAN err = iiMake_proc(h,NULL,&args);
  args.CleanUp(currRing);
  if (err || iiRETURNEXPR.Typ()!=MATRIX_CMD)
  {
    WerrorS("calling tropicalPointNewton failed");
    return gfan::QVector();
  }

  matrix w = (matrix) iiRETURNEXPR.data;
  iiRETURNEXPR.data=NULL;
  int n = w->cols();
  gfan::QVector v(n);
  for (int i=1; i<=n; i++)
  {
    poly cp = MATELEM(w,1,i);
    if (cp!=NULL)
    {
      number cn = p_GetCoeff(cp,r);
      v[i-1] = numberToRational(cn,r->cf);
    }
    else
    {
      v[i-1] = gfan::Rational((long) 0);
    }
  }
  mp_Delete(&w,r);

  if (origin != r)
    rChangeCurrRing(origin);
  return v;
}


gfan::ZMatrix callTropicalLinkNewton(ideal I, ring r)
{
  idhdl originhdl = currRingHdl;
  ring origin = currRing;
  idhdl rhdl;
  if (origin != r)
  {
    rhdl= enterid("rInterpreter",0,RING_CMD,&IDROOT,FALSE);
    IDRING(rhdl) = rCopy(r);
    rSetHdl(rhdl);
  }

  idhdl h = ggetid("tropicalLinkNewton");
  if ((h==NULL) || (h->typ!=PROC_CMD))
  {
    WerrorS("procedure tropicalLinkNewton not found");
    return gfan::ZMatrix();
  }
  sleftv args;
  args.Init();
  args.rtyp=IDEAL_CMD;
  args.data=(void*) id_Copy(I,r);
  args.next=NULL;

  BOOLEAN err = iiMake_proc(h,NULL,&args);
  args.CleanUp(currRing);
  if (err || iiRETURNEXPR.Typ()!=LIST_CMD)
  {
    WerrorS("calling tropicalLinkNewton failed");
    return gfan::ZMatrix();
  }

  lists l = (lists) iiRETURNEXPR.data;
  iiRETURNEXPR.data=NULL;
  gfan::ZMatrix V(0,rVar(r));
  for (int i=0; i<=lSize(l); i++)
  {
    matrix v0 = (matrix) l->m[i].Data();
    gfan::QMatrix v1 = matrixToQMatrix(v0,r);
    gfan::ZVector v = QToZVectorPrimitive(v1[0].toVector());
    V.appendRow(v);
  }
  l->Clean(currRing);

  if (origin != r)
  {
    rSetHdl(originhdl);
    killhdl(rhdl);
  }
  return V;
}


ideal tropical_kStd_wrapper(ideal I, ring r, tHomog h=testHomog)
{
  ring origin = currRing;
  if (origin != r)
    rChangeCurrRing(r);

  ideal stdI = kStd(I,currRing->qideal,h,NULL);
  idSkipZeroes(stdI);

  if (origin != r)
    rChangeCurrRing(origin);
  return stdI;
}


ideal tropical_kNF_wrapper(ideal dividend, ring dividendRing, ideal divisor, ring divisorRing)
{
  ring origin = currRing;
  if (origin != divisorRing)
    rChangeCurrRing(divisorRing);

  nMapFunc identity = n_SetMap(dividendRing->cf,divisorRing->cf);
  int k = IDELEMS(dividend);
  ideal dividendInDivisorRing = idInit(k);
  for (int l=0; l<k; l++)
    dividendInDivisorRing->m[l] = p_PermPoly(dividend->m[l],NULL,dividendRing,divisorRing,identity,NULL,0);

  ideal residueInDivisorRing = kNF(divisor,divisorRing->qideal,dividendInDivisorRing);

  identity = n_SetMap(divisorRing->cf,dividendRing->cf);
  ideal residue = idInit(k);
  for (int l=0; l<k; l++)
    residue->m[l] = p_PermPoly(residueInDivisorRing->m[l],NULL,divisorRing,dividendRing,identity,NULL,0);

  if (origin != divisorRing)
    rChangeCurrRing(origin);
  return residue;
}
