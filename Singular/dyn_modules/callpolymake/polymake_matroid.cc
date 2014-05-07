#include <polymake_conversion.h>

#include <Singular/blackbox.h>
#include <Singular/ipid.h>
#include <Singular/ipshell.h>

#include <Singular/dyn_modules/callgfanlib/bbfan.h>

int MATROID_CMD;

void* matroid_Init(blackbox* /*b*/)
{
  polymake::perl::Object* defaultMatroid = new polymake::perl::Object("Matroid");
  return (void*) defaultMatroid;
}

void matroid_Destroy(blackbox* /*b*/, void* toBeDeleted)
{
  if (toBeDeleted!=NULL)
  {
    polymake::perl::Object* deleteThis = (polymake::perl::Object*) toBeDeleted;
    delete deleteThis;
  }
}

void* matroid_Copy(blackbox* /*b*/, void* toBeCopied)
{
  if (toBeCopied!=NULL)
  {
    polymake::perl::Object* copyThis = (polymake::perl::Object*) toBeCopied;
    return new polymake::perl::Object(*copyThis);
  }
  return NULL;
}

BOOLEAN matroid_Assign(leftv leftSide, leftv rightSide)
{
  if (leftSide->Data()!=NULL)
  {
    polymake::perl::Object* oldDataOnLeftSide = (polymake::perl::Object*) leftSide->Data();
    delete oldDataOnLeftSide;
  }

  polymake::perl::Object* toBeAssigned;
  if (rightSide==NULL)
    toBeAssigned = new polymake::perl::Object("Matroid");
  else if (rightSide->Typ()==leftSide->Typ())
  {
    polymake::perl::Object* dataOnRightSide = (polymake::perl::Object*) rightSide->Data();
    toBeAssigned = new polymake::perl::Object(*dataOnRightSide);
  }
  else
  {
    WerrorS("matroid: unsupported assignment");
    toBeAssigned = new polymake::perl::Object("Matroid");
  }

  if (leftSide->rtyp==IDHDL)
    IDDATA((idhdl)leftSide->data) = (char*) toBeAssigned;
  else
    leftSide->data = (void*) toBeAssigned;
  return FALSE;
}


static std::string toString(polymake::Array<polymake::Set<int> > arrayOfSetInt)
{
  std::stringstream s;
  for (int i=0; i<arrayOfSetInt.size(); i++)
  {
    polymake::Vector<int> setInt(arrayOfSetInt[i]);
    for (int j=0; j<setInt.size()-1; j++)
      s << setInt[j] << ",";
    s << setInt[setInt.size()-1];
    s << std::endl;
  }
  return s.str();
}


std::string matroidToString(polymake::perl::Object* toBePrinted)
{
  std::stringstream output;

  int n_elements = toBePrinted->give("N_ELEMENTS");
  output << "N_ELEMENTS" << std::endl;
  output << n_elements << std::endl;

  polymake::Array<polymake::Set<int> > bases = toBePrinted->give("BASES");
  output << "BASES" << std::endl;
  output << toString(bases) << std::endl;

  return output.str();
}

char* matroid_String(blackbox* /*b*/, void* toBePrinted)
{
  if (toBePrinted==NULL) return omStrDup("invalid object");
  std::string printThis = matroidToString((polymake::perl::Object*) toBePrinted);
  return omStrDup(printThis.c_str());
}

static bool isListOfIntvecs(lists L)
{
  for (int i=0; i<=lSize(L); i++)
    if (L->m[i].Typ()!=INTVEC_CMD) return false;
  return true;
}

BOOLEAN matroidViaCircuits(leftv res, leftv args)
{
  leftv u = args;
  if ((u!=NULL) && (u->Typ()==LIST_CMD))
  {
    lists L = (lists) u->Data();
    leftv v = u->next;
    if ((v!=NULL) && (v->Typ()==INT_CMD))
    {
      if (v->next==NULL && isListOfIntvecs(L))
      {
        int n_elements = (int) (long) v->Data();
        polymake::Array<polymake::Set<int> > circuits = listOfIntvecsToArrayOfSetInt(L);
        polymake::perl::Object* matroid = new polymake::perl::Object("Matroid");
        matroid->take("CIRCUITS") << circuits;
        matroid->take("N_ELEMENTS") << n_elements;
        res->rtyp = MATROID_CMD;
        res->data = (polymake::perl::Object*) matroid;
        return FALSE;
      }
    }
  }
  WerrorS("matroidViaCircuits: unexpected parameters");
  return TRUE;
}

BOOLEAN bergmanFanMatroid(leftv res, leftv args)
{
	leftv u = args;
	if ((u != NULL) && (u->Typ() == MATROID_CMD))
	{
		leftv v = u->next;
		if ((v != NULL) && (v->Typ() == INT_CMD))
		{
			leftv w = v->next;
			if ((w != NULL) && (w->Typ() == INT_CMD))
			{
				polymake::perl::Object* matroid = (polymake::perl::Object*) u->Data();
				int modOutLin = (int) (long) v->Data();
				int coord = (int) (long) w->Data();
				gfan::ZFan* BFan;
				try
				{
					polymake::perl::Object Fan;
					CallPolymakeFunction("bergman_fan_matroid",*matroid,modOutLin,coord) >> Fan;
					BFan = PmFan2ZFan(&Fan);
				}
				catch (const std::exception& ex)
				{
					WerrorS("ERROR: "); WerrorS(ex.what()); WerrorS("\n");
					return TRUE;
				}
				res->rtyp = fanID;
				res->data = (char*) BFan;
				return FALSE;
			}
		}
		else if ( v == NULL )
		{
			polymake::perl::Object* matroid = (polymake::perl::Object*) u->Data();
			gfan::ZFan* BFan;
			try
			{
				polymake::perl::Object Fan;
				CallPolymakeFunction("bergman_fan_matroid", *matroid) >> Fan;
				BFan = PmFan2ZFan(&Fan);
			}
			catch (const std::exception& ex)
			{
				WerrorS("ERROR: "); WerrorS(ex.what()); WerrorS("\n");
				return TRUE;
			}
			res->rtyp = fanID;
			res->data = (char*) BFan;
			return FALSE;
		}
	}
	WerrorS("bergmanFanMatroid: unexpected parameters");
	return TRUE;
}


void matroid_setup(SModulFunctions* p)
{
  blackbox* b = (blackbox*) omAlloc0(sizeof(blackbox));
  // all undefined entries will be set to default in setBlackboxStuff
  // the default Print is quite usefull,
  // all other are simply error messages
  b->blackbox_Init=matroid_Init;
  b->blackbox_destroy=matroid_Destroy;
  b->blackbox_Copy=matroid_Copy;
  b->blackbox_Assign=matroid_Assign;
  b->blackbox_String=matroid_String;

  p->iiAddCproc("polymake.so","matroidViaCircuits",FALSE,matroidViaCircuits);
	p->iiAddCproc("polymake.so","bergmanFanMatroid",FALSE,bergmanFanMatroid);

  MATROID_CMD = setBlackboxStuff(b,"matroid");
}


