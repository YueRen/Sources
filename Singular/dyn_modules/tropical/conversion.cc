#include <vector>
#include <set>

#include <libpolys/misc/intvec.h>
#include <Singular/lists.h>

std::set<std::vector<int> > convertPermutations(lists generators0)
{
  std::set<std::vector<int> > generators1;
  for (int i=0; i<=lSize(generators0); i++)
  {
    if (generators0->m[i].Typ() != INTVEC_CMD)
    {
      WerrorS("generatePermutationGroup: illegal entries in list");
      return std::set<std::vector<int> >();
    }
    intvec* generator0 = (intvec*) generators0->m[i].Data();

    std::vector<int> generator1(generator0->rows());
    for (int j=0; j<generator0->rows(); j++)
      generator1[j] = (*generator0)[j]-1;

    generators1.insert(generator1);
  }
  return generators1;
}


lists convertPermutations(std::set<std::vector<int> > permutationGroup1)
{
  lists permutationGroup0 = (lists) omAllocBin(slists_bin);
  permutationGroup0->Init(permutationGroup1.size()); int i=0;
  for (std::set<std::vector<int> >::iterator it = permutationGroup1.begin(); it!=permutationGroup1.end(); ++it, ++i)
  {
    intvec* permutation0 = new intvec(it->size());
    for (std::size_t j=0; j<it->size(); j++)
      (*permutation0)[j] = (*it)[j]+1;

    permutationGroup0->m[i].rtyp = INTVEC_CMD;
    permutationGroup0->m[i].data = (void*) permutation0;
  }
  return permutationGroup0;
}
