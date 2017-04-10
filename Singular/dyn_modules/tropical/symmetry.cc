#include <map>
#include <set>
#include <vector>
#include <iostream>

#include <gfanlib/gfanlib.h>

std::vector<int> composePermutations(const std::vector<int> &sigma, const std::vector<int> &tau)
{
  std::vector<int> sigmaTau(sigma.size());

  for (std::size_t i=0; i<sigma.size(); i++)
    sigmaTau[i] = sigma[tau[i]];

  return sigmaTau;
}


std::set<std::vector<int> > permutationGroup(const std::set<std::vector<int> > &generators)
{
  std::map<std::vector<int>,bool> generatingSubgroup;

  for (std::set<std::vector<int> >::iterator it = generators.begin(); it!=generators.end(); ++it)
    generatingSubgroup.insert(std::make_pair(*it,true));

  for (std::map<std::vector<int>,bool>::iterator it = generatingSubgroup.begin(); it!=generatingSubgroup.end(); ++it)
  {
    if (it->second)
    {
      it->second = false;
      for (std::set<std::vector<int> >::iterator jt = generators.begin(); jt!=generators.end(); ++jt)
        generatingSubgroup.insert(std::make_pair(composePermutations(it->first,*jt),true));
      it = generatingSubgroup.begin();
      std::cout << "permutations found: " << generatingSubgroup.size() << std::endl;
    }
  }
  std::set<std::vector<int> > subgroup;
  for (std::map<std::vector<int>,bool>::iterator it = generatingSubgroup.begin(); it!=generatingSubgroup.end(); ++it)
    subgroup.insert(it->first);

  return subgroup;
}


gfan::ZVector applyPermutation(const gfan::ZVector &w, const std::vector<int> &g)
{
  gfan::ZVector gw(w.size());

  for (std::size_t i=0; i<gw.size(); i++)
    gw[i] = w[g[i]];

  return gw;
}

gfan::ZVector minimalRepresentative(const gfan::ZVector &w, const std::set<std::vector<int> > &G)
{
  gfan::ZVector wmin = w;

  for (std::set<std::vector<int> >::iterator g=G.begin(); g!=G.end(); ++g)
  {
    gfan::ZVector gw = applyPermutation(w,*g);
    if (gw < wmin)
      wmin = gw;
  }

  return wmin;
}

std::set<gfan::ZVector> orbit(const gfan::ZVector &w, const std::set<std::vector<int> > &G)
{
  std::set<gfan::ZVector> Gw;

  for (std::set<std::vector<int> >::iterator g=G.begin(); g!=G.end(); ++g)
    Gw.insert(applyPermutation(w,*g));

  return Gw;
}
