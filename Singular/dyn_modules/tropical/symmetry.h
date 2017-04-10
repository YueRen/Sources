#ifndef TROPICAL_SYMMETRY_H
#define TROPICAL_SYMMETRY_H

#include <vector>
#include <set>

#include <gfanlib/gfanlib.h>

std::vector<int> composePermutations(const std::vector<int> &sigma, const std::vector<int> &tau);
std::set<std::vector<int> > permutationGroup(const std::set<std::vector<int> > &generators);

gfan::ZVector applyPermutation(const gfan::ZVector &w, const std::vector<int> &g);
gfan::ZVector minimalRepresentative(const gfan::ZVector &w, const std::set<std::vector<int> > &G);
std::set<gfan::ZVector> orbit(const gfan::ZVector &w, const std::set<std::vector<int> > &G);

#endif
