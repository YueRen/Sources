#ifndef TROPICAL_CONVERSION_H
#define TROPICAL_CONVERSION_H

#include <vector>
#include <set>

#include <Singular/lists.h>


std::set<std::vector<int> > convertPermutations(lists generators0);
lists convertPermutations(std::set<std::vector<int> > permutationsGroup1);

int rationalToInt(const gfan::Rational &r);
gfan::Rational numberToRational(number c, coeffs cf);

gfan::ZVector intvecStarToZVector(intvec* v);
gfan::QMatrix matrixToQMatrix(matrix M, ring r);

#endif
