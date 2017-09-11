#ifndef TROPICAL_GROEBNERCONE_H
#define TROPICAL_GROEBNERCONE_H

#include <Singular/libsingular.h>
#include <gfanlib/gfanlib.h>

#include <set>
#include <vector>

namespace tropical {

class groebnerCone
{
private:
  ring polynomialRing;
  ideal polynomialIdeal;
  gfan::ZCone polyhedralCone;
  gfan::ZVector uniquePoint;

public:
  groebnerCone();
  groebnerCone(const groebnerCone& sigma);
  ~groebnerCone();
  groebnerCone& operator=(const groebnerCone& sigma);

	std::pair<gfan::ZVector,gfan::ZVector> facetPointingTo(const gfan::ZMatrix &P) const;

	groebnerCone(ideal I, ring r);
  groebnerCone(ideal I, ring r, gfan::ZVector interiorPoint, std::set<std::vector<int> > symmetryGroup);
  groebnerCone(ideal I, ideal inI, ring r, std::set<std::vector<int> > symmetryGroup=std::set<std::vector<int> >());

  ideal getPolynomialIdeal() const { return polynomialIdeal; };
  ring getPolynomialRing() const { return polynomialRing; };
  gfan::ZCone getPolyhedralCone() const { return polyhedralCone; };
  gfan::ZVector getUniquePoint() const { return uniquePoint; };

  bool operator<(const groebnerCone& sigma) const;

  void deletePolynomialIdealAndRing();
  void deletePolyhedralCone();
};

gfan::ZFan* groebnerConesToZFanStar(std::set<groebnerCone>& groebnerCones);
lists groebnerConesToListOfZCones(std::set<groebnerCone>& groebnerCones);
}

#endif
