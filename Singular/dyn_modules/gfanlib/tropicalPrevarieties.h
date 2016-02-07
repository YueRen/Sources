#ifndef GFANLIB_TROPICALPREVARIETY_H
#define GFANLIB_TROPICALPREVARIETY_H

#include <gfanlib/gfanlib.h>
#include <polys/monomials/p_polys.h>
#include <set>

#include <Singular/subexpr.h> // for leftv

class markedCone;
struct markedCone_compare;
typedef std::set<markedCone,markedCone_compare> markedCones;
class tropicalPrevariety;


class markedCone
{
 private:
  gfan::ZCone polyhedralCone;
  gfan::ZVector interiorPoint;

 public:
  markedCone(gfan::ZCone &zc);
  markedCone(const gfan::ZCone &zc, const gfan::ZVector &p);

  gfan::ZCone getPolyhedralCone() const
  {
    return polyhedralCone;
  };

  gfan::ZVector getInteriorPoint() const
  {
    return interiorPoint;
  };

  friend struct markedCone_compare;
};

struct markedCone_compare
{
  bool operator()(const markedCone &sigma, const markedCone &delta) const
  {
    int n = sigma.polyhedralCone.dimension();
    int m = delta.polyhedralCone.dimension();
    if (n!=m)
      return n<m;
    else
      return sigma.interiorPoint<delta.interiorPoint;
  }
};

class tropicalPrevariety
{
 public:
  markedCones maximalCones;
  int ambientDimension;
  int dimension;
  int expectedDimension;

  tropicalPrevariety();
  tropicalPrevariety(int ad, int d, int ed);
  tropicalPrevariety(gfan::ZCone &fundamentalDomain);
  tropicalPrevariety(const poly g, const ring r, const gfan::ZCone* fundamentalDomain);

  void deleteNonMaximalCones();
};

tropicalPrevariety intersectTropicalPrevarieties(const tropicalPrevariety &Sigma, const tropicalPrevariety &Delta);

BOOLEAN gfanlib_tropicalPrevariety(leftv res, leftv args);
BOOLEAN gfanlib_intersectTropicalPrevarieties(leftv res, leftv args);

#endif
