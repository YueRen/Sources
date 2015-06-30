#ifndef CALLGFANLIB_PREFAN_H
#define CALLGFANLIB_PREFAN_H

#include <Singular/blackbox.h>
#include <Singular/ipid.h>

#include <set>
#include <gfanlib/gfanlib.h>

extern int markedConeCollection_CMD;

class markedCone
{
private:
  gfan::ZCone polyhedralCone;
  gfan::ZVector interiorPoint;

public:
  friend struct markedCone_compare;

  markedCone(gfan::ZCone& zc);
  bool contains(const gfan::ZVector& p) const;

  gfan::ZCone getPolyhedralCone() const
  {
    return polyhedralCone;
  };
};

struct markedCone_compare
{
  bool operator()(const markedCone &sigma, const markedCone &theta) const
  {
    const gfan::ZVector w = sigma.interiorPoint;
    const gfan::ZVector v = theta.interiorPoint;
    return v < w;
  }
};



/**
 * implementation of a type representing the generating set of polyhedral fans.
 * consisting of a set of marked cones, this type does not offer the plethory of features of gfan::ZFans
 * but is lighter to deal with.
 * this type was created to be used in various fan traversals.
 */

class markedConeCollection
{
private:
  /**
   * collection of maximal cones that will generate the fan of interest
   * each cone comes equipped with a unique interior point to speed up comparisons
   * IMPORTANT: the interior pont is only unique, if the cones form a fan!!!
   */
  std::set<markedCone,markedCone_compare> setOfMarkedCones;

public:
  bool isEmpty() const;
  int size() const;
  lists getListOfCones() const;
  void insert(gfan::ZCone& zc);
  void remove(gfan::ZCone& zc);
  bool hasConeContaining(const gfan::ZVector& p) const;
  markedCone getConeAndDelete();
};


void bbmarkedConeCollection_setup(SModulFunctions* p);

#endif
