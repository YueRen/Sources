#ifndef TROPICAL_TRAVERSAL_H
#define TROPICAL_TRAVERSAL_H

#include <Singular/libsingular.h>
#include <gfanlib/gfanlib.h>

std::set<tropical::groebnerCone> tropicalTraversal(const tropical::groebnerCone &startingCone, const gfan::ZVector &w,
                                                   const std::set<std::vector<int> > &symmetryGroup = std::set<std::vector<int> >());
tropical::groebnerCone groebnerWalkTraversal(const tropical::groebnerCone &startingCone, const gfan::ZMatrix &P);

#endif
