#ifndef TROPICAL_TRAVERSAL_H
#define TROPICAL_TRAVERSAL_H

#include <Singular/libsingular.h>
#include <gfanlib/gfanlib.h>

std::set<tropical::groebnerCone> tropicalTraversal(const tropical::groebnerCone &startingCone, const std::set<std::vector<int> > &symmetryGroup = std::set<std::vector<int> >());

#endif
