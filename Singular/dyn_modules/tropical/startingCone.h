#ifndef TROPICAL_STARTINGCONE_H
#define TROPICAL_STARTINGCONE_H

#include <Singular/libsingular.h>
#include <gfanlib/gfanlib.h>

#include <groebnerCone.h>

tropical::groebnerCone tropicalStartingCone(const ideal I, const ring r, const gfan::ZVector &w,
                                            const std::set<std::vector<int> > &symmetryGroup = std::set<std::vector<int> >());

#endif
