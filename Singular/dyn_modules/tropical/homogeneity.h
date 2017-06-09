#ifndef TROPICAL_HOMOGENEITY_H
#define TROPICAL_HOMOGENEITY_H

#include <Singular/libsingular.h>
#include <gfanlib/gfanlib.h>

/**
 * returns a strictly positive vector with respect to whom I is weighted homogeneous
 * returns a vector of size 0 if none such vector exists
 */
gfan::ZVector homogeneityVector(const ideal I, const ring r);

#endif
