#ifndef TROPICAL_SINGULARWRAPPERS_H
#define TROPICAL_SINGULARWRAPPERS_H

#include <Singular/libsingular.h>
#include <gfanlib/gfanlib.h>

gfan::QVector callTropicalPointNewton(ideal I, ring r);
gfan::ZMatrix callTropicalLinkNewton(ideal I, ring r);

ideal tropical_kStd_wrapper(ideal I, ring r, tHomog h=testHomog);
ideal tropical_kNF_wrapper(ideal dividend, ring dividendRing, ideal divisor, ring divisorRing);
ideal tropical_kInterRed_wrapper(ideal I, ring r);

#endif
