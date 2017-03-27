#ifndef TROPICAL_VALUATION_H
#define TROPICAL_VALUATION_H

#include <Singular/libsingular.h>
#include <gfanlib/gfanlib.h>

/***
 * computes the valuation of c in cf as puiseux series
 * ASSUMES: - cf transcendental extention
 *          - puiseux variable = first transcendental variable
 **/
gfan::Rational puiseuxValuation(number c, const coeffs cf);

#endif
