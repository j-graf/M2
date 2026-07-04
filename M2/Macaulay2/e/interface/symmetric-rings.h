// Copyright 2026

#ifndef M2_INTERFACE_SYMMETRIC_RINGS_H_
#define M2_INTERFACE_SYMMETRIC_RINGS_H_

#include "engine-includes.hpp"

#if defined(__cplusplus)
class Ring;
class RingElement;
#else
typedef struct Ring Ring;
typedef struct RingElement RingElement;
#endif

#if defined(__cplusplus)
extern "C" {
#endif

const Ring *rawSymmetricRing(const Ring *A);
bool rawSymmetricRingsSetHallLittlewoodParameter(const Ring *R,
                                                 const RingElement *t);
bool rawSymmetricRingsRememberBasis(const Ring *R,
                                    int basisId,
                                    M2_string displaySymbol,
                                    int displayOrder,
                                    bool isMultiplicative);
const RingElement *rawSymmetricRingsBasisElement(const Ring *R,
                                                int basisId,
                                                M2_string displaySymbol,
                                                int displayOrder,
                                                bool isMultiplicative,
                                                int innerLength,
                                                M2_arrayint index);
const RingElement *rawSymmetricRingsSum(const Ring *R,
                                        engine_RawRingElementArray elements);
const RingElement *rawSymmetricRingsProduct(const Ring *R,
                                            engine_RawRingElementArray elements);
const RingElement *rawSymmetricRingsJacobiTrudi(const Ring *R,
                                                int basisId,
                                                M2_string displaySymbol,
                                                int displayOrder,
                                                bool isMultiplicative,
                                                M2_arrayint outer,
                                                M2_arrayint inner);
const RingElement *rawSymmetricRingsToBasis(const RingElement *f,
                                            int powerSumBasisId,
                                            M2_string powerSumDisplaySymbol,
                                            int powerSumDisplayOrder,
                                            bool powerSumIsMultiplicative,
                                            int targetBasisId,
                                            M2_string targetDisplaySymbol,
                                            int targetDisplayOrder,
                                            bool targetIsMultiplicative);
const RingElement *rawSymmetricRingsToSchurViaHRecursive(const RingElement *f,
                                                         int hBasisId,
                                                         M2_string hDisplaySymbol,
                                                         int hDisplayOrder,
                                                         bool hIsMultiplicative,
                                                         int schurBasisId,
                                                         M2_string schurDisplaySymbol,
                                                         int schurDisplayOrder);
const RingElement *rawSymmetricRingsToSchurFast(const RingElement *f,
                                                int powerSumBasisId,
                                                M2_string powerSumDisplaySymbol,
                                                int powerSumDisplayOrder,
                                                bool powerSumIsMultiplicative,
                                                int schurBasisId,
                                                M2_string schurDisplaySymbol,
                                                int schurDisplayOrder);
const RingElement *rawSymmetricRingsMultiplyToSchurFast(
    const RingElement *f,
    const RingElement *g,
    int powerSumBasisId,
    M2_string powerSumDisplaySymbol,
    int powerSumDisplayOrder,
    bool powerSumIsMultiplicative,
    int schurBasisId,
    M2_string schurDisplaySymbol,
    int schurDisplayOrder);
const RingElement *rawSymmetricRingsMultiplyToBasisFast(
    const RingElement *f,
    const RingElement *g,
    int powerSumBasisId,
    M2_string powerSumDisplaySymbol,
    int powerSumDisplayOrder,
    bool powerSumIsMultiplicative,
    int targetBasisId,
    M2_string targetDisplaySymbol,
    int targetDisplayOrder,
    bool targetIsMultiplicative);
const RingElement *rawSymmetricRingsPlethysm(const RingElement *f,
                                             const RingElement *g,
                                             int powerSumBasisId,
                                             M2_string powerSumDisplaySymbol,
                                             int powerSumDisplayOrder,
                                             bool powerSumIsMultiplicative);
const RingElement *rawSymmetricRingsPlethysmToBasis(const RingElement *f,
                                                    const RingElement *g,
                                                    int powerSumBasisId,
                                                    M2_string powerSumDisplaySymbol,
                                                    int powerSumDisplayOrder,
                                                    bool powerSumIsMultiplicative,
                                                    int targetBasisId,
                                                    M2_string targetDisplaySymbol,
                                                    int targetDisplayOrder,
                                                    bool targetIsMultiplicative);
int rawSymmetricRingsSingleBasisId(const RingElement *f);
const RingElement *rawSymmetricRingsOmega(const RingElement *f,
                                          M2_arrayint omegaMap,
                                          bool useSomega);
const RingElement *rawSymmetricRingsStraighten(const RingElement *f);
const RingElement *rawSymmetricRingsHallInnerProduct(const RingElement *f,
                                                    const RingElement *g,
                                                    M2_arrayint innerProductMap);
int rawSymmetricRingsTermCount(const RingElement *f);
const RingElement *rawSymmetricRingsTermCoefficient(const RingElement *f, int i);
M2_arrayint rawSymmetricRingsTermMonomial(const RingElement *f, int i);
M2_string rawSymmetricRingsElementToString(const RingElement *f);
M2_string rawSymmetricRingsElementToStringLimited(const RingElement *f,
                                                  int maxTerms);
int rawSymmetricRingsElementWeight(const RingElement *f);

#if defined(__cplusplus)
}
#endif

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
