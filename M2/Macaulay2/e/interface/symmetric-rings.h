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
const RingElement *rawSymmetricRingsToBasis(const RingElement *f,
                                            int powerSumBasisId,
                                            M2_string powerSumDisplaySymbol,
                                            int powerSumDisplayOrder,
                                            bool powerSumIsMultiplicative,
                                            int targetBasisId,
                                            M2_string targetDisplaySymbol,
                                            int targetDisplayOrder,
                                            bool targetIsMultiplicative);
M2_string rawSymmetricRingsElementToString(const RingElement *f);
int rawSymmetricRingsElementWeight(const RingElement *f);

#if defined(__cplusplus)
}
#endif

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
