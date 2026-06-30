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
                                                M2_arrayint index);
const RingElement *rawSymmetricRingsSum(const Ring *R,
                                        engine_RawRingElementArray elements);
const RingElement *rawSymmetricRingsProduct(const Ring *R,
                                            engine_RawRingElementArray elements);
M2_string rawSymmetricRingsElementToString(const RingElement *f);
int rawSymmetricRingsElementWeight(const RingElement *f);

#if defined(__cplusplus)
}
#endif

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
