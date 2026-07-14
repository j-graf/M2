// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_RAW_INTERFACE_HPP_
#define M2_SYMMETRIC_RINGS_RAW_INTERFACE_HPP_

#include "interface/symmetric-rings.h"

namespace symmetric_rings {

// Ring, basis, and arithmetic interface.

const Ring *rawSymmetricRing(const Ring *A);
bool rawSymmetricRingsSetHallLittlewoodParameter(const Ring *R,
                                                 const RingElement *t);
bool rawSymmetricRingsRememberBasis(const Ring *R,
                                    int basisId,
                                    M2_string canonicalBasisKey,
                                    M2_string displaySymbol,
                                    int displayOrder,
                                    bool isMultiplicative);
const RingElement *rawSymmetricRingsBasisElement(const Ring *R,
                                                int basisId,
                                                int innerLength,
                                                M2_arrayint index);
const RingElement *rawSymmetricRingsSum(const Ring *R,
                                        engine_RawRingElementArray elements);
const RingElement *rawSymmetricRingsPromoteCollected(const Ring *R,
                                                     const RingElement *f);
const RingElement *rawSymmetricRingsLiftCollected(const Ring *R,
                                                  const RingElement *f);
const RingElement *rawSymmetricRingsProduct(const Ring *R,
                                            engine_RawRingElementArray elements);
const RingElement *rawSymmetricRingsJacobiTrudi(const Ring *R,
                                                int basisId,
                                                M2_arrayint outer,
                                                M2_arrayint inner);
// Conversion, product, and plethysm interface.
const RingElement *rawSymmetricRingsToBasis(const RingElement *f,
                                            int targetBasisId);
const RingElement *rawSymmetricRingsProductToBasisDispatch(
    const RingElement *f,
    const RingElement *g,
    int targetBasisId);
const RingElement *rawSymmetricRingsPlethysm(const RingElement *f,
                                             const RingElement *g);
const RingElement *rawSymmetricRingsPlethysmToBasis(const RingElement *f,
                                                    const RingElement *g,
                                                    int targetBasisId);
// Conversion metadata interface.
int rawSymmetricRingsSingleBasisId(const RingElement *f);
bool rawSymmetricRingsHasPlethysmProvenance(const RingElement *f);
bool rawSymmetricRingsCopyConversionMetadata(const RingElement *source,
                                             const RingElement *target);
// Omega, straightening, and pairing interface.
const RingElement *rawSymmetricRingsOmega(const RingElement *f,
                                          M2_arrayint omegaMap,
                                          bool useSomega);
const RingElement *rawSymmetricRingsStraighten(const RingElement *f);
const RingElement *rawSymmetricRingsHallInnerProduct(const RingElement *f,
                                                    const RingElement *g,
                                                    int innerProductKind,
                                                    M2_arrayint innerProductMap);
const RingElement *rawSymmetricRingsBasisCoefficient(
    const RingElement *f,
    const RingElement *targetBasisElement);
// Term and presentation introspection.
int rawSymmetricRingsTermCount(const RingElement *f);
const RingElement *rawSymmetricRingsTermCoefficient(const RingElement *f, int i);
M2_arrayint rawSymmetricRingsTermMonomial(const RingElement *f, int i);
M2_arrayint rawSymmetricRingsPresentationTermIndices(const RingElement *f,
                                                     int maxTerms);
M2_arrayint rawSymmetricRingsPresentationTermMonomial(const RingElement *f,
                                                      int i);
M2_string rawSymmetricRingsElementToString(const RingElement *f);
int rawSymmetricRingsElementWeight(const RingElement *f);

} // namespace symmetric_rings

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
