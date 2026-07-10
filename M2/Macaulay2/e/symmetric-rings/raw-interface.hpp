// Copyright 2026

#ifndef M2_SYMMETRIC_RINGS_RAW_INTERFACE_HPP_
#define M2_SYMMETRIC_RINGS_RAW_INTERFACE_HPP_

#include "interface/symmetric-rings.h"

namespace symmetric_rings {

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
const RingElement *rawSymmetricRingsPromoteCollected(const Ring *R,
                                                     const RingElement *f);
const RingElement *rawSymmetricRingsLiftCollected(const Ring *R,
                                                  const RingElement *f);
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
const RingElement *rawSymmetricRingsProductToBasisDispatch(
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
bool rawSymmetricRingsHasPlethysmProvenance(const RingElement *f);
bool rawSymmetricRingsCopyConversionMetadata(const RingElement *source,
                                             const RingElement *target);
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

} // namespace symmetric_rings

#endif

// Local Variables:
// indent-tabs-mode: nil
// End:
