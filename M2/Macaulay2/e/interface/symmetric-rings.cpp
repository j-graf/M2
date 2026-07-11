// Copyright 2026

#include "interface/symmetric-rings.h"

#include "symmetric-rings/raw-interface.hpp"

const Ring *rawSymmetricRing(const Ring *A)
{
  return symmetric_rings::rawSymmetricRing(A);
}

bool rawSymmetricRingsSetHallLittlewoodParameter(const Ring *R,
                                                 const RingElement *t)
{
  return symmetric_rings::rawSymmetricRingsSetHallLittlewoodParameter(R, t);
}

bool rawSymmetricRingsRememberBasis(const Ring *R,
                                    int basisId,
                                    M2_string displaySymbol,
                                    int displayOrder,
                                    bool isMultiplicative)
{
  return symmetric_rings::rawSymmetricRingsRememberBasis(R,
                                                         basisId,
                                                         displaySymbol,
                                                         displayOrder,
                                                         isMultiplicative);
}

const RingElement *rawSymmetricRingsBasisElement(const Ring *R,
                                                int basisId,
                                                M2_string displaySymbol,
                                                int displayOrder,
                                                bool isMultiplicative,
                                                int innerLength,
                                                M2_arrayint index)
{
  return symmetric_rings::rawSymmetricRingsBasisElement(R,
                                                        basisId,
                                                        displaySymbol,
                                                        displayOrder,
                                                        isMultiplicative,
                                                        innerLength,
                                                        index);
}

const RingElement *rawSymmetricRingsSum(const Ring *R,
                                        engine_RawRingElementArray elements)
{
  return symmetric_rings::rawSymmetricRingsSum(R, elements);
}

const RingElement *rawSymmetricRingsPromoteCollected(const Ring *R,
                                                     const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsPromoteCollected(R, f);
}

const RingElement *rawSymmetricRingsLiftCollected(const Ring *R,
                                                  const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsLiftCollected(R, f);
}

const RingElement *rawSymmetricRingsProduct(const Ring *R,
                                            engine_RawRingElementArray elements)
{
  return symmetric_rings::rawSymmetricRingsProduct(R, elements);
}

const RingElement *rawSymmetricRingsJacobiTrudi(const Ring *R,
                                                int basisId,
                                                M2_string displaySymbol,
                                                int displayOrder,
                                                bool isMultiplicative,
                                                M2_arrayint outer,
                                                M2_arrayint inner)
{
  return symmetric_rings::rawSymmetricRingsJacobiTrudi(R,
                                                       basisId,
                                                       displaySymbol,
                                                       displayOrder,
                                                       isMultiplicative,
                                                       outer,
                                                       inner);
}

const RingElement *rawSymmetricRingsToBasis(const RingElement *f,
                                            int powerSumBasisId,
                                            M2_string powerSumDisplaySymbol,
                                            int powerSumDisplayOrder,
                                            bool powerSumIsMultiplicative,
                                            int targetBasisId,
                                            M2_string targetDisplaySymbol,
                                            int targetDisplayOrder,
                                            bool targetIsMultiplicative)
{
  return symmetric_rings::rawSymmetricRingsToBasis(f,
                                                   powerSumBasisId,
                                                   powerSumDisplaySymbol,
                                                   powerSumDisplayOrder,
                                                   powerSumIsMultiplicative,
                                                   targetBasisId,
                                                   targetDisplaySymbol,
                                                   targetDisplayOrder,
                                                   targetIsMultiplicative);
}

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
    bool targetIsMultiplicative)
{
  return symmetric_rings::rawSymmetricRingsProductToBasisDispatch(
      f,
      g,
      powerSumBasisId,
      powerSumDisplaySymbol,
      powerSumDisplayOrder,
      powerSumIsMultiplicative,
      targetBasisId,
      targetDisplaySymbol,
      targetDisplayOrder,
      targetIsMultiplicative);
}

const RingElement *rawSymmetricRingsPlethysm(const RingElement *f,
                                             const RingElement *g,
                                             int powerSumBasisId,
                                             M2_string powerSumDisplaySymbol,
                                             int powerSumDisplayOrder,
                                             bool powerSumIsMultiplicative)
{
  return symmetric_rings::rawSymmetricRingsPlethysm(f,
                                                    g,
                                                    powerSumBasisId,
                                                    powerSumDisplaySymbol,
                                                    powerSumDisplayOrder,
                                                    powerSumIsMultiplicative);
}

const RingElement *rawSymmetricRingsPlethysmToBasis(const RingElement *f,
                                                    const RingElement *g,
                                                    int powerSumBasisId,
                                                    M2_string powerSumDisplaySymbol,
                                                    int powerSumDisplayOrder,
                                                    bool powerSumIsMultiplicative,
                                                    int targetBasisId,
                                                    M2_string targetDisplaySymbol,
                                                    int targetDisplayOrder,
                                                    bool targetIsMultiplicative)
{
  return symmetric_rings::rawSymmetricRingsPlethysmToBasis(f,
                                                           g,
                                                           powerSumBasisId,
                                                           powerSumDisplaySymbol,
                                                           powerSumDisplayOrder,
                                                           powerSumIsMultiplicative,
                                                           targetBasisId,
                                                           targetDisplaySymbol,
                                                           targetDisplayOrder,
                                                           targetIsMultiplicative);
}

int rawSymmetricRingsSingleBasisId(const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsSingleBasisId(f);
}

bool rawSymmetricRingsHasPlethysmProvenance(const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsHasPlethysmProvenance(f);
}

bool rawSymmetricRingsCopyConversionMetadata(const RingElement *source,
                                             const RingElement *target)
{
  return symmetric_rings::rawSymmetricRingsCopyConversionMetadata(source, target);
}

const RingElement *rawSymmetricRingsOmega(const RingElement *f,
                                          M2_arrayint omegaMap,
                                          bool useSomega)
{
  return symmetric_rings::rawSymmetricRingsOmega(f, omegaMap, useSomega);
}

const RingElement *rawSymmetricRingsStraighten(const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsStraighten(f);
}

const RingElement *rawSymmetricRingsHallInnerProduct(const RingElement *f,
                                                    const RingElement *g,
                                                    int innerProductKind,
                                                    M2_arrayint innerProductMap)
{
  return symmetric_rings::rawSymmetricRingsHallInnerProduct(f,
                                                            g,
                                                            innerProductKind,
                                                            innerProductMap);
}

const RingElement *rawSymmetricRingsBasisCoefficient(
    const RingElement *f,
    const RingElement *targetBasisElement)
{
  return symmetric_rings::rawSymmetricRingsBasisCoefficient(
      f, targetBasisElement);
}

int rawSymmetricRingsTermCount(const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsTermCount(f);
}

const RingElement *rawSymmetricRingsTermCoefficient(const RingElement *f, int i)
{
  return symmetric_rings::rawSymmetricRingsTermCoefficient(f, i);
}

M2_arrayint rawSymmetricRingsTermMonomial(const RingElement *f, int i)
{
  return symmetric_rings::rawSymmetricRingsTermMonomial(f, i);
}

M2_string rawSymmetricRingsElementToString(const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsElementToString(f);
}

M2_string rawSymmetricRingsElementToStringLimited(const RingElement *f,
                                                  int maxTerms)
{
  return symmetric_rings::rawSymmetricRingsElementToStringLimited(f, maxTerms);
}

int rawSymmetricRingsElementWeight(const RingElement *f)
{
  return symmetric_rings::rawSymmetricRingsElementWeight(f);
}

// Local Variables:
// indent-tabs-mode: nil
// End:
