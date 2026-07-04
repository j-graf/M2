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

const RingElement *rawSymmetricRingsToSchurViaHRecursive(
    const RingElement *f,
    int hBasisId,
    M2_string hDisplaySymbol,
    int hDisplayOrder,
    bool hIsMultiplicative,
    int schurBasisId,
    M2_string schurDisplaySymbol,
    int schurDisplayOrder)
{
  return symmetric_rings::rawSymmetricRingsToSchurViaHRecursive(f,
                                                                hBasisId,
                                                                hDisplaySymbol,
                                                                hDisplayOrder,
                                                                hIsMultiplicative,
                                                                schurBasisId,
                                                            schurDisplaySymbol,
                                                            schurDisplayOrder);
}

const RingElement *rawSymmetricRingsToSchurFast(
    const RingElement *f,
    int powerSumBasisId,
    M2_string powerSumDisplaySymbol,
    int powerSumDisplayOrder,
    bool powerSumIsMultiplicative,
    int schurBasisId,
    M2_string schurDisplaySymbol,
    int schurDisplayOrder)
{
  return symmetric_rings::rawSymmetricRingsToSchurFast(f,
                                                       powerSumBasisId,
                                                       powerSumDisplaySymbol,
                                                       powerSumDisplayOrder,
                                                       powerSumIsMultiplicative,
                                                       schurBasisId,
                                                       schurDisplaySymbol,
                                                       schurDisplayOrder);
}

const RingElement *rawSymmetricRingsMultiplyToSchurFast(
    const RingElement *f,
    const RingElement *g,
    int powerSumBasisId,
    M2_string powerSumDisplaySymbol,
    int powerSumDisplayOrder,
    bool powerSumIsMultiplicative,
    int schurBasisId,
    M2_string schurDisplaySymbol,
    int schurDisplayOrder)
{
  return symmetric_rings::rawSymmetricRingsMultiplyToSchurFast(
      f,
      g,
      powerSumBasisId,
      powerSumDisplaySymbol,
      powerSumDisplayOrder,
      powerSumIsMultiplicative,
      schurBasisId,
      schurDisplaySymbol,
      schurDisplayOrder);
}

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
    bool targetIsMultiplicative)
{
  return symmetric_rings::rawSymmetricRingsMultiplyToBasisFast(
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
                                                    M2_arrayint innerProductMap)
{
  return symmetric_rings::rawSymmetricRingsHallInnerProduct(f,
                                                            g,
                                                            innerProductMap);
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
