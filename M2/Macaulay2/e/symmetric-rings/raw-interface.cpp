// Copyright 2026

#include "symmetric-rings/raw-interface.hpp"

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"
#include "exceptions.hpp"
#include "ring-elements/ring-element.hpp"

namespace symmetric_rings {

// ============================================================================
// Engine Interface Methods
// ============================================================================

const Ring *rawSymmetricRing(const Ring *A)
{
  try
    {
      return SymmetricEngineRing::create(A);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

bool rawSymmetricRingsSetHallLittlewoodParameter(const Ring *R,
                                                 const RingElement *t)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return false;
      bool ok = S->setHallLittlewoodParameter(t);
      if (!ok) ERROR("expected a Hall-Littlewood parameter in the coefficient ring");
      return ok;
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return false;
    }
}

bool rawSymmetricRingsRememberBasis(const Ring *R,
                                    int basisId,
                                    M2_string displaySymbol,
                                    int displayOrder,
                                    bool isMultiplicative)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return false;
      S->rememberBasisMetadata(basisId,
                               fromM2String(displaySymbol),
                               displayOrder,
                               isMultiplicative);
      return true;
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return false;
    }
}

const RingElement *rawSymmetricRingsBasisElement(const Ring *R,
                                                int basisId,
                                                M2_string displaySymbol,
                                                int displayOrder,
                                                bool isMultiplicative,
                                                int innerLength,
                                                M2_arrayint index)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      ring_elem result = S->basisElement(basisId,
                                         fromM2String(displaySymbol),
                                         displayOrder,
                                         isMultiplicative,
                                         innerLength,
                                         index);
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsSum(const Ring *R,
                                        engine_RawRingElementArray elements)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      return RingElement::make_raw(S, S->batchSum(elements));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsProduct(const Ring *R,
                                            engine_RawRingElementArray elements)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      return RingElement::make_raw(S, S->batchProduct(elements));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsJacobiTrudi(const Ring *R,
                                                int basisId,
                                                M2_string displaySymbol,
                                                int displayOrder,
                                                bool isMultiplicative,
                                                M2_arrayint outer,
                                                M2_arrayint inner)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      ring_elem result = S->jacobiTrudiBasis(basisId,
                                             fromM2String(displaySymbol),
                                             displayOrder,
                                             isMultiplicative,
                                             partitionFromM2Array(outer),
                                             partitionFromM2Array(inner));
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
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
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result = S->toBasis(f->get_value(),
                                    powerSumBasisId,
                                    fromM2String(powerSumDisplaySymbol),
                                    powerSumDisplayOrder,
                                    powerSumIsMultiplicative,
                                    targetBasisId,
                                    fromM2String(targetDisplaySymbol),
                                    targetDisplayOrder,
                                    targetIsMultiplicative);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
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
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result = S->hToSchurViaRecTrans(f->get_value(),
                                                hBasisId,
                                                fromM2String(hDisplaySymbol),
                                                hDisplayOrder,
                                                hIsMultiplicative,
                                                schurBasisId,
                                                fromM2String(schurDisplaySymbol),
                                                schurDisplayOrder);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsPlethysm(const RingElement *f,
                                             const RingElement *g,
                                             int powerSumBasisId,
                                             M2_string powerSumDisplaySymbol,
                                             int powerSumDisplayOrder,
                                             bool powerSumIsMultiplicative)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      if (g->get_ring() != S)
        {
          ERROR("expected elements in the same symmetric ring");
          return nullptr;
        }
      ring_elem result = S->plethysm(f->get_value(),
                                     g->get_value(),
                                     powerSumBasisId,
                                     fromM2String(powerSumDisplaySymbol),
                                     powerSumDisplayOrder,
                                     powerSumIsMultiplicative);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
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
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      if (g->get_ring() != S)
        {
          ERROR("expected elements in the same symmetric ring");
          return nullptr;
        }
      ring_elem result = S->plethysmToBasis(f->get_value(),
                                            g->get_value(),
                                            powerSumBasisId,
                                            fromM2String(powerSumDisplaySymbol),
                                            powerSumDisplayOrder,
                                            powerSumIsMultiplicative,
                                            targetBasisId,
                                            fromM2String(targetDisplaySymbol),
                                            targetDisplayOrder,
                                            targetIsMultiplicative);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

int rawSymmetricRingsSingleBasisId(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return -1;
      return S->uniformBasisId(f->get_value());
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return -1;
    }
}

const RingElement *rawSymmetricRingsOmega(const RingElement *f,
                                          M2_arrayint omegaMap,
                                          bool useSomega)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result = S->omegaInvolution(f->get_value(), omegaMap, useSomega);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsStraighten(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result = S->straighten(f->get_value());
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsHallInnerProduct(const RingElement *f,
                                                    const RingElement *g,
                                                    M2_arrayint innerProductMap)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      if (g->get_ring() != S)
        {
          ERROR("expected elements in the same symmetric ring");
          return nullptr;
        }
      ring_elem result = S->hallInnerProduct(f->get_value(),
                                             g->get_value(),
                                             innerProductMap);
      if (error()) return nullptr;
      return RingElement::make_raw(S->getCoefficientRing(), result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

int rawSymmetricRingsTermCount(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return 0;
      (void) S;
      return static_cast<int>(polyValue(f->get_value())->terms.size());
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return 0;
    }
}

const RingElement *rawSymmetricRingsTermCoefficient(const RingElement *f, int i)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      const auto *poly = polyValue(f->get_value());
      if (i < 0 || static_cast<size_t>(i) >= poly->terms.size())
        {
          ERROR("symmetric-ring term index out of range");
          return nullptr;
        }
      const Ring *A = S->getCoefficientRing();
      return RingElement::make_raw(A, A->copy(poly->terms[i].coeff));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

M2_arrayint rawSymmetricRingsTermMonomial(const RingElement *f, int i)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      (void) S;
      const auto *poly = polyValue(f->get_value());
      if (i < 0 || static_cast<size_t>(i) >= poly->terms.size())
        {
          ERROR("symmetric-ring term index out of range");
          return nullptr;
        }
      const auto& data = poly->terms[i].monomial.data;
      M2_arrayint result = M2_makearrayint(static_cast<int>(data.size()));
      for (size_t j = 0; j < data.size(); ++j)
        result->array[j] = data[j];
      return result;
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

M2_string rawSymmetricRingsElementToString(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      return toM2String(S->elementString(f->get_value()));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

M2_string rawSymmetricRingsElementToStringLimited(const RingElement *f,
                                                  int maxTerms)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      return toM2String(S->elementString(f->get_value(), maxTerms));
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

int rawSymmetricRingsElementWeight(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return -1;
      return S->elementWeight(f->get_value());
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return -1;
    }
}


} // namespace symmetric_rings
