// Copyright 2026

#include "symmetric-rings/raw-interface.hpp"

#include "symmetric-rings/symmetric-engine-ring.hpp"

#include "error.h"
#include "exceptions.hpp"
#include "ring-elements/ring-element.hpp"

namespace symmetric_rings {

namespace {

Partition partitionFromM2Array(M2_arrayint a)
{
    Partition result;
    if (a == nullptr) return result;
    result.reserve(a->len);
    for (int i = 0; i < a->len; ++i) result.push_back(a->array[i]);
    return result;
}

} // namespace

// ============================================================================
// Ring, Basis, And Arithmetic Interface
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
                                    M2_string canonicalBasisKey,
                                    M2_string displaySymbol,
                                    int displayOrder,
                                    bool isMultiplicative)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return false;
      S->rememberBasisMetadata(basisId,
                               fromM2String(canonicalBasisKey),
                               fromM2String(displaySymbol),
                               displayOrder,
                               isMultiplicative);
      if (error()) return false;
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
                                                int innerLength,
                                                M2_arrayint index)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      ring_elem result = S->basisElement(basisId, innerLength, index);
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

const RingElement *rawSymmetricRingsPromoteCollected(const Ring *R,
                                                     const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      ring_elem result;
      if (!S->promoteCollectedExpansion(f, result)) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsLiftCollected(const Ring *R,
                                                  const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      ring_elem result;
      if (!S->liftCollectedExpansion(f, result)) return nullptr;
      return RingElement::make_raw(S, result);
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
                                                M2_arrayint outer,
                                                M2_arrayint inner)
{
  try
    {
      const auto *S = symmetricRingFromRing(R);
      if (error()) return nullptr;
      ring_elem result = S->jacobiTrudiBasis(
          basisId, partitionFromM2Array(outer), partitionFromM2Array(inner));
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

// ============================================================================
// Conversion, Product, And Plethysm Interface
// ============================================================================

const RingElement *rawSymmetricRingsToBasis(const RingElement *f,
                                            int targetBasisId)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      ring_elem result = S->toBasis(f->get_value(), targetBasisId);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

const RingElement *rawSymmetricRingsProductToBasisDispatch(
    const RingElement *f,
    const RingElement *g,
    int targetBasisId)
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
      ring_elem result = S->productToBasisDispatch(
          f->get_value(),
          g->get_value(),
          targetBasisId);
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
                                             const RingElement *g)
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
      ring_elem result = S->plethysm(f->get_value(), g->get_value());
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
                                                    int targetBasisId)
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
      ring_elem result = S->plethysmToBasisDispatch(f->get_value(),
                                                     g->get_value(),
                                                     targetBasisId);
      if (error()) return nullptr;
      return RingElement::make_raw(S, result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

// ============================================================================
// Conversion Metadata Interface
// ============================================================================

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

bool rawSymmetricRingsHasPlethysmProvenance(const RingElement *f)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return false;
      (void) S;
      return hasCombinatorialTag(
          polyValue(f->get_value())->combinatorialTags,
          CombinatorialTag::Plethysm);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return false;
    }
  }

bool rawSymmetricRingsCopyConversionMetadata(const RingElement *source,
                                             const RingElement *target)
{
  try
    {
      const auto *sourceRing = symmetricRingFromElement(source);
      if (error()) return false;
      const auto *targetRing = symmetricRingFromElement(target);
      if (error()) return false;
      (void) sourceRing;
      (void) targetRing;
      mutablePolyValue(target->get_value())->conversionMetadata =
          polyValue(source->get_value())->conversionMetadata;
      mutablePolyValue(target->get_value())->combinatorialTags =
          polyValue(source->get_value())->combinatorialTags;
      return true;
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return false;
    }
}

// ============================================================================
// Omega, Straightening, And Pairing Interface
// ============================================================================

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
                                                    int innerProductKind,
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
                                             innerProductKind,
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

const RingElement *rawSymmetricRingsBasisCoefficient(
    const RingElement *f,
    const RingElement *targetBasisElement)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      if (targetBasisElement->get_ring() != S)
        {
          ERROR("expected elements in the same symmetric ring");
          return nullptr;
        }
      ring_elem result = S->basisCoefficient(
          f->get_value(), targetBasisElement->get_value());
      if (error()) return nullptr;
      return RingElement::make_raw(S->getCoefficientRing(), result);
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

// ============================================================================
// Term And Presentation Introspection
// ============================================================================

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

M2_arrayint rawSymmetricRingsPresentationTermIndices(const RingElement *f,
                                                     int maxTerms)
{
  try
    {
      const auto *S = symmetricRingFromElement(f);
      if (error()) return nullptr;
      std::vector<size_t> order = S->presentationTermOrder(f->get_value());
      // Presentation always sorts the complete view before applying its limit.
      size_t count = order.size();
      if (maxTerms >= 0)
        count = std::min(count, static_cast<size_t>(maxTerms));
      M2_arrayint result = M2_makearrayint(static_cast<int>(count));
      for (size_t i = 0; i < count; ++i)
        result->array[i] = static_cast<int>(order[i]);
      return result;
    }
  catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return nullptr;
    }
}

M2_arrayint rawSymmetricRingsPresentationTermMonomial(const RingElement *f,
                                                      int i)
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
      std::vector<int> data =
          S->presentationMonomialData(poly->terms[i].monomial);
      M2_arrayint result = M2_makearrayint(static_cast<int>(data.size()));
      for (size_t j = 0; j < data.size(); ++j) result->array[j] = data[j];
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
