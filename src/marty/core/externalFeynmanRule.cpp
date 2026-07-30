// This file is part of MARTY.
//
// MARTY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

#include "externalFeynmanRule.h"

#include "diracology.h"
#include "graph.h"
#include "interactionTerm.h"

#include <algorithm>
#include <numeric>
#include <tuple>
#include <utility>

namespace mty {
namespace {


bool containsMajoranaFermion(
    std::vector<mty::QuantumField> const &fields)
{
    return std::any_of(fields.begin(), fields.end(), [](auto const &field) {
        return field.isFermionic() && field.isSelfConjugate();
    });
}

bool isDiracTrilinearFermionRule(
    std::vector<mty::QuantumField> const &fields)
{
    if (fields.size() != 3)
        return false;

    std::size_t fermionCount = 0;
    for (auto const &field : fields) {
        if (!field.isFermionic())
            continue;
        if (field.isSelfConjugate())
            return false;
        ++fermionCount;
    }
    return fermionCount == 2;
}

bool sameDiracPhysicalOrientation(
    std::vector<mty::QuantumField> const &storedFields,
    std::vector<mty::QuantumField> const &physicalFields)
{
    if (storedFields.size() != physicalFields.size())
        return false;

    std::vector<bool> used(physicalFields.size(), false);
    for (auto const &stored : storedFields) {
        auto physical = stored.getConjugatedField();
        std::size_t match = physicalFields.size();
        for (std::size_t i = 0; i < physicalFields.size(); ++i) {
            if (used[i])
                continue;
            if (physical.getQuantumParent()
                != physicalFields[i].getQuantumParent())
                continue;
            if (!physical.isSelfConjugate()
                && physical.isComplexConjugate()
                       != physicalFields[i].isComplexConjugate())
                continue;
            match = i;
            break;
        }
        if (match == physicalFields.size())
            return false;
        used[match] = true;
    }
    return true;
}

std::vector<std::size_t> mapDiracPhysicalFieldOrder(
    std::vector<mty::QuantumField> const &physicalFields,
    std::vector<mty::QuantumField> const &storedFields)
{
    std::vector<std::size_t> identity(physicalFields.size());
    std::iota(identity.begin(), identity.end(), std::size_t{0});
    if (storedFields.size() != physicalFields.size())
        return identity;

    std::vector<std::size_t> order;
    order.reserve(physicalFields.size());
    std::vector<bool> used(physicalFields.size(), false);

    for (auto const &stored : storedFields) {
        auto physical = stored.getConjugatedField();
        std::size_t match = physicalFields.size();
        for (std::size_t i = 0; i < physicalFields.size(); ++i) {
            if (used[i])
                continue;
            if (physical.getQuantumParent()
                != physicalFields[i].getQuantumParent())
                continue;
            if (!physical.isSelfConjugate()
                && physical.isComplexConjugate()
                       != physicalFields[i].isComplexConjugate())
                continue;
            match = i;
            break;
        }
        if (match == physicalFields.size())
            return identity;
        used[match] = true;
        order.push_back(match);
    }

    return order;
}

std::vector<std::size_t>
mapCanonicalFieldOrder(
    std::vector<mty::QuantumField> const &fields,
    std::vector<mty::QuantumField> const &nativeFields)
{
    std::vector<std::size_t> identity(fields.size());
    std::iota(identity.begin(), identity.end(), std::size_t{0});
    if (nativeFields.size() != fields.size())
        return identity;

    std::vector<std::size_t> order;
    order.reserve(fields.size());
    std::vector<bool> used(fields.size(), false);

    for (auto const &nativeStored : nativeFields) {
        std::size_t match = fields.size();

        // Conjugation never changes the QuantumFieldParent. Match the particle
        // identity first; using getConjugatedField() here is fragile for
        // Majorana fermions because MARTY may flip their fermion-flow flag
        // without changing the physical particle.
        for (std::size_t i = 0; i < fields.size(); ++i) {
            if (used[i])
                continue;
            if (fields[i].getQuantumParent()
                == nativeStored.getQuantumParent()) {
                match = i;
                break;
            }
        }
        if (match == fields.size())
            return identity;

        used[match] = true;
        order.push_back(match);
    }

    return order;
}

csl::Expr remapExpressionToNativeFields(
    csl::Expr const &expression,
    std::vector<mty::QuantumField> const &fields,
    std::vector<mty::QuantumField> const &nativeFields,
    std::vector<std::size_t> const &order)
{
    if (nativeFields.size() != fields.size() || order.size() != fields.size())
        return csl::DeepRefreshed(expression);

    csl::Expr remapped = csl::DeepRefreshed(expression);

    // Use a two-stage replacement so that an index that happens to be a
    // destination for one leg can never be consumed prematurely while
    // remapping another leg.
    struct IndexMap {
        csl::Index source;
        csl::Index temporary;
        csl::Index destination;
    };
    std::vector<IndexMap> mappings;

    for (std::size_t nativePos = 0; nativePos < nativeFields.size();
         ++nativePos) {
        const auto sourcePos = order[nativePos];
        if (sourcePos >= fields.size())
            return csl::DeepRefreshed(expression);

        auto const &sourceIndices = fields[sourcePos].getIndexStructureView();
        auto const &nativeIndices
            = nativeFields[nativePos].getIndexStructureView();
        if (sourceIndices.size() != nativeIndices.size())
            return csl::DeepRefreshed(expression);

        std::vector<bool> usedNative(nativeIndices.size(), false);
        for (auto const &sourceIndex : sourceIndices) {
            std::size_t nativeIndex = nativeIndices.size();
            for (std::size_t j = 0; j < nativeIndices.size(); ++j) {
                if (usedNative[j])
                    continue;
                if (sourceIndex.getSpace()->isIndexCompatibleWith(
                        nativeIndices[j].getSpace())) {
                    nativeIndex = j;
                    break;
                }
            }
            if (nativeIndex == nativeIndices.size())
                return csl::DeepRefreshed(expression);

            usedNative[nativeIndex] = true;
            mappings.push_back({sourceIndex,
                                sourceIndex.rename(),
                                nativeIndices[nativeIndex]});
        }
    }

    for (auto const &mapping : mappings)
        csl::Replace(remapped, mapping.source, mapping.temporary);
    for (auto const &mapping : mappings)
        csl::Replace(remapped, mapping.temporary, mapping.destination);

    return remapped;
}

mty::Lagrangian::TermType makeBookkeepingTerm(
    std::vector<mty::QuantumField> const &fields,
    csl::Expr const &expression)
{
    // InteractionTerm must represent a Lorentz/gauge scalar.  A bare product
    // of fields is not sufficient for indexed fields (e.g. FFV vertices),
    // because its Dirac/Lorentz/colour indices remain free and the
    // InteractionTerm copy constructor subsequently rejects it.
    //
    // The external rule expression already contains exactly the tensors that
    // contract the free field indices.  Build a bookkeeping scalar by
    // multiplying it by copies of the fields with flipped indices, mirroring
    // FeynmanRule::getFieldProduct().  This term is metadata only: the actual
    // vertex used in amplitudes remains `expression` below.
    std::vector<csl::Expr> factors;
    factors.reserve(fields.size() + 1);

    for (auto field : fields) {
        field.setPoint(mty::defaultSpaceTimePoint);
        field.setExternal(false);
        field.setIncoming(false);

        csl::Expr fieldExpr = field.copy();
        for (auto &index : fieldExpr->getIndexStructureView())
            index = index.getFlipped();
        factors.push_back(std::move(fieldExpr));
    }

    factors.push_back(csl::DeepCopy(expression));
    csl::Expr scalar = csl::prod_s(std::move(factors), true);

    HEPAssert(csl::Abbrev::getFreeStructure(scalar).empty(),
              mty::error::IndexError,
              "External Feynman rule does not contract all field indices.");

    return std::make_shared<mty::InteractionTerm>(scalar);
}

} // namespace

void AddExternalFeynmanRule(
    mty::Model &model,
    std::vector<mty::QuantumField> const &fields,
    csl::Expr const &expression,
    mty::Lagrangian::TermType bookkeepingTerm,
    std::vector<std::size_t> *canonicalOrder)
{
    HEPAssert(not fields.empty(),
              mty::error::ValueError,
              "An external Feynman rule must contain at least one field.");
    HEPAssert(expression != CSL_0,
              mty::error::ValueError,
              "An external Feynman rule cannot have a zero expression.");

    const bool majoranaRule = containsMajoranaFermion(fields);
    const bool diracTrilinearRule = isDiracTrilinearFermionRule(fields);
    if (not bookkeepingTerm) {
        // External formats such as UFO encode the Feynman-rule factor i in
        // the vertex expression.  A native MARTY FeynmanRule, however, is
        // derived from a Lagrangian interaction and supplies that factor
        // itself.  Majorana vertices and ordinary trilinear Dirac vertices
        // are reconstructed through native MARTY interaction terms, so remove
        // the external vertex i for those paths.
        bookkeepingTerm = makeBookkeepingTerm(
            fields,
            (majoranaRule || diracTrilinearRule)
                ? (-CSL_I * expression)
                : expression);
    }

    std::vector<std::size_t> order(fields.size());
    std::iota(order.begin(), order.end(), std::size_t{0});

    if (majoranaRule) {
        // Majorana vertices must be canonicalized as Hermitian pairs, not as
        // independent UFO vertices.  In a native MARTY model an interaction
        // T and T^\dagger are both present while FeynmanRule derives either
        // member of the pair.  Deriving each UFO vertex independently can
        // leave one crossing represented with explicit C matrices even though
        // the conjugate crossing is canonical.
        //
        // First see whether a previously imported Majorana pair already
        // contains the physical orientation requested by this UFO vertex.  A
        // UFO Vertex normally aggregates all Lorentz/coupling components for
        // a fixed particle tuple, so the physical external-field orientation
        // is a suitable partner key here.
        auto samePhysicalOrientation = [](std::vector<mty::QuantumField> const &lhs,
                                          std::vector<mty::QuantumField> const &rhs) {
            if (lhs.size() != rhs.size())
                return false;
            std::vector<bool> used(rhs.size(), false);
            for (auto const &stored : lhs) {
                auto physical = stored.getConjugatedField();
                std::size_t match = rhs.size();
                for (std::size_t i = 0; i < rhs.size(); ++i) {
                    if (used[i])
                        continue;
                    if (physical.getQuantumParent() != rhs[i].getQuantumParent())
                        continue;
                    if (!physical.isSelfConjugate()
                        && physical.isComplexConjugate()
                               != rhs[i].isComplexConjugate())
                        continue;
                    match = i;
                    break;
                }
                if (match == rhs.size())
                    return false;
                used[match] = true;
            }
            return true;
        };

        for (auto const &existing : model.getFeynmanRules()) {
            if (!containsMajoranaFermion(existing.getFieldProduct()))
                continue;
            if (!samePhysicalOrientation(existing.getFieldProduct(), fields))
                continue;
            order = mapCanonicalFieldOrder(fields, existing.getFieldProduct());
            if (canonicalOrder)
                *canonicalOrder = order;
            return;
        }

        auto nativeTerm = std::make_shared<mty::InteractionTerm>(
            *bookkeepingTerm);
        auto conjugatedTerms = mty::InteractionTerm::createAndDispatch(
            csl::GetHermitianConjugate(nativeTerm->getTerm(), &mty::dirac4));
        HEPAssert(!conjugatedTerms.empty(),
                  mty::error::RuntimeError,
                  "Cannot build Hermitian-conjugate Majorana bookkeeping "
                  "context.");

        std::vector<mty::Lagrangian::TermType> pairTerms;
        pairTerms.reserve(1 + conjugatedTerms.size());
        pairTerms.push_back(nativeTerm);
        for (auto const &conjugated : conjugatedTerms) {
            if (*conjugated == *nativeTerm)
                continue;
            pairTerms.push_back(
                std::make_shared<mty::InteractionTerm>(*conjugated));
        }

        auto savedInteraction = std::move(model.L.interaction);
        model.L.interaction = pairTerms;

        std::vector<mty::FeynmanRule> pairRules;
        pairRules.reserve(pairTerms.size());
        try {
            for (auto const &term : pairTerms) {
                mty::FeynmanRule candidate(model, term);
                if (candidate.isEmpty() || candidate.isZero())
                    continue;
                candidate.renameIndices();
                bool duplicate = std::any_of(
                    pairRules.begin(), pairRules.end(),
                    [&](auto const &other) { return candidate == other; });
                if (!duplicate)
                    pairRules.push_back(std::move(candidate));
            }
        }
        catch (...) {
            model.L.interaction = std::move(savedInteraction);
            throw;
        }
        model.L.interaction = std::move(savedInteraction);

        HEPAssert(!pairRules.empty(),
                  mty::error::RuntimeError,
                  "MARTY could not canonicalize an external Majorana "
                  "interaction pair.");

        // Put first the rule whose displayed/physical field orientation
        // matches the UFO vertex currently being processed.  This keeps the
        // generated UFO metadata aligned with the external-rule vector; the
        // Hermitian partner follows immediately and will be recognized (and
        // not re-added) when its UFO Vertex is visited later.
        auto matching = std::find_if(
            pairRules.begin(), pairRules.end(), [&](auto const &candidate) {
                return samePhysicalOrientation(candidate.getFieldProduct(), fields);
            });
        if (matching != pairRules.end() && matching != pairRules.begin())
            std::iter_swap(pairRules.begin(), matching);

        order = mapCanonicalFieldOrder(fields, pairRules.front().getFieldProduct());
        if (canonicalOrder)
            *canonicalOrder = order;

        for (auto &candidate : pairRules)
            model.addFeynmanRule(std::move(candidate));
        return;
    }

    if (diracTrilinearRule) {
        // A manually assembled external rule has the correct analytical
        // vertex but not necessarily MARTY's canonical fermion ordering and
        // Wick graph.  For an ordinary trilinear Dirac interaction there are
        // no identical-boson combinatorial factors, so reconstruct the rule
        // through the native Lagrangian-to-Feynman-rule path.  This keeps the
        // UFO coupling and Lorentz structure while making all downstream
        // operations (Dirac traces, squared amplitudes and numerical library
        // generation) identical to a native MARTY model.
        auto nativeTerm = std::make_shared<mty::InteractionTerm>(
            *bookkeepingTerm);
        auto conjugatedTerms = mty::InteractionTerm::createAndDispatch(
            csl::GetHermitianConjugate(nativeTerm->getTerm(), &mty::dirac4));
        HEPAssert(!conjugatedTerms.empty(),
                  mty::error::RuntimeError,
                  "Cannot build Hermitian-conjugate Dirac bookkeeping "
                  "context.");

        std::vector<mty::Lagrangian::TermType> pairTerms;
        pairTerms.reserve(1 + conjugatedTerms.size());
        pairTerms.push_back(nativeTerm);
        for (auto const &conjugated : conjugatedTerms) {
            if (*conjugated == *nativeTerm)
                continue;
            pairTerms.push_back(
                std::make_shared<mty::InteractionTerm>(*conjugated));
        }

        auto savedInteraction = std::move(model.L.interaction);
        model.L.interaction = pairTerms;

        std::vector<mty::FeynmanRule> candidates;
        candidates.reserve(pairTerms.size());
        try {
            for (auto const &term : pairTerms) {
                mty::FeynmanRule candidate(model, term);
                if (candidate.isEmpty() || candidate.isZero())
                    continue;
                candidate.renameIndices();
                bool duplicate = std::any_of(
                    candidates.begin(), candidates.end(),
                    [&](auto const &other) { return candidate == other; });
                if (!duplicate)
                    candidates.push_back(std::move(candidate));
            }
        }
        catch (...) {
            model.L.interaction = std::move(savedInteraction);
            throw;
        }
        model.L.interaction = std::move(savedInteraction);

        auto matching = std::find_if(
            candidates.begin(), candidates.end(), [&](auto const &candidate) {
                return sameDiracPhysicalOrientation(
                    candidate.getFieldProduct(), fields);
            });
        HEPAssert(matching != candidates.end(),
                  mty::error::RuntimeError,
                  "MARTY could not canonicalize an external trilinear "
                  "Dirac rule in the requested physical orientation.");

        order = mapDiracPhysicalFieldOrder(
            fields, matching->getFieldProduct());
        if (canonicalOrder)
            *canonicalOrder = order;

        model.addFeynmanRule(std::move(*matching));
        return;
    }

    if (canonicalOrder)
        *canonicalOrder = order;

    mty::FeynmanRule rule;
    rule.setInteractionTerm(bookkeepingTerm);
    rule.setFieldProduct(fields);
    rule.setExpr(csl::DeepRefreshed(expression));
    rule.setDiagram(std::make_shared<mty::wick::Graph>());
    model.addFeynmanRule(std::move(rule));
}

} // namespace mty
