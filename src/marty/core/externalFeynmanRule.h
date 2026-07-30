// This file is part of MARTY.
//
// MARTY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

#pragma once

#include "feynmanRule.h"
#include "model.h"

#include <cstddef>
#include <vector>

namespace mty {

/**
 * @brief Registers a vertex supplied by an external model format.
 *
 * @param model Model receiving the rule.
 * @param fields Ordered external fields. Their points and free indices must
 *        match those used in @p expression.
 * @param expression Vertex expression in MARTY conventions.
 * @param bookkeepingTerm Optional interaction term used by diagram expansion
 *        and filters. When omitted, a scalar bookkeeping term is built from
 *        the supplied fields and the external vertex expression.
 *
 * @details This helper deliberately bypasses the Lagrangian-to-rule
 * derivation. It is therefore appropriate for formats such as UFO, which
 * encode Feynman vertices directly. The interaction term is only used for
 * particle counting, diagram expansion and filters; all coupling and tensor
 * information comes from @p expression.
 */
void AddExternalFeynmanRule(
    mty::Model &model,
    std::vector<mty::QuantumField> const &fields,
    csl::Expr const &expression,
    mty::Lagrangian::TermType bookkeepingTerm = nullptr,
    std::vector<std::size_t> *canonicalOrder = nullptr);

} // namespace mty
