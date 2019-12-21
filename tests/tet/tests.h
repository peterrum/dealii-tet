// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include "../tests.h"

using namespace dealii;

#include <deal.II/fe/fe.h>
#include <deal.II/base/quadrature.h>

namespace
{
  template <int dim, int spacedim = dim>
  class FEValues_
  {
  public:
    FEValues_(const FiniteElement<dim, spacedim> &fe,
              const Quadrature<dim> &             quad)
      : fe(fe)
      , quad(quad)
    {}

    double
    JxW(const unsigned int quadrature_point) const
    {
      return quad.weight(quadrature_point);
    }

    Tensor<1, spacedim>
    shape_grad(const unsigned int function_no,
               const unsigned int quadrature_point) const
    {
      return fe.shape_grad(function_no, quad.point(quadrature_point));
    }

  private:
    const FiniteElement<dim, spacedim> &fe;
    const Quadrature<dim> &             quad;
  };

} // namespace