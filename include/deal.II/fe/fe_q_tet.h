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

#ifndef dealii_fe_q_tet_h
#define dealii_fe_q_tet_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_tet.h>

#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN

template <int dim>
class FE_QTet : public FE_Poly<PolynomialsTet<dim>, dim>
{
public:
  FE_QTet(const unsigned int degree);

  std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

  std::string
  get_name() const override;

private:
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int deg)
  {
    AssertDimension(deg, 1);
    AssertDimension(dim, 2);

    std::vector<unsigned int> dpo(dim + 1, 0U);
    dpo[0] = 1;
    return dpo;
  }
};

DEAL_II_NAMESPACE_CLOSE

#endif
