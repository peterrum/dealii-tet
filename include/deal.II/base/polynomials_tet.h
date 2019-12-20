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


#ifndef dealii_polynomials_tet_h
#define dealii_polynomials_tet_h

#include <deal.II/base/config.h>

#include <deal.II/base/scalar_polynomials_base.h>

DEAL_II_NAMESPACE_OPEN


template <int dim>
class PolynomialsTet : public ScalarPolynomialsBase<dim>
{
public:
  static const unsigned int dimension = dim;

  PolynomialsTet(const unsigned int degree)
    : ScalarPolynomialsBase<dim>(1 /*degree*/, 3 /*n_polynomials*/)
  {
    AssertDimension(degree, 1);
    AssertDimension(dimension, 2);
  }

  void
  evaluate(const Point<dim> &           unit_point,
           std::vector<double> &        values,
           std::vector<Tensor<1, dim>> &grads,
           std::vector<Tensor<2, dim>> &grad_grads,
           std::vector<Tensor<3, dim>> &third_derivatives,
           std::vector<Tensor<4, dim>> &fourth_derivatives) const override
  {
    (void)grads;
    (void)grad_grads;
    (void)third_derivatives;
    (void)fourth_derivatives;

    if (values.size() == 0)
      return;

    AssertDimension(values.size(), 3);

    values[0] = unit_point[0];
    values[1] = unit_point[1];
    values[2] = 1.0 - unit_point[0] - unit_point[1];
  }

  std::string
  name() const override
  {
    return "Tet";
  }

  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override
  {}
};


DEAL_II_NAMESPACE_CLOSE

#endif
