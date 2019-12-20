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

  PolynomialsTet(const unsigned int degree);

  double
  compute_value(const unsigned int i, const Point<dim> &p) const;

  template <int order>
  Tensor<order, dim>
  compute_derivative(const unsigned int i, const Point<dim> &p) const;

  Tensor<1, dim>
  compute_grad(const unsigned int i, const Point<dim> &p) const;

  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const;

  void
  evaluate(const Point<dim> &           unit_point,
           std::vector<double> &        values,
           std::vector<Tensor<1, dim>> &grads,
           std::vector<Tensor<2, dim>> &grad_grads,
           std::vector<Tensor<3, dim>> &third_derivatives,
           std::vector<Tensor<4, dim>> &fourth_derivatives) const override;

  std::string
  name() const override;

  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override;
};

// template functions
template <int dim>
template <int order>
Tensor<order, dim>
PolynomialsTet<dim>::compute_derivative(const unsigned int i,
                                        const Point<dim> & p) const
{
  Tensor<order, dim> der;

  AssertDimension(order, 1);
  const auto grad = compute_grad(i, p);

  for (unsigned int i = 0; i < dim; i++)
    der[i] = grad[i];

  return der;
}


DEAL_II_NAMESPACE_CLOSE

#endif
