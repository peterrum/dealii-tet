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


#include <deal.II/base/polynomials_tet.h>

DEAL_II_NAMESPACE_OPEN

template <int dim>
PolynomialsTet<dim>::PolynomialsTet(const unsigned int degree)
  : ScalarPolynomialsBase<dim>(1 /*degree*/, 3 /*n_polynomials*/)
{
  AssertDimension(degree, 1);
  AssertDimension(dimension, 2);
}
template <int dim>
double
PolynomialsTet<dim>::compute_value(const unsigned int i,
                                   const Point<dim> & p) const
{
  if (i == 0)
    return p[0];
  else if (i == 1)
    return p[1];
  else if (i == 2)
    return 1.0 - p[0] - p[1];

  AssertIndexRange(i, 3);
  return 0;
}

template <int dim>
Tensor<1, dim>
PolynomialsTet<dim>::compute_grad(const unsigned int i,
                                  const Point<dim> & p) const
{
  (void)p;

  Tensor<1, dim> grad;

  if (i == 0)
    {
      grad[0] = +1.0;
      grad[1] = +0.0;
    }
  else if (i == 1)
    {
      grad[0] = +0.0;
      grad[1] = +1.0;
    }
  else if (i == 2)
    {
      grad[0] = -1.0;
      grad[1] = -1.0;
    }

  return grad;
}



template <int dim>
Tensor<2, dim>
PolynomialsTet<dim>::compute_grad_grad(const unsigned int i,
                                       const Point<dim> & p) const
{
  (void)i;
  (void)p;

  Assert(false, ExcNotImplemented());
  return Tensor<2, dim>();
}

template <int dim>
void
PolynomialsTet<dim>::evaluate(
  const Point<dim> &           unit_point,
  std::vector<double> &        values,
  std::vector<Tensor<1, dim>> &grads,
  std::vector<Tensor<2, dim>> &grad_grads,
  std::vector<Tensor<3, dim>> &third_derivatives,
  std::vector<Tensor<4, dim>> &fourth_derivatives) const
{
  (void)grads;
  (void)grad_grads;
  (void)third_derivatives;
  (void)fourth_derivatives;

  if (values.size() == 0)
    return;

  AssertDimension(values.size(), 3);

  for (unsigned int i = 0; i < 3; i++)
    values[i] = compute_value(i, unit_point);
}

template <int dim>
std::string
PolynomialsTet<dim>::name() const
{
  return "Tet";
}

template <int dim>
std::unique_ptr<ScalarPolynomialsBase<dim>>
PolynomialsTet<dim>::clone() const
{
  Assert(false, ExcNotImplemented());
}


template class PolynomialsTet<1>;
template class PolynomialsTet<2>;
template class PolynomialsTet<3>;

DEAL_II_NAMESPACE_CLOSE