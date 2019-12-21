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


// Create a serial triangulation and copy it.


#include <deal.II/base/polynomials_tet.h>

#include "./tests.h"

using namespace dealii;

template <int dim>
void
test(const unsigned int degree = 1)
{
  PolynomialsTet<dim> poly(degree);

  const Point<dim>            unit_point(1.0 / 3.0, 1.0 / 3.0);
  std::vector<double>         values(poly.n());
  std::vector<Tensor<1, dim>> grads;
  std::vector<Tensor<2, dim>> grad_grads;
  std::vector<Tensor<3, dim>> third_derivatives;
  std::vector<Tensor<4, dim>> fourth_derivatives;

  poly.evaluate(unit_point,
                values,
                grads,
                grad_grads,
                third_derivatives,
                fourth_derivatives);

  for (auto v : values)
    deallog << v << " ";
  deallog << std::endl;
}

int
main(int argc, char *argv[])
{
  initlog();

  {
    deallog.push("2d");
    test<2>();
    deallog.pop();
  }
}
