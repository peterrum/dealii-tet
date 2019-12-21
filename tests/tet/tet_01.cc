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


// Test PolynomialsTet on quadrature points returned by QGaussTet.


#include <deal.II/base/polynomials_tet.h>
#include <deal.II/base/quadrature_lib.h>

#include "./tests.h"

using namespace dealii;

template <int dim>
void
test(const unsigned int degree, const unsigned int n_points)
{
  PolynomialsTet<dim> poly(degree);
  QGaussTet<dim>      quad(n_points);

  std::vector<double>         values(poly.n());
  std::vector<Tensor<1, dim>> grads;
  std::vector<Tensor<2, dim>> grad_grads;
  std::vector<Tensor<3, dim>> third_derivatives;
  std::vector<Tensor<4, dim>> fourth_derivatives;

  for (unsigned int i = 0; i < n_points; i++)
    {
      poly.evaluate(quad.point(i),
                    values,
                    grads,
                    grad_grads,
                    third_derivatives,
                    fourth_derivatives);

      deallog << quad.point(i) << " ";
      deallog << quad.weight(i) << " ";

      for (auto v : values)
        deallog << v << " ";
      deallog << std::endl;
    }
}

int
main()
{
  initlog();

  {
    deallog.push("2d-1");
    test<2>(1 /*degree*/, 1 /*n_points*/);
    deallog.pop();
  }
  {
    deallog.push("2d-3");
    test<2>(1 /*degree*/, 3 /*n_points*/);
    deallog.pop();
  }
  {
    deallog.push("3d-1");
    test<3>(1 /*degree*/, 1 /*n_points*/);
    deallog.pop();
  }
  {
    deallog.push("3d-1");
    test<3>(1 /*degree*/, 4 /*n_points*/);
    deallog.pop();
  }
}
