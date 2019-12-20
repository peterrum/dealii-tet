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


#include <deal.II/fe/fe_q_tet.h>

#include "../tests.h"

using namespace dealii;

template <int dim>
void
test(const unsigned int degree = 1)
{
  FE_QTet<dim> fe(degree);

  const Point<dim> unit_point(1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0));

  std::cout << fe.dofs_per_cell << std::endl;

  for (unsigned int i = 0; i < 3; i++)
    deallog << fe.shape_value(i, unit_point) << " ";
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
