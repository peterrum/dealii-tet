// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


#include <deal.II/grid/tria_tet.h>

#include "./tests.h"

using namespace dealii;

template <int dim>
void
test()
{
  TetTriangulation<dim> tria;

  for (auto cell : tria.cell_iterators())
    {
    }
}

int
main()
{
  initlog();

  {
    deallog.push("2d-1");
    test<2>();
    deallog.pop();
  }
}
