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


#include <deal.II/tet/tria_cell_type.h>
#include <deal.II/tet/tria_connectivity.h>

#include "./tests.h"

using namespace dealii;

int
main()
{
  initlog();

  const int dim = 2;

  // clang-format off
  const std::vector<Tet::CellTypeEnum>     cell_types{Tet::CellTypeEnum::quad, Tet::CellTypeEnum::quad, Tet::CellTypeEnum::tet};
  const std::vector<unsigned int> cell_vertices{0, 1, 2, 3, 1, 4, 5, 2, 4, 6, 5};
  // clang-format on

  Tet::Connectivity<dim> connectivity;

  connectivity.build(cell_types, cell_vertices);
  connectivity.print(deallog.get_file_stream());
}
