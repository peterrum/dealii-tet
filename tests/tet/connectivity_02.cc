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


template <int dim>
void
test(const std::vector<Tet::CellTypeEnum> &cell_types,
     const std::vector<unsigned int> &     cell_vertices)
{
  Tet::Connectivity<dim> connectivity;

  connectivity.build(cell_types, cell_vertices);
  connectivity.print(deallog.get_file_stream());
}

int
main()
{
  initlog();

  test<3>({Tet::CellTypeEnum::tet, Tet::CellTypeEnum::tet},
          {0, 1, 2, 4, 1, 3, 5});

  test<3>({Tet::CellTypeEnum::tet,
           Tet::CellTypeEnum::tet,
           Tet::CellTypeEnum::tet,
           Tet::CellTypeEnum::tet,
           Tet::CellTypeEnum::tet},
          {0, 1, 2, 4, 1, 2, 3, 7, 1, 4, 5, 7, 2, 4, 6, 7, 1, 2, 4, 7});
}
