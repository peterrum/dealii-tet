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


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/tria_tet_cell_type.h>

#include "./tests.h"

using namespace dealii;


template <int dim>
void
test_cell_type(const std::shared_ptr<Tet::CellTypeBase<dim>> &cell)
{
  deallog.push(cell->get_name());
  for (unsigned int d = 1; d <= dim; d++)
    {
      deallog << "  " << cell->n_entities(d) << std::endl;
      for (unsigned int e = 0; e < cell->n_entities(d); e++)
        {
          for (const unsigned int v : cell->vertices_of_entity(d, e))
            deallog << v << " ";
          deallog << std::endl;
        }
      deallog << std::endl;
    }
  deallog.pop();
}

int
main()
{
  initlog();

  test_cell_type(Tet::CellTypeFactory<2>::build(Tet::CellTypeEnum::tet));
  test_cell_type(Tet::CellTypeFactory<2>::build(Tet::CellTypeEnum::quad));
}
