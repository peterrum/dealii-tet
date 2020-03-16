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


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_tet.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_tet.h>

#include <deal.II/grid/tria.h>

#include "./tests.h"

using namespace dealii;

template <int dim>
void
test()
{
  Tet::Triangulation<dim> tria;

  // TODO: create triangulation via GridIn
  {
    std::vector<Point<dim>>         vertices;
    std::vector<Tet::CellData<dim>> cells;

    vertices.push_back(Point<dim>(1, 0));
    vertices.push_back(Point<dim>(0, 1));
    vertices.push_back(Point<dim>(0, 0));
    vertices.push_back(Point<dim>(1, 1));

    Tet::CellData<dim> cell_1;
    cell_1.vertices = {0, 1, 2};
    cells.push_back(cell_1);

    Tet::CellData<dim> cell_2;
    cell_2.vertices = {1, 2, 3};
    cells.push_back(cell_2);

    tria.create_triangulation_tet(vertices, cells);
  }

  // print cell-by-cell the vertices
  for (const auto &cell : tria.cell_iterators())
    {
      for (unsigned int i = 0; i < cell->n_vertices(); i++)
        deallog << "  " << i << " " << cell->vertex_index(i) << " "
                << cell->vertex(i) << std::endl;
      deallog << std::endl;
    }

  // TODO: write file for Paraview via GridOut
}

int
main()
{
  initlog();

  {
    deallog.push("2d");
    test<2>();
    deallog.pop();
  }
}
