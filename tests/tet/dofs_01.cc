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

#include <deal.II/dofs/dof_handler.h>

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
  // 1) Create triangulation (TODO: this should be general mesh - not only for
  // TET)
  std::vector<Point<dim>>         vertices;
  std::vector<Tet::CellData<dim>> cells;

  vertices.push_back(Point<dim>(1, 0));
  vertices.push_back(Point<dim>(0, 1));
  vertices.push_back(Point<dim>(0, 0));
  vertices.push_back(Point<dim>(1, 1));

  Tet::CellData<dim> cell_1;
  cell_1.type     = Tet::CellTypeEnum::tet;
  cell_1.vertices = {0, 1, 2};
  cells.push_back(cell_1);

  Tet::CellData<dim> cell_2;
  cell_2.type     = Tet::CellTypeEnum::tet;
  cell_2.vertices = {1, 2, 3};
  cells.push_back(cell_2);

  Tet::Triangulation<dim> tria;
  tria.create_triangulation_tet(vertices, cells); // TODO: load mesh

  // 2) Create finite element (for TET)
  Tet::FE_Q<dim> fe(1);

  // 3) Create quadrature rule (for TET)
  Tet::QGauss<dim> quad(3);

  // 4) Create mapping (for TET)
  Tet::MappingQ<dim> mapping(1);

  // 5) Create FEValues (for a single set of FiniteElement, Quadrature, Mapping)
  FEValues<dim> fe_values(mapping,
                          fe,
                          quad,
                          update_quadrature_points | update_JxW_values |
                            update_contravariant_transformation |
                            update_covariant_transformation | update_gradients);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // 6) Loop over all cells of triangulation
  for (auto &cell : tria.cell_iterators())
    {
      fe_values.reinit(cell);

      typename DoFHandler<dim>::cell_iterator cell_dh(*cell, &dof_handler);

      std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);
      cell_dh->set_dof_indices(dof_indices);
    }
  for (auto &cell : dof_handler.cell_iterators())
    {
      fe_values.reinit(cell);

      std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);
      cell->get_dof_indices(dof_indices);
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
