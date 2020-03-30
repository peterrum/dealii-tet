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


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/tet/fe_q.h>
#include <deal.II/tet/mapping_q.h>
#include <deal.II/tet/quadrature_lib.h>

#include "./tests.h"

using namespace dealii;

template <int dim>
void
test(const unsigned int degree)
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
  cell_2.vertices = {1, 0, 3};
  cells.push_back(cell_2);

  Tet::Triangulation<dim> tria;
  tria.create_triangulation_tet(vertices, cells); // TODO: load mesh

  // 2) Create finite element (for TET)
  Tet::FE_Q<dim> fe(degree);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  deallog << "n_dofs: " << dof_handler.n_dofs() << std::endl;

  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);
  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);
  sparsity_pattern.print(deallog.get_file_stream());

  std::ofstream out("sparsity_pattern.svg");
  sparsity_pattern.print_svg(out);
}

int
main()
{
  initlog();

  {
    deallog.push("2d-1");
    test<2>(1);
    deallog.pop();
  }
  {
    deallog.push("2d-2");
    test<2>(2);
    deallog.pop();
  }
}
