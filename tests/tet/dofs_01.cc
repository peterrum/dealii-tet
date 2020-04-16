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

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/tria.h>

#include <deal.II/tet/fe_q.h>
#include <deal.II/tet/mapping_q.h>
#include <deal.II/tet/quadrature_lib.h>

#include "./tests.h"

using namespace dealii;

template <int dim>
void
test(const unsigned int degree)
{
  // 1) Create mesh
  Tet::Triangulation<dim> tria;
  create_mesh_1(tria);

  // 2) Create finite element (for TET)
  Tet::FE_Q<dim> fe(degree);

  // 3) create DoFHandler
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // 4) Loop over all cells of triangulation and print dofs
  for (const auto &cell : dof_handler.cell_iterators())
    {
      std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);
      cell->get_dof_indices(dof_indices);

      for (auto &dof_index : dof_indices)
        deallog << dof_index << " ";
      deallog << std::endl;
    }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
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
