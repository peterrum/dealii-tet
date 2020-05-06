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
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/tria.h>

#include <deal.II/tet/fe_q.h>

#include "./tests.h"

using namespace dealii;

template <int dim>
void
test(const unsigned int degree)
{
  // 1) Create mesh
  Tet::Triangulation<dim> tria(MPI_COMM_WORLD);
  create_mesh_1(tria);

  // 2) Create finite element (for TET)
  Tet::FE_Q<dim> fe(degree);

  // 3) create DoFHandler
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints constraint_matrix;
  DoFTools::make_zero_boundary_constraints(dof_handler, constraint_matrix);
  constraint_matrix.close();

  constraint_matrix.print(deallog.get_file_stream());
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
