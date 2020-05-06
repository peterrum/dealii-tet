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


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/tet/data_out.h>
#include <deal.II/tet/fe_q.h>
#include <deal.II/tet/grid_generator.h>
#include <deal.II/tet/mapping_q.h>
#include <deal.II/tet/quadrature_lib.h>

#include "./tests.h"

using namespace dealii;

static unsigned int counter = 0;

template <int dim, int spacedim>
void
test()
{
  const MPI_Comm comm = MPI_COMM_WORLD;

  Tet::Triangulation<dim, spacedim> tria(comm, false);

  Tet::GridGenerator::subdivided_hyper_rectangle(tria,
                                                 std::vector<unsigned int>{1,
                                                                           1,
                                                                           1},
                                                 Point<dim>(0, 0, 0),
                                                 Point<dim>(1, 1, 1),
                                                 false);

  Tet::FE_Q<dim> fe(2);


  DoFHandler<dim, spacedim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  using VectorType = LinearAlgebra::distributed::Vector<double>;
  VectorType solution;

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  solution.reinit(dof_handler.locally_owned_dofs(),
                  locally_relevant_dofs,
                  comm);


  std::ofstream output("solution_tet." +
                       std::to_string(Utilities::MPI::this_mpi_process(comm)) +
                       ".vtk");
  Tet::data_out(dof_handler, solution, "solution", output);
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  initlog();

  test<3, 3>();
}