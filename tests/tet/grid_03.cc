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


// Test GridIn and GridOut for TET meshes.


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/tet/fe_q.h>
#include <deal.II/tet/mapping_q.h>

#include "./tests.h"

using namespace dealii;

template <int dim, int spacedim = dim>
void
test(const MPI_Comm &   comm,
     const std::string  file_name_in,
     const unsigned int degree)
{
  Tet::Triangulation<dim> tria(MPI_COMM_WORLD);

  GridIn<dim, spacedim> grid_in;
  grid_in.attach_triangulation(tria);

  std::ifstream input_file(file_name_in);

  grid_in.read_ucd(input_file);

  unsigned int counter = 0;
  for (const auto &cell : tria.cell_iterators())
    cell->set_subdomain_id(
      std::min(counter++, Utilities::MPI::n_mpi_processes(comm) - 1));

  GridOut       grid_out;
  std::ofstream out(
    "grid." + std::to_string(Utilities::MPI::this_mpi_process(comm)) + ".vtk");
  grid_out.write_vtk(tria, out);

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
  MPILogInitAll                    all;

  const MPI_Comm comm = MPI_COMM_WORLD;

  {
    deallog.push("2d");
    test<2>(comm, SOURCE_DIR "/grid/tri_2elements.inp", 1);
    deallog.pop();
  }
}
