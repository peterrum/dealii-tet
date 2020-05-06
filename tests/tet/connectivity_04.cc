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


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/tet/grid_generator.h>

#include "./tests.h"

using namespace dealii;

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  const MPI_Comm                   comm = MPI_COMM_WORLD;

  initlog();

  std::vector<unsigned int> repetitions{1, 1, 1};
  Point<3>                  p1(0, 0, 0);
  Point<3>                  p2(1, 1, 1);

  Tet::Triangulation<3, 3> tria(comm, false);
  Tet::GridGenerator::subdivided_hyper_rectangle(
    tria, repetitions, p1, p2, false);

  const auto &table = tria.get_entity_table();

  tria.get_connectivity().print(deallog.get_file_stream());

  // loop over all cells
  for (unsigned int i = 0; i < table[3][3].ptr.size() - 1; ++i)
    {
      // print lines of cells
      for (unsigned j = table[3][2].ptr[i]; j < table[3][2].ptr[i + 1]; ++j)
        deallog << std::setw(3) << table[3][2].col[j] << " ";
      deallog << std::endl;

      // print lines of cells
      for (unsigned j = table[3][1].ptr[i]; j < table[3][1].ptr[i + 1]; ++j)
        deallog << std::setw(3) << table[3][1].col[j] << " ";
      deallog << std::endl;

      // loop over all faces of cell
      for (unsigned k = table[3][2].ptr[i]; k < table[3][2].ptr[i + 1]; ++k)
        {
          // print lines of faces
          for (unsigned j = table[2][1].ptr[table[3][2].col[k]];
               j < table[2][1].ptr[table[3][2].col[k] + 1];
               ++j)
            deallog << std::setw(3) << table[2][1].col[j] << " ";
          deallog << std::endl;
        }
      deallog << std::endl;
    }
}