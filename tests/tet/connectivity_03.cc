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

  std::vector<unsigned int> repetitions{5, 5, 5};
  Point<3>                  p1(0, 0, 0);
  Point<3>                  p2(1, 1, 1);

  Triangulation<3, 3> tria_hex;
  GridGenerator::subdivided_hyper_rectangle(
    tria_hex, repetitions, p1, p2, false);

  Tet::Triangulation<3, 3> tria(comm, false);
  Tet::GridGenerator::hex_to_tet_grid(tria_hex, tria);

  auto cell_hex = tria_hex.begin();
  auto cell_tet = tria.begin();

  deallog << numbers::internal_face_boundary_id << std::endl;

  const auto fu = [&](const auto &cell, const auto i) {
    const auto bid = cell->face(i)->boundary_id();

    if (bid == numbers::internal_face_boundary_id)
      deallog << std::setw(3) << "-"
              << " ";
    else
      deallog << std::setw(3) << bid << " ";
  };

  for (; cell_hex != tria_hex.end(); ++cell_hex)
    {
      for (unsigned int i = 0; i < 6; ++i)
        fu(cell_hex, i);
      deallog << std::endl;

      for (unsigned int i = 0; i < 5; ++i, ++cell_tet)
        {
          for (unsigned int i = 0; i < 4; ++i)
            fu(cell_tet, i);
          deallog << std::endl;
        }
      deallog << std::endl << std::endl << std::endl;
    }
}
