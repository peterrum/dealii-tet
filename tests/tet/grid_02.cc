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

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/tet/fe_q.h>
#include <deal.II/tet/mapping_q.h>

#include "./tests.h"

using namespace dealii;

template <int dim, int spacedim = dim>
void
test(const std::string file_name)
{
  Tet::Triangulation<dim> tria;

  GridIn<dim, spacedim> grid_in;
  grid_in.attach_triangulation(tria);

  std::ifstream input_file(file_name);

  grid_in.read_ucd(input_file);

  // print cell-by-cell the vertices
  for (const auto &cell : tria.cell_iterators())
    {
      for (unsigned int i = 0; i < cell->n_vertices(); i++)
        deallog << "  " << i << " " << cell->vertex_index(i) << " "
                << cell->vertex(i) << std::endl;
      deallog << std::endl;
    }

  GridOut       grid_out;
  std::ofstream out("grid_tri.vtk");
  grid_out.write_vtk(tria, out);
}

int
main()
{
  initlog();

  {
    deallog.push("2d");
    test<2>(SOURCE_DIR "/grid/tri_1element.inp");
    test<2>(SOURCE_DIR "/grid/tri_2elements.inp");
    deallog.pop();
  }
}
