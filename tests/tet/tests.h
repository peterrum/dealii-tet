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

#include "../tests.h"

using namespace dealii;

#include <deal.II/base/quadrature.h>

#include <deal.II/fe/fe.h>

namespace
{
  template <int dim, int spacedim>
  void
  create_mesh_0(Tet::Triangulation<dim, spacedim> &tria)
  {
    std::vector<Point<spacedim>>    vertices;
    std::vector<Tet::CellData<dim>> cells;

    vertices.push_back(Point<dim>(1, 0));
    vertices.push_back(Point<dim>(0, 1));
    vertices.push_back(Point<dim>(0, 0));

    Tet::CellData<dim> cell_1;
    cell_1.type     = Tet::CellTypeEnum::tet;
    cell_1.vertices = {0, 1, 2};
    cells.push_back(cell_1);

    tria.create_triangulation_tet(vertices, cells);
  }

  template <int dim, int spacedim>
  void
  create_mesh_1(Tet::Triangulation<dim, spacedim> &tria)
  {
    std::vector<Point<spacedim>>    vertices;
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

    tria.create_triangulation_tet(vertices, cells);
  }

  template <int dim, int spacedim>
  void
  create_mesh_2(Tet::Triangulation<dim, spacedim> &tria)
  {
    std::vector<Point<spacedim>>    vertices;
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
    cell_2.vertices = {1, 2, 3}; // TODO: why different?
    cells.push_back(cell_2);

    tria.create_triangulation_tet(vertices, cells);
  }

} // namespace