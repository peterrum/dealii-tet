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


#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/tet/fe_q.h>
#include <deal.II/tet/grid_generator.h>
#include <deal.II/tet/mapping_q.h>
#include <deal.II/tet/quadrature_lib.h>

#include "./tests.h"

using namespace dealii;

static unsigned int counter = 0;

void
test(const Point<3> &p1,
     const Point<3> &p2,
     const Point<3> &p3,
     const Point<3> &p4)
{
  const MPI_Comm comm = MPI_COMM_WORLD;

  Tet::Triangulation<3, 3> tria(comm, false);

  std::vector<Point<3>>         vertices;
  std::vector<Tet::CellData<3>> cells;

  vertices.push_back(p1);
  vertices.push_back(p2);
  vertices.push_back(p3);
  vertices.push_back(p4);

  Tet::CellData<3> cell_1;
  cell_1.type     = Tet::CellTypeEnum::tet;
  cell_1.vertices = {0, 1, 2, 3};
  cells.push_back(cell_1);

  tria.create_triangulation_tet(vertices, cells);

  const Tet::FE_Q<3> fe(2);

  Tet::QGauss<3>         quad(4);
  const Tet::MappingQ<3> mapping(1);


  FEValues<3, 3> fe_values(
    mapping, fe, quad, update_JxW_values | update_values | update_gradients);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quad.size();

  Vector<double> cell_rhs(dofs_per_cell);

  for (const auto &cell : tria.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      cell_rhs = 0;

      for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                            1.0 *                               // 1.0
                            fe_values.JxW(q_index));            // dx
          }

      for (auto i : cell_rhs)
        deallog << i << " ";
      deallog << std::endl;
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  initlog();

  const std::vector<Point<3>> p{{0, 0, 0},
                                {1, 0, 0},
                                {0, 1, 0},
                                {1, 1, 0},
                                {0, 0, 1},
                                {1, 0, 1},
                                {0, 1, 1},
                                {1, 1, 1}};

  test(p[0], p[1], p[2], p[4]);
  test(p[0], p[1], p[3], p[5]);

  test(p[1], p[3], p[2], p[7]);
  test(p[0], p[3], p[2], p[6]);

  test(p[1], p[4], p[5], p[7]);
  test(p[0], p[4], p[5], p[6]);

  test(p[2], p[4], p[7], p[6]);
  test(p[3], p[5], p[7], p[6]);

  test(p[1], p[2], p[4], p[7]);
  test(p[0], p[3], p[6], p[5]);
}