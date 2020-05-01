// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2019 by the deal.II authors
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

#ifndef dealii_tet_grid_generator_h
#define dealii_tet_grid_generator_h


#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <deal.II/grid/tria.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace Tet
{
  namespace GridGenerator
  {
    template <int dim, int spacedim>
    void
    subdivided_hyper_rectangle(Tet::Triangulation<dim, spacedim> &tria,
                               const std::vector<unsigned int> &  repetitions,
                               const Point<dim> &                 p1,
                               const Point<dim> &                 p2,
                               const bool colorize = false)
    {
      AssertDimension(dim, spacedim);

      (void)colorize;

      std::vector<Point<spacedim>>    vertices;
      std::vector<Tet::CellData<dim>> cells;

      if (dim == 2)
        {
          Point<dim> dx((p2[0] - p1[0]) / repetitions[0],
                        (p2[1] - p1[1]) / repetitions[1]);

          for (unsigned int j = 0; j <= repetitions[1]; ++j)
            for (unsigned int i = 0; i <= repetitions[0]; ++i)
              vertices.push_back(
                Point<spacedim>(p1[0] + dx[0] * i, p1[1] + dx[1] * j));

          for (unsigned int j = 0; j < repetitions[1]; ++j)
            for (unsigned int i = 0; i < repetitions[0]; ++i)
              {
                std::array<unsigned int, 4> quad{
                  (j + 0) * (repetitions[0] + 1) + i + 0, //
                  (j + 0) * (repetitions[0] + 1) + i + 1, //
                  (j + 1) * (repetitions[0] + 1) + i + 0, //
                  (j + 1) * (repetitions[0] + 1) + i + 1  //
                };                                        //

                {
                  CellData<dim> tri1;
                  tri1.type = Tet::CellTypeEnum::tet;
                  if (true || (i % 2 == 0) == (j % 2 == 0))
                    tri1.vertices = {quad[1], quad[2], quad[0]};
                  else
                    tri1.vertices = {quad[0], quad[3], quad[2]};

                  cells.push_back(tri1);
                }

                {
                  CellData<dim> tri1;
                  tri1.type = Tet::CellTypeEnum::tet;
                  if (true || (i % 2 == 0) == (j % 2 == 0))
                    tri1.vertices = {quad[2], quad[1], quad[3]};
                  else
                    tri1.vertices = {quad[3], quad[0], quad[1]};

                  cells.push_back(tri1);
                }
              }
        }
      else
        {
          Point<dim> dx((p2[0] - p1[0]) / repetitions[0],
                        (p2[1] - p1[1]) / repetitions[1],
                        (p2[2] - p1[2]) / repetitions[1]);

          for (unsigned int k = 0; k <= repetitions[2]; ++k)
            for (unsigned int j = 0; j <= repetitions[1]; ++j)
              for (unsigned int i = 0; i <= repetitions[0]; ++i)
                vertices.push_back(Point<spacedim>(p1[0] + dx[0] * i,
                                                   p1[1] + dx[1] * j,
                                                   p1[2] + dx[2] * k));

          for (unsigned int k = 0; k < repetitions[2]; ++k)
            for (unsigned int j = 0; j < repetitions[1]; ++j)
              for (unsigned int i = 0; i < repetitions[0]; ++i)
                {
                  std::array<unsigned int, 8> quad{
                    (k + 0) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                      (j + 0) * (repetitions[0] + 1) + i + 0, //
                    (k + 0) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                      (j + 0) * (repetitions[0] + 1) + i + 1, //
                    (k + 0) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                      (j + 1) * (repetitions[0] + 1) + i + 0, //
                    (k + 0) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                      (j + 1) * (repetitions[0] + 1) + i + 1, //
                    (k + 1) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                      (j + 0) * (repetitions[0] + 1) + i + 0, //
                    (k + 1) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                      (j + 0) * (repetitions[0] + 1) + i + 1, //
                    (k + 1) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                      (j + 1) * (repetitions[0] + 1) + i + 0, //
                    (k + 1) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                      (j + 1) * (repetitions[0] + 1) + i + 1 //
                  };                                         //


                  {
                    Tet::CellData<3> cell;
                    cell.type     = Tet::CellTypeEnum::tet;
                    cell.vertices = {quad[0], quad[1], quad[2], quad[4]};
                    cells.push_back(cell);
                  }

                  {
                    Tet::CellData<3> cell;
                    cell.type     = Tet::CellTypeEnum::tet;
                    cell.vertices = {quad[1], quad[2], quad[3], quad[7]};
                    cells.push_back(cell);
                  }

                  {
                    Tet::CellData<3> cell;
                    cell.type     = Tet::CellTypeEnum::tet;
                    cell.vertices = {quad[1], quad[4], quad[5], quad[7]};
                    cells.push_back(cell);
                  }

                  {
                    Tet::CellData<3> cell;
                    cell.type     = Tet::CellTypeEnum::tet;
                    cell.vertices = {quad[2], quad[4], quad[6], quad[7]};
                    cells.push_back(cell);
                  }

                  {
                    Tet::CellData<3> cell;
                    cell.type     = Tet::CellTypeEnum::tet;
                    cell.vertices = {quad[1], quad[2], quad[4], quad[7]};
                    cells.push_back(cell);
                  }
                }
        }

      tria.create_triangulation_tet(vertices, cells);
    }

  } // namespace GridGenerator
} // namespace Tet



DEAL_II_NAMESPACE_CLOSE

#endif
