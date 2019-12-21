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


// Create a serial triangulation and copy it.


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q_tet.h>
#include <deal.II/fe/fe_system.h>

#include "./tests.h"

using namespace dealii;

template <int dim>
void
test(const unsigned int degree = 1)
{
  FESystem<dim> fe(FE_QTet<dim>(degree), dim);
  QDuffy        quad(2, 1.0); // TODO: only working for 2D

  FEValues_<dim> fe_values(fe, quad); // TODO: use actual FEValues

  const unsigned int dofs_per_cell = 3 * dim; // TODO: fe.dofs_per_cell
  const unsigned int n_q_points    = quad.size();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

  for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      const unsigned int component_i = fe.system_to_component_index(i).first;
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          const unsigned int component_j =
            fe.system_to_component_index(j).first;
          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            cell_matrix(i, j) +=
              ((fe_values.shape_grad(i, q_point)[component_i] * //
                fe_values.shape_grad(j, q_point)[component_j])  //
               +                                                //
               (fe_values.shape_grad(i, q_point)[component_j] * //
                fe_values.shape_grad(j, q_point)[component_i])  //
               +                                                //
               ((component_i == component_j) ?                  //
                  (fe_values.shape_grad(i, q_point) *           //
                   fe_values.shape_grad(j, q_point)) :          //
                  0)                                            //
               ) *                                              //
              fe_values.JxW(q_point);
        }
    }

  cell_matrix.print_formatted(deallog.get_file_stream(), 3, false, 10);
}

int
main()
{
  initlog();

  {
    deallog.push("2d");
    test<2>();
    deallog.pop();
  }
}
