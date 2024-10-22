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


// Compute element stiffness matrix with system of FE_QTet, QGaussTet, and
// FEValues.

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/tria.h>

#include <deal.II/tet/fe_q.h>
#include <deal.II/tet/mapping_q.h>
#include <deal.II/tet/quadrature_lib.h>

#include "./tests.h"

using namespace dealii;

template <int dim>
void
test(const unsigned int degree = 1)
{
  // 1) Create mesh
  Tet::Triangulation<dim> tria(MPI_COMM_WORLD);
  create_mesh_0(tria);

  // 2) Create finite element (for TET)
  FESystem<dim> fe(Tet::FE_Q<dim>(degree), dim);

  // 3) Create quadrature rule (for TET)
  Tet::QGauss<dim> quad(3);

  // 4) Create mapping (for TET)
  Tet::MappingQ<dim> mapping(1);

  // 5) Create FEValues (for a single set of FiniteElement, Quadrature, Mapping)
  FEValues<dim> fe_values(mapping,
                          fe,
                          quad,
                          update_quadrature_points | update_JxW_values |
                            update_contravariant_transformation |
                            update_covariant_transformation | update_gradients);

  fe_values.reinit(tria.begin());

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
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
