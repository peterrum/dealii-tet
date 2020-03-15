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


// Test PolynomialsTet on quadrature points returned by QGaussTet.


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_tet.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_tet.h>

#include <deal.II/grid/tria.h>

#include "./tests.h"

using namespace dealii;

template <int dim>
void
test()
{
  // 1) Create triangulation (TODO: this should be general mesh - not only for
  // TET)
  TetTriangulation<dim> tria;
  tria.setup(); // TODO: load mesh

  // 2) Create finite element (for TET)
  FE_QTet<dim> fe(1);

  // 3) Create quadrature rule (for TET)
  QGaussTet<dim> quad(3);

  // 4) Create mapping (for TET)
  MappingTet<dim> mapping(1);

  // 5) Create FEValues (for a single set of FiniteElement, Quadrature, Mapping)
  FEValues<dim> fe_values(mapping,
                          fe,
                          quad,
                          update_quadrature_points | update_JxW_values |
                            update_contravariant_transformation |
                            update_covariant_transformation | update_gradients);

  // 6) Loop over all cells of triangulation
  for (auto &cell : tria.cell_iterators())
    {
      // 6a) get different cell quantities
      deallog << "level              = " << cell->level() << " " << std::endl;
      deallog << "index              = " << cell->index() << " " << std::endl;
      deallog << "id                 = " << cell->id() << " " << std::endl;
      deallog << "active             = " << cell->active() << " " << std::endl;
      deallog << "is_locally_owned   = " << cell->is_locally_owned() << " "
              << std::endl;
      deallog << "is_ghost           = " << cell->is_ghost() << " "
              << std::endl;
      deallog << "material_id        = " << cell->material_id() << " "
              << std::endl;
      deallog << "manifold_id        = " << cell->manifold_id() << " "
              << std::endl;
      deallog << "subdomain_id       = " << cell->subdomain_id() << " "
              << std::endl;
      deallog << "level_subdomain_id = " << cell->level_subdomain_id() << " "
              << std::endl;

      deallog << "vertices:" << std::endl;
      for (unsigned int i = 0; i < cell->n_vertices(); i++)
        deallog << "  " << i << " " << cell->vertex_index(i) << " "
                << cell->vertex(i) << std::endl;
      deallog << std::endl;

      // 6b) test FEValues::reinit
      fe_values.reinit(cell);

      // 6c) get quadrature points
      for (unsigned int i = 0; i < quad.size(); i++)
        deallog << fe_values.quadrature_point(i) << std::endl;

      deallog << std::endl;

      // 6d) compute element-stiffness matrix
      {
        const unsigned int dofs_per_cell = 3; // TODO: fe.dofs_per_cell
        const unsigned int n_q_points    = quad.size();
        FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

        for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                 fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                 fe_values.JxW(q_index));           // dx

        cell_matrix.print_formatted(deallog.get_file_stream(), 3, false, 10);
      }
    }
}

int
main()
{
  initlog();

  {
    deallog.push("2d-1");
    test<2>();
    deallog.pop();
  }
}
