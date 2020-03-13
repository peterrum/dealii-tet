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
  TetTriangulation<dim> tria;
  tria.setup();

  MappingTet<dim> mapping(1);
  FE_QTet<dim>    fe(1);
  QGaussTet<dim>  quad(3);

  FEValues<dim> fe_values(mapping, fe, quad, update_quadrature_points);

  for (auto &cell : tria.cell_iterators())
    {
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


      fe_values.reinit(cell);


      for (unsigned int i = 0; i < quad.size(); i++)
        deallog << fe_values.quadrature_point(i) << std::endl;

      deallog << std::endl;
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
