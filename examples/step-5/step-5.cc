/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, University of Heidelberg, 1999
 */


// @sect3{Include files}

// Again, the first few include files are already known, so we won't comment
// on them:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

// This one is new. We want to read a triangulation from disk, and the class
// which does this is declared in the following file:
#include <deal.II/grid/grid_in.h>

// We will use a circular domain, and the object describing the boundary of it
// comes from this file:
#include <deal.II/grid/manifold_lib.h>

// This is C++ ...
#include <fstream>
#include <iostream>


// Finally, this has been discussed in previous tutorial programs before:
using namespace dealii;

template <int dim>
class Step5
{
public:
  Step5();
  void run();

private:
  void output_results() const;

  Triangulation<dim> triangulation;

};


template <int dim>
Step5<dim>::Step5()
{}


template <int dim>
void Step5<dim>::output_results() const
{
//  DataOut<dim> data_out;
//
//  data_out.attach_dof_handler(dof_handler);
//  data_out.add_data_vector(solution, "solution");
//
//  data_out.build_patches();
//
//  std::ofstream output("solution-" + std::to_string(cycle) + ".vtu");
//  data_out.write_vtu(output);
}




template <int dim>
void Step5<dim>::run()
{
	std::cout << "Hello!" << std::endl;

	Tet::Triangulation<dim> tria;

	  // TODO: create triangulation via GridIn
	  {
	    std::vector<Point<dim>>         vertices;
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
	    cell_2.vertices = {1, 2, 3};
	    cells.push_back(cell_2);

	    tria.create_triangulation_tet(vertices, cells);
	  }

	  // print cell-by-cell the vertices
	  for (const auto &cell : tria.cell_iterators())
	    {
	      for (unsigned int i = 0; i < cell->n_vertices(); i++)
	        std::cout << "  " << i << " " << cell->vertex_index(i) << " "
	                << cell->vertex(i) << std::endl;
	    }

	  // TODO: write file for Paraview via GridOut


//  GridIn<dim> grid_in;
//  grid_in.attach_triangulation(triangulation);
//  std::ifstream input_file("circle-grid.inp");
//
//  Assert(dim == 2, ExcInternalError());
//
//  grid_in.read_ucd(input_file);
//
//  const SphericalManifold<dim> boundary;
//  triangulation.set_all_manifold_ids_on_boundary(0);
//  triangulation.set_manifold(0, boundary);
//
//  for (unsigned int cycle = 0; cycle < 6; ++cycle)
//    {
//      std::cout << "Cycle " << cycle << ':' << std::endl;
//
//      if (cycle != 0)
//        triangulation.refine_global(1);
//
//      // Now that we have a mesh for sure, we write some output and do all the
//      // things that we have already seen in the previous examples.
//      std::cout << "   Number of active cells: "  //
//                << triangulation.n_active_cells() //
//                << std::endl                      //
//                << "   Total number of cells: "   //
//                << triangulation.n_cells()        //
//                << std::endl;
//
//
//    }
}


// @sect3{The <code>main</code> function}

// The main function looks mostly like the one in the previous example, so we
// won't comment on it further:
int main()
{
  Step5<2> laplace_problem_2d;
  laplace_problem_2d.run();
  return 0;
}
