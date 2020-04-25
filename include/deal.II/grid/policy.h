// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2019 by the deal.II authors
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

#ifndef dealii_tria_policy_h
#define dealii_tria_policy_h


#include <deal.II/base/config.h>

#include <deal.II/grid/tria.h>


DEAL_II_NAMESPACE_OPEN


namespace TriangulationPolicy
{
  template <int dim, int spacedim>
  class Base
  {
  public:
    using T = dealii::Triangulation<dim, spacedim>;

    Base(Triangulation<dim, spacedim> &tria)
      : tria(tria)
    {}

    virtual void
    create_triangulation(
      const TriangulationDescription::Description<dim, spacedim>
        &construction_data)
    {
      tria.T::create_triangulation(construction_data);
    }

    virtual void
    create_triangulation(const std::vector<Point<spacedim>> &      vertices,
                         const std::vector<dealii::CellData<dim>> &cells,
                         const SubCellData &                       subcelldata)
    {
      tria.T::create_triangulation(vertices, cells, subcelldata);
    }

    virtual void
    copy_triangulation(const dealii::Triangulation<dim, spacedim> &other_tria)
    {
      tria.T::copy_triangulation(other_tria);
    }

    virtual bool
    prepare_coarsening_and_refinement()
    {
      return tria.T::prepare_coarsening_and_refinement();
    }

    virtual void
    execute_coarsening_and_refinement()
    {
      tria.T::execute_coarsening_and_refinement();
    }

    virtual unsigned int
    coarse_cell_id_to_coarse_cell_index(
      const types::coarse_cell_id coarse_cell_id) const
    {
      return tria.T::coarse_cell_id_to_coarse_cell_index(coarse_cell_id);
    }

    virtual types::coarse_cell_id
    coarse_cell_index_to_coarse_cell_id(
      const unsigned int coarse_cell_index) const
    {
      return tria.T::coarse_cell_index_to_coarse_cell_id(coarse_cell_index);
    }

    virtual std::size_t
    memory_consumption() const
    {
      return 0;
    }

    T &
    get_triangulation()
    {
      return tria;
    }

    const T &
    get_triangulation() const
    {
      return tria;
    }

  protected:
    T &tria;
  };
} // namespace TriangulationPolicy


DEAL_II_NAMESPACE_CLOSE


#endif