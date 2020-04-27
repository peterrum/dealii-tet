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

#ifndef dealii_tria_parallel_policy_h
#define dealii_tria_parallel_policy_h


#include <deal.II/base/config.h>

#include <deal.II/grid/policy.h>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace TriangulationPolicy
  {
    template <int dim, int spacedim>
    class Base : public dealii::TriangulationPolicy::Base<dim, spacedim>
    {
    public:
      /**
       * A structure that contains information about the distributed
       * triangulation.
       */
      struct NumberCache
      {
        /**
         * Number of locally owned active cells of this MPI rank.
         */
        unsigned int n_locally_owned_active_cells;
        /**
         * The total number of active cells (sum of @p
         * n_locally_owned_active_cells).
         */
        types::global_cell_index n_global_active_cells;
        /**
         * The global number of levels computed as the maximum number of levels
         * taken over all MPI ranks, so <tt>n_levels()<=n_global_levels =
         * max(n_levels() on proc i)</tt>.
         */
        unsigned int n_global_levels;
        /**
         * A set containing the subdomain_id (MPI rank) of the owners of the
         * ghost cells on this processor.
         */
        std::set<types::subdomain_id> ghost_owners;
        /**
         * A set containing the MPI ranks of the owners of the level ghost cells
         * on this processor (for all levels).
         */
        std::set<types::subdomain_id> level_ghost_owners;

        NumberCache();
      };

      Base(dealii::Triangulation<dim, spacedim> &tria,
           MPI_Comm                              mpi_communicator);

      const NumberCache &
      get_number_cache() const;

      virtual void
      copy_triangulation(
        const dealii::Triangulation<dim, spacedim> &other_tria) override;

      virtual bool
      is_multilevel_hierarchy_constructed() const;

    protected:
      const MPI_Comm     mpi_communicator;
      const unsigned int my_subdomain;
      const unsigned int n_subdomains;
      NumberCache        number_cache;

      void
      update_number_cache();
    };
  } // namespace TriangulationPolicy
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE


#endif
