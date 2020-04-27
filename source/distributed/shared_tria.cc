// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/sparsity_tools.h>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace shared
  {
    template <int dim, int spacedim>
    Policy<dim, spacedim>::Policy(
      dealii::parallel::TriangulationBase<dim, spacedim> &tria,
      MPI_Comm                                            mpi_communicator,
      const bool     allow_artificial_cells,
      const Settings settings)
      : dealii::parallel::TriangulationPolicy::Base<dim, spacedim>(
          tria,
          mpi_communicator)
      , settings(settings)
      , allow_artificial_cells(allow_artificial_cells)
    {
      const auto partition_settings =
        (Settings::partition_zoltan | Settings::partition_metis |
         Settings::partition_zorder | Settings::partition_custom_signal) &
        settings;
      (void)partition_settings;
      Assert((partition_settings == Settings::partition_auto) ||
               (partition_settings == Settings::partition_metis) ||
               (partition_settings == Settings::partition_zoltan) ||
               (partition_settings == Settings::partition_zorder) ||
               (partition_settings == Settings::partition_custom_signal),
             ExcMessage("Settings must contain exactly one type of the active "
                        "cell partitioning scheme."));

      if (settings & Settings::construct_multigrid_hierarchy)
        Assert(allow_artificial_cells,
               ExcMessage("construct_multigrid_hierarchy requires "
                          "allow_artificial_cells to be set to true."));
    }



    template <int dim, int spacedim>
    bool
    Policy<dim, spacedim>::is_multilevel_hierarchy_constructed() const
    {
      return (settings & Settings::construct_multigrid_hierarchy);
    }



    template <int dim, int spacedim>
    void
    Policy<dim, spacedim>::partition()
    {
      const auto n_subdomains =
        Utilities::MPI::n_mpi_processes(this->mpi_communicator);
      const auto my_subdomain =
        Utilities::MPI::this_mpi_process(this->mpi_communicator);

#ifdef DEBUG
      // Check that all meshes are the same (or at least have the same
      // total number of active cells):
      const unsigned int max_active_cells =
        Utilities::MPI::max(this->tria.n_active_cells(),
                            this->mpi_communicator);
      Assert(
        max_active_cells == this->tria.n_active_cells(),
        ExcMessage(
          "A parallel::shared::Triangulation needs to be refined in the same "
          "way on all processors, but the participating processors don't "
          "agree on the number of active cells."));
#endif

      auto partition_settings =
        (Settings::partition_zoltan | Settings::partition_metis |
         Settings::partition_zorder | Settings::partition_custom_signal) &
        settings;
      if (partition_settings == Settings::partition_auto)
#ifdef DEAL_II_TRILINOS_WITH_ZOLTAN
        partition_settings = Settings::partition_zoltan;
#elif defined DEAL_II_WITH_METIS
        partition_settings = Settings::partition_metis;
#else
        partition_settings = Settings::partition_zorder;
#endif

      if (partition_settings == Settings::partition_zoltan)
        {
#ifndef DEAL_II_TRILINOS_WITH_ZOLTAN
          AssertThrow(false,
                      ExcMessage(
                        "Choosing 'partition_zoltan' requires the library "
                        "to be compiled with support for Zoltan! "
                        "Instead, you might use 'partition_auto' to select "
                        "a partitioning algorithm that is supported "
                        "by your current configuration."));
#else
          GridTools::partition_triangulation(
            n_subdomains, this->tria, SparsityTools::Partitioner::zoltan);
#endif
        }
      else if (partition_settings == Settings::partition_metis)
        {
#ifndef DEAL_II_WITH_METIS
          AssertThrow(false,
                      ExcMessage(
                        "Choosing 'partition_metis' requires the library "
                        "to be compiled with support for METIS! "
                        "Instead, you might use 'partition_auto' to select "
                        "a partitioning algorithm that is supported "
                        "by your current configuration."));
#else
          GridTools::partition_triangulation(n_subdomains,
                                             this->tria,
                                             SparsityTools::Partitioner::metis);
#endif
        }
      else if (partition_settings == Settings::partition_zorder)
        {
          GridTools::partition_triangulation_zorder(n_subdomains, this->tria);
        }
      else if (partition_settings == Settings::partition_custom_signal)
        {
          // User partitions mesh manually
        }
      else
        {
          AssertThrow(false, ExcInternalError())
        }

      // do not partition multigrid levels if user is
      // defining a custom partition
      if ((settings &
           Triangulation<dim,
                         spacedim>::Settings::construct_multigrid_hierarchy) &&
          !(settings & Settings::partition_custom_signal))
        dealii::GridTools::partition_multigrid_levels(this->tria);

      true_subdomain_ids_of_cells.resize(this->tria.n_active_cells());

      // loop over all cells and mark artificial:
      auto cell = this->tria.begin_active();
      auto endc = this->tria.end();

      if (allow_artificial_cells)
        {
          // get active halo layer of (ghost) cells
          // parallel::shared::Triangulation<dim>::
          auto predicate = IteratorFilters::SubdomainEqualTo(my_subdomain);

          auto active_halo_layer_vector =
            dealii::GridTools::compute_active_cell_halo_layer(this->tria,
                                                              predicate);
          std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>
            active_halo_layer(active_halo_layer_vector.begin(),
                              active_halo_layer_vector.end());

          for (unsigned int index = 0; cell != endc; cell++, index++)
            {
              // store original/true subdomain ids:
              true_subdomain_ids_of_cells[index] = cell->subdomain_id();

              if (cell->is_locally_owned() == false &&
                  active_halo_layer.find(cell) == active_halo_layer.end())
                cell->set_subdomain_id(numbers::artificial_subdomain_id);
            }

          // loop over all cells in multigrid hierarchy and mark artificial:
          if (settings &
              Triangulation<dim,
                            spacedim>::Settings::construct_multigrid_hierarchy)
            {
              true_level_subdomain_ids_of_cells.resize(this->tria.n_levels());

              auto predicate = IteratorFilters::LocallyOwnedLevelCell();
              for (unsigned int lvl = 0; lvl < this->tria.n_levels(); ++lvl)
                {
                  true_level_subdomain_ids_of_cells[lvl].resize(
                    this->tria.n_cells(lvl));

                  const auto level_halo_layer_vector =
                    dealii::GridTools::compute_cell_halo_layer_on_level(
                      this->tria, predicate, lvl);
                  std::set<typename Triangulation<dim, spacedim>::cell_iterator>
                    level_halo_layer(level_halo_layer_vector.begin(),
                                     level_halo_layer_vector.end());

                  auto cell = this->tria.begin(lvl);
                  auto endc = this->tria.end(lvl);
                  for (unsigned int index = 0; cell != endc; cell++, index++)
                    {
                      // Store true level subdomain IDs before setting
                      // artificial
                      true_level_subdomain_ids_of_cells[lvl][index] =
                        cell->level_subdomain_id();

                      // for active cells, we must have knowledge of level
                      // subdomain ids of all neighbors to our subdomain, not
                      // just neighbors on the same level. if the cells
                      // subdomain id was not set to artitficial above, we will
                      // also keep its level subdomain id since it is either
                      // owned by this processor or in the ghost layer of the
                      // active mesh.
                      if (cell->is_active() &&
                          cell->subdomain_id() !=
                            numbers::artificial_subdomain_id)
                        continue;

                      // we must have knowledge of our parent in the hierarchy
                      if (cell->has_children())
                        {
                          bool keep_cell = false;
                          for (unsigned int c = 0;
                               c < GeometryInfo<dim>::max_children_per_cell;
                               ++c)
                            if (cell->child(c)->level_subdomain_id() ==
                                my_subdomain)
                              {
                                keep_cell = true;
                                break;
                              }
                          if (keep_cell)
                            continue;
                        }

                      // we must have knowledge of our neighbors on the same
                      // level
                      if (!cell->is_locally_owned_on_level() &&
                          level_halo_layer.find(cell) != level_halo_layer.end())
                        continue;

                      // mark all other cells to artificial
                      cell->set_level_subdomain_id(
                        numbers::artificial_subdomain_id);
                    }
                }
            }
        }
      else
        {
          // just store true subdomain ids
          for (unsigned int index = 0; cell != endc; cell++, index++)
            true_subdomain_ids_of_cells[index] = cell->subdomain_id();
        }

#ifdef DEBUG
      {
        // Assert that each cell is owned by a processor
        unsigned int n_my_cells = 0;
        auto         cell       = this->tria.begin_active();
        auto         endc       = this->tria.end();
        for (; cell != endc; ++cell)
          if (cell->is_locally_owned())
            n_my_cells += 1;

        const unsigned int total_cells =
          Utilities::MPI::sum(n_my_cells, this->mpi_communicator);
        Assert(total_cells == this->tria.n_active_cells(),
               ExcMessage("Not all cells are assigned to a processor."));
      }

      // If running with multigrid, assert that each level
      // cell is owned by a processor
      if (settings & Settings::construct_multigrid_hierarchy)
        {
          unsigned int n_my_cells = 0;
          auto         cell       = this->tria.begin();
          auto         endc       = this->tria.end();
          for (; cell != endc; ++cell)
            if (cell->is_locally_owned_on_level())
              n_my_cells += 1;

          const unsigned int total_cells =
            Utilities::MPI::sum(n_my_cells, this->mpi_communicator);
          Assert(total_cells == this->tria.n_cells(),
                 ExcMessage("Not all cells are assigned to a processor."));
        }
#endif
    }



    template <int dim, int spacedim>
    bool
    Policy<dim, spacedim>::with_artificial_cells() const
    {
      return allow_artificial_cells;
    }



    template <int dim, int spacedim>
    void
    Policy<dim, spacedim>::execute_coarsening_and_refinement()
    {
      this->tria
        .dealii::Triangulation<dim,
                               spacedim>::execute_coarsening_and_refinement();
      this->partition();
      this->update_number_cache();
    }



    template <int dim, int spacedim>
    const std::vector<types::subdomain_id> &
    Policy<dim, spacedim>::get_true_level_subdomain_ids_of_cells(
      const unsigned int level) const
    {
      Assert(level < true_level_subdomain_ids_of_cells.size(),
             ExcInternalError());
      Assert(true_level_subdomain_ids_of_cells[level].size() ==
               this->tria.n_cells(level),
             ExcInternalError());
      return true_level_subdomain_ids_of_cells[level];
    }



    template <int dim, int spacedim>
    const std::vector<types::subdomain_id> &
    Policy<dim, spacedim>::get_true_subdomain_ids_of_cells() const
    {
      return true_subdomain_ids_of_cells;
    }



    template <int dim, int spacedim>
    void
    Policy<dim, spacedim>::create_triangulation(
      const std::vector<Point<spacedim>> &vertices,
      const std::vector<CellData<dim>> &  cells,
      const SubCellData &                 subcelldata)
    {
      try
        {
          this->tria.dealii::Triangulation<dim, spacedim>::create_triangulation(
            vertices, cells, subcelldata);
        }
      catch (
        const typename dealii::Triangulation<dim, spacedim>::DistortedCellList
          &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          Assert(false, ExcInternalError());
        }
      this->partition();
      this->update_number_cache();
    }



    template <int dim, int spacedim>
    void
    Policy<dim, spacedim>::create_triangulation(
      const TriangulationDescription::Description<dim, spacedim>
        &construction_data)
    {
      (void)construction_data;

      Assert(false, ExcInternalError());
    }



    template <int dim, int spacedim>
    void
    Policy<dim, spacedim>::copy_triangulation(
      const dealii::Triangulation<dim, spacedim> &other_tria)
    {
      Assert(
        (dynamic_cast<
           const dealii::parallel::distributed::Triangulation<dim, spacedim> *>(
           &other_tria) == nullptr),
        ExcMessage(
          "Cannot use this function on parallel::distributed::Triangulation."));

      this->dealii::parallel::TriangulationPolicy::Base<dim, spacedim>::
        copy_triangulation(other_tria);
      this->partition();
      this->update_number_cache();
    }



    template <int dim, int spacedim>
    Triangulation<dim, spacedim>::Triangulation(
      MPI_Comm mpi_communicator,
      const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
                     smooth_grid,
      const bool     allow_artificial_cells,
      const Settings settings)
      : dealii::parallel::TriangulationBase<dim, spacedim>(
          mpi_communicator,
          std::shared_ptr<TriangulationPolicy::Base<dim, spacedim>>(
            new Policy<dim, spacedim>(*this,
                                      mpi_communicator,
                                      allow_artificial_cells,
                                      settings)),
          smooth_grid,
          false)
    {}



    template <int dim, int spacedim>
    bool
    Triangulation<dim, spacedim>::is_multilevel_hierarchy_constructed() const
    {
      return this->policy->is_multilevel_hierarchy_constructed();
    }



    template <int dim, int spacedim>
    bool
    Triangulation<dim, spacedim>::with_artificial_cells() const
    {
      return static_cast<Policy<dim, spacedim> *>(this->policy.get())
        ->with_artificial_cells();
    }



    template <int dim, int spacedim>
    const std::vector<types::subdomain_id> &
    Triangulation<dim, spacedim>::get_true_subdomain_ids_of_cells() const
    {
      return static_cast<Policy<dim, spacedim> *>(this->policy.get())
        ->get_true_subdomain_ids_of_cells();
    }



    template <int dim, int spacedim>
    const std::vector<types::subdomain_id> &
    Triangulation<dim, spacedim>::get_true_level_subdomain_ids_of_cells(
      const unsigned int level) const
    {
      return static_cast<Policy<dim, spacedim> *>(this->policy.get())
        ->get_true_level_subdomain_ids_of_cells(level);
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::execute_coarsening_and_refinement()
    {
      this->policy->execute_coarsening_and_refinement();
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::create_triangulation(
      const std::vector<Point<spacedim>> &vertices,
      const std::vector<CellData<dim>> &  cells,
      const SubCellData &                 subcelldata)
    {
      this->policy->create_triangulation(vertices, cells, subcelldata);
    }



    template <int dim, int spacedim>
    void
    Triangulation<dim, spacedim>::create_triangulation(
      const TriangulationDescription::Description<dim, spacedim>
        &construction_data)
    {
      this->policy->create_triangulation(construction_data);
    }

  } // namespace shared
} // namespace parallel



/*-------------- Explicit Instantiations -------------------------------*/
#include "shared_tria.inst"

DEAL_II_NAMESPACE_CLOSE
