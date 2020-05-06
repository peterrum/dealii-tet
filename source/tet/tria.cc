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

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>


DEAL_II_NAMESPACE_OPEN

namespace Tet
{
  template <int dim, int spacedim>
  const std::vector<types::subdomain_id> &
  Triangulation<dim, spacedim>::get_true_subdomain_ids_of_cells() const
  {
    true_subdomain_ids_of_cells.clear();
    true_subdomain_ids_of_cells.resize(
      std::distance(this->begin(), this->end()));

    // loop over all cells and mark artificial:
    auto cell = this->begin_active();
    auto endc = this->end();

    bool allow_artificial_cells = true;

    const unsigned int my_subdomain =
      Utilities::MPI::this_mpi_process(this->get_communicator());

    if (allow_artificial_cells)
      {
        // get active halo layer of (ghost) cells
        // parallel::shared::Triangulation<dim>::
        auto predicate = IteratorFilters::SubdomainEqualTo(my_subdomain);

        const auto active_halo_layer_vector =
          dealii::GridTools::compute_active_cell_halo_layer(
            *dynamic_cast<const dealii::Triangulation<dim, spacedim> *>(this),
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
      }
    else
      {
        for (unsigned int index = 0; cell != endc; cell++, index++)
          true_subdomain_ids_of_cells[index] = cell->subdomain_id();
      }

    return true_subdomain_ids_of_cells;
  }

// explicit instantiations
#include "tria.inst"

} // namespace Tet

DEAL_II_NAMESPACE_CLOSE
