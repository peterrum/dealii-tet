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

#include <deal.II/distributed/tria_base.h>

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
      Base(dealii::parallel::TriangulationBase<dim, spacedim> &tria_parallel,
           MPI_Comm                                            mpi_communicator)
        : dealii::TriangulationPolicy::Base<dim, spacedim>(tria_parallel)
        , tria_parallel(tria_parallel)
        , mpi_communicator(mpi_communicator)
        , my_subdomain(Utilities::MPI::this_mpi_process(mpi_communicator))
        , n_subdomains(Utilities::MPI::n_mpi_processes(mpi_communicator))
      {}


      virtual bool
      is_multilevel_hierarchy_constructed() const = 0;

    protected:
      dealii::parallel::TriangulationBase<dim, spacedim> &tria_parallel;
      const MPI_Comm                                      mpi_communicator;
      const unsigned int                                  my_subdomain;
      const unsigned int                                  n_subdomains;
    };
  } // namespace TriangulationPolicy
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE


#endif
