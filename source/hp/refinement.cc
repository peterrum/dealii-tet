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


#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/grid/grid_refinement.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/refinement.h>

#include <deal.II/lac/vector.h>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  namespace Refinement
  {
    /**
     * Setting p adaptivity flags
     */
    template <int dim, int spacedim>
    void
    full_p_adaptivity(const dealii::DoFHandler<dim, spacedim> &dof_handler)
    {
      std::vector<bool> p_flags(
        dof_handler.get_triangulation().n_active_cells(), true);

      p_adaptivity_from_flags(dof_handler, p_flags);
    }



    template <int dim, int spacedim>
    void
    p_adaptivity_from_flags(
      const dealii::DoFHandler<dim, spacedim> &dof_handler,
      const std::vector<bool> &                p_flags)
    {
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      p_flags.size());

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned() && p_flags[cell->active_cell_index()])
          {
            if (cell->refine_flag_set())
              {
                const unsigned int super_fe_index =
                  dof_handler.get_fe_collection().next_in_hierarchy(
                    cell->active_fe_index());

                // Reject update if already most superordinate element.
                if (super_fe_index != cell->active_fe_index())
                  cell->set_future_fe_index(super_fe_index);
              }
            else if (cell->coarsen_flag_set())
              {
                const unsigned int sub_fe_index =
                  dof_handler.get_fe_collection().previous_in_hierarchy(
                    cell->active_fe_index());

                // Reject update if already least subordinate element.
                if (sub_fe_index != cell->active_fe_index())
                  cell->set_future_fe_index(sub_fe_index);
              }
          }
    }



    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_absolute_threshold(
      const dealii::DoFHandler<dim, spacedim> &dof_handler,
      const Vector<Number> &                   criteria,
      const Number                             p_refine_threshold,
      const Number                             p_coarsen_threshold,
      const ComparisonFunction<typename identity<Number>::type> &compare_refine,
      const ComparisonFunction<typename identity<Number>::type>
        &compare_coarsen)
    {
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      criteria.size());

      std::vector<bool> p_flags(
        dof_handler.get_triangulation().n_active_cells(), false);

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned() &&
            ((cell->refine_flag_set() &&
              compare_refine(criteria[cell->active_cell_index()],
                             p_refine_threshold)) ||
             (cell->coarsen_flag_set() &&
              compare_coarsen(criteria[cell->active_cell_index()],
                              p_coarsen_threshold))))
          p_flags[cell->active_cell_index()] = true;

      p_adaptivity_from_flags(dof_handler, p_flags);
    }



    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_relative_threshold(
      const dealii::DoFHandler<dim, spacedim> &dof_handler,
      const Vector<Number> &                   criteria,
      const double                             p_refine_fraction,
      const double                             p_coarsen_fraction,
      const ComparisonFunction<typename identity<Number>::type> &compare_refine,
      const ComparisonFunction<typename identity<Number>::type>
        &compare_coarsen)
    {
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      criteria.size());
      Assert((p_refine_fraction >= 0) && (p_refine_fraction <= 1),
             dealii::GridRefinement::ExcInvalidParameterValue());
      Assert((p_coarsen_fraction >= 0) && (p_coarsen_fraction <= 1),
             dealii::GridRefinement::ExcInvalidParameterValue());

      // We first have to determine the maximal and minimal values of the
      // criteria of all flagged cells. We start with the minimal and maximal
      // values of all cells, a range within which the minimal and maximal
      // values on cells flagged for refinement must surely lie.
      Number max_criterion_refine =
               *std::min_element(criteria.begin(), criteria.end()),
             min_criterion_refine =
               *std::max_element(criteria.begin(), criteria.end());
      Number max_criterion_coarsen = max_criterion_refine,
             min_criterion_coarsen = min_criterion_refine;

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            if (cell->refine_flag_set())
              {
                max_criterion_refine =
                  std::max(max_criterion_refine,
                           criteria(cell->active_cell_index()));
                min_criterion_refine =
                  std::min(min_criterion_refine,
                           criteria(cell->active_cell_index()));
              }
            else if (cell->coarsen_flag_set())
              {
                max_criterion_coarsen =
                  std::max(max_criterion_coarsen,
                           criteria(cell->active_cell_index()));
                min_criterion_coarsen =
                  std::min(min_criterion_coarsen,
                           criteria(cell->active_cell_index()));
              }
          }

      const parallel::TriangulationBase<dim, spacedim> *parallel_tria =
        dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
          &dof_handler.get_triangulation());
      if (parallel_tria != nullptr &&
          dynamic_cast<const parallel::shared::Triangulation<dim, spacedim> *>(
            &dof_handler.get_triangulation()) == nullptr)
        {
          max_criterion_refine =
            Utilities::MPI::max(max_criterion_refine,
                                parallel_tria->get_communicator());
          min_criterion_refine =
            Utilities::MPI::min(min_criterion_refine,
                                parallel_tria->get_communicator());
          max_criterion_coarsen =
            Utilities::MPI::max(max_criterion_coarsen,
                                parallel_tria->get_communicator());
          min_criterion_coarsen =
            Utilities::MPI::min(min_criterion_coarsen,
                                parallel_tria->get_communicator());
        }

      // Absent any better strategies, we will set the threshold by linear
      // interpolation for both classes of cells individually.
      const Number threshold_refine =
                     min_criterion_refine +
                     p_refine_fraction *
                       (max_criterion_refine - min_criterion_refine),
                   threshold_coarsen =
                     min_criterion_coarsen +
                     p_coarsen_fraction *
                       (max_criterion_coarsen - min_criterion_coarsen);

      p_adaptivity_from_absolute_threshold(dof_handler,
                                           criteria,
                                           threshold_refine,
                                           threshold_coarsen,
                                           compare_refine,
                                           compare_coarsen);
    }



    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_fixed_number(
      const dealii::DoFHandler<dim, spacedim> &dof_handler,
      const Vector<Number> &                   criteria,
      const double                             p_refine_fraction,
      const double                             p_coarsen_fraction,
      const ComparisonFunction<typename identity<Number>::type> &compare_refine,
      const ComparisonFunction<typename identity<Number>::type>
        &compare_coarsen)
    {
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      criteria.size());
      Assert((p_refine_fraction >= 0) && (p_refine_fraction <= 1),
             dealii::GridRefinement::ExcInvalidParameterValue());
      Assert((p_coarsen_fraction >= 0) && (p_coarsen_fraction <= 1),
             dealii::GridRefinement::ExcInvalidParameterValue());

      // ComparisonFunction returning 'true' or 'false' for any set of
      // parameters. These will be used to overwrite user-provided comparison
      // functions whenever no actual comparison is required in the decision
      // process, i.e. when no or all cells will be refined or coarsened.
      const ComparisonFunction<Number> compare_false =
        [](const Number &, const Number &) { return false; };
      const ComparisonFunction<Number> compare_true =
        [](const Number &, const Number &) { return true; };

      // 1.) First extract from the vector of indicators the ones that
      //     correspond to cells that we locally own.
      unsigned int   n_flags_refinement = 0;
      unsigned int   n_flags_coarsening = 0;
      Vector<Number> indicators_refinement(
        dof_handler.get_triangulation().n_active_cells());
      Vector<Number> indicators_coarsening(
        dof_handler.get_triangulation().n_active_cells());
      for (const auto &cell :
           dof_handler.get_triangulation().active_cell_iterators())
        if (!cell->is_artificial() && cell->is_locally_owned())
          {
            if (cell->refine_flag_set())
              indicators_refinement(n_flags_refinement++) =
                criteria(cell->active_cell_index());
            else if (cell->coarsen_flag_set())
              indicators_coarsening(n_flags_coarsening++) =
                criteria(cell->active_cell_index());
          }
      indicators_refinement.grow_or_shrink(n_flags_refinement);
      indicators_coarsening.grow_or_shrink(n_flags_coarsening);

      // 2.) Determine the number of cells for p-refinement and p-coarsening on
      //     basis of the flagged cells.
      //
      // 3.) Find thresholds for p-refinment and p-coarsening on only those
      //     cells flagged for adaptation.
      //
      //     For cases in which no or all cells flagged for refinement and/or
      //     coarsening are subject to p-adaptation, we usually pick thresholds
      //     that apply to all or none of the cells at once. However here, we
      //     do not know which threshold would suffice for this task because the
      //     user could provide any comparison function. Thus if necessary, we
      //     overwrite the user's choice with suitable functions simplying
      //     returning 'true' and 'false' for any cell with reference wrappers.
      //     Thus, no function object copies are stored.
      //
      // 4.) Perform p-adaptation with absolute thresholds.
      Number threshold_refinement      = 0.;
      Number threshold_coarsening      = 0.;
      auto   reference_compare_refine  = std::cref(compare_refine);
      auto   reference_compare_coarsen = std::cref(compare_coarsen);

      const parallel::TriangulationBase<dim, spacedim> *parallel_tria =
        dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
          &dof_handler.get_triangulation());
      if (parallel_tria != nullptr &&
          dynamic_cast<const parallel::shared::Triangulation<dim, spacedim> *>(
            &dof_handler.get_triangulation()) == nullptr)
        {
#ifndef DEAL_II_WITH_P4EST
          Assert(false, ExcInternalError());
#else
          //
          // parallel implementation with distributed memory
          //

          MPI_Comm mpi_communicator = parallel_tria->get_communicator();

          // 2.) Communicate the number of cells scheduled for p-adaptation
          //     globally.
          const unsigned int n_global_flags_refinement =
            Utilities::MPI::sum(n_flags_refinement, mpi_communicator);
          const unsigned int n_global_flags_coarsening =
            Utilities::MPI::sum(n_flags_coarsening, mpi_communicator);

          const unsigned int target_index_refinement =
            static_cast<unsigned int>(
              std::floor(p_refine_fraction * n_global_flags_refinement));
          const unsigned int target_index_coarsening =
            static_cast<unsigned int>(
              std::ceil((1 - p_coarsen_fraction) * n_global_flags_coarsening));

          // 3.) Figure out the global max and min of the criteria. We don't
          //     need it here, but it's a collective communication call.
          const std::pair<Number, Number> global_min_max_refinement =
            dealii::internal::parallel::distributed::GridRefinement::
              compute_global_min_and_max_at_root(indicators_refinement,
                                                 mpi_communicator);

          const std::pair<Number, Number> global_min_max_coarsening =
            dealii::internal::parallel::distributed::GridRefinement::
              compute_global_min_and_max_at_root(indicators_coarsening,
                                                 mpi_communicator);

          // 3.) Compute thresholds if necessary.
          if (target_index_refinement == 0)
            reference_compare_refine = std::cref(compare_false);
          else if (target_index_refinement == n_global_flags_refinement)
            reference_compare_refine = std::cref(compare_true);
          else
            threshold_refinement = dealii::internal::parallel::distributed::
              GridRefinement::RefineAndCoarsenFixedNumber::compute_threshold(
                indicators_refinement,
                global_min_max_refinement,
                target_index_refinement,
                mpi_communicator);

          if (target_index_coarsening == n_global_flags_coarsening)
            reference_compare_coarsen = std::cref(compare_false);
          else if (target_index_coarsening == 0)
            reference_compare_coarsen = std::cref(compare_true);
          else
            threshold_coarsening = dealii::internal::parallel::distributed::
              GridRefinement::RefineAndCoarsenFixedNumber::compute_threshold(
                indicators_coarsening,
                global_min_max_coarsening,
                target_index_coarsening,
                mpi_communicator);
#endif
        }
      else
        {
          //
          // serial implementation (and parallel::shared implementation)
          //

          // 2.) Determine the number of cells scheduled for p-adaptation.
          const unsigned int n_p_refine_cells = static_cast<unsigned int>(
            std::floor(p_refine_fraction * n_flags_refinement));
          const unsigned int n_p_coarsen_cells = static_cast<unsigned int>(
            std::floor(p_coarsen_fraction * n_flags_coarsening));

          // 3.) Compute thresholds if necessary.
          if (n_p_refine_cells == 0)
            reference_compare_refine = std::cref(compare_false);
          else if (n_p_refine_cells == n_flags_refinement)
            reference_compare_refine = std::cref(compare_true);
          else
            {
              std::nth_element(indicators_refinement.begin(),
                               indicators_refinement.begin() +
                                 n_p_refine_cells - 1,
                               indicators_refinement.end(),
                               std::greater<Number>());
              threshold_refinement =
                *(indicators_refinement.begin() + n_p_refine_cells - 1);
            }

          if (n_p_coarsen_cells == 0)
            reference_compare_coarsen = std::cref(compare_false);
          else if (n_p_coarsen_cells == n_flags_coarsening)
            reference_compare_coarsen = std::cref(compare_true);
          else
            {
              std::nth_element(indicators_coarsening.begin(),
                               indicators_coarsening.begin() +
                                 n_p_coarsen_cells - 1,
                               indicators_coarsening.end(),
                               std::less<Number>());
              threshold_coarsening =
                *(indicators_coarsening.begin() + n_p_coarsen_cells - 1);
            }
        }

      // 4.) Finally perform adaptation.
      p_adaptivity_from_absolute_threshold(dof_handler,
                                           criteria,
                                           threshold_refinement,
                                           threshold_coarsening,
                                           std::cref(reference_compare_refine),
                                           std::cref(
                                             reference_compare_coarsen));
    }



    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_regularity(
      const dealii::DoFHandler<dim, spacedim> &dof_handler,
      const Vector<Number> &                   sobolev_indices)
    {
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      sobolev_indices.size());

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            if (cell->refine_flag_set())
              {
                const unsigned int super_fe_index =
                  dof_handler.get_fe_collection().next_in_hierarchy(
                    cell->active_fe_index());

                // Reject update if already most superordinate element.
                if (super_fe_index != cell->active_fe_index())
                  {
                    const unsigned int super_fe_degree =
                      dof_handler.get_fe_collection()[super_fe_index].degree;

                    if (sobolev_indices[cell->active_cell_index()] >
                        super_fe_degree)
                      cell->set_future_fe_index(super_fe_index);
                  }
              }
            else if (cell->coarsen_flag_set())
              {
                const unsigned int sub_fe_index =
                  dof_handler.get_fe_collection().previous_in_hierarchy(
                    cell->active_fe_index());

                // Reject update if already least subordinate element.
                if (sub_fe_index != cell->active_fe_index())
                  {
                    const unsigned int sub_fe_degree =
                      dof_handler.get_fe_collection()[sub_fe_index].degree;

                    if (sobolev_indices[cell->active_cell_index()] <
                        sub_fe_degree)
                      cell->set_future_fe_index(sub_fe_index);
                  }
              }
          }
    }



    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_reference(
      const dealii::DoFHandler<dim, spacedim> &                  dof_handler,
      const Vector<Number> &                                     criteria,
      const Vector<Number> &                                     references,
      const ComparisonFunction<typename identity<Number>::type> &compare_refine,
      const ComparisonFunction<typename identity<Number>::type>
        &compare_coarsen)
    {
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      criteria.size());
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      references.size());

      std::vector<bool> p_flags(
        dof_handler.get_triangulation().n_active_cells(), false);

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned() &&
            ((cell->refine_flag_set() &&
              compare_refine(criteria[cell->active_cell_index()],
                             references[cell->active_cell_index()])) ||
             (cell->coarsen_flag_set() &&
              compare_coarsen(criteria[cell->active_cell_index()],
                              references[cell->active_cell_index()]))))
          p_flags[cell->active_cell_index()] = true;

      p_adaptivity_from_flags(dof_handler, p_flags);
    }



    /**
     * Error prediction
     */
    template <int dim, typename Number, int spacedim>
    void
    predict_error(const dealii::DoFHandler<dim, spacedim> &dof_handler,
                  const Vector<Number> &                   error_indicators,
                  Vector<Number> &                         predicted_errors,
                  const double                             gamma_p,
                  const double                             gamma_h,
                  const double                             gamma_n)
    {
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      error_indicators.size());
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      predicted_errors.size());
      Assert(0 < gamma_p && gamma_p < 1,
             dealii::GridRefinement::ExcInvalidParameterValue());
      Assert(0 < gamma_h, dealii::GridRefinement::ExcInvalidParameterValue());
      Assert(0 < gamma_n, dealii::GridRefinement::ExcInvalidParameterValue());

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            if (cell->future_fe_index_set()) // p adaptation
              {
                Assert(cell->future_fe_index_set() && !cell->refine_flag_set(),
                       ExcMessage(
                         "For error prediction, a cell marked for p-adaptation "
                         "should not also be flagged for h-refinement!"));
                Assert(cell->future_fe_index_set() && !cell->coarsen_flag_set(),
                       ExcMessage(
                         "For error prediction, a cell marked for p-adaptation "
                         "should not also be flagged for h-coarsening!"));

                const int degree_difference =
                  dof_handler.get_fe_collection()[cell->future_fe_index()]
                    .degree -
                  cell->get_fe().degree;

                predicted_errors[cell->active_cell_index()] =
                  error_indicators[cell->active_cell_index()] *
                  std::pow(gamma_p, degree_difference);
              }
            else if (cell->refine_flag_set()) // h refinement
              {
                Assert(
                  cell->refine_flag_set() ==
                    RefinementCase<dim>::isotropic_refinement,
                  ExcMessage(
                    "Error prediction is only valid for isotropic refinement!"));

                predicted_errors[cell->active_cell_index()] =
                  error_indicators[cell->active_cell_index()] *
                  (gamma_h * std::pow(.5, dim + cell->get_fe().degree));
              }
            else if (cell->coarsen_flag_set()) // h coarsening
              {
                predicted_errors[cell->active_cell_index()] =
                  error_indicators[cell->active_cell_index()] /
                  (gamma_h * std::pow(.5, cell->get_fe().degree));
              }
            else // no changes
              {
                predicted_errors[cell->active_cell_index()] =
                  error_indicators[cell->active_cell_index()] * gamma_n;
              }
          }
    }



    /**
     * Decide between h and p adaptivity
     */
    template <int dim, int spacedim>
    void
    force_p_over_h(const dealii::DoFHandler<dim, spacedim> &dof_handler)
    {
      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned() && cell->future_fe_index_set())
          {
            cell->clear_refine_flag();
            cell->clear_coarsen_flag();
          }
    }



    template <int dim, int spacedim>
    void
    choose_p_over_h(const dealii::DoFHandler<dim, spacedim> &dof_handler)
    {
      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned() && cell->future_fe_index_set())
          {
            cell->clear_refine_flag();

            // A cell will only be coarsened into its parent if all of its
            // siblings are flagged for h coarsening as well. We must take this
            // into account for our decision whether we would like to impose h
            // or p adaptivity.
            if (cell->coarsen_flag_set())
              {
                const auto &       parent     = cell->parent();
                const unsigned int n_children = parent->n_children();

                unsigned int h_flagged_children = 0, p_flagged_children = 0;
                for (unsigned int child_index = 0; child_index < n_children;
                     ++child_index)
                  {
                    const auto &child = parent->child(child_index);
                    if (child->is_active())
                      {
                        if (child->coarsen_flag_set())
                          ++h_flagged_children;
                        if (child->future_fe_index_set())
                          ++p_flagged_children;
                      }
                  }

                if (h_flagged_children == n_children &&
                    p_flagged_children != n_children)
                  // Perform pure h coarsening and
                  // drop all p adaptation flags.
                  for (unsigned int child_index = 0; child_index < n_children;
                       ++child_index)
                    parent->child(child_index)->clear_future_fe_index();
                else
                  // Perform p adaptation on all children and
                  // drop all h coarsening flags.
                  for (unsigned int child_index = 0; child_index < n_children;
                       ++child_index)
                    if (parent->child(child_index)->is_active())
                      parent->child(child_index)->clear_coarsen_flag();
              }
          }
    }
  } // namespace Refinement
} // namespace hp


// explicit instantiations
#include "refinement.inst"

DEAL_II_NAMESPACE_CLOSE
