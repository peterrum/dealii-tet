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

#include <deal.II/base/config.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler_policy.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_levels.h>

#include <deal.II/tet/fe_q.h>

#include <algorithm>
#include <set>
#include <unordered_set>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
const typename DoFHandler<dim, spacedim>::active_fe_index_type
  DoFHandler<dim, spacedim>::invalid_active_fe_index;

namespace internal
{
  template <int dim, int spacedim>
  std::string
  policy_to_string(const dealii::internal::DoFHandlerImplementation::Policy::
                     PolicyBase<dim, spacedim> &policy)
  {
    std::string policy_name;
    if (dynamic_cast<const typename dealii::internal::DoFHandlerImplementation::
                       Policy::Sequential<dealii::DoFHandler<dim, spacedim>> *>(
          &policy) ||
        dynamic_cast<const typename dealii::internal::DoFHandlerImplementation::
                       Policy::Sequential<dealii::DoFHandler<dim, spacedim>> *>(
          &policy))
      policy_name = "Policy::Sequential<";
    else if (dynamic_cast<
               const typename dealii::internal::DoFHandlerImplementation::
                 Policy::ParallelDistributed<dealii::DoFHandler<dim, spacedim>>
                   *>(&policy) ||
             dynamic_cast<
               const typename dealii::internal::DoFHandlerImplementation::
                 Policy::ParallelDistributed<dealii::DoFHandler<dim, spacedim>>
                   *>(&policy))
      policy_name = "Policy::ParallelDistributed<";
    else if (dynamic_cast<
               const typename dealii::internal::DoFHandlerImplementation::
                 Policy::ParallelShared<dealii::DoFHandler<dim, spacedim>> *>(
               &policy) ||
             dynamic_cast<
               const typename dealii::internal::DoFHandlerImplementation::
                 Policy::ParallelShared<dealii::DoFHandler<dim, spacedim>> *>(
               &policy))
      policy_name = "Policy::ParallelShared<";
    else
      AssertThrow(false, ExcNotImplemented());
    policy_name += Utilities::int_to_string(dim) + "," +
                   Utilities::int_to_string(spacedim) + ">";
    return policy_name;
  }


  namespace DoFHandlerImplementation
  {
    /**
     * A class with the same purpose as the similarly named class of the
     * Triangulation class. See there for more information.
     */
    struct Implementation
    {
      /**
       * Implement the function of same name in
       * the mother class.
       */
      template <int spacedim>
      static unsigned int
      max_couplings_between_dofs(const DoFHandler<1, spacedim> &dof_handler)
      {
        return std::min(static_cast<types::global_dof_index>(
                          3 * dof_handler.fe_collection.max_dofs_per_vertex() +
                          2 * dof_handler.fe_collection.max_dofs_per_line()),
                        dof_handler.n_dofs());
      }

      template <int spacedim>
      static unsigned int
      max_couplings_between_dofs(const DoFHandler<2, spacedim> &dof_handler)
      {
        // get these numbers by drawing pictures
        // and counting...
        // example:
        //   |     |     |
        // --x-----x--x--X--
        //   |     |  |  |
        //   |     x--x--x
        //   |     |  |  |
        // --x--x--*--x--x--
        //   |  |  |     |
        //   x--x--x     |
        //   |  |  |     |
        // --X--x--x-----x--
        //   |     |     |
        // x = vertices connected with center vertex *;
        //   = total of 19
        // (the X vertices are connected with * if
        // the vertices adjacent to X are hanging
        // nodes)
        // count lines -> 28 (don't forget to count
        // mother and children separately!)
        types::global_dof_index max_couplings;
        switch (dof_handler.tria->max_adjacent_cells())
          {
            case 4:
              max_couplings =
                19 * dof_handler.fe_collection.max_dofs_per_vertex() +
                28 * dof_handler.fe_collection.max_dofs_per_line() +
                8 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 5:
              max_couplings =
                21 * dof_handler.fe_collection.max_dofs_per_vertex() +
                31 * dof_handler.fe_collection.max_dofs_per_line() +
                9 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 6:
              max_couplings =
                28 * dof_handler.fe_collection.max_dofs_per_vertex() +
                42 * dof_handler.fe_collection.max_dofs_per_line() +
                12 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 7:
              max_couplings =
                30 * dof_handler.fe_collection.max_dofs_per_vertex() +
                45 * dof_handler.fe_collection.max_dofs_per_line() +
                13 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 8:
              max_couplings =
                37 * dof_handler.fe_collection.max_dofs_per_vertex() +
                56 * dof_handler.fe_collection.max_dofs_per_line() +
                16 * dof_handler.fe_collection.max_dofs_per_quad();
              break;

            // the following numbers are not based on actual counting but by
            // extrapolating the number sequences from the previous ones (for
            // example, for dofs_per_vertex, the sequence above is 19, 21, 28,
            // 30, 37, and is continued as follows):
            case 9:
              max_couplings =
                39 * dof_handler.fe_collection.max_dofs_per_vertex() +
                59 * dof_handler.fe_collection.max_dofs_per_line() +
                17 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 10:
              max_couplings =
                46 * dof_handler.fe_collection.max_dofs_per_vertex() +
                70 * dof_handler.fe_collection.max_dofs_per_line() +
                20 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 11:
              max_couplings =
                48 * dof_handler.fe_collection.max_dofs_per_vertex() +
                73 * dof_handler.fe_collection.max_dofs_per_line() +
                21 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 12:
              max_couplings =
                55 * dof_handler.fe_collection.max_dofs_per_vertex() +
                84 * dof_handler.fe_collection.max_dofs_per_line() +
                24 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 13:
              max_couplings =
                57 * dof_handler.fe_collection.max_dofs_per_vertex() +
                87 * dof_handler.fe_collection.max_dofs_per_line() +
                25 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 14:
              max_couplings =
                63 * dof_handler.fe_collection.max_dofs_per_vertex() +
                98 * dof_handler.fe_collection.max_dofs_per_line() +
                28 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 15:
              max_couplings =
                65 * dof_handler.fe_collection.max_dofs_per_vertex() +
                103 * dof_handler.fe_collection.max_dofs_per_line() +
                29 * dof_handler.fe_collection.max_dofs_per_quad();
              break;
            case 16:
              max_couplings =
                72 * dof_handler.fe_collection.max_dofs_per_vertex() +
                114 * dof_handler.fe_collection.max_dofs_per_line() +
                32 * dof_handler.fe_collection.max_dofs_per_quad();
              break;

            default:
              Assert(false, ExcNotImplemented());
              max_couplings = 0;
          }
        return std::min(max_couplings, dof_handler.n_dofs());
      }

      template <int spacedim>
      static unsigned int
      max_couplings_between_dofs(const DoFHandler<3, spacedim> &dof_handler)
      {
        // TODO:[?] Invent significantly better estimates than the ones in this
        // function

        // doing the same thing here is a rather complicated thing, compared
        // to the 2d case, since it is hard to draw pictures with several
        // refined hexahedra :-) so I presently only give a coarse
        // estimate for the case that at most 8 hexes meet at each vertex
        //
        // can anyone give better estimate here?
        const unsigned int max_adjacent_cells =
          dof_handler.tria->max_adjacent_cells();

        types::global_dof_index max_couplings;
        if (max_adjacent_cells <= 8)
          max_couplings =
            7 * 7 * 7 * dof_handler.fe_collection.max_dofs_per_vertex() +
            7 * 6 * 7 * 3 * dof_handler.fe_collection.max_dofs_per_line() +
            9 * 4 * 7 * 3 * dof_handler.fe_collection.max_dofs_per_quad() +
            27 * dof_handler.fe_collection.max_dofs_per_hex();
        else
          {
            Assert(false, ExcNotImplemented());
            max_couplings = 0;
          }

        return std::min(max_couplings, dof_handler.n_dofs());
      }

      /**
       * Reserve enough space in the <tt>levels[]</tt> objects to store the
       * numbers of the degrees of freedom needed for the given element. The
       * given element is that one which was selected when calling
       * @p distribute_dofs the last time.
       */
      template <int spacedim>
      static void reserve_space(DoFHandler<1, spacedim> &dof_handler)
      {
        dof_handler.new_dofs[0][0].resize(
          dof_handler.tria->n_vertices() * dof_handler.get_fe().dofs_per_vertex,
          numbers::invalid_dof_index);

        for (unsigned int i = 0; i < dof_handler.tria->n_levels(); ++i)
          {
            dof_handler.new_dofs[i][1].resize(
              dof_handler.tria->n_raw_cells(i) *
                dof_handler.get_fe().dofs_per_line,
              numbers::invalid_dof_index);

            dof_handler.new_dofs_ptr[i][1].reserve(
              dof_handler.tria->n_raw_cells(i) + 1);
            for (unsigned int j = 0; j < dof_handler.tria->n_raw_cells(i) + 1;
                 j++)
              dof_handler.new_dofs_ptr[i][1].push_back(
                j * dof_handler.get_fe().dofs_per_line);

            dof_handler.new_cell_dofs_cache[i].resize(
              dof_handler.tria->n_raw_cells(i) *
                dof_handler.get_fe().dofs_per_cell,
              numbers::invalid_dof_index);

            dof_handler.new_cell_dofs_cache_ptr[i].reserve(
              dof_handler.tria->n_raw_cells(i) + 1);
            for (unsigned int j = 0; j < dof_handler.tria->n_raw_cells(i) + 1;
                 j++)
              dof_handler.new_cell_dofs_cache_ptr[i].push_back(
                j * dof_handler.get_fe().dofs_per_cell);
          }

        dof_handler.new_dofs[0][0].resize(
          dof_handler.tria->n_vertices() * dof_handler.get_fe().dofs_per_vertex,
          numbers::invalid_dof_index);
      }

      template <int spacedim>
      static void reserve_space(DoFHandler<2, spacedim> &dof_handler)
      {
        dof_handler.new_dofs[0][0].resize(
          dof_handler.tria->n_vertices() * dof_handler.get_fe().dofs_per_vertex,
          numbers::invalid_dof_index);

        for (unsigned int i = 0; i < dof_handler.tria->n_levels(); ++i)
          {
            dof_handler.new_dofs[i][2].resize(
              dof_handler.tria->n_raw_cells(i) *
                dof_handler.get_fe().dofs_per_quad,
              numbers::invalid_dof_index);

            dof_handler.new_dofs_ptr[i][2].reserve(
              dof_handler.tria->n_raw_cells(i) + 1);
            for (unsigned int j = 0; j < dof_handler.tria->n_raw_cells(i) + 1;
                 j++)
              dof_handler.new_dofs_ptr[i][2].push_back(
                j * dof_handler.get_fe().dofs_per_quad);

            dof_handler.new_cell_dofs_cache[i].resize(
              dof_handler.tria->n_raw_cells(i) *
                dof_handler.get_fe().dofs_per_cell,
              numbers::invalid_dof_index);

            dof_handler.new_cell_dofs_cache_ptr[i].reserve(
              dof_handler.tria->n_raw_cells(i) + 1);
            for (unsigned int j = 0; j < dof_handler.tria->n_raw_cells(i) + 1;
                 j++)
              dof_handler.new_cell_dofs_cache_ptr[i].push_back(
                j * dof_handler.get_fe().dofs_per_cell);
          }

        dof_handler.new_dofs[0][0].resize(
          dof_handler.tria->n_vertices() * dof_handler.get_fe().dofs_per_vertex,
          numbers::invalid_dof_index);

        if (dof_handler.tria->n_cells() > 0)
          {
            // line
            dof_handler.new_dofs_ptr[0][1].reserve(
              dof_handler.tria->n_raw_lines() + 1);
            for (unsigned int i = 0; i < dof_handler.tria->n_raw_lines() + 1;
                 i++)
              dof_handler.new_dofs_ptr[0][1].push_back(
                i * dof_handler.get_fe().dofs_per_line);

            dof_handler.new_dofs[0][1].resize(
              dof_handler.tria->n_raw_lines() *
                dof_handler.get_fe().dofs_per_line,
              numbers::invalid_dof_index);
          }
      }

      template <int spacedim>
      static void reserve_space(DoFHandler<3, spacedim> &dof_handler)
      {
        dof_handler.new_dofs[0][0].resize(
          dof_handler.tria->n_vertices() * dof_handler.get_fe().dofs_per_vertex,
          numbers::invalid_dof_index);

        for (unsigned int i = 0; i < dof_handler.tria->n_levels(); ++i)
          {
            dof_handler.new_dofs[i][3].resize(
              dof_handler.tria->n_raw_cells(i) *
                dof_handler.get_fe().dofs_per_hex,
              numbers::invalid_dof_index);

            dof_handler.new_dofs_ptr[i][3].reserve(
              dof_handler.tria->n_raw_cells(i) + 1);
            for (unsigned int j = 0; j < dof_handler.tria->n_raw_cells(i) + 1;
                 j++)
              dof_handler.new_dofs_ptr[i][3].push_back(
                j * dof_handler.get_fe().dofs_per_hex);

            dof_handler.new_cell_dofs_cache[i].resize(
              dof_handler.tria->n_raw_cells(i) *
                dof_handler.get_fe().dofs_per_cell,
              numbers::invalid_dof_index);

            dof_handler.new_cell_dofs_cache_ptr[i].reserve(
              dof_handler.tria->n_raw_cells(i) + 1);
            for (unsigned int j = 0; j < dof_handler.tria->n_raw_cells(i) + 1;
                 j++)
              dof_handler.new_cell_dofs_cache_ptr[i].push_back(
                j * dof_handler.get_fe().dofs_per_cell);
          }

        dof_handler.new_dofs[0][0].resize(
          dof_handler.tria->n_vertices() * dof_handler.get_fe().dofs_per_vertex,
          numbers::invalid_dof_index);

        if (dof_handler.tria->n_cells() > 0)
          {
            // lines
            dof_handler.new_dofs_ptr[0][1].reserve(
              dof_handler.tria->n_raw_lines() + 1);
            for (unsigned int i = 0; i < dof_handler.tria->n_raw_lines() + 1;
                 i++)
              dof_handler.new_dofs_ptr[0][1].push_back(
                i * dof_handler.get_fe().dofs_per_line);

            dof_handler.new_dofs[0][1].resize(
              dof_handler.tria->n_raw_lines() *
                dof_handler.get_fe().dofs_per_line,
              numbers::invalid_dof_index);

            // faces
            dof_handler.new_dofs_ptr[0][2].reserve(
              dof_handler.tria->n_raw_quads() + 1);
            for (unsigned int i = 0; i < dof_handler.tria->n_raw_quads() + 1;
                 i++)
              dof_handler.new_dofs_ptr[0][2].push_back(
                i * dof_handler.get_fe().dofs_per_quad);

            dof_handler.new_dofs[0][2].resize(
              dof_handler.tria->n_raw_quads() *
                dof_handler.get_fe().dofs_per_quad,
              numbers::invalid_dof_index);
          }
      }

      template <int spacedim>
      static void reserve_space_mg(DoFHandler<1, spacedim> &dof_handler)
      {
        Assert(dof_handler.get_triangulation().n_levels() > 0,
               ExcMessage("Invalid triangulation"));
        dof_handler.clear_mg_space();

        const dealii::Triangulation<1, spacedim> &tria =
          dof_handler.get_triangulation();
        const unsigned int dofs_per_line = dof_handler.get_fe().dofs_per_line;
        const unsigned int n_levels      = tria.n_levels();

        for (unsigned int i = 0; i < n_levels; ++i)
          {
            dof_handler.mg_levels.emplace_back(
              new internal::DoFHandlerImplementation::DoFLevel<1>);
            dof_handler.mg_levels.back()->dof_object.dofs =
              std::vector<types::global_dof_index>(tria.n_raw_lines(i) *
                                                     dofs_per_line,
                                                   numbers::invalid_dof_index);
          }

        const unsigned int n_vertices = tria.n_vertices();

        dof_handler.mg_vertex_dofs.resize(n_vertices);

        std::vector<unsigned int> max_level(n_vertices, 0);
        std::vector<unsigned int> min_level(n_vertices, n_levels);

        for (typename dealii::Triangulation<1, spacedim>::cell_iterator cell =
               tria.begin();
             cell != tria.end();
             ++cell)
          {
            const unsigned int level = cell->level();

            for (const unsigned int vertex : GeometryInfo<1>::vertex_indices())
              {
                const unsigned int vertex_index = cell->vertex_index(vertex);

                if (min_level[vertex_index] > level)
                  min_level[vertex_index] = level;

                if (max_level[vertex_index] < level)
                  max_level[vertex_index] = level;
              }
          }

        for (unsigned int vertex = 0; vertex < n_vertices; ++vertex)
          if (tria.vertex_used(vertex))
            {
              Assert(min_level[vertex] < n_levels, ExcInternalError());
              Assert(max_level[vertex] >= min_level[vertex],
                     ExcInternalError());
              dof_handler.mg_vertex_dofs[vertex].init(
                min_level[vertex],
                max_level[vertex],
                dof_handler.get_fe().dofs_per_vertex);
            }

          else
            {
              Assert(min_level[vertex] == n_levels, ExcInternalError());
              Assert(max_level[vertex] == 0, ExcInternalError());
              dof_handler.mg_vertex_dofs[vertex].init(1, 0, 0);
            }
      }

      template <int spacedim>
      static void reserve_space_mg(DoFHandler<2, spacedim> &dof_handler)
      {
        Assert(dof_handler.get_triangulation().n_levels() > 0,
               ExcMessage("Invalid triangulation"));
        dof_handler.clear_mg_space();

        const dealii::FiniteElement<2, spacedim> &fe = dof_handler.get_fe();
        const dealii::Triangulation<2, spacedim> &tria =
          dof_handler.get_triangulation();
        const unsigned int n_levels = tria.n_levels();

        for (unsigned int i = 0; i < n_levels; ++i)
          {
            dof_handler.mg_levels.emplace_back(
              std_cxx14::make_unique<
                internal::DoFHandlerImplementation::DoFLevel<2>>());
            dof_handler.mg_levels.back()->dof_object.dofs =
              std::vector<types::global_dof_index>(tria.n_raw_quads(i) *
                                                     fe.dofs_per_quad,
                                                   numbers::invalid_dof_index);
          }

        dof_handler.mg_faces = std_cxx14::make_unique<
          internal::DoFHandlerImplementation::DoFFaces<2>>();
        dof_handler.mg_faces->lines.dofs = std::vector<types::global_dof_index>(
          tria.n_raw_lines() * fe.dofs_per_line, numbers::invalid_dof_index);

        const unsigned int n_vertices = tria.n_vertices();

        dof_handler.mg_vertex_dofs.resize(n_vertices);

        std::vector<unsigned int> max_level(n_vertices, 0);
        std::vector<unsigned int> min_level(n_vertices, n_levels);

        for (typename dealii::Triangulation<2, spacedim>::cell_iterator cell =
               tria.begin();
             cell != tria.end();
             ++cell)
          {
            const unsigned int level = cell->level();

            for (const unsigned int vertex : GeometryInfo<2>::vertex_indices())
              {
                const unsigned int vertex_index = cell->vertex_index(vertex);

                if (min_level[vertex_index] > level)
                  min_level[vertex_index] = level;

                if (max_level[vertex_index] < level)
                  max_level[vertex_index] = level;
              }
          }

        for (unsigned int vertex = 0; vertex < n_vertices; ++vertex)
          if (tria.vertex_used(vertex))
            {
              Assert(min_level[vertex] < n_levels, ExcInternalError());
              Assert(max_level[vertex] >= min_level[vertex],
                     ExcInternalError());
              dof_handler.mg_vertex_dofs[vertex].init(min_level[vertex],
                                                      max_level[vertex],
                                                      fe.dofs_per_vertex);
            }

          else
            {
              Assert(min_level[vertex] == n_levels, ExcInternalError());
              Assert(max_level[vertex] == 0, ExcInternalError());
              dof_handler.mg_vertex_dofs[vertex].init(1, 0, 0);
            }
      }

      template <int spacedim>
      static void reserve_space_mg(DoFHandler<3, spacedim> &dof_handler)
      {
        Assert(dof_handler.get_triangulation().n_levels() > 0,
               ExcMessage("Invalid triangulation"));
        dof_handler.clear_mg_space();

        const dealii::FiniteElement<3, spacedim> &fe = dof_handler.get_fe();
        const dealii::Triangulation<3, spacedim> &tria =
          dof_handler.get_triangulation();
        const unsigned int n_levels = tria.n_levels();

        for (unsigned int i = 0; i < n_levels; ++i)
          {
            dof_handler.mg_levels.emplace_back(
              std_cxx14::make_unique<
                internal::DoFHandlerImplementation::DoFLevel<3>>());
            dof_handler.mg_levels.back()->dof_object.dofs =
              std::vector<types::global_dof_index>(tria.n_raw_hexs(i) *
                                                     fe.dofs_per_hex,
                                                   numbers::invalid_dof_index);
          }

        dof_handler.mg_faces = std_cxx14::make_unique<
          internal::DoFHandlerImplementation::DoFFaces<3>>();
        dof_handler.mg_faces->lines.dofs = std::vector<types::global_dof_index>(
          tria.n_raw_lines() * fe.dofs_per_line, numbers::invalid_dof_index);
        dof_handler.mg_faces->quads.dofs = std::vector<types::global_dof_index>(
          tria.n_raw_quads() * fe.dofs_per_quad, numbers::invalid_dof_index);

        const unsigned int n_vertices = tria.n_vertices();

        dof_handler.mg_vertex_dofs.resize(n_vertices);

        std::vector<unsigned int> max_level(n_vertices, 0);
        std::vector<unsigned int> min_level(n_vertices, n_levels);

        for (typename dealii::Triangulation<3, spacedim>::cell_iterator cell =
               tria.begin();
             cell != tria.end();
             ++cell)
          {
            const unsigned int level = cell->level();

            for (const unsigned int vertex : GeometryInfo<3>::vertex_indices())
              {
                const unsigned int vertex_index = cell->vertex_index(vertex);

                if (min_level[vertex_index] > level)
                  min_level[vertex_index] = level;

                if (max_level[vertex_index] < level)
                  max_level[vertex_index] = level;
              }
          }

        for (unsigned int vertex = 0; vertex < n_vertices; ++vertex)
          if (tria.vertex_used(vertex))
            {
              Assert(min_level[vertex] < n_levels, ExcInternalError());
              Assert(max_level[vertex] >= min_level[vertex],
                     ExcInternalError());
              dof_handler.mg_vertex_dofs[vertex].init(min_level[vertex],
                                                      max_level[vertex],
                                                      fe.dofs_per_vertex);
            }

          else
            {
              Assert(min_level[vertex] == n_levels, ExcInternalError());
              Assert(max_level[vertex] == 0, ExcInternalError());
              dof_handler.mg_vertex_dofs[vertex].init(1, 0, 0);
            }
      }

      template <int spacedim>
      static types::global_dof_index
      get_dof_index(
        const DoFHandler<1, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<1>>
          &mg_level,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<1>>
          &,
        const unsigned int obj_index,
        const unsigned int fe_index,
        const unsigned int local_index,
        const std::integral_constant<int, 1>)
      {
        Assert(dof_handler.is_hp_dof_handler == false, ExcNotImplemented());

        return mg_level->dof_object.get_dof_index(
          static_cast<const DoFHandler<1, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index);
      }

      template <int spacedim>
      static types::global_dof_index
      get_dof_index(
        const DoFHandler<2, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<2>>
          &,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<2>>
          &                mg_faces,
        const unsigned int obj_index,
        const unsigned int fe_index,
        const unsigned int local_index,
        const std::integral_constant<int, 1>)
      {
        return mg_faces->lines.get_dof_index(
          static_cast<const DoFHandler<2, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index);
      }

      template <int spacedim>
      static types::global_dof_index
      get_dof_index(
        const DoFHandler<2, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<2>>
          &mg_level,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<2>>
          &,
        const unsigned int obj_index,
        const unsigned int fe_index,
        const unsigned int local_index,
        const std::integral_constant<int, 2>)
      {
        Assert(dof_handler.is_hp_dof_handler == false, ExcNotImplemented());
        return mg_level->dof_object.get_dof_index(
          static_cast<const DoFHandler<2, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index);
      }

      template <int spacedim>
      static types::global_dof_index
      get_dof_index(
        const DoFHandler<3, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<3>>
          &,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<3>>
          &                mg_faces,
        const unsigned int obj_index,
        const unsigned int fe_index,
        const unsigned int local_index,
        const std::integral_constant<int, 1>)
      {
        Assert(dof_handler.is_hp_dof_handler == false, ExcNotImplemented());
        return mg_faces->lines.get_dof_index(
          static_cast<const DoFHandler<3, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index);
      }

      template <int spacedim>
      static types::global_dof_index
      get_dof_index(
        const DoFHandler<3, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<3>>
          &,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<3>>
          &                mg_faces,
        const unsigned int obj_index,
        const unsigned int fe_index,
        const unsigned int local_index,
        const std::integral_constant<int, 2>)
      {
        Assert(dof_handler.is_hp_dof_handler == false, ExcNotImplemented());
        return mg_faces->quads.get_dof_index(
          static_cast<const DoFHandler<3, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index);
      }

      template <int spacedim>
      static types::global_dof_index
      get_dof_index(
        const DoFHandler<3, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<3>>
          &mg_level,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<3>>
          &,
        const unsigned int obj_index,
        const unsigned int fe_index,
        const unsigned int local_index,
        const std::integral_constant<int, 3>)
      {
        Assert(dof_handler.is_hp_dof_handler == false, ExcNotImplemented());
        return mg_level->dof_object.get_dof_index(
          static_cast<const DoFHandler<3, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index);
      }

      template <int spacedim>
      static void
      set_dof_index(
        const DoFHandler<1, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<1>>
          &mg_level,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<1>>
          &,
        const unsigned int            obj_index,
        const unsigned int            fe_index,
        const unsigned int            local_index,
        const types::global_dof_index global_index,
        const std::integral_constant<int, 1>)
      {
        Assert(dof_handler.is_hp_dof_handler == false, ExcNotImplemented());
        mg_level->dof_object.set_dof_index(
          static_cast<const DoFHandler<1, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index,
          global_index);
      }

      template <int spacedim>
      static void
      set_dof_index(
        const DoFHandler<2, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<2>>
          &,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<2>>
          &                           mg_faces,
        const unsigned int            obj_index,
        const unsigned int            fe_index,
        const unsigned int            local_index,
        const types::global_dof_index global_index,
        const std::integral_constant<int, 1>)
      {
        Assert(dof_handler.is_hp_dof_handler == false, ExcNotImplemented());
        mg_faces->lines.set_dof_index(
          static_cast<const DoFHandler<2, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index,
          global_index);
      }

      template <int spacedim>
      static void
      set_dof_index(
        const DoFHandler<2, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<2>>
          &mg_level,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<2>>
          &,
        const unsigned int            obj_index,
        const unsigned int            fe_index,
        const unsigned int            local_index,
        const types::global_dof_index global_index,
        const std::integral_constant<int, 2>)
      {
        Assert(dof_handler.is_hp_dof_handler == false, ExcNotImplemented());
        mg_level->dof_object.set_dof_index(
          static_cast<const DoFHandler<2, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index,
          global_index);
      }

      template <int spacedim>
      static void
      set_dof_index(
        const DoFHandler<3, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<3>>
          &,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<3>>
          &                           mg_faces,
        const unsigned int            obj_index,
        const unsigned int            fe_index,
        const unsigned int            local_index,
        const types::global_dof_index global_index,
        const std::integral_constant<int, 1>)
      {
        Assert(dof_handler.is_hp_dof_handler == false, ExcNotImplemented());
        mg_faces->lines.set_dof_index(
          static_cast<const DoFHandler<3, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index,
          global_index);
      }

      template <int spacedim>
      static void
      set_dof_index(
        const DoFHandler<3, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<3>>
          &,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<3>>
          &                           mg_faces,
        const unsigned int            obj_index,
        const unsigned int            fe_index,
        const unsigned int            local_index,
        const types::global_dof_index global_index,
        const std::integral_constant<int, 2>)
      {
        Assert(dof_handler.is_hp_dof_handler == false, ExcNotImplemented());
        mg_faces->quads.set_dof_index(
          static_cast<const DoFHandler<3, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index,
          global_index);
      }

      template <int spacedim>
      static void
      set_dof_index(
        const DoFHandler<3, spacedim> &dof_handler,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFLevel<3>>
          &mg_level,
        const std::unique_ptr<internal::DoFHandlerImplementation::DoFFaces<3>>
          &,
        const unsigned int            obj_index,
        const unsigned int            fe_index,
        const unsigned int            local_index,
        const types::global_dof_index global_index,
        const std::integral_constant<int, 3>)
      {
        Assert(dof_handler.is_hp_dof_handler == false, ExcNotImplemented());
        mg_level->dof_object.set_dof_index(
          static_cast<const DoFHandler<3, spacedim> &>(dof_handler),
          obj_index,
          fe_index,
          local_index,
          global_index);
      }
    };
  } // namespace DoFHandlerImplementation

  namespace hp
  {
    namespace DoFHandlerImplementation
    {
      /**
       * A class with the same purpose as the similarly named class of the
       * Triangulation class. See there for more information.
       */
      struct Implementation
      {
        /**
         * No future_fe_indices should have been assigned when partitioning a
         * triangulation, since they are only available locally and will not be
         * communicated.
         */
        template <int dim, int spacedim>
        static void
        ensure_absence_of_future_fe_indices(
          DoFHandler<dim, spacedim> &dof_handler)
        {
          (void)dof_handler;
          for (const auto &cell : dof_handler.active_cell_iterators())
            if (cell->is_locally_owned())
              Assert(
                !cell->future_fe_index_set(),
                ExcMessage(
                  "There shouldn't be any cells flagged for p-adaptation when partitioning."));
        }



        /**
         * Do that part of reserving space that pertains to releasing
         * the previously used memory.
         */
        template <int dim, int spacedim>
        static void
        reserve_space_release_space(DoFHandler<dim, spacedim> &dof_handler)
        {
          // Release all space except the fields for active_fe_indices and
          // refinement flags which we have to back up before
          {
            std::vector<std::vector<
              typename DoFHandler<dim, spacedim>::active_fe_index_type>>
              active_fe_backup(dof_handler.new_active_fe_indices.size()),
              future_fe_backup(dof_handler.new_future_fe_indices.size());
            for (unsigned int level = 0;
                 level < dof_handler.new_future_fe_indices.size();
                 ++level)
              {
                active_fe_backup[level] =
                  std::move(dof_handler.new_active_fe_indices[level]);
                future_fe_backup[level] =
                  std::move(dof_handler.new_future_fe_indices[level]);
              }

            // delete all levels and set them up newly, since vectors
            // are troublesome if you want to change their size
            dof_handler.clear_space();

            dof_handler.new_dofs.resize(dof_handler.tria->n_levels());
            dof_handler.new_dofs_ptr.resize(dof_handler.tria->n_levels());
            dof_handler.new_cell_dofs_cache.resize(
              dof_handler.tria->n_levels());
            dof_handler.new_cell_dofs_cache_ptr.resize(
              dof_handler.tria->n_levels());
            dof_handler.new_active_fe_indices.resize(
              dof_handler.tria->n_levels());
            dof_handler.new_future_fe_indices.resize(
              dof_handler.tria->n_levels());

            for (unsigned int level = 0; level < dof_handler.tria->n_levels();
                 ++level)
              {
                // recover backups
                dof_handler.new_active_fe_indices[level] =
                  std::move(active_fe_backup[level]);
                dof_handler.new_future_fe_indices[level] =
                  std::move(future_fe_backup[level]);
              }
          }
        }



        /**
         * Do that part of reserving space that pertains to vertices,
         * since this is the same in all space dimensions.
         */
        template <int dim, int spacedim>
        static void
        reserve_space_vertices(DoFHandler<dim, spacedim> &dof_handler)
        {
          // The final step in all of the reserve_space() functions is to set
          // up vertex dof information. since vertices are sequentially
          // numbered, what we do first is to set up an array in which
          // we record whether a vertex is associated with any of the
          // given fe's, by setting a bit. in a later step, we then
          // actually allocate memory for the required dofs
          //
          // in the following, we only need to consider vertices that are
          // adjacent to either a locally owned or a ghost cell; we never
          // store anything on vertices that are only surrounded by
          // artificial cells. so figure out that subset of vertices
          // first
          std::vector<bool> locally_used_vertices(
            dof_handler.tria->n_vertices(), false);
          for (const auto &cell : dof_handler.active_cell_iterators())
            if (!cell->is_artificial())
              for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
                locally_used_vertices[cell->vertex_index(v)] = true;

          std::vector<std::vector<bool>> vertex_fe_association(
            dof_handler.fe_collection.size(),
            std::vector<bool>(dof_handler.tria->n_vertices(), false));

          for (const auto &cell : dof_handler.active_cell_iterators())
            if (!cell->is_artificial())
              for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
                vertex_fe_association[cell->active_fe_index()]
                                     [cell->vertex_index(v)] = true;

                // in debug mode, make sure that each vertex is associated
                // with at least one fe (note that except for unused
                // vertices, all vertices are actually active). this is of
                // course only true for vertices that are part of either
                // ghost or locally owned cells
#ifdef DEBUG
          for (unsigned int v = 0; v < dof_handler.tria->n_vertices(); ++v)
            if (locally_used_vertices[v] == true)
              if (dof_handler.tria->vertex_used(v) == true)
                {
                  unsigned int fe = 0;
                  for (; fe < dof_handler.fe_collection.size(); ++fe)
                    if (vertex_fe_association[fe][v] == true)
                      break;
                  Assert(fe != dof_handler.fe_collection.size(),
                         ExcInternalError());
                }
#endif

          const unsigned int d = 0;
          const unsigned int l = 0;

          dof_handler.new_hp_ptr[d].clear();
          dof_handler.new_hp_fe[d].clear();
          dof_handler.new_dofs_ptr[l][d].clear();
          dof_handler.new_dofs[l][d].clear();

          dof_handler.new_hp_ptr[d].reserve(dof_handler.tria->n_vertices() + 1);

          unsigned int vertex_slots_needed = 0;
          unsigned int fe_slots_needed     = 0;

          for (unsigned int v = 0; v < dof_handler.tria->n_vertices(); ++v)
            {
              dof_handler.new_hp_ptr[d].push_back(fe_slots_needed);

              if (dof_handler.tria->vertex_used(v) && locally_used_vertices[v])
                {
                  for (unsigned int fe = 0;
                       fe < dof_handler.fe_collection.size();
                       ++fe)
                    if (vertex_fe_association[fe][v] == true)
                      {
                        fe_slots_needed++;
                        vertex_slots_needed +=
                          dof_handler.get_fe(fe).dofs_per_vertex;
                      }
                }
            }

          dof_handler.new_hp_ptr[d].push_back(fe_slots_needed);

          dof_handler.new_hp_fe[d].reserve(fe_slots_needed);
          dof_handler.new_dofs_ptr[l][d].reserve(fe_slots_needed + 1);

          dof_handler.new_dofs[l][d].reserve(vertex_slots_needed);

          for (unsigned int v = 0; v < dof_handler.tria->n_vertices(); ++v)
            if (dof_handler.tria->vertex_used(v) && locally_used_vertices[v])
              {
                for (unsigned int fe = 0; fe < dof_handler.fe_collection.size();
                     ++fe)
                  if (vertex_fe_association[fe][v] == true)
                    {
                      dof_handler.new_hp_fe[d].push_back(fe);
                      dof_handler.new_dofs_ptr[l][d].push_back(
                        dof_handler.new_dofs[l][d].size());

                      for (unsigned int i = 0;
                           i < dof_handler.get_fe(fe).dofs_per_vertex;
                           i++)
                        dof_handler.new_dofs[l][d].push_back(
                          numbers::invalid_dof_index);
                    }
              }


          dof_handler.new_dofs_ptr[l][d].push_back(
            dof_handler.new_dofs[l][d].size());

          AssertDimension(vertex_slots_needed,
                          dof_handler.new_dofs[l][d].size());
          AssertDimension(fe_slots_needed, dof_handler.new_hp_fe[d].size());
          AssertDimension(fe_slots_needed + 1,
                          dof_handler.new_dofs_ptr[l][d].size());
          AssertDimension(dof_handler.tria->n_vertices() + 1,
                          dof_handler.new_hp_ptr[d].size());

          dof_handler.new_dofs[l][d].assign(vertex_slots_needed,
                                            numbers::invalid_dof_index);
        }



        /**
         * Do that part of reserving space that pertains to cells,
         * since this is the same in all space dimensions.
         */
        template <int dim, int spacedim>
        static void
        reserve_space_cells(DoFHandler<dim, spacedim> &dof_handler)
        {
          (void)dof_handler;
          // count how much space we need on each level for the cell
          // dofs and set the dof_*_offsets data. initially set the
          // latter to an invalid index, and only later set it to
          // something reasonable for active dof_handler.cells
          //
          // note that for dof_handler.cells, the situation is simpler
          // than for other (lower dimensional) objects since exactly
          // one finite element is used for it
          for (unsigned int level = 0; level < dof_handler.tria->n_levels();
               ++level)
            {
              dof_handler.new_dofs_ptr[level][dim] =
                std::vector<typename DoFHandler<dim, spacedim>::offset_type>(
                  dof_handler.tria->n_raw_cells(level),
                  static_cast<typename DoFHandler<dim, spacedim>::offset_type>(
                    -1));
              dof_handler.new_cell_dofs_cache_ptr[level] =
                std::vector<typename DoFHandler<dim, spacedim>::offset_type>(
                  dof_handler.tria->n_raw_cells(level),
                  static_cast<typename DoFHandler<dim, spacedim>::offset_type>(
                    -1));

              types::global_dof_index next_free_dof = 0;
              types::global_dof_index cache_size    = 0;

              for (auto cell :
                   dof_handler.active_cell_iterators_on_level(level))
                if (cell->is_active() && !cell->is_artificial())
                  {
                    dof_handler.new_dofs_ptr[level][dim][cell->index()] =
                      next_free_dof;
                    next_free_dof +=
                      cell->get_fe().template n_dofs_per_object<dim>();

                    dof_handler.new_cell_dofs_cache_ptr[level][cell->index()] =
                      cache_size;
                    cache_size += cell->get_fe().dofs_per_cell;
                  }

              dof_handler.new_dofs[level][dim] =
                std::vector<types::global_dof_index>(
                  next_free_dof, numbers::invalid_dof_index);
              dof_handler.new_cell_dofs_cache[level] =
                std::vector<types::global_dof_index>(
                  cache_size, numbers::invalid_dof_index);
            }
        }



        /**
         * Do that part of reserving space that pertains to faces,
         * since this is the same in all space dimensions.
         */
        template <int dim, int spacedim>
        static void
        reserve_space_faces(DoFHandler<dim, spacedim> &dof_handler)
        {
          // FACE DOFS
          //
          // Count face dofs, then allocate as much space
          // as we need and prime the linked list for faces (see the
          // description in hp::DoFLevel) with the indices we will
          // need. Note that our task is more complicated than for the
          // cell case above since two adjacent cells may have different
          // active_fe_indices, in which case we need to allocate
          // *two* sets of face dofs for the same face. But they don't
          // *have* to be different, and so we need to prepare for this
          // as well.
          //
          // The way we do things is that we loop over all active
          // cells (these are the only ones that have DoFs
          // anyway) and all their faces. We note in the
          // user flags whether we have previously visited a face and
          // if so skip it (consequently, we have to save and later
          // restore the face flags)
          {
            std::vector<bool> saved_face_user_flags;
            switch (dim)
              {
                case 2:
                  {
                    const_cast<dealii::Triangulation<dim, spacedim> &>(
                      *dof_handler.tria)
                      .save_user_flags_line(saved_face_user_flags);
                    const_cast<dealii::Triangulation<dim, spacedim> &>(
                      *dof_handler.tria)
                      .clear_user_flags_line();

                    break;
                  }

                case 3:
                  {
                    const_cast<dealii::Triangulation<dim, spacedim> &>(
                      *dof_handler.tria)
                      .save_user_flags_quad(saved_face_user_flags);
                    const_cast<dealii::Triangulation<dim, spacedim> &>(
                      *dof_handler.tria)
                      .clear_user_flags_quad();

                    break;
                  }

                default:
                  Assert(false, ExcNotImplemented());
              }

            const unsigned int d = dim - 1;
            const unsigned int l = 0;

            dof_handler.new_hp_ptr[d].clear();
            dof_handler.new_hp_fe[d].clear();
            dof_handler.new_dofs_ptr[l][d].clear();
            dof_handler.new_dofs[l][d].clear();

            dof_handler.new_hp_ptr[d].resize(dof_handler.tria->n_raw_faces() +
                                             1);

            // An array to hold how many slots (see the hp::DoFLevel
            // class) we will have to store on each level
            unsigned int n_face_slots = 0;

            for (const auto &cell : dof_handler.active_cell_iterators())
              if (!cell->is_artificial())
                for (const unsigned int face :
                     GeometryInfo<dim>::face_indices())
                  if (cell->face(face)->user_flag_set() == false)
                    {
                      unsigned int fe_slots_needed = 0;

                      if (cell->at_boundary(face) ||
                          cell->face(face)->has_children() ||
                          cell->neighbor_is_coarser(face) ||
                          (!cell->at_boundary(face) &&
                           cell->neighbor(face)->is_artificial()) ||
                          (!cell->at_boundary(face) &&
                           !cell->neighbor(face)->is_artificial() &&
                           (cell->active_fe_index() ==
                            cell->neighbor(face)->active_fe_index())))
                        {
                          fe_slots_needed = 1;
                          n_face_slots +=
                            dof_handler.get_fe(cell->active_fe_index())
                              .template n_dofs_per_object<dim - 1>();
                        }
                      else
                        {
                          fe_slots_needed = 2;
                          n_face_slots +=
                            dof_handler.get_fe(cell->active_fe_index())
                              .template n_dofs_per_object<dim - 1>() +
                            dof_handler
                              .get_fe(cell->neighbor(face)->active_fe_index())
                              .template n_dofs_per_object<dim - 1>();
                        }

                      // mark this face as visited
                      cell->face(face)->set_user_flag();

                      dof_handler.new_hp_ptr[d][cell->face(face)->index() + 1] =
                        fe_slots_needed;
                    }

            for (unsigned int i = 1; i < dof_handler.new_hp_ptr[d].size(); i++)
              dof_handler.new_hp_ptr[d][i] += dof_handler.new_hp_ptr[d][i - 1];


            dof_handler.new_hp_fe[d].resize(dof_handler.new_hp_ptr[d].back());
            dof_handler.new_dofs_ptr[l][d].resize(
              dof_handler.new_hp_ptr[d].back() + 1);

            dof_handler.new_dofs[l][d].reserve(n_face_slots);


            // With the memory now allocated, loop over the
            // dof_handler cells again and prime the _offset values as
            // well as the fe_index fields
            switch (dim)
              {
                case 2:
                  {
                    const_cast<dealii::Triangulation<dim, spacedim> &>(
                      *dof_handler.tria)
                      .clear_user_flags_line();

                    break;
                  }

                case 3:
                  {
                    const_cast<dealii::Triangulation<dim, spacedim> &>(
                      *dof_handler.tria)
                      .clear_user_flags_quad();

                    break;
                  }

                default:
                  Assert(false, ExcNotImplemented());
              }

            for (const auto &cell : dof_handler.active_cell_iterators())
              if (!cell->is_artificial())
                for (const unsigned int face :
                     GeometryInfo<dim>::face_indices())
                  if (!cell->face(face)->user_flag_set())
                    {
                      // Same decision tree as before
                      if (cell->at_boundary(face) ||
                          cell->face(face)->has_children() ||
                          cell->neighbor_is_coarser(face) ||
                          (!cell->at_boundary(face) &&
                           cell->neighbor(face)->is_artificial()) ||
                          (!cell->at_boundary(face) &&
                           !cell->neighbor(face)->is_artificial() &&
                           (cell->active_fe_index() ==
                            cell->neighbor(face)->active_fe_index())))
                        {
                          const unsigned int fe = cell->active_fe_index();
                          const unsigned int n_dofs =
                            dof_handler.get_fe(fe)
                              .template n_dofs_per_object<dim - 1>();
                          const unsigned int offset =
                            dof_handler
                              .new_hp_ptr[d][cell->face(face)->index()];

                          dof_handler.new_hp_fe[d][offset]           = fe;
                          dof_handler.new_dofs_ptr[l][d][offset + 1] = n_dofs;

                          for (unsigned int i = 0; i < n_dofs; i++)
                            dof_handler.new_dofs[l][d].push_back(
                              numbers::invalid_dof_index);
                        }
                      else
                        {
                          unsigned int fe_1 = cell->active_fe_index();
                          unsigned int fe_2 =
                            cell->neighbor(face)->active_fe_index();

                          if (fe_2 < fe_1)
                            std::swap(fe_1, fe_2);

                          const unsigned int n_dofs_1 =
                            dof_handler.get_fe(fe_1)
                              .template n_dofs_per_object<dim - 1>();

                          const unsigned int n_dofs_2 =
                            dof_handler.get_fe(fe_2)
                              .template n_dofs_per_object<dim - 1>();

                          const unsigned int offset =
                            dof_handler
                              .new_hp_ptr[d][cell->face(face)->index()];

                          dof_handler.new_hp_fe[d].push_back(
                            cell->active_fe_index());
                          dof_handler.new_dofs_ptr[l][d].push_back(
                            dof_handler.new_dofs[l][d].size());

                          for (unsigned int i = 0;
                               i < dof_handler.get_fe(cell->active_fe_index())
                                     .template n_dofs_per_object<dim - 1>();
                               i++)
                            dof_handler.new_dofs[l][d].push_back(
                              numbers::invalid_dof_index);

                          dof_handler.new_hp_fe[d][offset + 0]       = fe_1;
                          dof_handler.new_hp_fe[d][offset + 1]       = fe_2;
                          dof_handler.new_dofs_ptr[l][d][offset + 1] = n_dofs_1;
                          dof_handler.new_dofs_ptr[l][d][offset + 2] = n_dofs_2;


                          for (unsigned int i = 0; i < n_dofs_1 + n_dofs_2; i++)
                            dof_handler.new_dofs[l][d].push_back(
                              numbers::invalid_dof_index);
                        }

                      // mark this face as visited
                      cell->face(face)->set_user_flag();
                    }

            for (unsigned int i = 1; i < dof_handler.new_dofs_ptr[l][d].size();
                 i++)
              dof_handler.new_dofs_ptr[l][d][i] +=
                dof_handler.new_dofs_ptr[l][d][i - 1];

            // at the end, restore the user flags for the faces
            switch (dim)
              {
                case 2:
                  {
                    const_cast<dealii::Triangulation<dim, spacedim> &>(
                      *dof_handler.tria)
                      .load_user_flags_line(saved_face_user_flags);

                    break;
                  }

                case 3:
                  {
                    const_cast<dealii::Triangulation<dim, spacedim> &>(
                      *dof_handler.tria)
                      .load_user_flags_quad(saved_face_user_flags);

                    break;
                  }

                default:
                  Assert(false, ExcNotImplemented());
              }
          }
        }



        /**
         * Reserve enough space in the <tt>levels[]</tt> objects to
         * store the numbers of the degrees of freedom needed for the
         * given element. The given element is that one which was
         * selected when calling @p distribute_dofs the last time.
         */
        template <int spacedim>
        static void reserve_space(dealii::DoFHandler<1, spacedim> &dof_handler)
        {
          Assert(dof_handler.fe_collection.size() > 0,
                 (typename dealii::DoFHandler<1, spacedim>::ExcNoFESelected()));
          Assert(dof_handler.tria->n_levels() > 0,
                 ExcMessage("The current Triangulation must not be empty."));
          Assert(dof_handler.tria->n_levels() ==
                   dof_handler.new_future_fe_indices.size(),
                 ExcInternalError());

          reserve_space_release_space(dof_handler);

          Threads::TaskGroup<> tasks;
          tasks +=
            Threads::new_task(&reserve_space_cells<1, spacedim>, dof_handler);
          tasks += Threads::new_task(&reserve_space_vertices<1, spacedim>,
                                     dof_handler);
          tasks.join_all();
        }



        template <int spacedim>
        static void reserve_space(dealii::DoFHandler<2, spacedim> &dof_handler)
        {
          Assert(dof_handler.fe_collection.size() > 0,
                 (typename dealii::DoFHandler<1, spacedim>::ExcNoFESelected()));
          Assert(dof_handler.tria->n_levels() > 0,
                 ExcMessage("The current Triangulation must not be empty."));
          Assert(dof_handler.tria->n_levels() ==
                   dof_handler.new_future_fe_indices.size(),
                 ExcInternalError());

          reserve_space_release_space(dof_handler);

          Threads::TaskGroup<> tasks;
          tasks +=
            Threads::new_task(&reserve_space_cells<2, spacedim>, dof_handler);
          tasks +=
            Threads::new_task(&reserve_space_faces<2, spacedim>, dof_handler);
          tasks += Threads::new_task(&reserve_space_vertices<2, spacedim>,
                                     dof_handler);
          tasks.join_all();
        }



        template <int spacedim>
        static void reserve_space(dealii::DoFHandler<3, spacedim> &dof_handler)
        {
          const unsigned int dim = 3;

          Assert(dof_handler.fe_collection.size() > 0,
                 (typename dealii::DoFHandler<1, spacedim>::ExcNoFESelected()));
          Assert(dof_handler.tria->n_levels() > 0,
                 ExcMessage("The current Triangulation must not be empty."));
          Assert(dof_handler.tria->n_levels() ==
                   dof_handler.new_future_fe_indices.size(),
                 ExcInternalError());

          reserve_space_release_space(dof_handler);

          Threads::TaskGroup<> tasks;
          tasks +=
            Threads::new_task(&reserve_space_cells<3, spacedim>, dof_handler);
          tasks +=
            Threads::new_task(&reserve_space_faces<3, spacedim>, dof_handler);
          tasks += Threads::new_task(&reserve_space_vertices<3, spacedim>,
                                     dof_handler);

          // While the tasks above are running, we can turn to line dofs

          // the situation here is pretty much like with vertices:
          // there can be an arbitrary number of finite elements
          // associated with each line.
          //
          // the algorithm we use is somewhat similar to what we do in
          // reserve_space_vertices()
          {
            // what we do first is to set up an array in which we
            // record whether a line is associated with any of the
            // given fe's, by setting a bit. in a later step, we
            // then actually allocate memory for the required dofs
            std::vector<std::vector<bool>> line_fe_association(
              dof_handler.fe_collection.size(),
              std::vector<bool>(dof_handler.tria->n_raw_lines(), false));

            for (const auto &cell : dof_handler.active_cell_iterators())
              if (!cell->is_artificial())
                for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell;
                     ++l)
                  line_fe_association[cell->active_fe_index()]
                                     [cell->line_index(l)] = true;

            // first check which of the lines is used at all,
            // i.e. is associated with a finite element. we do this
            // since not all lines may actually be used, in which
            // case we do not have to allocate any memory at all
            std::vector<bool> line_is_used(dof_handler.tria->n_raw_lines(),
                                           false);
            for (unsigned int line = 0; line < dof_handler.tria->n_raw_lines();
                 ++line)
              for (unsigned int fe = 0; fe < dof_handler.fe_collection.size();
                   ++fe)
                if (line_fe_association[fe][line] == true)
                  {
                    line_is_used[line] = true;
                    break;
                  }



            const unsigned int d = 1;
            const unsigned int l = 0;

            dof_handler.new_hp_ptr[d].clear();
            dof_handler.new_hp_fe[d].clear();
            dof_handler.new_dofs_ptr[l][d].clear();
            dof_handler.new_dofs[l][d].clear();

            dof_handler.new_hp_ptr[d].reserve(dof_handler.tria->n_raw_lines() +
                                              1);

            unsigned int line_slots_needed = 0;
            unsigned int fe_slots_needed   = 0;

            for (unsigned int line = 0; line < dof_handler.tria->n_raw_lines();
                 ++line)
              {
                dof_handler.new_hp_ptr[d].push_back(fe_slots_needed);

                if (line_is_used[line] == true)
                  {
                    for (unsigned int fe = 0;
                         fe < dof_handler.fe_collection.size();
                         ++fe)
                      if (line_fe_association[fe][line] == true)
                        {
                          fe_slots_needed++;
                          line_slots_needed +=
                            dof_handler.get_fe(fe).dofs_per_line;
                        }
                  }
              }

            dof_handler.new_hp_ptr[d].push_back(fe_slots_needed);

            dof_handler.new_hp_fe[d].reserve(fe_slots_needed);
            dof_handler.new_dofs_ptr[l][d].reserve(fe_slots_needed + 1);

            dof_handler.new_dofs[l][d].reserve(line_slots_needed);

            for (unsigned int line = 0; line < dof_handler.tria->n_raw_lines();
                 ++line)
              if (line_is_used[line] == true)
                {
                  for (unsigned int fe = 0;
                       fe < dof_handler.fe_collection.size();
                       ++fe)
                    if (line_fe_association[fe][line] == true)
                      {
                        dof_handler.new_hp_fe[d].push_back(fe);
                        dof_handler.new_dofs_ptr[l][d].push_back(
                          dof_handler.new_dofs[l][d].size());

                        for (unsigned int i = 0;
                             i < dof_handler.get_fe(fe).dofs_per_line;
                             i++)
                          dof_handler.new_dofs[l][d].push_back(
                            numbers::invalid_dof_index);
                      }
                }
          }

          // Ensure that everything is done at this point.
          tasks.join_all();
        }



        /**
         * Given a hp::DoFHandler object, make sure that the active_fe_indices
         * that a user has set for locally owned cells are communicated to all
         * other relevant cells as well.
         *
         * For parallel::shared::Triangulation objects,
         * this information is distributed on both ghost and artificial cells.
         *
         * In case a parallel::distributed::Triangulation is used,
         * indices are communicated only to ghost cells.
         */
        template <int dim, int spacedim>
        static void
        communicate_active_fe_indices(DoFHandler<dim, spacedim> &dof_handler)
        {
          Assert(dof_handler.is_hp_dof_handler == true, ExcNotImplemented());

          if (const dealii::parallel::shared::Triangulation<dim, spacedim> *tr =
                dynamic_cast<
                  const dealii::parallel::shared::Triangulation<dim, spacedim>
                    *>(&dof_handler.get_triangulation()))
            {
              // we have a shared triangulation. in this case, every processor
              // knows about all cells, but every processor only has knowledge
              // about the active_fe_index on the cells it owns.
              //
              // we can create a complete set of active_fe_indices by letting
              // every processor create a vector of indices for all cells,
              // filling only those on the cells it owns and setting the indices
              // on the other cells to zero. then we add all of these vectors
              // up, and because every vector entry has exactly one processor
              // that owns it, the sum is correct
              std::vector<unsigned int> active_fe_indices(tr->n_active_cells(),
                                                          0u);
              for (const auto &cell : dof_handler.active_cell_iterators())
                if (cell->is_locally_owned())
                  active_fe_indices[cell->active_cell_index()] =
                    cell->active_fe_index();

              Utilities::MPI::sum(active_fe_indices,
                                  tr->get_communicator(),
                                  active_fe_indices);

              // now go back and fill the active_fe_index on all other
              // cells. we would like to call cell->set_active_fe_index(),
              // but that function does not allow setting these indices on
              // non-locally_owned cells. so we have to work around the
              // issue a little bit by accessing the underlying data
              // structures directly
              for (const auto &cell : dof_handler.active_cell_iterators())
                if (!cell->is_locally_owned())
                  dof_handler
                    .new_active_fe_indices[cell->level()][cell->index()] =
                    active_fe_indices[cell->active_cell_index()];
            }
          else if (const dealii::parallel::distributed::Triangulation<dim,
                                                                      spacedim>
                     *tr = dynamic_cast<const dealii::parallel::distributed::
                                          Triangulation<dim, spacedim> *>(
                       &dof_handler.get_triangulation()))
            {
              // For completely distributed meshes, use the function that is
              // able to move data from locally owned cells on one processor to
              // the corresponding ghost cells on others. To this end, we need
              // to have functions that can pack and unpack the data we want to
              // transport -- namely, the single unsigned int active_fe_index
              // objects
              auto pack = [](const typename dealii::DoFHandler<dim, spacedim>::
                               active_cell_iterator &cell) -> unsigned int {
                return cell->active_fe_index();
              };

              auto unpack = [&dof_handler](
                              const typename dealii::DoFHandler<dim, spacedim>::
                                active_cell_iterator &cell,
                              const unsigned int      active_fe_index) -> void {
                // we would like to say
                //   cell->set_active_fe_index(active_fe_index);
                // but this is not allowed on cells that are not
                // locally owned, and we are on a ghost cell
                dof_handler
                  .new_active_fe_indices[cell->level()][cell->index()] =
                  active_fe_index;
              };

              GridTools::exchange_cell_data_to_ghosts<
                unsigned int,
                dealii::DoFHandler<dim, spacedim>>(
                static_cast<dealii::DoFHandler<dim, spacedim> &>(dof_handler),
                pack,
                unpack);
            }
          else
            {
              // a sequential triangulation. there is nothing we need to do here
              Assert(
                (dynamic_cast<
                   const dealii::parallel::TriangulationBase<dim, spacedim> *>(
                   &dof_handler.get_triangulation()) == nullptr),
                ExcInternalError());
            }
        }



        /**
         * Collect all finite element indices on cells that will be affected by
         * future refinement and coarsening. Further, prepare those indices to
         * be distributed on on the updated triangulation later.
         *
         * On cells to be refined, the active_fe_index will be inherited to
         * their children and thus will be stored as such.
         *
         * On cells to be coarsened, we choose the finite element on the parent
         * cell from those assigned to their children to be the one that is
         * dominated by all children. If none was found, we pick the most
         * dominant element in the whole collection that is dominated by all
         * children. See documentation of
         * hp::FECollection::find_dominated_fe_extended() for further
         * information.
         *
         * On cells intended for p-refinement or p-coarsening, those
         * active_fe_indices will be determined by the corresponding flags that
         * have been set on the relevant cells.
         */
        template <int dim, int spacedim>
        static void
        collect_fe_indices_on_cells_to_be_refined(
          DoFHandler<dim, spacedim> &dof_handler)
        {
          const auto &fe_transfer = dof_handler.active_fe_index_transfer;

          for (const auto &cell : dof_handler.active_cell_iterators())
            if (cell->is_locally_owned())
              {
                if (cell->refine_flag_set())
                  {
                    // Store the active_fe_index of each cell that will be
                    // refined to and distribute it later on its children.
                    // Pick their future index if flagged for p-refinement.
                    fe_transfer->refined_cells_fe_index.insert(
                      {cell, cell->future_fe_index()});
                  }
                else if (cell->coarsen_flag_set())
                  {
                    // From all cells that will be coarsened, determine their
                    // parent and calculate its proper active_fe_index, so that
                    // it can be set after refinement. But first, check if that
                    // particular cell has a parent at all.
                    Assert(cell->level() > 0, ExcInternalError());
                    const auto &parent = cell->parent();

                    // Check if the active_fe_index for the current cell has
                    // been determined already.
                    if (fe_transfer->coarsened_cells_fe_index.find(parent) ==
                        fe_transfer->coarsened_cells_fe_index.end())
                      {
                        // Find a suitable active_fe_index for the parent cell
                        // based on the 'least dominant finite element' of its
                        // children. Consider the childrens' hypothetical future
                        // index when they have been flagged for p-refinement.
                        std::set<unsigned int> fe_indices_children;
                        for (unsigned int child_index = 0;
                             child_index < parent->n_children();
                             ++child_index)
                          {
                            const auto sibling = parent->child(child_index);
                            Assert(sibling->is_active() &&
                                     sibling->coarsen_flag_set(),
                                   typename dealii::Triangulation<
                                     dim>::ExcInconsistentCoarseningFlags());

                            fe_indices_children.insert(
                              sibling->future_fe_index());
                          }
                        Assert(!fe_indices_children.empty(),
                               ExcInternalError());

                        const unsigned int fe_index =
                          dof_handler.fe_collection.find_dominated_fe_extended(
                            fe_indices_children, /*codim=*/0);

                        Assert(fe_index != numbers::invalid_unsigned_int,
                               typename dealii::hp::FECollection<dim>::
                                 ExcNoDominatedFiniteElementAmongstChildren());

                        fe_transfer->coarsened_cells_fe_index.insert(
                          {parent, fe_index});
                      }
                  }
                else
                  {
                    // No h-refinement is scheduled for this cell.
                    // However, it may have p-refinement indicators, so we
                    // choose a new active_fe_index based on its flags.
                    if (cell->future_fe_index_set() == true)
                      fe_transfer->persisting_cells_fe_index.insert(
                        {cell, cell->future_fe_index()});
                  }
              }
        }



        /**
         * Distribute active finite element indices that have been previously
         * prepared in collect_fe_indices_on_cells_to_be_refined().
         */
        template <int dim, int spacedim>
        static void
        distribute_fe_indices_on_refined_cells(
          DoFHandler<dim, spacedim> &dof_handler)
        {
          const auto &fe_transfer = dof_handler.active_fe_index_transfer;

          // Set active_fe_indices on persisting cells.
          for (const auto &persist : fe_transfer->persisting_cells_fe_index)
            {
              const auto &cell = persist.first;

              if (cell->is_locally_owned())
                {
                  Assert(cell->is_active(), ExcInternalError());
                  cell->set_active_fe_index(persist.second);
                }
            }

          // Distribute active_fe_indices from all refined cells on their
          // respective children.
          for (const auto &refine : fe_transfer->refined_cells_fe_index)
            {
              const auto &parent = refine.first;

              for (unsigned int child_index = 0;
                   child_index < parent->n_children();
                   ++child_index)
                {
                  const auto &child = parent->child(child_index);
                  Assert(child->is_locally_owned() && child->is_active(),
                         ExcInternalError());
                  child->set_active_fe_index(refine.second);
                }
            }

          // Set active_fe_indices on coarsened cells that have been determined
          // before the actual coarsening happened.
          for (const auto &coarsen : fe_transfer->coarsened_cells_fe_index)
            {
              const auto &cell = coarsen.first;
              Assert(cell->is_locally_owned() && cell->is_active(),
                     ExcInternalError());
              cell->set_active_fe_index(coarsen.second);
            }
        }


        /**
         * Coarsening strategy for the CellDataTransfer object responsible for
         * tranferring the active_fe_index of each cell on
         * parallel::distributed::Triangulation objects that have been refined.
         *
         * A finite element index needs to be determined for the (not yet
         * active) parent cell from its (still active) children.  Out of the set
         * of elements previously assigned to the former children, we choose the
         * one dominated by all children for the parent cell.
         */
        template <int dim, int spacedim>
        static unsigned int
        determine_fe_from_children(
          const std::vector<unsigned int> &        children_fe_indices,
          dealii::hp::FECollection<dim, spacedim> &fe_collection)
        {
          Assert(!children_fe_indices.empty(), ExcInternalError());

          // convert vector to set
          const std::set<unsigned int> children_fe_indices_set(
            children_fe_indices.begin(), children_fe_indices.end());

          const unsigned int dominated_fe_index =
            fe_collection.find_dominated_fe_extended(children_fe_indices_set,
                                                     /*codim=*/0);

          Assert(dominated_fe_index != numbers::invalid_unsigned_int,
                 typename dealii::hp::FECollection<
                   dim>::ExcNoDominatedFiniteElementAmongstChildren());

          return dominated_fe_index;
        }
      };
    } // namespace DoFHandlerImplementation
  }   // namespace hp
} // namespace internal



template <int dim, int spacedim>
DoFHandler<dim, spacedim>::DoFHandler(const bool is_hp_dof_handler)
  : is_hp_dof_handler(is_hp_dof_handler)
  , tria(nullptr, typeid(*this).name())
  , mg_faces(nullptr)
{}



template <int dim, int spacedim>
DoFHandler<dim, spacedim>::DoFHandler(const Triangulation<dim, spacedim> &tria,
                                      const bool is_hp_dof_handler)
  : is_hp_dof_handler(is_hp_dof_handler)
  , tria(&tria, typeid(*this).name())
  , mg_faces(nullptr)
{
  if (is_hp_dof_handler)
    {
      this->setup_policy_and_listeners();
      this->create_active_fe_table();
    }
  else
    {
      this->setup_policy();
    }
}

template <int dim, int spacedim>
DoFHandler<dim, spacedim>::~DoFHandler()
{
  if (is_hp_dof_handler)
    {
      // unsubscribe as a listener to refinement of the underlying
      // triangulation
      for (auto &connection : this->tria_listeners)
        connection.disconnect();
      this->tria_listeners.clear();

      // ...and release allocated memory
      // virtual functions called in constructors and destructors never use the
      // override in a derived class
      // for clarity be explicit on which function is called
      DoFHandler<dim, spacedim>::clear();
    }
  else
    {
      // release allocated memory
      // virtual functions called in constructors and destructors never use the
      // override in a derived class
      // for clarity be explicit on which function is called
      DoFHandler<dim, spacedim>::clear();

      // also release the policy. this needs to happen before the
      // current object disappears because the policy objects
      // store references to the DoFhandler object they work on
      this->policy.reset();
    }
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::initialize(const Triangulation<dim, spacedim> &tria,
                                      const FiniteElement<dim, spacedim> &fe)
{
  this->initialize(tria, hp::FECollection<dim, spacedim>(fe));
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::initialize(const Triangulation<dim, spacedim> &tria,
                                      const hp::FECollection<dim, spacedim> &fe)
{
  if (is_hp_dof_handler)
    {
      this->clear();

      if (this->tria != &tria)
        {
          for (auto &connection : this->tria_listeners)
            connection.disconnect();
          this->tria_listeners.clear();

          this->tria = &tria;

          this->setup_policy_and_listeners();
        }

      this->create_active_fe_table();

      this->distribute_dofs(fe);
    }
  else
    {
      this->tria = &tria;
      // this->faces                      = nullptr;
      this->number_cache.n_global_dofs = 0;

      this->setup_policy();

      this->distribute_dofs(fe);
    }
}



/*------------------------ Cell iterator functions ------------------------*/

template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::cell_iterator
DoFHandler<dim, spacedim>::begin(const unsigned int level) const
{
  typename Triangulation<dim, spacedim>::cell_iterator cell =
    this->get_triangulation().begin(level);
  if (cell == this->get_triangulation().end(level))
    return end(level);
  return cell_iterator(*cell, this);
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::active_cell_iterator
DoFHandler<dim, spacedim>::begin_active(const unsigned int level) const
{
  // level is checked in begin
  cell_iterator i = begin(level);
  if (i.state() != IteratorState::valid)
    return i;
  while (i->has_children())
    if ((++i).state() != IteratorState::valid)
      return i;
  return i;
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::cell_iterator
DoFHandler<dim, spacedim>::end() const
{
  return cell_iterator(&this->get_triangulation(), -1, -1, this);
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::cell_iterator
DoFHandler<dim, spacedim>::end(const unsigned int level) const
{
  typename Triangulation<dim, spacedim>::cell_iterator cell =
    this->get_triangulation().end(level);
  if (cell.state() != IteratorState::valid)
    return end();
  return cell_iterator(*cell, this);
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::active_cell_iterator
DoFHandler<dim, spacedim>::end_active(const unsigned int level) const
{
  typename Triangulation<dim, spacedim>::cell_iterator cell =
    this->get_triangulation().end_active(level);
  if (cell.state() != IteratorState::valid)
    return active_cell_iterator(end());
  return active_cell_iterator(*cell, this);
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::level_cell_iterator
DoFHandler<dim, spacedim>::begin_mg(const unsigned int level) const
{
  // Assert(this->has_level_dofs(), ExcMessage("You can only iterate over mg "
  //     "levels if mg dofs got distributed."));
  typename Triangulation<dim, spacedim>::cell_iterator cell =
    this->get_triangulation().begin(level);
  if (cell == this->get_triangulation().end(level))
    return end_mg(level);
  return level_cell_iterator(*cell, this);
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::level_cell_iterator
DoFHandler<dim, spacedim>::end_mg(const unsigned int level) const
{
  // Assert(this->has_level_dofs(), ExcMessage("You can only iterate over mg "
  //     "levels if mg dofs got distributed."));
  typename Triangulation<dim, spacedim>::cell_iterator cell =
    this->get_triangulation().end(level);
  if (cell.state() != IteratorState::valid)
    return end();
  return level_cell_iterator(*cell, this);
}



template <int dim, int spacedim>
typename DoFHandler<dim, spacedim>::level_cell_iterator
DoFHandler<dim, spacedim>::end_mg() const
{
  return level_cell_iterator(&this->get_triangulation(), -1, -1, this);
}



template <int dim, int spacedim>
IteratorRange<typename DoFHandler<dim, spacedim>::cell_iterator>
DoFHandler<dim, spacedim>::cell_iterators() const
{
  return IteratorRange<typename DoFHandler<dim, spacedim>::cell_iterator>(
    begin(), end());
}



template <int dim, int spacedim>
IteratorRange<typename DoFHandler<dim, spacedim>::active_cell_iterator>
DoFHandler<dim, spacedim>::active_cell_iterators() const
{
  return IteratorRange<
    typename DoFHandler<dim, spacedim>::active_cell_iterator>(begin_active(),
                                                              end());
}



template <int dim, int spacedim>
IteratorRange<typename DoFHandler<dim, spacedim>::level_cell_iterator>
DoFHandler<dim, spacedim>::mg_cell_iterators() const
{
  return IteratorRange<typename DoFHandler<dim, spacedim>::level_cell_iterator>(
    begin_mg(), end_mg());
}



template <int dim, int spacedim>
IteratorRange<typename DoFHandler<dim, spacedim>::cell_iterator>
DoFHandler<dim, spacedim>::cell_iterators_on_level(
  const unsigned int level) const
{
  return IteratorRange<typename DoFHandler<dim, spacedim>::cell_iterator>(
    begin(level), end(level));
}



template <int dim, int spacedim>
IteratorRange<typename DoFHandler<dim, spacedim>::active_cell_iterator>
DoFHandler<dim, spacedim>::active_cell_iterators_on_level(
  const unsigned int level) const
{
  return IteratorRange<
    typename DoFHandler<dim, spacedim>::active_cell_iterator>(
    begin_active(level), end_active(level));
}



template <int dim, int spacedim>
IteratorRange<typename DoFHandler<dim, spacedim>::level_cell_iterator>
DoFHandler<dim, spacedim>::mg_cell_iterators_on_level(
  const unsigned int level) const
{
  return IteratorRange<typename DoFHandler<dim, spacedim>::level_cell_iterator>(
    begin_mg(level), end_mg(level));
}



//---------------------------------------------------------------------------



template <int dim, int spacedim>
types::global_dof_index
DoFHandler<dim, spacedim>::n_boundary_dofs() const
{
  Assert(!(dim == 2 && spacedim == 3) || is_hp_dof_handler == false,
         ExcNotImplemented());

  Assert(this->fe_collection.size() > 0, ExcNoFESelected());

  std::unordered_set<types::global_dof_index> boundary_dofs;
  std::vector<types::global_dof_index>        dofs_on_face;
  dofs_on_face.reserve(this->get_fe_collection().max_dofs_per_face());

  const IndexSet &owned_dofs = locally_owned_dofs();

  // loop over all faces to check whether they are at a
  // boundary. note that we need not take special care of single
  // lines in 3d (using @p{cell->has_boundary_lines}), since we do
  // not support boundaries of dimension dim-2, and so every
  // boundary line is also part of a boundary face.
  for (const auto &cell : this->active_cell_iterators())
    if (cell->is_locally_owned() && cell->at_boundary())
      {
        for (auto f : GeometryInfo<dim>::face_indices())
          if (cell->at_boundary(f))
            {
              const unsigned int dofs_per_face = cell->get_fe().dofs_per_face;
              dofs_on_face.resize(dofs_per_face);

              cell->face(f)->get_dof_indices(dofs_on_face,
                                             cell->active_fe_index());
              for (unsigned int i = 0; i < dofs_per_face; ++i)
                {
                  const unsigned int global_idof_index = dofs_on_face[i];
                  if (owned_dofs.is_element(global_idof_index))
                    {
                      boundary_dofs.insert(global_idof_index);
                    }
                }
            }
      }
  return boundary_dofs.size();
}



template <int dim, int spacedim>
types::global_dof_index
DoFHandler<dim, spacedim>::n_boundary_dofs(
  const std::set<types::boundary_id> &boundary_ids) const
{
  Assert(!(dim == 2 && spacedim == 3) || is_hp_dof_handler == false,
         ExcNotImplemented());

  Assert(this->fe_collection.size() > 0, ExcNoFESelected());
  Assert(boundary_ids.find(numbers::internal_face_boundary_id) ==
           boundary_ids.end(),
         ExcInvalidBoundaryIndicator());

  // same as above, but with additional checks for set of boundary
  // indicators
  std::unordered_set<types::global_dof_index> boundary_dofs;
  std::vector<types::global_dof_index>        dofs_on_face;
  dofs_on_face.reserve(this->get_fe_collection().max_dofs_per_face());

  const IndexSet &owned_dofs = locally_owned_dofs();

  for (const auto &cell : this->active_cell_iterators())
    if (cell->is_locally_owned() && cell->at_boundary())
      {
        for (auto f : GeometryInfo<dim>::face_indices())
          if (cell->at_boundary(f) &&
              (boundary_ids.find(cell->face(f)->boundary_id()) !=
               boundary_ids.end()))
            {
              const unsigned int dofs_per_face = cell->get_fe().dofs_per_face;
              dofs_on_face.resize(dofs_per_face);

              cell->face(f)->get_dof_indices(dofs_on_face,
                                             cell->active_fe_index());
              for (unsigned int i = 0; i < dofs_per_face; ++i)
                {
                  const unsigned int global_idof_index = dofs_on_face[i];
                  if (owned_dofs.is_element(global_idof_index))
                    {
                      boundary_dofs.insert(global_idof_index);
                    }
                }
            }
      }
  return boundary_dofs.size();
}



template <int dim, int spacedim>
std::size_t
DoFHandler<dim, spacedim>::memory_consumption() const
{
  if (is_hp_dof_handler)
    {
      std::size_t mem =
        (MemoryConsumption::memory_consumption(this->tria) +
         MemoryConsumption::memory_consumption(this->fe_collection) +
         MemoryConsumption::memory_consumption(this->tria) +
         //         MemoryConsumption::memory_consumption(this->levels_hp) +
         // MemoryConsumption::memory_consumption(*this->faces_hp) +
         MemoryConsumption::memory_consumption(this->number_cache) // +
         // MemoryConsumption::memory_consumption(this->vertex_dofs) +
         // MemoryConsumption::memory_consumption(this->vertex_dof_offsets)
        );
      //      for (unsigned int i = 0; i < this->levels_hp.size(); ++i)
      //        mem +=
      //        MemoryConsumption::memory_consumption(*this->levels_hp[i]);
      // mem += MemoryConsumption::memory_consumption(*this->faces_hp);

      return mem;
    }
  else
    {
      std::size_t mem =
        (MemoryConsumption::memory_consumption(this->tria) +
         MemoryConsumption::memory_consumption(this->fe_collection) +
         MemoryConsumption::memory_consumption(this->block_info_object) +
         // MemoryConsumption::memory_consumption(this->levels) +
         // MemoryConsumption::memory_consumption(*this->faces) +
         // MemoryConsumption::memory_consumption(this->faces) +
         sizeof(this->number_cache) +
         MemoryConsumption::memory_consumption(this->n_dofs())
         // + MemoryConsumption::memory_consumption(this->vertex_dofs)
        );
      // for (unsigned int i = 0; i < this->levels.size(); ++i)
      //  mem += MemoryConsumption::memory_consumption(*this->levels[i]);

      for (unsigned int level = 0; level < this->mg_levels.size(); ++level)
        mem += this->mg_levels[level]->memory_consumption();

      if (this->mg_faces != nullptr)
        mem += MemoryConsumption::memory_consumption(*this->mg_faces);

      for (unsigned int i = 0; i < this->mg_vertex_dofs.size(); ++i)
        mem += sizeof(MGVertexDoFs) +
               (1 + this->mg_vertex_dofs[i].get_finest_level() -
                this->mg_vertex_dofs[i].get_coarsest_level()) *
                 sizeof(types::global_dof_index);

      return mem;
    }
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::set_fe(const FiniteElement<dim, spacedim> &fe)
{
  this->set_fe(hp::FECollection<dim, spacedim>(fe));
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::set_fe(const hp::FECollection<dim, spacedim> &ff)
{
  Assert(
    this->tria != nullptr,
    ExcMessage(
      "You need to set the Triangulation in the DoFHandler using initialize() or "
      "in the constructor before you can distribute DoFs."));
  Assert(this->tria->n_levels() > 0,
         ExcMessage("The Triangulation you are using is empty!"));
  Assert(ff.size() > 0, ExcMessage("The hp::FECollection given is empty!"));

  // don't create a new object if the one we have is already appropriate
  if (this->fe_collection != ff)
    this->fe_collection = hp::FECollection<dim, spacedim>(ff);

  if (is_hp_dof_handler)
    {
      // ensure that the active_fe_indices vectors are initialized correctly
      this->create_active_fe_table();

      // make sure every processor knows the active_fe_indices
      // on both its own cells and all ghost cells
      dealii::internal::hp::DoFHandlerImplementation::Implementation::
        communicate_active_fe_indices(*this);

      // make sure that the fe collection is large enough to
      // cover all fe indices presently in use on the mesh
      for (const auto &cell : this->active_cell_iterators())
        if (!cell->is_artificial())
          Assert(cell->active_fe_index() < this->fe_collection.size(),
                 ExcInvalidFEIndex(cell->active_fe_index(),
                                   this->fe_collection.size()));
    }
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::distribute_dofs(
  const FiniteElement<dim, spacedim> &fe)
{
  this->distribute_dofs(hp::FECollection<dim, spacedim>(fe));
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::distribute_dofs(
  const hp::FECollection<dim, spacedim> &ff)
{
  if (is_hp_dof_handler)
    {
      new_dofs.resize(this->tria->n_levels());
      new_dofs_ptr.resize(this->tria->n_levels());
      new_cell_dofs_cache.resize(this->tria->n_levels());
      new_cell_dofs_cache_ptr.resize(this->tria->n_levels());
      new_active_fe_indices.resize(this->tria->n_levels());
      new_future_fe_indices.resize(this->tria->n_levels());
      // assign the fe_collection and initialize all active_fe_indices
      this->set_fe(ff);

      // If an underlying shared::Tria allows artificial cells,
      // then save the current set of subdomain ids, and set
      // subdomain ids to the "true" owner of each cell. we later
      // restore these flags
      std::vector<types::subdomain_id>                      saved_subdomain_ids;
      const parallel::shared::Triangulation<dim, spacedim> *shared_tria =
        (dynamic_cast<const parallel::shared::Triangulation<dim, spacedim> *>(
          &this->get_triangulation()));
      if (shared_tria != nullptr && shared_tria->with_artificial_cells())
        {
          saved_subdomain_ids.resize(shared_tria->n_active_cells());

          const std::vector<types::subdomain_id> &true_subdomain_ids =
            shared_tria->get_true_subdomain_ids_of_cells();

          for (const auto &cell : shared_tria->active_cell_iterators())
            {
              const unsigned int index   = cell->active_cell_index();
              saved_subdomain_ids[index] = cell->subdomain_id();
              cell->set_subdomain_id(true_subdomain_ids[index]);
            }
        }

      // then allocate space for all the other tables
      dealii::internal::hp::DoFHandlerImplementation::Implementation::
        reserve_space(*this);

      // now undo the subdomain modification
      if (shared_tria != nullptr && shared_tria->with_artificial_cells())
        for (const auto &cell : shared_tria->active_cell_iterators())
          cell->set_subdomain_id(
            saved_subdomain_ids[cell->active_cell_index()]);


      // Clear user flags because we will need them. But first we save
      // them and make sure that we restore them later such that at the
      // end of this function the Triangulation will be in the same
      // state as it was at the beginning of this function.
      std::vector<bool> user_flags;
      this->tria->save_user_flags(user_flags);
      const_cast<Triangulation<dim, spacedim> &>(*this->tria)
        .clear_user_flags();


      /////////////////////////////////

      // Now for the real work:
      this->number_cache = this->policy->distribute_dofs();

      /////////////////////////////////

      // do some housekeeping: compress indices
      //{
      //  Threads::TaskGroup<> tg;
      //  for (int level = this->levels_hp.size() - 1; level >= 0; --level)
      //    tg += Threads::new_task(
      //      &dealii::internal::hp::DoFLevel::compress_data<dim, spacedim>,
      //      *this->levels_hp[level],
      //      this->fe_collection);
      //  tg.join_all();
      //}

      // finally restore the user flags
      const_cast<Triangulation<dim, spacedim> &>(*this->tria)
        .load_user_flags(user_flags);
    }
  else
    {
      // first, assign the finite_element
      this->set_fe(ff);

      if (auto tria = dynamic_cast<const Tet::Triangulation<dim, spacedim> *>(
            this->tria.operator->()))
        {
          // 1) clear old state
          this->clear_space();
          new_dofs.resize(1);
          new_dofs_ptr.resize(1);

          // 2) get number of dofs per entity of a single FE
          //    (TODO: should be part of Finite Element)
          const auto n_dofs_per_object =
            [](const FiniteElement<dim, spacedim> &fe) {
              std::array<unsigned int, dim + 1> n_dofs_per_object;

              if (dim >= 0)
                n_dofs_per_object[0] = fe.template n_dofs_per_object<0>();
              if (dim >= 1)
                n_dofs_per_object[1] = fe.template n_dofs_per_object<1>();
              if (dim >= 2)
                n_dofs_per_object[2] = fe.template n_dofs_per_object<2>();
              if (dim >= 3)
                n_dofs_per_object[3] = fe.template n_dofs_per_object<3>();

              return n_dofs_per_object;
            }(ff[0]);

          // 3) reserve space for each entity of triangulation
          const auto &entity_table = tria->get_entity_table();

          for (int d = 0; d <= dim; d++)
            {
              new_dofs[0][d].clear();
              new_dofs_ptr[0][d].clear();
              new_dofs_ptr[0][d].push_back(0);

              for (unsigned int i = 0;
                   i < entity_table[d][d == 0 ? 1 : 0].ptr.size() - 1;
                   i++)
                {
                  for (unsigned int j = 0; j < n_dofs_per_object[d]; j++)
                    new_dofs[0][d].push_back(numbers::invalid_dof_index);
                  new_dofs_ptr[0][d].push_back(new_dofs[0][d].size());
                }
            }

          // 4) distribute dofs
          number_cache = policy->distribute_dofs();

          return;
        }

      // delete all levels and set them up newly. note that we still have to
      // allocate space for all degrees of freedom on this mesh (including ghost
      // and cells that are entirely stored on different processors), though we
      // may not assign numbers to some of them (i.e. they will remain at
      // invalid_dof_index). We need to allocate the space because we will want
      // to be able to query the dof_indices on each cell, and simply be told
      // that we don't know them on some cell (i.e. get back invalid_dof_index)
      this->clear_space();
      new_dofs.resize(this->tria->n_levels());
      new_dofs_ptr.resize(this->tria->n_levels());
      new_cell_dofs_cache.resize(this->tria->n_levels());
      new_cell_dofs_cache_ptr.resize(this->tria->n_levels());
      internal::DoFHandlerImplementation::Implementation::reserve_space(*this);

      // hand things off to the policy
      this->number_cache = this->policy->distribute_dofs();

      // initialize the block info object only if this is a sequential
      // triangulation. it doesn't work correctly yet if it is parallel
      if (dynamic_cast<
            const parallel::DistributedTriangulationBase<dim, spacedim> *>(
            &*this->tria) == nullptr)
        this->block_info_object.initialize(*this, false, true);
    }
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::set_entity_dofs(unsigned int d,
                                           unsigned int index,
                                           const types::global_dof_index *&ptr)
{
  AssertIndexRange(d, new_dofs_ptr[0].size());
  AssertIndexRange(index + 1, new_dofs_ptr[0][d].size());

  const unsigned int n_dofs =
    new_dofs_ptr[0][d][index + 1] - new_dofs_ptr[0][d][index];

  std::memcpy(new_dofs[0][d].data() + new_dofs_ptr[0][d][index],
              ptr,
              n_dofs * sizeof(types::global_dof_index));

  ptr += n_dofs;
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::get_entity_dofs(unsigned int              d,
                                           unsigned int              index,
                                           types::global_dof_index *&ptr) const
{
  AssertIndexRange(index + 1, new_dofs_ptr[0][d].size());

  const unsigned int n_dofs =
    new_dofs_ptr[0][d][index + 1] - new_dofs_ptr[0][d][index];

  std::memcpy(ptr,
              new_dofs[0][d].data() + new_dofs_ptr[0][d][index],
              n_dofs * sizeof(types::global_dof_index));

  ptr += n_dofs;
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::distribute_mg_dofs()
{
  AssertThrow(is_hp_dof_handler == false, ExcNotImplemented());

  Assert(
    this->new_dofs.size() > 0,
    ExcMessage(
      "Distribute active DoFs using distribute_dofs() before calling distribute_mg_dofs()."));

  Assert(
    ((this->tria->get_mesh_smoothing() &
      Triangulation<dim, spacedim>::limit_level_difference_at_vertices) !=
     Triangulation<dim, spacedim>::none),
    ExcMessage(
      "The mesh smoothing requirement 'limit_level_difference_at_vertices' has to be set for using multigrid!"));

  this->clear_mg_space();

  internal::DoFHandlerImplementation::Implementation::reserve_space_mg(*this);
  this->mg_number_cache = this->policy->distribute_mg_dofs();

  // initialize the block info object
  // only if this is a sequential
  // triangulation. it doesn't work
  // correctly yet if it is parallel
  if (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &*this->tria) == nullptr)
    this->block_info_object.initialize(*this, true, false);
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::initialize_local_block_info()
{
  AssertThrow(is_hp_dof_handler == false, ExcNotImplemented());

  this->block_info_object.initialize_local(*this);
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::setup_policy()
{
  // decide whether we need a sequential or a parallel distributed policy
  if (dynamic_cast<const dealii::parallel::shared::Triangulation<dim, spacedim>
                     *>(&this->get_triangulation()) != nullptr ||
      (dynamic_cast<const dealii::Tet::Triangulation<dim, spacedim> *>(
         &this->get_triangulation()) != nullptr &&
       dynamic_cast<const dealii::Tet::Triangulation<dim, spacedim> *>(
         &this->get_triangulation())
           ->is_distributed() == false))
    this->policy =
      std_cxx14::make_unique<internal::DoFHandlerImplementation::Policy::
                               ParallelShared<DoFHandler<dim, spacedim>>>(
        *this);
  else if (dynamic_cast<
             const dealii::parallel::DistributedTriangulationBase<dim, spacedim>
               *>(&this->get_triangulation()) != nullptr ||
           dynamic_cast<const dealii::Tet::Triangulation<dim, spacedim> *>(
             &this->get_triangulation()) != nullptr)
    this->policy =
      std_cxx14::make_unique<internal::DoFHandlerImplementation::Policy::
                               ParallelDistributed<DoFHandler<dim, spacedim>>>(
        *this);
  else
    this->policy =
      std_cxx14::make_unique<internal::DoFHandlerImplementation::Policy::
                               Sequential<DoFHandler<dim, spacedim>>>(*this);
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::clear()
{
  if (is_hp_dof_handler)
    {
      // release memory
      this->clear_space();
    }
  else
    {
      // release memory
      this->clear_space();
      this->clear_mg_space();
    }
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::clear_space()
{
  new_cell_dofs_cache.clear();

  new_cell_dofs_cache_ptr.clear();

  new_dofs.clear();

  new_dofs_ptr.clear();

  if (is_hp_dof_handler)
    {
      //      this->levels_hp.clear();
      this->new_active_fe_indices.clear();
      this->new_future_fe_indices.clear();

      // this->faces_hp.reset();

      // this->vertex_dofs        = std::vector<types::global_dof_index>();
      new_dofs.clear();
      // this->vertex_dof_offsets = std::vector<unsigned int>();
    }
  else
    {
      // std::vector<types::global_dof_index> tmp;
      // std::swap(this->vertex_dofs, tmp);

      this->number_cache.clear();
    }
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::clear_mg_space()
{
  this->mg_levels.clear();
  this->mg_faces.reset();

  std::vector<MGVertexDoFs> tmp;

  std::swap(this->mg_vertex_dofs, tmp);

  this->mg_number_cache.clear();
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::renumber_dofs(
  const std::vector<types::global_dof_index> &new_numbers)
{
  if (is_hp_dof_handler)
    {
      Assert(this->new_future_fe_indices.size() > 0,
             ExcMessage(
               "You need to distribute DoFs before you can renumber them."));

      AssertDimension(new_numbers.size(), this->n_locally_owned_dofs());

#ifdef DEBUG
      // assert that the new indices are
      // consecutively numbered if we are
      // working on a single
      // processor. this doesn't need to
      // hold in the case of a parallel
      // mesh since we map the interval
      // [0...n_dofs()) into itself but
      // only globally, not on each
      // processor
      if (this->n_locally_owned_dofs() == this->n_dofs())
        {
          std::vector<types::global_dof_index> tmp(new_numbers);
          std::sort(tmp.begin(), tmp.end());
          std::vector<types::global_dof_index>::const_iterator p = tmp.begin();
          types::global_dof_index                              i = 0;
          for (; p != tmp.end(); ++p, ++i)
            Assert(*p == i, ExcNewNumbersNotConsecutive(i));
        }
      else
        for (const auto new_number : new_numbers)
          Assert(new_number < this->n_dofs(),
                 ExcMessage(
                   "New DoF index is not less than the total number of dofs."));
#endif

      // uncompress the internal storage scheme of dofs on cells so that
      // we can access dofs in turns. uncompress in parallel, starting
      // with the most expensive levels (the highest ones)
      //{
      //  Threads::TaskGroup<> tg;
      //  for (int level = this->levels_hp.size() - 1; level >= 0; --level)
      //    tg += Threads::new_task(
      //      &dealii::internal::hp::DoFLevel::uncompress_data<dim, spacedim>,
      //      *this->levels_hp[level],
      //      this->fe_collection);
      //  tg.join_all();
      //}

      // do the renumbering
      this->number_cache = this->policy->renumber_dofs(new_numbers);

      // now re-compress the dof indices
      //{
      //  Threads::TaskGroup<> tg;
      //  for (int level = this->levels_hp.size() - 1; level >= 0; --level)
      //    tg += Threads::new_task(
      //      &dealii::internal::hp::DoFLevel::compress_data<dim, spacedim>,
      //     *this->levels_hp[level],
      //      this->fe_collection);
      //  tg.join_all();
      //}
    }
  else
    {
      Assert(this->new_dofs.size() > 0,
             ExcMessage(
               "You need to distribute DoFs before you can renumber them."));

#ifdef DEBUG
      if (dynamic_cast<const parallel::shared::Triangulation<dim, spacedim> *>(
            &*this->tria) != nullptr)
        {
          Assert(new_numbers.size() == this->n_dofs() ||
                   new_numbers.size() == this->n_locally_owned_dofs(),
                 ExcMessage("Incorrect size of the input array."));
        }
      else if (dynamic_cast<
                 const parallel::DistributedTriangulationBase<dim, spacedim> *>(
                 &*this->tria) != nullptr)
        {
          AssertDimension(new_numbers.size(), this->n_locally_owned_dofs());
        }
      else
        {
          AssertDimension(new_numbers.size(), this->n_dofs());
        }

      // assert that the new indices are
      // consecutively numbered if we are
      // working on a single
      // processor. this doesn't need to
      // hold in the case of a parallel
      // mesh since we map the interval
      // [0...n_dofs()) into itself but
      // only globally, not on each
      // processor
      if (this->n_locally_owned_dofs() == this->n_dofs())
        {
          std::vector<types::global_dof_index> tmp(new_numbers);
          std::sort(tmp.begin(), tmp.end());
          std::vector<types::global_dof_index>::const_iterator p = tmp.begin();
          types::global_dof_index                              i = 0;
          for (; p != tmp.end(); ++p, ++i)
            Assert(*p == i, ExcNewNumbersNotConsecutive(i));
        }
      else
        for (const auto new_number : new_numbers)
          Assert(new_number < this->n_dofs(),
                 ExcMessage(
                   "New DoF index is not less than the total number of dofs."));
#endif

      this->number_cache = this->policy->renumber_dofs(new_numbers);
    }
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::renumber_dofs(
  const unsigned int                          level,
  const std::vector<types::global_dof_index> &new_numbers)
{
  AssertThrow(is_hp_dof_handler == false, ExcNotImplemented());

  Assert(
    this->mg_levels.size() > 0 && this->new_dofs.size() > 0,
    ExcMessage(
      "You need to distribute active and level DoFs before you can renumber level DoFs."));
  AssertIndexRange(level, this->get_triangulation().n_global_levels());
  AssertDimension(new_numbers.size(),
                  this->locally_owned_mg_dofs(level).n_elements());

#ifdef DEBUG
  // assert that the new indices are consecutively numbered if we are working
  // on a single processor. this doesn't need to hold in the case of a
  // parallel mesh since we map the interval [0...n_dofs(level)) into itself
  // but only globally, not on each processor
  if (this->n_locally_owned_dofs() == this->n_dofs())
    {
      std::vector<types::global_dof_index> tmp(new_numbers);
      std::sort(tmp.begin(), tmp.end());
      std::vector<types::global_dof_index>::const_iterator p = tmp.begin();
      types::global_dof_index                              i = 0;
      for (; p != tmp.end(); ++p, ++i)
        Assert(*p == i, ExcNewNumbersNotConsecutive(i));
    }
  else
    for (const auto new_number : new_numbers)
      Assert(new_number < this->n_dofs(level),
             ExcMessage(
               "New DoF index is not less than the total number of dofs."));
#endif

  this->mg_number_cache[level] =
    this->policy->renumber_mg_dofs(level, new_numbers);
}



template <int dim, int spacedim>
unsigned int
DoFHandler<dim, spacedim>::max_couplings_between_boundary_dofs() const
{
  Assert(this->fe_collection.size() > 0, ExcNoFESelected());

  switch (dim)
    {
      case 1:
        return this->fe_collection.max_dofs_per_vertex();
      case 2:
        return (3 * this->fe_collection.max_dofs_per_vertex() +
                2 * this->fe_collection.max_dofs_per_line());
      case 3:
        // we need to take refinement of one boundary face into
        // consideration here; in fact, this function returns what
        // #max_coupling_between_dofs<2> returns
        //
        // we assume here, that only four faces meet at the boundary;
        // this assumption is not justified and needs to be fixed some
        // time. fortunately, omitting it for now does no harm since
        // the matrix will cry foul if its requirements are not
        // satisfied
        return (19 * this->fe_collection.max_dofs_per_vertex() +
                28 * this->fe_collection.max_dofs_per_line() +
                8 * this->fe_collection.max_dofs_per_quad());
      default:
        Assert(false, ExcNotImplemented());
        return 0;
    }
}



template <int dim, int spacedim>
unsigned int
DoFHandler<dim, spacedim>::max_couplings_between_dofs() const
{
  Assert(this->fe_collection.size() > 0, ExcNoFESelected());
  return internal::DoFHandlerImplementation::Implementation::
    max_couplings_between_dofs(*this);
}



template <int dim, int spacedim>
template <int structdim>
types::global_dof_index
DoFHandler<dim, spacedim>::get_dof_index(const unsigned int obj_level,
                                         const unsigned int obj_index,
                                         const unsigned int fe_index,
                                         const unsigned int local_index) const
{
  if (is_hp_dof_handler)
    {
      Assert(false, ExcNotImplemented());
      return numbers::invalid_dof_index;
    }
  else
    {
      return internal::DoFHandlerImplementation::Implementation::get_dof_index(
        *this,
        this->mg_levels[obj_level],
        this->mg_faces,
        obj_index,
        fe_index,
        local_index,
        std::integral_constant<int, structdim>());
    }
}



template <int dim, int spacedim>
template <int structdim>
void
DoFHandler<dim, spacedim>::set_dof_index(
  const unsigned int            obj_level,
  const unsigned int            obj_index,
  const unsigned int            fe_index,
  const unsigned int            local_index,
  const types::global_dof_index global_index) const
{
  if (is_hp_dof_handler)
    {
      Assert(false, ExcNotImplemented());
      return;
    }
  else
    {
      internal::DoFHandlerImplementation::Implementation::set_dof_index(
        *this,
        this->mg_levels[obj_level],
        this->mg_faces,
        obj_index,
        fe_index,
        local_index,
        global_index,
        std::integral_constant<int, structdim>());
    }
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::set_active_fe_indices(
  const std::vector<unsigned int> &active_fe_indices)
{
  Assert(active_fe_indices.size() == this->get_triangulation().n_active_cells(),
         ExcDimensionMismatch(active_fe_indices.size(),
                              this->get_triangulation().n_active_cells()));

  this->create_active_fe_table();
  // we could set the values directly, since they are stored as
  // protected data of this object, but for simplicity we use the
  // cell-wise access. this way we also have to pass some debug-mode
  // tests which we would have to duplicate ourselves otherwise
  for (const auto &cell : this->active_cell_iterators())
    if (cell->is_locally_owned())
      cell->set_active_fe_index(active_fe_indices[cell->active_cell_index()]);
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::get_active_fe_indices(
  std::vector<unsigned int> &active_fe_indices) const
{
  active_fe_indices.resize(this->get_triangulation().n_active_cells());

  // we could try to extract the values directly, since they are
  // stored as protected data of this object, but for simplicity we
  // use the cell-wise access.
  for (const auto &cell : this->active_cell_iterators())
    if (!cell->is_artificial())
      active_fe_indices[cell->active_cell_index()] = cell->active_fe_index();
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::setup_policy_and_listeners()
{
  // connect functions to signals of the underlying triangulation
  this->tria_listeners.push_back(this->tria->signals.pre_refinement.connect(
    [this]() { this->pre_refinement_action(); }));
  this->tria_listeners.push_back(this->tria->signals.post_refinement.connect(
    [this]() { this->post_refinement_action(); }));
  this->tria_listeners.push_back(this->tria->signals.create.connect(
    [this]() { this->post_refinement_action(); }));

  // decide whether we need a sequential or a parallel shared/distributed
  // policy and attach corresponding callback functions dealing with the
  // transfer of active_fe_indices
  if (dynamic_cast<
        const dealii::parallel::distributed::Triangulation<dim, spacedim> *>(
        &this->get_triangulation()))
    {
      this->policy = std_cxx14::make_unique<
        internal::DoFHandlerImplementation::Policy::ParallelDistributed<
          DoFHandler<dim, spacedim>>>(*this);

      // repartitioning signals
      this->tria_listeners.push_back(
        this->tria->signals.pre_distributed_repartition.connect([this]() {
          internal::hp::DoFHandlerImplementation::Implementation::
            ensure_absence_of_future_fe_indices<dim, spacedim>(*this);
        }));
      this->tria_listeners.push_back(
        this->tria->signals.pre_distributed_repartition.connect(
          [this]() { this->pre_distributed_active_fe_index_transfer(); }));
      this->tria_listeners.push_back(
        this->tria->signals.post_distributed_repartition.connect(
          [this] { this->post_distributed_active_fe_index_transfer(); }));

      // refinement signals
      this->tria_listeners.push_back(
        this->tria->signals.pre_distributed_refinement.connect(
          [this]() { this->pre_distributed_active_fe_index_transfer(); }));
      this->tria_listeners.push_back(
        this->tria->signals.post_distributed_refinement.connect(
          [this]() { this->post_distributed_active_fe_index_transfer(); }));

      // serialization signals
      this->tria_listeners.push_back(
        this->tria->signals.post_distributed_save.connect([this]() {
          this->post_distributed_serialization_of_active_fe_indices();
        }));
    }
  else if (dynamic_cast<
             const dealii::parallel::shared::Triangulation<dim, spacedim> *>(
             &this->get_triangulation()) != nullptr)
    {
      this->policy =
        std_cxx14::make_unique<internal::DoFHandlerImplementation::Policy::
                                 ParallelShared<DoFHandler<dim, spacedim>>>(
          *this);

      // partitioning signals
      this->tria_listeners.push_back(
        this->tria->signals.pre_partition.connect([this]() {
          internal::hp::DoFHandlerImplementation::Implementation::
            ensure_absence_of_future_fe_indices(*this);
        }));

      // refinement signals
      this->tria_listeners.push_back(this->tria->signals.pre_refinement.connect(
        [this] { this->pre_active_fe_index_transfer(); }));
      this->tria_listeners.push_back(
        this->tria->signals.post_refinement.connect(
          [this] { this->post_active_fe_index_transfer(); }));
    }
  else
    {
      this->policy =
        std_cxx14::make_unique<internal::DoFHandlerImplementation::Policy::
                                 Sequential<DoFHandler<dim, spacedim>>>(*this);

      // refinement signals
      this->tria_listeners.push_back(this->tria->signals.pre_refinement.connect(
        [this] { this->pre_active_fe_index_transfer(); }));
      this->tria_listeners.push_back(
        this->tria->signals.post_refinement.connect(
          [this] { this->post_active_fe_index_transfer(); }));
    }
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::create_active_fe_table()
{
  AssertThrow(is_hp_dof_handler == true, ExcNotImplemented());


  // Create sufficiently many hp::DoFLevels.
  //  while (this->levels_hp.size() < this->tria->n_levels())
  //    this->levels_hp.emplace_back(new dealii::internal::hp::DoFLevel);

  while (this->new_active_fe_indices.size() < this->tria->n_levels())
    this->new_active_fe_indices.push_back({});

  while (this->new_future_fe_indices.size() < this->tria->n_levels())
    this->new_future_fe_indices.push_back({});

  // then make sure that on each level we have the appropriate size
  // of active_fe_indices; preset them to zero, i.e. the default FE
  for (unsigned int level = 0; level < this->new_future_fe_indices.size();
       ++level)
    {
      if (this->new_active_fe_indices[level].size() == 0 &&
          this->new_future_fe_indices[level].size() == 0)
        {
          this->new_active_fe_indices[level].resize(
            this->tria->n_raw_cells(level), 0);
          this->new_future_fe_indices[level].resize(
            this->tria->n_raw_cells(level), invalid_active_fe_index);
        }
      else
        {
          // Either the active_fe_indices have size zero because
          // they were just created, or the correct size. Other
          // sizes indicate that something went wrong.
          Assert(this->new_active_fe_indices[level].size() ==
                     this->tria->n_raw_cells(level) &&
                   this->new_future_fe_indices[level].size() ==
                     this->tria->n_raw_cells(level),
                 ExcInternalError());
        }

      // it may be that the previous table was compressed; in that
      // case, restore the correct active_fe_index. the fact that
      // this no longer matches the indices in the table is of no
      // importance because the current function is called at a
      // point where we have to recreate the dof_indices tables in
      // the levels anyway
      // this->levels_hp[level]->normalize_active_fe_indices();
    }
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::pre_refinement_action()
{
  create_active_fe_table();
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::post_refinement_action()
{
  //  // Normally only one level is added, but if this Triangulation
  //  // is created by copy_triangulation, it can be more than one level.
  //  while (this->levels_hp.size() < this->tria->n_levels())
  //    this->levels_hp.emplace_back(new dealii::internal::hp::DoFLevel);
  //
  //  // Coarsening can lead to the loss of levels. Hence remove them.
  //  while (this->levels_hp.size() > this->tria->n_levels())
  //    {
  //      // drop the last element. that also releases the memory pointed to
  //      this->levels_hp.pop_back();
  //    }

  while (this->new_active_fe_indices.size() < this->tria->n_levels())
    this->new_active_fe_indices.push_back({});

  while (this->new_active_fe_indices.size() > this->tria->n_levels())
    this->new_active_fe_indices.pop_back();

  while (this->new_future_fe_indices.size() < this->tria->n_levels())
    this->new_future_fe_indices.push_back({});

  while (this->new_future_fe_indices.size() > this->tria->n_levels())
    this->new_future_fe_indices.pop_back();



  Assert(this->new_future_fe_indices.size() == this->tria->n_levels(),
         ExcInternalError());
  for (unsigned int i = 0; i < this->new_future_fe_indices.size(); ++i)
    {
      // Resize active_fe_indices vectors. Use zero indicator to extend.
      this->new_active_fe_indices[i].resize(this->tria->n_raw_cells(i), 0);

      // Resize future_fe_indices vectors. Make sure that all
      // future_fe_indices have been cleared after refinement happened.
      //
      // We have used future_fe_indices to update all active_fe_indices
      // before refinement happened, thus we are safe to clear them now.
      this->new_future_fe_indices[i].assign(this->tria->n_raw_cells(i),
                                            invalid_active_fe_index);
    }
}


template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::pre_active_fe_index_transfer()
{
  // Finite elements need to be assigned to each cell by calling
  // distribute_dofs() first to make this functionality available.
  if (this->fe_collection.size() > 0)
    {
      Assert(this->active_fe_index_transfer == nullptr, ExcInternalError());

      this->active_fe_index_transfer =
        std_cxx14::make_unique<ActiveFEIndexTransfer>();

      dealii::internal::hp::DoFHandlerImplementation::Implementation::
        collect_fe_indices_on_cells_to_be_refined(*this);
    }
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::pre_distributed_active_fe_index_transfer()
{
#ifndef DEAL_II_WITH_P4EST
  Assert(false, ExcInternalError());
#else
  // Finite elements need to be assigned to each cell by calling
  // distribute_dofs() first to make this functionality available.
  if (this->fe_collection.size() > 0)
    {
      Assert(this->active_fe_index_transfer == nullptr, ExcInternalError());

      this->active_fe_index_transfer =
        std_cxx14::make_unique<ActiveFEIndexTransfer>();

      // If we work on a p::d::Triangulation, we have to transfer all
      // active_fe_indices since ownership of cells may change. We will
      // use our p::d::CellDataTransfer member to achieve this. Further,
      // we prepare the values in such a way that they will correspond to
      // the active_fe_indices on the new mesh.

      // Gather all current future_fe_indices.
      this->active_fe_index_transfer->active_fe_indices.resize(
        this->get_triangulation().n_active_cells(),
        numbers::invalid_unsigned_int);

      for (const auto &cell : this->active_cell_iterators())
        if (cell->is_locally_owned())
          this->active_fe_index_transfer
            ->active_fe_indices[cell->active_cell_index()] =
            cell->future_fe_index();

      // Create transfer object and attach to it.
      const auto *distributed_tria = dynamic_cast<
        const parallel::distributed::Triangulation<dim, spacedim> *>(
        &this->get_triangulation());

      this->active_fe_index_transfer->cell_data_transfer =
        std_cxx14::make_unique<
          parallel::distributed::
            CellDataTransfer<dim, spacedim, std::vector<unsigned int>>>(
          *distributed_tria,
          /*transfer_variable_size_data=*/false,
          [this](const std::vector<unsigned int> &children_fe_indices) {
            return dealii::internal::hp::DoFHandlerImplementation::
              Implementation::determine_fe_from_children<dim, spacedim>(
                children_fe_indices, this->fe_collection);
          });

      this->active_fe_index_transfer->cell_data_transfer
        ->prepare_for_coarsening_and_refinement(
          this->active_fe_index_transfer->active_fe_indices);
    }
#endif
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::post_active_fe_index_transfer()
{
  // Finite elements need to be assigned to each cell by calling
  // distribute_dofs() first to make this functionality available.
  if (this->fe_collection.size() > 0)
    {
      Assert(this->active_fe_index_transfer != nullptr, ExcInternalError());

      dealii::internal::hp::DoFHandlerImplementation::Implementation::
        distribute_fe_indices_on_refined_cells(*this);

      // We have to distribute the information about active_fe_indices
      // of all cells (including the artificial ones) on all processors,
      // if a parallel::shared::Triangulation has been used.
      dealii::internal::hp::DoFHandlerImplementation::Implementation::
        communicate_active_fe_indices(*this);

      // Free memory.
      this->active_fe_index_transfer.reset();
    }
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::post_distributed_active_fe_index_transfer()
{
#ifndef DEAL_II_WITH_P4EST
  Assert(false, ExcInternalError());
#else
  // Finite elements need to be assigned to each cell by calling
  // distribute_dofs() first to make this functionality available.
  if (this->fe_collection.size() > 0)
    {
      Assert(this->active_fe_index_transfer != nullptr, ExcInternalError());

      // Unpack active_fe_indices.
      this->active_fe_index_transfer->active_fe_indices.resize(
        this->get_triangulation().n_active_cells(),
        numbers::invalid_unsigned_int);
      this->active_fe_index_transfer->cell_data_transfer->unpack(
        this->active_fe_index_transfer->active_fe_indices);

      // Update all locally owned active_fe_indices.
      this->set_active_fe_indices(
        this->active_fe_index_transfer->active_fe_indices);

      // Update active_fe_indices on ghost cells.
      dealii::internal::hp::DoFHandlerImplementation::Implementation::
        communicate_active_fe_indices(*this);

      // Free memory.
      this->active_fe_index_transfer.reset();
    }
#endif
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::prepare_for_serialization_of_active_fe_indices()
{
#ifndef DEAL_II_WITH_P4EST
  Assert(false,
         ExcMessage(
           "You are attempting to use a functionality that is only available "
           "if deal.II was configured to use p4est, but cmake did not find a "
           "valid p4est library."));
#else
  Assert(
    (dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim> *>(
       &this->get_triangulation()) != nullptr),
    ExcMessage(
      "This functionality requires a parallel::distributed::Triangulation object."));

  // Finite elements need to be assigned to each cell by calling
  // distribute_dofs() first to make this functionality available.
  if (this->fe_collection.size() > 0)
    {
      Assert(active_fe_index_transfer == nullptr, ExcInternalError());

      active_fe_index_transfer =
        std_cxx14::make_unique<ActiveFEIndexTransfer>();

      // Create transfer object and attach to it.
      const auto *distributed_tria = dynamic_cast<
        const parallel::distributed::Triangulation<dim, spacedim> *>(
        &this->get_triangulation());

      active_fe_index_transfer->cell_data_transfer = std_cxx14::make_unique<
        parallel::distributed::
          CellDataTransfer<dim, spacedim, std::vector<unsigned int>>>(
        *distributed_tria,
        /*transfer_variable_size_data=*/false,
        [this](const std::vector<unsigned int> &children_fe_indices) {
          return dealii::internal::hp::DoFHandlerImplementation::
            Implementation::determine_fe_from_children<dim, spacedim>(
              children_fe_indices, this->fe_collection);
        });

      // If we work on a p::d::Triangulation, we have to transfer all
      // active fe indices since ownership of cells may change.

      // Gather all current active_fe_indices
      this->get_active_fe_indices(active_fe_index_transfer->active_fe_indices);

      // Attach to transfer object
      active_fe_index_transfer->cell_data_transfer->prepare_for_serialization(
        active_fe_index_transfer->active_fe_indices);
    }
#endif
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::post_distributed_serialization_of_active_fe_indices()
{
#ifndef DEAL_II_WITH_P4EST
  Assert(false,
         ExcMessage(
           "You are attempting to use a functionality that is only available "
           "if deal.II was configured to use p4est, but cmake did not find a "
           "valid p4est library."));
#else
  if (this->fe_collection.size() > 0)
    {
      Assert(this->active_fe_index_transfer != nullptr, ExcInternalError());

      // Free memory.
      this->active_fe_index_transfer.reset();
    }
#endif
}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::deserialize_active_fe_indices()

{
#ifndef DEAL_II_WITH_P4EST
  Assert(false,
         ExcMessage(
           "You are attempting to use a functionality that is only available "
           "if deal.II was configured to use p4est, but cmake did not find a "
           "valid p4est library."));
#else
  Assert(
    (dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim> *>(
       &this->get_triangulation()) != nullptr),
    ExcMessage(
      "This functionality requires a parallel::distributed::Triangulation object."));

  // Finite elements need to be assigned to each cell by calling
  // distribute_dofs() first to make this functionality available.
  if (this->fe_collection.size() > 0)
    {
      Assert(active_fe_index_transfer == nullptr, ExcInternalError());

      active_fe_index_transfer =
        std_cxx14::make_unique<ActiveFEIndexTransfer>();

      // Create transfer object and attach to it.
      const auto *distributed_tria = dynamic_cast<
        const parallel::distributed::Triangulation<dim, spacedim> *>(
        &this->get_triangulation());

      active_fe_index_transfer->cell_data_transfer = std_cxx14::make_unique<
        parallel::distributed::
          CellDataTransfer<dim, spacedim, std::vector<unsigned int>>>(
        *distributed_tria,
        /*transfer_variable_size_data=*/false,
        [this](const std::vector<unsigned int> &children_fe_indices) {
          return dealii::internal::hp::DoFHandlerImplementation::
            Implementation::determine_fe_from_children<dim, spacedim>(
              children_fe_indices, this->fe_collection);
        });

      // Unpack active_fe_indices.
      active_fe_index_transfer->active_fe_indices.resize(
        this->get_triangulation().n_active_cells(),
        numbers::invalid_unsigned_int);
      active_fe_index_transfer->cell_data_transfer->deserialize(
        active_fe_index_transfer->active_fe_indices);

      // Update all locally owned active_fe_indices.
      this->set_active_fe_indices(active_fe_index_transfer->active_fe_indices);

      // Update active_fe_indices on ghost cells.
      dealii::internal::hp::DoFHandlerImplementation::Implementation::
        communicate_active_fe_indices(*this);

      // Free memory.
      active_fe_index_transfer.reset();
    }
#endif
}



template <int dim, int spacedim>
DoFHandler<dim, spacedim>::MGVertexDoFs::MGVertexDoFs()
  : coarsest_level(numbers::invalid_unsigned_int)
  , finest_level(0)
{}



template <int dim, int spacedim>
void
DoFHandler<dim, spacedim>::MGVertexDoFs::init(
  const unsigned int cl,
  const unsigned int fl,
  const unsigned int dofs_per_vertex)
{
  coarsest_level = cl;
  finest_level   = fl;

  if (coarsest_level <= finest_level)
    {
      const unsigned int n_levels  = finest_level - coarsest_level + 1;
      const unsigned int n_indices = n_levels * dofs_per_vertex;

      indices = std_cxx14::make_unique<types::global_dof_index[]>(n_indices);
      std::fill(indices.get(),
                indices.get() + n_indices,
                numbers::invalid_dof_index);
    }
  else
    indices.reset();
}



template <int dim, int spacedim>
unsigned int
DoFHandler<dim, spacedim>::MGVertexDoFs::get_coarsest_level() const
{
  return coarsest_level;
}



template <int dim, int spacedim>
unsigned int
DoFHandler<dim, spacedim>::MGVertexDoFs::get_finest_level() const
{
  return finest_level;
}

/*-------------- Explicit Instantiations -------------------------------*/
#include "dof_handler.inst"


DEAL_II_NAMESPACE_CLOSE
