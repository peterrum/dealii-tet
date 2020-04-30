// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2018 by the deal.II authors
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

#ifndef dealii_tria_faces_h
#define dealii_tria_faces_h

#include <deal.II/base/config.h>

#include <deal.II/grid/tria_object.h>
#include <deal.II/grid/tria_objects.h>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace TriangulationImplementation
  {
    /**
     * General template for information belonging to the faces of a
     * triangulation. These classes are similar to the TriaLevel classes. As
     * cells are organised in a hierarchical structure of levels, each
     * triangulation consists of several such TriaLevels. However the faces of
     * a triangulation, lower dimensional objects like lines in 2D or lines
     * and quads in 3D, do not have to be based on such a hierarchical
     * structure. In fact we have to organise them in only one object if we
     * want to enable anisotropic refinement. Therefore the TriaFaces classes
     * store the information belonging to the faces of a triangulation
     * separately from the TriaLevel classes.
     *
     * @author Tobias Leicht, 2006
     */
    class TriaFaces
    {
    public:
      // for Boost::serialization
      TriaFaces() = default;

      TriaFaces(unsigned int dim)
        : dim(dim){};

      unsigned int dim;

      /**
       * The TriaObject containing the data of quads.
       */

      TriaObjects<TriaObject<2>> quads;

      /**
       * In effect, this field has <code>4*n_quads</code> elements, being the
       * number of quads times the four lines each has.
       */
      std::vector<bool> line_orientations;

      /**
       * The TriaObject containing the data of lines.
       */
      TriaObjects<TriaObject<1>> lines;

      /**
       * Assert that enough space is allocated to accommodate
       * <code>new_quads_in_pairs</code> new quads, stored in pairs, plus
       * <code>new_quads_single</code> stored individually. This function does
       * not only call <code>vector::reserve()</code>, but does really append
       * the needed elements.
       */
      void
      reserve_space(const unsigned int new_quads_in_pairs,
                    const unsigned int new_quads_single = 0)
      {
        AssertDimension(this->dim, 3);

        Assert(new_quads_in_pairs % 2 == 0, ExcInternalError());

        unsigned int next_free_single = 0;
        unsigned int next_free_pair   = 0;

        // count the number of objects, of unused single objects and of
        // unused pairs of objects
        unsigned int n_quads          = 0;
        unsigned int n_unused_pairs   = 0;
        unsigned int n_unused_singles = 0;
        for (unsigned int i = 0; i < quads.used.size(); ++i)
          {
            if (quads.used[i])
              ++n_quads;
            else if (i + 1 < quads.used.size())
              {
                if (quads.used[i + 1])
                  {
                    ++n_unused_singles;
                    if (next_free_single == 0)
                      next_free_single = i;
                  }
                else
                  {
                    ++n_unused_pairs;
                    if (next_free_pair == 0)
                      next_free_pair = i;
                    ++i;
                  }
              }
            else
              ++n_unused_singles;
          }
        Assert(n_quads + 2 * n_unused_pairs + n_unused_singles ==
                 quads.used.size(),
               ExcInternalError());

        // how many single quads are needed in addition to n_unused_quads?
        const int additional_single_quads = new_quads_single - n_unused_singles;

        unsigned int new_size =
          quads.used.size() + new_quads_in_pairs - 2 * n_unused_pairs;
        if (additional_single_quads > 0)
          new_size += additional_single_quads;

        // see above...
        if (new_size > quads.cells.size())
          {
            // reserve the field of the derived class
            line_orientations.reserve(new_size *
                                      GeometryInfo<2>::lines_per_cell);
            line_orientations.insert(line_orientations.end(),
                                     new_size *
                                         GeometryInfo<2>::lines_per_cell -
                                       line_orientations.size(),
                                     true);
          }
      }

      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Read or write the data of this object to or from a stream for the
       * purpose of serialization
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };



    template <class Archive>
    void
    TriaFaces::serialize(Archive &ar, const unsigned int)
    {
      ar &dim;

      if (dim == 2)
        ar &lines;

      if (dim == 3)
        ar &quads &lines &line_orientations;
    }
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
