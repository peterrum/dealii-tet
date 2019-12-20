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

#include <deal.II/grid/tria_tet_cell_type.h>

DEAL_II_NAMESPACE_OPEN


namespace Tet
{
  namespace
  {
    template <int dim>
    static std::array<CellTypeEntities, dim + 1>
    convert_to_crs(
      const std::vector<std::vector<std::vector<unsigned int>>> &entities_in)
    {
      std::array<CellTypeEntities, dim + 1> result;

      AssertDimension(dim, entities_in.size());

      for (int d = 0; d < dim; d++)
        {
          const auto &entities = entities_in[d];

          for (auto entity : entities)
            {
              for (auto vertex : entity)
                result[dim - d].vertices.push_back(vertex);

              result[dim - d].vertices_ptr.push_back(
                result[dim - d].vertices.size());
            }
        }

      return result;
    }
  } // namespace

  template <int dim>
  CellTypeBase<dim>::CellTypeBase(
    const std::string                                          name,
    const std::vector<std::vector<std::vector<unsigned int>>> &entities_in)
    : name(name)
    , entities(convert_to_crs<dim>(entities_in))
    , n_vertices_(entities[dim].vertices.size())
  {}

  template <int dim>
  unsigned int
  CellTypeBase<dim>::n_entities(const unsigned int d) const
  {
    AssertIndexRange(d, dim + 1);
    return entities[d].vertices_ptr.size() - 1;
  }

  template <int dim>
  unsigned int
  CellTypeBase<dim>::n_vertices()
  {
    return n_vertices_;
  }

  template <int dim>
  dealii::ArrayView<const unsigned int>
  CellTypeBase<dim>::vertices_of_entity(const unsigned int d,
                                        const unsigned int e) const
  {
    AssertIndexRange(d, dim + 1);
    return dealii::ArrayView<const unsigned int>(
      entities[d].vertices.data() + entities[d].vertices_ptr[e],
      entities[d].vertices_ptr[e + 1] - entities[d].vertices_ptr[e]);
  }

  template <int dim>
  std::string
  CellTypeBase<dim>::get_name() const
  {
    return name;
  }

  template <int dim>
  struct CellTypeTet : public CellTypeBase<dim>
  {
    CellTypeTet()
      : CellTypeBase<dim>("tet",
                          {
                            {{0, 1, 2}},             // cell
                            {{0, 1}, {1, 2}, {2, 0}} // edges
                          })
    {
      AssertDimension(dim, 2);
    }
  };

  template <int dim>
  struct CellTypeQuad : public CellTypeBase<dim>
  {
    CellTypeQuad()
      : CellTypeBase<dim>("quad",
                          {
                            {{0, 1, 2, 3}},                  // cell
                            {{0, 1}, {1, 2}, {2, 3}, {3, 0}} // edges
                          })
    {
      AssertDimension(dim, 2);
    }
  };

  template <int dim>
  std::shared_ptr<CellTypeBase<dim>>
  CellTypeFactory<dim>::build(const CellTypeEnum type)
  {
    std::shared_ptr<CellTypeBase<dim>> result;

    switch (type)
      {
        // clang-format off
        case CellTypeEnum::tet:  result.reset(new CellTypeTet<dim>());  break;
        case CellTypeEnum::quad: result.reset(new CellTypeQuad<dim>()); break;
        // clang-format on
        default:
          Assert(false,
                 dealii::StandardExceptions::ExcMessage(
                   "Element type is not supported!"));
      }

    return result;
  }


} // namespace Tet

template class Tet::CellTypeBase<1>;
template class Tet::CellTypeBase<2>;
template class Tet::CellTypeBase<3>;

template class Tet::CellTypeFactory<1>;
template class Tet::CellTypeFactory<2>;
template class Tet::CellTypeFactory<3>;

DEAL_II_NAMESPACE_CLOSE
