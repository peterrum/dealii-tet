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

#ifndef dealii_tria_tet_cell_type_h
#define dealii_tria_tet_cell_type_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>

DEAL_II_NAMESPACE_OPEN

namespace Tet
{
  enum class CellTypeEnum
  {
    quad,
    tet,
  };

  struct CellTypeEntities
  {
    std::vector<unsigned int> vertices     = {};
    std::vector<unsigned int> vertices_ptr = {0};
  };

  template <int dim>
  struct CellTypeBase
  {
    CellTypeBase(
      const std::string                                          name,
      const std::vector<std::vector<std::vector<unsigned int>>> &entities_in);

    unsigned int
    n_entities(const unsigned int d) const;

    unsigned int
    n_vertices();

    dealii::ArrayView<const unsigned int>
    vertices_of_entity(const unsigned int d, const unsigned int e) const;

    std::string
    get_name() const;

  private:
    const std::string name;

    const std::array<CellTypeEntities, dim + 1> entities;

    const unsigned int n_vertices_;
  };

  template <int dim>
  class CellTypeFactory
  {
  public:
    static std::shared_ptr<CellTypeBase<dim>>
    build(const CellTypeEnum type);
  };

} // namespace Tet

DEAL_II_NAMESPACE_CLOSE

#endif