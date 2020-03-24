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

#ifndef dealii_tria_tet_connectivity_h
#define dealii_tria_tet_connectivity_h

#include <deal.II/base/config.h>

#include <deal.II/grid/tria_tet_cell_type.h>
#include <deal.II/grid/tria_tet_crs.h>

DEAL_II_NAMESPACE_OPEN

namespace Tet
{
  template <int dim>
  class Connectivity
  {
  public:
    void
    build(const std::vector<CellTypeEnum> &cell_types,
          const std::vector<unsigned int> &cell_vertices);

    void
    print(std::ostream &out) const;

    // private:
    std::array<std::array<CRS<unsigned int>, dim + 1>, dim + 1> table;
  };

} // namespace Tet

DEAL_II_NAMESPACE_CLOSE

#endif