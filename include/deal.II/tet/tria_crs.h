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

#ifndef dealii_tria_tet_crs_h
#define dealii_tria_tet_crs_h

#include <deal.II/base/config.h>

#include <deal.II/tet/tria_cell_type.h>

DEAL_II_NAMESPACE_OPEN

namespace Tet
{
  template <typename T = unsigned int>
  struct CRS
  {
    std::vector<std::size_t> ptr = {0};
    std::vector<T>           col;
  };

} // namespace Tet

DEAL_II_NAMESPACE_CLOSE

#endif
