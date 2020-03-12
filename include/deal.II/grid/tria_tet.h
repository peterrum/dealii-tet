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

#ifndef dealii_tria_tet_h
#define dealii_tria_tet_h


#include <deal.II/base/config.h>

#include <deal.II/grid/tria.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim = dim>
class TetTriangulation : public Triangulation<dim, spacedim>
{
public:
  using cell_iterator =
    typename dealii::Triangulation<dim, spacedim>::cell_iterator;

  using active_cell_iterator =
    typename dealii::Triangulation<dim, spacedim>::active_cell_iterator;

  using CellStatus = typename dealii::Triangulation<dim, spacedim>::CellStatus;

  virtual IteratorRange<cell_iterator>
  cell_iterators() const
  {
    Assert(false, ExcNotImplemented());
  }
};

DEAL_II_NAMESPACE_CLOSE

#endif
