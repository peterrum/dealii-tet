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

#ifndef dealii_tet_quadrature_lib_h
#define dealii_tet_quadrature_lib_h


#include <deal.II/base/config.h>

#include <deal.II/base/quadrature.h>

DEAL_II_NAMESPACE_OPEN

namespace Tet
{
  template <int dim>
  class QGauss : public Quadrature<dim>
  {
  public:
    QGauss(const unsigned int n_points);
  };
} // namespace Tet

DEAL_II_NAMESPACE_CLOSE

#endif
