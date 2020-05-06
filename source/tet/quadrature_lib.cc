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

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>

#include <deal.II/tet/quadrature_lib.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>


DEAL_II_NAMESPACE_OPEN

namespace Tet
{
  template <int dim>
  QGauss<dim>::QGauss(const unsigned int n_points)
  {
    const double Q14  = 1.0 / 4.0;
    const double Q16  = 1.0 / 6.0;
    const double Q124 = 1.0 / 6.0 / 4.0;

    if (dim == 2)
      {
        if (n_points == 1) // intrule_tri_1point
          {
            const double p = 1.0 / 3.0;
            this->quadrature_points.emplace_back(p, p);
            this->weights.emplace_back(1.0);
          }
        else if (n_points == 3) // intrule_tri_3point
          {
            const double Q23 = 2.0 / 3.0;
            const double Q16 = 1.0 / 6.0;

            this->quadrature_points.emplace_back(Q23, Q16);
            this->quadrature_points.emplace_back(Q16, Q23);
            this->quadrature_points.emplace_back(Q16, Q16);
            this->weights.emplace_back(Q16);
            this->weights.emplace_back(Q16);
            this->weights.emplace_back(Q16);
          }
        else if (n_points == 7) // intrule_tri_7point
          {
            const double q12 = 0.5;

            // clang-format off
            this->quadrature_points.emplace_back(0.3333333333330, 0.3333333333330);
            this->quadrature_points.emplace_back(0.7974269853530, 0.1012865073230);
            this->quadrature_points.emplace_back(0.1012865073230, 0.7974269853530);
            this->quadrature_points.emplace_back(0.1012865073230, 0.1012865073230);
            this->quadrature_points.emplace_back(0.0597158717898, 0.4701420641050);
            this->quadrature_points.emplace_back(0.4701420641050, 0.0597158717898);
            this->quadrature_points.emplace_back(0.4701420641050, 0.4701420641050);
            // clang-format on

            this->weights.emplace_back(q12 * 0.225);
            this->weights.emplace_back(q12 * 0.125939180545);
            this->weights.emplace_back(q12 * 0.125939180545);
            this->weights.emplace_back(q12 * 0.125939180545);
            this->weights.emplace_back(q12 * 0.132394152789);
            this->weights.emplace_back(q12 * 0.132394152789);
            this->weights.emplace_back(q12 * 0.132394152789);
          }
      }
    else if (dim == 3)
      {
        if (n_points == 1) // intrule_tet_1point
          {
            this->quadrature_points.emplace_back(Q14, Q14, Q14);
            this->weights.emplace_back(Q16);
          }
        else if (n_points == 4) // intrule_tet_4point
          {
            const double palpha = (5.0 + 3.0 * sqrt(5.0)) / 20.0;
            const double pbeta  = (5.0 - sqrt(5.0)) / 20.0;
            this->quadrature_points.emplace_back(pbeta, pbeta, pbeta);
            this->quadrature_points.emplace_back(palpha, pbeta, pbeta);
            this->quadrature_points.emplace_back(pbeta, palpha, pbeta);
            this->quadrature_points.emplace_back(pbeta, pbeta, palpha);
            this->weights.emplace_back(Q124);
            this->weights.emplace_back(Q124);
            this->weights.emplace_back(Q124);
            this->weights.emplace_back(Q124);
          }
        else if (n_points == 10) // intrule_tet_10point
          {
            const double Q16 = 1.0 / 6.0;

            // clang-format off
            this->quadrature_points.emplace_back(0.5684305841968444, 0.1438564719343852, 0.1438564719343852);
            this->quadrature_points.emplace_back(0.1438564719343852, 0.1438564719343852, 0.1438564719343852);
            this->quadrature_points.emplace_back(0.1438564719343852, 0.1438564719343852, 0.5684305841968444);
            this->quadrature_points.emplace_back(0.1438564719343852, 0.5684305841968444, 0.1438564719343852);
            this->quadrature_points.emplace_back(0.0000000000000000, 0.5000000000000000, 0.5000000000000000);
            this->quadrature_points.emplace_back(0.5000000000000000, 0.0000000000000000, 0.5000000000000000);
            this->quadrature_points.emplace_back(0.5000000000000000, 0.5000000000000000, 0.0000000000000000);
            this->quadrature_points.emplace_back(0.5000000000000000, 0.0000000000000000, 0.0000000000000000);
            this->quadrature_points.emplace_back(0.0000000000000000, 0.5000000000000000, 0.0000000000000000);
            this->quadrature_points.emplace_back(0.0000000000000000, 0.0000000000000000, 0.5000000000000000);
            // clang-format on

            this->weights.emplace_back(0.2177650698804054 * Q16);
            this->weights.emplace_back(0.2177650698804054 * Q16);
            this->weights.emplace_back(0.2177650698804054 * Q16);
            this->weights.emplace_back(0.2177650698804054 * Q16);
            this->weights.emplace_back(0.0214899534130631 * Q16);
            this->weights.emplace_back(0.0214899534130631 * Q16);
            this->weights.emplace_back(0.0214899534130631 * Q16);
            this->weights.emplace_back(0.0214899534130631 * Q16);
            this->weights.emplace_back(0.0214899534130631 * Q16);
            this->weights.emplace_back(0.0214899534130631 * Q16);
          }
      }

    AssertDimension(this->quadrature_points.size(), this->weights.size());
    Assert(this->quadrature_points.size() > 0,
           ExcMessage("No valid quadrature points!"));
  }

} // namespace Tet


template class Tet::QGauss<2>;
template class Tet::QGauss<3>;

DEAL_II_NAMESPACE_CLOSE
