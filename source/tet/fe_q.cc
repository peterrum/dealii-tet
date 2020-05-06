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

#ifndef dealii_fe_dgp_monomial_h
#define dealii_fe_dgp_monomial_h

#include <deal.II/base/config.h>

#include <deal.II/tet/fe_q.h>

DEAL_II_NAMESPACE_OPEN

namespace Tet
{
  template <int dim, int spacedim>
  FE_Q<dim, spacedim>::FE_Q(const unsigned int degree)
    : FE_Poly<Tet::ScalarPolynomial<dim>, dim, spacedim>(
        Tet::ScalarPolynomial<dim>(degree),
        FiniteElementData<dim>(get_dpo_vector(degree),
                               dim == 2 ?
                                 GeometryInfoWrapper<dim>(3, 3, 1, 0, 2, 1, 0) :
                                 GeometryInfoWrapper<dim>(4, 6, 4, 1, 3, 3, 1),
                               1,
                               degree,
                               FiniteElementData<dim>::L2),
        std::vector<bool>(FiniteElementData<dim>(
                            get_dpo_vector(degree),
                            dim == 2 ?
                              GeometryInfoWrapper<dim>(3, 3, 1, 0, 2, 1, 0) :
                              GeometryInfoWrapper<dim>(4, 6, 4, 1, 3, 3, 1),
                            1,
                            degree)
                            .dofs_per_cell,
                          true),
        std::vector<ComponentMask>(
          FiniteElementData<dim>(
            get_dpo_vector(degree),
            dim == 2 ? GeometryInfoWrapper<dim>(3, 3, 1, 0, 2, 1, 0) :
                       GeometryInfoWrapper<dim>(4, 6, 4, 1, 3, 3, 1),
            1,
            degree)
            .dofs_per_cell,
          std::vector<bool>(1, true)))
  {
    this->unit_support_points.clear();

    if (dim == 2)
      {
        if (degree == 1) // DRT::Element::tri3 (TODO: change order)
          {
            this->unit_support_points.emplace_back(1.0, 0.0);
            this->unit_support_points.emplace_back(0.0, 1.0);
            this->unit_support_points.emplace_back(0.0, 0.0);
          }
        else if (degree == 2) // DRT::Element::tri6
          {
            this->unit_support_points.emplace_back(1.0, 0.0);
            this->unit_support_points.emplace_back(0.0, 1.0);
            this->unit_support_points.emplace_back(0.0, 0.0);
            this->unit_support_points.emplace_back(0.5, 0.5);
            this->unit_support_points.emplace_back(0.0, 0.5);
            this->unit_support_points.emplace_back(0.5, 0.0);
          }
        else
          {
            Assert(false, ExcNotImplemented());
          }
      }
    else if (dim == 3)
      {
        if (degree == 1)
          {
            this->unit_support_points.emplace_back(0.0, 0.0, 0.0);
            this->unit_support_points.emplace_back(1.0, 0.0, 0.0);
            this->unit_support_points.emplace_back(0.0, 1.0, 0.0);
            this->unit_support_points.emplace_back(0.0, 0.0, 1.0);
          }
        else if (degree == 2)
          {
            this->unit_support_points.emplace_back(0.0, 0.0, 0.0);
            this->unit_support_points.emplace_back(1.0, 0.0, 0.0);
            this->unit_support_points.emplace_back(0.0, 1.0, 0.0);
            this->unit_support_points.emplace_back(0.0, 0.0, 1.0);

            this->unit_support_points.emplace_back(0.5, 0.0, 0.0);
            this->unit_support_points.emplace_back(0.5, 0.5, 0.0);
            this->unit_support_points.emplace_back(0.0, 0.5, 0.0);
            this->unit_support_points.emplace_back(0.0, 0.0, 0.5);
            this->unit_support_points.emplace_back(0.5, 0.0, 0.5);
            this->unit_support_points.emplace_back(0.0, 0.5, 0.5);
          }
        else
          {
            Assert(false, ExcNotImplemented());
          }
      }
    else
      {
        Assert(false, ExcNotImplemented());
      }
  }

  template <int dim, int spacedim>
  std::unique_ptr<FiniteElement<dim, dim>>
  FE_Q<dim, spacedim>::clone() const
  {
    return std_cxx14::make_unique<FE_Q<dim, spacedim>>(*this);
  }

  template <int dim, int spacedim>
  std::string
  FE_Q<dim, spacedim>::get_name() const
  {
    std::ostringstream namebuf;
    namebuf << "FE_Q<" << dim << ">(" << this->degree << ")";

    return namebuf.str();
  }

} // namespace Tet

// explicit instantiations
#include "fe_q.inst"

DEAL_II_NAMESPACE_CLOSE

#endif
