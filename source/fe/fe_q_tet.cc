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

#include <deal.II/fe/fe_q_tet.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
FE_QTet<dim, spacedim>::FE_QTet(const unsigned int degree)
  : FE_Poly<PolynomialsTet<dim>, dim, spacedim>(
      PolynomialsTet<dim>(degree),
      FiniteElementData<dim>(get_dpo_vector(degree),
                             1,
                             degree,
                             FiniteElementData<dim>::L2),
      std::vector<bool>(
        FiniteElementData<dim>(get_dpo_vector(degree), 1, degree).dofs_per_cell,
        true),
      std::vector<ComponentMask>(
        FiniteElementData<dim>(get_dpo_vector(degree), 1, degree).dofs_per_cell,
        std::vector<bool>(1, true)))
{}

template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, dim>>
FE_QTet<dim, spacedim>::clone() const
{
  return std_cxx14::make_unique<FE_QTet<dim, spacedim>>(*this);
}

template <int dim, int spacedim>
std::string
FE_QTet<dim, spacedim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_QTet<" << dim << ">(" << this->degree << ")";

  return namebuf.str();
}


// explicit instantiations
#include "fe_q_tet.inst"

DEAL_II_NAMESPACE_CLOSE

#endif
