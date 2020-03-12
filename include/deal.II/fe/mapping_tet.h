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

#ifndef dealii_mapping_tet_h
#define dealii_mapping_tet_h

#include <deal.II/base/config.h>

#include <deal.II/base/std_cxx14/memory.h>

#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/mapping.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim = dim>
class MappingTet : public Mapping<dim, spacedim>
{
public:
  class InternalData : public Mapping<dim, spacedim>::InternalDataBase
  {
  public:
    InternalData(const unsigned int polynomial_degree)
      : polynomial_degree(polynomial_degree)
    {}

    void
    initialize(const UpdateFlags      update_flags,
               const Quadrature<dim> &quadrature,
               const unsigned int     n_original_q_points)
    {
      (void)update_flags;
      (void)quadrature;
      (void)n_original_q_points;
    }

    const unsigned int polynomial_degree;

    mutable std::vector<Point<spacedim>> mapping_support_points;

    mutable typename Triangulation<dim, spacedim>::cell_iterator
      cell_of_current_support_points;
  };

  MappingTet(const unsigned int degree)
    : polynomial_degree(polynomial_degree)
  {}

  const unsigned int polynomial_degree;

  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual std::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
  get_vertices(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual Point<spacedim>
  get_center(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
             const bool map_center_of_reference_cell = true) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual BoundingBox<spacedim>
  get_bounding_box(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual bool
  preserves_vertex_locations() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual Point<spacedim>
  transform_unit_to_real_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<dim> &                                          p) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual Point<dim>
  transform_real_to_unit_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<spacedim> &                                     p) const
  {
    Assert(false, ExcNotImplemented());
  }

  Point<dim - 1>
  project_real_point_to_unit_point_on_face(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const Point<spacedim> &                                     p) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  transform(const ArrayView<const Tensor<1, dim>> &                  input,
            const MappingKind                                        type,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<1, spacedim>> &output) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  transform(const ArrayView<const DerivativeForm<1, dim, spacedim>> &input,
            const MappingKind                                        type,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<2, spacedim>> &output) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  transform(const ArrayView<const Tensor<2, dim>> &                  input,
            const MappingKind                                        type,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<2, spacedim>> &output) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  transform(const ArrayView<const DerivativeForm<2, dim, spacedim>> &input,
            const MappingKind                                        type,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<3, spacedim>> &output) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  transform(const ArrayView<const Tensor<3, dim>> &                  input,
            const MappingKind                                        type,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<3, spacedim>> &output) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const
  {
    return update_flags;
  }

  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_data(const UpdateFlags update_flags, const Quadrature<dim> &q) const
  {
    std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
          data_ptr = std_cxx14::make_unique<InternalData>(polynomial_degree);
    auto &data     = dynamic_cast<InternalData &>(*data_ptr);
    data.initialize(this->requires_update_flags(update_flags), q, q.size());

    return data_ptr;
  }

  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_face_data(const UpdateFlags          update_flags,
                const Quadrature<dim - 1> &quadrature) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_subface_data(const UpdateFlags          update_flags,
                   const Quadrature<dim - 1> &quadrature) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual CellSimilarity::Similarity
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const
  {
    // ensure that the following static_cast is really correct:
    Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
           ExcInternalError());
    const InternalData &data = static_cast<const InternalData &>(internal_data);

    const unsigned int n_q_points = quadrature.size();

    data.mapping_support_points = this->compute_mapping_support_points(cell);
    data.cell_of_current_support_points = cell;

    Assert(false, ExcNotImplemented());
  }

  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    dealii ::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const
  {
    Assert(false, ExcNotImplemented());
  }


  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          subface_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual std::vector<Point<spacedim>>
  compute_mapping_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
  {
    std::vector<Point<spacedim>> a;
    Assert(false, ExcNotImplemented());

    return a;
  }
};

DEAL_II_NAMESPACE_CLOSE

#endif