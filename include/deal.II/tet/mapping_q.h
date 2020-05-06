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

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/std_cxx14/memory.h>

#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/tet/fe_q.h>

DEAL_II_NAMESPACE_OPEN

namespace Tet
{
  template <int dim, int spacedim = dim>
  class MappingQ : public dealii::Mapping<dim, spacedim>
  {
  public:
    class InternalData : public dealii::Mapping<dim, spacedim>::InternalDataBase
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

        this->update_each = update_flags;

        Tet::FE_Q<dim> fe(polynomial_degree);
        n_shape_functions = fe.dofs_per_cell;

        quadrature_points.resize(quadrature.size());


        // step 1: shape functions
        shape_.clear();
        shape_.reserve(n_shape_functions * quadrature.size());
        shape_derivatives.reserve(n_shape_functions * quadrature.size());

        for (unsigned int j = 0; j < quadrature.size(); j++)
          for (unsigned int i = 0; i < n_shape_functions; i++)
            {
              shape_.push_back(fe.shape_value(i, quadrature.point(j)));
              shape_derivatives.push_back(
                fe.shape_grad(i, quadrature.point(j)));
            }
      }

      const double &
      shape(const unsigned int qpoint, const unsigned int shape_nr) const
      {
        return shape_[qpoint * n_shape_functions + shape_nr];
      }

      const Tensor<1, dim> &
      derivative(const unsigned int qpoint, const unsigned int shape_nr) const
      {
        return shape_derivatives[qpoint * n_shape_functions + shape_nr];
      }


      const unsigned int polynomial_degree;

      unsigned int n_shape_functions;

      mutable std::vector<Point<spacedim>> mapping_support_points;

      mutable typename Triangulation<dim, spacedim>::cell_iterator
        cell_of_current_support_points;

      std::vector<Point<spacedim>> quadrature_points;

      std::vector<double> shape_;

      std::vector<Tensor<1, dim>> shape_derivatives;

      mutable std::vector<DerivativeForm<1, dim, spacedim>> covariant;

      mutable std::vector<DerivativeForm<1, dim, spacedim>> contravariant;
    };

    MappingQ(const unsigned int degree)
      : polynomial_degree(degree)
    {}

    const unsigned int polynomial_degree;

    virtual std::unique_ptr<dealii::Mapping<dim, spacedim>>
    clone() const
    {
      return std_cxx14::make_unique<MappingQ<dim, spacedim>>(
        this->polynomial_degree);
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
    transform(
      const ArrayView<const Tensor<1, dim>> &                          input,
      const MappingKind                                                type,
      const typename dealii::Mapping<dim, spacedim>::InternalDataBase &internal,
      const ArrayView<Tensor<1, spacedim>> &output) const
    {
      AssertDimension(input.size(), output.size());
      Assert(mapping_covariant == type, ExcNotImplemented());

      const typename dealii::Tet::MappingQ<dim, spacedim>::InternalData &data =
        static_cast<
          const typename dealii::Tet::MappingQ<dim, spacedim>::InternalData &>(
          internal);

      for (unsigned int i = 0; i < output.size(); ++i)
        output[i] = apply_transformation(data.covariant[i], input[i]);
    }

    virtual void
    transform(
      const ArrayView<const DerivativeForm<1, dim, spacedim>> &        input,
      const MappingKind                                                type,
      const typename dealii::Mapping<dim, spacedim>::InternalDataBase &internal,
      const ArrayView<Tensor<2, spacedim>> &output) const
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    transform(
      const ArrayView<const Tensor<2, dim>> &                          input,
      const MappingKind                                                type,
      const typename dealii::Mapping<dim, spacedim>::InternalDataBase &internal,
      const ArrayView<Tensor<2, spacedim>> &output) const
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    transform(
      const ArrayView<const DerivativeForm<2, dim, spacedim>> &        input,
      const MappingKind                                                type,
      const typename dealii::Mapping<dim, spacedim>::InternalDataBase &internal,
      const ArrayView<Tensor<3, spacedim>> &output) const
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    transform(
      const ArrayView<const Tensor<3, dim>> &                          input,
      const MappingKind                                                type,
      const typename dealii::Mapping<dim, spacedim>::InternalDataBase &internal,
      const ArrayView<Tensor<3, spacedim>> &output) const
    {
      Assert(false, ExcNotImplemented());
    }

    virtual UpdateFlags
    requires_update_flags(const UpdateFlags in) const
    {
      // add flags if the respective quantities are necessary to compute
      // what we need. note that some flags appear in both the conditions
      // and in subsequent set operations. this leads to some circular
      // logic. the only way to treat this is to iterate. since there are
      // 5 if-clauses in the loop, it will take at most 5 iterations to
      // converge. do them:
      UpdateFlags out = in;
      for (unsigned int i = 0; i < 5; ++i)
        {
          // The following is a little incorrect:
          // If not applied on a face,
          // update_boundary_forms does not
          // make sense. On the other hand,
          // it is necessary on a
          // face. Currently,
          // update_boundary_forms is simply
          // ignored for the interior of a
          // cell.
          if (out & (update_JxW_values | update_normal_vectors))
            out |= update_boundary_forms;

          if (out & (update_covariant_transformation | update_JxW_values |
                     update_jacobians | update_jacobian_grads |
                     update_boundary_forms | update_normal_vectors))
            out |= update_contravariant_transformation;

          if (out &
              (update_inverse_jacobians | update_jacobian_pushed_forward_grads |
               update_jacobian_pushed_forward_2nd_derivatives |
               update_jacobian_pushed_forward_3rd_derivatives))
            out |= update_covariant_transformation;

          // The contravariant transformation is used in the Piola
          // transformation, which requires the determinant of the Jacobi
          // matrix of the transformation.  Because we have no way of
          // knowing here whether the finite element wants to use the
          // contravariant or the Piola transforms, we add the JxW values
          // to the list of flags to be updated for each cell.
          if (out & update_contravariant_transformation)
            out |= update_volume_elements;

          // the same is true when computing normal vectors: they require
          // the determinant of the Jacobian
          if (out & update_normal_vectors)
            out |= update_volume_elements;
        }

      return out;
    }

    virtual std::unique_ptr<
      typename dealii::Mapping<dim, spacedim>::InternalDataBase>
    get_data(const UpdateFlags update_flags, const Quadrature<dim> &q) const
    {
      std::unique_ptr<typename dealii::Mapping<dim, spacedim>::InternalDataBase>
            data_ptr = std_cxx14::make_unique<InternalData>(polynomial_degree);
      auto &data     = dynamic_cast<InternalData &>(*data_ptr);
      data.initialize(this->requires_update_flags(update_flags), q, q.size());

      return data_ptr;
    }

    virtual std::unique_ptr<
      typename dealii::Mapping<dim, spacedim>::InternalDataBase>
    get_face_data(const UpdateFlags          update_flags,
                  const Quadrature<dim - 1> &quadrature) const
    {
      Assert(false, ExcNotImplemented());
    }

    virtual std::unique_ptr<
      typename dealii::Mapping<dim, spacedim>::InternalDataBase>
    get_subface_data(const UpdateFlags          update_flags,
                     const Quadrature<dim - 1> &quadrature) const
    {
      Assert(false, ExcNotImplemented());
    }

    virtual CellSimilarity::Similarity
    fill_fe_values(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      const CellSimilarity::Similarity cell_similarity,
      const Quadrature<dim> &          quadrature,
      const typename dealii::Mapping<dim, spacedim>::InternalDataBase
        &internal_data,
      dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                   spacedim>
        &output_data) const;

    virtual void
    fill_fe_face_values(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      const unsigned int                                          face_no,
      const Quadrature<dim - 1> &                                 quadrature,
      const typename dealii::Mapping<dim, spacedim>::InternalDataBase
        &internal_data,
      dealii ::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                    spacedim>
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
      const typename dealii::Mapping<dim, spacedim>::InternalDataBase
        &internal_data,
      dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                   spacedim>
        &output_data) const
    {
      Assert(false, ExcNotImplemented());
    }

    virtual std::vector<Point<spacedim>>
    compute_mapping_support_points(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
    {
      AssertDimension(polynomial_degree, 1);

      std::vector<Point<spacedim>> a(cell->n_vertices());

      for (unsigned int i = 0; i < a.size(); i++)
        a[i] = cell->vertex(i);

      return a;
    }
  };



  namespace internal
  {
    namespace MappingTet
    {
      template <int dim, int spacedim>
      void
      maybe_compute_q_points(
        const typename dealii::Tet::MappingQ<dim, spacedim>::InternalData &data,
        std::vector<Point<spacedim>> &quadrature_points)
      {
        const UpdateFlags update_flags = data.update_each;

        if (update_flags & update_quadrature_points)
          {
            for (unsigned int point = 0; point < quadrature_points.size();
                 ++point)
              {
                const double *  shape = &data.shape(point, 0);
                Point<spacedim> result =
                  (shape[0] * data.mapping_support_points[0]);
                for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                  for (unsigned int i = 0; i < spacedim; ++i)
                    result[i] += shape[k] * data.mapping_support_points[k][i];
                quadrature_points[point] = result;
              }
          }
      }

      template <int dim, int spacedim>
      void
      maybe_update_Jacobians(
        const typename dealii::Tet::MappingQ<dim, spacedim>::InternalData &data)
      {
        const UpdateFlags update_flags = data.update_each;

        if (update_flags & update_contravariant_transformation)
          {
            const unsigned int n_q_points = data.contravariant.size();

            std::fill(data.contravariant.begin(),
                      data.contravariant.end(),
                      DerivativeForm<1, dim, spacedim>());

            Assert(data.n_shape_functions > 0, ExcInternalError());
            const Tensor<1, spacedim> *supp_pts =
              &data.mapping_support_points[0];

            for (unsigned int point = 0; point < n_q_points; ++point)
              {
                const Tensor<1, dim> *data_derv = &data.derivative(point, 0);

                double result[spacedim][dim];

                // peel away part of sum to avoid zeroing the
                // entries and adding for the first time
                for (unsigned int i = 0; i < spacedim; ++i)
                  for (unsigned int j = 0; j < dim; ++j)
                    result[i][j] = data_derv[0][j] * supp_pts[0][i];
                for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim; ++j)
                      result[i][j] += data_derv[k][j] * supp_pts[k][i];

                // write result into contravariant data. for
                // j=dim in the case dim<spacedim, there will
                // never be any nonzero data that arrives in
                // here, so it is ok anyway because it was
                // initialized to zero at the initialization
                for (unsigned int i = 0; i < spacedim; ++i)
                  for (unsigned int j = 0; j < dim; ++j)
                    data.contravariant[point][i][j] = result[i][j];
              }
          }

        if (update_flags & update_covariant_transformation)
          {
            const unsigned int n_q_points = data.contravariant.size();
            for (unsigned int point = 0; point < n_q_points; ++point)
              data.covariant[point] =
                (data.contravariant[point]).covariant_form();
          }
      }
    } // namespace MappingTet
  }   // namespace internal

  template <int dim, int spacedim>
  CellSimilarity::Similarity
  MappingQ<dim, spacedim>::fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const typename dealii::Mapping<dim, spacedim>::InternalDataBase
      &internal_data,
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

    // TODO
    const CellSimilarity::Similarity computed_cell_similarity =
      CellSimilarity::none;

    internal::MappingTet::maybe_compute_q_points<dim, spacedim>(
      data, output_data.quadrature_points);

    data.contravariant.resize(n_q_points);
    data.covariant.resize(n_q_points);


    internal::MappingTet::maybe_update_Jacobians<dim, spacedim>(data);

    const UpdateFlags          update_flags = data.update_each;
    const std::vector<double> &weights      = quadrature.get_weights();

    if (update_flags & update_JxW_values)
      {
        output_data.JxW_values.resize(n_q_points);

        AssertDimension(output_data.JxW_values.size(), n_q_points);
        AssertDimension(dim, spacedim);

        for (unsigned int point = 0; point < n_q_points; ++point)
          output_data.JxW_values[point] =
            weights[point] * data.contravariant[point].determinant();
      }



    // Assert(false, ExcNotImplemented());
  }

} // namespace Tet

DEAL_II_NAMESPACE_CLOSE

#endif