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

#ifndef dealii_fe_q_tet_h
#define dealii_fe_q_tet_h

#include <deal.II/base/config.h>

#include <deal.II/fe/fe_poly.h>

#include <deal.II/tet/polynomials.h>

DEAL_II_NAMESPACE_OPEN

namespace Tet
{
  template <int dim, int spacedim = dim>
  class FE_Q : public FE_Poly<Tet::ScalarPolynomial<dim>, dim>
  {
  public:
    FE_Q(const unsigned int degree);

    std::unique_ptr<FiniteElement<dim, dim>>
    clone() const override;

    std::string
    get_name() const override;

  private:
    static std::vector<unsigned int>
    get_dpo_vector(const unsigned int deg)
    {
      // AssertDimension(deg, 1);
      (void)deg;
      AssertDimension(dim, 2);

      std::vector<unsigned int> dpo(dim + 1, 0U);
      dpo[0] = 1;
      return dpo;
    }

  protected:
    virtual std::unique_ptr<
      typename FiniteElement<dim, spacedim>::InternalDataBase>
    get_data(const UpdateFlags update_flags,
             const Mapping<dim, spacedim> & /*mapping*/,
             const Quadrature<dim> &quadrature,
             dealii::internal::FEValuesImplementation::FiniteElementRelatedData<
               dim,
               spacedim> &output_data) const override
    {
      // generate a new data object and
      // initialize some fields
      std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
            data_ptr = std_cxx14::make_unique<InternalData>();
      auto &data     = dynamic_cast<InternalData &>(*data_ptr);
      // data.update_each = requires_update_flags(update_flags); // TODO
      data.update_each = update_flags;

      const unsigned int n_q_points    = quadrature.size();
      const unsigned int dofs_per_cell = dim == 2 ? 3 : 4; // TODO

      // initialize some scratch arrays. we need them for the underlying
      // polynomial to put the values and derivatives of shape functions
      // to put there, depending on what the user requested
      std::vector<double> values(update_flags & update_values ? dofs_per_cell :
                                                                0);
      std::vector<Tensor<1, dim>> grads(
        update_flags & update_gradients ? dofs_per_cell : 0);
      std::vector<Tensor<2, dim>> grad_grads(
        update_flags & update_hessians ? dofs_per_cell : 0);
      std::vector<Tensor<3, dim>> third_derivatives(
        update_flags & update_3rd_derivatives ? dofs_per_cell : 0);
      std::vector<Tensor<4, dim>>
        fourth_derivatives; // won't be needed, so leave empty

      // now also initialize fields the fields of this class's own
      // temporary storage, depending on what we need for the given
      // update flags.
      //
      // there is one exception from the rule: if we are dealing with
      // cells (i.e., if this function is not called via
      // get_(sub)face_data()), then we can already store things in the
      // final location where FEValues::reinit() later wants to see
      // things. we then don't need the intermediate space. we determine
      // whether we are on a cell by asking whether the number of
      // elements in the output array equals the number of quadrature
      // points (yes, it's a cell) or not (because in that case the
      // number of quadrature points we use here equals the number of
      // quadrature points summed over *all* faces or subfaces, whereas
      // the number of output slots equals the number of quadrature
      // points on only *one* face)
      if ((update_flags & update_values) &&
          !((output_data.shape_values.n_rows() > 0) &&
            (output_data.shape_values.n_cols() == n_q_points)))
        data.shape_values.reinit(dofs_per_cell, n_q_points);

      if (update_flags & update_gradients)
        data.shape_gradients.reinit(dofs_per_cell, n_q_points);

      if (update_flags & update_hessians)
        data.shape_hessians.reinit(dofs_per_cell, n_q_points);

      if (update_flags & update_3rd_derivatives)
        data.shape_3rd_derivatives.reinit(dofs_per_cell, n_q_points);

      // next already fill those fields of which we have information by
      // now. note that the shape gradients are only those on the unit
      // cell, and need to be transformed when visiting an actual cell
      if (update_flags & (update_values | update_gradients | update_hessians |
                          update_3rd_derivatives))
        for (unsigned int i = 0; i < n_q_points; ++i)
          {
            this->poly_space.evaluate(quadrature.point(i),
                                      values,
                                      grads,
                                      grad_grads,
                                      third_derivatives,
                                      fourth_derivatives);

            // the values of shape functions at quadrature points don't change.
            // consequently, write these values right into the output array if
            // we can, i.e., if the output array has the correct size. this is
            // the case on cells. on faces, we already precompute data on *all*
            // faces and subfaces, but we later on copy only a portion of it
            // into the output object; in that case, copy the data from all
            // faces into the scratch object
            if (update_flags & update_values)
              if (output_data.shape_values.n_rows() > 0)
                {
                  if (output_data.shape_values.n_cols() == n_q_points)
                    for (unsigned int k = 0; k < dofs_per_cell; ++k)
                      output_data.shape_values[k][i] = values[k];
                  else
                    for (unsigned int k = 0; k < dofs_per_cell; ++k)
                      data.shape_values[k][i] = values[k];
                }

            // for everything else, derivatives need to be transformed,
            // so we write them into our scratch space and only later
            // copy stuff into where FEValues wants it
            if (update_flags & update_gradients)
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  data.shape_gradients[k][i] = grads[k];
                }

            if (update_flags & update_hessians)
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                data.shape_hessians[k][i] = grad_grads[k];

            if (update_flags & update_3rd_derivatives)
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                data.shape_3rd_derivatives[k][i] = third_derivatives[k];
          }
      return data_ptr;
    }

    virtual void
    fill_fe_values(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      const CellSimilarity::Similarity                         cell_similarity,
      const Quadrature<dim> &                                  quadrature,
      const Mapping<dim, spacedim> &                           mapping,
      const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
      const dealii::internal::FEValuesImplementation::
        MappingRelatedData<dim, spacedim> &mapping_data,
      const typename FiniteElement<dim, spacedim>::InternalDataBase
        &fe_internal,
      dealii::internal::FEValuesImplementation::
        FiniteElementRelatedData<dim, spacedim> &output_data) const override
    {
      (void)cell;

      // convert data object to internal data for this class. fails with an
      // exception if that is not possible
      Assert(dynamic_cast<const InternalData *>(&fe_internal) != nullptr,
             ExcInternalError());
      const InternalData &fe_data =
        static_cast<const InternalData &>(fe_internal); // NOLINT


      const unsigned int dofs_per_cell = dim == 2 ? 3 : 4; // TODO

      // transform gradients and higher derivatives. there is nothing to do
      // for values since we already emplaced them into output_data when
      // we were in get_data()
      if (fe_data.update_each & update_gradients &&
          cell_similarity != CellSimilarity::translation)
        for (unsigned int k = 0; k < dofs_per_cell; ++k)
          mapping.transform(make_array_view(fe_data.shape_gradients, k),
                            mapping_covariant,
                            mapping_internal,
                            make_array_view(output_data.shape_gradients, k));

      if (fe_data.update_each & update_hessians &&
          cell_similarity != CellSimilarity::translation)
        {
          for (unsigned int k = 0; k < dofs_per_cell; ++k)
            mapping.transform(make_array_view(fe_data.shape_hessians, k),
                              mapping_covariant_gradient,
                              mapping_internal,
                              make_array_view(output_data.shape_hessians, k));

          for (unsigned int k = 0; k < dofs_per_cell; ++k)
            for (unsigned int i = 0; i < quadrature.size(); ++i)
              for (unsigned int j = 0; j < 2; ++j)
                output_data.shape_hessians[k][i] -=
                  mapping_data.jacobian_pushed_forward_grads[i][j] *
                  output_data.shape_gradients[k][i][j];
        }

      if (fe_data.update_each & update_3rd_derivatives &&
          cell_similarity != CellSimilarity::translation)
        {
          for (unsigned int k = 0; k < dofs_per_cell; ++k)
            mapping.transform(make_array_view(fe_data.shape_3rd_derivatives, k),
                              mapping_covariant_hessian,
                              mapping_internal,
                              make_array_view(output_data.shape_3rd_derivatives,
                                              k));

          // for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
          //  correct_third_derivatives(output_data,
          //                            mapping_data,
          //                            quadrature.size(),
          //                            k);
        }
    }


    class InternalData : public FiniteElement<dim, dim>::InternalDataBase
    {
    public:
      /**
       * Array with shape function values in quadrature points. There is one row
       * for each shape function, containing values for each quadrature point.
       *
       * In this array, we store the values of the shape function in the
       * quadrature points on the unit cell. Since these values do not change
       * under transformation to the real cell, we only need to copy them over
       * when visiting a concrete cell.
       */
      Table<2, double> shape_values;

      /**
       * Array with shape function gradients in quadrature points. There is one
       * row for each shape function, containing values for each quadrature
       * point.
       *
       * We store the gradients in the quadrature points on the unit cell. We
       * then only have to apply the transformation (which is a matrix-vector
       * multiplication) when visiting an actual cell.
       */
      Table<2, Tensor<1, dim>> shape_gradients;

      /**
       * Array with shape function hessians in quadrature points. There is one
       * row for each shape function, containing values for each quadrature
       * point.
       *
       * We store the hessians in the quadrature points on the unit cell. We
       * then only have to apply the transformation when visiting an actual
       * cell.
       */
      Table<2, Tensor<2, dim>> shape_hessians;

      /**
       * Array with shape function third derivatives in quadrature points. There
       * is one row for each shape function, containing values for each
       * quadrature point.
       *
       * We store the third derivatives in the quadrature points on the unit
       * cell. We then only have to apply the transformation when visiting an
       * actual cell.
       */
      Table<2, Tensor<3, dim>> shape_3rd_derivatives;
    };
  };

} // namespace Tet

DEAL_II_NAMESPACE_CLOSE

#endif
