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


//#include <deal.II/base/config.h>

//#include <deal.II/grid/tria.h>

#include <deal.II/grid/tria_tet_connectivity.h>

DEAL_II_NAMESPACE_OPEN

namespace Tet
{
  template <int dim>
  struct CellData
  {
    CellTypeEnum type;

    std::vector<unsigned int> vertices;
  };

  template <int dim, int spacedim = dim>
  class Triangulation : public dealii::Triangulation<dim, spacedim>
  {
  public:
    using cell_iterator =
      typename dealii::Triangulation<dim, spacedim>::cell_iterator;

    using active_cell_iterator =
      typename dealii::Triangulation<dim, spacedim>::active_cell_iterator;

    using face_iterator =
      typename dealii::Triangulation<dim, spacedim>::face_iterator;

    using active_face_iterator =
      typename dealii::Triangulation<dim, spacedim>::active_face_iterator;

    using vertex_iterator =
      typename dealii::Triangulation<dim, spacedim>::vertex_iterator;

    using active_vertex_iterator =
      typename dealii::Triangulation<dim, spacedim>::active_vertex_iterator;

    using line_iterator =
      typename dealii::Triangulation<dim, spacedim>::line_iterator;

    using active_line_iterator =
      typename dealii::Triangulation<dim, spacedim>::active_line_iterator;

    using quad_iterator =
      typename dealii::Triangulation<dim, spacedim>::quad_iterator;

    using active_quad_iterator =
      typename dealii::Triangulation<dim, spacedim>::active_quad_iterator;

    using hex_iterator =
      typename dealii::Triangulation<dim, spacedim>::hex_iterator;

    using active_hex_iterator =
      typename dealii::Triangulation<dim, spacedim>::active_hex_iterator;

    using CellStatus =
      typename dealii::Triangulation<dim, spacedim>::CellStatus;

    std::vector<types::global_dof_index> vertex_indices;
    std::vector<types::global_dof_index> vertex_indices_ptr;

    Connectivity<dim> connectivity;

    void
    create_triangulation_tet(const std::vector<Point<spacedim>> &vertices,
                             const std::vector<CellData<dim>> &  cells)
    {
      {
        std::vector<CellTypeEnum> cell_types;
        std::vector<unsigned int> cell_vertices;

        for (const auto &cell : cells)
          {
            cell_types.push_back(cell.type);
            for (const auto &vertex : cell.vertices)
              cell_vertices.push_back(vertex);
          }

        connectivity.build(cell_types, cell_vertices);
      }

      this->levels.clear();
      this->levels.push_back(
        std::make_unique<
          internal::TriangulationImplementation::TriaLevel<dim>>());


      const unsigned int n_cell = cells.size();

      this->levels[0]->reserve_space(n_cell, dim, spacedim);

      // step 1: used
      this->levels[0]->cells.used.clear();
      this->levels[0]->cells.used.resize(n_cell);

      for (unsigned int i = 0; i < n_cell; i++)
        this->levels[0]->cells.used[i] = true;

      // if removed obly one cell is generated. 2nd argument n_cell?
      this->levels[0]->cells.reserve_space(0, 0); // TODO: what is happening here?

      // step 2: material id
      this->levels[0]->cells.boundary_or_material_id.clear();
      this->levels[0]->cells.boundary_or_material_id.resize(n_cell);
      for (unsigned int i = 0; i < n_cell; i++)
        this->levels[0]->cells.boundary_or_material_id[i].material_id = 0;

      // step 3: manifold id
      this->levels[0]->cells.manifold_id.clear();
      this->levels[0]->cells.manifold_id.resize(n_cell);
      for (unsigned int i = 0; i < n_cell; i++)
        this->levels[0]->cells.manifold_id[i] = 0;

      // step 4: subdomain id
      this->levels[0]->subdomain_ids.clear();
      this->levels[0]->subdomain_ids.resize(n_cell);
      for (unsigned int i = 0; i < n_cell; i++)
        this->levels[0]->subdomain_ids[i] = 0;

      // step 5: level_subdomain id
      this->levels[0]->level_subdomain_ids.clear();
      this->levels[0]->level_subdomain_ids.resize(n_cell);
      for (unsigned int i = 0; i < n_cell; i++)
        this->levels[0]->level_subdomain_ids[i] = 0;

      // step 6: no children...
      this->levels[0]->cells.children.clear();
      this->levels[0]->cells.children.resize(
        GeometryInfo<dim>::max_children_per_cell / 2 * n_cell);
      for (unsigned int i = 0;
           i < GeometryInfo<dim>::max_children_per_cell / 2 * n_cell;
           i++)
        this->levels[0]->cells.children[i] = -1;

      // faces [TODO]
      this->faces = std::make_unique<
        internal::TriangulationImplementation::TriaFaces<dim>>();

      // vertices
      vertex_indices_ptr.push_back(vertex_indices.size());

      for (const auto &cell : cells)
        {
          for (const auto vertex : cell.vertices)
            vertex_indices.push_back(vertex);

          vertex_indices_ptr.push_back(vertex_indices.size());
        }

      this->vertices = vertices;
    }

    virtual void
    clear() override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    set_mesh_smoothing(
      const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
        mesh_smoothing) override
    {
      Assert(false, ExcNotImplemented());
      (void)mesh_smoothing;
    }
    virtual const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing &
    get_mesh_smoothing() const override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    set_manifold(const types::manifold_id       number,
                 const Manifold<dim, spacedim> &manifold_object) override
    {
      Assert(false, ExcNotImplemented());
      (void)number;
      (void)manifold_object;
    }

    DEAL_II_DEPRECATED
    virtual void
    set_manifold(const types::manifold_id number) override
    {
      Assert(false, ExcNotImplemented());
      (void)number;
    }

    virtual void
    reset_manifold(const types::manifold_id manifold_number) override
    {
      Assert(false, ExcNotImplemented());
      (void)manifold_number;
    }

    virtual void
    reset_all_manifolds() override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    set_all_manifold_ids(const types::manifold_id number) override
    {
      Assert(false, ExcNotImplemented());
      (void)number;
    }

    virtual void
    set_all_manifold_ids_on_boundary(const types::manifold_id number) override
    {
      Assert(false, ExcNotImplemented());
      (void)number;
    }

    virtual void
    set_all_manifold_ids_on_boundary(const types::boundary_id b_id,
                                     const types::manifold_id number) override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual const Manifold<dim, spacedim> &
    get_manifold(const types::manifold_id number) const override
    {
      Assert(false, ExcNotImplemented());
      (void)number;
    }

    virtual std::vector<types::boundary_id>
    get_boundary_ids() const override
    {
      Assert(false, ExcNotImplemented());
      return std::vector<types::boundary_id>();
    }

    virtual std::vector<types::manifold_id>
    get_manifold_ids() const override
    {
      Assert(false, ExcNotImplemented());
      return std::vector<types::manifold_id>();
    }

    virtual void
    copy_triangulation(
      const dealii::Triangulation<dim, spacedim> &other_tria) override
    {
      Assert(false, ExcNotImplemented());
      (void)other_tria;
    }

    virtual void
    create_triangulation(const std::vector<Point<spacedim>> &      vertices,
                         const std::vector<dealii::CellData<dim>> &cells,
                         const SubCellData &subcelldata) override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    create_triangulation_compatibility(
      const std::vector<Point<spacedim>> &      vertices,
      const std::vector<dealii::CellData<dim>> &cells,
      const SubCellData &                       subcelldata) override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    flip_all_direction_flags() override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    set_all_refine_flags() override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    refine_global(const unsigned int times = 1) override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    execute_coarsening_and_refinement() override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual bool
    prepare_coarsening_and_refinement() override
    {
      Assert(false, ExcNotImplemented());

      return false;
    }

    virtual void
    save_refine_flags(std::ostream &out) const override
    {
      Assert(false, ExcNotImplemented());
      (void)out;
    }

    virtual void
    save_refine_flags(std::vector<bool> &v) const override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    load_refine_flags(std::istream &in) override
    {
      Assert(false, ExcNotImplemented());
      (void)in;
    }

    virtual void
    load_refine_flags(const std::vector<bool> &v) override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    save_coarsen_flags(std::ostream &out) const override
    {
      Assert(false, ExcNotImplemented());
      (void)out;
    }

    virtual void
    save_coarsen_flags(std::vector<bool> &v) const override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    load_coarsen_flags(std::istream &out) override
    {
      Assert(false, ExcNotImplemented());
      (void)out;
    }

    virtual void
    load_coarsen_flags(const std::vector<bool> &v) override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual bool
    get_anisotropic_refinement_flag() const override
    {
      Assert(false, ExcNotImplemented());

      return false;
    }

    virtual void
    clear_user_flags() override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    save_user_flags(std::ostream &out) const override
    {
      Assert(false, ExcNotImplemented());
      (void)out;
    }

    virtual void
    save_user_flags(std::vector<bool> &v) const override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    load_user_flags(std::istream &in) override
    {
      Assert(false, ExcNotImplemented());
      (void)in;
    }

    virtual void
    load_user_flags(const std::vector<bool> &v) override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    clear_user_flags_line() override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    save_user_flags_line(std::ostream &out) const override
    {
      Assert(false, ExcNotImplemented());
      (void)out;
    }

    virtual void
    save_user_flags_line(std::vector<bool> &v) const override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    load_user_flags_line(std::istream &in) override
    {
      Assert(false, ExcNotImplemented());
      (void)in;
    }

    virtual void
    load_user_flags_line(const std::vector<bool> &v) override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }


    virtual void
    clear_user_flags_quad() override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    save_user_flags_quad(std::ostream &out) const override
    {
      Assert(false, ExcNotImplemented());
      (void)out;
    }

    virtual void
    save_user_flags_quad(std::vector<bool> &v) const override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    load_user_flags_quad(std::istream &in) override
    {
      Assert(false, ExcNotImplemented());
      (void)in;
    }

    virtual void
    load_user_flags_quad(const std::vector<bool> &v) override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }


    virtual void
    clear_user_flags_hex() override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    save_user_flags_hex(std::ostream &out) const override
    {
      Assert(false, ExcNotImplemented());
      (void)out;
    }

    virtual void
    save_user_flags_hex(std::vector<bool> &v) const override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    load_user_flags_hex(std::istream &in) override
    {
      Assert(false, ExcNotImplemented());
      (void)in;
    }

    virtual void
    load_user_flags_hex(const std::vector<bool> &v) override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    clear_user_data() override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual void
    save_user_indices(std::vector<unsigned int> &v) const override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    load_user_indices(const std::vector<unsigned int> &v) override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    save_user_pointers(std::vector<void *> &v) const override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    load_user_pointers(const std::vector<void *> &v) override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    save_user_indices_line(std::vector<unsigned int> &v) const override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    load_user_indices_line(const std::vector<unsigned int> &v) override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    save_user_indices_quad(std::vector<unsigned int> &v) const override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    load_user_indices_quad(const std::vector<unsigned int> &v) override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    save_user_indices_hex(std::vector<unsigned int> &v) const override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    load_user_indices_hex(const std::vector<unsigned int> &v) override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    save_user_pointers_line(std::vector<void *> &v) const override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    load_user_pointers_line(const std::vector<void *> &v) override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    save_user_pointers_quad(std::vector<void *> &v) const override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    load_user_pointers_quad(const std::vector<void *> &v) override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    save_user_pointers_hex(std::vector<void *> &v) const override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual void
    load_user_pointers_hex(const std::vector<void *> &v) override
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }

    virtual cell_iterator
    begin(const unsigned int level = 0) const override
    {
      AssertDimension(level, 0);
      return cell_iterator(
        const_cast<dealii::Triangulation<dim, spacedim> *>(
          dynamic_cast<const dealii::Triangulation<dim, spacedim> *>(this)),
        0,
        0);
    }

    virtual active_cell_iterator
    begin_active(const unsigned int level = 0) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
    }

    virtual cell_iterator
    end() const override
    {
      return cell_iterator(
        const_cast<dealii::Triangulation<dim, spacedim> *>(
          dynamic_cast<const dealii::Triangulation<dim, spacedim> *>(this)),
        -1,
        -1);
    }

    virtual cell_iterator
    end(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
    }

    virtual active_cell_iterator
    end_active(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
    }

    virtual cell_iterator
    last() const override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual active_cell_iterator
    last_active() const override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual IteratorRange<cell_iterator>
    cell_iterators() const override
    {
      return IteratorRange<
        typename Triangulation<dim, spacedim>::cell_iterator>(begin(), end());
    }

    virtual IteratorRange<active_cell_iterator>
    active_cell_iterators() const override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual IteratorRange<cell_iterator>
    cell_iterators_on_level(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
    }

    virtual IteratorRange<active_cell_iterator>
    active_cell_iterators_on_level(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
    }

    virtual face_iterator
    begin_face() const override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual active_face_iterator
    begin_active_face() const override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual face_iterator
    end_face() const override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual IteratorRange<active_face_iterator>
    active_face_iterators() const override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual vertex_iterator
    begin_vertex() const override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual active_vertex_iterator
    begin_active_vertex() const override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual vertex_iterator
    end_vertex() const override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual unsigned int
    n_lines() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual unsigned int
    n_lines(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
      return 0;
    }

    virtual unsigned int
    n_active_lines() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual unsigned int
    n_active_lines(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
      return 0;
    }

    virtual unsigned int
    n_quads() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual unsigned int
    n_quads(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
      return 0;
    }

    virtual unsigned int
    n_active_quads() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual unsigned int
    n_active_quads(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
      return 0;
    }

    virtual unsigned int
    n_hexs() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual unsigned int
    n_hexs(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
      return 0;
    }

    virtual unsigned int
    n_active_hexs() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual unsigned int
    n_active_hexs(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
      return 0;
    }

    virtual unsigned int
    n_cells() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual unsigned int
    n_cells(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
      return 0;
    }

    virtual unsigned int
    n_active_cells() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual types::global_dof_index
    n_global_active_cells() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual unsigned int
    n_active_cells(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
      return 0;
    }

    virtual unsigned int
    n_faces() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual unsigned int
    n_active_faces() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual unsigned int
    n_levels() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual unsigned int
    n_global_levels() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual bool
    has_hanging_nodes() const override
    {
      Assert(false, ExcNotImplemented());
      return false;
    }

    virtual unsigned int
    n_vertices() const override
    {
      return this->vertices.size();
    }

    virtual const std::vector<Point<spacedim>> &
    get_vertices() const override
    {
      return this->vertices;
    }

    virtual unsigned int
    n_used_vertices() const override
    {
      Assert(false, ExcNotImplemented());
      return this->vertices.size();
    }

    virtual bool
    vertex_used(const unsigned int index) const override
    {
      Assert(false, ExcNotImplemented());
      (void)index;
      return true;
    }

    virtual const std::vector<bool> &
    get_used_vertices() const override
    {
      Assert(false, ExcNotImplemented());
      return this->vertices_used;
    }

    virtual unsigned int
    max_adjacent_cells() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual types::subdomain_id
    locally_owned_subdomain() const override
    {
      return 0; // [TODO]: make parallel
    }

    virtual Triangulation<dim, spacedim> &
    get_triangulation() override
    {
      Assert(false, ExcNotImplemented());
      return *this;
    }

    virtual const Triangulation<dim, spacedim> &
    get_triangulation() const override
    {
      Assert(false, ExcNotImplemented());
      return *this;
    }

    virtual unsigned int
    n_raw_lines() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual unsigned int
    n_raw_lines(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
      return 0;
    }

    virtual unsigned int
    n_raw_quads() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual unsigned int
    n_raw_quads(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
      return 0;
    }

    virtual unsigned int
    n_raw_hexs(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
      return 0;
    }

    virtual unsigned int
    n_raw_cells(const unsigned int level) const override
    {
      Assert(false, ExcNotImplemented());
      (void)level;
      return 0;
    }

    virtual unsigned int
    n_raw_faces() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual std::size_t
    memory_consumption() const override
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    virtual void
    add_periodicity(
      const std::vector<GridTools::PeriodicFacePair<cell_iterator>> &) override
    {
      Assert(false, ExcNotImplemented());
    }

    virtual const std::map<
      std::pair<cell_iterator, unsigned int>,
      std::pair<std::pair<cell_iterator, unsigned int>, std::bitset<3>>> &
    get_periodic_face_map() const override
    {
      Assert(false, ExcNotImplemented());
      return this->periodic_face_map;
    }

    virtual unsigned int
    coarse_cell_id_to_coarse_cell_index(
      const types::coarse_cell_id coarse_cell_id) const override
    {
      Assert(false, ExcNotImplemented());
      return coarse_cell_id;
    }

    virtual types::coarse_cell_id
    coarse_cell_index_to_coarse_cell_id(
      const unsigned int coarse_cell_index) const override
    {
      return coarse_cell_index;
    }
  };
} // namespace Tet

DEAL_II_NAMESPACE_CLOSE

#endif
