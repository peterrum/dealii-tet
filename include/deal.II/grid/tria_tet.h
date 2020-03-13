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

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim = dim>
class TetTriangulation : public Triangulation<dim, spacedim>
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

  using CellStatus = typename dealii::Triangulation<dim, spacedim>::CellStatus;

  std::vector<types::global_dof_index> vertex_indices;
  std::vector<types::global_dof_index> vertex_indices_ptr;

  void
  setup()
  {
    this->levels.clear();
    this->levels.push_back(
      std::make_unique<
        internal::TriangulationImplementation::TriaLevel<dim>>());


    const unsigned int n_cell = 1;

    this->levels[0]->reserve_space(n_cell, dim, spacedim);

    // step 1: used
    this->levels[0]->cells.used.clear();
    this->levels[0]->cells.used.resize(1);

    for (unsigned int i = 0; i < n_cell; i++)
      this->levels[0]->cells.used[i] = true;

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
    this->faces =
      std::make_unique<internal::TriangulationImplementation::TriaFaces<dim>>();

    // vertices
    vertex_indices_ptr.push_back(vertex_indices.size());
    vertex_indices.push_back(0);
    vertex_indices.push_back(1);
    vertex_indices.push_back(2);

    vertex_indices_ptr.push_back(vertex_indices.size());

    this->vertices.emplace_back(Point<spacedim>(1, 0));
    this->vertices.emplace_back(Point<spacedim>(0, 1));
    this->vertices.emplace_back(Point<spacedim>(0, 0));

    // this->vertices.resize(3);
  }

  virtual void
  clear()
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  set_mesh_smoothing(
    const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
      mesh_smoothing)
  {
    Assert(false, ExcNotImplemented());
    (void)mesh_smoothing;
  }
  virtual const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing &
  get_mesh_smoothing() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  set_manifold(const types::manifold_id       number,
               const Manifold<dim, spacedim> &manifold_object)
  {
    Assert(false, ExcNotImplemented());
  }

  DEAL_II_DEPRECATED
  virtual void
  set_manifold(const types::manifold_id number)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  reset_manifold(const types::manifold_id manifold_number)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  reset_all_manifolds()
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  set_all_manifold_ids(const types::manifold_id number)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  set_all_manifold_ids_on_boundary(const types::manifold_id number)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  set_all_manifold_ids_on_boundary(const types::boundary_id b_id,
                                   const types::manifold_id number)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual const Manifold<dim, spacedim> &
  get_manifold(const types::manifold_id number) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual std::vector<types::boundary_id>
  get_boundary_ids() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual std::vector<types::manifold_id>
  get_manifold_ids() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  copy_triangulation(const Triangulation<dim, spacedim> &other_tria)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  create_triangulation(const std::vector<Point<spacedim>> &vertices,
                       const std::vector<CellData<dim>> &  cells,
                       const SubCellData &                 subcelldata)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  create_triangulation_compatibility(
    const std::vector<Point<spacedim>> &vertices,
    const std::vector<CellData<dim>> &  cells,
    const SubCellData &                 subcelldata)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  flip_all_direction_flags()
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  set_all_refine_flags()
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  refine_global(const unsigned int times = 1)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  execute_coarsening_and_refinement()
  {
    Assert(false, ExcNotImplemented());
  }

  virtual bool
  prepare_coarsening_and_refinement()
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_refine_flags(std::ostream &out) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_refine_flags(std::vector<bool> &v) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_refine_flags(std::istream &in)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_refine_flags(const std::vector<bool> &v)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_coarsen_flags(std::ostream &out) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_coarsen_flags(std::vector<bool> &v) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_coarsen_flags(std::istream &out)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_coarsen_flags(const std::vector<bool> &v)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual bool
  get_anisotropic_refinement_flag() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  clear_user_flags()
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_flags(std::ostream &out) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_flags(std::vector<bool> &v) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_flags(std::istream &in)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_flags(const std::vector<bool> &v)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  clear_user_flags_line()
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_flags_line(std::ostream &out) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_flags_line(std::vector<bool> &v) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_flags_line(std::istream &in)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_flags_line(const std::vector<bool> &v)
  {
    Assert(false, ExcNotImplemented());
  }


  virtual void
  clear_user_flags_quad()
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_flags_quad(std::ostream &out) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_flags_quad(std::vector<bool> &v) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_flags_quad(std::istream &in)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_flags_quad(const std::vector<bool> &v)
  {
    Assert(false, ExcNotImplemented());
  }


  virtual void
  clear_user_flags_hex()
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_flags_hex(std::ostream &out) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_flags_hex(std::vector<bool> &v) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_flags_hex(std::istream &in)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_flags_hex(const std::vector<bool> &v)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  clear_user_data()
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_indices(std::vector<unsigned int> &v) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_indices(const std::vector<unsigned int> &v)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_pointers(std::vector<void *> &v) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_pointers(const std::vector<void *> &v)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_indices_line(std::vector<unsigned int> &v) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_indices_line(const std::vector<unsigned int> &v)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_indices_quad(std::vector<unsigned int> &v) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_indices_quad(const std::vector<unsigned int> &v)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_indices_hex(std::vector<unsigned int> &v) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_indices_hex(const std::vector<unsigned int> &v)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_pointers_line(std::vector<void *> &v) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_pointers_line(const std::vector<void *> &v)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_pointers_quad(std::vector<void *> &v) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_pointers_quad(const std::vector<void *> &v)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  save_user_pointers_hex(std::vector<void *> &v) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  load_user_pointers_hex(const std::vector<void *> &v)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual cell_iterator
  begin(const unsigned int level = 0) const
  {
    AssertDimension(level, 0);
    return cell_iterator(
      const_cast<dealii::Triangulation<dim, spacedim> *>(
        dynamic_cast<const dealii::Triangulation<dim, spacedim> *>(this)),
      0,
      0);
  }

  virtual active_cell_iterator
  begin_active(const unsigned int level = 0) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual cell_iterator
  end() const
  {
    return cell_iterator(
      const_cast<dealii::Triangulation<dim, spacedim> *>(
        dynamic_cast<const dealii::Triangulation<dim, spacedim> *>(this)),
      -1,
      -1);
  }

  virtual cell_iterator
  end(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual active_cell_iterator
  end_active(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual cell_iterator
  last() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual active_cell_iterator
  last_active() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual IteratorRange<cell_iterator>
  cell_iterators() const
  {
    return IteratorRange<typename Triangulation<dim, spacedim>::cell_iterator>(
      begin(), end());
  }

  virtual IteratorRange<active_cell_iterator>
  active_cell_iterators() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual IteratorRange<cell_iterator>
  cell_iterators_on_level(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual IteratorRange<active_cell_iterator>
  active_cell_iterators_on_level(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual face_iterator
  begin_face() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual active_face_iterator
  begin_active_face() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual face_iterator
  end_face() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual IteratorRange<active_face_iterator>
  active_face_iterators() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual vertex_iterator
  begin_vertex() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual active_vertex_iterator
  begin_active_vertex() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual vertex_iterator
  end_vertex() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_lines() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_lines(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_active_lines() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_active_lines(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_quads() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_quads(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_active_quads() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_active_quads(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_hexs() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_hexs(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_active_hexs() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_active_hexs(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_cells() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_cells(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_active_cells() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual types::global_dof_index
  n_global_active_cells() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_active_cells(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_faces() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_active_faces() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_levels() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_global_levels() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual bool
  has_hanging_nodes() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_vertices() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual const std::vector<Point<spacedim>> &
  get_vertices() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_used_vertices() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual bool
  vertex_used(const unsigned int index) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual const std::vector<bool> &
  get_used_vertices() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  max_adjacent_cells() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual types::subdomain_id
  locally_owned_subdomain() const
  {
    return 0; // [TODO]: make parallel
  }

  virtual Triangulation<dim, spacedim> &
  get_triangulation()
  {
    Assert(false, ExcNotImplemented());
  }

  virtual const Triangulation<dim, spacedim> &
  get_triangulation() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_raw_lines() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_raw_lines(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_raw_quads() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_raw_quads(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_raw_hexs(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_raw_cells(const unsigned int level) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  n_raw_faces() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual std::size_t
  memory_consumption() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual void
  add_periodicity(
    const std::vector<GridTools::PeriodicFacePair<cell_iterator>> &)
  {
    Assert(false, ExcNotImplemented());
  }

  virtual const std::map<
    std::pair<cell_iterator, unsigned int>,
    std::pair<std::pair<cell_iterator, unsigned int>, std::bitset<3>>> &
  get_periodic_face_map() const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual unsigned int
  coarse_cell_id_to_coarse_cell_index(
    const types::coarse_cell_id coarse_cell_id) const
  {
    Assert(false, ExcNotImplemented());
  }

  virtual types::coarse_cell_id
  coarse_cell_index_to_coarse_cell_id(
    const unsigned int coarse_cell_index) const
  {
    return coarse_cell_index;
  }
};

DEAL_II_NAMESPACE_CLOSE

#endif
