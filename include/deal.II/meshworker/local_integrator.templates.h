// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2017 by the deal.II authors
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
#ifndef dealii_local_integrator_templates_h
#define dealii_local_integrator_templates_h

#include <deal.II/lac/block_indices.h>

#include <deal.II/meshworker/local_integrator.h>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  template <int dim, int spacedim, typename number, typename CellIterator, typename FaceIterator>
  LocalIntegrator<dim, spacedim, number, CellIterator, FaceIterator>::LocalIntegrator()
    : use_cell(true)
    , use_boundary(true)
    , use_face(true)
  {}


  template <int dim, int spacedim, typename number, typename CellIterator, typename FaceIterator>
  LocalIntegrator<dim, spacedim, number, CellIterator, FaceIterator>::LocalIntegrator(bool c,
                                                          bool b,
                                                          bool f)
    : use_cell(c)
    , use_boundary(b)
    , use_face(f)
  {}



  template <int dim, int spacedim, typename number, typename CellIterator, typename FaceIterator>
  void
  LocalIntegrator<dim, spacedim, number, CellIterator, FaceIterator>::cell(
    DoFInfo<dim, spacedim, number, CellIterator, FaceIterator> &,
    IntegrationInfo<dim, spacedim> &) const
  {
    Assert(false, ExcPureFunction());
  }


  template <int dim, int spacedim, typename number, typename CellIterator, typename FaceIterator>
  void
  LocalIntegrator<dim, spacedim, number, CellIterator, FaceIterator>::boundary(
    DoFInfo<dim, spacedim, number, CellIterator, FaceIterator> &,
    IntegrationInfo<dim, spacedim> &) const
  {
    Assert(false, ExcPureFunction());
  }


  template <int dim, int spacedim, typename number, typename CellIterator, typename FaceIterator>
  void
  LocalIntegrator<dim, spacedim, number, CellIterator, FaceIterator>::face(
    DoFInfo<dim, spacedim, number, CellIterator, FaceIterator> &,
    DoFInfo<dim, spacedim, number, CellIterator, FaceIterator> &,
    IntegrationInfo<dim, spacedim> &,
    IntegrationInfo<dim, spacedim> &) const
  {
    Assert(false, ExcPureFunction());
  }
} // namespace MeshWorker


DEAL_II_NAMESPACE_CLOSE

#endif
