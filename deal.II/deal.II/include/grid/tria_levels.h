//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__tria_levels_h
#define __deal2__tria_levels_h


#include <base/config.h>
#include <vector>
#include <grid/tria_object.h>
#include <base/point.h>
#include <grid/tria_objects.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace Triangulation
  {

/**
 * General template for information belonging to one level of a multilevel
 * hierarchy of a triangulation. However, only the cells are really
 * level-based data if anisotropic rfeinement is allowed. Therefore, only the
 * information concerning cells is considered here, the information concerning
 * cell-faces can be found in the <tt>TriaFaces</tt> classes.
 *
 * This template is only declared to allow
 * specializations for different dimensions.
 *
 * @ingroup grid
 */
    template <int dim>
    class TriaLevel
    {
    };


/**
 * Store all information which belongs to one level of the multilevel
 * hierarchy.
 *
 * In TriaLevel<0> all data is stored which is not dependent on the
 * dimension, e.g. a field to store the refinement flag for the cells
 * (what a cell actually is is declared elsewhere), etc. Actually, it
 * is only cell-based data, like neighborship info or refinement
 * flags. There is another field, which may fit in here, namely the
 * material data (for cells) or the boundary indicators (for faces),
 * but since we need for a line or quad either boundary information or
 * material data, we store them with the lines and quads rather than
 * with the common data. Likewise, in 3d, we need boundary indicators
 * for lines and quads (we need to know how to refine a line if the
 * two adjacent faces have different boundary indicators), and
 * material data for cells.
 *
 * @ingroup grid
 * @author Wolfgang Bangerth, 1998
 */
    template <>
    class TriaLevel<0>
    {
      public:
                                         /**
                                          *  Flags for the cells whether they are
                                          *  to be refined or not. The meaning
                                          *  what a cell is, is dimension specific,
                                          *  therefore also the length of this
                                          *  vector depends on the dimension: in
                                          *  one dimension, the length of this
                                          *  vector equals the length of the
                                          *  @p lines vector, in two dimensions
                                          *  that of the @p quads vector, etc.
                                          */
        std::vector<bool> refine_flags;

                                         /**
                                          * Same meaning as the one above, but
                                          * specifies whether a cell must be
                                          * coarsened.
                                          */
        std::vector<bool> coarsen_flags;
    
                                         /**
                                          *  Levels and indices of the neighbors
                                          *  of the cells. Convention is, that the
                                          *  neighbors of the cell with index @p i
                                          *  are stored in the fields following
                                          *  $i*(2*real_space_dimension)$, e.g. in
                                          *  one spatial dimension, the neighbors
                                          *  of cell 0 are stored in <tt>neighbors[0]</tt>
                                          *  and <tt>neighbors[1]</tt>, the neighbors of
                                          *  cell 1 are stored in <tt>neighbors[2]</tt>
                                          *  and <tt>neighbors[3]</tt>, and so on.
                                          *
                                          *  In neighbors, <tt>neighbors[i].first</tt> is
                                          *  the level, while <tt>neighbors[i].first</tt>
                                          *  is the index of the neighbor.
                                          *
                                          *  If a neighbor does not exist (cell is
                                          *  at the boundary), <tt>level=index=-1</tt>
                                          *  is set.
                                          *
                                          *  <em>Conventions:</em> The
                                          *  @p ith neighbor of a cell is
                                          *  the one which shares the
                                          *  @p ith face (@p Line in 2D,
                                          *  @p Quad in 3D) of this cell.
                                          *
                                          *  The neighbor of a cell has at most the
                                          *  same level as this cell, i.e. it may
                                          *  or may not be refined.
                                          *
                                          *  In one dimension, a neighbor may have
                                          *  any level less or equal the level of
                                          *  this cell. If it has the same level,
                                          *  it may be refined an arbitrary number
                                          *  of times, but the neighbor pointer
                                          *  still points to the cell on the same
                                          *  level, while the neighbors of the
                                          *  childs of the neighbor may point to
                                          *  this cell or its children.
                                          *
                                          *  In two and more dimensions, the
                                          *  neighbor is either on the same level
                                          *  and refined (in which case its children
                                          *  have neighbor pointers to this cell or
                                          *  its direct children), unrefined on
                                          *  the same level or one level down (in
                                          *  which case its neighbor pointer points
                                          *  to the mother cell of this cell).
                                          */
        std::vector<std::pair<int,int> > neighbors;

                                         /**
                                          * One integer per cell to store
                                          * which subdomain it belongs
                                          * to. This field is most often
                                          * used in parallel computations,
                                          * where it denotes which
                                          * processor shall work on the
                                          * cells with a given subdomain
                                          * number.
                                          */
        std::vector<unsigned int> subdomain_ids;
    
                                         /**
                                          *  Reserve enough space to accomodate
                                          *  @p total_cells cells on this level.
                                          *  Since there are no @p used flags on this
                                          *  level, you have to give the total number
                                          *  of cells, not only the number of newly
                                          *  to accomodate ones, like in the
                                          *  <tt>TriaLevel<N>::reserve_space</tt>
                                          *  functions, with <tt>N>0</tt>.
                                          *
                                          *  Since the
                                          *  number of neighbors per cell depends
                                          *  on the dimensions, you have to pass
                                          *  that additionally.
                                          */
        void reserve_space (const unsigned int total_cells,
                            const unsigned int dimension);

                                         /**
                                          *  Check the memory consistency of the
                                          *  different containers. Should only be
                                          *  called with the prepro flag @p DEBUG
                                          *  set. The function should be called from
                                          *  the functions of the higher
                                          *  TriaLevel classes.
                                          */
        void monitor_memory (const unsigned int true_dimension) const;

                                         /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        unsigned int memory_consumption () const;

                                         /**
                                          *  Exception
                                          */
        DeclException3 (ExcMemoryWasted,
                        char*, int, int,
                        << "The container " << arg1 << " contains "
                        << arg2 << " elements, but it`s capacity is "
                        << arg3 << ".");
                                         /**
                                          *  Exception
                                          */
        DeclException2 (ExcMemoryInexact,
                        int, int,
                        << "The containers have sizes " << arg1 << " and "
                        << arg2 << ", which is not as expected.");
    };



/**
 *  Store all information which belongs to one level of the multilevel hierarchy.
 *  
 *  In one dimension, this contains the lines on this level.
 *
 * @ingroup grid
 *  @author Wolfgang Bangerth, 1998
 */
    template <>
    class TriaLevel<1> : public TriaLevel<0>
			 
    {
      public:
					 /**
					  * The object containing the data on lines and
					  * related functions
					  */
	TriaObjects<Line> lines;

      public:
					 /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        unsigned int memory_consumption () const;
    };



/**
 *  Store all information which belongs to one level of the multilevel hierarchy.
 *
 *  In 2D this is a vector ocontains the quads on this level.
 *
 * @ingroup grid
 *  @author Wolfgang Bangerth, 1998
 */
    template <>
    class TriaLevel<2> : public TriaLevel<0>
    {
      public:
					 /**
					  * The object containing the data on quads and
					  * related functions
					  */
	TriaObjects<Quad> quads;

      public:
					 /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        unsigned int memory_consumption () const;
    };



/**
 *  Store all information which belongs to one level of the multilevel hierarchy.
 *
 *  In 3D this contains the hexahedra on this levels, 
 *
  * @ingroup grid
 *  @author Wolfgang Bangerth, 1998, 2003
 */
    template <>
    class TriaLevel<3> :  public TriaLevel<0>
    {
      public:
					 /**
					  * The object containing the data on hexes and
					  * related functions
					  */
	TriaObjectsHex hexes;

      public:
					 /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        unsigned int memory_consumption () const;
    };
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
