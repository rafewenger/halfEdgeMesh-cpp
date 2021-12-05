/// \file half_edge_mesh_DCMT.hpp
/// template classes for half edge mesh supporting decimation (DCMT) operations.

/*
  Copyright (C) 2021 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include<sstream>

#include "half_edge_mesh.hpp"
#include "half_edge_mesh_coord.hpp"

// *** DEBUG ***
// Uncomment if you want to include print statements 
//   for debugging in this file.
#include<iostream>


#ifndef _HALF_EDGE_MESH_DCMT_
#define _HALF_EDGE_MESH_DCMT_

namespace HMESH {

  /// Vertex supporting decimation.
  template <const int _DIM, typename HALF_EDGE_PTR, typename CTYPE>
  class VERTEX_DCMT_BASE: public VERTEX_BASE<_DIM,HALF_EDGE_PTR,CTYPE> 
  {
  protected:

    /// Internal flag used for detecting vertex adjacencies.
    mutable bool visited_flag;

    /// Return visited_flag;
    bool _IsVisited() const { return visited_flag; }

    /// Set visited_flag to flag in all neighbors of this.
    void _SetVisitedFlagsInAdjacentVertices(const bool flag) const;

    /// Set visited_flag to false in all neighbors of this.
    void _ClearVisitedFlagsInAdjacentVertices() const
    { _SetVisitedFlagsInAdjacentVertices(false); }

    /// Process first and last edges in half_edges_from[].
    /// - Swap if first edge is not boundary and last edge is boundary.
    /// - Swap if first edge is not boundary and 
    ///   first edge->PrevHalfEdgeInCell() is not boundary,
    ///   but last edge->PrevHalfEdgeInCell() is boundary.
    void _ProcessFirstLastHalfEdgesFrom();

  public:
    
    /// Return true if vertex is incident on more than two edges.
    bool IsIncidentOnMoreThanTwoEdges() const;

    /// Allow HALF_EDGE_MESH_DCMT_BASE to access all vertex fields.
    template <const int D, typename VTYPE, typename HTYPE, 
              typename CELL_TYPE>
    friend class HALF_EDGE_MESH_DCMT_BASE;

    /// Allow CELL_DCMT_BASE to access all vertex fields.
    template <const int D, typename HTYPE>
    friend class CELL_DCMT_BASE;
  };

  /// 3D version of VERTEX_DCMT_BASE
  template <typename HALF_EDGE_PTR, typename CTYPE>
  class VERTEX3D_DCMT_BASE:public VERTEX_DCMT_BASE<3,HALF_EDGE_PTR,CTYPE>
  {};


  /// Half edge supporting decimation.
  /// @tparam DIM Dimension of coordinates. 
  ///   Must equal dimension of *VERTEX_PTR.
  template <const int DIM, typename VERTEX_PTR, typename HALF_EDGE_PTR, 
            typename CELL_PTR>
  class HALF_EDGE_DCMT_BASE: 
    public HALF_EDGE_BASE<VERTEX_PTR, HALF_EDGE_PTR, CELL_PTR>
  {
  public:

    /// Compute square of edge length.
    double ComputeLengthSquared() const
    { return(compute_squared_distance<DIM>
             (this->FromVertex()->coord, this->ToVertex()->coord)); }

    /// Compute cosine of angle at this->FromVertex() in triangle
    ///   formed by PreviousHalfEdge()->FromVertex(),
    ///   this->FromVertex() and this->ToVertex().
    /// - Sets flag_zero to true if middle vertex has same coordinates 
    ///   as one of the other two (or middle vertex is very, very
    ///   close to one of the other two.)
    double ComputeCosAngleAtFromVertex(bool & flag_zero) const;

    /// Allow HALF_EDGE_MESH_DCMT_BASE to access all half edge fields.
    template <const int D, typename VTYPE, typename HTYPE, 
              typename CELL_TYPE>
    friend class HALF_EDGE_MESH_DCMT_BASE;
  };

  /// Cell supporting decimation.
  /// @tparam DIM Dimension of coordinates. 
  ///   Must equal dimension of *VERTEX_PTR.
  template <const int DIM, typename HALF_EDGE_PTR>
  class CELL_DCMT_BASE:public CELL_BASE<HALF_EDGE_PTR>
  {
  protected:
    /// Set visited_flag to flag in all cell vertices.
    void _SetVisitedFlagsInAllVertices(const bool flag) const;

    /// Set visited_flag to false in all cell vertices.
    void _ClearVisitedFlagsInAllVertices() const
    { _SetVisitedFlagsInAllVertices(false); }


  public:

    /// Compute min and max squared edge lengths in a cell.
    /// @param[out] min_edge_length_squared Square of minimum edge length.
    /// @param[out] max_edge_length_squared Square of maximum edge length.
    /// @param[out] ihalf_edge_min Index of half edge with squared length
    ///    equal to min_edge_length_squared.
    /// @param[out] ihalf_edge_max Index of half edge with squared length
    ///    equal to max_edge_length_squared.
    template <typename CTYPE>
    void ComputeMinMaxEdgeLengthSquared
    (CTYPE & min_edge_length_squared, CTYPE & max_edge_length_squared,
     int & ihalf_edge_min, int & ihalf_edge_max) const;

    /// Compute min and max squared edge lengths in a cell.
    /// - Version that only returns min and max edge length squared.
    template <typename CTYPE>
    void ComputeMinMaxEdgeLengthSquared
    (CTYPE & min_edge_length_squared, 
     CTYPE & max_edge_length_squared) const
    { 
      int ihalf_edge_min, ihalf_edge_max;
      ComputeMinMaxEdgeLengthSquared
        (min_edge_length_squared, max_edge_length_squared,
         ihalf_edge_min, ihalf_edge_max);
    }

    /// Compute cosine of the min or max angle between consecutive
    ///   cell edges.
    /// - Note: cos_min_angle >= cos_max_angle.
    /// @param[out] cos_min_angle Cosine of the min angle.
    /// - Note: cosine of the min angle =  max cosine A, where A is any
    ///     angle between consecutive cell edges.
    /// - The smallest angle is 0 and cos(0) = 1.
    /// @param[out] cos_max_angle Cosine of the max angle.
    /// - cosine of the max angle =  min cosine A, where A is any
    ///     angle between consecutive cell edges.
    /// - The largest angle is pi and cos(pi) = -1.
    /// @param[out] ihalf_edge_min Half edge whose from vertex
    ///   has the min angle.
    /// @param[out] ihalf_edge_max Half edge whose from vertex
    ///   has the max angle.
    template <typename CTYPE>
    void ComputeCosMinMaxAngle
    (CTYPE & cos_min_angle, CTYPE & cos_max_angle,
     int & ihalf_edge_min, int & ihalf_edge_max) const;


    /// Allow HALF_EDGE_MESH_DCMT_BASE to access all cell fields.
    template <const int D, typename VTYPE, typename HTYPE, 
              typename CELL_TYPE>
    friend class HALF_EDGE_MESH_DCMT_BASE;
  };

  // Forward declaration of simple vertex, half_edge and cell types.
  class VERTEX_DCMT_A;
  class HALF_EDGE_DCMT_A;
  class CELL_DCMT_A;

  typedef VERTEX_DCMT_A * VERTEX_DCMT_A_PTR;
  typedef HALF_EDGE_DCMT_A * HALF_EDGE_DCMT_A_PTR;
  typedef CELL_DCMT_A * CELL_DCMT_A_PTR;

  // Declaration of VERTEX_DCMT_A, HALF_EDGE_DCMT_A, CELL_DCMT_A.
  class VERTEX_DCMT_A:public VERTEX3D_DCMT_BASE<HALF_EDGE_DCMT_A_PTR,float>
  {};

  class HALF_EDGE_DCMT_A:public HALF_EDGE_DCMT_BASE
  <3, VERTEX_DCMT_A_PTR,HALF_EDGE_DCMT_A_PTR,CELL_DCMT_A_PTR>
  {};

  class CELL_DCMT_A:public CELL_DCMT_BASE<3, HALF_EDGE_DCMT_A_PTR>
  {};


  /// Half edge mesh with join, split, edge collapse, etc operations.
  /// @tparam DIM Dimension of coordinates. 
  ///   Must equal dimension of _VERTEX_TYPE.
  /// - Note: _VERTEX_TYPE, _HALF_EDGE_TYPE and _CELL_TYPE must
  ///   be derived from VERTEX_DCMT_BASE, HALF_EDGE_DCMT_BASE and
  ///   CELL_DCMT_BASE.
  template < const int DIM, typename _VERTEX_TYPE, 
             typename _HALF_EDGE_TYPE, typename _CELL_TYPE >
  class HALF_EDGE_MESH_DCMT_BASE:
    public HALF_EDGE_MESH_BASE<_VERTEX_TYPE,_HALF_EDGE_TYPE,_CELL_TYPE>
  {

  public:
    typedef _VERTEX_TYPE VERTEX_TYPE;
    typedef _HALF_EDGE_TYPE HALF_EDGE_TYPE;
    typedef _CELL_TYPE CELL_TYPE;
    typedef typename VERTEX_TYPE::COORD_TYPE COORD_TYPE;

  protected:

    // *** Internal split routines ***

    /// Split an internal edge.
    /// - Returns new vertex.
    /// @pre half_edgeA is an internal edge.
    const VERTEX_TYPE * _SplitInternalEdge(HALF_EDGE_TYPE * half_edgeA);

    /// Split a boundary edge.
    /// - Returns new vertex.
    /// @pre half_edgeA is a boundary edge.
    const VERTEX_TYPE * _SplitBoundaryEdge(HALF_EDGE_TYPE * half_edgeA);


    // *** Internal routines ***

    /// Remove half edge from the half_edge_from list of its from_vertex.
    /// - Does not ensure that first element is a boundary half edge.
    /// - Call _MoveBoundaryHalfEdgeToHalfEdgeFrom0() to ensure
    ///   that first half edge in vertex list is a boundary edge.
    void _RemoveHalfEdgeFromVertexList
    (const HALF_EDGE_TYPE * half_edge0);

    /// Move all half edges from vA to be from vB.
    /// - Moves vA->half_edge_list to vB->half_edge_list.
    /// - Clears vA->half_edge_lis[].
    void _MoveVertexHalfEdgeFromList
    (VERTEX_TYPE * vA, VERTEX_TYPE *vB);

    /// Swap next_half_edge_around_edge.
    void _SwapNextHalfEdgeAroundEdge
    (HALF_EDGE_TYPE * half_edgeA, HALF_EDGE_TYPE * half_edgeB);

    /// Return non constant pointer to previous half edge around edge.
    HALF_EDGE_TYPE * 
    _FindPrevHalfEdgeAroundEdgeNC(HALF_EDGE_TYPE * half_edge0);

    /// Return half edge (v0,v1) or (v1,v0), if it exists.
    /// - Return NULL if no edge found.
    /// - Internal version that returns non-const half_edge.
    HALF_EDGE_TYPE *
    _FindEdge(const VERTEX_TYPE * v0, const VERTEX_TYPE * v1);

    /// Find some half edge (v0,v1) or (v1,v0) and link with half_edgeA
    ///   in half edge around edge cycle.
    /// - If half edge not found, then do nothing.
    void _FindAndLinkHalfEdgeAroundEdge
    (VERTEX_TYPE * v0, VERTEX_TYPE * v1, HALF_EDGE_TYPE * half_edgeA);

    /// Link half edges around edges merged by merging v0 and v1.
    /// - Searches for all possible merged edges, not just triangles
    ///   containing A and B, to handle non-manifolds.
    /// - Running time: O(v0->NumHalfEdgesFrom() + v1->NumHalfEdgesFrom()).
    void _LinkHalfEdgesAroundMergedEdges
    (VERTEX_TYPE * v0, VERTEX_TYPE * v1);

    /// Relink half edges in cell.
    /// - Overwrites previous links.
    void _RelinkHalfEdgesInCell(HALF_EDGE_TYPE * hprev,
                                HALF_EDGE_TYPE * hnext);

    /// Delete vertex.
    void _DeleteVertex(VERTEX_TYPE * v);

    /// Delete half edge.
    void _DeleteHalfEdge(HALF_EDGE_TYPE * half_edge);

    /// Delete half edges around edge.
    /// @param max_numh Upper bound on the number of edges
    ///    in half edges around edge cycle containing half_edge0.
    /// - max_numh bounds number of iterations in function as 
    ///   a guard against infinite loops if the data structure
    ///   is corrupted.
    void _DeleteHalfEdgesAroundEdge(HALF_EDGE_TYPE * half_edge0,
                                    const int max_numh);

    /// Delete duplicate half edges from cycle
    ///   of half edges around edge.
    /// - Triangle collapses can cause two half edges 
    ///     in the same cell to merge.
    /// - Such edges need to be deleted.
    void _RemoveDuplicatesFromHalfEdgeAroundEdgeCycle
    (HALF_EDGE_TYPE * half_edge0);


    /// Delete cell.
    void _DeleteCell(CELL_TYPE * cell);


  public:


    // *** Get functions ***

    /// Dimension of coordinates.
    int Dimension() const { return(DIM); }    


    // *** Find functions ***

    /// Return half edge (v0,v1) or (v1,v0), if it exists.
    /// - Return NULL if no edge found.
    const HALF_EDGE_TYPE *
    FindEdge(const VERTEX_TYPE * v0, const VERTEX_TYPE * v1) const;

    // *** Functions to modify the mesh ***

    /// Collapse edge mapping two vertices to a single vertex.
    /// - New vertex is at edge midpoint.
    /// - Returns pointer to the new vertex.
    /// - Returns NULL if collapse is illegal. 
    ///   (See IsIllegalEdgeCollapse().)
    /// - Note: This routine DOES NOT check/ensure that mesh
    ///   will be a manifold after the edge collapse.
    const VERTEX_TYPE * CollapseEdge(const int ihalf_edge0);

    /// Split cell with diagonal connecting the two from vertices.
    /// - Diagonal (half_edgeA->FromVertex(), half_edgeB->FromVertex()).
    /// - Returns the half edge h forming the split diagonal
    ///   where h->NextEdgeInCell()->ToVertex() equals 
    ///   half_edgeA->ToVertex().
    /// - Check if diagonal is already a mesh edge, before splitting cell.
    ///   If diagonal is already a mesh edge, splitting will create an
    ///   edge incident on three or more cells.
    /// @pre half_edgeA->Cell() == half_edgeB->Cell().
    /// @pre half_edgeA->ToVertex() != half_edgeB->FromVertex() and
    ///      half_edgeA->FromVertex() != half_edgeB->ToVertex().
    const HALF_EDGE_TYPE * 
    SplitCell(const int ihalf_edgeA, const int ihalf_edgeB);

    /// Join two cells sharing an edge.
    /// - Returns edge incident on the joined cell.
    /// @pre Exactly two cells share the edge indicated by half_edgeA.
    /// - Returns NULL if join fails.
    const HALF_EDGE_TYPE * JoinTwoCells(const int ihalf_edgeA);

    /// Split an edge.
    /// - Returns new vertex.
    const VERTEX_TYPE * SplitEdge(const int ihalf_edgeA);


    // *** Functions to check potential edge collapses. ***

    /// Return true if half edge endpoints and v are in a mesh triangle.
    bool IsInTriangle(const HALF_EDGE_TYPE * half_edge0,
                      const VERTEX_TYPE * v) const;

    /// Return true if cell icell is in the boundary of a tetrahedron.
    bool IsInTetrahedron(const int icell0) const;

    /// Return true if both endpoints (vfrom, vto) of half_edge
    ///   are neighbors of some vertex vC, but (vfrom, vto, vC)
    ///   is not a mesh triangle.
    /// - @param[out] ivC Index of vertex vC, where vC is a neighbor
    ///   of vfrom and vto, but (vftom, vto, vC) is not a mesh triangle.
    /// - Running time is: 
    ///     O(vfrom->NumHalfEdgesFrom() + vto->NumHalfEdgesTo()).
    bool FindTriangleHole
    (const HALF_EDGE_TYPE * half_edge, int & ivC) const;


    /// Return true if cell icell is a triangle whose 3 edges
    ///   are boundary edges.
    bool IsIsolatedTriangle(const int icell) const;

    /// Count number of vertices shared by two cells.
    int CountNumVerticesSharedByTwoCells
    (const CELL_TYPE * cellA, const CELL_TYPE * cellB) const;

    /// Return true if edge collapse is illegal.
    /// - Edge collapse (vA,vB) is illegal if some cell contains
    ///   both vA and vB but not edge (vA,vB).
    /// - CollapseEdge will not perform an illegal collapse.
    /// - Note: This is only possible if the cell has four or more edges.
    bool IsIllegalEdgeCollapse(const VERTEX_TYPE * vA,
                               const VERTEX_TYPE * vB) const;

    /// Return true if edge collapse is illegal.
    bool IsIllegalEdgeCollapse
    (const HALF_EDGE_TYPE * half_edge) const
    { return IsIllegalEdgeCollapse(half_edge->FromVertex(),
                                   half_edge->ToVertex()); }

    /// Return true if split cell is illegal.
    /// - Split cell is illegal 
    ///   if half_edgeA and half_edgeB are in different cells or 
    ///   if half_edgeA->FromVertex and half_edgeB->FromVertex
    ///     are adjacent vertices.
    bool IsIllegalSplitCell
    (const HALF_EDGE_TYPE * half_edgeA,
     const HALF_EDGE_TYPE * half_edgeB) const;

    /// Return true if join cells is illegal.
    /// - Join cells is illegal if half_edge is a boundary half edge
    ///   or some endpoints of half edge has degree 2.
    bool IsIllegalJoinCells(const HALF_EDGE_TYPE * half_edge) const;


    // *** Functions to compute mesh information. ***

    /// Compute min and max squared edge lengths in the mesh.
    /// @param[out] min_edge_length_squared Square of minimum edge length.
    /// @param[out] max_edge_length_squared Square of maximum edge length.
    /// @param[out] ihalf_edge_min Index of half edge with squared length
    ///    equal to min_edge_length_squared.
    /// @param[out] ihalf_edge_max Index of half edge with squared length
    ///    equal to max_edge_length_squared.
    template <typename CTYPE>
    void ComputeMinMaxEdgeLengthSquared
    (CTYPE & min_edge_length_squared, CTYPE & max_edge_length_squared,
     int & ihalf_edge_min, int & ihalf_edge_max) const;

    /// Compute min and max squared edge lengths in the mesh.
    /// - Version that only returns min and max edge length squared.
    template <typename CTYPE>
    void ComputeMinMaxEdgeLengthSquared
    (CTYPE & min_edge_length_squared, 
     CTYPE & max_edge_length_squared) const
    { 
      int ihalf_edge_min, ihalf_edge_max;
      ComputeMinMaxEdgeLengthSquared
        (min_edge_length_squared, max_edge_length_squared,
         ihalf_edge_min, ihalf_edge_max);
    };


    /// Compute min squared ratio of the min to max edge in any cell.
    /// - Ignores cells with all edge lengths 0.
    /// - Returns 1.0 if there are no cells or all edges are length 0.
    template <typename CTYPE>
    void ComputeMinCellEdgeLengthRatioSquared
    (CTYPE & min_edge_length_ratio_squared,
     int & icell_min_ratio,
     CTYPE & min_edge_length_squared, CTYPE & max_edge_length_squared,
     int & ihalf_edge_min, int & ihalf_edge_max) const;

    /// Compute min squared ratio of the min to max edge in any cell.
    /// - Version that only returns min squared ratio and icell.
    template <typename CTYPE>
    void ComputeMinCellEdgeLengthRatioSquared
    (CTYPE & min_edge_length_ratio_squared,
     int & icell_min_ratio) const;

    /// Compute min squared ratio of the min to max edge in any cell.
    /// - Version that only returns min squared ratio.
    template <typename CTYPE>
    void ComputeMinCellEdgeLengthRatioSquared
    (CTYPE & min_edge_length_ratio_squared) const;

    /// Compute cos angle at v1 of triangle (v0,v1,v2).
    /// - Sets flag_zero to true if v1 has same coordinates as v1 or v2
    ///   (or v1 is very, very close to v1 or v2.)
    double ComputeCosTriangleAngle
    (const VERTEX_TYPE * v0, const VERTEX_TYPE * v1, 
     const VERTEX_TYPE * v2, bool & flag_zero) const;

    /// Compute cosine of min and max cell angles.
    /// @param[out] cos_min_angle Cosine of the min angle.
    /// @param[out] cos_max_angle Cosine of the max angle.
    /// @param[out] ihalf_edge_min Index of half edge whose from vertex
    ///    has the min angle.
    /// @param[out] ihalf_edge_max Index of half edge whose from vertex
    ///    has the max angle.
    template <typename CTYPE>
    void ComputeCosMinMaxAngle
    (CTYPE & cos_min_angle, CTYPE & cos_max_angle,
     int & ihalf_edge_min, int & ihalf_edge_max) const;

    /// Compute cosine of min and max cell angles.
    /// - Version that only returns cosine of min and max angles.
    template <typename CTYPE>
    void ComputeCosMinMaxAngle
    (CTYPE & cos_min_angle, CTYPE & cos_max_angle) const;
  };


  /// Simple HALF_EDGE_MESH_DCMT_BASE.
  class HALF_EDGE_MESH_DCMT_A:
    public HALF_EDGE_MESH_DCMT_BASE
  <3, VERTEX_DCMT_A, HALF_EDGE_DCMT_A, CELL_DCMT_A>
  {
  };


  // ******************************************************************
  // Member functions for join/split/collapse/etc.
  // ******************************************************************

  // Return non constant pointer to previous half edge around edge.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  _HALF_EDGE_TYPE * 
  HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  _FindPrevHalfEdgeAroundEdgeNC(HALF_EDGE_TYPE * half_edge0)
  {
    // Cannot have more than max_num half edges around an edge.
    const int max_numh = 
      half_edge0->FromVertex()->NumHalfEdgesFrom() +
      half_edge0->ToVertex()->NumHalfEdgesFrom();
      
    HALF_EDGE_TYPE * half_edge = half_edge0->next_half_edge_around_edge;
    for (int k = 0; k < max_numh; k++) {

      if (half_edge0 == half_edge->next_half_edge_around_edge)
        { return half_edge; }

      half_edge = half_edge->next_half_edge_around_edge;
    }

    // Should never reach here.  Data structure inconsistency.
    throw SIMPLE_EXCEPTION
      ("Programming error. Unable to find previous half edge around half edge.");
  }


  // Collapse edge mapping two vertices to a single vertex.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  const _VERTEX_TYPE * 
  HALF_EDGE_MESH_DCMT_BASE<DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  CollapseEdge(const int ihalf_edge0)
  {
    HALF_EDGE_TYPE * half_edge0 = this->HalfEdgeNC(ihalf_edge0);
    if (half_edge0 == NULL) {
      // Can't collapse a half edge that doesn't exist.
      return NULL;
    }

    // Don't collapse half_edge0 if it's two endpoints (vA,vB) are in some cell
    //   but edge (vA,vB) is not in the cell.
    if (IsIllegalEdgeCollapse(half_edge0)) 
      { return NULL; }

    const int max_num_half_edges_around_edge = 
      half_edge0->FromVertex()->NumHalfEdgesFrom() +
      half_edge0->ToVertex()->NumHalfEdgesFrom();

    VERTEX_TYPE * v0 = half_edge0->FromVertex();
    VERTEX_TYPE * v1 = half_edge0->ToVertex();

    // Update *.next_half_edge_around_edge links.
    _LinkHalfEdgesAroundMergedEdges(v0, v1);
    
    // Update *.prev_half_edge_in_cell and *.next_half_edge_in_cell links.
    // Delete triangle cells containing edge (v0,v1).
    HALF_EDGE_TYPE * half_edge = half_edge0;
    int k = 0;
    do {
      CELL_TYPE * cell = half_edge->Cell();
      HALF_EDGE_TYPE * prev_half_edge = 
        half_edge->prev_half_edge_in_cell;
      HALF_EDGE_TYPE * next_half_edge = 
        half_edge->next_half_edge_in_cell;

      if (cell->HalfEdge() == half_edge) 
        { cell->half_edge = next_half_edge; }

      if (cell->NumVertices() == 3) {
        // Cell is a triangle.
        // Collapsing half edge removes cell.
        VERTEX_TYPE * v2 = prev_half_edge->from_vertex;
        _DeleteHalfEdge(prev_half_edge);
        _DeleteHalfEdge(next_half_edge);
        _DeleteCell(cell);
        v2->_MoveBoundaryHalfEdgeToHalfEdgeFrom0();
      }
      else {
        // Link previous and next half edges.
        prev_half_edge->next_half_edge_in_cell = next_half_edge;
        next_half_edge->prev_half_edge_in_cell = prev_half_edge;

        cell->num_vertices--;
      }

      half_edge = half_edge->next_half_edge_around_edge;

      k++;
    } while ((k < max_num_half_edges_around_edge) &&
             (half_edge != half_edge0));

    // Compute an upper bound on the number of half edges 
    // with endpoints (v0,v1).
    const int max_numh0 =
      v0->NumHalfEdgesFrom() + v1->NumHalfEdgesFrom();

    // Move all half edges from v0 to be from v1.
    // - Moves v0->half_edge_list to v1->half_edge_list.
    _MoveVertexHalfEdgeFromList(v0, v1);

    // Set v1 to midpoint of (v0, v1).
    compute_midpoint<DIM>(v0->coord, v1->coord, v1->coord);

    // Delete half edges around edge (v0,v1).
    _DeleteHalfEdgesAroundEdge(half_edge0, max_numh0);

    _DeleteVertex(v0);

    v1->_MoveBoundaryHalfEdgeToHalfEdgeFrom0();

    return v1;
  }


  // Split cell with diagonal connecting the two from vertices.
  // - Diagonal (half_edgeA->FromVertex(), half_edgeB->FromVertex()).
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  const _HALF_EDGE_TYPE * HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  SplitCell(const int ihalf_edgeA, const int ihalf_edgeB)
  {
    HALF_EDGE_TYPE * half_edgeA = this->half_edge_list[ihalf_edgeA];
    HALF_EDGE_TYPE * half_edgeB = this->half_edge_list[ihalf_edgeB];

    if (half_edgeA == NULL || half_edgeB == NULL) {
      throw SIMPLE_EXCEPTION
        ("Programming error. Arguments to HALF_EDGE_MESH_DCMT_BASE::SplitCell are not half edge indices.");
    }

    if (IsIllegalSplitCell(half_edgeA, half_edgeB))
      { return NULL; }

    VERTEX_TYPE * vA = half_edgeA->FromVertex();
    VERTEX_TYPE * vB = half_edgeB->FromVertex();
    CELL_TYPE * cellA = half_edgeA->Cell();
    CELL_TYPE * cellB = this->_AddCell();
    HALF_EDGE_TYPE * diagA = this->_AddHalfEdge();
    HALF_EDGE_TYPE * diagB = this->_AddHalfEdge();

    HALF_EDGE_TYPE * half_edgeC = _FindEdge(vA, vB);

    // Link diagA and diagB around edge.
    diagA->next_half_edge_around_edge = diagB;
    diagB->next_half_edge_around_edge = diagA;

    if (half_edgeC != NULL) {
      // Link half_edge_around_edge cycle of half_edgeC and diagA/diagB.
      _SwapNextHalfEdgeAroundEdge(half_edgeC, diagA);
    }

    // Set cells containing diagA and diagB.
    diagA->cell = cellA;
    diagB->cell = cellB;

    // Set diagA and diagB from vertices.
    diagA->from_vertex = vB;
    diagB->from_vertex = vA;

    // Add diagA and diagB to vertex half_edge_from[] lists.
    diagA->from_vertex->half_edge_from.push_back(diagA);
    diagB->from_vertex->half_edge_from.push_back(diagB);

    // Change cell of half edges from half_edgeB to half_edgeA.
    int numvA = cellA->NumVertices();
    HALF_EDGE_TYPE * half_edge = half_edgeB;
    int k = 0;
    while ((k < numvA) && (half_edge != half_edgeA)) {
      half_edge->cell = cellB;
      half_edge = half_edge->next_half_edge_in_cell;
      k++;
    }

    // Set num_vertices in cellA and cell B.
    cellB->num_vertices = k+1;
    cellA->num_vertices = numvA+1-k;

    // cellB->half_edge.
    cellB->half_edge = half_edgeB;

    // Change cellA->half_edge, if necessary.
    if (cellA->HalfEdge()->Cell() != cellA) 
      { cellA->half_edge = half_edgeA; }

    HALF_EDGE_TYPE * hprevA = half_edgeA->prev_half_edge_in_cell;
    HALF_EDGE_TYPE * hprevB = half_edgeB->prev_half_edge_in_cell;

    // Link half edges in cell.
    _RelinkHalfEdgesInCell(hprevB, diagA);
    _RelinkHalfEdgesInCell(diagA, half_edgeA);
    _RelinkHalfEdgesInCell(hprevA, diagB);
    _RelinkHalfEdgesInCell(diagB, half_edgeB);

    // Swap first and last edges in half_edge_list[], if necessary.
    // diagA and diagB are not boundary edges, but
    //   diagA->PrevEdgeInCell() or diagB->PrevEdgeInCell() could
    //   be boundary edges.
    diagA->from_vertex->_ProcessFirstLastHalfEdgesFrom();
    diagB->from_vertex->_ProcessFirstLastHalfEdgesFrom();

    return diagA;
  }


  // Join two cells sharing an edge.
  // - Returns edge incident on the joined cell.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  const _HALF_EDGE_TYPE * HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  JoinTwoCells(const int ihalf_edgeA)
  {
    HALF_EDGE_TYPE * half_edgeA = this->half_edge_list[ihalf_edgeA];

    if (half_edgeA == NULL) {
      throw SIMPLE_EXCEPTION
        ("Programming error. Argument to HALF_EDGE_MESH_DCMT_BASE::JoinTwoCells is not a half edge index.");
    }

    if (IsIllegalJoinCells(half_edgeA)) 
      { return NULL; }

    HALF_EDGE_TYPE * half_edgeB = half_edgeA->NextHalfEdgeAroundEdge();

    if (half_edgeB->NextHalfEdgeAroundEdge() != half_edgeA) {
      throw SIMPLE_EXCEPTION
        ("Programming error. Half edge passed to HALF_EDGE_MESH_DCMT_BASE::JoinTwoCells is in an edge shared by three or more cells.");
    }

    VERTEX_TYPE * vA = half_edgeA->FromVertex();
    VERTEX_TYPE * vB = half_edgeB->FromVertex();
    CELL_TYPE * cellA = half_edgeA->Cell();
    CELL_TYPE * cellB = half_edgeB->Cell();
    const int numvA = cellA->NumVertices();
    const int numvB = cellB->NumVertices();
    HALF_EDGE_TYPE * prevA = half_edgeA->PrevHalfEdgeInCell();
    HALF_EDGE_TYPE * nextA = half_edgeA->NextHalfEdgeInCell();
    HALF_EDGE_TYPE * prevB = half_edgeB->PrevHalfEdgeInCell();
    HALF_EDGE_TYPE * nextB = half_edgeB->NextHalfEdgeInCell();

    if (!vA->IsIncidentOnMoreThanTwoEdges() ||
        !vB->IsIncidentOnMoreThanTwoEdges()) {
      // Can't remove an edge if some endpoint only has degree 2.
      return NULL; 
    }

    // Change cellA->HalfEdge() if necessary.
    if (cellA->HalfEdge() == half_edgeA)
      { cellA->half_edge = half_edgeA->NextHalfEdgeInCell(); }

    // Change edges in cellB to be in cellA.
    HALF_EDGE_TYPE * half_edge = half_edgeB->NextHalfEdgeInCell();
    for (int k = 0; k+1 < numvB; k++) {
      half_edge->cell = cellA; 
      half_edge = half_edge->NextHalfEdgeInCell();
    }

    // Set number of vertices in cellA.
    cellA->num_vertices = numvA + numvB - 2;

    // Relink half edges in cell.
    _RelinkHalfEdgesInCell(prevA, nextB);
    _RelinkHalfEdgesInCell(prevB, nextA);

    // Delete cellB and half_edgeA and half_edgeB.
    _DeleteHalfEdge(half_edgeA);
    _DeleteHalfEdge(half_edgeB);
    _DeleteCell(cellB);

    // half_edgeA and half_edgeB are not boundary edges,
    //   but the previous edges in the cell might be boundary edges.
    vA->_MoveBoundaryHalfEdgeToHalfEdgeFrom0();
    vB->_MoveBoundaryHalfEdgeToHalfEdgeFrom0();

    return nextA;
  }


  // Split an edge.
  // - Returns new vertex.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  const _VERTEX_TYPE * HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  SplitEdge(const int ihalf_edgeA)
  {
    HALF_EDGE_TYPE * half_edgeA = this->half_edge_list[ihalf_edgeA];

    if (half_edgeA == NULL) {
      throw SIMPLE_EXCEPTION
        ("Programming error. Argument to HALF_EDGE_MESH_DCMT_BASE::SplitEdge is not a half edge index.");
    }

    if (half_edgeA->IsBoundary()) {
      return _SplitBoundaryEdge(half_edgeA);
    }
    else {
      return _SplitInternalEdge(half_edgeA);
    }

  }


  // ******************************************************************
  // Find member functions
  // ******************************************************************

  // Return half edge (v0,v1) or (v1,v0), if it exists.
  // - Return NULL if no edge found.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  const _HALF_EDGE_TYPE * HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  FindEdge(const VERTEX_TYPE * v0, const VERTEX_TYPE * v1) const
  {
    const HALF_EDGE_TYPE * half_edge =
      v0->FindIncidentHalfEdge(v1->Index());

    if (half_edge != NULL) { return half_edge; }

    half_edge = v1->FindIncidentHalfEdge(v0->Index());

    return half_edge;
  }


  // ******************************************************************
  // Member functions for checking potential edge collapses.
  // ******************************************************************

  // Return true if half edge endpoints and v are in a mesh triangle.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  bool HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  IsInTriangle(const HALF_EDGE_TYPE * half_edge0,
               const VERTEX_TYPE * v) const
  {
    // Cannot have more than max_num half edges around an edge.
    const int max_numh = 
      half_edge0->FromVertex()->NumHalfEdgesFrom() +
      half_edge0->ToVertex()->NumHalfEdgesFrom();

    const HALF_EDGE_TYPE * half_edge = half_edge0;
    int k = 0;
    do {
      if (half_edge->Cell()->IsTriangle()) {
        const HALF_EDGE_TYPE * prev_half_edge =
          half_edge->PrevHalfEdgeInCell();

        if (prev_half_edge->FromVertex() == v)
          { return true; }
      }

      half_edge = half_edge->NextHalfEdgeAroundEdge();
      k++;
    } while ((k < max_numh) && (half_edge != half_edge0));

    return false;
  }


  // Return true if both endpoints (vfrom, vto) of half_edge
  //   are neighbors of some vertex vC, but (vfrom, vto, vC)
  //   is not a mesh triangle.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  bool HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  FindTriangleHole
  (const HALF_EDGE_TYPE * half_edge, int & ivC) const
  {
    // Initialize.
    ivC = 0;

    const VERTEX_TYPE * vfrom = half_edge->FromVertex();
    const VERTEX_TYPE * vto = half_edge->ToVertex();
    vto->_ClearVisitedFlagsInAdjacentVertices();
    vfrom->_SetVisitedFlagsInAdjacentVertices(true);
    
    for (int k = 0; k < vto->NumHalfEdgesFrom(); k++) {
      const HALF_EDGE_TYPE * half_edgeA = vto->KthHalfEdgeFrom(k);
      const VERTEX_TYPE * vtoA = half_edgeA->ToVertex();

      if (vtoA->visited_flag) {
        // vtoA is a neighbor of vfrom and vto.
        if (!IsInTriangle(half_edge, half_edgeA->ToVertex())) {
          ivC = half_edgeA->ToVertex()->Index();
          return true;
        }
      }

      const HALF_EDGE_TYPE * half_edgeB = 
        half_edgeA = half_edgeA->PrevHalfEdgeInCell();
      const VERTEX_TYPE * vfromB = half_edgeB->FromVertex();

      // Check vfromB to handle boundary edges and/or cells 
      //   with arbitrary orientations.
      if (vfromB->visited_flag) {
        // vfromB is a neighbor of vfrom and vto.
        if (!IsInTriangle(half_edge, vfromB)) {
          ivC = half_edgeB->FromVertex()->Index();
          return true;
        }
      }
    }

    return false;
  }


  /// Return true if cell icell is a triangle whose 3 edges
  ///   are boundary edges.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  bool HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  IsIsolatedTriangle(const int icell) const
  {
    const int THREE(3);

    const CELL_TYPE * cell = this->Cell(icell);
    if (cell == NULL) { return false; }

    if (!cell->IsTriangle())
      { return false; }

    const HALF_EDGE_TYPE * half_edge = cell->HalfEdge();
    for (int i = 0; i < THREE; i++) {
      if (!half_edge->IsBoundary()) 
        { return false; }

      half_edge = half_edge->NextHalfEdgeInCell();
    }

    // Cell has three vertices (and three edges) and all edges
    //   are boundary edges.
    return true;
  }


  /// Return true if cell icell is in the boundary of a tetrahedron.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  bool HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  IsInTetrahedron(const int icell0) const
  {
    const CELL_TYPE * cell0 = this->Cell(icell0);
    if (cell0 == NULL) { return false; }

    if (!cell0->IsTriangle())
      { return false; }

    const HALF_EDGE_TYPE * half_edge0 = cell0->HalfEdge();
    const VERTEX_TYPE * iv2 = 
      half_edge0->PrevHalfEdgeInCell()->FromVertex();

    // Cannot have more than max_num half edges around an edge.
    const int max_numh = 
      half_edge0->FromVertex()->NumHalfEdgesFrom() +
      half_edge0->ToVertex()->NumHalfEdgesFrom();
      
    HALF_EDGE_TYPE * half_edge = half_edge0->NextHalfEdgeAroundEdge();
    int k = 0;
    while (k < max_numh && half_edge != half_edge0) {
      const CELL_TYPE * cell = half_edge->Cell();

      if (cell->IsTriangle()) {

        const HALF_EDGE_TYPE * prev_half_edge = 
          half_edge->PrevHalfEdgeInCell();

        const HALF_EDGE_TYPE * next_half_edge = 
          half_edge->NextHalfEdgeInCell();

        if (IsInTriangle(prev_half_edge, iv2) &&
            IsInTriangle(next_half_edge, iv2)) {
          // cell0, cell, and two triangles form a tetrahedron.
          return true;
        }
      }

      k = k+1;
      half_edge = half_edge->NextHalfEdgeAroundEdge();
    }

    return false;
  }


  // Count number of vertices shared by two cells.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  int HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  CountNumVerticesSharedByTwoCells
  (const CELL_TYPE * cellA, const CELL_TYPE * cellB) const
  {
    int num_shared_vertices = 0;

    cellB->_ClearVisitedFlagsInAllVertices();
    cellA->_SetVisitedFlagsInAllVertices(true);

    const HALF_EDGE_TYPE * half_edgeB = cellB->HalfEdge();
    for (int kB = 0; kB < cellB->NumVertices(); kB++) {

      const VERTEX_TYPE * v = half_edgeB->FromVertex();
      if (v->_IsVisited())
        { num_shared_vertices++; }

      half_edgeB = half_edgeB->NextHalfEdgeInCell();
    }

    return num_shared_vertices;
  }


  // Return true if edge collapse is illegal.
  // - Edge collapse (vA,vB) is illegal if some cell contains
  //   both vA and vB but not edge (vA,vB).
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  bool HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  IsIllegalEdgeCollapse(const VERTEX_TYPE * vA,
                        const VERTEX_TYPE * vB) const
  {
    if (vA->NumHalfEdgesFrom() > vB->NumHalfEdgesFrom()) {
      // Swap vA and vB to reduce number of cells processed.
      return IsIllegalEdgeCollapse(vB, vA); 
    }
    else {

      for (int k = 0; k < vA->NumHalfEdgesFrom(); k++) {
        const HALF_EDGE_TYPE * half_edge0 = vA->KthHalfEdgeFrom(k);
        const CELL_TYPE * cell = half_edge0->Cell();
        if (cell->NumVertices() < 4) { 
          // All pairs of cell vertices form an edge.
          continue;
        }

        const HALF_EDGE_TYPE * half_edge = 
          (half_edge0->NextHalfEdgeInCell())->NextHalfEdgeInCell();
        for (int i = 2; i <= cell->NumVertices()-2; i++) {
          if (half_edge->FromVertex() == vB)
            { return true; }
        }
      }

      return false;
    }

  }


  // Return true if split cell is illegal.
  // - Split cell is illegal 
  //   if half_edgeA and half_edgeB are in different cells or 
  //   if half_edgeA->FromVertex and half_edgeB->FromVertex
  //     are adjacent vertices.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  bool HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  IsIllegalSplitCell
  (const HALF_EDGE_TYPE * half_edgeA,
   const HALF_EDGE_TYPE * half_edgeB) const
  {
    if (half_edgeA->Cell() != half_edgeB->Cell())
      { return true; }

    if (half_edgeA == half_edgeB)
      { return true; }

    if (half_edgeA->FromVertex() == half_edgeB->ToVertex())
      { return true; }

    if (half_edgeA->ToVertex() == half_edgeB->FromVertex())
      { return true; }

    return false;
  }


  // Return true if join cells is illegal.
  // - Join cells is illegal if half_edge is a boundary half edge
  //   or more than two cells are incident on the edge
  //   or some endpoints of half edge has degree 2.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  bool HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  IsIllegalJoinCells(const HALF_EDGE_TYPE * half_edge) const
  {
    const int TWO(2);

    if (half_edge->IsBoundary()) { return true; }

    if (!half_edge->FromVertex()->IsIncidentOnMoreThanTwoEdges())
      { return true; }

    if (!half_edge->ToVertex()->IsIncidentOnMoreThanTwoEdges())
      { return true; }

    const HALF_EDGE_TYPE * half_edgeX = 
      half_edge->NextHalfEdgeAroundEdge();

    if (half_edge != half_edgeX->NextHalfEdgeAroundEdge()) {
      // More than two cells are incident on edge
      //   (half_edge->FromVertex(), half_edge->ToVertex()).
      return true;
    }

    if (CountNumVerticesSharedByTwoCells
        (half_edge->Cell(), half_edgeX->Cell()) > TWO) {
      // Cells share more than two vertices.
      return true;
    }

    return false;
  }


  // ******************************************************************
  // Member functions for computing mesh information.
  // ******************************************************************

  // Compute min and max squared edge lengths in the mesh.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  template <typename CTYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  ComputeMinMaxEdgeLengthSquared
  (CTYPE & min_edge_length_squared, CTYPE & max_edge_length_squared,
   int & ihalf_edge_min, int & ihalf_edge_max) const
  {
    bool flag_found = false;

    // Initialize
    min_edge_length_squared = 0.0;
    max_edge_length_squared = 0.0;
    ihalf_edge_min = 0;
    ihalf_edge_max = 0;

    for (int ihalf_edge = 0; ihalf_edge < this->HalfEdgeListLength();
         ihalf_edge++) {
      const HALF_EDGE_TYPE * half_edge = this->HalfEdge(ihalf_edge);
      if (half_edge == NULL) { continue; }

      const CTYPE length_squared = half_edge->ComputeLengthSquared();
      if (!flag_found || (length_squared < min_edge_length_squared)) {
        min_edge_length_squared = length_squared;
        ihalf_edge_min = half_edge->Index();
      }

      if (!flag_found || (length_squared > max_edge_length_squared)) {
        max_edge_length_squared = length_squared;
        ihalf_edge_max = half_edge->Index();
      }

      flag_found = true;
    }
  }


  // Compute min squared ratio of the min to max edge in any cell.
  // - Ignores cells with all edge lengths 0.
  // - Returns 1.0 if there are no cells or all edges are length 0.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  template <typename CTYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  ComputeMinCellEdgeLengthRatioSquared
  (CTYPE & min_edge_length_ratio_squared,
   int & icell_min_ratio,
   CTYPE & min_edge_length_squared, CTYPE & max_edge_length_squared,
   int & ihalf_edge_min, int & ihalf_edge_max) const
  {
    CTYPE min_Lsquared, max_Lsquared;
    int ihalf_min, ihalf_max;

    // Initialize
    min_edge_length_ratio_squared = 1.0;
    icell_min_ratio = 0;
    min_edge_length_squared = 0.0;
    max_edge_length_squared = 0.0;
    ihalf_edge_min = 0;
    ihalf_edge_max = 0;

    for (int icell = 0; icell < this->CellListLength(); icell++) {
      const CELL_TYPE * cell = this->Cell(icell);
      if (cell == NULL) { continue; }

      cell->ComputeMinMaxEdgeLengthSquared
        (min_Lsquared, max_Lsquared, ihalf_min, ihalf_max);

      if (max_Lsquared == 0) {
        // Skip.
        continue;
      }

      CTYPE ratio = min_Lsquared/max_Lsquared;
      if (ratio < min_edge_length_ratio_squared) {
        min_edge_length_ratio_squared = ratio;
        icell_min_ratio = icell;
        min_edge_length_squared = min_Lsquared;
        max_edge_length_squared = max_Lsquared;
        ihalf_edge_min = ihalf_min;
        ihalf_edge_max = ihalf_max;
      }
      
    }
  }


  // Compute min squared ratio of the min to max edge in any cell.
  // - Version that only returns min squared ratio and icell.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  template <typename CTYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  ComputeMinCellEdgeLengthRatioSquared
  (CTYPE & min_edge_length_ratio_squared,
   int & icell_min_ratio) const
  {
    CTYPE min_edge_length_squared, max_edge_length_squared;
    int ihalf_edge_min, ihalf_edge_max;

    ComputeMinCellEdgeLengthRatioSquared
      (min_edge_length_ratio_squared, icell_min_ratio,
       min_edge_length_squared, max_edge_length_squared,
       ihalf_edge_min, ihalf_edge_max);
  }


  // Compute min squared ratio of the min to max edge in any cell.
  // - Version that only returns min squared ratio.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  template <typename CTYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  ComputeMinCellEdgeLengthRatioSquared
  (CTYPE & min_edge_length_ratio_squared) const
  {
    int icell;

    ComputeMinCellEdgeLengthRatioSquared
      (min_edge_length_ratio_squared, icell);
  }


  // Compute cos angle at v1 of triangle (v0,v1,v2).
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  double HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  ComputeCosTriangleAngle
  (const VERTEX_TYPE * v0, const VERTEX_TYPE * v1, 
   const VERTEX_TYPE * v2, bool & flag_zero) const
  { 
    return compute_cos_triangle_angle
      (v0->coord, v1->coord, v2->coord, flag_zero);
  }


  // Compute cosine of min and max cell angles.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  template <typename CTYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  ComputeCosMinMaxAngle
  (CTYPE & cos_min_angle, CTYPE & cos_max_angle,
   int & ihalf_edge_min, int & ihalf_edge_max) const
  {
    CTYPE cos_minA, cos_maxA;
    int ihalf_min, ihalf_max;

    // Initialize
    cos_min_angle = 0.0;
    cos_max_angle = 0.0;
    ihalf_edge_min = 0;
    ihalf_edge_max = 0;

    bool flag_found = false;
    for (int icell = 0; icell < this->CellListLength(); icell++) {
      const CELL_TYPE * cell = this->Cell(icell);
      if (cell == NULL) { continue; }

      cell->ComputeCosMinMaxAngle(cos_minA, cos_maxA, ihalf_min, ihalf_max);

      if (!flag_found || cos_minA > cos_min_angle) {
        cos_min_angle = cos_minA;
        ihalf_edge_min = ihalf_min;
      }

      if (!flag_found || cos_maxA < cos_max_angle) {
        cos_max_angle = cos_maxA;
        ihalf_edge_max = ihalf_max;
      }

      flag_found = true;
    }
  }


  // Compute cosine of min and max cell angles.
  // - Version that only returns cosine of min and max angles.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  template <typename CTYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  ComputeCosMinMaxAngle
  (CTYPE & cos_min_angle, CTYPE & cos_max_angle) const
  {
    int ihalf_edge_min, ihalf_edge_max;

    ComputeCosMinMaxAngle
      (cos_min_angle, cos_max_angle, ihalf_edge_min, ihalf_edge_max);
  }


  // ******************************************************************
  // Internal split functions.
  // ******************************************************************

  // Split an internal edge.
  // - Returns new vertex.
  // @pre half_edgeA is an internal edge.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  const _VERTEX_TYPE * HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  _SplitInternalEdge(HALF_EDGE_TYPE * half_edgeA)
  {
    if (half_edgeA == NULL) {
      throw SIMPLE_EXCEPTION
        ("Programming error. Argument to HALF_EDGE_MESH_DCMT_BASE::_SplitInternalEdge is not a half edge index.");
    }

    HALF_EDGE_TYPE * half_edgeB = half_edgeA->NextHalfEdgeAroundEdge();

    if (half_edgeB->NextHalfEdgeAroundEdge() != half_edgeA) {
      throw SIMPLE_EXCEPTION
        ("Programming error. Half edge passed to HALF_EDGE_MESH_DCMT_BASE::_SplitInternalEdge is in an edge shared by three or more cells.");
    }

    if (half_edgeB == half_edgeA) {
      throw SIMPLE_EXCEPTION
        ("Programming error. Half edge passed to HALF_EDGE_MESH_DCMT_BASE::_SplitInternalEdge is a boundary edge. Call _SplitBoundaryEdge().");
    }

    VERTEX_TYPE * vA = half_edgeA->FromVertex();
    VERTEX_TYPE * vB = half_edgeB->FromVertex();
    CELL_TYPE * cellA = half_edgeA->Cell();
    CELL_TYPE * cellB = half_edgeB->Cell();
    const int numvA = cellA->NumVertices();
    const int numvB = cellB->NumVertices();
    HALF_EDGE_TYPE * nextA = half_edgeA->NextHalfEdgeInCell();
    HALF_EDGE_TYPE * nextB = half_edgeB->NextHalfEdgeInCell();

    // Create new vertex
    const int ivnew = this->AddVertex();
    VERTEX_TYPE * newv = this->vertex_list[ivnew];

    if (newv == NULL) {
      throw SIMPLE_EXCEPTION
        ("Error creating new vertex. Out of memory?");
    }

    // Set newv to midpoint of (vA, vB).
    compute_midpoint<DIM>(vA->coord, vB->coord, newv->coord);

    // Create two new half edges.
    HALF_EDGE_TYPE * new_half_edgeA = this->_AddHalfEdge();
    HALF_EDGE_TYPE * new_half_edgeB = this->_AddHalfEdge();
    new_half_edgeA->cell = cellA;
    new_half_edgeB->cell = cellB;
    new_half_edgeA->from_vertex = newv;
    new_half_edgeB->from_vertex = newv;
    newv->half_edge_from.push_back(new_half_edgeA);
    newv->half_edge_from.push_back(new_half_edgeB);

    // Relink half edges in cell.
    _RelinkHalfEdgesInCell(half_edgeA, new_half_edgeA);
    _RelinkHalfEdgesInCell(half_edgeB, new_half_edgeB);
    _RelinkHalfEdgesInCell(new_half_edgeA, nextA);
    _RelinkHalfEdgesInCell(new_half_edgeB, nextB);

    // Increase number of vertices in each cell.
    cellA->num_vertices++;
    cellB->num_vertices++;

    // Unlink half_edgeA->next_half_edge_around_edge and
    // half_edgeB->next_half_edge_around_edge.
    half_edgeA->next_half_edge_around_edge = half_edgeA;
    half_edgeB->next_half_edge_around_edge = half_edgeB;

    // Link half edges around edge.
    this->_LinkHalfEdgesAroundEdge(half_edgeA, new_half_edgeB);
    this->_LinkHalfEdgesAroundEdge(half_edgeB, new_half_edgeA);

    // half_edgeA and half_edgeB are not boundary edges,
    //   but the previous edges in the cell might be boundary edges.
    vA->_MoveBoundaryHalfEdgeToHalfEdgeFrom0();
    vB->_MoveBoundaryHalfEdgeToHalfEdgeFrom0();

    return newv;
  }


  // Split a boundary edge.
  // - Returns new vertex.
  // @pre half_edgeA is a boundary edge.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  const _VERTEX_TYPE * HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  _SplitBoundaryEdge(HALF_EDGE_TYPE * half_edgeA)
  {
    if (half_edgeA == NULL) {
      throw SIMPLE_EXCEPTION
        ("Programming error. Argument to HALF_EDGE_MESH_DCMT_BASE::_SplitBoundaryEdge is not a half edge index.");
    }

    if (!half_edgeA->IsBoundary()) {
      throw SIMPLE_EXCEPTION
        ("Programming error. Half edge passed to HALF_EDGE_MESH_DCMT_BASE::_SplitBoundaryEdge is an internal edge. Call _SplitInternalEdge().");
    }

    VERTEX_TYPE * vA = half_edgeA->FromVertex();
    VERTEX_TYPE * vB = half_edgeA->ToVertex();
    CELL_TYPE * cellA = half_edgeA->Cell();
    const int numvA = cellA->NumVertices();
    HALF_EDGE_TYPE * nextA = half_edgeA->NextHalfEdgeInCell();

    // Create new vertex
    const int ivnew = this->AddVertex();
    VERTEX_TYPE * newv = this->vertex_list[ivnew];

    if (newv == NULL) {
      throw SIMPLE_EXCEPTION
        ("Error creating new vertex. Out of memory?");
    }

    // Set newv to midpoint of (vA, vB).
    compute_midpoint<DIM>(vA->coord, vB->coord, newv->coord);

    // Create two new half edges.
    HALF_EDGE_TYPE * new_half_edgeA = this->_AddHalfEdge();
    new_half_edgeA->cell = cellA;
    new_half_edgeA->from_vertex = newv;
    newv->half_edge_from.push_back(new_half_edgeA);

    // Relink half edges in cell.
    _RelinkHalfEdgesInCell(half_edgeA, new_half_edgeA);
    _RelinkHalfEdgesInCell(new_half_edgeA, nextA);

    // Increase number of vertices in cell.
    cellA->num_vertices++;

    // No need to move edges in half_edge_from[] lists.

    return newv;
  }

  // ******************************************************************
  // Member functions for internal use.
  // ******************************************************************

  // Remove half edge from the half_edge_from list of its from_vertex.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  _RemoveHalfEdgeFromVertexList
  (const HALF_EDGE_TYPE * half_edge0)
  {
    VERTEX_TYPE * v0 = half_edge0->FromVertex();
    const int list_length = v0->NumHalfEdgesFrom();
    const int ilast = list_length-1;

    for (int k = 0; k < list_length; k++) {

      const HALF_EDGE_TYPE * half_edge = v0->KthHalfEdgeFrom(k);

      if (half_edge0 == half_edge) {
        if (k != ilast) {
          // Move half_edge0 to back of lis.
          std::swap(v0->half_edge_from[k], v0->half_edge_from[ilast]);
        }

        v0->half_edge_from.pop_back();

        return;
      }
    }
  }


  // Move half edges in vA->half_edge_from[] to
  //   vB->half_edge_from[].
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  _MoveVertexHalfEdgeFromList
  (VERTEX_TYPE * vA, VERTEX_TYPE *vB)
  {
    for (int k = 0; k < vA->NumHalfEdgesFrom(); k++) {
      HALF_EDGE_TYPE * half_edge = vA->KthHalfEdgeFrom(k);
      half_edge->from_vertex = vB;
    }

    // Append vA->half_edge_from[] to vB->half_edge_from[].
    vB->half_edge_from.insert(vB->half_edge_from.end(),
                              vA->half_edge_from.begin(),
                              vA->half_edge_from.end());

    vA->half_edge_from.clear();
  }


  /// Swap next_half_edge_around_edge.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  _SwapNextHalfEdgeAroundEdge
  (HALF_EDGE_TYPE * half_edgeA, HALF_EDGE_TYPE * half_edgeB)
  {
    HALF_EDGE_TYPE * tempA = half_edgeA->next_half_edge_around_edge;
    HALF_EDGE_TYPE * tempB = half_edgeB->next_half_edge_around_edge;

    half_edgeA->next_half_edge_around_edge = tempB;
    half_edgeB->next_half_edge_around_edge = tempA;
  }


  // Delete vertex.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  _DeleteVertex(VERTEX_TYPE * v)
  {
    if (v == NULL) {
      // Can't delete a NULL vertex.
      return;
    }

    const int iv = v->Index();
    this->vertex_list[iv] = NULL;
    delete v;
  }


  // Delete half edge.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  _DeleteHalfEdge(HALF_EDGE_TYPE * half_edge0)
  {
    if (!half_edge0->IsBoundary()) {
      HALF_EDGE_TYPE * next_half_edge_around_edge =
        half_edge0->next_half_edge_around_edge;

      HALF_EDGE_TYPE * prev_half_edge_around_edge =
        _FindPrevHalfEdgeAroundEdgeNC(half_edge0);

      prev_half_edge_around_edge->next_half_edge_around_edge =
        next_half_edge_around_edge;
    }

    _RemoveHalfEdgeFromVertexList(half_edge0);

    const int ihalf_edge0 = half_edge0->Index();
    this->half_edge_list[ihalf_edge0] = NULL;

    delete half_edge0;
  }


  // Return half edge (v0,v1) or (v1,v0), if it exists.
  // - Return NULL if no edge found.
  // - Internal version that returns non-const half_edge.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  _HALF_EDGE_TYPE * HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  _FindEdge(const VERTEX_TYPE * v0, const VERTEX_TYPE * v1)
  {
    HALF_EDGE_TYPE * half_edge =
      v0->FindIncidentHalfEdge(v1->Index());

    if (half_edge != NULL) { return half_edge; }

    half_edge = v1->FindIncidentHalfEdge(v0->Index());

    return half_edge;
  }


  // Find some half edge (v0,v1) or (v1,v0) and link with half_edgeA
  //   in half edge around edge cycle.
  // - If half edge not found, then do nothing.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  _FindAndLinkHalfEdgeAroundEdge
  (VERTEX_TYPE * v0, VERTEX_TYPE * v1, HALF_EDGE_TYPE * half_edgeA)
  {
    HALF_EDGE_TYPE * half_edgeB = _FindEdge(v0, v1);
    if (half_edgeB != NULL) {

      // Determine if half_edgeB is a boundary edge before link.
      const bool flag_boundary = half_edgeB->IsBoundary();
              
      _SwapNextHalfEdgeAroundEdge(half_edgeA, half_edgeB);

      // half_edgeA is no longer a boundary edge.
      v1->_MoveBoundaryHalfEdgeToHalfEdgeFrom0();
    }
  }


  // Link half edges around edges merged by merging v0 and v1.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  _LinkHalfEdgesAroundMergedEdges(VERTEX_TYPE * v0, VERTEX_TYPE * v1)
  {
    v1->_ClearVisitedFlagsInAdjacentVertices();
    v0->_SetVisitedFlagsInAdjacentVertices(true);
    
    for (int k = 0; k < v1->NumHalfEdgesFrom(); k++) {
      HALF_EDGE_TYPE * half_edgeA = v1->half_edge_from[k];
      VERTEX_TYPE * vtoA = half_edgeA->ToVertex();

      if (vtoA->visited_flag) {
        // vtoA is a neighbor of v0 and v1.

        _FindAndLinkHalfEdgeAroundEdge(v0, vtoA, half_edgeA);

        // Set vtoA->visited_flag to false so that vtoA
        //   will not be processed twice.
        vtoA->visited_flag = false;
      }

      HALF_EDGE_TYPE * half_edgeB =
        half_edgeA->prev_half_edge_in_cell;
      VERTEX_TYPE * vfromB = half_edgeB->FromVertex();

      // Check vfromB to handle boundary edges and/or cells 
      //   with arbitrary orientations.
      if (vfromB->visited_flag) {
        // vfromY is a neighbor of v0 and v1.

        _FindAndLinkHalfEdgeAroundEdge(v0, vfromB, half_edgeB);

        // Set vfromB->visited_flag to false so that vfromB
        //   will not be processed twice.
        vfromB->visited_flag = false;
      }
    }

  }


  // Relink half edges in cell.
  // - Overwrites previous links.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  _RelinkHalfEdgesInCell(HALF_EDGE_TYPE * hprev,
                         HALF_EDGE_TYPE * hnext)
  {
    hprev->next_half_edge_in_cell = hnext;
    hnext->prev_half_edge_in_cell = hprev;
  }


  // Delete half edges around edge.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  _DeleteHalfEdgesAroundEdge(HALF_EDGE_TYPE * half_edge0,
                             const int max_numh)
  {
    HALF_EDGE_TYPE * half_edge = half_edge0;

    for (int k = 0; k < max_numh; k++) {

      HALF_EDGE_TYPE * next_half_edge_around_edge =
        half_edge->next_half_edge_around_edge;
      
      if (next_half_edge_around_edge == half_edge) {
        // Delete half edge.
        const int ihalf_edge = half_edge->Index();
        this->half_edge_list[ihalf_edge] = NULL;
        _RemoveHalfEdgeFromVertexList(half_edge);
        delete half_edge;

        return;
      }
      else {
        // Delete next_half_edge_around_edge.
        half_edge->next_half_edge_around_edge =
          next_half_edge_around_edge->next_half_edge_around_edge;
        const int inext_half_edge = next_half_edge_around_edge->Index();
        this->half_edge_list[inext_half_edge] = NULL;
        _RemoveHalfEdgeFromVertexList(next_half_edge_around_edge);
        delete next_half_edge_around_edge;
      }
    }
  }


  /// Delete cell.
  template <const int DIM, typename _VERTEX_TYPE, 
            typename _HALF_EDGE_TYPE, typename _CELL_TYPE>
  void HALF_EDGE_MESH_DCMT_BASE
  <DIM, _VERTEX_TYPE, _HALF_EDGE_TYPE, _CELL_TYPE>::
  _DeleteCell(CELL_TYPE * cell)
  {
    const int icell = cell->Index();
    this->cell_list[icell] = NULL;

    delete cell;
  }


  // *****************************************************************
  // Member functions of VERTEX_DCMT_BASE.
  // *****************************************************************


  // Return true if vertex is incident on more than three edges.
  template <const int DIM, typename HALF_EDGE_PTR, typename CTYPE>
  bool VERTEX_DCMT_BASE<DIM,HALF_EDGE_PTR,CTYPE>::
  IsIncidentOnMoreThanTwoEdges() const
  {
    const int TWO(2);

    const int num_half_edges_from = this->NumHalfEdgesFrom();

    if (num_half_edges_from > TWO) { return true; }

    if (!this->IsBoundary()) { return false; }

    if (num_half_edges_from == TWO) { 
      // Boundary vertex in at two cells must have at least three
      //   incident edges.
      return true;
    }
    else {
      // Boundary vertex is in just one cell, and has exactly two
      //   incident edges.
      return false;
    }
  }


  // Set visited_flag to flag in all neighbors of this.
  template <const int DIM, typename HALF_EDGE_PTR, typename CTYPE>
  void VERTEX_DCMT_BASE<DIM,HALF_EDGE_PTR,CTYPE>::
  _SetVisitedFlagsInAdjacentVertices(const bool flag) const
  {
    for (int k = 0; k < this->NumHalfEdgesFrom(); k++) {
      const HALF_EDGE_PTR half_edgeA = this->KthHalfEdgeFrom(k);

      half_edgeA->ToVertex()->visited_flag = flag;

      const HALF_EDGE_PTR half_edgeB = 
        half_edgeA->PrevHalfEdgeInCell();

      // Set half_edgeB->FromVertex()->visited_flag in case of
      ///  boundary edges or cells with arbitrary orientations.
      half_edgeB->FromVertex()->visited_flag = flag;
    }
  }


  // Process first and last half edges in half_edges_from[].
  template <const int DIM, typename HALF_EDGE_PTR, typename CTYPE>
  void VERTEX_DCMT_BASE<DIM,HALF_EDGE_PTR,CTYPE>::
  _ProcessFirstLastHalfEdgesFrom()
  {
    if (this->NumHalfEdgesFrom() < 2) {
      // No swap
      return;
    }

    if (this->KthHalfEdgeFrom(0)->IsBoundary()) {
      // No swap
      return;
    }

    const int ilast = this->NumHalfEdgesFrom()-1;
    if (this->KthHalfEdgeFrom(ilast)->IsBoundary()) {
      std::swap(this->half_edge_from[0], this->half_edge_from[ilast]);
      return;
    }

    if (this->KthHalfEdgeFrom(0)->PrevHalfEdgeInCell()->IsBoundary()) {
      // No swap.
      return;
    }

    if (this->KthHalfEdgeFrom(ilast)->PrevHalfEdgeInCell()->IsBoundary()) {
      std::swap(this->half_edge_from[0], this->half_edge_from[ilast]);
      return;
    }
    
  }


  // *****************************************************************
  // Member functions of HALF_EDGE_DCMT_BASE.
  // *****************************************************************

  // Compute cosine of angle between PreviousHalfEdge()->FromVertex(),
  //   this->FromVertex() and this->ToVertex().
  // - Sets flag_zero to true if v1 has same coordinates as v0 or v2
  //   (or v1 is very, very close to v0 or v2.)
  template <const int DIM, typename VERTEX_PTR, typename HALF_EDGE_PTR, 
            typename CELL_PTR>
  double HALF_EDGE_DCMT_BASE<DIM,VERTEX_PTR,HALF_EDGE_PTR,CELL_PTR>::
  ComputeCosAngleAtFromVertex(bool & flag_zero) const
  {
    const HALF_EDGE_PTR prev_half_edge =
      this->PrevHalfEdgeInCell();
    const VERTEX_PTR v0 = prev_half_edge->FromVertex();
    const VERTEX_PTR v1 = this->FromVertex();
    const VERTEX_PTR v2 = this->ToVertex();

    return compute_cos_triangle_angle<DIM>
      (v0->coord, v1->coord, v2->coord, flag_zero);
  }


  // *****************************************************************
  // Member functions of CELL_DCMT_BASE.
  // *****************************************************************

  /// Compute max and min squared edge lengths in a cell.
  template <const int DIM, typename HALF_EDGE_PTR>
  template <typename CTYPE>
  void CELL_DCMT_BASE<DIM, HALF_EDGE_PTR>::
  ComputeMinMaxEdgeLengthSquared
  (CTYPE & min_edge_length_squared, CTYPE & max_edge_length_squared,
   int & ihalf_edge_min, int & ihalf_edge_max) const
  {
    if (this->NumVertices() == 0) {
      // Empty cell.  Set values to defaults and return.
      min_edge_length_squared = 0.0;
      max_edge_length_squared = 0.0;
      ihalf_edge_min = 0;
      ihalf_edge_max = 0;
    }

    HALF_EDGE_PTR half_edge = this->HalfEdge();
    min_edge_length_squared = half_edge->ComputeLengthSquared();
    max_edge_length_squared = min_edge_length_squared;
    ihalf_edge_min = half_edge->Index();
    ihalf_edge_max = ihalf_edge_min;

    for (int i = 1; i < this->NumVertices(); i++) {
      half_edge = half_edge->NextHalfEdgeInCell();

      const CTYPE length_squared = half_edge->ComputeLengthSquared();
      if (length_squared < min_edge_length_squared) {
        min_edge_length_squared = length_squared;
        ihalf_edge_min = half_edge->Index();
      }

      if (length_squared > max_edge_length_squared) {
        max_edge_length_squared = length_squared;
        ihalf_edge_max = half_edge->Index();
      }
    }
  }


  // Compute cosine of the min or max angle between consecutive
  //   cell edges.
  // - Note: cos_min_angle >= cos_max_angle.
  // - The smallest angle is 0 and cos(0) = 1.
  // - The largest angle is pi and cos(pi) = -1.
  template <const int DIM, typename HALF_EDGE_PTR>
  template <typename CTYPE>
  void CELL_DCMT_BASE<DIM, HALF_EDGE_PTR>::
  ComputeCosMinMaxAngle
  (CTYPE & cos_min_angle, CTYPE & cos_max_angle,
   int & ihalf_edge_min, int & ihalf_edge_max) const
  {
    bool flag_zero;

    if (this->NumVertices() == 0) {
      // Empty cell.  Set values to defaults and return.
      cos_min_angle = 0.0;
      cos_max_angle = 0.0;
      ihalf_edge_min = 0;
      ihalf_edge_max = 0;
    }

    HALF_EDGE_PTR half_edge = this->HalfEdge();
    
    bool flag_found = false;
    for (int i = 0; i < this->NumVertices(); i++) {

      double cos_angle = half_edge->ComputeCosAngleAtFromVertex(flag_zero);

      if (!flag_zero) {

        if (!flag_found) {
          cos_min_angle = cos_angle;
          cos_max_angle = cos_angle;
          ihalf_edge_min = half_edge->Index();
          ihalf_edge_max = ihalf_edge_min;
        }
        else if (cos_angle > cos_min_angle) {
          // Remember: Small angles have large cos values.
          cos_min_angle = cos_angle;
          ihalf_edge_min = half_edge->Index();
        }
        // Note: cos_min_angle >= cos_max_angle, so 
        //   if cos_angle > cos_min_angle, then 
        //   (cos_angle < cos_max_angle) is false.
        else if (cos_angle < cos_max_angle) {
          // Remember: Large angles have small (possibly negative) 
          //   cos values.
          cos_max_angle = cos_angle;
          ihalf_edge_max = half_edge->Index();
        }

        flag_found = true;
      }
      
      half_edge = half_edge->NextHalfEdgeInCell();
    }
    
    if (!flag_found) {
      // Set to default values.
      cos_min_angle = 0.0;
      cos_max_angle = 0.0;
      ihalf_edge_min = this->HalfEdge()->Index();
      ihalf_edge_max = ihalf_edge_min;
    }
  }


  // Set visited_flag to flag in all cell vertices.
  template <const int DIM, typename HALF_EDGE_PTR>
  void CELL_DCMT_BASE<DIM, HALF_EDGE_PTR>::
  _SetVisitedFlagsInAllVertices(const bool flag) const
  {
    HALF_EDGE_PTR half_edge = this->HalfEdge();
    for (int k = 0; k < this->NumVertices(); k++) {
      half_edge->FromVertex()->visited_flag = flag;
      half_edge = half_edge->NextHalfEdgeInCell();
    }
  }

}

#endif
