/// \file half_edge_mesh.hpp
/// template classes for 2D half edge mesh data structures.
/*! \mainpage Half edge mesh:
 *  This is an implementation of a half edge mesh.
 *
 *  The mesh is stored in HALF_EDGE_MESH_BASE, including lists
 *  of all the vertices, half edges and cells in the mesh.
 *  - Each vertex, half edge and cell is in its own class.
 *  - All allocations of vertices, half edges and cells should be done 
 *    in HALF_EDGE_MESH_BASE or a subclass of HALF_EDGE_MESH_BASE.
 *  - Each vertex, half edge and cell can be identified by a pointer
 *    to the object containing the vertex, half edge or cell, or
 *    by an integer index (identifier) of the vertex, half edge or cell.
 *  - This is NOT a very efficient/compact implementation of half edges. 
 *  - This implementation is meant to be simple and (hopefully) robust
 *    for use in OSU CSE homeworks and prototypes.
 *  - Note: Many of the simpler get functions do not check their arguments,
 *    e.g. NULL pointers or indices in range.
 *     Such checks would be too time consuming for large meshes. 
 *     The calling function is responsible to ensure that pointers are
 *     not NULL and indices are in a specified range.
 * - Requires C++-11 or later.
 * - Version 0.0.1
 */

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

#include<algorithm>
#include<exception>
#include<string>
#include<sstream>
#include<vector>

// *** DEBUG ***
// Uncomment if you want to include print statements 
//   for debugging in this file.
#include<iostream>


#ifndef _HALF_EDGE_MESH_
#define _HALF_EDGE_MESH_

namespace HMESH {

  template <const int DIM, typename HALF_EDGE_PTR, typename CTYPE>
  class VERTEX_BASE {

  protected:
    void Init();

    /// index: Unique non-negative integer identifying the vertex.
    int index;

    /// Store list of all half edges originating at vertex.
    /// - Useful if mesh is not a manifold and cells incident
    ///   on vertex do not form a fan.
    std::vector<HALF_EDGE_PTR> half_edge_from;

    /// Move boundary half edge to half_edge_from[0].
    /// - If there are no boundary half edges in half_edge_from[],
    ///   but half_edge_from[k]->PreviousHalfEdgeInCell() is a
    ///   boundary half edge, move half_edge_from[k] to half_edge_from[0].
    /// - Does nothing if half_edge_from[0] is a boundary half edge
    ///   or if vertex is not incident on any boundary half edges (from or to).
    /// - Revised: 11-24-2021 - RW
    void _MoveBoundaryHalfEdgeToHalfEdgeFrom0();

  public:
    typedef CTYPE COORD_TYPE;

  public:

    /// Vertex coordinates.
    /// - public.  Can be modified outside of HALF_EDGE_MESH 
    ///   data structure.
    CTYPE coord[DIM];

  public:
    /// Constructor.
    VERTEX_BASE() { Init(); }

    int Dimension() const { return(DIM); }

    /// Return coord[ic].
    /// @pre ic is in range [0..(Dimension()-1)].
    COORD_TYPE Coord(const int ic) const { return(coord[ic]); }

    /// Return vertex index.
    int Index() const { return(index); }

    /// Return number of half edges with from vertex 
    ///   equal to this.
    int NumHalfEdgesFrom() const
    { return(half_edge_from.size()); }

    /// Return k'th incident half edge with from vertex equal
    /// @pre k < NumHalfEdgesFrom() and k >= 0.
    const HALF_EDGE_PTR KthHalfEdgeFrom(const int k) const
    { return(half_edge_from[k]); }

    /// Return true if vertex is on the boundary or vertex is not
    ///   incident on any cells.
    bool IsBoundary() const
    { 
      if (NumHalfEdgesFrom() == 0) {
        // Vertex is not incident on any cells.
        return true;
      }
      // Vertex is on the boundary iff the half_edge_from[0].IsBoundary() and
      //   or half_edge_from[0].PrevHalfEdgeInCell().IsBoundary().
      return (KthHalfEdgeFrom(0)->IsBoundary() || 
              KthHalfEdgeFrom(0)->PrevHalfEdgeInCell()->IsBoundary());
    }

    /// Return incident half edge whose FromVertex() is current vertex
    ///   and whose ToVertexIndex() is iv.
    /// - Return NULL if no half edge found.
    const HALF_EDGE_PTR FindIncidentHalfEdge(const int iv) const;

    /// Count number of half edges in incidence list 
    ///   whose from vertex is current verex and whose to vertex is iv.
    int CountNumIncidentHalfEdges(const int iv) const;

    /// Print vertex coordinates.  Separate coordinates with {separator}.
    template <typename OSTREAM_TYPE>
    void PrintCoord(OSTREAM_TYPE & out, 
                    const std::string & separator) const;

    template <typename VTYPE, typename HTYPE, typename CELL_TYPE>
    friend class HALF_EDGE_MESH_BASE;
  };


  template <typename HALF_EDGE_PTR, typename CTYPE>
  class VERTEX3D_BASE:public VERTEX_BASE<3,HALF_EDGE_PTR,CTYPE>
  {};


  template <typename VERTEX_PTR, typename HALF_EDGE_PTR, typename CELL_PTR>
  class HALF_EDGE_BASE {

  protected:

    void Init();

    /// index: Unique non-negative integer identifying the half edge.
    int index;

    /// Pointer to next half edge in cell.
    HALF_EDGE_PTR next_half_edge_in_cell;

    /// Pointer to previous half edge in cell.
    HALF_EDGE_PTR prev_half_edge_in_cell;

    /// Pointer to next half edge around the edge.
    /// - Equivalent to opposite when at most two cells share an edge
    ///   and cells are consistently oriented.
    /// - Equals itself (this) if half edge is a boundary edge.
    HALF_EDGE_PTR next_half_edge_around_edge;

    /// Pointer to from vertex.
    VERTEX_PTR from_vertex;

    /// Pointer to cell containing half edge.
    CELL_PTR cell;


  public:

    /// Constructor.
    HALF_EDGE_BASE() { Init(); }

    /// Return half edge index.
    int Index() const { return(index); }

    /// Return pointer to next half edge in cell.
    const HALF_EDGE_PTR NextHalfEdgeInCell() const
    { return(next_half_edge_in_cell); }

    /// Return pointer to previous half edge in cell.
    const HALF_EDGE_PTR PrevHalfEdgeInCell() const
    { return(prev_half_edge_in_cell); }

    /// Return pointer to next half edge around edge.
    const HALF_EDGE_PTR NextHalfEdgeAroundEdge() const
    { return(next_half_edge_around_edge); }

    /// Return pointer to previous edge around from vertex.
    /// - Returns PrevHalfEdgeInCell() if
    ///   PrevHalfEdgeInCell() is boundary half edge.
    /// - NextHalfEdgeAroundFromVertex() is not defined, since
    ///   PrevHalfEdgeAroundFromVertex() should always be used in moving
    ///   around a vertex.
    const HALF_EDGE_PTR PrevHalfEdgeAroundFromVertex() const
    { return(PrevHalfEdgeInCell()->NextHalfEdgeAroundEdge()); }

    /// Return pointer to previous edge around vertex iv.
    /// - Returns PrevHalfEdgeInCell() if
    ///   PrevHalfEdgeInCell() is boundary half edge.
    /// - Note: If iv == ToVertex(), returns 
    ///    NextHalfEdgeInCell()->NextHalfEdgeAroundEdge() so that repeated calls
    ///    to PrevHalfEdgeAroundVertex() move in a consistent direction. 
    const HALF_EDGE_PTR PrevHalfEdgeAroundVertex(const int iv) const;

    /// Return pointer to cell containing half edge.
    const CELL_PTR Cell() const { return(cell); }

    /// Return index of cell.
    /// - Added: 11-25-2021 - RW
    const int CellIndex() const { return Cell()->Index(); }

    /// Return pointer to from vertex.
    const VERTEX_PTR FromVertex() const
    { return(from_vertex); }

    /// Return index of "from vertex".
    const int FromVertexIndex() const
    { return(FromVertex()->Index()); }

    /// Return pointer to "to vertex".
    const VERTEX_PTR ToVertex() const
    { return(NextHalfEdgeInCell()->FromVertex()); }

    /// Return index of "to vertex".
    const int ToVertexIndex() const
    { return(ToVertex()->Index()); }

    /// Return true if half edge is a boundary half edge.
    bool IsBoundary() const
    { return(this == (this->NextHalfEdgeAroundEdge())); }

    /// Count number of half edges around edge.
    int CountNumHalfEdgesAroundEdge() const;

    /// Return pointer to half edge with minimum index
    ///   in cycle of half edges around edge.
    /// - Useful in getting a single half edge representing an edge.
    /// - Could return pointer to self.
    /// - Added: 11-25-2021 - RW
    const HALF_EDGE_PTR MinIndexHalfEdgeAroundEdge() const;

    /// Return true if half_edgeB has same endpoints as {this}.
    /// - half_edgeB could be oriented in same direction as {this}
    ///   or in opposite direction as {this}.
    bool SameEndpoints(const HALF_EDGE_PTR half_edgeB) const;

    /// Return string of endpoints of half edge.
    /// Separate coordinates with {separator}.
    /// - Added: 11-25-2021 - RW
    std::string EndpointsStr(const std::string & separator) const;

    /// Return string of half edge index and endpoints.
    /// Separate coordinates with {separator}.
    /// - Added: 11-25-2021 - RW
    std::string IndexAndEndpointsStr
    (const std::string & separator) const;

    /// Print endpoints of half edge.
    template <typename OSTREAM_TYPE>
    void PrintEndpoints(OSTREAM_TYPE & out, 
                        const std::string & separator) const;

    /// Print half edge index and endpoints
    template <typename OSTREAM_TYPE>
    void PrintIndexAndEndpoints(OSTREAM_TYPE & out, 
                                const std::string & separator) const;

    template <typename VTYPE, typename HTYPE, typename CELL_TYPE>
    friend class HALF_EDGE_MESH_BASE;
  };


  template <typename HALF_EDGE_PTR>
  class CELL_BASE {

  protected:
    void Init() { half_edge = NULL; num_vertices = 0; };

    /// index: Unique non-negative integer identifying the cell.
    int index;

    /// Some half edge in the cell.
    HALF_EDGE_PTR half_edge;

    /// Number of cell vertices.
    int num_vertices;

  public:
    
    /// Constructor.
    CELL_BASE() { Init(); }

    /// Return cell index
    int Index() const { return(index); }

    /// Return half edge.
    const HALF_EDGE_PTR HalfEdge() const
    { return(half_edge); }

    /// Return number of cell vertices.
    int NumVertices() const
    { return(num_vertices); }

    /// Return true if cell has exactly 3 vertices.
    /// - Added: 11-24-2021 - RW
    bool IsTriangle() const
    { return (NumVertices() == 3); }

    template <typename VTYPE, typename HTYPE, typename CELL_TYPE>
    friend class HALF_EDGE_MESH_BASE;
  };


  // Forward declaration of simple vertex, half_edge and cell types.
  class VERTEX3D_A;
  class HALF_EDGE_A;
  class CELL_A;

  typedef VERTEX3D_A * VERTEX3D_A_PTR;
  typedef HALF_EDGE_A * HALF_EDGE_A_PTR;
  typedef CELL_A * CELL_A_PTR;

  // Declaration of VERTEX3D_A, HALF_EDGE_A, CELL_A.
  class VERTEX3D_A:public VERTEX3D_BASE<HALF_EDGE_A_PTR,float>
  {};

  class HALF_EDGE_A:public HALF_EDGE_BASE
  <VERTEX3D_A_PTR,HALF_EDGE_A_PTR,CELL_A_PTR>
  {};

  class CELL_A:public CELL_BASE<HALF_EDGE_A_PTR>
  {};


  /// Simple exception template
  class SIMPLE_EXCEPTION:public std::exception {

  public:
    const char * message;

  public:
    SIMPLE_EXCEPTION(const char * msg):message(msg){};

    virtual const char * what() const noexcept
    { return message; }
  };


  template <typename _VERTEX_TYPE, typename _HALF_EDGE_TYPE,
            typename _CELL_TYPE>
  class HALF_EDGE_MESH_BASE {

  public:
    typedef _VERTEX_TYPE VERTEX_TYPE;
    typedef _HALF_EDGE_TYPE HALF_EDGE_TYPE;
    typedef _CELL_TYPE CELL_TYPE;
    typedef typename VERTEX_TYPE::COORD_TYPE COORD_TYPE;

  protected:
    /// List of pointers to vertices.  Some elements could be NULL.
    /// - Should be using unique pointers, but...
    std::vector<VERTEX_TYPE *> vertex_list;

    /// List of pointers to hald edges.  Some elements could be NULL.
    /// - Should be using unique pointers, but...
    std::vector<HALF_EDGE_TYPE *> half_edge_list;

    /// List of pointers to cells.  Some elements could be NULL.
    /// - Should be using unique pointers, but...
    std::vector<CELL_TYPE *> cell_list;

    /// Create vertex with index iv, if it does not yet exist.
    /// - Returns pointer to vertex.
    /// - Returns pointer to vertex, if vertex already exists.
    VERTEX_TYPE * _CreateVertex(const int iv);

    /// Add half edge to half_edge_list[].
    /// - Returns pointer to new half edge.
    HALF_EDGE_TYPE * _AddHalfEdge();

    /// Add a half edge to half_edge_list().
    /// @param cell Cell containing the half edge.
    /// @param vfrom Half edge from vertex.
    /// @param vto Half edge to vertex
    /// @param hprev Previous half edge. Could be NULL.
    HALF_EDGE_TYPE * 
    _AddHalfEdge(CELL_TYPE * cell, 
                 VERTEX_TYPE * vfrom,
                 VERTEX_TYPE * vto,
                 HALF_EDGE_TYPE * hprev);

    /// Add a half edge to half_edge_list().
    /// @param icell Cell containing the half edge.
    /// @param vfrom Half edge from vertex.
    /// @param vto Half edge to vertex
    HALF_EDGE_TYPE *
    _AddHalfEdge(CELL_TYPE * cell, 
                 VERTEX_TYPE * vfrom,
                 VERTEX_TYPE * vto)
    { return _AddHalfEdge(cell, vfrom, vto, NULL); }

    /// Link half edge hprev to half edge hnext.
    /// @pre Both half edges are in the same polygon.
    void _LinkHalfEdgesInCell(HALF_EDGE_TYPE * hprev,
                              HALF_EDGE_TYPE * hnext);

    /// Add half_edgeB after half_edgeA to cyclic list
    ///   of half edges around edge.
    void _LinkHalfEdgesAroundEdge(HALF_EDGE_TYPE * half_edgeA,
                                  HALF_EDGE_TYPE * half_edgeB);

    /// Move boundary half edge to half_edge_from[0] for each vertex
    ///   in cell_vertex[].
    /// @param cell_vertex[] List of cell vertex indices.
    /// @pre Each vertex should already have been created.
    /// - Added: 11-24-2021 - RW
    void _MoveBoundaryHalfEdgeToHalfEdgeFrom0
    (const std::vector<int> & cell_vertex);

    /// Add cell to cell_list[].
    /// - Returns pointer to new cell.
    CELL_TYPE * _AddCell();

    /// Free list.
    template <typename ELEMENT_TYPE>
    void FreeList(std::vector<ELEMENT_TYPE> & v);

    /// Free all memory.
    void FreeAll();


  public:

    /// Destructor.
    ~HALF_EDGE_MESH_BASE() { FreeAll(); };


    // Get functions

    /// Return pointer to vertex with index iv.
    /// - Could return NULL if no vertex associated with iv,
    ///   e.g. vertex iv was deleted.
    const VERTEX_TYPE * Vertex(const int iv) const
    { return(vertex_list[iv]); }

    /// Return non-constant (NC) pointer to vertex with index iv.
    /// - Could return NULL if no vertex associated with iv,
    ///   e.g. vertex iv was deleted.
    /// - Added: 11-16-2021 - RW
    VERTEX_TYPE * VertexNC(const int iv)
    { return(vertex_list[iv]); }

    /// Return length (.size()) of vertex_list.
    /// - Note vertex_list could contain some NULL elements.
    int VertexListLength() const
    { return(vertex_list.size()); }

    /// Return pointer to half edge with index ihalf.
    /// - Could return NULL if no half edge associated with ihalf,
    ///   e.g. half edge ihalf was deleted.
    const HALF_EDGE_TYPE * HalfEdge(const int ihalf) const
    { return(half_edge_list[ihalf]); }

    /// Return non-constant (NC) pointer to half edge with index ihalf.
    /// - Could return NULL if no half edge associated with ihalf,
    ///   e.g. half edge ihalf was deleted.
    /// - Added: 11-16-2021 - RW
    HALF_EDGE_TYPE * HalfEdgeNC(const int ihalf)
    { return(half_edge_list[ihalf]); }

    /// Return length (.size()) of half_edge_list.
    /// - Note half_edge_list could contain some NULL elements.
    int HalfEdgeListLength() const
    { return(half_edge_list.size()); }

    /// Return pointer to cell with index icell.
    /// - Could return NULL if no cell associated with icell,
    ///   e.g. cell icell was deleted.
    const CELL_TYPE * Cell(const int icell) const
    { return(cell_list[icell]); }

    /// Return non-constant (NC) pointer to cell with index icell.
    /// - Could return NULL if no cell associated with icell,
    ///   e.g. cell icell was deleted.
    /// - Added: 11-16-2021 - RW
    CELL_TYPE * CellNC(const int icell)
    { return(cell_list[icell]); }

    /// Return length (.size()) of cell_list.
    /// - Note cell_list could contain some NULL elements.
    int CellListLength() const
    { return(cell_list.size()); }

    /// Return index of vertex v.
    int VertexIndex(const VERTEX_TYPE * v)
    { return(v->Index()); }

    /// Return index of half edge.
    int HalfEdgeIndex(const HALF_EDGE_TYPE * half_edge) const
    { return(half_edge->Index()); }

    /// Return index of cell.
    int CellIndex(const CELL_TYPE * cell) const
    { return(cell->Index()); }

    /// Return true if iv is the index of some vertex.
    bool IsVertexIndex(const int iv) const;

    /// Return true if ihalf is the index of some half edge.
    bool IsHalfEdgeIndex(const int ihalf) const;

    /// Return true if icell is the index of some half edge.
    bool IsCellIndex(const int icell) const;

    /// Count number of vertices.
    /// - Added: 11-24-2021 - RW
    int CountNumVertices() const;

    /// Count number of isolated vertices.
    /// - Isolated vertices are not in any mesh cell.
    /// - Added: 11-24-2021 - RW
    int CountNumIsolatedVertices() const;

    /// Count number of edges.
    /// - Added: 11-25-2021 - RW
    int CountNumEdges() const;

    /// Count number of boundary edges.
    /// - Added: 11-25-2021 - RW
    int CountNumBoundaryEdges() const;

    /// Count number of cells.
    /// - Added: 11-23-2021 - RW
    int CountNumCells() const;

    /// Count number of cells with a given number of vertices.
    /// - Added: 11-25-2021 - RW
    int CountNumCellsOfSize(const int numv) const;

    /// Count number of cells with number of vertices 
    ///  greater than or equal to.
    // - Added: 11-25-2021 - RW
    int CountNumCellsOfSizeGE(const int numv) const;

    /// Count number of triangles.
    /// - Added: 11-25-2021 - RW
    int CountNumTriangles() const
    { return CountNumCellsOfSize(3); }

    /// Count number of quadrilaterals.
    /// - Added: 11-25-2021 - RW
    int CountNumQuads() const
    { return CountNumCellsOfSize(4); }

    /// Count number of pentagons.
    /// - Added: 11-25-2021 - RW
    int CountNumPentagons() const
    { return CountNumCellsOfSize(5); }


    // Set functions.

    /// Add new vertex with index iv.
    /// - Returns a refence to the new vertex.
    /// @pre iv is not the index of any existing vertex.
    VERTEX_TYPE * AddVertex(const int iv);

    /// Add new vertex.
    /// - Return new vertex index.
    int AddVertex();

    /// Add numv vertices to vertex_list().
    /// - May only be called once.
    void AddVertices(const int numv);

    /// Add a cell.
    /// - Cell vertices are integers from 0 to VertexListLength().
    /// - Usually, max vertex index should be n-1, where n is the 
    ///   number of vertices.
    /// - Cells must have at least 3 vertices.
    /// - Returns index of new cell.
    int AddCell(const std::vector<int> & cell_vertex);

    /// Set coordinate.
    /// - Adds vertices (iv-vertex_list.size()) 
    ///     if vertex iv is not in vertex_list().
    template <typename CTYPE>
      void SetCoord(const int iv, const int ic, const CTYPE c);

    /// Set coordinate.
    /// - Adds vertices (iv-vertex_list.size()) 
    ///     if vertex iv is not in vertex_list().
    template <typename CTYPE>
      void SetCoord(const int iv, const CTYPE * c);


    // ***** Check routines ***

    /// Check data structure vertices.
    /// - Returns index of problem vertex.
    bool CheckVertices(int & iv, std::string & error_msg) const;

    /// Check data structure half edges.
    /// - Returns index of problem half edge.
    bool CheckHalfEdges(int & ihalf_edge, 
                        std::string & error_msg) const;

    /// Check data structure cells.
    /// - Returns index of problem cell.
    bool CheckCells(int & icell,
                    std::string & error_msg) const;

    /// Check vertices, half edges and cells.
    /// - Does not check manifold properties.
    bool CheckAll(std::string & error_msg) const;

    /// Check if mesh cells are consistenly oriented.
    /// - Return true if all adjacent cells are consistely oriented.
    /// - Returns index of half_edge where half_edge->Cell()
    ///   and half_edge->NextHalfEdgeAroundEdge()->Cell() have
    ///   opposite orientations.
    bool CheckOrientation(int & ihalf_edge) const;

    /// Check manifold edge property.
    /// - Return true if all edges have 2 or fewer incident cells.
    /// - Returns index of non-manifold edge.
    bool CheckManifoldEdges(int & ihalf_edge) const;

    /// Check manifold vertex property.
    /// - Return true if the cells on each vertex form a fan and
    ///   are consistently oriented around the vertex.
    /// - Returns index of non-manifold vertex.
    bool CheckManifoldVertices(int & iv) const;

    /// Check manifold edge and vertex properties.
    /// - Return true if all vertices and edges have manifold properties.
    /// @param flag_non_manifold_vertex 
    ///   True if cell incident on vertex iv do not form a fan.
    /// @param flag_non_manifold_edge
    ///   True if edge ihalf_edge has 3 or more incident cells.
    bool CheckManifold
    (int & iv, int & ihalf_edge, 
     bool & flag_non_manifold_vertex,
     bool & flag_non_manifold_edge) const;

    /// Check vertex index.
    /// - Return true if iv is an index of some mesh vertex.
    /// - Otherwise set error_msg and return false.
    /// - Added: 11-25-2021 - RW
    bool CheckVertexIndex(const int iv, std::string & error_msg) const;

  };


  /// Simple half edge mesh.
  typedef HALF_EDGE_MESH_BASE<VERTEX3D_A, HALF_EDGE_A, CELL_A>
  HALF_EDGE_MESH_A;


  // *****************************************************************
  // Member functions of VERTEX_BASE.
  // *****************************************************************

  // Initialize function.
  template <const int DIM, typename HALF_EDGE_PTR, typename CTYPE>
  void VERTEX_BASE<DIM,HALF_EDGE_PTR,CTYPE>::
  VERTEX_BASE::Init()
  {
    for (int ic = 0; ic < Dimension(); ic++)
      { coord[ic] = 0; }
  }


  // Return incident half edge whose from vertex is the current vertex
  //   and whose ToVertexIndex() is iv.
  // - Return NULL if no half edge found.
  template <const int DIM, typename HALF_EDGE_PTR, typename CTYPE>
  const HALF_EDGE_PTR
  VERTEX_BASE<DIM,HALF_EDGE_PTR,CTYPE>::
  FindIncidentHalfEdge(const int iv) const
  {
    for (int k = 0; k < NumHalfEdgesFrom(); k++) {
      const HALF_EDGE_PTR half_edge = KthHalfEdgeFrom(k);

      if (half_edge->ToVertexIndex() == iv)
        { return(half_edge); }
    }

    // No half edge found.
    return(NULL);
  }


  // Count number of half edges whose from vertex is the current vertex
  //   and whose to vertex is iv.
  template <const int DIM, typename HALF_EDGE_PTR, typename CTYPE>
  int VERTEX_BASE<DIM,HALF_EDGE_PTR,CTYPE>::
  CountNumIncidentHalfEdges(const int iv) const
  {
    int num = 0;
    for (int k = 0; k < NumHalfEdgesFrom(); k++) {
      const HALF_EDGE_PTR half_edge = KthHalfEdgeFrom(k);
      if (half_edge->ToVertexIndex() == iv)
        { num++; }
    }

    return(num);
  }


  // Move boundary half edge to half_edge_from[0].
  // - If there are no boundary half edges in half_edge_from[],
  //   but half_edge_from[k]->PreviousHalfEdgeInCell() is a boundary half edge,
  //   move half_edge_from[k] to half_edge_from[0].
  // - Revised: 11-24-2021 - RW
  template <const int DIM, typename HALF_EDGE_PTR, typename CTYPE>
  void VERTEX_BASE<DIM,HALF_EDGE_PTR,CTYPE>::
  _MoveBoundaryHalfEdgeToHalfEdgeFrom0()
  {
    if (NumHalfEdgesFrom() < 1) { return; }
    if (half_edge_from[0]->IsBoundary()) { 
      // half_edge_from[0] is already a boundary half edge.
      // Do nothing.
      return;
    }

    for (int k = 1; k < NumHalfEdgesFrom(); k++) {
      const HALF_EDGE_PTR half_edge = KthHalfEdgeFrom(k);
      if (half_edge->IsBoundary()) {
        std::swap(half_edge_from[0], half_edge_from[k]);
        return;
      }
    }

    // No boundary half edges found.

    // Extra processing in case cells are inconsistently oriented.
    // Check if half_edge_from[k]->PreviousHalfEdgeInCell()
    //   is a boundary half edge for some k.

    const HALF_EDGE_PTR prev_half_edge0 = 
      half_edge_from[0]->PrevHalfEdgeInCell();

    if (prev_half_edge0->IsBoundary()) { 
      // prev_half_edge0 is already a boundary half edge.
      // Do nothing.
      return;
    }

    for (int k = 1; k < NumHalfEdgesFrom(); k++) {
      const HALF_EDGE_PTR half_edge = KthHalfEdgeFrom(k);

      if (half_edge->PrevHalfEdgeInCell()->IsBoundary()) {
        std::swap(half_edge_from[0], half_edge_from[k]);
        return;
      }
    }

    return;
  }


  // Print coordinates.
  template <const int DIM, typename HALF_EDGE_PTR, typename CTYPE>
  template <typename OSTREAM_TYPE>
  void VERTEX_BASE<DIM,HALF_EDGE_PTR,CTYPE>::
  PrintCoord(OSTREAM_TYPE & out, const std::string & separator) const
  {
    if (Dimension() < 1) { return; }
    out << coord[0];

    for (int i = 1; i < Dimension(); i++) {
      out << separator;
      out << coord[i];
    }
  }

  // *****************************************************************
  // Member functions of HALF_EDGE_BASE.
  // *****************************************************************

  // Initialize function.
  template <typename VERTEX_PTR, typename HALF_EDGE_PTR, typename CELL_PTR>
  void HALF_EDGE_BASE<VERTEX_PTR,HALF_EDGE_PTR,CELL_PTR>::
  HALF_EDGE_BASE::Init()
  {
    next_half_edge_in_cell = NULL;
    prev_half_edge_in_cell = NULL;
    next_half_edge_around_edge = HALF_EDGE_PTR(this); 
    from_vertex = NULL;
    cell = NULL;
  }


  // Count number of half edges around edge.
  template <typename VERTEX_PTR, typename HALF_EDGE_PTR, typename CELL_PTR>
  int HALF_EDGE_BASE<VERTEX_PTR,HALF_EDGE_PTR,CELL_PTR>::
  CountNumHalfEdgesAroundEdge() const
  {
    // Cannot have more than max_num half edges around an edge.
    const int max_num = 
      FromVertex()->NumHalfEdgesFrom() +
      ToVertex()->NumHalfEdgesFrom();

    int num = 1;
    HALF_EDGE_PTR half_edge = this->NextHalfEdgeAroundEdge();
    // Check that num <= max_num to avoid infinite loop in case
    //   data structures is corrupted.
    while (half_edge != this &&
           num <= max_num) {
      half_edge = half_edge->NextHalfEdgeAroundEdge();
      num++;
    }

    return(num);
  }


  // Return pointer to half edge with minimum index
  //   in cycle of half edges around edge.
  // - Added: 11-25-2021 - RW
  template <typename VERTEX_PTR, typename HALF_EDGE_PTR, typename CELL_PTR>
  const HALF_EDGE_PTR 
  HALF_EDGE_BASE<VERTEX_PTR,HALF_EDGE_PTR,CELL_PTR>::
  MinIndexHalfEdgeAroundEdge() const
  {
    // Cannot have more than max_num half edges around an edge.
    const int max_num = 
      FromVertex()->NumHalfEdgesFrom() +
      ToVertex()->NumHalfEdgesFrom();

    HALF_EDGE_PTR half_edge = this->NextHalfEdgeAroundEdge();

    // Initialize.
    HALF_EDGE_PTR min_index_half_edge = half_edge;
    int min_index = half_edge->Index();

    // Check that num <= max_num to avoid infinite loop in case
    //   data structures is corrupted.
    int k = 0;
    half_edge = half_edge->NextHalfEdgeAroundEdge();
    do {
      if (half_edge->Index() < min_index) {
        min_index_half_edge = half_edge;
        min_index = half_edge->Index();
      }
      half_edge = half_edge->NextHalfEdgeAroundEdge();
      k++;
    } while (half_edge != this && k <= max_num);

    return min_index_half_edge;
  }


  // Return true if half_edgeB has same endpoints as this.
  template <typename VERTEX_PTR, typename HALF_EDGE_PTR, typename CELL_PTR>
  bool HALF_EDGE_BASE<VERTEX_PTR,HALF_EDGE_PTR,CELL_PTR>::
  SameEndpoints(const HALF_EDGE_PTR half_edgeB) const
  {
    if ((this->FromVertex() == half_edgeB->ToVertex()) &&
        (this->ToVertex() == half_edgeB->FromVertex()))
      { return(true); }

    if ((this->FromVertex() == half_edgeB->FromVertex()) &&
        (this->ToVertex() == half_edgeB->ToVertex()))
      { return(true); }

    return(false);
  }


  // Return pointer to previous edge around vertex iv.
  template <typename VERTEX_PTR, typename HALF_EDGE_PTR, typename CELL_PTR>
  const HALF_EDGE_PTR HALF_EDGE_BASE<VERTEX_PTR,HALF_EDGE_PTR,CELL_PTR>::
  PrevHalfEdgeAroundVertex(const int iv) const
  {
    if (FromVertexIndex() == iv) 
      { return(PrevHalfEdgeAroundFromVertex()); }
    else {
      // Return NextHalfEdgeInCell()->NextHalfEdgeAroundEdge() 
      //   so that repeated calls to PrevHalfEdgeAroundVertex()
      //   move in a consistent direction around vertex iv.
      return(NextHalfEdgeInCell()->NextHalfEdgeAroundEdge());
    }
  }


  // Return string of endpoints of half edge.
  // Separate coordinates with {separator}.
  // - Added: 11-25-2021 - RW
  template <typename VERTEX_PTR, typename HALF_EDGE_PTR, typename CELL_PTR>
  std::string HALF_EDGE_BASE<VERTEX_PTR,HALF_EDGE_PTR,CELL_PTR>::
  EndpointsStr(const std::string & separator) const
  {
    std::stringstream ss;

    if (FromVertex() == NULL) 
      { ss << "NULL" << separator; }
    else
      { ss << FromVertexIndex() << separator; }

    if (NextHalfEdgeInCell() == NULL)
      { ss << "NULL"; }
    else if (ToVertex() == NULL)
      { ss << "NULL"; }
    else
      { ss << ToVertexIndex(); }

    return ss.str();
  }


  // Return string of half edge index and endpoints.
  // Separate coordinates with {separator}.
  // - Added: 11-25-2021 - RW
  template <typename VERTEX_PTR, typename HALF_EDGE_PTR, typename CELL_PTR>
  std::string HALF_EDGE_BASE<VERTEX_PTR,HALF_EDGE_PTR,CELL_PTR>::
  IndexAndEndpointsStr(const std::string & separator) const
  {
    std::stringstream ss;

    ss << Index() << " (" << EndpointsStr(separator) << ")";

    return ss.str();
  }


  // Print endpoints of half edge.
  template <typename VERTEX_PTR, typename HALF_EDGE_PTR, typename CELL_PTR>
  template <typename OSTREAM_TYPE>
  void HALF_EDGE_BASE<VERTEX_PTR,HALF_EDGE_PTR,CELL_PTR>::
  PrintEndpoints(OSTREAM_TYPE & out, 
                 const std::string & separator) const
  {
    out << FromVertexIndex() << separator;

    // Check that next_half_edge_in_cell is not NULL.
    if (NextHalfEdgeInCell() == NULL) {
      // Data structure problem. Cannot determine ToVertex(). Print a "*".
      out << "*";
    }
    else {
      // Print ToVertex().
      out << ToVertexIndex();
    }
  }


  // Print index and endpoints of half edge.
  template <typename VERTEX_PTR, typename HALF_EDGE_PTR, typename CELL_PTR>
  template <typename OSTREAM_TYPE>
  void HALF_EDGE_BASE<VERTEX_PTR,HALF_EDGE_PTR,CELL_PTR>::
  PrintIndexAndEndpoints(OSTREAM_TYPE & out, 
                         const std::string & separator) const
  {
    out << Index() << " (";
    PrintEndpoints(out, separator);
    out << ")";
  }


  // *****************************************************************
  // Member functions of HALF_EDGE_MESH_BASE.
  // *****************************************************************

  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  VERTEX_TYPE *  
  HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  _CreateVertex(const int iv)
  {
    if (iv < 0) { 
      throw SIMPLE_EXCEPTION
        ("Illegal argument to HALF_EDGE_MESH_BASE::_CreateVertex(). Vertex index iv must be non-negative.");
    }

    if (iv >= VertexListLength()) {
      vertex_list.resize(iv+1, NULL); 
    }

    // Check that resize worked.
    if (iv >= VertexListLength()) {
      throw SIMPLE_EXCEPTION
        ("Resize of vertex_list failed in HALF_EDGE_MESH_BASE::_CreateVertex(). Probably out of memory."); 
    }

    if (Vertex(iv) == NULL) {
      vertex_list[iv] = new VERTEX_TYPE;
      vertex_list[iv]->index = iv;
    }

    return(vertex_list[iv]);
  }


  // Return true if iv is the index of some vertex.
  // - Added: 11-25-2021 - RW
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  bool HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  IsVertexIndex(const int iv) const
  {
    if (iv < 0) { return false; }
    if (iv >= VertexListLength()) { return false; }
    if (Vertex(iv) == NULL) { return false; }

    return(true);
  }


  // Return true if ihalf_edge is the index of some half edge.
  // - Added: 11-23-2021 - RW
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  bool HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  IsHalfEdgeIndex(const int ihalf_edge) const
  {
    if (ihalf_edge < 0) { return false; }
    if (ihalf_edge >= HalfEdgeListLength()) { return false; }
    if (HalfEdge(ihalf_edge) == NULL) { return false; }

    return(true);
  }


  // Return true if icell is the index of some cell.
  // - Added: 11-23-2021 - RW
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  bool HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  IsCellIndex(const int icell) const
  {
    if (icell < 0) { return false; }
    if (icell >= CellListLength()) { return false; }
    if (Cell(icell) == NULL) { return false; }

    return(true);
  }


  // Count number of vertices.
  // - Added: 11-24-2021 - RW
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  int HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CountNumVertices() const
  {
    int num_vertices = 0;
    for (int iv = 0; iv < VertexListLength(); iv++) {
      const VERTEX_TYPE * v = Vertex(iv);
      if (v != NULL)
        { num_vertices++; }
    }

    return num_vertices;
  }


  // Count number of isolated vertices.
  // - Added: 11-24-2021 - RW
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  int HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CountNumIsolatedVertices() const
  {
    int num_isolated_vertices = 0;
    for (int iv = 0; iv < VertexListLength(); iv++) {
      const VERTEX_TYPE * v = Vertex(iv);
      if (v != NULL) {
        if (v->NumHalfEdgesFrom() == 0) 
          { num_isolated_vertices++; }
      }
    }

    return num_isolated_vertices;
  }


  // Count number of edges.
  // - Added: 11-25-2021 - RW
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  int HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CountNumEdges() const
  {
    int num_edges = 0;
    for (int ihalf_edge = 0; ihalf_edge < HalfEdgeListLength(); 
         ihalf_edge++) {
      const HALF_EDGE_TYPE * half_edge = HalfEdge(ihalf_edge);
      if (half_edge != NULL) {
        const HALF_EDGE_TYPE * min_index_half_edge =
          half_edge->MinIndexHalfEdgeAroundEdge();
        if (half_edge == min_index_half_edge) 
          { num_edges++; }
      }
    }

    return num_edges;
  }

  // Count number of boundary edges
  // - Added: 11-25-2021 - RW
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  int HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CountNumBoundaryEdges() const
  {
    int num_boundary_edges = 0;
    for (int ihalf_edge = 0; ihalf_edge < HalfEdgeListLength(); 
         ihalf_edge++) {
      const HALF_EDGE_TYPE * half_edge = HalfEdge(ihalf_edge);
      if (half_edge != NULL) {
        if (half_edge->IsBoundary())
          { num_boundary_edges++; }
      }
    }

    return num_boundary_edges;
  }


  // Count number of cells.
  // - Added: 11-23-2021 - RW
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  int HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CountNumCells() const
  {
    int num_cells = 0;
    for (int icell = 0; icell < CellListLength(); icell++) {
      const CELL_TYPE * cell = Cell(icell);
      if (cell != NULL)
        { num_cells++; }
    }

    return num_cells;
  }


  // Count number of cells with a given number of vertices.
  // - Added: 11-25-2021 - RW
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  int HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CountNumCellsOfSize(const int numv) const
  {
    int num_cells = 0;
    for (int icell = 0; icell < CellListLength(); icell++) {
      const CELL_TYPE * cell = Cell(icell);
      if (cell != NULL) {
        if (cell->NumVertices() == numv)
          { num_cells++; }
      }
    }

    return num_cells;
  }


  // Count number of cells with number of vertices 
  //  greater than or equal to.
  // - Added: 11-25-2021 - RW
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  int HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CountNumCellsOfSizeGE(const int numv) const
  {
    int num_cells = 0;
    for (int icell = 0; icell < CellListLength(); icell++) {
      const CELL_TYPE * cell = Cell(icell);
      if (cell != NULL) {
        if (cell->NumVertices() >= numv)
          { num_cells++; }
      }
    }

    return num_cells;
  }


  // Add new vertex with index iv.
  // - Returns a refence to the new vertex.
  // @pre iv is not the index of any existing vertex.
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  VERTEX_TYPE *
  HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  AddVertex(const int iv)
  {
    if (iv < vertex_list.size()) {
      if (vertex_list[iv] != NULL) {
        throw SIMPLE_EXCEPTION
          ("Illegal vertex identifier in call to HALF_EDGE_MESH_BASE::AddVertex(). Mesh already has vertex with given identifier.");
      }
    }
    else {
      vertex_list.resize(iv+1, NULL);
    }

    return(_CreateVertex(iv));
  }


  // Add new vertex.
  // - Return new vertex index.
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  int HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  AddVertex()
  {
    // CORRECTION: 11-16-2021 - RW
    // INCORRECT: const int iv = vertex_list.size() + 1;
    const int iv = vertex_list.size();
    AddVertex(iv);
    return(iv);
  }


  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  void HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  AddVertices(const int numv)
  {
    if (numv < 1) { return; }

    // Check that vertex table is empty.
    if (vertex_list.size() != 0) {
      throw SIMPLE_EXCEPTION
        ("Illegal call to HALF_EDGE_MESH_BASE::AddVertices(). Mesh already has vertices.");
    }

    vertex_list.resize(numv, NULL);

    // Check that resize worked.
    if (numv != VertexListLength()) {
      // Resize failed. Probably out of memory.
      throw SIMPLE_EXCEPTION
        ("Resize of vertex_list failed in HALF_EDGE_MESH_BASE::AddVertices(). Probably out of memory.");
    }

    for (int i = 0; i < vertex_list.size(); i++) {
      vertex_list[i] = new VERTEX_TYPE;
      vertex_list[i]->index = i;
    }
  }


  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  HALF_EDGE_TYPE *  
  HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  _AddHalfEdge()
  {
    const int ihalf_edge  = half_edge_list.size();
    half_edge_list.resize(ihalf_edge+1);

    // Check that resize worked.
    if (ihalf_edge >= HalfEdgeListLength()) {
      // Resize failed. Probably out of memory.
      throw SIMPLE_EXCEPTION
        ("Resize of half_edge_list failed in HALF_EDGE_MESH_BASE::_AddHalfEdge(). Probably out of memory.");
    }

    half_edge_list[ihalf_edge] = new HALF_EDGE_TYPE;
    HALF_EDGE_TYPE * half_edge_ptr = half_edge_list[ihalf_edge];
    half_edge_ptr->index = ihalf_edge;

    return(half_edge_ptr);
  }


  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  HALF_EDGE_TYPE * 
  HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  _AddHalfEdge(CELL_TYPE * cell, 
               VERTEX_TYPE * vfrom,
               VERTEX_TYPE * vto,
               HALF_EDGE_TYPE * hprev)
  {
    if (cell == NULL || vfrom == NULL || vto == NULL) {
      throw SIMPLE_EXCEPTION
        ("Illegal argument to HALF_EDGE_MESH_BASE::_AddHalfEdge(). Arguments cell, vfrom, and vto cannot be NULL."); 
    }

    HALF_EDGE_TYPE * half_edge = _AddHalfEdge();

    half_edge->cell = cell;
    cell->num_vertices++;
    half_edge->from_vertex = vfrom;
    half_edge->prev_half_edge_in_cell = hprev;
    if (hprev != NULL) 
      { hprev->next_half_edge_in_cell = half_edge; }

    // Link half_edge with other half edges around edge.
    HALF_EDGE_TYPE * half_edgeB =
      vto->FindIncidentHalfEdge(vfrom->Index());
    if (half_edgeB == NULL) {
      // Check whether there is a half edge around edge
      //   with the same orientation as half_edge.
      HALF_EDGE_TYPE * half_edgeC =
        vfrom->FindIncidentHalfEdge(vto->Index());
      if (half_edgeC == NULL)
        { half_edge->next_half_edge_around_edge = half_edge; }
      else {
        // Link half_edge with half_edgeC, even though they
        //   have the same orientation.
        _LinkHalfEdgesAroundEdge(half_edgeC, half_edge);
      }
    }
    else {
      _LinkHalfEdgesAroundEdge(half_edgeB, half_edge);
    }

    vfrom->half_edge_from.push_back(half_edge);

    // Deleted: 11-24-2021 - RW
    // Revised version of _MoveBoundaryHalfEdgeToHalfEdgFrom(),
    //   requires prev_half_edge_in_cell to be set before call.
    // Calls have been moved to _AddCell().
    // Making these calls here will cause a segmentation fault.
    // vfrom->_MoveBoundaryHalfEdgeToHalfEdgeFrom0();
    // vto->_MoveBoundaryHalfEdgeToHalfEdgeFrom0();

    return(half_edge);
  }


  // Link half edge hprev to half edge hnext.
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  void HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  _LinkHalfEdgesInCell(HALF_EDGE_TYPE * hprev,
                       HALF_EDGE_TYPE * hnext)
  {
    if (hprev == NULL || hnext == NULL) {
      throw SIMPLE_EXCEPTION
        ("Illegal argument to HALF_EDGE_MESH_BASE::_LinkHalfEdgesInCell(). Arguments hprev and hnext cannot be NULL.");
    }

    if (hprev->Cell() != hnext->Cell()) {
      throw SIMPLE_EXCEPTION
        ("Link of half edges failed in HALF_EDGE_MESH_BASE::_LinkHalfEdgesInCell(). hprev and hnext are in different cells.");
    }

    if ((hprev->NextHalfEdgeInCell() != NULL) && 
        (hprev->NextHalfEdgeInCell() != hnext)) {
      throw SIMPLE_EXCEPTION
        ("Link of half edges failed in HALF_EDGE_MESH_BASE::_LinkHalfEdgesInCell(). hprev->next_half_edge_in_cell is already defined.");
    }

    if ((hnext->PrevHalfEdgeInCell() != NULL) &&
        (hnext->PrevHalfEdgeInCell() != hprev)) {
      throw SIMPLE_EXCEPTION
        ("Link of half edges failed in HALF_EDGE_MESH_BASE::_LinkHalfEdgesInCell(). hnext->prev_half_edge_in_cell is already defined.");
    }

    hprev->next_half_edge_in_cell = hnext;
    hnext->prev_half_edge_in_cell = hprev;
  }


  // Add half_edgeB after half_edgeA to cyclic list
  //   of half edges around edge.
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  void HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  _LinkHalfEdgesAroundEdge(HALF_EDGE_TYPE * half_edgeA,
                           HALF_EDGE_TYPE * half_edgeB)
  {
    half_edgeB->next_half_edge_around_edge = half_edgeB;
    std::swap(half_edgeA->next_half_edge_around_edge,
              half_edgeB->next_half_edge_around_edge);
  }

  // Move boundary half edge to half_edge_from[0] for each vertex
  //   in cell_vertex[].
  // - Added: 11-24-2021 - RW
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  void HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  _MoveBoundaryHalfEdgeToHalfEdgeFrom0
  (const std::vector<int> & cell_vertex)
  {
    for (int i = 0; i < cell_vertex.size(); i++) {
      const int iv = cell_vertex[i];
      VERTEX_TYPE * v = vertex_list[iv];
      if (v == NULL) {
        throw SIMPLE_EXCEPTION
          ("Programming error. Attempt to access non-existant vertex in _MoveBoundaryHalfEdgeToHalfEdgeFrom0().");
      }

      v->_MoveBoundaryHalfEdgeToHalfEdgeFrom0();
    }
  }




  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  CELL_TYPE *  
  HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  _AddCell()
  {
    const int icell = cell_list.size();
    cell_list.resize(icell+1);

    // Check that resize worked.
    if (icell >= CellListLength()) {
      // Resize failed. Probably out of memory.
      throw SIMPLE_EXCEPTION
        ("Resize of cell_list failed in HALF_EDGE_MESH_BASE::_AddCell(). Probably out of memory.");
    }

    cell_list[icell] = new CELL_TYPE;
    CELL_TYPE * cell_ptr = cell_list[icell];
    cell_ptr->index = icell;

    return(cell_ptr);
  }


  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  int HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  AddCell(const std::vector<int> & cell_vertex)
  {
    if (cell_vertex.size() < 3) {
      throw SIMPLE_EXCEPTION
        ("Illegal argument to HALF_EDGE_MESH_BASE::AddCell(). Vector cell_vertex[] must have 3 or more vertices.");
    }

    CELL_TYPE * cell = _AddCell();

    // Create first half edge;
    const int iv0 = cell_vertex[0];
    const int iv1 = cell_vertex[1];
    VERTEX_TYPE * v0 = _CreateVertex(iv0);
    VERTEX_TYPE * v1 = _CreateVertex(iv1);

    HALF_EDGE_TYPE * half_edge0 = _AddHalfEdge(cell, v0, v1);
    cell->half_edge = half_edge0;

    HALF_EDGE_TYPE * hprev = half_edge0;

    for (int i0 = 1; i0 < cell_vertex.size(); i0++) {
      const int i1 = (i0+1) % cell_vertex.size();
      const int iv0 = cell_vertex[i0];
      const int iv1 = cell_vertex[i1];
      VERTEX_TYPE * v0 = vertex_list[iv0];
      VERTEX_TYPE * v1 = _CreateVertex(iv1);
      
      HALF_EDGE_TYPE * half_edge = _AddHalfEdge(cell, v0, v1, hprev);

      hprev = half_edge;
    }

    // Link last half edge (hprev) and first half edge (half_edge0)
    _LinkHalfEdgesInCell(hprev, half_edge0);

    // - Added: 11-24-2021 - RW
    // - This call must be AFTER _LinkHalfEdgesInCell(hprev, half_edge0).
    _MoveBoundaryHalfEdgeToHalfEdgeFrom0(cell_vertex);

    if (cell_vertex.size() != cell->NumVertices()) {
      throw SIMPLE_EXCEPTION
        ("Error in HALF_EDGE_MESH_BASE::AddCell(). Incorrect number of vertices in cell.");

    }

    return(CellIndex(cell));
  }


  // Set coordinate.
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  template <typename CTYPE>
  void HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  SetCoord(const int iv, const int ic, const CTYPE c)
  {
    VERTEX_TYPE * v = _CreateVertex(iv);

    if (ic < 0 || ic >= v->Dimension()) {
      throw SIMPLE_EXCEPTION
        ("Illegal argument to HALF_EDGE_MESH_BASE::SetCoord(). Coordinate ic out of bounds.");
    }

    v->coord[ic] = c;
  }


  // Set coordinate.
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  template <typename CTYPE>
  void HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  SetCoord(const int iv, const CTYPE * c)
  {
    VERTEX_TYPE * v = _CreateVertex(iv);

    for (int ic = 0; ic < v->Dimension(); ic++) 
      { v->coord[ic] = c[ic]; }
  }


  // Free list.
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  template <typename ELEMENT_TYPE>
  void HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  FreeList(std::vector<ELEMENT_TYPE> & list)
  {
    for (int i = 0; i < list.size(); i++) {
      if (list[i] != NULL) {
        delete list[i]; 
        list[i] = NULL;
      }
    }
  }


  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  void HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  FreeAll()
  {
    FreeList(vertex_list);
    FreeList(half_edge_list);
    FreeList(cell_list);
  }


  // *****************************************************************
  // Member check functions of HALF_EDGE_MESH_BASE.
  // *****************************************************************

  // Check vertices.
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  bool HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CheckVertices(int & iv,
                 std::string & error_msg) const
  {
    std::stringstream ss;

    // Initialize.
    iv = 0;

    for (int jv = 0; jv < vertex_list.size(); jv++) {
      iv = jv;

      const VERTEX_TYPE * v = Vertex(jv);
      if (v == NULL) { continue; }

      if (v->Index() != jv) {
        ss << "Incorrect vertex index for vertex " << jv << ".";
        error_msg = ss.str();
        return(false);
      }

      bool flag_boundary = false;
      const HALF_EDGE_TYPE * boundary_half_edge = NULL;
      for (int k = 0; k < v->NumHalfEdgesFrom(); k++) {
        const HALF_EDGE_TYPE * half_edge = v->KthHalfEdgeFrom(k);

        if (half_edge == NULL) {
          ss << "Vertex " << jv << " list half_edge_from[] contains a NULL pointer.";
          error_msg = ss.str();
          return(false);
        }

        if (half_edge->FromVertex() != v) {
          ss << "Error in list half_edge_from[] for vertex " << jv 
             << ".  List incorrectly includes half edge " << half_edge->Index()
             << ".";
          error_msg = ss.str();
          return(false);
        }

        if (half_edge->IsBoundary()) {
          // Vertex is on the boundary.
          flag_boundary = true;
          boundary_half_edge = half_edge;
        }
      }

      if (flag_boundary) {
        // Must be at least one half edge for flag_boundary to be true.
        const HALF_EDGE_TYPE * half_edge =
          v->KthHalfEdgeFrom(0);
        if (!half_edge->IsBoundary()) {
          // Error. First half edge should be a boundary half edge.
          ss << "Vertex " << jv << " is on boundary half edge ";
          boundary_half_edge->PrintIndexAndEndpoints(ss, ",");
          ss << " but first incident half edge ";
          half_edge->PrintIndexAndEndpoints(ss, ",");
          ss << " is not a boundary half edge.";
          error_msg = ss.str();
          return(false);
        }
      }
    }

    return(true);
  }


  // Check half edges.
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  bool HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CheckHalfEdges(int & ihalf_edge,
                 std::string & error_msg) const
  {
    std::stringstream ss;

    for (int j = 0; j < HalfEdgeListLength(); j++) {

      ihalf_edge = j;

      const HALF_EDGE_TYPE * half_edge = half_edge_list[j];
      if (half_edge == NULL) { continue; }
      if (half_edge->Index() != j) {
        ss << "Incorrect half edge index for half edge " << j << ".";
        error_msg = ss.str();
        return(false);
      }

      const VERTEX_TYPE * v = half_edge->FromVertex();
      if (v == NULL) {
        ss << "Missing (NULL) from vertex in half edge " << j << ".";
        error_msg = ss.str();
        return(false);
      }

      const int num_match =
        std::count(v->half_edge_from.begin(), 
                   v->half_edge_from.end(),
                   half_edge);
      if (num_match < 1) {
        // Half edge not in incident half edge list of vertex v.
        ss << "Half edge does not appear in half_edge_from[] list for vertex " << v->Index() << ".";
        error_msg = ss.str();
        return(false); 
      }
      else if (num_match > 1) {
        // Half edge appears more than once in incident half edge
        //   list of vertex v.
        ss << "Half edge appears more than once in half_edge_from[] list for vertex " << v->Index() << ".";
        error_msg = ss.str();
        return(false);
      }

      if (half_edge->IsBoundary() &&
          !(v->KthHalfEdgeFrom(0))->IsBoundary()) {
        ss << "Half edge " << half_edge->Index()
           << " is a boundary half edge but v->KthHalfEdgeFrom(0) is not a boundary half edge.";
        error_msg = ss.str();
        return(false);
      }

      const CELL_TYPE * cell = half_edge->Cell();
      const HALF_EDGE_TYPE * next_half_edge =
        half_edge->NextHalfEdgeInCell();
      const HALF_EDGE_TYPE * prev_half_edge =
        half_edge->PrevHalfEdgeInCell();

      if (next_half_edge == NULL) { return(false); }
      if (prev_half_edge == NULL) { return(false); }

      if (next_half_edge->PrevHalfEdgeInCell() != half_edge) {
        ss << "Error. half_edge->NextHalfEdgeInCell()->PrevHalfEdgeInCell() is not half_edge.\n";
        ss << "  Half edge: " << half_edge->Index()
           << "  Next half edge: " << next_half_edge->Index()
           << "\n";
        error_msg = ss.str();
        return(false); 
      }

      if (prev_half_edge->NextHalfEdgeInCell() != half_edge) 
        { return(false); }

      if (next_half_edge->Cell() != cell) {
        // next_half_edge should be in same cell as half_edge.
        ss << "Error.  Consecutive half edges, "
           << half_edge->Index() << " and "
           << next_half_edge->Index() << ", are in different cells.";
        error_msg = ss.str();
        return(false); 
      }
          
      // Check half edges around edge.
      HALF_EDGE_TYPE * half_edgeX = 
        half_edge->next_half_edge_around_edge;

      if (half_edgeX != half_edge) {

        if (half_edgeX == NULL) { 
          ss << "Error for half edge " << half_edge->Index()
             << ". half_edge->next_half_edge_around_edge == NULL.";
          error_msg = ss.str();
          return(false); 
        }

        if (!half_edge->SameEndpoints(half_edgeX)) {
          ss << "Error. Two half edges around edge have different endpoints. Half edges: " << half_edge->Index() << " " << half_edgeX->Index() << ".";
          error_msg = ss.str();
          return(false); 
        }

        if (half_edge->Cell() == half_edgeX->Cell()) {
          // Two half edges around the same edge cannot be in the same cell.
          ss << "Error.  Two half edges around edge are in the same cell. Half edges: " << half_edge->Index() << " " << half_edgeX->Index() << ".";
          error_msg = ss.str();
          return(false); 
        }
      }
    }

    // Check that number of half edges around edge match
    //   number of corresponding half edges in vectors
    //   ivfrom->half_edge_from and ivto->half_edge_from.

    // Vector keeps from revisiting half edges.
    std::vector<bool> is_visited(HalfEdgeListLength(), false);

    for (int j = 0; j < HalfEdgeListLength(); j++) {
      ihalf_edge = j;

      const HALF_EDGE_TYPE * half_edge = half_edge_list[j];
      if (half_edge == NULL) { continue; }

      const int ihalf_edge = half_edge->Index();
      if (is_visited[ihalf_edge]) { 
        // Skip half_edge.  Already visited.
        continue;
      }

      const int numh = half_edge->CountNumHalfEdgesAroundEdge();

      const VERTEX_TYPE * vto = half_edge->ToVertex();
      const VERTEX_TYPE * vfrom = half_edge->FromVertex();
      const int ivto = vto->Index();
      const int ivfrom = vfrom->Index();
      const int numh2 = 
        vto->CountNumIncidentHalfEdges(ivfrom) +
        vfrom->CountNumIncidentHalfEdges(ivto);

      if (numh != numh2) {
        ss << "Inconsistency between half edges around edge and vertex incident lists for edge (" << ivfrom << "," << ivto << ").";
        error_msg = ss.str();
        return(false); 
      }

      // Mark all visited half edges so that they are not processed again.
      // Reduces time spent checking half edges around edge.
      const HALF_EDGE_TYPE * half_edgeX = half_edge;
      for (int k = 0; k < numh; k++) {
        const int ihalf_edgeX = half_edgeX->Index();
        is_visited[ihalf_edgeX] = true;
        half_edgeX = half_edgeX->NextHalfEdgeAroundEdge();
      }
    }

    return(true);
  }


  // Check cells.
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  bool HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CheckCells(int & icell, std::string & error_msg) const
  {
    std::stringstream ss;

    for (int j = 0; j < CellListLength(); j++) {

      icell = j;

      const CELL_TYPE * cell = cell_list[j];
      if (cell == NULL) { continue; }
      if (cell->Index() != j) {
        ss << "Incorrect cell index for cell " << j << ".";
        error_msg = ss.str();
        return(false);
      }

      const HALF_EDGE_TYPE * half_edge0 = cell->HalfEdge();
      if (half_edge0->Cell() != cell) {
        ss << "Incorrect half edge stored in cell " << icell << ".";
        error_msg = ss.str();
        return(false);
      }

      int k = 1;
      const HALF_EDGE_TYPE * half_edge = half_edge0;

      while (k < cell->NumVertices()) {
        half_edge = half_edge->NextHalfEdgeInCell();

        if (half_edge == half_edge0) {
          // Fewer than expected number of vertices in cell.
          ss << "Incorrect number of vertices (" << cell->NumVertices()
             << ") stored in cell " << icell << "."
             << " Counted " << k << " vertices.";
          error_msg = ss.str();
          return(false);
        }

        k = k+1;
      }

      if (half_edge->NextHalfEdgeInCell() != half_edge0) {
        // More than expected number of vertices in cell.
        ss << "Incorrect number of vertices (" << cell->NumVertices()
           << ") stored in cell " << icell << "."
           << " Cell has more than " << cell->NumVertices() 
           << " vertices.";
        error_msg = ss.str();
        return(false);
      }
    }

    return(true);
  }

  // Check vertices, half edges and cells.
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  bool HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CheckAll(std::string & error_msg) const
  {
    int iv, ihalf_edge, icell;
    std::stringstream ss;

    if (!CheckVertices(iv, error_msg)) { 
      if (error_msg == "") {
        ss << "Error related to vertex " << iv << ".";
        error_msg = ss.str();
      }
      return(false);
    }

    if (!CheckHalfEdges(ihalf_edge, error_msg)) { 
      if (error_msg == "") {
        ss << "Error related to half edge " << ihalf_edge << ".";
        error_msg = ss.str();
      }
      return(false); 
    }

    if (!CheckCells(icell, error_msg)) { 
      if (error_msg == "") {
        ss << "Error related to cell " << icell << ".";
        error_msg = ss.str();
      }
      return(false); 
    }

    return(true);
  }


  // Check if mesh cells are consistenly oriented.
  // - Return true if all adjacent cells are consistely oriented.
  // - Returns index of half_edge where half_edge->Cell()
  //   and half_edge->NextHalfEdgeAroundEdge()->Cell() have
  //   opposite orientations.
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  bool HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CheckOrientation(int & ihalf_edge) const
  {
    for (int i = 0; i < HalfEdgeListLength(); i++) {
      ihalf_edge = i;

      const HALF_EDGE_TYPE * half_edge = HalfEdge(i);
      if (half_edge == NULL) { continue; }

      const HALF_EDGE_TYPE * half_edgeB = 
        half_edge->NextHalfEdgeAroundEdge();

      if (half_edge == half_edgeB) {
        // Skip boundary half edge.
        continue;
      }
      
      if (half_edge->FromVertex() == half_edgeB->FromVertex()) {
        // Cells half_edge->Cell() and half_edgeB->Cell()
        //   have inconsistent orientations.
        return(false);
      }
    }

    return(true);
  }


  // Check manifold edge property.
  // - Returns index of non-manifold edge.
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  bool HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CheckManifoldEdges(int & ihalf_edge) const
  {
    for (int i = 0; i < HalfEdgeListLength(); i++) {
      ihalf_edge = i;

      const HALF_EDGE_TYPE * half_edge = HalfEdge(i);
      if (half_edge == NULL) { continue; }

      const int numh = half_edge->CountNumHalfEdgesAroundEdge();
      if (numh >= 3) { return(false); }
    }

    return(true);
  }


  // Check manifold vertex property.
  // - Return true if the cells on each vertex form a fan and
  //   are consistently oriented around the vertex.
  // - Returns index of non-manifold vertex.
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  bool HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CheckManifoldVertices(int & iv) const
  {
    for (int j = 0; j < VertexListLength(); j++) {
      iv = j;

      const VERTEX_TYPE * v = Vertex(j);
      if (v == NULL) { continue; }

      const int numh = v->NumHalfEdgesFrom();

      if (numh == 0) { continue; }

      const HALF_EDGE_TYPE * half_edge0 =
        v->half_edge_from[0];
      const HALF_EDGE_TYPE * half_edge =
        half_edge0->PrevHalfEdgeAroundVertex(iv);

      int num_cells = 1;
      while (half_edge != half_edge0 &&
             !(half_edge->IsBoundary()) &&
             num_cells < numh) {
        num_cells++;
        half_edge = half_edge->PrevHalfEdgeAroundVertex(iv);
      }

      if (num_cells != numh) { return(false); }
    }

    return(true);
  }


  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  bool HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CheckManifold
  (int & iv, int & ihalf_edge, 
   bool & flag_non_manifold_vertex,
   bool & flag_non_manifold_edge) const
  {
    flag_non_manifold_vertex = !CheckManifoldVertices(iv);
    flag_non_manifold_edge = !CheckManifoldEdges(ihalf_edge);

    return(!(flag_non_manifold_vertex || flag_non_manifold_edge));
  }


  // Check vertex index.
  // - Return true if iv is an index of some mesh vertex.
  // - Added: 11-25-2021 - RW
  template <typename VERTEX_TYPE, typename HALF_EDGE_TYPE,
            typename CELL_TYPE>
  bool HALF_EDGE_MESH_BASE<VERTEX_TYPE, HALF_EDGE_TYPE, CELL_TYPE>::
  CheckVertexIndex(const int iv, std::string & error_msg) const
  {
    std::stringstream ss;

    if (this->VertexListLength() <= 0) {
      ss << "Mesh has no vertices.\n";
      error_msg = ss.str();
      return false;
    }

    if (iv < 0) {
      ss << "Illegal negative vertex index: " << iv << "\n";
      ss << "  Vertex indices cannot be negative.\n";
      error_msg = ss.str();
      return false;
    }

    if (iv >= this->VertexListLength()) {
      ss << "Illegal vertex index: " << iv << "\n";
      ss << "  Maximum vertex index: " 
         << this->VertexListLength()-1 << "\n";
      error_msg = ss.str();
      return false;
    }

    if (Vertex(iv) == NULL) {
      ss << "No vertex has index: " << iv << "\n";
      error_msg = ss.str();
      return false;
    }

    return true;
  }

}

#endif
