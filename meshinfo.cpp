/// \file meshinfo.cpp
/// Print mesh information.

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "half_edge_mesh.hpp"
#include "half_edge_mesh_DCMT.hpp"
#include "half_edge_mesh_IO.hpp"

using namespace HMESH;
using namespace std;

// Global variables
const char * input_filename(NULL);
bool flag_more_info(false);

// Functions
void parse_command_line(int argc, char ** argv);
void usage_error();
void help_msg();

// Print functions.
template <typename MESH_TYPE>
void print_mesh_size(const MESH_TYPE & mesh, bool flag_more_info);
template <typename MESH_TYPE>
void print_min_max_edge_lengths
(const MESH_TYPE & mesh, const bool flag_more_info);
template <typename MESH_TYPE>
void print_min_cell_edge_length_ratio
(const MESH_TYPE & mesh, bool flag_more_info);
template <typename MESH_TYPE>
void print_manifold_info
(const MESH_TYPE & mesh, bool flag_more_info);
template <typename MESH_TYPE>
void print_min_max_angles
(const MESH_TYPE & mesh, bool flag_more_info);

// Define M_PI, if necessary.
#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif


int main(int argc, char ** argv)
{
  HALF_EDGE_MESH_DCMT_A mesh;

  parse_command_line(argc, argv);

  open_and_read_off_file(input_filename, mesh);

  try {

    print_mesh_size(mesh, flag_more_info);
    print_min_max_edge_lengths(mesh, flag_more_info);
    print_min_cell_edge_length_ratio(mesh, flag_more_info);
    print_min_max_angles(mesh, flag_more_info);
    print_manifold_info(mesh, flag_more_info);
  }
  catch (SIMPLE_EXCEPTION & e) {
    cerr << e.what() << endl;
    cerr << "Exiting." << endl;
    exit(-1);
  }

}


// *** PRINT MESH INFORMATION ***

template <typename MESH_TYPE>
void print_mesh_size(const MESH_TYPE & mesh, const bool flag_more_info)
{
  const int num_vertices = mesh.CountNumVertices();
  const int num_isolated_vertices = 
    mesh.CountNumIsolatedVertices();
  const int num_edges = mesh.CountNumEdges();
  const int num_boundary_edges = mesh.CountNumBoundaryEdges();
  const int num_cells = mesh.CountNumCells();

  cout << "Number of mesh vertices: " 
       << num_vertices - num_isolated_vertices << endl;
  if (num_isolated_vertices > 0) {
    cout << "Total number of vertices in input file: " 
         << num_vertices << endl;
  }
  cout << "Number of mesh edges: " << num_edges << endl;
  cout << "Number of boundary mesh edges: " 
       << num_boundary_edges << endl;

  cout << "Number of mesh cells: " << num_cells << endl;

  if (flag_more_info) {
    const int num_triangles = mesh.CountNumTriangles();
    cout << "  Number of mesh triangles: " << num_triangles << endl;
    const int num_quads = mesh.CountNumQuads();
    cout << "  Number of mesh quadrilaterals: " << num_quads << endl;
    const int num_pentagons = mesh.CountNumPentagons();
    const int num_large_cells = mesh.CountNumCellsOfSizeGE(6);
    if (num_pentagons > 0) {
      cout << "  Number of mesh pentagons: " << num_pentagons << endl; 
      cout << "  Number of cells with > 5 vertices: " 
           << num_large_cells << endl;
    }
    else {
      cout << "  Number of cells with > 4 vertices: " 
           << num_large_cells << endl;
    }
  }
}


template <typename MESH_TYPE>
void print_min_max_edge_lengths
(const MESH_TYPE & mesh, bool flag_more_info)
{
  typedef typename MESH_TYPE::HALF_EDGE_TYPE HALF_EDGE_TYPE;

  double min_edge_length_squared, max_edge_length_squared;
  int ihalf_edge_min, ihalf_edge_max;


  if (flag_more_info) {
    mesh.ComputeMinMaxEdgeLengthSquared
      (min_edge_length_squared, max_edge_length_squared,
       ihalf_edge_min, ihalf_edge_max);

    cout << "Min edge length: " << sqrt(min_edge_length_squared) << endl;

    if (mesh.IsHalfEdgeIndex(ihalf_edge_min) &&
        mesh.IsHalfEdgeIndex(ihalf_edge_max)) {

      const HALF_EDGE_TYPE * half_edge_min = 
        mesh.HalfEdge(ihalf_edge_min);
      const HALF_EDGE_TYPE * half_edge_max = 
        mesh.HalfEdge(ihalf_edge_max);

      cout << "  Length of edge: (";
      half_edge_min->PrintEndpoints(cout, ",");
      cout << ").  In cell: " << half_edge_min->Cell()->Index()
           << "." << endl;

      cout << "Max edge length: " << sqrt(max_edge_length_squared) << endl;
      cout << "  Length of edge: (";
      half_edge_max->PrintEndpoints(cout, ",");
      cout << ").  In cell: " << half_edge_max->Cell()->Index()
           << "." << endl;
    }
  }
  else {
    mesh.ComputeMinMaxEdgeLengthSquared(min_edge_length_squared,
                                        max_edge_length_squared);

    cout << "Min edge length: " << sqrt(min_edge_length_squared) << endl;
    cout << "Max edge length: " << sqrt(max_edge_length_squared) << endl;
  }
}


template <typename MESH_TYPE>
void print_min_cell_edge_length_ratio
(const MESH_TYPE & mesh, bool flag_more_info)
{
  typedef typename MESH_TYPE::HALF_EDGE_TYPE HALF_EDGE_TYPE;

  double min_edge_length_ratio_squared;
  double min_edge_length_squared, max_edge_length_squared;
  int icell_min_ratio;
  int ihalf_edge_min, ihalf_edge_max;

  if (flag_more_info) {
    mesh.ComputeMinCellEdgeLengthRatioSquared
      (min_edge_length_ratio_squared, icell_min_ratio,
       min_edge_length_squared, max_edge_length_squared,
       ihalf_edge_min, ihalf_edge_max);

    cout << "Min cell edge length ratio: " 
         << sqrt(min_edge_length_ratio_squared) << endl;
    if (mesh.IsCellIndex(icell_min_ratio)) {
      cout << "  In cell: " << icell_min_ratio << "." << endl;

      if (mesh.IsHalfEdgeIndex(ihalf_edge_min) &&
          mesh.IsHalfEdgeIndex(ihalf_edge_max)) {

        const HALF_EDGE_TYPE * half_edge_min =
          mesh.HalfEdge(ihalf_edge_min);
        const HALF_EDGE_TYPE * half_edge_max =
          mesh.HalfEdge(ihalf_edge_max);

        cout << "  Min cell edge length: "
             << sqrt(min_edge_length_squared) << ".";
        cout << "  Edge: (";
        half_edge_min->PrintEndpoints(cout, ",");
        cout << ")." << endl;

        cout << "  Max cell edge length: "
             << sqrt(max_edge_length_squared) << ".";
        cout << "  Edge: (";
        half_edge_max->PrintEndpoints(cout, ",");
        cout << ")." << endl;
      }
    }
  }
  else {
    mesh.ComputeMinCellEdgeLengthRatioSquared
      (min_edge_length_ratio_squared);

    cout << "Min cell edge length ratio: " 
         << sqrt(min_edge_length_ratio_squared) << endl;
  }
}


template <typename MESH_TYPE>
void print_min_max_angles
(const MESH_TYPE & mesh, bool flag_more_info)
{
  typedef typename MESH_TYPE::HALF_EDGE_TYPE HALF_EDGE_TYPE;

  double cos_min_angle, cos_max_angle;
  int ihalf_edge_min, ihalf_edge_max;

  if (flag_more_info) {
    mesh.ComputeCosMinMaxAngle
      (cos_min_angle, cos_max_angle,
       ihalf_edge_min, ihalf_edge_max);

    cout << "Minimum cell angle: "
         << acos(cos_min_angle)*180.0/M_PI << endl;

    if (mesh.IsHalfEdgeIndex(ihalf_edge_min)) {
      const HALF_EDGE_TYPE * half_edge_min = 
        mesh.HalfEdge(ihalf_edge_min);
      cout << "  At vertex " << half_edge_min->FromVertexIndex()
           << " in cell " << half_edge_min->Cell()->Index()
           << "." << endl;
    }

    cout << "Maximum cell angle: "
         << acos(cos_max_angle)*180.0/M_PI << endl;

    if (mesh.IsHalfEdgeIndex(ihalf_edge_max)) {
      const HALF_EDGE_TYPE * half_edge_max = 
        mesh.HalfEdge(ihalf_edge_max);
      cout << "  At vertex " << half_edge_max->FromVertexIndex()
           << " in cell " << half_edge_max->Cell()->Index()
           << "." << endl;
    }

  }
  else {
    mesh.ComputeCosMinMaxAngle(cos_min_angle, cos_max_angle);
    cout << "Minimum cell angle: "
         << acos(cos_min_angle)*180.0/M_PI << endl;
    cout << "Maximum cell angle: "
         << acos(cos_max_angle)*180.0/M_PI << endl;
  }
}



template <typename MESH_TYPE>
void print_manifold_info
(const MESH_TYPE & mesh, bool flag_more_info)
{
  typedef typename MESH_TYPE::HALF_EDGE_TYPE HALF_EDGE_TYPE;

  int iv, ihalf_edgeA, ihalf_edgeB;
  bool flag_non_manifold_vertex, flag_non_manifold_edge;

  const bool flag_manifold =
    mesh.CheckManifold
    (iv, ihalf_edgeA, flag_non_manifold_vertex, 
     flag_non_manifold_edge);
    
  const bool is_oriented =
    mesh.CheckOrientation(ihalf_edgeB);

  if (flag_manifold && is_oriented) {
    cout << "Mesh is an oriented manifold." << endl;
  }
  else if (flag_non_manifold_edge) {
    cout << "Mesh has a non-manifold edge";
    
    if (mesh.IsHalfEdgeIndex(ihalf_edgeA)) {
      const HALF_EDGE_TYPE * half_edgeA = mesh.HalfEdge(ihalf_edgeA);
      cout << " (";
      half_edgeA->PrintEndpoints(cerr, ",");
      cout << ")";
    }
    cout << "." << endl;
  }
  else if (flag_non_manifold_vertex && is_oriented) {
    cout << "Mesh has a non-manifold vertex " << iv << "." << endl;
  }
  else if (flag_non_manifold_vertex) {
    cout << "Non-manifold or inconsistent orientations at vertex "
         << iv << "." << endl;
  }
  else if (!is_oriented) {
    cout << "Mesh is a manifold." << endl;
    const HALF_EDGE_TYPE * half_edgeB = mesh.HalfEdge(ihalf_edgeB);
    const HALF_EDGE_TYPE * half_edgeBX = half_edgeB->NextHalfEdgeAroundEdge();
    cout << "Inconsistent orientations of cells "
         << half_edgeB->Cell()->Index() << " and "
         << half_edgeBX->Cell()->Index() << "." << endl;
  }

}


// *** SUBROUTINES ***

void parse_command_line(int argc, char ** argv)
{
  using namespace std;

  int iarg = 1;
  while (iarg < argc && argv[iarg][0] == '-') {

    std::string s = argv[iarg];
    if (s == "-more") 
      { flag_more_info = true; }
    else if (s == "-h")
      { help_msg(); }    
    else {
      cerr << "Usage error. Option " + s + " is undefined." << endl;
      usage_error();
    }
    iarg = iarg+1;
  }


  if (iarg >= argc || iarg+1 < argc)
    { usage_error(); }

  input_filename = argv[iarg];
}


void usage_msg(std::ostream & out)
{
  out << "Usage: meshinfo [OPTIONS] <input filename>" 
      << endl;
  cout << "OPTIONS: " << endl;
  cout << "  [-more] [-h]" << endl;
}


void help_msg()
{
  usage_msg(cout);
  cout << endl;
  cout << "meshinfo - Print mesh information." << endl;
  cout << "Options:" << endl;
  cout << "-more:    Print additional information." << endl;
  cout << "-h:       Output this help message and exit." << endl;
  exit(0);
}


void usage_error()
{
  usage_msg(cerr);
  exit(-1);
}
