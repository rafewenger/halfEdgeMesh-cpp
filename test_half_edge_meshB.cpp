/// \file test_half_edge_meshB.hpp
/// Test HALF_EDGE_MESH data structure using classes derived from
/// VERTEX_BASE, HALF_EDGE_BASE, and CELL_BASE.
///
/// - Reads in a mesh.
/// - Stores "new values" in the new mesh vertices, half edges and cells.
/// - Creates a new mesh whose vertex and cell identifiers equal
///   the new values stored in the old mesh.
/// - Writes the new mesh to a .off file.
/// - The new .off file will have 5 vertices with coordinates (0 0 0)
///   since the new mesh has no vertices equal to 0, 1, 2, 3, or 4.

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>

#include "half_edge_mesh.hpp"
#include "half_edge_mesh_IO.hpp"

using namespace HMESH;
using namespace std;

// Global variables
const char * input_filename(NULL);
const char * output_filename(NULL);
bool flag_silent(false);
bool flag_no_warn(false);
bool flag_time(false);

// New classes.
class VERTEX_B;
class HALF_EDGE_B;
class CELL_B;
typedef VERTEX_B * VERTEX_B_PTR;
typedef HALF_EDGE_B * HALF_EDGE_B_PTR;
typedef CELL_B * CELL_B_PTR;

class VERTEX_B:public VERTEX3D_BASE<HALF_EDGE_B_PTR,float>
{
  // New field.
public:
  float new_val;

  //  typedef float COORD_TYPE;
};

class HALF_EDGE_B:public HALF_EDGE_BASE<VERTEX_B_PTR,HALF_EDGE_B_PTR,CELL_B_PTR>
{
  // New field.
public: 
  float new_val;
};

class CELL_B:public CELL_BASE<HALF_EDGE_B_PTR>
{
  // New field.
public:
  float new_val;
};

typedef HALF_EDGE_MESH_BASE<VERTEX_B, HALF_EDGE_B, CELL_B>
HALF_EDGE_MESH_B;



// Functions
void parse_command_line(int argc, char ** argv);
void usage_error();
void help_msg();
void print_time(const char * label, const double time);
bool warn_non_manifold_or_not_oriented(const HALF_EDGE_MESH_B & mesh);
bool check_mesh(const HALF_EDGE_MESH_B & mesh, const bool flag_no_warn);
void test_new_classes(HALF_EDGE_MESH_B & mesh);
void write_mesh_info(const HALF_EDGE_MESH_B & mesh);
void create_new_mesh(const HALF_EDGE_MESH_B & mesh,
                     HALF_EDGE_MESH_B & new_mesh);


int main(int argc, char ** argv)
{
  HALF_EDGE_MESH_B mesh, meshB;

  clock_t begin_time = clock();

  parse_command_line(argc, argv);

  open_and_read_off_file(input_filename, mesh);

  clock_t time2 = clock();

  check_mesh(mesh, flag_no_warn);

  test_new_classes(mesh);
  check_mesh(mesh, flag_no_warn);
  write_mesh_info(mesh);

  create_new_mesh(mesh, meshB);

  bool passed_check = check_mesh(meshB, flag_no_warn);

  if (!flag_silent) {
    if (passed_check) 
      { cout << "Mesh data structure passed check." << endl; }
  }

  clock_t time3 = clock();

  if (output_filename == NULL) 
    { output_filename = "out.off"; }
  if (std::string(output_filename) == std::string(input_filename))
    { output_filename = "out2.off"; }

  if (!flag_silent) {
    cout << "Writing file: " << output_filename << endl;
  }

  open_and_write_off_file(output_filename, meshB);

  clock_t end_time = clock();

  if (flag_time) {
    print_time("Time to read file:  ", time2 - begin_time);
    print_time("Time to check mesh: ", time3 - time2);
    print_time("Time to write file: ", end_time - time3);
    print_time("Total time:         ", end_time - begin_time);
  }
}


// *** SUBROUTINES ***

// Test new classes.
void test_new_classes(HALF_EDGE_MESH_B & mesh)
{
  for (int iv = 0; iv < mesh.VertexListLength(); iv++) {
    VERTEX_B_PTR v = mesh.VertexNC(iv);
    if (v != NULL) {
      v->new_val = v->Index() + 5;
    }
  }

  for (int ihalf_edge = 0; ihalf_edge < mesh.HalfEdgeListLength(); 
       ihalf_edge++) {
    HALF_EDGE_B_PTR half_edge = mesh.HalfEdgeNC(ihalf_edge);
    if (half_edge != NULL) {
      half_edge->new_val = half_edge->Index() + 5;
    }
  }

  for (int icell = 0; icell < mesh.CellListLength(); 
       icell++) {
    CELL_B_PTR cell = mesh.CellNC(icell);
    if (cell != NULL) {
      cell->new_val = cell->Index() + 5;
    }
  }

}


void write_mesh_info(const HALF_EDGE_MESH_B & mesh)
{
  const int iv = mesh.VertexListLength()/2;
  const VERTEX_B * v = mesh.Vertex(iv);
  if (v != NULL) {
    cout << "Vertex(" << iv << ").new_val = " << v->new_val << endl;
  }

  const int ihalf_edge = mesh.HalfEdgeListLength()/2;
  const HALF_EDGE_B * half_edge = mesh.HalfEdge(ihalf_edge);
  if (half_edge != NULL) {
    cout << "Half edge(" << ihalf_edge << ").new_val = " 
         << half_edge->new_val << endl;
  }

  const int icell = mesh.CellListLength()/2;
  const CELL_B * cell = mesh.Cell(icell);
  if (cell != NULL) {
    cout << "Cell(" << icell << ").new_val = " 
         << cell->new_val << endl;
  }
}

// Create new mesh using indices at new_val.
// @pre new_mesh is empty.
void create_new_mesh(const HALF_EDGE_MESH_B & mesh,
                     HALF_EDGE_MESH_B & new_mesh)
{
  for (int iv = 0; iv < mesh.VertexListLength(); iv++) {
    const VERTEX_B * v = mesh.Vertex(iv);
    const int ivnew = v->new_val;
    new_mesh.AddVertex(ivnew);
    new_mesh.SetCoord(ivnew, v->coord);
  }

  for (int icell = 0; icell < mesh.CellListLength(); icell++) {
    const CELL_B * cell = mesh.Cell(icell);
    const HALF_EDGE_B * half_edge = cell->HalfEdge();
    std::vector<int> cell_vlist;
    for (int j = 0; j < cell->NumVertices(); j++)  {
      const VERTEX_B * v = half_edge->FromVertex();
      const int iv_new = v->new_val;
      cell_vlist.push_back(iv_new);
      half_edge = half_edge->NextHalfEdgeInCell();
    }
    new_mesh.AddCell(cell_vlist);
  }
}


void parse_command_line(int argc, char ** argv)
{
  using namespace std;

  int iarg = 1;
  while (iarg < argc && argv[iarg][0] == '-') {

    std::string s = argv[iarg];
    if (s == "-s")
      { flag_silent = true; }
    else if (s == "-no_warn") 
      { flag_no_warn = true; }
    else if (s == "-time")
      { flag_time = true; }
    else if (s == "-h")
      { help_msg(); }    
    else {
      cerr << "Usage error. Option " + s + " is undefined." << endl;
      usage_error();
    }
    iarg = iarg+1;
  }


  if (iarg >= argc || iarg+2 < argc)
    { usage_error(); }

  input_filename = argv[iarg];

  if (iarg+1 < argc) 
    { output_filename = argv[iarg+1]; }
}


/// Print a warning message if the mesh is not a manifold or not oriented.
/// - Returns true if mesh is not a manifold or not oriented
bool warn_non_manifold_or_not_oriented(const HALF_EDGE_MESH_B & mesh)
{
    int iv, ihalf_edge;
    bool flag_non_manifold_vertex, flag_non_manifold_edge;
    
    if (!mesh.CheckManifold
        (iv, ihalf_edge, flag_non_manifold_vertex, 
         flag_non_manifold_edge)) {

      if (flag_non_manifold_vertex) {
        cerr << "Warning: Non-manifold vertex " << iv << "." << endl;
      }

      if (flag_non_manifold_edge) {
        cerr << "Warning: Non-manifold edge (";
        mesh.HalfEdge(ihalf_edge)->PrintEndpoints(cerr, ",");
        cerr << ")." << endl;
      }

      return true;
    }

    if (!flag_non_manifold_edge) {
      if (!mesh.CheckOrientation(ihalf_edge)) {
        cerr << "Warning: Inconsistent orientation of cells incident on edge (";
        mesh.HalfEdge(ihalf_edge)->PrintEndpoints(cerr, ",");
        cerr << ")." << endl;

        return true;
      }
    }

    return false;
}


// Return true if mesh passed mesh check, manifold check 
//   and orientation check.
bool check_mesh(const HALF_EDGE_MESH_B & mesh, const bool flag_no_warn)
{
  string error_msg;

  if (!mesh.CheckAll(error_msg)) {
    cerr << "Error detected in mesh data structure." << endl;
    if (error_msg != "") { cerr << error_msg << endl; }
    exit(-1);
  }

  bool flag_warning = false;
  if (!flag_no_warn) {
    flag_warning = warn_non_manifold_or_not_oriented(mesh); 
  }

  return (!flag_warning);
}


void print_time(const char * label, const double time)
{
  cout << label << time/CLOCKS_PER_SEC << " seconds." << endl;
}

void usage_msg(std::ostream & out)
{
  out << "Usage: test_half_edge_mesh [-s] [-no_warn] [-time] [-h] <input filename> [<output filename>]" 
      << endl;
}


void help_msg()
{
  usage_msg(cout);
  cout << endl;
  cout << "test_half_edge_meshB - Test derived vertex, half edge and cell classes in the half edge mesh data structure." << endl;
  cout << "  Reads in a mesh." << endl;
  cout << "  Stores \"new values\" in the new mesh vertices, half edges and cells." << endl;
  cout << "  Creates a new mesh whose vertex and cell identifiers equal"
       << endl;
  cout << "    the new values stored in the old mesh." << endl;
  cout << "  Writes the new mesh to a .off file." << endl;
  cout << "  Note: The new .off file will have 5 vertices with coordintes (0, 0, 0)" << endl;
  cout << "    since the new mesh has no vertices equal to 0, 1, 2, 3, or 4." << endl;
  cout << "Options:" << endl;
  cout << "-s:       Silent. Output only warnings and error messages." << endl;
  cout << "-no_warn: Do not output non-manifold or inconsistent orientation warnings." << endl;
  cout << "-time:    Report run time." << endl;
  cout << "-h:       Output this help message and exit." << endl;
  exit(0);
}


void usage_error()
{
  usage_msg(cerr);
  exit(-1);
}

