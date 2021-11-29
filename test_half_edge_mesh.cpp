/// \file test_half_edge_mesh.hpp
/// Test HALF_EDGE_MESH data structure.

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

// Functions
void parse_command_line(int argc, char ** argv);
void usage_error();
void help_msg();
bool warn_non_manifold_or_not_oriented(const HALF_EDGE_MESH_A & mesh);
void print_time(const char * label, const double time);


int main(int argc, char ** argv)
{
  HALF_EDGE_MESH_A mesh;
  string error_msg;

  clock_t begin_time = clock();

  parse_command_line(argc, argv);

  open_and_read_off_file(input_filename, mesh);

  clock_t time2 = clock();

  if (!mesh.CheckAll(error_msg)) {
    cerr << "Error detected in mesh data structure." << endl;
    if (error_msg != "") { cerr << error_msg << endl; }
    exit(-1);
  }

  bool flag_warning = false;
  if (!flag_no_warn) 
    { flag_warning = warn_non_manifold_or_not_oriented(mesh); }

  if (!flag_silent) {
    if (!flag_warning) 
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

  open_and_write_off_file(output_filename, mesh);

  clock_t end_time = clock();

  if (flag_time) {
    print_time("Time to read file:  ", time2 - begin_time);
    print_time("Time to check mesh: ", time3 - time2);
    print_time("Time to write file: ", end_time - time3);
    print_time("Total time:         ", end_time - begin_time);
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
bool warn_non_manifold_or_not_oriented(const HALF_EDGE_MESH_A & mesh)
{
  int iv, ihalf_edgeA, ihalf_edgeB;
  bool flag_non_manifold_vertex, flag_non_manifold_edge;

    const bool flag_manifold =
      mesh.CheckManifold
      (iv, ihalf_edgeA, flag_non_manifold_vertex, 
       flag_non_manifold_edge);
    
    if (flag_non_manifold_edge) {
      cerr << "Warning: Non-manifold edge (";
      mesh.HalfEdge(ihalf_edgeA)->PrintEndpoints(cerr, ",");
      cerr << ")." << endl;

      // Non-manifold edge automatically implies inconsistent orientations.
      return true;
    }

    const bool is_oriented =
      mesh.CheckOrientation(ihalf_edgeB);

    if (is_oriented) {
      if (flag_non_manifold_vertex) {
        cerr << "Warning: Non-manifold vertex " << iv << "." << endl;

        return true;
      }
    }
    else {
      cerr << "Warning: Inconsistent orientation of cells incident on edge (";
      mesh.HalfEdge(ihalf_edgeB)->PrintEndpoints(cerr, ",");
      cerr << ")." << endl;

      if (flag_non_manifold_vertex) {
        cerr << "Warning: Non-manifold vertex or inconsistent orientations in cells incident on vertex " << iv << "." << endl;
      }

      return true;
    }

    return false;
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
  cout << "test_half_edge_mesh - Test the HALF_EDGE_MESH and associated classes" << endl;
  cout << "  and I/O routines by reading a .off file to the mesh," << endl;
  cout << "  running check mesh, check manifold and check orientation routines" << endl;
  cout << "  and then writing the mesh to a .off file." << endl;
  cout << endl;
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

