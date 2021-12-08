/// \file decimate_mesh.hpp
/// Some simple mesh decimation routines.
/// - Uses data structure HALF_EDGE_MESH_DCMT_BASE (DCMT = decimate).

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>

#include "half_edge_mesh.hpp"
#include "half_edge_mesh_DCMT.hpp"
#include "half_edge_mesh_IO.hpp"

using namespace HMESH;
using namespace std;

// Global filename variables
const char * input_filename(NULL);
const char * output_filename(NULL);

// Global variables controlling output/checks.
bool flag_terse(false);
bool flag_silent(false);
bool flag_no_warn(false);
bool flag_time(false);
bool flag_reduce_checks(false);

// Global variables controlling decimation.
bool flag_collapse_edges(false);
bool flag_collapse_short_edges(false);
bool flag_split_cells(false);
bool flag_split_all_cells(false);
bool flag_join_cells(false);
bool flag_join_each_cell(false);
bool flag_split_edges(false);
bool flag_split_long_edges(false);
bool flag_allow_non_manifold(false);
bool flag_fail_on_non_manifold(false);

// Number of cells for data set to be considered large.
const int LARGE_DATA_NUM_CELLS(10000);

// Parse/print/prompt/check functions.
void parse_command_line(int argc, char ** argv);
const HALF_EDGE_DCMT_A *
prompt_for_mesh_edge(const HALF_EDGE_MESH_DCMT_A & mesh,
                     const bool flag_only_internal);
void usage_error();
void help_msg();
void print_time(const char * label, const double time);
void print_mesh_info(const HALF_EDGE_MESH_DCMT_A & mesh);
template <typename MESH_TYPE>
bool check_mesh(const MESH_TYPE & mesh, const bool flag_no_warn);

// Decimation functions.
void prompt_and_split_edges
(HALF_EDGE_MESH_DCMT_A & mesh, 
 const bool flag_terse, const bool flag_no_warn);
void split_longest_edge_in_each_cell
(HALF_EDGE_MESH_DCMT_A & mesh, const bool flag_terse, 
 const bool flag_no_warn);
void prompt_and_collapse_edges
(HALF_EDGE_MESH_DCMT_A & mesh, 
 const bool flag_terse, const bool flag_no_warn);
void collapse_shortest_edge_in_each_cell
(HALF_EDGE_MESH_DCMT_A & mesh, const bool flag_terse, 
 const bool flag_no_warn);
void prompt_and_split_cells
(HALF_EDGE_MESH_DCMT_A & mesh, const bool flag_terse,
 const bool flag_no_warn);
void split_all_cells
(HALF_EDGE_MESH_DCMT_A & mesh, const bool flag_terse, 
 const bool flag_no_warn);
void prompt_and_join_cells
(HALF_EDGE_MESH_DCMT_A & mesh, const bool flag_terse, 
 const bool flag_no_warn);
void join_each_cell
(HALF_EDGE_MESH_DCMT_A & mesh, const bool flag_terse, 
 const bool flag_no_warn);


int main(int argc, char ** argv)
{
  HALF_EDGE_MESH_DCMT_A mesh;

  clock_t begin_time = clock();

  parse_command_line(argc, argv);

  open_and_read_off_file(input_filename, mesh);

  clock_t time2 = clock();

  try {

    if (flag_split_edges) 
      { prompt_and_split_edges(mesh, flag_terse, flag_no_warn); }

    if (flag_collapse_edges) 
      { prompt_and_collapse_edges(mesh, flag_terse, flag_no_warn); }

    if (flag_split_cells) 
      { prompt_and_split_cells(mesh, flag_terse, flag_no_warn); }

    if (flag_join_cells) 
      { prompt_and_join_cells(mesh, flag_terse, flag_no_warn); }

    if (flag_collapse_short_edges) {
      collapse_shortest_edge_in_each_cell
        (mesh, flag_terse, flag_no_warn); 
    }
    
    if (flag_split_long_edges)
      { split_longest_edge_in_each_cell(mesh, flag_terse, flag_no_warn); }

    if (flag_split_all_cells)
      { split_all_cells(mesh, flag_terse, flag_no_warn); }

    if (flag_join_each_cell)
      { join_each_cell(mesh, flag_terse, flag_no_warn); }
  }
  catch (SIMPLE_EXCEPTION & e) {
    cerr << e.what() << endl;
    cerr << "Exiting." << endl;
    exit(-1);
  }


  const bool passed_check = check_mesh(mesh, flag_silent && flag_no_warn);

  if (!flag_silent) {
    if (passed_check) 
      { cout << "Mesh data structure passed check." << endl; }
  }

  if (!flag_silent) {
    cout << endl;
    print_mesh_info(mesh); 
  }
    
  clock_t time3 = clock();

  if (output_filename == NULL) 
    { output_filename = "out.off"; }
  if (std::string(output_filename) == std::string(input_filename))
    { output_filename = "out2.off"; }

  if (!flag_silent) {
    cout << endl;
    cout << "Writing file: " << output_filename << endl;
  }

  open_and_write_off_file(output_filename, mesh);

  clock_t end_time = clock();

  if (flag_time) {
    print_time("Time to read file:  ", time2 - begin_time);
    print_time("Time to process mesh: ", time3 - time2);
    print_time("Time to write file: ", end_time - time3);
    print_time("Total time:         ", end_time - begin_time);
  }
}


// *** Collapse edge routines ***

bool check_edge_collapse(const HALF_EDGE_MESH_DCMT_A & mesh, 
                         const HALF_EDGE_DCMT_A * half_edge,
                         const bool flag_no_warn);
template <typename MESH_TYPE>
void reduce_checks_on_large_datasets
(const MESH_TYPE &  mesh, const bool flag_no_warn);


// Collapse edge.
void collapse_edge
(HALF_EDGE_MESH_DCMT_A & mesh, const HALF_EDGE_DCMT_A * half_edge,
 const bool flag_terse, const bool flag_no_warn,
 const bool flag_check)
{
  const bool flag = check_edge_collapse(mesh, half_edge, flag_no_warn);

  if (mesh.IsIllegalEdgeCollapse(half_edge)) { return; }

  if (flag || flag_allow_non_manifold) {
    if (!flag_terse) {
      cout << "Collapsing edge (" << half_edge->EndpointsStr(",")
           << ")." << endl;
    }

    const VERTEX_DCMT_A * vnew = 
      mesh.CollapseEdge(half_edge->Index());
    if (vnew == NULL) {
      cout << "Skipped illegal collapse of edge (" 
           << half_edge->EndpointsStr(",") << ")." << endl;
    }

    if (flag_check)
      { check_mesh(mesh, flag_no_warn); }
  }
  else {
    if (!flag_terse) {
      cout << "Skipped collapse of edge (" 
           << half_edge->EndpointsStr(",") << ")." << endl;
    }
  }
}


// Prompt and collapse edges.
void prompt_and_collapse_edges
(HALF_EDGE_MESH_DCMT_A & mesh, 
 const bool flag_terse, const bool flag_no_warn)
{
  std::string error_msg;

  while (true) {
    const HALF_EDGE_DCMT_A * half_edge0 =
      prompt_for_mesh_edge(mesh, false);

    if (half_edge0 == NULL) {
      // End.
      cout << endl;
      return;
    }

    collapse_edge
      (mesh, half_edge0, flag_terse, flag_no_warn, true);

    cout << endl;
  }
}


// Collapse shortest cell edge.
void collapse_shortest_cell_edge
(HALF_EDGE_MESH_DCMT_A & mesh, const int icell,
 const bool flag_terse, const bool flag_no_warn, 
 const bool flag_check)
{
  int ihalf_edge_min, ihalf_edge_max;
  double minL, maxL;

  const CELL_DCMT_A * cell = mesh.Cell(icell);
  if (cell == NULL) { return; }

  cell->ComputeMinMaxEdgeLengthSquared
    (minL, maxL, ihalf_edge_min, ihalf_edge_max);

  const HALF_EDGE_DCMT_A * half_edge_min = 
    mesh.HalfEdge(ihalf_edge_min);

  collapse_edge
    (mesh, half_edge_min, flag_terse, flag_no_warn, flag_check);
}


// Collapse shortest edge in each cell.
void collapse_shortest_edge_in_each_cell
(HALF_EDGE_MESH_DCMT_A & mesh, const bool flag_terse, const bool flag_no_warn)
{
  const int n = mesh.CellListLength();

  reduce_checks_on_large_datasets(mesh, flag_no_warn);
  bool flag_check = !flag_reduce_checks;

  for (int icell = 0; icell < n; icell++) {
    const CELL_DCMT_A * cell = mesh.Cell(icell);
    if (cell == NULL) { continue; }

    collapse_shortest_cell_edge
      (mesh, icell, flag_terse, flag_no_warn, flag_check);

    if (flag_reduce_checks) {
      // Check mesh halfway through.
      if (icell == int(n/2)) 
        { check_mesh(mesh, flag_no_warn); }
    }
  }

}


// *** Split edge routines ***

// Split edge.
void split_edge
(HALF_EDGE_MESH_DCMT_A & mesh,
 const HALF_EDGE_DCMT_A * half_edge,
 const bool flag_terse, const bool flag_no_warn,
 const bool flag_check)
{
  if (!flag_terse) {
    cout << "Splitting edge (" << half_edge->EndpointsStr(",")
         << ")." << endl;
  }

  const VERTEX_DCMT_A * vnew = 
    mesh.SplitEdge(half_edge->Index());

  if (vnew == NULL) {
    cout << "Split of edge (" << half_edge->EndpointsStr(",")
         << " failed." << endl;
  }

  if (flag_check)
    { check_mesh(mesh, flag_no_warn); }
}


// Split edges.
void prompt_and_split_edges
(HALF_EDGE_MESH_DCMT_A & mesh, 
 const bool flag_terse, const bool flag_no_warn)
{
  while (true) {
    const HALF_EDGE_DCMT_A * half_edge0 =
      prompt_for_mesh_edge(mesh, false);

    if (half_edge0 == NULL) {
      // End.
      cout << endl;
      return;
    }

    split_edge(mesh, half_edge0, flag_terse, flag_no_warn, true);

    cout << endl;
  }
}


// Split longest cell edge.
void split_longest_cell_edge
(HALF_EDGE_MESH_DCMT_A & mesh, const int icell,
 const bool flag_terse, const bool flag_no_warn, 
 const bool flag_check)
{
  int ihalf_edge_min, ihalf_edge_max;
  double minL, maxL;

  const CELL_DCMT_A * cell = mesh.Cell(icell);
  if (cell == NULL) { return; }

  cell->ComputeMinMaxEdgeLengthSquared
    (minL, maxL, ihalf_edge_min, ihalf_edge_max);

  const HALF_EDGE_DCMT_A * half_edge_max =
    mesh.HalfEdge(ihalf_edge_max);

  split_edge
    (mesh, half_edge_max, flag_terse, flag_no_warn, flag_check);
}


// Split longest edge in each cell.
void split_longest_edge_in_each_cell
(HALF_EDGE_MESH_DCMT_A & mesh, const bool flag_terse, const bool flag_no_warn)
{
  const int n = mesh.CellListLength();
  const int original_numv = mesh.VertexListLength();

  reduce_checks_on_large_datasets(mesh, flag_no_warn);
  bool flag_check = !flag_reduce_checks;

  for (int icell = 0; icell < n; icell++) {
    const CELL_DCMT_A * cell = mesh.Cell(icell);
    if (cell == NULL) { continue; }

    split_longest_cell_edge
      (mesh, icell, flag_terse, flag_no_warn, flag_check);

    if (flag_reduce_checks) {
      // Check mesh halfway through.
      if (icell == int(n/2)) 
        { check_mesh(mesh, flag_no_warn); }
    }
  }
}


// *** Split cell routines ***

bool check_split_cell(const HALF_EDGE_MESH_DCMT_A & mesh, 
                      const HALF_EDGE_DCMT_A * half_edgeA,
                      const HALF_EDGE_DCMT_A * half_edgeB,
                      const bool flag_no_warn);
void print_cells_with_more_than_three_vertices
(const int max_num, const std::vector<int> & cell_list);



// Split cell with diagonal 
//   (half_edgeA->FromVertex()), half_edgeB->FromVertex()).
// - Returns split edge.
// - Returns NULL if split fails.
const HALF_EDGE_DCMT_A * split_cell
(HALF_EDGE_MESH_DCMT_A & mesh, 
 const HALF_EDGE_DCMT_A * half_edgeA,
 const HALF_EDGE_DCMT_A * half_edgeB,
 const bool flag_terse, const bool flag_no_warn,
 const bool flag_check)
{
  const int ihalf_edgeA = half_edgeA->Index();
  const int ihalf_edgeB = half_edgeB->Index();
  const VERTEX_DCMT_A * vA = half_edgeA->FromVertex();
  const VERTEX_DCMT_A * vB = half_edgeB->FromVertex();
  const int ivA = vA->Index();
  const int ivB = vB->Index();
  const int icell = half_edgeA->CellIndex();

  const bool flag = 
    check_split_cell(mesh, half_edgeA, half_edgeB, flag_no_warn);

  if (mesh.IsIllegalSplitCell(half_edgeA, half_edgeB)) 
    { return NULL; }

  if (flag || flag_allow_non_manifold) {
    if (!flag_terse) {
      cout << "Splitting cell " << icell << " at diagonal (" 
           << ivA << "," << ivB << ")." << endl;
    }

    const HALF_EDGE_DCMT_A * split_edge =
      mesh.SplitCell(ihalf_edgeA, ihalf_edgeB);
    if (split_edge == NULL) {
      cout << "Split of cell " << icell << " at diagonal (" 
           << vA->Index() << "," << vB->Index() 
           << ") failed." << endl;
      return NULL;
    }

    if (flag_check)
      { check_mesh(mesh, flag_no_warn); }

    return split_edge;
  }
  else {
    if (!flag_terse) {
      cout << "Skipping split of cell " << icell
           << " at diagonal (" << ivA << "," << ivB 
           << ")." << endl;
    }
    return NULL;
  }
}


// Get at most max_num cells with more than three vertices.
void get_cells_with_more_than_three_vertices
(const HALF_EDGE_MESH_DCMT_A & mesh,
 const int max_num, std::vector<int> & cell_list)
{
  const int THREE(3);

  cell_list.clear();
  if (max_num < 1) { return; }

  for (int icell = 0; icell < mesh.CellListLength(); icell++) {
    const CELL_DCMT_A * cell = mesh.Cell(icell);
    if (cell == NULL) { continue; }

    if (cell->NumVertices() > THREE) 
      { cell_list.push_back(icell); }

    if (cell_list.size() >= max_num)
      { return; }
  }
}


// Prompt and split cells.
void prompt_and_split_cells
(HALF_EDGE_MESH_DCMT_A & mesh, 
 const bool flag_terse, const bool flag_no_warn)
{
  const int THREE(3);
  const int MAX_NUM(10);
  int icell;
  std::vector<int> cell_list;

  get_cells_with_more_than_three_vertices
    (mesh, MAX_NUM, cell_list);

  if (cell_list.size() == 0) {

    if (!flag_no_warn) {
      cout << "All cells are triangles. No cells can be split." << endl;
    }

    return;
  }

  print_cells_with_more_than_three_vertices(MAX_NUM, cell_list);
  cout << endl;

  while (true) {

    cout << "Enter cell (-1 to end): ";
    cin >> icell;

    if (icell < 0) { return; }
    if (icell >= mesh.CellListLength()) {
      cout << "No cell has index " << icell << "." << endl;
      cout << "Maximum cell index: " << mesh.CellListLength()-1 << endl;
      continue;
    }

    const CELL_DCMT_A * cell = mesh.Cell(icell);
    if (cell == NULL) {
      cout << "No cell has index " << icell << "." << endl;
      cout << endl;
      continue;
    }

    if (cell->NumVertices() <= THREE) {
      cout << "Cell " << icell << " has fewer than four vertices." << endl;
      cout << endl;
      continue;
    }

    cout << "Vertices in cell " << icell << ":";
    const HALF_EDGE_DCMT_A * half_edge = cell->HalfEdge();
    for (int k = 0; k < cell->NumVertices(); k++) {
      cout << "  " << half_edge->FromVertexIndex();
      half_edge = half_edge->NextHalfEdgeInCell();
    }
    cout << endl;

    int ivA, ivB;
    cout << "Enter two distinct vertex indices (-1 to end): ";
    cin >> ivA;
    if (ivA < 0) { return; }
    cin >> ivB;
    if (ivB < 0) { return; }
    const HALF_EDGE_DCMT_A * half_edgeA = NULL;
    const HALF_EDGE_DCMT_A * half_edgeB = NULL;

    if (ivA == ivB) {
      cout << endl;
      cout << "Vertices are not dictinct. Start again." << endl;
      cout << endl;
      continue;
    }

    half_edge =  cell->HalfEdge();
    for (int k = 0; k < cell->NumVertices(); k++) {
      if (half_edge->FromVertexIndex() == ivA)
        { half_edgeA = half_edge; }
      if (half_edge->FromVertexIndex() == ivB)
        { half_edgeB = half_edge; }

      half_edge = half_edge->NextHalfEdgeInCell();
    }

    if (half_edgeA == NULL || half_edgeB == NULL) {
      cout << endl;
      cout << "Vertices are not in cell " << icell << "." << endl;
      cout << "Start again." << endl;
      cout << endl;
      continue;
    }

    if (half_edgeA->ToVertexIndex() == ivB ||
        half_edgeB->ToVertexIndex() == ivA) {
      cout << endl;
      cout << "(" << ivA << "," << ivB 
           << ") is a cell edge, not a cell diagonal."
           << endl;
      cout << "  Vertices must not be adjacent." << endl;
      cout << "Start again." << endl;
      cout << endl;
      continue;
    }

    split_cell(mesh, half_edgeA, half_edgeB, 
               flag_terse, flag_no_warn, true);

    cout << endl;
  }
}


// Split cell at largest angle.
// - Split cell at vertex forming the largest angle.
// - Split as many times as necessary to triangulate.
void split_cell_at_largest_angle
(HALF_EDGE_MESH_DCMT_A & mesh, const CELL_DCMT_A * cell, 
 const bool flag_terse, const bool flag_no_warn,
 const bool flag_check)
{
  int ihalf_edge_min, ihalf_edge_max;
  double cos_min_angle, cos_max_angle;

  cell->ComputeCosMinMaxAngle(cos_min_angle, cos_max_angle,
                              ihalf_edge_min, ihalf_edge_max);
  const HALF_EDGE_DCMT_A * half_edgeA = mesh.HalfEdge(ihalf_edge_max);

  while (!half_edgeA->Cell()->IsTriangle()) {

    const HALF_EDGE_DCMT_A * half_edgeB =
      half_edgeA->PrevHalfEdgeInCell()->PrevHalfEdgeInCell();
    const VERTEX_DCMT_A * vA = half_edgeA->FromVertex();
    const VERTEX_DCMT_A * vB = half_edgeB->FromVertex();
    const int ivA = vA->Index();
    const int ivB = vB->Index();

    const HALF_EDGE_DCMT_A * split_edge =
      split_cell(mesh, half_edgeA, half_edgeB, 
                 flag_terse, flag_no_warn, flag_check);
    if (split_edge == NULL) {
      // Cannot split cell at largest angle.
      return;
    }

    if (flag_check)
      { check_mesh(mesh, flag_no_warn); }

    // Get largest angle in remaining cell.
    const CELL_DCMT_A * cellA = split_edge->Cell();
    cell->ComputeCosMinMaxAngle(cos_min_angle, cos_max_angle,
                                ihalf_edge_min, ihalf_edge_max);
    half_edgeA = mesh.HalfEdge(ihalf_edge_max);
  }

}


// Split all cells.
void split_all_cells
(HALF_EDGE_MESH_DCMT_A & mesh, const bool flag_terse, 
 const bool flag_no_warn)
{
  const int n = mesh.CellListLength();

  reduce_checks_on_large_datasets(mesh, flag_no_warn);
  bool flag_check = !flag_reduce_checks;

  for (int icell = 0; icell < n; icell++) {
    const CELL_DCMT_A * cell = mesh.Cell(icell);
    if (cell == NULL) { continue; }
    if (cell->IsTriangle()) { continue; }

    split_cell_at_largest_angle
      (mesh, cell, flag_terse, flag_no_warn, flag_check);

    if (flag_reduce_checks) {
      // Check mesh halfway through.
      if (icell == int(n/2)) 
        { check_mesh(mesh, flag_no_warn); }
    }
  }
}


// *** Join cell routines ***

bool check_join_cell(const HALF_EDGE_MESH_DCMT_A & mesh, 
                     const HALF_EDGE_DCMT_A * half_edge,
                     const bool flag_no_warn);


// Join two cells separated half_edge.
void join_two_cells(HALF_EDGE_MESH_DCMT_A & mesh, 
                    const HALF_EDGE_DCMT_A * half_edge, 
                    const bool flag_terse, const bool flag_no_warn,
                    const bool flag_check)
{
  const HALF_EDGE_DCMT_A * half_edgeX =
    half_edge->NextHalfEdgeAroundEdge();

  const bool flag = 
    check_join_cell(mesh, half_edge, flag_no_warn);

  if (mesh.IsIllegalJoinCells(half_edge)) { return; }

  if (flag) {
    if (!flag_terse) {
      cout << "Joining cell " << half_edge->CellIndex()
           << " to cell " << half_edgeX->CellIndex() 
           << " by deleting edge (" 
           << half_edge->EndpointsStr(",") << ")." << endl;
    }

    const HALF_EDGE_DCMT_A * half_edgeB =
      mesh.JoinTwoCells(half_edge->Index());

    if (half_edgeB == NULL) {
      cout << "Join of cell " << half_edge->CellIndex()
           << " to cell " << half_edgeX->CellIndex() 
           << " failed." << endl;
    }
    else {
      if (flag_check)
        { check_mesh(mesh, flag_no_warn); }

      return;
    }
  }
  else {
    if (!flag_terse) {
      cout << "Skipping join of cell " << half_edge->CellIndex()
           << " with cell " << half_edgeX->CellIndex() << "." << endl;
    }
  }

}


// Prompt and join cells.
void prompt_and_join_cells
(HALF_EDGE_MESH_DCMT_A & mesh, 
 const bool flag_terse, const bool flag_no_warn)
{
  while (true) {
    const HALF_EDGE_DCMT_A * half_edge0 =
      prompt_for_mesh_edge(mesh, true);

    if (half_edge0 == NULL) {
      // End.
      cout << endl;
      return;
    }

    join_two_cells(mesh, half_edge0, flag_terse, flag_no_warn, true);

    cout << endl;
  }
}


// Attempt to join each cell by deleting longest edge.
void join_each_cell
(HALF_EDGE_MESH_DCMT_A & mesh, const bool flag_terse, const bool flag_no_warn)
{
  const int MAX_NUMV(6);
  const int n = mesh.CellListLength();

  reduce_checks_on_large_datasets(mesh, flag_no_warn);
  bool flag_check = !flag_reduce_checks;

  for (int icell = 0; icell < n; icell++) {
    const CELL_DCMT_A * cell = mesh.Cell(icell);
    if (cell == NULL) { continue; }

    if (cell->NumVertices() >= MAX_NUMV) {
      // Don't let cells get too large.
      continue;
    }

    double Lmin, Lmax;
    int ihalf_edge_min, ihalf_edge_max;
    cell->ComputeMinMaxEdgeLengthSquared
      (Lmin, Lmax, ihalf_edge_min, ihalf_edge_max);

    // CORRECTED: 12-07-2021 - RW
    //const HALF_EDGE_DCMT_A * half_edge = mesh.HalfEdge(ihalf_edge_min);
 
    const HALF_EDGE_DCMT_A * half_edge = mesh.HalfEdge(ihalf_edge_max);
    
    const HALF_EDGE_DCMT_A * half_edgeX = 
      half_edge->NextHalfEdgeAroundEdge();

    if (half_edgeX->Cell()->NumVertices() >= MAX_NUMV) {
      // Don't let cells get too large.
      continue;
    }
    
    join_two_cells
      (mesh, half_edge, flag_terse, flag_no_warn, flag_check);


    if (flag_reduce_checks) {
      // Check mesh halfway through.
      if (icell == int(n/2)) 
        { check_mesh(mesh, flag_no_warn); }
    }
  }
}


// *** Check routines ***

template <typename MESH_TYPE>
bool check_oriented_manifold
(const MESH_TYPE & mesh, const bool flag_no_warn);


// Return true if mesh passed mesh check, manifold check 
//   and orientation check.
template <typename MESH_TYPE>
bool check_mesh(const MESH_TYPE & mesh, const bool flag_no_warn)
{
  string error_msg;

  if (!mesh.CheckAll(error_msg)) {
    cerr << "Error detected in mesh data structure." << endl;
    if (error_msg != "") { cerr << error_msg << endl; }
    exit(-1);
  }

  if (!flag_no_warn || flag_fail_on_non_manifold) {
    const bool flag_oriented_manifold = 
      check_oriented_manifold(mesh, flag_no_warn);

    if (flag_fail_on_non_manifold && !flag_oriented_manifold) {
      if (!flag_no_warn) {
        cerr << "Detected non-manifold or inconsistent orientations."
             << endl;
        cerr << "Exiting." << endl;
      }

      exit(-1);
    }

    return flag_oriented_manifold;
  }
  else {
    return true;
  }
}


/// Check whether mesh is an oriented manifold.
/// - Returns true if mesh is an oriented manifold.
/// @param flag_no_warn If true, do not print warning messages.
template <typename MESH_TYPE>
bool check_oriented_manifold
(const MESH_TYPE & mesh, const bool flag_no_warn)
{
  int iv, ihalf_edgeA, ihalf_edgeB;
  bool flag_non_manifold_vertex, flag_non_manifold_edge;

  const bool flag_manifold =
    mesh.CheckManifold
    (iv, ihalf_edgeA, flag_non_manifold_vertex, 
     flag_non_manifold_edge);

  if (flag_non_manifold_edge) {
    if (!flag_no_warn) {
      cerr << "Warning: Non-manifold edge ("
           << mesh.HalfEdge(ihalf_edgeA)->EndpointsStr(",")
           << ")." << endl;
    }

    // Non-manifold edge automatically implies inconsistent orientations.
    return false;
  }

  const bool is_oriented =
    mesh.CheckOrientation(ihalf_edgeB);

  if (is_oriented) {
    if (flag_non_manifold_vertex) {
      if (!flag_no_warn) {
        cerr << "Warning: Non-manifold vertex " << iv << "." << endl;
      }
      return false;
    }
  }
  else {
    if (flag_non_manifold_vertex) {
      if (!flag_no_warn) {
        cerr << "Warning: Non-manifold vertex or inconsistent orientations in cells incident on vertex " << iv << "." << endl;
      }
    }
    else {
      if (!flag_no_warn) {
        cerr << "Warning: Inconsistent orientation of cells incident on edge ("
             << mesh.HalfEdge(ihalf_edgeB)->EndpointsStr(",")
             << ")." << endl;
      }
    }

    return false;
  }

  return true;
 }


 // Print a warning message if collapsing ihalf_edge is illegal or
 //   will change mesh topology.
 // - Return true if collapse is not illegal and does not change
 //   manifold topology.
 bool check_edge_collapse(const HALF_EDGE_MESH_DCMT_A & mesh, 
                          const HALF_EDGE_DCMT_A * half_edge,
                          const bool flag_no_warn)
 {
     const int icell = half_edge->Cell()->Index();
     bool return_flag = true;

     int ivC;
     if (mesh.IsIllegalEdgeCollapse(half_edge)) {
       if (!flag_no_warn) {
         cout << "Collapse of edge (" 
              << half_edge->EndpointsStr(",") << ") is illegal."
              << endl;
         cout << "  Some cell contains vertices "
              << half_edge->EndpointsStr(" and ")
              << " but not edge "
              << half_edge->EndpointsStr(",") << ")." << endl;
       }

       return_flag = false;
     }

     if (mesh.FindTriangleHole(half_edge, ivC)) {
       if (!flag_no_warn) {
         cout << "Collapsing edge ("
              << half_edge->EndpointsStr(",")
              << ") will change the mesh topology." << endl;
         cout << "  Vertices ("
              << half_edge->EndpointsStr(", ")
              << ", " << ivC << ") form a triangle hole." << endl;
       }

       return_flag = false;
     }

     if (!half_edge->IsBoundary()) {
       if (half_edge->FromVertex()->IsBoundary() &&
           half_edge->ToVertex()->IsBoundary()) {
         if (!flag_no_warn) {
           cout << "Collapsing edge ("
                << half_edge->EndpointsStr(",")
                << ") merges two non-adjacent boundary vertices." 
                << endl;
         }

         return_flag = false;
       }
     }

     if (mesh.IsIsolatedTriangle(icell)) {
       if (!flag_no_warn) {
         cout << "Collapsing edge ("
              << half_edge->EndpointsStr(",")
              << ") will delete isolated cell " << icell << endl;
       }

       return_flag = false;
     }

     if (mesh.IsInTetrahedron(icell)) {
       if (!flag_no_warn) {
         cout << "Collapsing edge ("
              << half_edge->EndpointsStr(",")
             << ") will collapse a tetrahedron." << endl;
      }

      return_flag = false;
    }

    return return_flag;
}


// Print a warning message if splitting cell at diagonal
//   (half_edgeA->FromVertex(), half_edgeB->FromVertex())
//   will change mesh topology.
// - Return true if split does not change manifold topology.
bool check_split_cell(const HALF_EDGE_MESH_DCMT_A & mesh, 
                         const HALF_EDGE_DCMT_A * half_edgeA,
                         const HALF_EDGE_DCMT_A * half_edgeB,
                         const bool flag_no_warn)
{
  const VERTEX_DCMT_A * vA = half_edgeA->FromVertex();
  const VERTEX_DCMT_A * vB = half_edgeB->FromVertex();
  const HALF_EDGE_DCMT_A * half_edgeC = mesh.FindEdge(vA, vB);

  bool flag_cell_edge = false;
  bool return_flag = true;

  if (mesh.IsIllegalSplitCell(half_edgeA, half_edgeB)) {
    const VERTEX_DCMT_A * vA = half_edgeA->FromVertex();
    const VERTEX_DCMT_A * vB = half_edgeB->FromVertex();

    if ((vA == half_edgeB->ToVertex()) ||
        (vB == half_edgeA->ToVertex())) 
      { flag_cell_edge = true; }

    if (!flag_no_warn) {
      if (flag_cell_edge) {
        cout << "(" << vA->Index() << "," << vB->Index() 
             << ") is a cell edge, not a cell diagonal."
             << endl;
      }
      else {
        cout << "Illegal split of cell " << half_edgeA->CellIndex()
             << " with diagonal (" << vA->Index() << ","
             << vB->Index() << ")." << endl;
      }
    }
    return_flag = false;
  }

  if (half_edgeC != NULL && !flag_cell_edge) {
    if (!flag_no_warn) {
      cout << "Splitting cell " << half_edgeA->CellIndex()
           << " with diagonal ("
           << vA->Index() << "," << vB->Index() 
           << ") creates an edge incident on three or more cells." << endl;
    }
    return_flag = false;
  }

  return return_flag;
}


// Print a warning message if joining cells separated by ihalf_edge 
//   is illegal.
// - Return true if join is not illegal.
bool check_join_cell(const HALF_EDGE_MESH_DCMT_A & mesh, 
                     const HALF_EDGE_DCMT_A * half_edge,
                     const bool flag_no_warn)
{
  const int TWO(2);
  bool return_flag = true;

  if (mesh.IsIllegalJoinCells(half_edge)) {

    const HALF_EDGE_DCMT_A * half_edgeX = 
      half_edge->NextHalfEdgeAroundEdge();

    if (!flag_no_warn) {
      if (half_edge->IsBoundary()) {
        cout << "Only one cell contains edge ("
             << half_edge->EndpointsStr(",") << ")." << endl;
      }
      else if (!half_edge->FromVertex()->IsIncidentOnMoreThanTwoEdges()) {
        cout << "Half edge endpoint " << half_edge->FromVertexIndex()
             << " is incident on only two edges." << endl;
      }
      else if (!half_edge->ToVertex()->IsIncidentOnMoreThanTwoEdges()) {
        cout << "Half edge endpoint " << half_edge->ToVertexIndex()
             << " is incident on only two edges." << endl;
      }
      else if (half_edge != half_edgeX->NextHalfEdgeAroundEdge()) {
        cout << "More than two cells are incident on edge ("
             << half_edge->EndpointsStr(",") << ")." << endl;
      }
      else {
        const int num_shared_vertices =
          mesh.CountNumVerticesSharedByTwoCells
          (half_edge->Cell(), half_edgeX->Cell());

        if (num_shared_vertices > TWO) {
          cout << "Cells " << half_edge->CellIndex()
               << " and " << half_edgeX->CellIndex()
               << " share " << num_shared_vertices 
               << " vertices." << endl;
        }
        else {
          cout << "Join of two cells incident on edge ("
               << half_edge->EndpointsStr(",") << ") is illegal." << endl;
        }
      }
    }

    return_flag = false;
  }

  return return_flag;
}


/// Reduce checks on large data sets.
template <typename MESH_TYPE>
void reduce_checks_on_large_datasets
(const MESH_TYPE &  mesh, const bool flag_no_warn)
{
  if (flag_reduce_checks) {
    // Already reducing checks.
    return;
  }

  const int num_cells = mesh.CountNumCells();
  if (num_cells >= LARGE_DATA_NUM_CELLS) {
    if (!flag_no_warn) {
      cout << "Warning: Large data set with " << num_cells 
           << " cells." << endl;
      cout << "  Reducing checks (using -flag_reduce_checks.)" << endl;
    }

    flag_reduce_checks = true;
  }
}


// *** Parse/print/prompt functions ***

void parse_command_line(int argc, char ** argv)
{
  using namespace std;

  int iarg = 1;
  while (iarg < argc && argv[iarg][0] == '-') {

    std::string s = argv[iarg];
    if (s == "-collapse_edges")
      { flag_collapse_edges = true; }
    else if (s == "-collapse_short_edges")
      { flag_collapse_short_edges = true; }
    else if (s == "-split_cells")
      { flag_split_cells = true; }
    else if (s == "-split_all_cells")
      { flag_split_all_cells = true; }
    else if (s == "-split_edges")
      { flag_split_edges = true; }
    else if (s == "-split_long_edges")
      { flag_split_long_edges = true; }
    else if (s == "-join_cells")
      { flag_join_cells = true; }
    else if (s == "-join_each_cell")
      { flag_join_each_cell = true; }
    else if (s == "-split_long_edges_cells") {
      flag_split_long_edges = true;
      flag_split_all_cells = true;
    }
    else if (s == "-allow_non_manifold")
      { flag_allow_non_manifold = true; }
    else if (s == "-fail_on_non_manifold")
      { flag_fail_on_non_manifold = true; }
    else if (s == "-reduce_checks")
      { flag_reduce_checks = true; }
    else if (s == "-terse")
      { flag_terse = true; }
    else if (s == "-s") {
      flag_terse = true;
      flag_silent = true; 
    }
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


// Prompt for mesh edge.
// - Return NULL if user enters -1.
const HALF_EDGE_DCMT_A *
prompt_for_mesh_edge(const HALF_EDGE_MESH_DCMT_A & mesh,
                     const bool flag_only_internal)
{
  int iv0, iv1, ihalf_edge;
  std::string error_msg;

  while (true) {
    cout << "Enter vertex (-1 to end): ";
    cin >> iv0;
    if (iv0 < 0) { return NULL; }

    if (!mesh.CheckVertexIndex(iv0, error_msg)) {
      cout << error_msg << endl;
      continue;
    }

    const VERTEX_DCMT_A * v0 = mesh.Vertex(iv0);

    if (v0->NumHalfEdgesFrom() == 0) {
      cout << "Vertex " << iv0 << " is not incident on any cell. " << endl;
      continue;
    }

    if (flag_only_internal) {
      int num_internal_half_edges_from = 0;
      cout << "Internal half edges from " << iv0 << ":";
      for (int k = 0; k < v0->NumHalfEdgesFrom(); k++) {
        const HALF_EDGE_DCMT_A * half_edge = v0->KthHalfEdgeFrom(k);
        if (!half_edge->IsBoundary()) {
          cout << "  (" << half_edge->EndpointsStr(",") << ")";
          num_internal_half_edges_from++;
        }
      }
      cout << endl;

      if (num_internal_half_edges_from == 0) {
        cout << "No internal half edges from " << iv0 << "." << endl;
        cout << "Start again." << endl;
        cout << endl;
        continue;
      }
    }
    else {
      cout << "Half edges from " << iv0 << ":";
      for (int k = 0; k < v0->NumHalfEdgesFrom(); k++) {
        const HALF_EDGE_DCMT_A * half_edge = v0->KthHalfEdgeFrom(k);
        cout << "  (" << half_edge->EndpointsStr(",") << ")";
      }
      cout << endl;
    }

    cout << "Enter vertex adjacent to vertex " << iv0 << " (-1 to end): ";
    cin >> iv1;
    if (iv1 < 0) { return NULL; }

    if (!mesh.CheckVertexIndex(iv1, error_msg)) {
      cout << error_msg << endl;
      continue;
    }

    HALF_EDGE_DCMT_A * half_edge0 = v0->FindIncidentHalfEdge(iv1);

    if (half_edge0 == NULL) {
      cout << "Mesh does not have a half edge ("
           << iv0 << "," << iv1 << ")." << endl << endl;;
      continue;
    }

    if (flag_only_internal && half_edge0->IsBoundary()) {
      cout << "Half edge (" << half_edge0->EndpointsStr(",")
           << ") is a boundary half edge." << endl;
      cout << "Start again." << endl;
      cout << endl;
      continue;
    }

    return half_edge0;
  }
}


void print_cells_with_more_than_three_vertices
(const int max_num, const std::vector<int> & cell_list)
{
  cout << "Cells with more than three vertices";
  if (cell_list.size() == max_num) 
      { cout << " (partial list)"; }
  cout << ": ";
  for (int i = 0; i < cell_list.size(); i++)
    { cout << " " << cell_list[i]; }
  cout << endl;
}


void print_time(const char * label, const double time)
{
  cout << label << time/CLOCKS_PER_SEC << " seconds." << endl;
}


void print_mesh_info(const HALF_EDGE_MESH_DCMT_A & mesh)
{
  const int FIVE(5);
  const int num_vertices = mesh.CountNumVertices();
  const int num_isolated_vertices = 
    mesh.CountNumIsolatedVertices();
  const int num_edges = mesh.CountNumEdges();
  const int num_boundary_edges = mesh.CountNumBoundaryEdges();
  const int num_cells = mesh.CountNumCells();
  const int num_triangles = mesh.CountNumTriangles();
  const int num_quads = mesh.CountNumQuads();
  const int num_large_cells = mesh.CountNumCellsOfSizeGE(FIVE);
  double min_edge_length_squared, max_edge_length_squared;
  double min_edge_length_ratio_squared;
  double cos_min_angle, cos_max_angle;
  int iv, ihalf_edgeA, ihalf_edgeB;
  bool flag_non_manifold_vertex, flag_non_manifold_edge;


  mesh.ComputeMinMaxEdgeLengthSquared
    (min_edge_length_squared, max_edge_length_squared);
  mesh.ComputeMinCellEdgeLengthRatioSquared
    (min_edge_length_ratio_squared);
  mesh.ComputeCosMinMaxAngle(cos_min_angle, cos_max_angle);

  const bool flag_manifold =
    mesh.CheckManifold
    (iv, ihalf_edgeA, flag_non_manifold_vertex, 
     flag_non_manifold_edge);
    
  const bool is_oriented =
    mesh.CheckOrientation(ihalf_edgeB);

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
  cout << "  Number of mesh triangles: " << num_triangles << endl;
  cout << "  Number of mesh quadrilaterals: " << num_quads << endl;
  cout << "  Number of cells with > 4 vertices: " 
       << num_large_cells << endl;
  cout << "Min edge length: " << sqrt(min_edge_length_squared) << endl;
  cout << "Max edge length: " << sqrt(max_edge_length_squared) << endl;
  cout << "Min cell edge length ratio: " 
       << sqrt(min_edge_length_ratio_squared) << endl;
  cout << "Minimum cell angle: "
       << acos(cos_min_angle)*180.0/M_PI << endl;
  cout << "Maximum cell angle: "
       << acos(cos_max_angle)*180.0/M_PI << endl;


  if (flag_manifold && is_oriented) {
    cout << "Mesh is an oriented manifold." << endl;
  }
  else {
    cout << "Mesh is non-manifold or has inconsistent cell orientations." 
         << endl;
  }

}

void usage_msg(std::ostream & out)
{
  out << "Usage: decimate_mesh [OPTIONS] <input filename> [<output filename>]" 
      << endl;
  out << "OPTIONS: " << endl;
  out << "  [-list_no_collapse] [-collapse_edges] [-collapse_short_edges]" 
       << endl;
  out << "  [-split_edges] [-split_long_edges]" << endl;
  out << "  [-split_cells] [-split_all_cells] [-split_all_edges_cells]" 
       << endl;
  out << "  [-join_cells] [-join_each_cell]" << endl;
  out << "  [-allow_non_manifold] [-fail_on_non_manifold]" << endl;
  out << "  [-s | -terse] [-no_warn] [-time] [-h]" << endl;
}


void help_msg()
{
  usage_msg(cout);
  cout << endl;
  cout << "decimate_mesh - Decimate mesh." << endl;
  cout << "   Collapse/split/join mesh edges or cells." << endl;

  cout << endl;
  cout << "Options:" << endl;
  cout << "-collapse_edges:  Prompt and collapse edges." << endl;
  cout << "-collapse_short_edges: Attempt to collapse shortest edge in each cell." << endl;
  cout << "-split_edges:     Prompt and split edges." << endl;  
  cout << "-split_long_edges: Split longest edge in each cell." << endl;
  cout << "-split_cells:     Prompt and split cells across diagonals." << endl;
  cout << "-split_all_cells: Attempt to split all cells." << endl;
  cout << "-split_long_edges_cells: Split long edges in each cell"
       << endl
       << "     and then split all cells."
       << endl;
  cout << "-join_cells:      Prompt and join cells sharing edges." << endl;
  cout << "-join_each_cell:  Attempt to join each cell with adjacent cell"
       << endl
       << "     sharing the largest cell edge." << endl;
  cout << "-allow_non_manifold:  Allow edge collapses or cell splits"
       << endl
       << "     that create non-manifold conditions." << endl;
  cout << "-fail_on_non_manifold: Exit with non-zero return code (fail)"
       << endl
       << "     if non-manifold or inconsistent orientations detected."
       << endl;
  cout << "-terse:   Terse output. Suppress messages output after each" 
       << endl;
  cout << "     collapse/join/split iteration." 
       << endl;
  cout << "   Does not supress warning messages at each iteration."
       << endl;
  cout << "   Does not supress final mesh information."
       << endl;
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
