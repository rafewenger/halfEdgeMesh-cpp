/// \file half_edge_mesh_IO.hpp
/// template classes for 2D half edge mesh read/write.
/// - Requires C++-11 or later.
/// - Version 0.0.1

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

#include<exception>
#include<fstream>
#include<iostream>
#include<string>
#include<sstream>
#include<vector>

#include "half_edge_mesh.hpp"


#ifndef _HALF_EDGE_MESH_IO_
#define _HALF_EDGE_MESH_IO_

namespace HMESH {

  // *****************************************************************
  // Read .off file
  // *****************************************************************

  /// Get first line that is not blank or a comment.
  /// - First non-blank character in a comment line is '#'.
  template <typename ISTREAM_TYPE>
  bool get_non_comment_line(ISTREAM_TYPE & in, std::string & line)
  {
    while (getline(in, line)) {

      const size_t i = line.find_first_not_of(' ');
      if (i == std::string::npos) {
        // Empty (blank) line.
        continue;
      }

      if (line[i] == '#') {
        // Comment line. Ignore.
        continue;
      }

      // Found non empty, non comment line.
      return(true);
    }

    return(false);
  }


  /// Read off file into MTYPE.
  /// MTYPE is assumed to be a class derived from HALF_EDGE_MESH_BASE.
  /// @pre Dimension of vertices is 3.
  template <typename ISTREAM_TYPE, typename MTYPE>
  void read_off_file(ISTREAM_TYPE & in, MTYPE & mesh)
  {
    const int DIM3(3);
    std::string line;
    std::stringstream str_stream;
    int numv;
    int numpoly;
    int numv_in_poly;
    typename MTYPE::COORD_TYPE coord[DIM3];

    getline(in, line);
    if (line != "OFF") {
      throw SIMPLE_EXCEPTION
        ("Read error. File does not begin with OFF.");
    }

    if (!get_non_comment_line(in, line)) {
      throw SIMPLE_EXCEPTION
        ("Read error. File does not contain line with number of vertices and polygons.");
    }

    str_stream.str(line);
    str_stream >> numv;
    str_stream >> numpoly;
    if (str_stream.bad() || str_stream.fail()) {
      throw SIMPLE_EXCEPTION
        ("Read error. Error reading number of vertices and number of polygons.");
    }

    mesh.AddVertices(numv);

    for (int iv = 0; iv < numv; iv++) {
      if (!get_non_comment_line(in, line)) {
        throw SIMPLE_EXCEPTION
          ("Read error. File is missing vertex coordinates.");
      }

      std::stringstream str_stream(line);

      str_stream >> coord[0];
      str_stream >> coord[1];
      str_stream >> coord[2];

      if (str_stream.bad() || str_stream.fail()) {
        throw SIMPLE_EXCEPTION
          ("Read error. Error reading vertex coordinates.");
      }

      mesh.SetCoord(iv, coord);
    }

    for (int ipoly = 0; ipoly < numpoly; ipoly++) {

      std::vector<int> cell_vlist;

      if (!get_non_comment_line(in, line)) {
        throw SIMPLE_EXCEPTION
          ("Read error. File is missing polygon vertex.");
      }

      std::stringstream str_stream(line);
      str_stream >> numv_in_poly;
      if (str_stream.bad() || str_stream.fail()) {
        throw SIMPLE_EXCEPTION
          ("Read error. Error reading polygon vertices.");
      }

      if (numv_in_poly < 3) {
        // Ignore polygons with less than 3 vertices.
        continue;
      }

      cell_vlist.resize(numv_in_poly);
      for (int k = 0; k < numv_in_poly; k++) {
        int iv;
        str_stream >> cell_vlist[k];
      }

      if (str_stream.bad() || str_stream.fail()) {
        throw SIMPLE_EXCEPTION
          ("Read error. Error reading polygon vertices.");
      }

      mesh.AddCell(cell_vlist);
    }

  }


  /// Open and read off file into MTYPE.
  /// - Version that opens file.
  /// MTYPE is assumed to be a class derived from HALF_EDGE_MESH_BASE.
  /// @pre Dimension of vertices is 3.
  template <typename MTYPE>
  void open_and_read_off_file(const char * input_filename, MTYPE & mesh)
  {
    using namespace std;

    ifstream in(input_filename, ios::in);
    if (!in.good()) {
      cerr << "Unable to open input file " << input_filename << "." << endl;
      exit(-1);
    }

    try {
      read_off_file(in, mesh);
    }
    catch (SIMPLE_EXCEPTION e) {
      cerr << e.what() << endl;
      cerr << "Exiting." << endl;
      in.close();
      exit(-1);
    }

    in.close();
  }


  // *****************************************************************
  // Write .off file
  // *****************************************************************

  /// Write mesh MTYPE into off file.
  /// MTYPE is assumed to be a class derived from HALF_EDGE_MESH_BASE.
  /// @pre Dimension of vertices is 3.
  template <typename OSTREAM_TYPE, typename MTYPE> 
  void write_off_file(OSTREAM_TYPE & out, MTYPE & mesh)
  {
    typedef typename MTYPE::VERTEX_TYPE VERTEX_TYPE;
    typedef typename MTYPE::HALF_EDGE_TYPE HALF_EDGE_TYPE;
    typedef typename MTYPE::CELL_TYPE CELL_TYPE;

    // Print OFF label.
    out << "OFF\n";

    // Print number of vertices and polygons.
    // CORRECTION: 11-21-2021 - RW
    // INCORRECT: '<< mesh.CellListLength() << " "'
    out << mesh.VertexListLength() << " "
        << mesh.CountNumCells() << " "
        << 0 << "\n";
    out << "\n";

    // Print vertex coordinates.
    for (int iv = 0; iv < mesh.VertexListLength(); iv++) {
      const VERTEX_TYPE * v = mesh.Vertex(iv);
      if (v == NULL) 
        { out << "0 0 0" << "\n"; }
      else {
        v->PrintCoord(out, " "); 
        out << "\n";
      }
    }
    out << "\n";

    // Print polygon vertices.
    for (int ipoly = 0; ipoly < mesh.CellListLength(); ipoly++) {
      const CELL_TYPE * cell = mesh.Cell(ipoly);

      if (cell == NULL) { continue; }

      out << cell->NumVertices() << " ";
      const HALF_EDGE_TYPE * half_edge = cell->HalfEdge();
      for (int j = 0; j < cell->NumVertices(); j++) {
        const VERTEX_TYPE * v = half_edge->FromVertex();
        const int iv = v->Index();
        out << " " << iv;
        half_edge = half_edge->NextHalfEdgeInCell();
      }
      out << "\n";
    }
  }

  /// Open and write mesh MTYPE into off file.
  /// MTYPE is assumed to be a class derived from HALF_EDGE_MESH_BASE.
  /// @pre Dimension of vertices is 3.
  template <typename MTYPE> 
  void open_and_write_off_file
  (const char * output_filename, MTYPE & mesh)
  {
    using namespace std;

    ofstream out(output_filename, ios::out);
    if (!out.good()) {
      cerr << "Unable to open output file " << output_filename << "." << endl;
      exit(-1);
    }

    try {
      write_off_file(out, mesh);
    }
    catch (SIMPLE_EXCEPTION e) {
      cerr << e.what() << endl;
      cerr << "Exiting." << endl;
      out.close();
      exit(-1);
    }

    out.close();
  }

}

#endif
