# Makefile for test_half_edge_mesh

TEMPLATE_FILESA = half_edge_mesh.hpp half_edge_mesh_IO.hpp
TEMPLATE_FILESB = $(TEMPLATE_FILESA) half_edge_mesh_DCMT.hpp half_edge_mesh_coord.hpp

all: test_half_edge_mesh test_half_edge_meshB\
			meshinfo decimate_mesh


test_half_edge_mesh: $(TEMPLATE_FILESA) test_half_edge_mesh.cpp
	g++ -o $@ test_half_edge_mesh.cpp	

test_half_edge_meshB: $(TEMPLATE_FILESA) test_half_edge_meshB.cpp
	g++ -o $@ test_half_edge_meshB.cpp	

meshinfo: $(TEMPLATE_FILESB) meshinfo.cpp
	g++ -o $@ meshinfo.cpp

decimate_mesh: $(TEMPLATE_FILESB) decimate_mesh.cpp
	g++ -o $@ decimate_mesh.cpp

doc: $(TEMPLATE_FILESB) half_edge_mesh_doxygen.config
	doxygen half_edge_mesh_doxygen.config
