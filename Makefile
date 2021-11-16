# Makefile for test_half_edge_mesh

TEMPLATE_FILES = half_edge_mesh.hpp half_edge_mesh_IO.hpp

all: test_half_edge_mesh test_half_edge_meshB


test_half_edge_mesh: $(TEMPLATE_FILES) test_half_edge_mesh.cpp
	g++ -o $@ test_half_edge_mesh.cpp	

test_half_edge_meshB: $(TEMPLATE_FILES) test_half_edge_meshB.cpp
	g++ -o $@ test_half_edge_meshB.cpp	

doc: $(TEMPLATE_FILES) half_edge_mesh_doxygen.config
	doxygen half_edge_mesh_doxygen.config
