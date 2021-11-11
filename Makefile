# Makefile for test_half_edge_mesh

TEMPLATE_FILES = half_edge_mesh.hpp half_edge_mesh_IO.hpp

test_half_edge_mesh: $(TEMPLATE_FILES) test_half_edge_mesh.cpp
	g++ -o $@ test_half_edge_mesh.cpp	

doc: $(TEMPLATE_FILES) half_edge_mesh_doxygen.config
	doxygen half_edge_mesh_doxygen.config
