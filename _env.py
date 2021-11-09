# -*- coding: utf-8 -*-


NODE_INTERIOR = -1  # nodes located in the interior
# nodes located in the complement of (interior + frontier)
NODE_COMPLEMENTARY = -2
NODE_DIRICHLET = 1  # nodes with dirichlet boundary condition
NODE_NEUMANN = 2  # nodes with neumann boundary condition
NODE_ROBIN = 3  # nodes with robin boundary condition
NODE_LINER_A = 4  # nodes with robin boundary condition along liner_a
NODE_LINER_B = 5  # nodes with robin boundary condition along liner_b
