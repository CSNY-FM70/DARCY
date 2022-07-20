clear; clc;
%Discretization Degree
ns=20;
%Global indexing
[xv,yv,elt2vert,nvtx,ne,h]=uniform_mesh_info(ns);
[uh,A_int,rhs]=twod_linear_FEM(ns,xv,yv,elt2vert,nvtx,ne,h,exp(ones(ne,1)),zeros(ne,1));

fileID = fopen('A_K20.dat','w');
fprintf(fileID,'%8.5f\n',A_int);
type('A_K20.dat')


fileID = fopen('b_K20.dat','w');
fprintf(fileID,'%8.5f\n',rhs);
type('b_K20.dat')


fileID = fopen('sol_K20.dat','w');
fprintf(fileID,'%8.5f\n',uh);
type('sol_K20.dat')

