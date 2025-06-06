function [nx,ny,nz,nt,ninshape,comp] = countersymtensor(N,shape)
% produces commonly used counters and index variables

[nx,ny,nz,nt,comp]=ndgrid(-N(1):N(1),-N(2):N(2),-N(3):N(3),-N(4):N(4),1:3);

% select the shape
ninshape = shapetensor(nx,ny,nz,nt,shape);



