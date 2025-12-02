function Z0=getZ0bound_Wave(A,J,symmetry,Edaggershape,eta,etaOmega)
% computes the Z0 bound for traveling waves

% index set needed
M=sizeshape_Wave(Edaggershape);
[nx,ny,nz,nt,ninshape,comp] = countersymtensor(M,Edaggershape);
[symvar,symindex,~,grouporder,multiplicity] = symmetryindicestensor_ext(nx,ny,nz,nt,comp,M,ninshape,symmetry);

% indices for the finite part
nfinite=shapetensor(nx,ny,nz,nt,Edaggershape);
finitesymvar=(symvar & nfinite);
symindexfinite=symindex(finitesymvar);

% weights in the symmetrized norm
orbits=grouporder./multiplicity; % orbit-stabilizer formula
weights=eta.^(abs(nx)+abs(ny)+abs(nz)).*orbits;
weights=weights(symvar);
finiteweightsetaOmega=[weights(symindexfinite);etaOmega];

O=abs(eye(size(A))-A*J);

% operator norm
v=finiteweightsetaOmega';
Z0=max((v*O)./v);

end
