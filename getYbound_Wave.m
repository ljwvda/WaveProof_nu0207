function Y=getYbound_Wave(A,Fext,Fphase,solshape,symmetry,Edaggershape,nu,c,eta,etac)
% computes the Y0 bound for traveling waves

% biggest index set needed
M2N=max(sizeshape_Wave(Edaggershape),2*sizeshape_Wave(solshape));
M2Nshape.type='rec';
M2Nshape.Nrec=M2N;
[nx,ny,nz,nt,ninshape,comp] = countersymtensor(M2N,M2Nshape);
[symvar,symindex,~,grouporder,multiplicity] = symmetryindicestensor_ext(nx,ny,nz,nt,comp,M2N,ninshape,symmetry);
tilden2 = tildentensor(nx,ny,nz,nu);

% indices for the finite and tail parts
nfinite=shapetensor(nx,ny,nz,nt,Edaggershape);
finitesymvar=(symvar & nfinite);
symindexjac=symindex(finitesymvar);
tailsymvar=(symvar & ~nfinite);
symindextail=symindex(tailsymvar);

% weights in the symmetrized norm
orbits=grouporder./multiplicity; % orbit-stabilizer formula
weights=eta.^(abs(nx)+abs(ny)+abs(nz)).*orbits;
weights=weights(symvar);
finiteweightsetaOmega=[weights(symindexjac);etac];

% the diagonal
lambda=reshape(1./abs(-c*1i*nz(:)+nu*tilden2(:)),size(tilden2));
lambda=lambda(symvar);

Fext=setsizetensor(Fext,M2N);
Fext=Fext(symvar);

% finite part
finiteresidue=A*[Fext(symindexjac);Fphase];
Y1=sum(abs(finiteresidue).*finiteweightsetaOmega);

% tail part
tailresidue=Fext(symindextail).*lambda(symindextail);
Y2=sum(abs(tailresidue).*weights(symindextail));

Y=Y1+Y2;

end