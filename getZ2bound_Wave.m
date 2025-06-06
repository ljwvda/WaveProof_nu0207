function Z2=getZ2bound_Wave(A,symmetry,Edaggershape,nu,c,eta,etaOmega,Ndagger)
% computes the Z2 bound for traveling wave

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
% weights=eta.^(abs(nx)+abs(ny)+abs(nz)+abs(nt)).*orbits;
weights=eta.^(abs(nx)+abs(ny)+abs(nz)).*orbits;
weights=weights(symvar);
finiteweightsetetaOmega=[weights(symindexfinite);etaOmega];

% operator norm in the rescaled symmetrized norm
% excluding the phase equation
% but including the frequency variable
v=finiteweightsetetaOmega'; 
vv=v(1:end-1);
%maxn=max(max(max(max(1,abs(nx)),abs(ny)),abs(nz)),abs(nt));
maxn=max(max(max(1,abs(nx)),abs(ny)),abs(nz));
vvmaxn=vv./maxn(finitesymvar)';
normAmaxn=max(v*abs(A(1:end,1:end-1))./vvmaxn);

if exist('intval','file') && isintval(nu)
    roottwo=sqrt(intval(2));
else
    roottwo=sqrt(2);
end

% including the formula for the tail estimate
% Z2=(4+roottwo)*max([normAmaxn,1/c,1/(sqrt(nu*Ndagger))]);  
Z2=(4+roottwo)*max([normAmaxn,1/(sqrt(nu*Ndagger))]);  





