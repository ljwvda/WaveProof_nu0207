function [success,rmin,rmax,bounds] = checkpolytensorsplit_ext_Wave(x0,nu,forcing,ref,solshape,symmetry,Ndagger,Ntilde,eta,etaOmega,setup,parallel)
% This computes the radii polynomial and checks if it is negative somewhere
% If nu is an intval it uses interval arithmetic 

% Adjusted for extra symmetry cases

if ~exist('parallel','var')
    parallel=false;
end

success=false;
rmin=NaN;
rmax=NaN;

x0(end)=real(x0(end)); % make sure c is real-valued

if exist('intval','file') && isintval(nu)
  disp('including interval arithmetic')
  x0=intval(x0);
  forcing=intval(forcing);
  ref=intval(ref);
  eta=intval(eta);
  etaOmega=intval(etaOmega);
else
  disp('no intervals yet')
end

% guarantee symmetries
ref=conjugatesymmetric_ext_Wave(ref,solshape,symmetry);
x0(1:end-1)=divergencefree_ext_Wave(x0(1:end-1),solshape,symmetry);
x0(1:end-1)=conjugatesymmetric_ext_Wave(x0(1:end-1),solshape,symmetry);

if Ntilde<Ndagger
    error('need to choose Ntilde >= Ndagger')
end

c=real(x0(end));
wsol=x0(1:end-1);
w=symmetrytofulltensor_ext_Wave(wsol,solshape,symmetry);

if exist('setup','var') && strcmp(setup,'2D')
    % special for two dimensional case
    Edaggersh='ell2D'; 
    Etildesh='ell2D';
    disp('assuming 2D solution to reduce tensor size')
else
    Edaggersh='ell';
    Etildesh='ell';
    setup='not2D';
end

% Define the shapes for the index sets Edagger and Etilde
N=sizeshape_Wave(solshape);
Edaggershape.type=Edaggersh;
Edaggershape.Nell=Ndagger;
Edaggershape.nu=nu;
Edaggershape.omega=c;
M=sizeshape_Wave(Edaggershape);
Etildeshape.type=Etildesh;
Etildeshape.Nell=Ntilde;
Etildeshape.nu=nu;
Etildeshape.omega=c;
Q=sizeshape_Wave(Etildeshape);
MNQ=max(M+2*N,Q+2*N);  %in fact Q>=M is assumed anyway
disp(['size of big tensor is ',int2str(MNQ)]);

%%%% Start of bounds

[~,J,~,~,Fext,Fphase,tensors]=FDFsymflextensor_ext_Wave([wsol;c],nu,forcing,ref,solshape,symmetry,Edaggershape);

if exist('intval','file') && isintval(nu)
    A=intval(inv(mid(J)));
else
    A=inv(J);
    disp(['norm of A is ',num2str(norm(A,1))]);
end 
disp(['the dimension of A is ',int2str(size(A,1))]);

%cond(A)

%%% Y bound %%%

Y=getYbound_Wave(A,Fext,Fphase,solshape,symmetry,Edaggershape,nu,c,eta,etaOmega);
disp(['Y bound is (wave) ',num2str(altsup(Y))]);

%%% Z2 bound %%%

Z2=getZ2bound_Wave(A,symmetry,Edaggershape,nu,c,eta,etaOmega,Ndagger);
disp(['Z2 bound is (wave) ',num2str(altsup(Z2))]);

%%% Z0 bound %%%

Z0=getZ0bound_Wave(A,J,symmetry,Edaggershape,eta,etaOmega);
clear('J'); % no longer needed
disp(['Z0 bound is (wave) ',num2str(altsup(Z0))]);

%%% Z1 bound %%%

%%% Tail part of Z1 %%%%

Z1tail=getZ1tailbound_Wave(wsol,tensors,solshape,symmetry,nu,eta,Ntilde,setup);
disp(['the Z1 tail term is (wave) ',num2str(altsup(Z1tail))]); 

%%% Finite part of Z1 %%%
%%% This contains the long loop %%%

if parallel
    Z1finite=getZ1finiteboundparallel_Wave(A,w,tensors,...
                solshape,Edaggershape,Etildeshape,symmetry,nu,c,eta,etaOmega,ref);    
else
    Z1finite=getZ1finitebound_Wave(A,w,tensors,...
                solshape,Edaggershape,Etildeshape,symmetry,nu,c,eta,etaOmega,ref);
end
disp(['the finite term of Z1 is (wave) ',num2str(altsup(Z1finite))]);

Z1=max([Z1finite,Z1tail]);
disp(['Z1 bound is (wave) ',num2str(altsup(Z1))]);

%%%%%%%%%%% Conclusion %%%%%%%%%%%%

discr=(1-Z1-Z0)^2-2*Y*Z2;
disp(['discriminant is ',num2str(altsup(discr))]);
if 1-Z1-Z0>0 && discr>0
   disp('SUCCESS')
   if exist('intval','file') && isintval(nu)
      disp('including full interval arithmetic')
   else
      disp('but no interval arithmetic')
   end
   rmin=altsup((1-Z1-Z0-sqrt(discr))/Z2);
   rmax=-altsup(-(1-Z1-Z0)/Z2);
   disp(['r_min is ',num2str(rmin)]);
   disp(['r_max is ',num2str(rmax)]);
   success=true;
else
   disp('FAILURE')
end

% store the bounds which have been obtained
bounds.Y=altsup(Y);
bounds.Z0=altsup(Z0);
bounds.Z1=altsup(Z1);
bounds.Z2=altsup(Z2);
bounds.Z1finite=altsup(Z1finite);
bounds.Z1tail=altsup(Z1tail);

end





