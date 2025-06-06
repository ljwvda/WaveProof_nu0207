function wout = divergencefree_ext_Wave(win,shape,symmetry)
% returns divergence free "projection" of the input
% works for both intvals and floats

% find nonzero indices
[nx,ny,nz] = countersymtensor(sizeshape(shape),shape);
nx=nx(:,:,:,:,1);
ny=ny(:,:,:,:,2);
nz=nz(:,:,:,:,3);
nxnonzero=(nx~=0);
nynonzero=(ny~=0);
nznonzero=(nz~=0);
sumnonzero=nxnonzero+nynonzero+nznonzero;

% split into components
wfull = symmetrytofulltensor_ext_Wave(win,shape,symmetry);
w1= wfull(:,:,:,:,1);
w2= wfull(:,:,:,:,2);
w3= wfull(:,:,:,:,3);

% divergence up to a factor i
divx=nx.*w1+ny.*w2+nz.*w3; 

% subtracts divergence appropriately from the components
w1(nxnonzero)=w1(nxnonzero)-divx(nxnonzero)./(sumnonzero(nxnonzero).*nx(nxnonzero));
w2(nynonzero)=w2(nynonzero)-divx(nynonzero)./(sumnonzero(nynonzero).*ny(nynonzero));
w3(nznonzero)=w3(nznonzero)-divx(nznonzero)./(sumnonzero(nznonzero).*nz(nznonzero));

wfull(:,:,:,:,1) = w1;
wfull(:,:,:,:,2) = w2;
wfull(:,:,:,:,3) = w3;

% average over group action
[~,wfull] = symmetrizetensor_ext(wfull,symmetry);

wout = fulltosymmetrytensor_ext(wfull,shape,symmetry);

end

