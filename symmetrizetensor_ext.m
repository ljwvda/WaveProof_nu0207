function [wout,average] = symmetrizetensor_ext(win,symmetry)
% uses group actions to "symmetrize" the tensor
% the second output normalizes with the order of the group (averaging)

% Extended with cases

sz=(size(win));
N=(sz(1:4)-1)/2;
shape.type='rec';
shape.Nrec=N;
[nx,ny,nz,nt] = countersymtensor(N,shape);  % needed for some symmetries
win=+win; % converts logicals to integers if needed

switch symmetry
    case 0 % no symmetry
        wout=win;
        average=wout;
    case 1 % Sz symmetry
        wtmp=win;
        wtmp2=wtmp(:,:,end:-1:1,:,:);   % Sz symmetry
        wtmp2(:,:,:,:,[1 2])=-wtmp2(:,:,:,:,[1 2]);
        wout=wtmp+wtmp2;
        average=wout/2;
    case 2 % Sz and DSx and P2D symmetry
        wout=win;
        wtmp1=wout(:,:,end:-1:1,:,:);   % Sz symmetry
        wtmp1(:,:,:,:,[1 2])=-wtmp1(:,:,:,:,[1 2]);
        wout=wout+wtmp1;
        wtmp2=wout(end:-1:1,:,:,:,:);   % DSx symmetry
        wtmp2(:,:,:,:,[2 3])=-wtmp2(:,:,:,:,[2 3]);
        wtmp2=(-1).^(nx+ny).*wtmp2;
        wout=wout+wtmp2;
        wtmp3=(-1).^(nx+ny+nt).*wout;  % P2D symmetry
        wout=wout+wtmp3;
        average=wout/8;
    case 4  % 2D flow and SxSy and DSx and P4SxR and P2D symmetry
        wout=win;
        wtmp1=wout(end:-1:1,end:-1:1,:,:,:);   % SxSy symmetry
        wtmp1(:,:,:,:,[1 2])=-wtmp1(:,:,:,:,[1 2]);
        wout=wout+wtmp1;
        wtmp2=wout(end:-1:1,:,:,:,:);   % DSx symmetry
        wtmp2(:,:,:,:,[2 3])=-wtmp2(:,:,:,:,[2 3]);
        wtmp2=(-1).^(nx+ny).*wtmp2;
        wout=wout+wtmp2;
        wtmp3=permute(wout,[2 1 3 4 5]);  % P4SxR symmetry
        wtmp3(:,:,:,:,[1 2])=wtmp3(:,:,:,:,[2 1]);
        wtmp3=-(-1).^(ny).*(1i).^(nt).*wtmp3;
        wout=wout+wtmp3;
        wtmp4=(-1).^(nx+ny+nt).*wout;  % P2D symmetry
        wout=wout+wtmp4;
        average=wout/16;
   
    case 20 % Third Hopf bifurcation, specific for traveling wave, tested for nu=0.21
        wout=win;
        % For comp3 we can mirror in x and y direction (simultaneously), see case 4
        % For comp1 and comp2 the same holds but we also need to multiply by -1
        wtmp1=wout(end:-1:1,end:-1:1,:,:,:);
        wtmp1(:,:,:,:,[1 2])=-wtmp1(:,:,:,:,[1 2]);
        wout=wout+wtmp1;

        %%
        % % Mirror in x direction (and multiply by -1):
        wtmp2=wout(end:-1:1,:,:,:,:);   % DSx symmetry, see case 4
        wtmp2(:,:,:,:,[2 3])=-wtmp2(:,:,:,:,[2 3]);
        wtmp2=(-1).^(nx+ny).*wtmp2;
        wout=wout+wtmp2;

        % We can swap x and y (and multiply by (-)i) for third
        % component
        % 2nd comp is 1st comp "transposed" (up to multiplication by (-)i)
        wtmp3=permute(wout,[2 1 3 4 5]);  % P4SxR symmetry
        wtmp3(:,:,:,:,[1 2])=wtmp3(:,:,:,:,[2 1]);
        wtmp3=-(-1).^(ny).*(1i).^(nz).*wtmp3;
        wout=wout+wtmp3;

        % We have zeros on diagonal and even off-diagonals:
        wtmp4=(-1).^(nx+ny+nz).*wout;  % P2D symmetry, see case 4
        wout=wout+wtmp4;

        average = wout/16;
    case 21 % Hopf bifurcation problem (symmetry for eigenvectors) --> main_v1_proof.m
        wout=win;
        wtmp1=wout(end:-1:1,end:-1:1,:,:,:);   % SxSy symmetry
        wtmp1(:,:,:,:,[1 2])=-wtmp1(:,:,:,:,[1 2]);
        wout=wout+wtmp1;
        wtmp2=wout(end:-1:1,:,:,:,:);   % DSx symmetry
        wtmp2(:,:,:,:,[2 3])=-wtmp2(:,:,:,:,[2 3]);
        wtmp2=(-1).^(nx+ny).*wtmp2;
        wout=wout+wtmp2;
        % wtmp3=permute(wout,[2 1 3 4 5]);  % P4SxR symmetry
        % wtmp3(:,:,:,:,[1 2])=wtmp3(:,:,:,:,[2 1]);
        % wtmp3=-(-1).^(ny).*(1i).^(nt).*wtmp3;
        % wout=wout+wtmp3;
        % wtmp4=(-1).^(nx+ny+nt).*wout;  % P2D symmetry
        % wout=wout+wtmp4;
        average=wout/4;
    case 22 % Hopf bifurcation problem (symmetry for eigenvectors) --> main_v1_proof.m
        wout=win;
        wtmp1=wout(end:-1:1,end:-1:1,:,:,:);   % SxSy symmetry
        wtmp1(:,:,:,:,[1 2])=-wtmp1(:,:,:,:,[1 2]);
        wout=wout+wtmp1;
        % wtmp2=wout(end:-1:1,:,:,:,:);   % DSx symmetry
        % wtmp2(:,:,:,:,[2 3])=-wtmp2(:,:,:,:,[2 3]);
        % wtmp2=(-1).^(nx+ny).*wtmp2;
        % wout=wout+wtmp2;
        % wtmp3=permute(wout,[2 1 3 4 5]);  % P4SxR symmetry
        % wtmp3(:,:,:,:,[1 2])=wtmp3(:,:,:,:,[2 1]);
        % wtmp3=-(-1).^(ny).*(1i).^(nt).*wtmp3;
        % wout=wout+wtmp3;
        % wtmp4=(-1).^(nx+ny+nt).*wout;  % P2D symmetry
        % wout=wout+wtmp4;
        average=wout/2;
    case 23 % Hopf bifurcation problem (symmetry for eigenvectors) --> main_v1_proof.m
        wout=win;
        % wtmp1=wout(end:-1:1,end:-1:1,:,:,:);   % SxSy symmetry
        % wtmp1(:,:,:,:,[1 2])=-wtmp1(:,:,:,:,[1 2]);
        % wout=wout+wtmp1;
        wtmp2=wout(end:-1:1,:,:,:,:);   % DSx symmetry
        wtmp2(:,:,:,:,[2 3])=-wtmp2(:,:,:,:,[2 3]);
        wtmp2=(-1).^(nx+ny).*wtmp2;
        wout=wout+wtmp2;
        % wtmp3=permute(wout,[2 1 3 4 5]);  % P4SxR symmetry
        % wtmp3(:,:,:,:,[1 2])=wtmp3(:,:,:,:,[2 1]);
        % wtmp3=-(-1).^(ny).*(1i).^(nt).*wtmp3;
        % wout=wout+wtmp3;
        % wtmp4=(-1).^(nx+ny+nt).*wout;  % P2D symmetry
        % wout=wout+wtmp4;
        average=wout/2;
end

end