function Z1finite=getZ1finitebound_Wave(A,w,tensors,solshape,Edaggershape,Etildeshape,symmetry,nu,c,eta,etaOmega,ref)
% computes the finite part of the Z1 bound for traveling waves

Dw=tensors.Dw;
Mw=tensors.Mw;
DMw=tensors.DMw;

% biggest index set needed
% sum of 2*solshape and Etildeshape
N=sizeshape_Wave(solshape);
Q2Nshape.type='otherplusrec';
Q2Nshape.other=Etildeshape;
Q2Nshape.Nrec=2*N;
Q2N=sizeshape_Wave(Q2Nshape);

[nx,ny,nz,nt,ninshape,comp] = countersymtensor(Q2N,Q2Nshape);
[symvar,symindex,symfactor,grouporder,multiplicity,symindexextra] = symmetryindicestensor_ext(nx,ny,nz,nt,comp,Q2N,ninshape,symmetry);
[tilden2,tilden2reci] = tildentensor(nx,ny,nz,nu);

% indices for the solution
nsol=shapetensor(nx,ny,nz,nt,solshape);
solsymvar=(symvar & nsol);
% indices for the finite part
nfinite=shapetensor(nx,ny,nz,nt,Edaggershape);
finitesymvar=(symvar & nfinite);
finitesymindex=symindex(finitesymvar);
% indices for the tail part
tailsymvar=(symvar & ~nfinite);
tailsymindex=symindex(tailsymvar);
% indices for the solution variables outside the finite (jacobian) part
soltailsymvar=(solsymvar & ~finitesymvar);

% the set Stilde= Ssol+Etilde
% assuming Ssol is rectangular
% (if not, then it is still rigorous, just computes too much)
QNshape.type='otherplusrec';
QNshape.other=Etildeshape;
QNshape.Nrec=N;
% the set Ssol+Edagger
MNshape.type='otherplusrec';
MNshape.other=Edaggershape;
MNshape.Nrec=N;

% indices for the set Ssol+Edagger
nfiniteplussol=shapetensor(nx,ny,nz,nt,MNshape);
% symmetry reduced variables in Sol+Edagger \ Edagger
finiteplussoltailsymvar=(symvar & nfiniteplussol & ~finitesymvar);
% indices for the set Stilde
ntildeS=shapetensor(nx,ny,nz,nt,QNshape);
tildeSsymvar=(symvar & ntildeS);

wfull=setsizetensor(w,Q2N);
reffull=setsizetensor(symmetrytofulltensor_ext(ref,solshape,symmetry),Q2N);

% weights in the symmetrized norm
orbits=grouporder./multiplicity; % orbit-stabilizer formula
weights=eta.^(abs(nx)+abs(ny)+abs(nz)).*orbits;
weightsetaOmega=[weights(finitesymvar);etaOmega];
weightsymvar=weights(symvar);

% diagonal part
lambda=reshape(1./abs(-1i*c*nz(:)+nu*tilden2(:)),size(tilden2)); 
% multiplied by weights
lambdaweights=lambda.*weights;
lambdaweightstail=lambdaweights(tailsymvar);
lambdaweightssoltail=lambdaweights(soltailsymvar);

% the variables in the finite part of the derivative
finitepartofZ1=find(tildeSsymvar);
nsteps=length(finitepartofZ1);
disp(['number of elements in the finite part of Z1 is ',int2str(nsteps)]);

% number of nontrivial variables in the finite part of the derivative
NQ=length(find(symvar));
zerosQ=altzeros([NQ,1],nu);
Qcolumnnorms=altzeros([1,nsteps],nu);

% initialize indices
[nsymvar,NN]=makeindices(N,nx,ny,nz,nt,symvar);
compsymvar=comp(symvar);
compind{1}=(compsymvar==1);
compind{2}=(compsymvar==2);
compind{3}=(compsymvar==3);

% operator norm matrix A -- This is added so we do not have to do the matrix
% multiplication in each loop
v=weightsetaOmega';
Anorm=max((v*abs(A))./v);

for counter=1:nsteps
    % This loop can be parallelized
    j=finitepartofZ1(counter);
    Q0=zerosQ;

    if mod(counter,100)==0
        % just progress report
        % doesn't make much sense when parallelized
        if mod(counter,1000)==0
            disp(['step ',int2str(counter),' out of ',int2str(nsteps),' and so far Z1finite is ', num2str(max(mid(Qcolumnnorms)))]);
        else
            disp(['step ',int2str(counter),' out of ',int2str(nsteps)]);
        end
    end

    % find the variables in the orbit of j
    jjjindex=unique(symindexextra(symindex(j),:));

    for jjj=jjjindex
        % for loop over the group orbit
        nn=[nx(jjj);ny(jjj);nz(jjj)];

        % identify the nonzero components of M_n^{l,m}
        % needed for the terms DPsi2 and DPsi3
        nnreci=1i*nn*tilden2reci(jjj);
        matrix=[0 -3 2; 3 0 -1; -2 1 0];
        column=matrix(:,comp(jjj)); %nonvanishing components l in M_n^{l,m}
        lvals=find(column~=0);
        nvals=column(lvals(:));
        Mvals=sign(nvals).*nnreci(abs(nvals));

        % shift of the indices (coming from derivative of convolution)
        nsx=nsymvar(:,1)-nx(jjj);
        nsy=nsymvar(:,2)-ny(jjj);
        nsz=nsymvar(:,3)-nz(jjj);
        nst=nsymvar(:,4)-nt(jjj);

        % only the shifted indices within the numerical solution rectangle
        % will lead to a nonvanishing contribution
        jj=(abs(nsx)<=N(1) & abs(nsy)<=N(2) & abs(nsz)<=N(3) & abs(nst)<=N(4));
        % Compute Mnlm only for relevant indices
        Mnlm_jj = Mvals(1) * (compind{lvals(1)}(jj)) + Mvals(2) * (compind{lvals(2)}(jj));

        % convert tensor indices into linear indices for nx,ny,nz dimensions of tensor
        Sall1=1+nsx(jj)+N(1)+NN(1)*(nsy(jj)+N(2))+NN(2)*(nsz(jj)+N(3))+NN(3)*(nst(jj)+N(4));
        % next, including the nt dimension
        Sall2=Sall1+NN(4)*(compsymvar(jj)-1);
        % including the component dimension, but selecting the component corresponding to jjj
        % as this is the one needed below in DPsi3
        Sall3=Sall2+NN(5)*(comp(jjj)-1);

        % Compute the four terms in DPsi

        % term sum_{p=1}^3 n_p(Mw)^p_{k-n}
        % note that sum_{p=1}^3 n_p(Mw)^p_{k-n} = sum_{p=1}^3 k_p(Mw)^p_{k-n}
        % since Mw is divergence free for any w
        DPsi1=Mw(Sall1)*nx(jjj)+Mw(Sall1+NN(4))*ny(jjj)+Mw(Sall1+2*NN(4))*nz(jjj);
        % the factor M_n^{l,m}
        DPsi1=DPsi1.*(compsymvar(jj)==comp(jjj));

        % term (D_m(Mw)^l)_{k-n}
        DPsi2=DMw(Sall3);

        % term sum_{p=1}^3 M_n^{p,m} (D_p w^l)_{k-n}
        DPsi3=Dw(Sall2+NN(5)*(lvals(1)-1))*Mvals(1)+Dw(Sall2+NN(5)*(lvals(2)-1))*Mvals(2);

        % term sum_{p=1}^3 n_p w^p_{k-n}
        DPsi4=w(Sall1)*nx(jjj)+w(Sall1+NN(4))*ny(jjj)+w(Sall1+2*NN(4))*nz(jjj);
        % the factor M_n^{l,m}
        %DPsi4=DPsi4.*Mnlm(jj);
        DPsi4=DPsi4.*Mnlm_jj;

        % add to the column with factor symfactor(j)=\tilde{\alpha}(jj,j)
        % no need for factor i
        Q0(jj)=Q0(jj)+symfactor(jjj)*(DPsi1-DPsi2+DPsi3-DPsi4);
    end
    if finiteplussoltailsymvar(j)
        % compute finite part
        if soltailsymvar(j)
            % add derivative phase condition instead of 0
            Q0finite=A*[Q0(finitesymindex);reffull(j)'*nz(j)];
        else
            Q0finite=A*[Q0(finitesymindex);0];
            if norm(Anorm*[Q0(finitesymindex);0])>0.03
                Q0finite=A*[Q0(finitesymindex);0];
            end
        end
        % norm of finite part
        Q0finite=sum(abs(Q0finite).*weightsetaOmega);
    else
        % finite part cancels (absorbed in Z0)
        % or we are far into the tail, where the finite part of A
        % does not hit the nonvanishing part of Q0
        Q0finite=0;
    end
    Q0tail=Q0(tailsymindex);
    Q0tail=sum(abs(Q0tail).*lambdaweightstail);
    % rescale by weight of j-the variable
    Qtotal=(Q0finite+Q0tail)/weightsymvar(symindex(j));
    Qcolumnnorms(counter)=Qtotal;
end

% the norm of the matrix
Z1finite=max(Qcolumnnorms);

% derivative w.r.t. the frequency
Z1freqterm=abs(nz(soltailsymvar).*wfull(soltailsymvar)); 
Z1freqterm=Z1freqterm.*lambdaweightssoltail;
Z1freqterm=sum(Z1freqterm)/etaOmega;
disp(['contribution of ',int2str(length(find(soltailsymvar))),' elements in frequency term to Z1 is ',num2str(altsup(Z1freqterm))]);

Z1finite=max([Z1finite,Z1freqterm]);

end