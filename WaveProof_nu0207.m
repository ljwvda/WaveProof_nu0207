close all

% Load the data
data = load('WaveProof_nu0207.mat');
phi0=data.x0;
nu=data.nu;
solshape=data.solshape;
symmetry=data.symmetry;

forcing=classicforcingtensor;

nu=intval(nu); % Comment this if you want to run the proof without
%interval arithmetic

tic

%% Proof
parallel=false; % to force serial execution of the loop 
%parallel=true; % uses parfor (also runs without the parallel toolbox)

ref = phi0(1:end-1); % Exclude wave speed

%  Ntilde >= Ndagger (to ensure finite part A does not influence tail
%  estimate)
Ndagger = 300; % Size of the computational/finite part A
Ntilde = 700; % To balance the quality of the estimates and the computational costs

eta=1+1e-7;
etac=1+1e-7;
setup='not2D';

% run the proof
[success,rmin,rmax,bounds]=checkpolytensorsplit_ext_Wave(phi0,nu,forcing,ref,...
          solshape,symmetry,Ndagger,Ntilde,eta,etac,setup,parallel); 

toc