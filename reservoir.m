function [para] = reservoir
 
para.NX = 54;     % number of blocks
para.NY = 22;
para.N = para.NX *para.NY;
para.h = 100; %ft
para.W = 5753; %ft
para.L = 7060.5; %ft
para.ct = 1e-6; %psi^-1
para.mu = 1; %cP
para.Bw = 1; %RB/STB
para.dx = para.L/para.NX; %ft
para.dy = para.W/para.NY; %ft
para.Ax = para.dy*para.h; %ft^2
para.Ay = para.dx*para.h; %ft^2

global k phi
% A = importdata(''); just another way to do this I guess
% phi = A(:,1)';
% k = A(:,2)';
k = load('Nechelik_perm.dat')
phi = load('Nechelik_poro.dat')
para.k = k;
para.phi = phi;

