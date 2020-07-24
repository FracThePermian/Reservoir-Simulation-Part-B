function T_frick = T_frac(cha,chacha,Ax,mu,Bw,dx)
global k
k_frac = 2*(1/k(cha)+1/k(chacha))^-1;
T_frick = 6.33e-3*k_frac*Ax/(Bw*mu*dx);
