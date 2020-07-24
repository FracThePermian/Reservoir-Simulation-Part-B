function [ T,B,Q,J ] = TBQ_box_f(BC,P_B,para_wells,q_wells)
global k phi
[para] = reservoir; N = para.N; NX = para.NX; NY = para.NY; mu = para.mu; Bw = para.Bw; dx = para.dx; dy = para.dy; Ax = para.Ax; Ay = para.Ay;
h = para.h; ct = para.ct; Bw = para.Bw
T = sparse(N,N); B = sparse(N,N); Q = sparse(N,1); Q(para_wells) = Q(para_wells) + q_wells; J = sparse(N,N);

for i = 1:N
    if i+NX <= N 
        T_diag = T_frac(i,i+NX,Ay,mu,Bw,dy);
        T(i,i+NX) = T(i,i+NX)-T_diag;
        T(i+NX,i) = T(i+NX,i)-T_diag;
    end
    
    % This is for diagonals toward the inside of sparsity
    if (mod(i,NX) ~= 0) && (i+1 <= N)  
        T_diag = T_frac(i,i+1,Ax,mu,Bw,dx);
        T(i,i+1) = T(i,i+1)-T_diag;
        T(i+1,i) = T(i+1,i)-T_diag;
    end
    
    T(i,i) = abs(sum(T(i,:)));%Sum up every row on every iteration
    
    if BC(i) == 1 %if and only if the boundary condition exists
        T(i,i) = T(i,i)+2*T_frac(i,i,Ax,mu,Bw,dx);
        Q(i) = Q(i) + 2*P_B(i)*T_frac(i,i,Ax,mu,Bw,dx);
    end
    
    if BC(i) == -1 %Because of the boundary conditions
        J(i,i) = 6.33e-3*2*pi*para.h*k(i)/(para.mu*para.Bw*log(4*(0.28*(para.dx^2+para.dy^2)^0.5/2)))
        Q(i) = Q(i) + P_B(i)*J(i,i)
    end
    
B(i,i) = dx*dy*h*phi(i)*ct/Bw; %phi as such for next project '2'

end
end




