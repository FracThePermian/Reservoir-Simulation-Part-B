# Reservoir-Simulation-Part-B

Compares analytical/numerical results for the drainage of multiple wells w/n the dimensions that you define (drainage area and depth). The Nechelik (.dat) files provide porosity and permeability for each unit block of reservoir. Modify the following variables to your liking.

Control the placement of wells and whether they are injecting and/or producing.

It is written in Matlab, is self-contained, and no external dependencys. 

### Assumptions
  * Incompressible Fluid
  * Homogeneous/Heterogeneous Reservoir
  * Neumann boundary (Reservoir side)
  * Dirichlet boundary (Well side)
### Modify Parameters
 * Reservoir Depth, Area, and Volume (ft ; ft^2 ; ft^3)
 * Reservoir Pressure (psi)
 * Time and/or Iterations (days)
 * Porosity (phi)
 * Permeability (k)
 * Compressibility (ct)
 * Flow Rate (+/- Q)
 * Damage (hk)

## Method and Materials
* The Finite Element Method (FEM) is an essential numerical tool for solving boundary value problems of PDE's. 

**Snippet of module that generates sparsity matrix**
```matlab
for i = 1:N  		%Generate sparsity matrix
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
    if BC(i) ~= 0 %if and only if the boundary condition exists
        T(i,i) = T(i,i)+BC(i)*2*T_frac(i,i,Ax,mu,Bw,dx);
    end
```

**Analytical Solution to Pressure Profile**
```matlab
function P = P_analytical(r,t)
Pi = 1000;q = 1000/5.615;mu = 1;k = 25;h = 10;phi = 0.2;ct = 1e-6
P = Pi -70.6*q*mu/(h*k)*expint((39.516*phi*mu*ct*(r).^2)/(k*t));
```

## Installation

Download package:
```git
git clone https://github.com/FracThePermian/Reservoir-Simulation-Part-B
```

Navigate to folder from Matlab (~/Reservoir-Simulation-Part-B/) and run the Project_1_Main.m script.


## Results

**3D Pressure Profile**

![Pressure Visualization](https://github.com/FracThePermian/Reservoir-Simulation-Part-B/blob/master/Graphical-Output/3d_contour.png?raw=true "Pressure Gain/Reduction Profile")

**Pressure topology**

![Pressure Contours](https://github.com/FracThePermian/Reservoir-Simulation-Part-B/blob/master/Graphical-Output/Pcontour1.png?raw=true "Pressure Contours")

**3D Pressure Profile**

![Pressure Contour 2](https://github.com/FracThePermian/Reservoir-Simulation-Part-B/blob/master/Graphical-Output/Pcontour2.png?raw=true "Pressure Contour 2")

**Flow Rate vs. Time**

![Rate vs. Time](https://github.com/FracThePermian/Reservoir-Simulation-Part-B/blob/master/Graphical-Output/RateVTime.png?raw=true "Q vs. Time")

**Cumulative Oil Production**

![Cumulative Oil Production](https://github.com/FracThePermian/Reservoir-Simulation-Part-B/blob/master/Graphical-Output/CumulativeOilProduction.png?raw=true "Cumulative Oil Production")


### License
[MIT](https://opensource.org/licenses/MIT "MIT")
