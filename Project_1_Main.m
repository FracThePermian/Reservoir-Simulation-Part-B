close all
clear all
clc

%1. When inserting a well, add the x,y location on line 14,15 with the
%BHP ratio inserted on line 16. Then change the rates on line 219.
global k phi
[para] = reservoir; %class to help clean up program
 x = linspace(para.dx,para.L-para.dx/2,para.NX); y = linspace(para.dy,para.W-para.dy/2,para.NY); %for graphing
P = 3700*ones(para.N,1);    BC = sparse(para.N,1);      P_B = zeros(para.N,1);
well_location = [2500,4050]; q_well = [1000]; %multiple constant 'rate-wells' locations.
inj_well_locations = ceil(well_location(1)/para.dx)+floor(well_location(2)/para.dy)*para.NX; %which is 5050

x_loc_constBHP_well = [5536,5474,3600,2352,2000]; %2000 for 6th well
y_loc_constBHP_well = [3500,4708,4937,3322,5000]; %5000 for 6th well
P_prod_well = 2000*[1,1,1,1,1.025]; %1.025 for 6th well
BHP_well_locations = ceil(x_loc_constBHP_well/para.dx)+floor(y_loc_constBHP_well/para.dy)*para.NX;
BC(BHP_well_locations) = -1;
P_B(BHP_well_locations) = P_prod_well;

[T, B, Q, J] =  TBQ_box_f(BC,P_B,inj_well_locations,q_well); %a box that spits out matrices T, B, Q
 
dt = 1; t = 0; t_end = 201; n = 1; inv = sparse(T+J+B/dt); %change
P_array = zeros(t_end/dt,para.N); elim_k = find(k<0.001); %No other option than to use find function
while t < t_end
    P_2 = P; P = inv\(B*P_2/dt+Q);
    P(elim_k,:) = NaN;
    P_array(n,:) = P;
    t = t + dt; n = n + 1;
end

[Xx,Yy] = meshgrid(x,y); %Makes the meshgrid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_days_plot2 = reshape(P_array(2,:),para.NX,para.NY); P_days_plot2 = P_days_plot2';
P_days_plot10 = reshape(P_array(10,:),para.NX,para.NY); P_days_plot10 = P_days_plot10';
P_plot20 = reshape(P_array(20,:),para.NX,para.NY); P_plot20 = P_plot20';
P_days_plot100 = reshape(P_array(100,:),para.NX,para.NY); P_days_plot100 = P_days_plot100';
%Figure 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1),surf(Xx,Yy,P_days_plot2)
set(gca,'Ydir','reverse') %I did this because others told me my res. looked backward
colorbar
colormap jet
shading interp
xlabel('Reservoir Length (ft)','FontSize',14)
ylabel('Reservoir Width (ft)','FontSize',14)
zlabel('Res. Pressure (psi)','FontSize',14)
title ('t = 2 days','FontSize',20)
view (-35,75)
subplot(2,2,2),surf(Xx,Yy,P_days_plot10)
set(gca,'Ydir','reverse')
colorbar
colormap jet
shading interp
xlabel('Reservoir Length (ft)','FontSize',14)
ylabel('Reservoir Width (ft)','FontSize',14)
zlabel('Res. Pressure (psi)','FontSize',14)
title ('t = 10 days','FontSize',20)
view (-35,75)
subplot(2,2,3),surf(Xx,Yy,P_plot20)
set(gca,'Ydir','reverse')
colorbar
colormap jet
shading interp
xlabel('Reservoir Length (ft)','FontSize',14)
ylabel('Reservoir Width (ft)','FontSize',14)
zlabel('Res. Pressure (psi)','FontSize',14)
title ('t = 20 days','FontSize',20)
view (-35,75)
subplot(2,2,4),surf(Xx,Yy,P_days_plot100)
set(gca,'Ydir','reverse')
colorbar
colormap jet
shading interp
xlabel('Reservoir Length (ft)','FontSize',14)
ylabel('Reservoir Width (ft)','FontSize',14)
zlabel('Res. Pressure (psi)','FontSize',14)
title ('t = 100 days','FontSize',20)
view (-35,75)

figure %Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1),surf(Xx,Yy,P_days_plot2)
set(gca,'Ydir','reverse')
colorbar
colormap jet
shading interp
xlabel('Reservoir Length (ft)','FontSize',14)
ylabel('Reservoir Width (ft)','FontSize',14)
zlabel('Res. Pressure (psi)','FontSize',14)
title ('t = 2 days','FontSize',20)
view (-35,15)
subplot(2,2,2),surf(Xx,Yy,P_days_plot10)
set(gca,'Ydir','reverse')
colorbar
colormap jet
shading interp
xlabel('Reservoir Length (ft)','FontSize',14)
ylabel('Reservoir Width (ft)','FontSize',14)
zlabel('Res. Pressure (psi)','FontSize',14)
title ('t = 10 days','FontSize',20)
view (-35,15)
subplot(2,2,3),surf(Xx,Yy,P_plot20)
set(gca,'Ydir','reverse')
colorbar
colormap jet
shading interp
xlabel('Reservoir Length (ft)','FontSize',14)
ylabel('Reservoir Width (ft)','FontSize',14)
zlabel('Res. Pressure (psi)','FontSize',14)
title ('t = 20 days','FontSize',20)
view (-35,15)
subplot(2,2,4),surf(Xx,Yy,P_days_plot100)
set(gca,'Ydir','reverse')
colorbar
colormap jet
shading interp
xlabel('Reservoir Length (ft)','FontSize',14)
ylabel('Reservoir Width (ft)','FontSize',14)
zlabel('Res. Pressure (psi)','FontSize',14)
title ('t = 100 days','FontSize',20)
view (-35,15)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%CONTOUR MAP SECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure %Figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1),contourf(Xx,Yy,P_days_plot2)
set(gca,'Ydir','reverse')
colorbar
colormap jet
shading interp
xlabel('Reservoir Length (ft)','FontSize',14)
ylabel('Reservoir Width (ft)','FontSize',14)
zlabel('Res. Pressure (psi)','FontSize',14)
title ('t = 2 days','FontSize',20)

subplot(2,2,2),contourf(Xx,Yy,P_days_plot10)
set(gca,'Ydir','reverse')
colorbar
colormap jet
shading interp
xlabel('Reservoir Length (ft)','FontSize',14)
ylabel('Reservoir Width (ft)','FontSize',14)
zlabel('Res. Pressure (psi)','FontSize',14)
title ('t = 10 days','FontSize',20)

subplot(2,2,3),contourf(Xx,Yy,P_plot20)
set(gca,'Ydir','reverse')
colorbar
colormap jet
shading interp
xlabel('Reservoir Length (ft)','FontSize',14)
ylabel('Reservoir Width (ft)','FontSize',14)
zlabel('Res. Pressure (psi)','FontSize',14)
title ('t = 20 days','FontSize',20)

subplot(2,2,4),contourf(Xx,Yy,P_days_plot100)
set(gca,'Ydir','reverse')
colorbar
colormap jet
shading interp
xlabel('Reservoir Length (ft)','FontSize',14)
ylabel('Reservoir Width (ft)','FontSize',14)
zlabel('Res. Pressure (psi)','FontSize',14)
title ('t = 100 days','FontSize',20)


J_well1 = full(J(745,745));
J_well2 = full(J(1014,1014));
J_well3 = full(J(1000,1000));
J_well4 = full(J(666,666));
J_well5 = full(J(1042,1042));

p1 = P_array(1:100,745);
p2 = P_array(1:100,1014);
p3 = P_array(1:100,1000);
p4 = P_array(1:100,666);
p5 = P_array(1:100,1042);
t = 1:100;
q1 = J_well1*(p1-P_prod_well(1))
q2 = J_well2*(p2-P_prod_well(2))
q3 = J_well3*(p3-P_prod_well(3))
q4 = J_well4*(p4-P_prod_well(4))
q5 = J_well5*(p5-P_prod_well(5))

figure %Figure 4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1),plot(t,q1)
xlabel('Time (days)','FontSize',14)
ylabel('Well Rate (STB/day)','FontSize',14)
title('Producing Well #1','FontSize',14)
subplot(2,2,2),plot(t,q2)
xlabel('Time (days)','FontSize',14)
ylabel('Well Rate (STB/day)','FontSize',14)
title('Producing Well #2','FontSize',14)
subplot(2,2,3),plot(t,q3)
xlabel('Time (days)','FontSize',14)
ylabel('Well Rate (STB/day)','FontSize',14)
title('Producing Well #3','FontSize',14)
subplot(2,2,4),plot(t,q4)
xlabel('Time (days)','FontSize',14)
ylabel('Well Rate (STB/day)','FontSize',14)
title('Producing Well #4','FontSize',14)
figure
plot(t,q5)
xlabel('Time (days)','FontSize',14)
ylabel('Well Rate (STB/day)','FontSize',14)
title('Producing Well #5','FontSize',14)

figure %Figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_q = q1+q2+q3+q4+q5
J0 = 6.33e-3*2*pi*para.h*k(830)/(para.mu*para.Bw*log(4*(0.28*(para.dx^2+para.dy^2)^0.5/2)))
P1_XXX = P_array(1:100,830)
P_wfXXX = P1_XXX+q_well/J0
plot(t,P_wfXXX)
xlabel('Time (days)','FontSize',14)
ylabel('BottomHole Pressure (psi)','FontSize',14)
title('Injector Well #1','FontSize',20)

figure %Cumulative Production
Np_before = cumsum(total_q)
plot(t,Np_before)
xlabel('Time (days)','FontSize',12)
ylabel('Cumulative Oil produced (STB)','FontSize',12)
title('Cumulative Oil Production vs. Time','FontSize',12)

%For project 1a plotting.
% P_plot50 = reshape(P_days(50,:),para.NX,para.NY);
% surf(I,J,P_plot50) %This is the 3-D plot of well pressure wrt res. dimensions
% colorbar
% xlabel('Reservoir Length (ft)','FontSize',14)
% ylabel('Reservoir Width (ft)','FontSize',14)
% zlabel('Reservoir Pressure(psi)','FontSize',14)
% title('3-D Pressure Visualization','FontSize',20)
% 
% figure %This prints our the contour map of pressures over the well 
% contourf(I,J,P_plot50)
% colorbar
% xlabel('Reservoir Length(ft)','FontSize',14)
% ylabel('Reservoir Width(ft)','FontSize',14)
% title('Reservoir Contour','FontSize',20)
% 
% figure %This plots the final analytical and numerical reservoir pressure wrt length 
% P_1 = P_days(1,:);  P_10 = P_days(10,:);  P_50 = P_days(50,:);
% P1=P_1(5001:5100);  P10=P_10(5001:5100);  P50=P_50(5001:5100);
% plot(x,P1,'r+',x,P_analytical(x-10000,1),'r-',x,P10,'b+',x,P_analytical(x-10000,10),'b-',x,P50,'k+',x,P_analytical(x-10000,50),'k-');
% xlabel('Reservoir Length(ft)','FontSize',14)
% ylabel('Reservoir Pressure(psi)','FontSize',14)
% title('Reservoir Pressure vs Reservoir Length (Numerical vs. Analytical)','FontSize',20)
% legend('t = 1 day','Analytical 1 day','t = 10 days','Analytical 10 days','t = 50days','Analytical 50 days')
% 
