function P = P_analytical(r,t)
Pi = 1000;q = 1000/5.615;mu = 1;k = 25;h = 10;phi = 0.2;ct = 1e-6
P = Pi -70.6*q*mu/(h*k)*expint((39.516*phi*mu*ct*(r).^2)/(k*t));
