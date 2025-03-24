clear
close all

param_init
syms p v z
h_x = z - Th*v - 0.5*(v-vl)^2/(cd*g);
dhx = [diff(h_x,p) diff(h_x,v) diff(h_x,z)];
hx = matlabFunction(h_x);
dhx = matlabFunction(dhx);
V = (v - vd)^2;
dV = [diff(V,p) diff(V,v) diff(V,z)];
mV = matlabFunction(V);
mdV = matlabFunction(dV);
observer = 0;
storage = [];

%% with observer
for t = 0:t_step:T

    z = x_c(3);
    v = x_c(2);
    p = x_c(1);

    V = mV(v);
    Fr = GetFr(f0,f1,f2,v);
    fx = [v;-Fr/m;vl - v];
    g1 = [0;1/m;0];
    g2 = [0;g/m;0];
    kdx = -g;

    dVdx = mdV(v);

    LfV = dVdx*fx;
    Lg1V = dVdx*g1;
    lg2V = dVdx*g1;
    Aclf = [Lg1V,-1];
    % Bclf = -alpha_v * V-LfV-lg2V*hat_d;
    Bclf = -alpha_v * V-LfV;

    h = hx(v,z);
    dhdx = dhx(v);
    Lfh = dhdx*fx;
    Lg1h = dhdx*g1;
    Lg2h = dhdx*g2;
    Acbf = [-Lg1h,0];
    if observer == 1
        Bcbf = beta_c * h + Lfh - (Lg2h^2)/(2*beta*(2*beta_e - beta_c));
    else
        Bcbf = beta_c * h + Lfh;
    end
    
    H0 = 2/m^2;
    mu = 2*(0.01);
    H = diag([2*H0,2*mu]);
    f = [-(4/m^2)*Fr;0];
    Au = [Aclf;Acbf];
    Bu = [Bclf;Bcbf];
    out = quadprog(H,f,Au,Bu);
    hat_u = out(1);
    d = 200*cos(2*pi*t);
    u = out(1) + kdx*hat_d;
    [hat_d,z_new] = DisturbanceObserver(A,C,g2,rx*v,sx,fx,g1,doz,u,t,t_step);
    if hat_d > 200
        hat_d = 200;
    elseif hat_d < -200
        hat_d = -200;
    end
    de = hat_d - d;
    he = h - beta*(0.5*de^2);

    doz = z_new;

    u = InputConstraints(u ,cd,ca,m,g);
    sigma = out(2);
    
    [~,Y] = ode45(@(t,y) Dynamics(t,y,u,d,vl,f0,f1,f2,m,g),[t,t+t_step],x_c);
    x_c = Y(end,:);

    storage.Dox(cnt,:) = x_c';
    storage.Dou(cnt,:) = u';
    storage.Dosigma(cnt,:) = sigma';
    storage.Doh(cnt,:) = h';
    storage.Dode(cnt,:) = de';
    storage.DoV(cnt,:) = V';
    storage.Dohe(cnt,:) = he';
    cnt = cnt+1;

end

%% without observer
param_init
for t = 0:t_step:T

    z = x_c(3);
    v = x_c(2);
    p = x_c(1);

    V = mV(v);
    Fr = GetFr(f0,f1,f2,v);
    fx = [v;-Fr/m;vl - v];
    g1 = [0;1/m;0];
    g2 = [0;g/m;0];
    kdx = -g;

    dVdx = mdV(v);

    LfV = dVdx*fx;
    Lg1V = dVdx*g1;
    lg2V = dVdx*g1;
    Aclf = [Lg1V,-1];
    % Bclf = -alpha_v * V-LfV-lg2V*hat_d;
    Bclf = -alpha_v * V-LfV;

    h = hx(v,z);
    dhdx = dhx(v);
    Lfh = dhdx*fx;
    Lg1h = dhdx*g1;
    Lg2h = dhdx*g2;
    Acbf = [-Lg1h,0];
    if observer == 1
        Bcbf = beta_c * h + Lfh - (Lg2h^2)/(2*beta*(2*beta_e - beta_c));
    else
        Bcbf = beta_c * h + Lfh;
    end
    
    H0 = 2/m^2;
    mu = 2*(0.01);
    H = diag([2*H0,2*mu]);
    f = [-(4/m^2)*Fr;0];
    Au = [Aclf;Acbf];
    Bu = [Bclf;Bcbf];
    out = quadprog(H,f,Au,Bu);
    hat_u = out(1);
    d = 200*cos(2*pi*t);
    u = out(1);
    u = InputConstraints(u ,cd,ca,m,g);
    sigma = out(2);
    
    [~,Y] = ode45(@(t,y) Dynamics(t,y,u,d,vl,f0,f1,f2,m,g),[t,t+t_step],x_c);
    x_c = Y(end,:);

    storage.x(cnt,:) = x_c';
    storage.u(cnt,:) = u';
    storage.sigma(cnt,:) = sigma';
    storage.h(cnt,:) = h';
    storage.de(cnt,:) = de';
    storage.V(cnt,:) = V';
    storage.he(cnt,:) = he';
    cnt = cnt+1;

end
close all
figure
t = 0:t_step:T;
subplot(2,2,1)
hold on
plot(t,vd*ones(size(t)),"Color",[1,0,0],LineStyle="--");
hold on
plot(t,vl*ones(size(t)),"Color",[0,1,0],LineStyle="--");
hold on
plot(t,storage.Dox(:,2));
hold on
plot(t,storage.x(:,2),"Color",[0,0,1],LineStyle="--");

subplot(2,2,2)
plot(t,storage.Dox(:,3));
hold on
plot(t,storage.x(:,3),"Color",[0,0,1],LineStyle="--");

subplot(2,2,3)
hold on
plot(t,cd*m*g*ones(size(t)),"Color",[1,0,0],LineStyle="--");
hold on
plot(t,-cd*m*g*ones(size(t)),"Color",[0,1,0],LineStyle="--");
hold on
plot(t,storage.Dou);
hold on
plot(t,storage.u,"Color",[0,0,1],LineStyle="--");

subplot(2,2,4)
plot(t,storage.Dosigma);
hold on
plot(t,storage.sigma,"Color",[0,0,1],LineStyle="--");

figure
subplot(2,1,1)
plot(t,storage.Doh);
hold on

plot(t,storage.h,"Color",[0,0,1],LineStyle="--");
subplot(2,1,2)
plot(t,storage.DoV);
hold on
plot(t,storage.V,"Color",[0,0,1],LineStyle="--");
figure
plot(t,storage.Dode);
function dy =  Dynamics(~,s,u,d,vl,f0,f1,f2,m,g)
    % p = s(1);
    v = s(2);
    dz = vl - v;
    dp = v;
    Fr = GetFr(f0,f1,f2,v);
    dv = - Fr/m + 1/m *u + 1/m *g*d;
    dy = [dp;dv;dz];
end

function Fr = GetFr(f0,f1,f2,v)
    Fr = f0 + f1*v + f2*v^2;
end

function u_out = InputConstraints(u,cd,ca,m,g)
    u_out = u;   
    if u > ca*m*g
        u_out = ca*m*g;
    elseif u < -cd*m*g
        u_out = -cd*m*g;
    end
end

% function h = GetHx(z,v,Th,vl,cd,g)
%     h = z - Th*v - 0.5*(v-vl)^2/(cd*g);
% end

function [hat_d,z_new] = DisturbanceObserver(A,C,g2,rx,sx,f,g1,z,u,t,t_step)
    [~,Y] = ode45(@(t,z) DoDyamic(t,z,A,C,g2,rx,sx,f,g1,u),[t,t+t_step],z);
    z_new = Y(end,:)';
    hat_xi = z_new + rx;
    hat_d = C*hat_xi;
    % dz = z;
    function dz = DoDyamic(~,z,A,C,g2,rx,sx,f,g1,u)
        dz = (A - sx*g2*C)*z + A*(rx) - sx*(g2*C*(rx) + f +g1*u);
    end
end
