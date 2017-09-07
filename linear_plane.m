clear all;
close all
clc

%% TUT IIA

%%
tic
p0 = 20e-6;
c = 343; %m/s
f_pw = 870; %Hz
pA = sqrt(2);
w_pw = 2*pi*f_pw;
k = w_pw/c;
y_ref = 1; 

d_XY = 0.05;
angle = [45 60 90];
d_x0 = [0.1 0.4]; %Abständ zwischen den Punktquellen

for i = 1:length(angle)
    for j = 1:length(d_x0)
% i =1;j=1; % hier dann die schleife
theta_PW = angle(i) * pi/180; 
x_0 = -50:d_x0(j):50;  %Quellverteilung

k_ypw = (w_pw*sin(theta_PW))/c;
k_xpw = (w_pw*cos(theta_PW))/c;
%1m referenzabstand
X = -5.5:d_XY:5.5;
Y = -1:d_XY:10;
X(find(X == 0)) = 0.01;
Y(find(Y == 0)) = 0.01;
[x,y] = meshgrid(X,Y);

%%

D25 = 4 * 1j * pA * exp(-1j .* k_ypw .* y_ref) .* exp(-1j .* k_xpw .* x_0) ./...
    besselh(0,2,k_ypw .* y_ref);


%%

for n = 1:length(x_0)
    dist(:,:,n) = sqrt((x-x_0(n)).^2 + y.^2);
end

%%
G = 1/(4*pi) * exp(-1j * k * dist) ./ dist;

%% 
P2D = zeros(size(x,1),size(x,1));

for a = 1:length(x_0)
    t = G(:,:,a) * D25(a);
    P2D = P2D + t;
end

%% plot

if j==1 
   figure;
end
subplot(1,2,j)
surf(x,y,real(P2D/10));
grid on;
view(0,90)
axis equal;
shading flat;
cstep = 20;
cmin = -3; cmax = 3;
cm = cbrewer('seq', 'YlOrRd', cstep);
colormap(cm);
set(gca, 'CLim', [cmin cmax]);
cbh = colorbar;
set(cbh, 'YTick');
xlabel('x  [m]');
ylabel('y  [m]');
axis([-5.5 5.5 -1 10]);
str = ['R[P(x,\omega_P_W)]/10  \theta_P_W =' num2str(angle(i)) '°  \Deltax_0 = ' num2str(d_x0(j)) 'm']; 
title(str);
set(gcf, 'Position', get(0, 'Screensize'));


    end
end

toc