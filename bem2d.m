clc;clear;close all

pier_geometry = 'pier_fine.nas';
[node,elem] = nastranmat_comsol(pier_geometry);

% Find boundary edges
bound = findboundary(elem);
bound = sortboundary(bound,false);

n = size(bound,1);
vertex = node(bound(:,1),:);
edge = [transpose(1:n) transpose([2:n 1])];

vertex = vertex./10;
xw = max(vertex(:,1)) - min(vertex(:,1));
yw = max(vertex(:,2)) - min(vertex(:,2));
vertex(:,1) = vertex(:,1)-xw/2;
vertex(:,2) = vertex(:,2)-yw/2;

figure

for tt=1:size(edge,1)

    colel = 'r';
	plot(vertex(edge(tt,1:2),1),vertex(edge(tt,1:2),2),[colel ':'])
    hold on
    
    nvec = cross([vertex(edge(tt,2),1)-vertex(edge(tt,1),1) vertex(edge(tt,2),2)-vertex(edge(tt,1),2) 0],[0 0 1]);
    nvec = nvec./norm(nvec);

    quiver(mean(vertex(edge(tt,1:2),1)),mean(vertex(edge(tt,1:2),2)),nvec(1),nvec(2),'r','AutoScale','off','MaxHeadSize',0.1)
    
end
xlabel('x/m','FontSize',14)
ylabel('y/m','FontSize',14)
axis equal
%%
freq = 1e5;
omega = 2*pi*freq;
%Iron
sigma = 4e6; mu_r = 5e3;
% Copper
% sigma = 57e6; mu_r  = 1;
mu0 = pi*4e-7;
mu = mu_r * mu0;
k = sqrt(-1j*omega*mu*sigma);

Je = 1;

int_2element = @(p,e,i_to,i_from,k) int_3points(...
[(p(e(i_to,1),1)+p(e(i_to,2),1))./2 (p(e(i_to,1),2)+p(e(i_to,2),2))./2], ...
[p(e(i_from,1),1) p(e(i_from,1),2)], ...
[p(e(i_from,2),1) p(e(i_from,2),2)], k);

H = zeros(n,n);
G = zeros(n,n);
U = zeros(n,n);
V = zeros(n,n);
F = zeros(n,n);
C = zeros(n,n);
F2 = zeros(n,n);
disp([' 2D BEM calculation, ' num2str(n) ' elements, k = ' num2str(k)])
idx = [transpose(1:n) transpose([2:n 1])];
% coefficient matrix calculation
for j=1:n
    for i=1:n
        [h,g,u,v,f1,f2] = int_2element(vertex,edge,j,i,k);

        H(j,idx(i,:)) = H(j,idx(i,:)) + h;
        G(j,idx(i,:)) = G(j,idx(i,:)) + g;
        U(j,idx(i,:)) = U(j,idx(i,:)) + u;
        V(j,idx(i,:)) = V(j,idx(i,:)) + v;
        F(j,i) = F(j,i) + mu*Je*f1;
	F2(j,i) = F2(j,i) + mu*Je*f2;
    end

    H(j,idx(j,:)) = H(j,idx(j,:)) + [1/4 1/4];
    U(j,idx(j,:)) = U(j,idx(j,:)) + [1/4 1/4];
    if mod(j,50)==0
       disp(['Iteration: ' num2str(j)]) 
    end
end

FI = F*ones(n,1);
FI2 = F2*ones(n,1);

A = pinv(H-G*pinv(V)*U)*FI;

figure
plot3([vertex(:,1);vertex(1,1)],[vertex(:,2);vertex(1,2)],abs([A;A(1,1)]),'.-')
xlabel('x/m','FontSize',14)
ylabel('y/m','FontSize',14)
zlabel('Magnetic vector potential','FontSize',14)

J = Je - 1j*sigma*omega*A;

figure
plot3([vertex(:,1);vertex(1,1)],[vertex(:,2);vertex(1,2)],abs([J;J(1,1)]),'.-','MarkerSize',8)
a = patch('XData',vertex(:,1),'YData',vertex(:,2),'FaceColor',[0.9 0.9 0.9],'LineStyle','-','FaceAlpha',1);
hh = hatchfill(a,'cross',45,10,[0.8 0.8 0.8]);
set(hh,'Color','w')
xdata = [transpose(vertex(edge(:,1),1)); transpose(vertex(edge(:,2),1)); transpose(vertex(edge(:,2),1)); transpose(vertex(edge(:,1),1))];
ydata = [transpose(vertex(edge(:,1),2)); transpose(vertex(edge(:,2),2)); transpose(vertex(edge(:,2),2)); transpose(vertex(edge(:,1),2))];
zdata = [zeros(1,n); zeros(1,n); transpose(J(edge(:,2))); transpose(J(edge(:,1)))]; 
patch(xdata,ydata,abs(zdata),[0.00,0.45,0.74],'LineStyle','none','FaceAlpha',0.5)
set(gca,'DataAspectRatio',[1 1 max(abs(J))/min([xw yw])])
view([-30 30])
box on
xl = xlabel('x/m','FontSize',14);
yl = ylabel('y/m','FontSize',14);
zlabel('Current density of boundary/A\cdotm^{-2}','FontSize',14)
set(xl,'Rotation',16);
set(yl,'Rotation',-49);
