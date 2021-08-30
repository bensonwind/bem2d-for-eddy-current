clc;clear;close all

pier_geometry = 'C:\Users\benson_wind\Documents\MATLAB\bem2d_pier\pier_fine.nas';
[node,elem] = nastranmat_comsol(pier_geometry);

figure

% [vertex,elem] = uniformbisect(vertex,elem);
% [vertex,elem] = uniformrefine(vertex,elem);
% showmesh(node,elem);
% findvertex(vertex)
% axis equal
% view([0 90])

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

%%
% n = 100;
% phi = transpose(linspace(2*pi,0,n+1));
% phi = phi(1:end-1);
% a = 2;
% vertex = [a.*cos(phi) a.*sin(phi)];
% edge = [transpose(1:n-1) transpose(2:n);n 1];
%%
% nn = 100;
% w = 3; h = 1;
% vertex1 = [transpose(linspace(-w/2,w/2,nn)) -h/2*ones(nn,1)];
% vertex2 = [w/2*ones(nn,1) transpose(linspace(-h/2,h/2,nn))];
% vertex3 = [transpose(linspace(w/2,-w/2,nn)) h/2*ones(nn,1)];
% vertex4 = [-w/2*ones(nn,1) transpose(linspace(h/2,-h/2,nn))];
% vertex = [vertex1;vertex2(2:end,:);vertex3(2:end,:);vertex4(2:end-1,:)];
% n = (nn-1)*4;
% edge = [transpose(1:n-1) transpose(2:n);n 1];
%%
% hold on
figure
% nn = zeros(size(edge,1),2);
for tt=1:size(edge,1)
    colel = 'r';
	plot(vertex(edge(tt,1:2),1),vertex(edge(tt,1:2),2),[colel ':'])
    hold on
    
    nvec = cross([vertex(edge(tt,2),1)-vertex(edge(tt,1),1) vertex(edge(tt,2),2)-vertex(edge(tt,1),2) 0],[0 0 1]);
    nvec = nvec./norm(nvec);
%     nn(tt,:) = [nvec(1) nvec(2)];
    quiver(mean(vertex(edge(tt,1:2),1)),mean(vertex(edge(tt,1:2),2)),nvec(1),nvec(2),'r','AutoScale','off','MaxHeadSize',0.1)
    
%     [~, ~, ~, nx, ny]=elemshape([vertex(edge(tt,1),:);vertex(edge(tt,2),:)],0);
%     quiver(mean(vertex(edge(tt,1:2),1)),mean(vertex(edge(tt,1:2),2)),nx,ny,'r','AutoScale','off','MaxHeadSize',0.5)
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
% k = (1+j)*sqrt(omega*mu*sigma/2);
Je = 1;

% midx = ( vertex(edge(:,1),1) + vertex(edge(:,2),1) )./2;
% midy = ( vertex(edge(:,1),2) + vertex(edge(:,2),2) )./2;
% [G0,dG0] = greendef1(k,vertex(100,:),midx,midy);
% figure
% plot3([midx;midx(1)],[midy;midy(1)],abs([G0;G0(1)]),'.-')

int_2element = @(p,e,i_to,i_from,k) int_3points_v2( ...
[(p(e(i_to,1),1)+p(e(i_to,2),1))./2 (p(e(i_to,1),2)+p(e(i_to,2),2))./2], ...
[p(e(i_from,1),1) p(e(i_from,1),2)], ...
[p(e(i_from,2),1) p(e(i_from,2),2)], k);

% int_2element = @(p,e,i_to,i_from,k) int_3points( ...
% [p(i_to,1,1) p(i_to,2)], ...
% [p(e(i_from,1),1) p(e(i_from,1),2)], ...
% [p(e(i_from,2),1) p(e(i_from,2),2)], k);

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
%         [h,g,u,v,f] = int_3points_v2([0 0], ...
%                         [vertex(edge(i,1),1) vertex(edge(i,1),2)], ...
%                         [vertex(edge(i,2),1) vertex(edge(i,2),2)], k);
        H(j,idx(i,:)) = H(j,idx(i,:)) + h;
        G(j,idx(i,:)) = G(j,idx(i,:)) + g;
        U(j,idx(i,:)) = U(j,idx(i,:)) + u;
        V(j,idx(i,:)) = V(j,idx(i,:)) + v;
        F(j,i) = F(j,i) + mu*Je*f1;
        F2(j,i) = F2(j,i) + mu*Je*f2;
    end
%     P0 = [(vertex(edge(j,1),1)+vertex(edge(j,2),1))./2 (vertex(edge(j,1),2)+vertex(edge(j,2),2))./2];
%     F2(j) = mu*Je*int_2d(node,elem,k,P0);
%     F(j,j) = mu*Je*1j/4/k*(2j/k);

    H(j,idx(j,:)) = H(j,idx(j,:)) + [1/4 1/4];
    U(j,idx(j,:)) = U(j,idx(j,:)) + [1/4 1/4];
    if mod(j,50)==0
       disp(['Iteration: ' num2str(j)]) 
    end
end

FI = F*ones(n,1);
FI2 = F2*ones(n,1);
% A = pinv(H)*F1;
% A = lsqr(H,F1);
% A = 4.6e-7*ones(n,1);
% % dAdn = -pinv(G)*FI;
% Q = pinv(V)*U*A;
% dAdn = Q;
% % A = pinv(H)*(G*ones(n,1)+FI);
A = pinv(H-G*pinv(V)*U)*FI;
% % Q = mu*ones(n,1);
% % A = inv(H)*(G*Q+F1);
figure
plot3([vertex(:,1);vertex(1,1)],[vertex(:,2);vertex(1,2)],abs([A;A(1,1)]),'.-')
xlabel('x/m','FontSize',14)
ylabel('y/m','FontSize',14)
zlabel('Magnetic vector potential','FontSize',14)
% 
J = Je - 1j*sigma*omega*A;
% figure
% % plot3([vertex(:,1);vertex(1,1)],[vertex(:,2);vertex(1,2)],abs([J;J(1,1)]),'.-')
% xlabel('x/m','FontSize',14)
% ylabel('y/m','FontSize',14)
% zlabel('Total current density/A\cdotm^{-2}','FontSize',14)
% grid on
% 
% % dAdn = @(r) mu*Je/k*besselj(1,k*r)/besselj(0,k*a);
% % dAdn(a)
% 
% newA = @(r) mu*Je/k^2-besselj(0,k*r)/k/besselj(1,k*a);
% newA(a)

% dAdn_from_A = @(P,r) -k*(P+mu*Je/k^2)/besselj(0,k*a)*besselj(1,k*r);
% 
% P = 1*ones(n,1);
% Q = pinv(G)*(H*P-FI);
% [mean(Q) std(Q)]
% Q_true = dAdn_from_A(P,a);
% [mean(Q_true) std(Q_true)]
% 
% A_from_dAdn = @(Q,r) -mu*Je/k^2 - Q/k/besselj(1,k*a)*besselj(0,k*r);
% 
% Q = 0*ones(n,1);
% P = pinv(H)*(G*Q+FI);
% [mean(P) std(P)]
% P_true = A_from_dAdn(Q,a);
% [mean(P_true) std(P_true)]
% 
% dAdn_from_A2 = @(r) -mu*Je/2*r;
% P2 = 0*ones(n,1);
% Q2 = pinv(V)*(U*P-FI2);
% [mean(Q2) std(Q2)]
% Q2_true = dAdn_from_A2(a);
% [mean(Q2_true) std(Q2_true)]
% 
% 
% % A_from_dAdn2 = @(Q,r) Q*a*log(r);
% % -mu*Je/k^2+(P+mu*Je/k^2)/besselj(0,k*a)*besselj(0,k*r)
% figure
% clf
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