function [h,g,u,v,f1,f2] = int_3points_v2(P0,P1,P2,k)
    
    persistent bpregular wfregular bpleft wfleft bpmiddle wfmiddle bpright wfright
    
    c = rcond([P1-P0;P2-P0]);
    eps = 1e-12;
    
    pxy = P0;
    elknxy = [P1;P2];

    % Gauss-Legendre order
    n=20;

    if isempty(bpregular)
       [bpregular,wfregular]=gaussrule(n);
       [bpleft,wfleft]=nsingrule(8,-1,0.002);
       [bpmiddle,wfmiddle]=nsingrule(8,0,0.002);
       [bpright,wfright]=nsingrule(8,1,0.002);
    end
    
    bp = bpmiddle;
    wf = wfmiddle;
%     bp = bpregular;
%     wf = wfregular;

    % obtain shape functions, normal vector and global coordinates of the integration points
    [psi, xq, yq, nx, ny]=elemshape(elknxy,bp);
    
%     jacobi=sqrt(nx.^2+ny.^2);    % jacobian = modulus of the normal vector
%     R = sqrt((xq-pxy(1)).^2+(yq-pxy(2)).^2); % Distances from IPs to collocation point
%     dRdn=((xq-pxy(1)).*nx+(yq-pxy(2)).*ny)./R;

    jacobi = norm(P2-P1)./2.*ones(size(psi,1),1);
    nx=nx./jacobi; ny=ny./jacobi; % normalize normal vector
    r = [xq-pxy(1) yq-pxy(2)];
    r = r./vecnorm(r,2,2);
    dRdn = dot(r,[nx ny],2);
    
	[G1,dG1dR]=greendef1(k,pxy,xq,yq);
       
	F1 = excitationdef1(k,pxy,xq,yq);
    if c < eps
        f1 = 1j/4/k*(2j/k);
        f1 = 0;
    else
        f1 = wf'*(F1.*dRdn.*jacobi);
    end

       %Calculates h1, h2 and h3 without the corrective term

	dG1dn=dG1dR.*dRdn;
    if c < eps
        h1=1j/4*1/2*k*(2j/k);
        h1=0;
        h2=h1;
    else
        h1=wf'*(psi(:,1).*dG1dn.*jacobi);
        h2=wf'*(psi(:,2).*dG1dn.*jacobi);
    end
	
%        h3=wf'*(psi(:,3).*dG0dn);

	g1=wf'*(psi(:,1).*G1.*jacobi);
	g2=wf'*(psi(:,2).*G1.*jacobi);
        
	h = [h1 h2];
	g = [g1 g2];
    
    [G2,dG2dR]=greendef2(pxy,xq,yq);
    
    F2 = excitationdef2(pxy,xq,yq);
    if c < eps
        f2 = 1j/4/k*(2j/k);
        f2 = 0;
    else
        f2 = wf'*(F2.*dRdn.*jacobi);
    end
    
    
% 	dG2dn=dG2dR.*(-dRdn);
    dG2dn=dG2dR.*dRdn;
    if c < eps
        u1=0;
        u2=u1;
    else
        u1=wf'*(psi(:,1).*dG2dn.*jacobi);
        u2=wf'*(psi(:,2).*dG2dn.*jacobi);
    end
	
    
    v1=wf'*(psi(:,1).*G2.*jacobi);
	v2=wf'*(psi(:,2).*G2.*jacobi);

	u = [-u1 -u2];
	v = [-v1 -v2];

end

