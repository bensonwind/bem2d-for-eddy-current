function [Gdir,dGdirdR]=greendef1(k,pxy,xq,yq)

% [G0dir,G0ref,dG0dirdR1,dG0refdR2,Pbeta,dPbetadx,dPbetady]=greendef(k,pxyb,betaP,xq,yq);
%
% Calculates the 2D Green's function for free field, rigid plane or
% plane with impedance.
%
% Note:
%   To obtain the incident pressure on the nodes from a line source in the absence
%   of the body using this function, do like this:
%
%    [G0dir,G0ref,dG0dirdR1,dG0refdR2,Pbeta,dPbetadx,dPbetady]=greendef(k,position,betaP,xyb(:,1),xyb(:,2));
%    inc_pressure=i/4*(G0dir+G0ref)+Pbeta;
%
%   where 'position' are the [x,y] coordinates of the line source and xyb is the node
%   positions geometry file. 
%
%
% Input variables:
%   -k:       wavenumber.
%   -pxyb:    real vector containing the (x,y,body) values for
%             the point 'P', source or collocation point.
%   -xq:      real column vector containing the global x-coordinates
%             for each point to calculate (e.g. integration points). 
%   -yq:      real column vector containing the global y-coordinates
%             for each point to calculate (e.g. integration points). 
%
%  Output variables:
%   -Gdir:    direct Green's function for free-field or rigid plane.
%   -dGdirdR:derivative of the direct Green's function.

%
% Reference:  -S.N. Chandler-Wilde and D.C. Hothersall: "Efficient calculation
%             of the Green function for acoustic propagation above an homogeneous
%             impedance plane", Journal of Sound and Vibration (1995) 180(5),
%             pp. 705-724.

% Vicente Cutanda 5-2001.

	R=sqrt((xq-pxy(1)).^2+(yq-pxy(2)).^2); % Distances from IPs to collocation/source point
    
% Free field Green function and its derivative
	Gdir=1j/4*besselh(0,1,k*R,1);
	dGdirdR=-1j/4*k*besselh(1,1,k*R,1);
