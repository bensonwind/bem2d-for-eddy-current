function [F0]=excitationdef2(pxy,xq,yq)

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
%   -G0dir:    direct Green's function for free-field or rigid plane.
%   -G0ref:    reflected Green's function for free-field or rigid plane.
%   -dG0dirdR1:derivative of the direct Green's function.
%   -dG0dirdR2:derivative of the reflected Green's function.
%   -Pbeta:   correction term for a plane with finite impedance.
%   -dPbetadx,dPbetady :derivatives of the correction term.
%
% Reference:  -S.N. Chandler-Wilde and D.C. Hothersall: "Efficient calculation
%             of the Green function for acoustic propagation above an homogeneous
%             impedance plane", Journal of Sound and Vibration (1995) 180(5),
%             pp. 705-724.

% Vicente Cutanda 5-2001.

	R=sqrt((xq-pxy(1)).^2+(yq-pxy(2)).^2); % Distances from IPs to collocation/source point


% Free field Green function and its derivative
    F0 = R/4/pi.*(log(1./R)+1/2);

