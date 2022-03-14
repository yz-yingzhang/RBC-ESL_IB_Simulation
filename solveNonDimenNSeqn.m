function[u,p] = solveNonDimenNSeqn( u0, p0, f, fT, Dpx, Dpy, G0x, G0y, ...
    L, IL1, IL2, ILp )
% Solve the non-dimensionalized incompressible unsteady Stokes Equation
%
%  [u, p] = solveNonDimenNSeqnStag2( u0, p0, f, fT, rho, D0x, D0y, Lp, IL1, IL2, dt )
%
%  Parameters:
%       u0    = current Eulerian velocity field
%       p0    = pressure calculated from the previoud time step
%       f     = current Eulerian force
%       fT    = Eulerian force calculated using X^*
%       Dpx, Dpy, G0x, G0y, L, IL1, IL2, ILp    
%             = discrete operators
%
%  Return:
%       u   = Eulerian velocity at time dt after u0 solution
%       p   = Eulerian pressure at time dt
%
%
%  Created on 13 Aug 2020
%          by Ying Zhang (yingzhang@brandeis.edu)
%
%

u = u0;
M = size(u0,1);
Ne = sqrt(M);
uM = u;

w = IL2 * u0 + ((1/2)*(f+fT) - [G0x*p0 G0y*p0]);

IL1(1:Ne,:) = 0;
for(i = 1:Ne)
    IL1(i,i) = 1;
end
w(1:Ne,:) = 0;
w((M-Ne+1):M,:) = 0;
uM = IL1\w;

% solve for phi using uM
L(1,:) = 1;
DDuM = (Dpx*uM(:,1)+Dpy*uM(:,2));
DDuM(1) = 0;
phi = L\DDuM;

% update the pressure for the next time step
p = p0 + ILp*phi;

% solve for u at t+dt by projection
GxPhi = G0x*phi;
GyPhi = G0y*phi;
u = uM - [GxPhi GyPhi];

