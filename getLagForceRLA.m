function[F] = getLagForceRLA( Xl, Nl, dtheta, Ks, Kb, Ka, DA )
% Calculates the force at Lagrangian markers for with resting length of
% elasticity and area-conservation
%
%  [F] = getLagForceRLA( Xl, Nl, dtheta, Ks, Kb, Ka, DA )
%
%  Parameters:
%       Xl     = input Lagrangian marker position
%       Nl     = number of Lagrangian markers
%       dtheta = Lagrangian marker spacing
%       Ks     = stretching coefficient
%       Kb     = bending rigidity
%       Ka     = area-preserving constant
%       DA     = difference in area: (A - A0)
%
%  Return:
%       F      = the Lagrangian force
%
%
%  Created on 30 Sept 2021
%          by Ying Zhang (yingzhang@brandeis.edu)
%
%

% spring force
e = ones(Nl,1);
Lap = spdiags([e -2*e e], [-1:1], Nl, Nl);

% periodic
Lap(1,Nl) = 1;
Lap(Nl,1) = 1;

Lap = Lap ./ (dtheta * dtheta);
LX = Lap * Xl;

idxUp = [2:Nl 1]';
idxDn = [Nl 1:(Nl-1)]';

tmp    = Xl(idxUp,:) - Xl;
LUp    = sqrt(dot(tmp,tmp,2));
TUp    = [(tmp(:,1) ./ LUp) (tmp(:,2) ./ LUp)];

tmp    = Xl - Xl(idxDn,:);
LDn    = sqrt(dot(tmp,tmp,2));
TDn    = [(tmp(:,1) ./ LDn) (tmp(:,2) ./ LDn)];

Fs      = zeros(Nl, 2);
Fs(:,1) = (1/(dtheta)) * ( -TUp(:,1) + TDn(:,1) );
Fs(:,2) = (1/(dtheta)) * ( -TUp(:,2) + TDn(:,2) );

% elastic force
Fs = Ks * (Fs + LX);

% bending force
Fb = -Kb * Lap * LX;

% area-preserving force
Fa = [(Xl(idxUp,2)-Xl(idxDn,2)) (-Xl(idxUp,1)+Xl(idxDn,1))];
Fa = -(Ka*DA/2) * Fa;

% the total force density
F = Fs + Fb + Fa;