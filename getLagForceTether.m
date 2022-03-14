function[F] = getLagForceTether( Xl, Tl, KT )
% Calculates the force at Lagrangian markers for bending filament
%
%  [F] = getLagForceTether( Xl, Tl, Nl, dtheta, K, isPer )
%
%  Parameters:
%       Xl     = input Lagrangian marker position
%       Tl     = input tether marker position
%       KT     = tether force coefficient
%
%  Return:
%       F      = the Lagrangian force
%
%
%  Created on 26 May 2020
%          by Ying Zhang (phzhang@bu.edu)
%
%

% tether force
F = -KT * (Xl - Tl);