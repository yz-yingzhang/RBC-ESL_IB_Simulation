function[f] = spreadLagForcePhysBCs( F, Xl, dtheta, Ne, hx, hy )
% Spreads a force on the Lagrangian markers, Xl, to the fluid 
%     at pts (X,Y)
%
%  [f] = spreadLagForcePhysBCs( F, Xl, dtheta, Ne, hx, hy )
%
%  Parameters:
%       F      = Lagrangian force
%       Xl     = input Lagrangian marker position
%       dtheta = Lagrangian marker spacing
%       Ne     = number of Eulerian mesh points in each direction
%       hx, hy = Eulerian mesh spacing
%
%  Return:
%       f      = the Eulerian force
%
%
%  Created on 19 June 2020
%          by Ying Zhang (yingzhang@brandeis.edu)
%
%

% evaluate the discrete delta function given the Lagrangian markers and the
% Eulerian grid
% returns indices of the Eulerian points where the weigths should be
% applied
[idxs, delta] = evalDeltaPhysBCs(Xl, Ne, hx, hy);

f = zeros(Ne*Ne,2);
XlLen = length(Xl(:,1));

% for each Lagrangian marker calculate the force to apply
for( i = 1:XlLen )
    
   f(idxs(i,:),1) = f(idxs(i,:),1) + F(i,1) * dtheta * delta(i,:)';
   f(idxs(i,:),2) = f(idxs(i,:),2) + F(i,2) * dtheta * delta(i,:)';
   
end