function[phi] = evalPhi(r)
% Evaluate one phi function given r
%
%  [phi] = evalPhi( r )
%
%  Parameters:
%       r    = the points to evaluate the discrete delta function at
%
%  Return:
%       phi  = value of unscaled discrete delta function at each r
%
%
%  Created on 26 May 2020
%          by Ying Zhang (phzhang@bu.edu)
%
%

phi = zeros(length(r),1);

idx1 = find( 0 <= r & r <= 1 );
idx2 = find( r > 1 & r <= 2 );
idx3 = find( -1 <= r & r < 0 );
idx4 = find( -2 <= r & r < -1 );

reval = r(idx1);
phi(idx1) = 3 - 2*reval + sqrt(1 + 4*reval - 4*reval.*reval);

reval = r(idx2);
phi(idx2) = 5 - 2*reval - sqrt(-7 + 12*reval - 4*reval.*reval);

reval = r(idx3);
phi(idx3) = 3 + 2*reval + sqrt(1 - 4*reval - 4*reval.*reval);

reval = r(idx4);
phi(idx4) = 5 + 2*reval - sqrt(-7 - 12*reval - 4*reval.*reval);

phi = phi/8;
