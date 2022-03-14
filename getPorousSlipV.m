function[Por_Mat,nX,nY] = getPorousSlipV(dq, xLag, yLag, F_Lag, Kp)
% Computes porous slip velocity based on Darcy's Law
%
%   [Por_Mat, nX, nY] = getPorousSlipV( dq, xLag, yLag, F_Lag, Kp )
%
%  Parameters:
%       xLag  = a vector of the x-coordinates of the Lagrangian points
%       yLag  = a vector of the y-coordinates of the Lagrangian points
%       F_Lag = force density associated with the Lagrangian points
%       Kp    = porous slip parameter/permeability coefficient
%
%  Return:
%       Por_Mat = the modified Lagrangian point velocity
%       nX, nY  = unit normal vectors in x and y
%
%
%  Created on 30 Sept 2021
%          by Ying Zhang (yingzhang@brandeis.edu)
%

% number of porous media pts.
Np = length( xLag );

% Compute Lagrangian Derivatives
[xL_q, yL_q] = getLagDs(dq, Np, xLag, yLag);

% Compute Normal Vector (unit normals)
[nX, nY, sqrtNorm] = getLagNVecs(xL_q, yL_q);

% Compute Porous Slip Velocity
Up_X = - Kp .* F_Lag(:,1).*nX ./ sqrtNorm;
Up_Y = - Kp .* F_Lag(:,2).*nY ./ sqrtNorm;
Por_Mat = [Up_X Up_Y];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HELPER: computes Lagrangian Derivatives
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xL_q, yL_q] = getLagDs(dq, Np, xL, yL)


xL_q = zeros(Np,1);
yL_q = zeros(Np,1);


for( i = 1:Np )
    c = i-Np;
   if( c == -(Np-1) )
       xL_q(i) = ( -25/12*xL(i) + 4*xL(i+1) - 3*xL(i+2) + 4/3*xL(i+3) - 1/4*xL(i+4) ) / dq;
       yL_q(i) = ( -25/12*yL(i) + 4*yL(i+1) - 3*yL(i+2) + 4/3*yL(i+3) - 1/4*yL(i+4) ) / dq;
       
   elseif( c == -(Np-2) )
       xL_q(i) = ( -0.25*xL(i-1) - 5/6*xL(i) + 1.5*xL(i+1) - 0.5*xL(i+2) + 1/12*xL(i+3) ) / dq;
       yL_q(i) = ( -0.25*yL(i-1) - 5/6*yL(i) + 1.5*yL(i+1) - 0.5*yL(i+2) + 1/12*yL(i+3) ) / dq;

   
   elseif( c >= -(Np-3) && c < -1 )
       xL_q(i) = ( 1/12*xL(i-2) - 2/3*xL(i-1) + 2/3*xL(i+1) - 1/12*xL(i+2) ) / dq;
       yL_q(i) = ( 1/12*yL(i-2) - 2/3*yL(i-1) + 2/3*yL(i+1) - 1/12*yL(i+2) ) / dq;

       
   elseif( c == -1 )
       xL_q(i) = ( -1/12*xL(i-3) + 0.5*xL(i-2) - 1.5*xL(i-1) + 5/6*xL(i) + 0.25*xL(i+1) ) / dq;
       yL_q(i) = ( -1/12*yL(i-3) + 0.5*yL(i-2) - 1.5*yL(i-1) + 5/6*yL(i) + 0.25*yL(i+1) ) / dq;
            
   elseif( c == 0 )
       xL_q(i) = ( 0.25*xL(i-4) - 4/3*xL(i-3) + 3*xL(i-2) - 4*xL(i-1) + 25/12*xL(i) ) / dq;
       yL_q(i) = ( 0.25*yL(i-4) - 4/3*yL(i-3) + 3*yL(i-2) - 4*yL(i-1) + 25/12*yL(i) ) / dq;
       
   else
      fprintf('\n\n');
      error('Error: Bad file format inside your .porous file!\n'); 
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HELPER: computes Lagrangian UNIT Normal Vectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nX,nY,sqrtN] = getLagNVecs(xL_q, yL_q)

sqrtN = sqrt( (xL_q).^2 + (yL_q).^2 );

nX = ( yL_q ) ./ sqrtN;
nY = ( -xL_q ) ./ sqrtN;