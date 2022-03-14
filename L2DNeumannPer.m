function [L] = L2DNeumannPer(Ne, hx, hy)
% Constructs the 2D Laplacian for a square mesh using the standard 5-point
% stencil with ghost cell at boundaries
%
% [L] = lap2DNeumann( Ne, h )
%
%  Parameters:
%       Ne    = the number of mesh points in each direction
%       h     = the mesh width
%
%
%  Return:
%       L = discrete Laplacian for zero Neumann BCs
%
%  Created on 23 June 2020
%          by Ying Zhang (yingzhang@brandeis.edu)
%
%

M = Ne * Ne;
e = ones(Ne,1);

% discrete Laplacian in y
MBlock = -2*speye(Ne);
OdBlock = speye(Ne);
Ly = blktridiag(MBlock, OdBlock, OdBlock, Ne);
Ly(1:Ne, Ne+1:2*Ne) = 2*Ly(1:Ne, Ne+1:2*Ne);
Ly((M-Ne+1):M, (M-2*Ne+1):(M-Ne)) = 2*Ly((M-Ne+1):M, (M-2*Ne+1):(M-Ne));

% discrete Laplacian in x
OdBlock = 0*OdBlock;
MBlock = spdiags([e -2*e e], [-1:1], Ne, Ne );
MBlock(1,Ne) = 1;
MBlock(Ne,1) = 1;
Lx = blktridiag(MBlock, OdBlock, OdBlock, Ne);

L = Lx/(hx*hx) + Ly/(hy*hy);

