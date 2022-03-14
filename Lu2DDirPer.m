function [L] = Lu2DDirPer(Ne, hx, hy)
%  Constructs the 2D Laplacian for a square mesh with zero Dirichlet 
%     boundary conditions using the standard 5-point stencil.
%
%  [L] = Lu2DDirPer( N, h )
%
%
%  Parameters:
%       Ne     = the number of mesh points in each direction
%       hx, hy = the mesh width in x and y
%
%  Return:
%       L = discrete Laplacian
%
%  Created on 23 June 2020
%          by Ying Zhang (phzhang@bu.edu)
%
%


M = Ne * Ne;
e = ones(Ne,1);

% discrete Laplacian in y
MBlock = -2*speye(Ne);
OdBlock = speye(Ne);
Ly = blktridiag(MBlock, OdBlock, OdBlock, Ne);
Ly(1:Ne, 1:Ne) = Ly(1:Ne, 1:Ne)+(8/3)*OdBlock;
Ly(1:Ne, Ne+1:2*Ne) = -2*OdBlock;
Ly(1:Ne, 2*Ne+1:3*Ne) = (1/3)*OdBlock;

% discrete Laplacian in x
OdBlock = 0*OdBlock;
MBlock = spdiags([e -2*e e], [-1:1], Ne, Ne );
MBlock(1,Ne) = 1;
MBlock(Ne,1) = 1;
Lx = blktridiag(MBlock, OdBlock, OdBlock, Ne);

L = Lx/(hx*hx) + Ly/(hy*hy);