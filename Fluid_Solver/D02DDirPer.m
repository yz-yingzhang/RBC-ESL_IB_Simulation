function [D0x, D0y] = D02DDirPer(Ne, hx, hy)
%  Computes the 2D centered difference operators on a domain
%   zero Dirichlet for y boundaries and periodic for x boundaries
%
%  [D0x, D0y] = D02DDirPer(N, h);
%     
%  Parameters:
%       Ne     = the number of mesh points in each direction
%       hx, hy = the mesh width in x and y
%
%
%  Return:
%       D0x = centered diff. op in the X direction
%       D0y = centered diff. op in the Y direction
%
%  Created on 26 May 2020
%          by Ying Zhang (phzhang@bu.edu)
%
%

M = Ne * Ne;
e = ones(Ne,1);
OdBlock = 0*speye(Ne);

% gradient in x direction
MBlock = spdiags([-e e], [0:1], Ne, Ne );
MBlock(Ne,1) = 1;
D0x = blktridiag(MBlock, OdBlock, MBlock, Ne);

% gradient in y direction
MBlock = spdiags([-e -e], [0:1], Ne, Ne );
MBlock(Ne,1) = -1;
D0y = blktridiag(MBlock, OdBlock, -MBlock, Ne);

D0x = D0x / (2*hx);
D0y = D0y / (2*hy);
