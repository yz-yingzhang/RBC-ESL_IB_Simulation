function [G0x, G0y] = D02DNeumann(Ne, hx, hy)
%  Computes the 2D centered difference operators on with zero Neumann BCs
%
%  [G0x, G0y] = D02DNeumann(N, h);
%     
%  Parameters:
%       Ne    = the number of mesh points in each direction
%       h     = the mesh width
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
MBlock = spdiags([-e e], [-1:0], Ne, Ne );
MBlock(1,Ne) = -1;
G0x = blktridiag(MBlock, MBlock, OdBlock, Ne);
G0x(1:Ne, 1:Ne) = G0x(1:Ne, 1:Ne)*2;

% gradient in y direction
MBlock = spdiags([e e], [-1:0], Ne, Ne );
MBlock(1,Ne) = 1;
G0y = blktridiag(MBlock, -MBlock, OdBlock, Ne);
G0y(1:Ne, (M-Ne+1):M) = -MBlock;
G0y(1:Ne,:)  = 0;

G0x = G0x / (2*hx);
G0y = G0y / (2*hy);