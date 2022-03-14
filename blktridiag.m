function [A] = blktridiag( Amd, Asub, Asup, n )
% Computes a sparse (block) tridiagonal matrix with n blocks
%
% usage1: identical blocks
%   A = blktridiag( Amd, Asub, Asup, n ) 
%
% usage2: a list of distinct blocks
%   A = blktridiag( Amd, Asub, Asup )
%
%  Parameters (usage1):
%       Amd  = a pxq array, for the main diagonal blocks
%       Asub = a pxq array, for sub diagonal block
%         Asub must be the same size and shape as Amd
%       Asup = a pxq array, super diagonal block
%         Asup must be the same size and shape as Amd
%       n    = a scalar integer, defines the number of blocks
%         When n == 1, only a single block will be formed, A == Amd
%
%  Parameters (usage2):
%       Amd = a pxqxn array, a list of n distinct pxq arrays
%         Each plane of Amd corresponds to a single block
%         on the main diagonal.
%       Asub = a pxqx(n-1) array, a list of n-1 distinct pxq arrays,
%         Each plane of Asub corresponds to a single block
%         on the sub-diagonal.
%       Asup = a pxqx(n-1) array, a list of n-1 distinct pxq arrays,
%         Each plane of Asup corresponds to a single block
%         on the super-diagonal.
%
%   NOTE: the sizes of Amd, Asub, and Asup must be consistent
%       with each other, or an error will be generated.
% 
%  Return:
%       A = (n*p by n*q) SPARSE block tridiagonal array
%
%  Created on 26 May 2020
%          by Ying Zhang (phzhang@bu.edu)
%
%


% Which mode of operation are we in?
if( nargin == 4 )
  % replicated block mode
  
  % verify the inputs in this mode are 2-d arrays.
  if( length(size(Amd) )~=2 ) || ...
     ( length(size(Asub) )~=2 ) || ...
     ( length(size(Asup) )~=2 ) 
    error 'Inputs must be 2d arrays if a replication factor is provided'
  end
  
  % get block sizes, check for consistency
  [p,q] = size(Amd);
  if( isempty(Amd) )
    error 'Blocks must be non-empty arrays or scalars'
  end
  if( any(size(Amd)~=size(Asub)) || any(size(Amd)~=size(Asup)) )
    error 'Amd, Asub, Asup are not identical in size'
  end

  if( isempty(n) || (length(n)>1) || (n<1) || (n~=floor(n)) )
    error 'n must be a positive scalar integer'
  end
  
  % scalar inputs?
  % since p and q are integers...
  if( (p*q) == 1 )
    if( n == 1 )
      A = Amd;
    else
      % faster as Jos points out
      A = spdiags(repmat([Asub Amd Asup],n,1),-1:1,n,n);
    end
    % no need to go any farther
    return
  end
  
  % use sparse. the main diagonal elements of each array are...
  v = repmat(Amd(:),n,1);
  % then the sub and super diagonal blocks.
  if( n > 1 )
    % sub-diagonal
    v = [v;repmat(Asub(:),n-1,1)];
    
    % super-diagonal
    v = [v;repmat(Asup(:),n-1,1)];
  end
  
elseif( nargin == 3 )
  % non-replicated blocks, supplied as planes of a 3-d array
  
  % get block sizes, check for consistency
  [p,q,n] = size(Amd);
  if( isempty(Amd) )
    error 'Blocks must be (non-empty) arrays or scalars'
  end
  
  if( (p~=size(Asub,1)) || ...
      (q~=size(Asub,2)) || ...
      (p~=size(Asup,1)) || ...
      (q~=size(Asup,2)) )
    error 'Amd, Asub, Asup do not have the same size blocks'
  end

  if( (n>1) && (((n-1) ~= size(Asub,3)) || ((n-1) ~= size(Asup,3))) )
    error 'Asub and Asup must each have one less block than Amd'
  end
  
  % scalar inputs?
  if( (p*q) == 1 )
    if( n == 1 )
      A = Amd(1);
    else
      % best to just use spdiags
      A = spdiags([[Asub(:);0], Amd(:), [0;Asup(:)]], -1:1, n, n);
    end
    % no need to go any farther
    return
  end
  
  % The main diagonal elements
  v = Amd(:);
  % then the sub and super diagonal blocks.
  if( n > 1 )
    % sub-diagonal
    v = [v; Asub(:)];

    % super-diagonal
    v = [v; Asup(:)];
  end
else
  % must have 3 or 4 arguments
  error 'Must have either 3 or 4 arguments to BLKTRIDIAG'
end

% now generate the index arrays
% first the main diagonal
[ind1,ind2,ind3] = ndgrid(0:p-1,0:q-1,0:n-1);
rind = 1 + ind1(:) + p*ind3(:);
cind = 1 + ind2(:) + q*ind3(:);

% then the sub and super diagonal blocks.
if( n > 1 )
  % sub-diagonal
  [ind1,ind2,ind3] = ndgrid(0:p-1, 0:q-1, 0:n-2);
  rind = [rind; 1+p+ind1(:)+p*ind3(:)];
  cind = [cind; 1+ind2(:)+q*ind3(:)];

  % super-diagonal
  rind = [rind; 1+ind1(:)+p*ind3(:)];
  cind = [cind; 1+q+ind2(:)+q*ind3(:)];
end

% build the final array all in one call to sparse
A = sparse(rind,cind,v,n*p,n*q);
