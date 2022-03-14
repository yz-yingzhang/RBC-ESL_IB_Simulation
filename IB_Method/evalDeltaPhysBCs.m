function[idxs, delta, xIdx, yIdx] = evalDeltaPhysBCs( Xl, Ne, hx, hy )
% Calculates the delta_h function at Eulerian grid locations idxs
%   using Donev's version
%
%  [idxs, delta] = evalDeltaPhysBCs( Xl, Ne, hx, hy )
%
%  Parameters:
%       Xl    = input Lagrangian marker position
%       Ne    = number of Eulerian grid points in each direction
%       hx    = Eulerian mesh spacing in x direction
%       hy    = Eulerian mesh spacing in y direction
%
%  Return:
%       idxs  = Nl by 16 matrix, where Nl is the number of Langrangian
%           markers, that gives the indices in the Euclidean mesh of where
%           the weight should be applied
%       delta = the weights (from evaluating delta_h) to apply at each
%           index above
%
%
%  Created on 17 July 2020
%          by Ying Zhang (yingzhang@brandeis.edu)
%
%

xL = -2:1;
xR = -1:2;
yD = -2:1;
yU = -1:2;

XlLen = length(Xl(:,1));
idxs  = zeros(XlLen, 16);
delta = zeros(XlLen, 16);

for( i = 1:XlLen )
    
    idxY = [];
    
    fLidx = [(Xl(i,1)/hx) (Xl(i,2)/hy)];
    Lidx = round(fLidx);
    
    if( Lidx(1) > fLidx(1) )
        xIdx = xL + Lidx(1);
    else
        xIdx = xR + Lidx(1);
    end
    
    if( Lidx(2) > fLidx(2) )
        yIdx = yD + Lidx(2);
    else
        yIdx = yU + Lidx(2);
    end
    
    % list the indices of total 16 pairs of Eulerian points
    [A,B] = meshgrid(xIdx,yIdx);
    totIdxVec = reshape(cat(2, A', B'), [], 2);
    
    % evaluate delta_h
    evalPhix = evalPhi(fLidx(1)-totIdxVec(:,1));
    evalPhiy = evalPhi(fLidx(2)-totIdxVec(:,2));
    
    % enforce no-slip BCs in y by flipping sign
    % i.e. u_{i,-1} = -u_{i,1} and v_{i,-1} = v_{i,1}
    idxY = find(yIdx < 0);
    for( k = 1:length(idxY) )
        evalPhiy((4*(idxY(k)-1)+1):(4*(idxY(k)-1)+1)+4-1) = -evalPhiy((4*(idxY(k)-1)+1):(4*(idxY(k)-1)+1)+4-1);
    end
    yIdx = abs(yIdx);
    delta(i,:) = (evalPhix .* evalPhiy) ./ (hx*hy);
    
    % update indices in x for periodic BCs
    xIdx = mod(xIdx, Ne) + 1;
    
    idxs(i,:) = [ (xIdx + Ne * yIdx(1)) (xIdx + Ne * yIdx(2)) (xIdx + Ne * yIdx(3)) (xIdx + Ne * yIdx(4))];
    
end
