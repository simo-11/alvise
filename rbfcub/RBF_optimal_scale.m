
function RBF_scale=RBF_optimal_scale(centers,f_centers,RBF_type,...
    minep,maxep)

%--------------------------------------------------------------------------
% Object:
%--------------------------------------------------------------------------
% As explained in [1], we determine an almost optimal parameter "RBF_scale"
% for RBF interpolation by means of a LOOCV type procedure.
%--------------------------------------------------------------------------
% Input:
%--------------------------------------------------------------------------
%
% centers: RBF centers (N x 2 matrix, where N is the number of centers);
%
% f_centers: values of the function at centers;
%
% RBF_type:
%      1: phi=@(r) (1+r.*r).^(1/2);             % Multiquadric
%      2: phi=@(r) exp(-r.*r);                  % Gaussian
%      3: phi=@(r) (1+r.*r).^(-1/2);            % Inverse Multiquadric
%      4: phi=@(r) (1+4*r).*(max(0,(1-r))).^4;  % Wendland 2
%      5: phi=@(r) r.*r.*log(r+eps*(1-sign(r))); % TPS
%      6: phi=@(r) r.^3;                        % polyharmonic spline
%      7: phi=@(r) r.^5;                        % polyharmonic spline
%      8: phi=@(r) r.^7;                        % polyharmonic spline
%      9: phi=@(r) (max(0,(1-r))).^2;             % Wendland W0
%      10: phi=@(r) (35*r.^2+18*r+3).*(max(0,(1-r))).^6;         % Wendland W4
%      11: phi=@(r) (32*r.^3+25*r.^2+8*r+1).*(max(0,(1-r))).^8;  % Wendland W6
%      12: phi=@(r) (sqrt(2)/(3*sqrt(pi))*(3*r^2*log(r/(1+sqrt(1-r.^2)))+...
%            (2*r^2+1).*sqrt(1-r^2)))*max(0,(1-r));  % Missing Wendland
%      13: phi=@(r) exp(-r);                     % Matern beta_1=3/2.
%      14: phi=@(r) (1+r).*exp(-r);              % Matern beta_2=5/2.
%
%      Note: Optimal scaling "RBF_scale" is sought for RBF_type not in
%            [5,8], otherwise we set "RBF_scale" to 1 (those RBF are scale 
%            independent).
%
% minep,maxep: define range for searching a (near) optimal RBF shape
%       parameter, e.g. epsilon in [minep,maxep]=[0.5,15]
%--------------------------------------------------------------------------
% Output:
%--------------------------------------------------------------------------
%
% RBF_scale: almost optimal RBF shape parameter.
%
%--------------------------------------------------------------------------
%% Copyright (C) 2021- R. Cavoretto, A. De Rossi, A. Sommariva, M. Vianello
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%%
%% Authors:
%%          Roberto Cavoretto <roberto.cavoretto@unito.it>
%%          Alessandra De Rossi   <?alessandra.derossi@unito.it>
%%          Alvise Sommariva <alvise@euler.math.unipd.it>
%%          Marco Vianello   <marcov@euler.math.unipd.it>
%%
%% Last Update: July 9, 2021.
%--------------------------------------------------------------------------

% defaults
if nargin < 4, minep=.5; end
if nargin < 5, maxep=15; end

Npts=size(centers,1);

if (RBF_type >= 5 & RBF_type <= 8)
    RBF_scale=ones(Ncenters,1);
else
    % .... define loocv ....
    switch RBF_type
        case 1, RBF_order=1;
        case 5, RBF_order=2;
        case 6, RBF_order=2;
        case 7, RBF_order=3;
        case 8, RBF_order=4;
        otherwise, RBF_order=0;
    end
    
    % ... adding polynomial components (if required) ...
    switch RBF_order
        case 1
            border_RBF_cubmat=ones(Npts,1);
        case 2
            border_RBF_cubmat=[ones(Npts,1)  centers];
        case 3
            c1=(centers(:,1)); c2=(centers(:,2));
            border_RBF_cubmat=[ones(Npts,1) centers c1.^2 c1.*c2 c2.^2];
        case 4
            c1=(centers(:,1)); c2=(centers(:,2));
            border_RBF_cubmat=[ones(Npts,1) centers c1.^2 c1.*c2 c2.^2 ...
                c1.^3 (c1.^2).*c2 c1.*(c2.^2) c2.^3];
        otherwise
            border_RBF_cubmat=[];
    end
    
    DM = DistanceMatrix(centers,centers);
    phi = RBF(RBF_type);
    
    optep = fminbnd(@(epx) CostEpsilonLOOCV(epx,...
        DM,phi,f_centers,border_RBF_cubmat,RBF_order),...
        minep,maxep);
    RBF_scale=1/optep*ones(Npts,1);
end











function DM = DistanceMatrix(dsites,ctrs)

%--------------------------------------------------------------------------
% Object:
%--------------------------------------------------------------------------
% DM = DistanceMatrix(dsites,ctrs)
% Forms the distance matrix of two sets of points in R^s,
% i.e., DM(i,j) = || datasite_i - center_j ||_2.
%--------------------------------------------------------------------------
% Input:
%--------------------------------------------------------------------------
%   dsites: Mxs matrix representing a set of M data sites in R^s
%              (i.e., each row contains one s-dimensional point)
%   ctrs:   Nxs matrix representing a set of N centers in R^s
%              (one center per row)
%--------------------------------------------------------------------------
% Output:
%--------------------------------------------------------------------------
%   DM:     MxN matrix whose i,j position contains the Euclidean
%              distance between the i-th data site and j-th center
%--------------------------------------------------------------------------
% Last Update: July 9, 2021.
%--------------------------------------------------------------------------

[M,~] = size(dsites); [N,s] = size(ctrs);
DM = zeros(M,N);

%--------------------------------------------------------------------------
% Accumulate sum of squares of coordinate differences
% The ndgrid command produces two MxN matrices:
%
%   dr, consisting of N identical columns (each containing
%       the d-th coordinate of the M data sites);
%
%   cc, consisting of M identical rows (each containing
%       the d-th coordinate of the N centers).
%--------------------------------------------------------------------------

for d=1:s
    [dr,cc] = ndgrid(dsites(:,d),ctrs(:,d));
    DM = DM + (dr-cc).^2;
end
DM = sqrt(DM);






function ceps = CostEpsilonLOOCV(ep,r,phi,rhs,P,RBF_order)

%--------------------------------------------------------------------------
% Object:
%--------------------------------------------------------------------------
% Implements cost function for optimization of shape parameter via Rippa's
% LOOCV algorithm-
%--------------------------------------------------------------------------
% Last Update: July 9, 2021.
%--------------------------------------------------------------------------

A = phi(r.*ep);
if RBF_order >= 1
    zeromat = zeros(size(P,2));
    A = [A P; ...
        P' zeromat];
    rhs = [rhs; zeros(size(zeromat,1),1)];
end
invA = pinv(A);
EF = (invA*rhs)./diag(invA);
EF = EF(1:end-RBF_order);
ceps = norm(EF(:),inf);

