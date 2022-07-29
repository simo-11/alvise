
function [tw,xw]=trigauss(n,alpha,beta,method)

%--------------------------------------------------------------------------
% AUTHORS.
%--------------------------------------------------------------------------
% Alvise Sommariva and Marco Vianello, University of Padova
% July 25, 2016.
%
% previous versions with the help of G. Da Fies.
%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% It computes the angles and weights of a trigonometric quadrature formula
% on [alpha,beta], 0<beta-alpha<=pi, matching the moments up to 10^(-14).
%
% Depending on the choice of the variable 'method', some methods are
% implemented.
%
% The formula integrates the canonical trigonometric basis with accuracy
% from about 10^(-15) (small omega) to about 10^(-13) (omega-->pi)
% up to n=300.
%--------------------------------------------------------------------------
% INPUT.
%--------------------------------------------------------------------------
% n: trigonometric degree of exactness.
%
% [alpha,beta]: angular interval, 0<(beta-alpha)/2<=pi.
%
% method:
%         1. 'classic' implements classical trigauss rule.
%         2. 'legendre' implements (shifted) Gauss-Legendre rule
%             matching the trigonometric moments up to 10^(-14). Few nodes
%             if angles are small otherwise the cardinality might be higher
%             than in 'classic'.
%         3. 'better': chooses 'classic' or 'legendre', so to have the
%             minimal number of points, still matching the moments up to
%             10^(-14).
%         4.  'antigauss': antigaussian formula based on Legendre rule of
%             degree M=M(n).
%             It provides a (M+1) x 4 matrix whose entries are
%                             tw=[tAGL wAGL tGLe wGLe]
%             where
%             a) (tAGL,wAGL) is the antigaussian rule;
%             b) (tGLe(1:end-1),wGLe(1:end-1)) is the gaussian rule;
%         5.  'kronrod': Gauss-Kronrod formula based on Legendre rule of
%             degree M=M(n).
%             It provides a (2*M+1) x 3 matrix whose entries are
%                             tw=[tGK wGK wGLf]
%             where
%             a) (tGK,wGK) is the Gauss-Kronrod rule;
%             b) (tGK(2:2:end),wGLf(2:2:end)) is the gaussian rule;
%         6.  'antitrigauss': antigaussian formula based on trigonometric
%             gaussian rule of degree n (whose cardinality is n+1).
%             It provides a (n+2) x 4 matrix whose entries are
%                             tw=[tAGL wAGL tGLe wGLe];
%             where
%             a) (tAGL,wAGL) is the antigaussian rule;
%             b) (tGLe(1:end-1),wGLe(1:end-1)) is the gaussian rule;
%         7.  'trig_kronrod': Gauss-Kronrod based on trigonometric gaussian
%             rule of degree n.
%             It provides a (2*(n+1)+1) x 3 matrix whose entries are
%                             tw=[tTK wTK wTf];
%             where
%             a) (tTK,wTK) is the Gauss-Kronrod rule;
%             b) (tTK(2:2:end),wTf(2:2:end)) is the trig. gaussian rule;
%         8.  'classic_fast': classic version of the trigauss (see paper 
%              Fast variants of the Golub and Welsch algorithm for symmetric
%              weight functions by G. Meurant and A. Sommariva (2012)). 
%         otherwise: if no string is given it chooses 'better' option by
%             default.
%--------------------------------------------------------------------------
% OUTPUT.
%--------------------------------------------------------------------------
% tw: 1) for 'classic', 'legendre', 'better', it is a
%     M x 2 array of (angles,weights) (M depends on the choosen rule)
%     2) for 'antigauss' and 'antitrigauss' it is a M x 4 matrix (see
%     explanation above at points 4. and 6.)
%     3) for 'kronrod' and 'trig_kronrod' it is a M x 3 matrix (see
%     explanation above at points 5. and 7.)
% xw: it is the M x 2 array of the nodes and weights array with respect to
%     the associated weight function of the form 1/sqrt{1-a^2 x^2} in the
%     interval [-1,1]. This matrix is non null only for the 'classic' and 
%     'antigaussian' preferences.
%--------------------------------------------------------------------------
% EXAMPLES.
%--------------------------------------------------------------------------
%
% >> format long e;
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'better');
% >> twB
%
% twB =
%
%     -1.508419779244729e-02     1.590094129717467e-03
%     -1.251400776401954e-02     3.493153120682099e-03
%     -8.255043791082401e-03     4.927692470361337e-03
%     -2.881384626391015e-03     5.697023547188071e-03
%      2.881384626391019e-03     5.697023547188069e-03
%      8.255043791082396e-03     4.927692470361324e-03
%      1.251400776401955e-02     3.493153120682099e-03
%      1.508419779244729e-02     1.590094129717463e-03
%
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'classic');
% >> twB
%
% twB =
%
%     -1.536597391335750e-02     8.744544324142629e-04
%     -1.393392018525980e-02     1.972635825710960e-03
%     -1.146915301262857e-02     2.926255469814080e-03
%     -8.153889677130131e-03     3.662992813352760e-03
%     -4.233938891704753e-03     4.128095270647346e-03
%      2.178667008161081e-18     4.287058912019128e-03
%      4.233938891704759e-03     4.128095270647326e-03
%      8.153889677130132e-03     3.662992813352761e-03
%      1.146915301262857e-02     2.926255469814094e-03
%      1.393392018525980e-02     1.972635825710951e-03
%      1.536597391335750e-02     8.744544324142609e-04
%
% >> format short e
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'antigauss');
% >> twB
%
% twB =
%
%   -1.5612e-02   5.4060e-04  -1.5084e-02   1.5901e-03
%   -1.4036e-02   2.5826e-03  -1.2514e-02   3.4932e-03
%   -1.0564e-02   4.2828e-03  -8.2550e-03   4.9277e-03
%   -5.6645e-03   5.4042e-03  -2.8814e-03   5.6970e-03
%   -7.0636e-19   5.7956e-03   2.8814e-03   5.6970e-03
%    5.6645e-03   5.4042e-03   8.2550e-03   4.9277e-03
%    1.0564e-02   4.2828e-03   1.2514e-02   3.4932e-03
%    1.4036e-02   2.5826e-03   1.5084e-02   1.5901e-03
%    1.5612e-02   5.4060e-04            0            0
%
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'kronrod');
% >> twB
%
% twB =
%
%   -1.5604e-02   2.7995e-04            0
%   -1.5084e-02   7.7659e-04   1.5901e-03
%   -1.4045e-02   1.2956e-03            0
%   -1.2514e-02   1.7537e-03   3.4932e-03
%   -1.0561e-02   2.1404e-03            0
%   -8.2550e-03   2.4607e-03   4.9277e-03
%   -5.6659e-03   2.7029e-03            0
%   -2.8814e-03   2.8494e-03   5.6970e-03
%    3.4551e-18   2.8973e-03            0
%    2.8814e-03   2.8494e-03   5.6970e-03
%    5.6659e-03   2.7029e-03            0
%    8.2550e-03   2.4607e-03   4.9277e-03
%    1.0561e-02   2.1404e-03            0
%    1.2514e-02   1.7537e-03   3.4932e-03
%    1.4045e-02   1.2956e-03            0
%    1.5084e-02   7.7659e-04   1.5901e-03
%    1.5604e-02   2.7995e-04            0
%
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'antitrigauss');
% >> twB
%
% twB =
%
%   -1.5366e-02   8.7445e-04  -1.5298e-02   1.0473e-03
%   -1.3934e-02   1.9726e-03  -1.3588e-02   2.3476e-03
%   -1.1469e-02   2.9263e-03  -1.0672e-02   3.4414e-03
%   -8.1539e-03   3.6630e-03  -6.8077e-03   4.2296e-03
%   -4.2339e-03   4.1281e-03  -2.3385e-03   4.6420e-03
%    2.1787e-18   4.2871e-03   2.3385e-03   4.6420e-03
%    4.2339e-03   4.1281e-03   6.8077e-03   4.2296e-03
%    8.1539e-03   3.6630e-03   1.0672e-02   3.4414e-03
%    1.1469e-02   2.9263e-03   1.3588e-02   2.3476e-03
%    1.3934e-02   1.9726e-03   1.5298e-02   1.0473e-03
%    1.5366e-02   8.7445e-04            0            0
%
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'trig_kronrod');
% >> twB
%
% twB =
%
%   -1.5640e-02   1.8370e-04            0
%   -1.5298e-02   5.1143e-04   1.0473e-03
%   -1.4611e-02   8.6012e-04            0
%   -1.3588e-02   1.1787e-03   2.3476e-03
%   -1.2265e-02   1.4628e-03            0
%   -1.0672e-02   1.7183e-03   3.4414e-03
%   -8.8397e-03   1.9398e-03            0
%   -6.8077e-03   2.1160e-03   4.2296e-03
%   -4.6243e-03   2.2427e-03            0
%   -2.3385e-03   2.3207e-03   4.6420e-03
%   -4.7437e-18   2.3475e-03            0
%    2.3385e-03   2.3207e-03   4.6420e-03
%    4.6243e-03   2.2427e-03            0
%    6.8077e-03   2.1160e-03   4.2296e-03
%    8.8397e-03   1.9398e-03            0
%    1.0672e-02   1.7183e-03   3.4414e-03
%    1.2265e-02   1.4628e-03            0
%    1.3588e-02   1.1787e-03   2.3476e-03
%    1.4611e-02   8.6012e-04            0
%    1.5298e-02   5.1143e-04   1.0473e-03
%    1.5640e-02   1.8370e-04            0
%
%--------------------------------------------------------------------------
% NOTE.
%--------------------------------------------------------------------------
% For examples about the usage of formulas of antigaussian or Kronrod type,
% see the file "demo_trigauss_error.m".
%--------------------------------------------------------------------------
%% Copyright (C) 2016
%% Gaspare Da Fies, Alvise Sommariva, Marco Vianello.
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
%% Gaspare Da Fies, Gerard Meurant, Alvise Sommariva, Marco Vianello.
%%
%% Date: JULY 25, 2016
%--------------------------------------------------------------------------

if nargin <= 3
    method = 'better';
end

validStrings = {'classic','legendre','better','antigauss',...
    'kronrod','antitrigauss','trig_kronrod','classic_fast'};

xw=[];

if ( any(strcmpi(method, validStrings)) )
    
    if ( strcmpi(method, 'classic') )
        [tw,xw]=trigauss_classic(n,alpha,beta);
    end
    
    if ( strcmpi(method, 'legendre') )
        if beta-alpha==2*pi
            tw=circle_trigquad(n,alpha,beta);
        else
            NN=degree_trigGL(n,alpha,beta);
            tw=gauss_legendre(NN,alpha,beta);
        end
    end
    
    if ( strcmpi(method, 'better') )
        if beta-alpha==2*pi
            tw=circle_trigquad(n,alpha,beta);
        else
            NN=degree_trigGL(n,alpha,beta);
            if NN <= n
                tw=gauss_legendre(NN,alpha,beta);
            else
                tw=trigauss_classic(n,alpha,beta);
            end
        end
    end
    
    if ( strcmpi(method, 'antigauss') )
        NN=degree_trigGL(n,alpha,beta);
        tw=antigauss_full(NN,alpha,beta);
    end
    
    if ( strcmpi(method, 'kronrod') )
        NN=degree_trigGL(n,alpha,beta);
        tw=gauss_kronrod_full(NN,alpha,beta);
    end
    
    if ( strcmpi(method, 'antitrigauss') )
        tw=antitrigauss_full(n,alpha,beta);
    end
    
    if ( strcmpi(method, 'trig_kronrod') )
        tw=trigauss_kronrod_full(n,alpha,beta);
    end
    
    if ( strcmpi(method, 'classic_fast') )
        omega=(beta-alpha)/2;
        tw=fast_trigauss(n+1,omega);
        tw(:,1)=((alpha+beta)/2)+tw(:,1);
    end
else
    warning('wrong string as method, choosen classical trigauss');
    tw=trigauss(n,alpha,beta,'classic');
end





function NN=degree_trigGL(n,alpha,beta)

% this routine computes the number of nodes that Gauss-Legendre needs so to
% match the trigonometric moments up to 10^(-14).

theta=(beta-alpha)/2;

if n < 500
    
    u=[1 (25:25:500)];
    
    L=[8 29    45    60    75    89   103   119   132   146   158   173  ...
        185   199 212   225   240   253   264   279   291];
    
    s=spline(u,L);
    
    NN=ceil(ppval(s,n*theta));
    
else
    
    u=n*theta;
    s=0.54*u+21;
    NN=ceil(s);
    
end



function tw=gauss_kronrod_full(n,alpha,beta)
n0=ceil(3*n/2)+1;
ab0=r_jacobi(n0,0,0);

% KRONROD.
xw=kronrod(n,ab0);
xGK=xw(:,1); wGK=xw(:,2);
tGK=(beta+alpha)/2+(beta-alpha)*xGK/2;
wGK=(beta-alpha)*wGK/2;


% GAUSS.
xwGL=gauss(n,ab0);
xGL=xwGL(:,1); wGL=xwGL(:,2);
tGL=(beta+alpha)/2+(beta-alpha)*xGL/2;

wGLl=(beta-alpha)*wGL/2;
wGLf=zeros(size(wGK));

wGLf(2:2:end)=wGLl;

tw=[tGK wGK wGLf];





function tw=antigauss_full(n,alpha,beta)

ab=r_jacobi(n+1,0,0);
ab(end,2)=2*ab(end,2);
xwAGL=gauss(n+1,ab);

xAGL=xwAGL(:,1); wAGL=xwAGL(:,2);
tAGL=(beta+alpha)/2+(beta-alpha)*xAGL/2;
wAGL=(beta-alpha)*wAGL/2;

xwGL=gauss(n,ab(1:end-1,:));
xGL=xwGL(:,1); wGL=xwGL(:,2);
tGL=(beta+alpha)/2+(beta-alpha)*xGL/2; tGLe=[tGL; 0];
wGL=(beta-alpha)*wGL/2; wGLe=[wGL; 0];

tw=[tAGL wAGL tGLe wGLe];







function tw=trigauss_kronrod_full(n,alpha,beta)

omega=(beta-alpha)/2;
n0=ceil(3*n/2)+1;
ab0=r_trigauss(n0,omega);

% KRONROD.
xw=kronrod(n,ab0);
xTK=xw(:,1); wTK=xw(:,2);
tTK(:,1)=2*asin(sin(omega/2)*xTK(:,1))+(beta+alpha)/2;


% GAUSS.
xwT=gauss(n,ab0);
xT=xwT(:,1); wT=xwT(:,2);
tT(:,1)=2*asin(sin(omega/2)*xT(:,1))+(beta+alpha)/2;

wTf=zeros(size(wTK));

wTf(2:2:end)=wT;

tw=[tTK wTK wTf];





function tw=antitrigauss_full(n,alpha,beta)

omega=(beta-alpha)/2;

ab=r_trigauss(n+1,omega);
ab(end,2)=2*ab(end,2);
xwAT=gauss(n+1,ab);

xAT=xwAT(:,1); wAT=xwAT(:,2);
tAT=2*asin(sin(omega/2)*xwAT(:,1))+(beta+alpha)/2;

xwT=gauss(n,ab(1:end-1,:));
xT=xwT(:,1); wT=xwT(:,2);
tT=2*asin(sin(omega/2)*xT)+(beta+alpha)/2;
tTe=[tT; 0];
wTe=[wT; 0];

tw=[tAT wAT tTe wTe];




function tw=gauss_legendre(n,alpha,beta)

% [xGL, wGL] = legpts(N); wGL=wGL';
ab=r_jacobi(n,0,0);
xw=gauss(n,ab); xGL=xw(:,1); wGL=xw(:,2);
t=(beta+alpha)/2+(beta-alpha)*xGL/2;
w=(beta-alpha)*wGL/2;
tw=[t w];




function [tw,xw]=trigauss_classic(n,alpha,beta)
% by Gaspare Da Fies and Marco Vianello, University of Padova
% 8 Nov 2011

% computes the n+1 angles and weights of a trigonometric gaussian
% quadrature formula on [alpha,beta], 0<beta-alpha<=pi

% uses the routines chebyshev.m, gauss.m from
% www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
% we suggest to put the following statements
% ab = zeros(N,2); sig = zeros(N+1,2*N);
% at the beginning of the body of chebyshev.m to speed-up execution

% input:
% n: trigonometric degree of exactness
% [alpha,beta]: angular interval, 0<beta-alpha<=pi

% output:
% tw: (n+1) x 2 array of (angles,weights)

% the formula integrates the canonical trigonometric basis with accuracy
% from about 10^(-15) (small omega) to about 10^(-13) (omega-->pi)
% up to n=300


% half-length of the angular interval
omega=(beta-alpha)/2;

if omega == pi
    tw=circle_trigquad(n,0,2*pi);
else
    
    ab=r_trigauss(n,omega);
    % format long e; ab
    
    % Gaussian formula for the weight function above
    xw=gauss(n+1,ab);
    
    % angles and weights for the trigonometric gaussian formula
    tw(:,1)=2*asin(sin(omega/2)*xw(:,1))+(beta+alpha)/2;
    tw(:,2)=xw(:,2);
    
end



function ab=r_trigauss(n,omega)

% modified Chebyshev moments by recurrence
z(1)=2*omega;

z(n+1)=quadgk(@(t)cos(2*n*acos(sin(t/2)/sin(omega/2))),...
    -omega,omega,'MaxIntervalCount',50000);
temp=(2:2:2*n-1);
dl=1/4-1./(4*(temp-1));
dc=1/2-1/sin(omega/2)^2-1./(2*(temp.^2-1));
du=1/4+1./(4*(temp+1));
d=4*cos(omega/2)/sin(omega/2)./(temp.^2-1)';
d(n-1)=d(n-1)-du(n-1)*z(n+1);
z(2:n)=tridisolve(dl(2:n-1),dc(1:n-1),du(1:n-2),d(1:n-1));
mom=zeros(1,2*n+2);
mom(1:2:2*n+1)=z(1:n+1);

% normalization of the moments (monic polynomials)
k=(3:length(mom));
mom(3:end)=exp((2-k)*log(2)).*mom(3:end);

% recurrence coeffs of the monic Chebyshev polynomials
abm(:,1)=zeros(2*n+1,1);
abm(:,2)=0.25*ones(2*n+1,1); abm(1,2)=pi; abm(2,2)=0.5;

% recurrence coeffs for the monic OPS w.r.t. the weight function
% w(x)=2*sin(omega/2)/sqrt(1-sin^2(omega/2)*x^2)
% by the modified Chebyshev algorithm
[ab,normsq]=chebyshev(n+1,mom,abm);






function x = tridisolve(a,b,c,d)
%   TRIDISOLVE  Solve tridiagonal system of equations.
% From Cleve Moler's Matlab suite
% http://www.mathworks.it/moler/ncmfilelist.html

%     x = TRIDISOLVE(a,b,c,d) solves the system of linear equations
%     b(1)*x(1) + c(1)*x(2) = d(1),
%     a(j-1)*x(j-1) + b(j)*x(j) + c(j)*x(j+1) = d(j), j = 2:n-1,
%     a(n-1)*x(n-1) + b(n)*x(n) = d(n).
%
%   The algorithm does not use pivoting, so the results might
%   be inaccurate if abs(b) is much smaller than abs(a)+abs(c).
%   More robust, but slower, alternatives with pivoting are:
%     x = T\d where T = diag(a,-1) + diag(b,0) + diag(c,1)
%     x = S\d where S = spdiags([[a; 0] b [0; c]],[-1 0 1],n,n)

x = d;
n = length(x);
for j = 1:n-1
    mu = a(j)/b(j);
    b(j+1) = b(j+1) - mu*c(j);
    x(j+1) = x(j+1) - mu*x(j);
end
x(n) = x(n)/b(n);
for j = n-1:-1:1
    x(j) = (x(j)-c(j)*x(j+1))/b(j);
end






function tw=circle_trigquad(n,alpha,beta)

if nargin == 1
    alpha=0;
    beta=2*pi;
end
N=n+1;
w=(2*pi/N)*ones(N,1);
t=linspace(pi/N,2*pi-pi/N,N); t=alpha+t';

tw=[t w];




























function [tw,xw]=anti_trigauss(n,alpha,beta)
% by Gaspare Da Fies and Marco Vianello, University of Padova
% 8 Nov 2011

% computes the n+1 angles and weights of a trigonometric gaussian
% quadrature formula on [alpha,beta], 0<beta-alpha<=pi

% uses the routines chebyshev.m, gauss.m from
% www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
% we suggest to put the following statements
% ab = zeros(N,2); sig = zeros(N+1,2*N);
% at the beginning of the body of chebyshev.m to speed-up execution

% input:
% n: trigonometric degree of exactness
% [alpha,beta]: angular interval, 0<beta-alpha<=pi

% output:
% tw: (n+1) x 2 array of (angles,weights)

% the formula integrates the canonical trigonometric basis with accuracy
% from about 10^(-15) (small omega) to about 10^(-13) (omega-->pi)
% up to n=300


% half-length of the angular interval
omega=(beta-alpha)/2;

% We're calculating Antigaussian points, so we want to shift from
% G(n) to AG(n+1)
n = n+1;

ab=r_trigauss(n,omega);

% Only other real modification from original gaussian nodes
ab(end,2) = ab(end,2)*2;

% Gaussian formula for the weight function above
xw=gauss(n+1,ab);

% angles and weights for the trigonometric gaussian formula
tw(:,1)=2*asin(sin(omega/2)*xw(:,1))+(beta+alpha)/2;
tw(:,2)=xw(:,2);


















function xw=fast_trigauss(N,omega)

if rem(N,2) == 0
    ab = r_subchebyshev(N,omega);
else
    ab = r_subchebyshev(N+1,omega);
end
[xw_symm_AS,nloopS] = SymmAGWo(N,ab);
xw=trigauss_conversion(xw_symm_AS,omega);




function [xw,nloop]=radau_gw(N,ab,end0)
%RADAUGW Gauss-Radau quadrature rule
% use Golub and Welsch instead of gauss
%

% Input
% N : cardinality of the rule
% ab: 3-term recurrence for the orthogonal polynomials
% same as in OPQ
% ab(1,2) is the 0th moment
% end0: prescribed node

% Output
% xw : xw(:,1) nodes, xw(:,2) weights of the quadrature rule

%
% Author A. Sommariva
% Adapted from OPQ from W. Gautschi
% June 2012
%

p0 = 0;
p1 = 1;

for n = 1:N
    pm1 = p0;
    p0 = p1;
    p1 = (end0 - ab(n,1)) * p0 - ab(n,2) * pm1;
end

ab(N+1,1) = end0 - ab(N+1,2) * p0  /p1;

a = ones(size(ab,1),1);
[x,w,nloop] = GWo(a,-ab(:,1),ab(:,2),ab(1,2),1);
xw = [x' w'];




function [xw,nloop]=SymmAGWo(N,ab)
%SYMMAGWo computation of the nodes and weights for a symmetric weight
%function
% this version uses the reduced polynomials and the optimized Golub and
% Welsch function
%

%
% see: Fast variants of the Golub and Welsch algorithm for symmetric
% weight functions by G. Meurant and A. Sommariva (2012)

% Input
% N : cardinality of the rule
% ab: 3-term recurrence for the orthogonal polynomials
% same as in OPQ
% ab(1,2) is the 0th moment

% Output
% xw : xw(:,1) nodes, xw(:,2) weights of the quadrature rule
% nloop: number of iterations in QR

%
% Authors G. Meurant and A. Sommariva
% June 2012
%

N0 = size(ab,1);
if N0 < N
    error('SymmAGWo: input array ab is too short')
end

na = norm(ab(:,1));
if na > 0
    error('SymmAGWo: the weight function must be symmetric')
end

% computation of the 3-term recurrence for the reduced weight function

nu = ab(:,2);
nu_odd = nu(1:2:end);
nu_even = nu(2:2:end);
D = [nu(2); nu_odd(2:end) + nu_even(2:end)];
E = [nu(1) / 2; nu_odd(2:end) .* nu_even(1:end-1)];
ab_symm = [D E];

if rem(N,2) == 0
    % N even
    a = ones(size(ab_symm,1),1);
    % use Golub and Welsch
    [x,w,nloop] = GWo(a,-D,E,E(1),1);
    % nodes and weights
    x = sqrt(x');
    w = w';
    xw = [-flipud(x) flipud(w); x w];
    
else
    % N odd
    N_symm = (N - 1) / 2;
    [xw_symm,nloop] = radau_gw(N_symm,ab_symm,0);
    w0 = 2 * xw_symm(1,2);
    % nodes and weights
    x = sqrt(xw_symm(2:end,1));
    w = xw_symm(2:end,2);
    xw = [-flipud(x) flipud(w); 0 w0; x w];
end




function tw=trigauss_conversion(xw,omega)

tw(:,1)=2*asin(sin(omega/2)*xw(:,1));
tw(:,2)=xw(:,2);



function ab=r_subchebyshev(n,omega)

%--------------------------------------------------------------------------
% FAST TRIGAUSS.
% FAST VARIANT OF A CODE BY DA FIES AND VIANELLO.
%--------------------------------------------------------------------------
% INPUTS:
% n    : NUMBER OF POINTS.
% omega: ARC ANGLE.
%
% OUTPUTS:
% ab   : THREE TERMS RECURSION.
%
%--------------------------------------------------------------------------
%% Copyright (C) 2007-2012
%% Gaspare Da Fies, Gerard Meurant, Alvise Sommariva, Marco Vianello.
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
%% Gaspare Da Fies, Gerard Meurant, Alvise Sommariva, Marco Vianello.
%%
%% Date: June 3, 2012
%--------------------------------------------------------------------------

N = n;
n = n - 1;


% modified Chebyshev moments by recurrence

if rem(N,2) == 1
    NN=N+1; nn=n+1;
else
    NN=N; nn=n;
end

mom=fast_moments_computation(omega,2*nn+1);

% recurrence coeffs of the monic Chebyshev polynomials
abm(:,1)=zeros(2*nn+1,1);
abm(:,2)=0.25*ones(2*nn+1,1); abm(1,2)=pi; abm(2,2)=0.5;

% recurrence coeffs for the monic OPS w.r.t. the weight function
% w(x)=2*sin(omega/2)/sqrt(1-sin^2(omega/2)*x^2)
% by the modified Chebyshev algorithm

% ab = chebyshev(NN+1,mom,abm);
ab = fast_chebyshev(NN,mom,abm);







function ab=fast_chebyshev(N,mom,abm);
%SUBP_MOD_CHEBYSHEV Modified Chebyshev algorithm
% this works only for the subperiodic weight function
%
% From Gautschi's code (simplified)
% Mar 2012
%

ab = zeros(N,2);
sig = zeros(N+1,2*N);

ab(1,2) = mom(1);

sig(1,1:2*N) = 0;
sig(2,:) = mom(1:2*N);

for n = 3:N+1
    for m = n-1:2*N-n+2
        sig(n,m) = sig(n-1,m+1) + abm(m,2) * sig(n-1,m-1) - ab(n-2,2) * sig(n-2,m);
    end
    
    ab(n-1,2) = sig(n,n-1) / sig(n-1,n-2);
end





function mom=fast_moments_computation(omega,n)

mom=zeros(1,n+1);
mom(1)=2*omega; % FIRST MOMENT.

if(n>=2)
    
    if(omega<=1/4*pi)
        l=10;
    elseif(omega<=1/2*pi)
        l=20;
    elseif(omega<=3/4*pi)
        l=40;
    else
        if omega == pi
            l=2*ceil(10*pi);
        else
            l=2*ceil(10*pi/(pi-omega));
        end
    end
    
    
    temp=(2:2:n+2*l-2); % AUXILIAR VECTORS.
    temp2=temp.^2-1;
    
    dl=1/4 -1./(4*(temp-1)); % DIAGONALS.
    dc=1/2 -1/sin(omega/2)^2 -1./(2*temp2);
    du=1/4 +1./(4*(temp+1));
    
    d=4*cos(omega/2)/sin(omega/2)./temp2'; % COMPUTING KNOWN TERM.
    d(end)=d(end);                         % PUT LAST MOMENT NULL.
    
    z=tridisolve(dl(2:end),dc,du(1:end-1),d); % SOLVE SYSTEM.
    mom(3:2:n+1)=z(1:floor(n/2)); % SET ODD MOMENTS.
    
end

mom=mom';

normalized = 0;

if normalized == 0
    M=length(mom);
    kk=2.^(-((1:2:M)-2))'; kk(1)=1;
    v=ones(M,1);
    v(1:2:M)=kk;
    mom=v.*mom;
end


function [t,w,nloop]=GWo(a,b,c,muzero,symm)
%GWo vectorized coding of the Golub and Welsch algorithm
%

% Input
% (a,b,c): coefficients of the orthogonal polynomials
%   p(k)=(a(k) x + b(k)) p(k-1) - c(k) p(k-2)
% a is the diagonal of the tridiagonal matrix, b and c the sub and super
% diags
% muzero: 0th moment
% symm: = 1 symmetrize the matrix
%
% In the symmetric case (symm=0), c is not used (given as [])
%
% Output
% t: nodes, w: weights
% nloop: number of QR iterations
%

%
% Author G. Meurant
% optimized version Apr 2012
%

n=length(a);
if symm == 1
    ac = cumprod(a);
    a = -b ./ a;
    b(1:n-1) = sqrt(c(2:n) ./ ac(2:n));
end

b = b(:)';
a = a(:)';

% computation of the norm of the tridiagonal matrix
abb = abs(b(1:n-1));
apb = [abs(a(1:n-1)) + abb, abs(a(n))];
apb(2:n) = apb(2:n) + abb;
normJ = max(apb);

w = zeros(1,n);
t = zeros(1,n);
w(1) = 1;
% relative zero tolerance
epss = eps * normJ;
lambda = normJ;
lambda1 = lambda;
lambda2 = lambda;
rho = lambda;
m = n;

inspect = 1;
nloop = 0;

while inspect == 1
    
    nloop = nloop + 1;
    
    m1 = m-1;
    k = m1;
    i = m1;
    
    if m1 == 0
        t(1) = a(1);
        w(1) = muzero * w(1)^2;
        % sort the abscissas and return
        [t,ind] = sort(t);
        w = w(ind);
        %   nloop
        return
    end
    
    if abs(b(m1)) <= epss
        t(m) = a(m);
        w(m) = muzero * w(m)^2;
        rho = min(lambda1,lambda2);
        m  = m1;
        continue
    end
    
    if abs(b(i)) <= epss
        k = i;
        break
    else
        k = 1;
    end
    
    % find the shift with the eigenvalues of lower 2x2 block
    b2 = b(m1)^2;
    det = sqrt((a(m1) - a(m))^2 + 4 * b2);
    aa = a(m1) + a(m);
    
    if aa > 0
        lambda2=(aa+det)/2;
    else
        lambda2=(aa-det)/2;
    end
    
    lambda1 = (a(m1) * a(m) - b2) / lambda2;
    eigmax = max(lambda1,lambda2);
    
    if abs(eigmax-rho) <= abs(eigmax)/8
        lambda = eigmax;
    end
    rho = eigmax;
    
    % transform block from k to m
    cj = b(k);
    % shift the diagonal
    bk1 = a(k) - lambda;
    
    if k == 1
        % this is handled separately to avoid a test for b(j-1) at each step in
        % the for loop on j
        kdeb = 2;
        r = sqrt(cj^2 + bk1^2);
        st = cj / r;
        ct = bk1 / r;
        aj = a(1);
        f = aj * ct + b(1) * st;
        q = b(1) * ct + a(2) * st;
        a(1) = f * ct + q * st;
        b(1) = f * st - q * ct;
        wj = w(1);
        a(2) = aj + a(2) - a(1);
        w(1) = wj * ct + w(2) * st;
        w(2) = wj * st - w(2) * ct;
        bk1 = b(1);
        cj = b(2) * st;
        b(2) = -b(2) * ct;
    else
        kdeb = k;
    end % if k
    
    for j = kdeb:m1
        
        % compute and apply rotations
        j1 = j + 1;
        r = sqrt(cj^2 + bk1^2);
        r1 = 1 / r;
        st = cj * r1;
        ct = bk1 * r1;
        aj = a(j);
        b(j-1) = r;
        f = aj * ct + b(j) * st;
        q = b(j) * ct + a(j1) * st;
        a(j) = f * ct + q * st;
        b(j) = f * st - q * ct;
        wj = w(j);
        a(j1) = aj + a(j1) - a(j);
        w(j) = wj * ct + w(j1) * st;
        w(j1) = wj * st - w(j1) * ct;
        bk1 = b(j);
        cj = b(j+1) * st;
        b(j1) = -b(j1) * ct;
        
    end % for
    
end % while

