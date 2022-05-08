
function demo_RBF_moms

%--------------------------------------------------------------------------
% * Object:
%--------------------------------------------------------------------------
% In this demo we test the routines computing RBF moments over
% polygons.
%--------------------------------------------------------------------------
% Reference paper:
% A. Sommariva and M. Vianello,
% RBF moment computation and meshless cubature on general polygonal regions.
%--------------------------------------------------------------------------
%% Copyright (C) 2019-2020 Alvise Sommariva, Marco Vianello.
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
%% Author:  Alvise Sommariva <alvise@euler.math.unipd.it>
%%          Marco Vianello   <marcov@euler.math.unipd.it>
%%
%% Date: February 27, 2020.
%% Last update: December 28, 2020.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% * RBF_type:
%--------------------------------------------------------------------------
% 1: phi=@(r) (1+r.*r).^(1/2);             % Multiquadric
% 2: phi=@(r) exp(-r.*r);                  % Gaussian
% 3: phi=@(r) (1+r.*r).^(-1/2);            % Inverse Multiquadric
% 4: phi=@(r) (1+4*r).*(max(0,(1-r))).^4;  % Wendland 2
% 5: phi=@(r) r.*r.*log(r+eps*(1-sign(r))); % TPS
% 6: phi=@(r) r.^3;                        % polyharmonic spline
% 7: phi=@(r) r.^5;                        % polyharmonic spline
% 8: phi=@(r) r.^7;                        % polyharmonic spline
% 9: phi=@(r) (max(0,(1-r))).^2;             % Wendland W0
% 10: phi=@(r) (35*r.^2+18*r+3).*(max(0,(1-r))).^6;         % Wendland W4
% 11: phi=@(r) (32*r.^3+25*r.^2+8*r+1).*(max(0,(1-r))).^8;  % Wendland W6
% 12: phi=@(r) (sqrt(2)/(3*sqrt(pi))*(3*r^2*log(r/(1+sqrt(1-r.^2)))+...
%            (2*r^2+1).*sqrt(1-r^2)))*max(0,(1-r));  % Missing Wendland
% 13 phi=@(r) exp(-r);                     % Matern beta_1=3/2.
% 14 phi=@(r) (1+r).*exp(-r);              % Matern beta_2=5/2.
%
% * Wendland W0 (9), Missing Wendland (12) and Matern (13) are under 
%   revision.
%--------------------------------------------------------------------------
% Note: In the tests of the reference paper: "polygon_typeV=[1:6]".
%--------------------------------------------------------------------------

RBF_typeV=1:14;

%--------------------------------------------------------------------------
% * polygon_typeV:
%--------------------------------------------------------------------------
% -6: polygon with 2 polygonal holes (non simply connected polygonal reg.)
%     and a polygon inside a hole, plus a disconnected polygon.
% -5: polygon with 2 polygonal holes (non simply connected polygonal reg.)
%     and a polygon inside a hole
% -4: polygon with 2 polygonal holes (non simply connected polygonal reg.)
% -3: square with triangular hole (non simply connected polygonal reg.)
% -2: polygonal doughnut (non simply connected polygonal reg.)
% -1: square doughnut (non simply connected polygonal reg.)
%  0: shifted square
%  1: square ([0,1]^2)
%  2: convex polygon (6 sides)
%  3: not convex polygon (9 sides) *
%  4: not convex polygon (9 sides) *
%  5: polygon (not simple domain, candy-like, 10 sides) *
%  6: square ([-1,1]^2)
%  7: simplex ([0 0; 1 0; 1 1])
%  8: polygon approximation of a circle.
%--------------------------------------------------------------------------
% Note: In the tests of the reference paper: "polygon_typeV=[-6 4]".
%--------------------------------------------------------------------------

polygon_typeV=[-6 4];      % see subroutine "define_polygon_polyshape" below
% ("polygon_typeV" may be a vector)
Ncenters=50;          % number of centers (choose less than 100)
RBF_scale=1;
doplot=1;               % domain pointset plot: 0: no 1: yes
strcub=''; % it can be 'triangulation' or 'polar' or '' (empty string).


% ........................... warning .....................................

if Ncenters > 50, warning('Routine can be time consuming'); end

% ...................... main code below ..................................

% .................. perform battery of experiments .......................

aeV=[];
reV=[];
cpus2V=[];
polygonV=[];
rbfV=[];
for jj=1:length(polygon_typeV)
    
    % .... define polygon and its vertices (polyshape object) ....
    polygon_type=polygon_typeV(jj);
    P=define_polygon_polyshape(polygon_type);
    
    % .... reference algebraic rule ....
    deg_pgs=1000; [xyw,xvc,yvc,P,tri]=polygauss_2018(deg_pgs,P);
    
    % .... define centers ....
    centers=generate_centers_polyshape(P,Ncenters);
    
    for ii=1:length(RBF_typeV)
        RBF_type=RBF_typeV(ii);
        
        tic;
        
        if strcmp(strcub,'') |  ~(strcmp(strcub,'triangulation') | ...
                strcmp(strcub,'polar'))
            if  RBF_type <= 3 | RBF_type >= 15
                strcub='triangulation';
            else
                strcub='polar';
            end
        end
        
        switch strcub
            case 'polar'
                moms_rbf=RBF_moms_polar(P,centers,RBF_type,RBF_scale);
            otherwise
                moms_rbf=RBF_moms_tri(P,centers,RBF_type,RBF_scale);
        end
        cpus_moms=toc;
        
        % ....................... reference cubature ......................
        
        phi = RBF(RBF_type);
        
        xx=xyw(:,1); yy=xyw(:,2); ww=xyw(:,3); % ref cubature rule: nodes, weights
        for ii=1:Ncenters
            PP=centers(ii,:);
            rr=sqrt((xx-PP(1)).^2+(yy-PP(2)).^2);
            zz=phi(rr/RBF_scale);
            moms_pgss(ii,1)=ww'*zz;
        end
        
        % format long e; moms_rbf-moms_pgss
        
        [ae,re]=statistics(RBF_type,RBF_scale,polygon_type,Ncenters,...
            moms_rbf,moms_pgss,cpus_moms,doplot,P,centers,strcub);
        aeV=[aeV; max(ae)]; reV=[reV; max(re)];
        cpus2V=[cpus2V;cpus_moms/Ncenters];
        polygonV=[polygonV; polygon_type];
        rbfV=[rbfV; RBF_type];
    end
end








% ........................ summary of experiments .........................

fprintf('\n \n \n \t ................... summary .............................');
fprintf('\n \t test | abs.err. | rel.err. | moms cpu | polyg.| RBF |');
fprintf('\n \t ......................................................... ');
for k=1:length(aeV)
    aeL=aeV(k); reL=reV(k); cpus2L=cpus2V(k);
    fprintf('\n \t %3.0f  |  %1.1e |  %1.1e |  %1.1e |   %2.0f  |  %2.0f |',...
        k,aeL,reL,cpus2L,polygonV(k),rbfV(k));
end
fprintf('\n \t ......................................................... \n \n');







function [ae,re]=statistics(RBF_type,RBF_scale,polygon_type,Ncenters,...
    moms_rbf,moms_pgss,cpus_moms,doplot,P,centers,strcub)


%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% Single moments cubature experiment, on a certain
% * polygonal domain, represented as polyshape object;
% * by a RBF.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% MAIN PROGRAM STARTS HERE.
%--------------------------------------------------------------------------


fprintf('\n \n \n \t ........ MOMENT COMPUTATION STATISTICS  ........ ');

% ... RBF ...

switch RBF_type
    case 1, RBF_string=' MULTIQUADRIC: (1+r^2)^(1/2)';
    case 2, RBF_string=' GAUSSIAN: exp(-r^2)';
    case 3, RBF_string=' INVERSE MULTIQUADRIC: (1+r^2)^(-1/2)';
    case 4, RBF_string=' WENDLAND, W2: (1+4*r)*(max(0,(1-r))^4)';
    case 5, RBF_string=' TPS: r^2*log(r))';
    case 6, RBF_string=' POLYHARMONIC: r^3';
    case 7, RBF_string=' POLYHARMONIC: r^5';
    case 8, RBF_string=' POLYHARMONIC: r^7';
    case 9, RBF_string=' WENDLAND, W0: max(0,(1-r)^2)';
    case 10, RBF_string=' WENDLAND, W4:(35*r^2+18*r+3)*max(0,(1-r)^6)';
    case 11, RBF_string=' WENDLAND, W6:(32*r^3+25*r^2+8*r+1)*max(0,(1-r)^8)';
    case 12, RBF_string=' MISSING WENDLAND: (32*r^3+25*r^2+8*r+1)*max(0,(1-r)^8)';
    case 13, RBF_string=' MATERN B1=(d+1)/2: exp(-r)';
    case 14, RBF_string=' MATERN B2=(d+3)/2: (1+r)*exp(-r)';
end
fprintf('\n \t RBF TYPE     : %2.0f',RBF_type);
fprintf('\n \t METHOD     : '); disp(strcub);
fprintf('\n \t *** '); disp(RBF_string);

% ... RBF scale ...

fprintf('\n \t RBF SCALE    : % 1.5e',RBF_scale);



% ... experiment data ...

fprintf('\n \t POLYGON TYPE : %2.0f',polygon_type);
fprintf('\n \t NODES CARDIN.: %4.0f',Ncenters);


% ... cubature errors ...

ae=abs(moms_rbf-moms_pgss);
re=abs(moms_rbf-moms_pgss)./(abs(moms_pgss)+...
    (1-sign(abs(moms_pgss)))*10^(-15));

fprintf('\n \t ** MAX.ABSOLUTE ERR.: %1.3e',max(ae));
fprintf('\n \t ** MAX.RELATIVE ERR.: %1.3e',max(re));
fprintf('\n \t ** AVER.CPU: %1.3e',cpus_moms/Ncenters);

% ............................... plots ...................................

if doplot == 1
    clf;
    hold on;
    axis square;
    % plot(xv,yv,'k-','LineWidth',4);
    plot(P);
    plot(centers(:,1),centers(:,2),'r.');
    hold off;
end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% MAIN PROGRAM ENDS HERE.
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





function P=define_polygon_polyshape(polygon_type)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% Define polygon through its vertices.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% polygon_type: parameter that determines the polygon.
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% polygon_sides: polygon described as a polyshape object.
%--------------------------------------------------------------------------


switch polygon_type
    case -1
        
        % outer boundary as polyshape object
        xv1=[0 3 3 0]; yv1=[0 0 3 3]; P1=polyshape(xv1,yv1);
        
        % inner boundary (hole) as polyshape object
        xv2=[1 2 2 1]; yv2=[1 1 2 2]; P2=polyshape(xv2,yv2);
        
        % domain as difference of polyshape objects
        P=subtract(P1,P2);
        
    case -2
        
        N=11; a=0; b=2*pi;
        h=(b-a)/N; theta=(a:h:b-h);
        x=cos(theta); y=sin(theta);
        
        % outer boundary as polyshape object
        xv1=0.5*x+0.5; yv1=0.5*y+0.5; P1=polyshape(xv1,yv1);
        
        % inner boundary (hole) as polyshape object
        xv2=0.15*x+0.5; yv2=0.15*y+0.5; P2=polyshape(xv2,yv2);
        
        % domain as difference of polyshape objects
        P=subtract(P1,P2);
        
        
    case -3
        % outer boundary as polyshape object
        xv1=[0 1 1 0]; yv1=[0 0 1 1]; P1=polyshape(xv1,yv1);
        
        % inner boundary as polyshape object
        xv2=[0.25 0.75 0.5]; yv2=[0.25 0.25 0.5]; P2=polyshape(xv2,yv2);
        
        % domain as difference of polyshape objects
        P=subtract(P1,P2);
        
    case -4
        
        % polygonal domain with two holes.
        
        N=11; a=0; b=2*pi;
        h=(b-a)/N; theta=(a:h:b-h);
        x=cos(theta); y=sin(theta);
        
        % outer boundary as polyshape object
        xv1=1*x+0; yv1=1*y+0; P1=polyshape(xv1,yv1);
        
        % inner boundary (hole) as polyshape object
        xv2=0.15*x+0.5; yv2=0.15*y+0.5; P2=polyshape(xv2,yv2);
        
        % inner boundary (hole) as polyshape object
        xv3=0.15*x-0.5; yv3=0.15*y-0.5; P3=polyshape(xv3,yv3);
        
        % domain as difference of polyshape objects
        P=subtract(P1,P2); P=subtract(P,P3);
        
        
    case -5
        
        % polygonal domain with two holes and a polygon inside a hole
        
        N=11; a=0; b=2*pi;
        h=(b-a)/N; theta=(a:h:b-h);
        x=cos(theta); y=sin(theta);
        
        % outer boundary as polyshape object
        xv1=1*x+0; yv1=1*y+0; P1=polyshape(xv1,yv1);
        
        % inner boundary (hole) as polyshape object
        xv2=0.15*x+0.5; yv2=0.15*y+0.5; P2=polyshape(xv2,yv2);
        
        % inner boundary (hole) as polyshape object
        xv3=0.15*x-0.5; yv3=0.15*y-0.5; P3=polyshape(xv3,yv3);
        
        % outer boundary of an inner polyshape object
        xv4=0.10*x-0.5; yv4=0.10*y-0.5; P4=polyshape(xv4,yv4);
        
        % domain as difference/union of polyshape objects
        P=subtract(P1,P2); P=subtract(P,P3); P=union(P,P4);
        
    case -6
        
        % polygonal domain with two holes, a polygon inside a hole, and a
        % disconnected outer square
        
        N=10; a=0; b=2*pi;
        h=(b-a)/N; theta=(a:h:b-h);
        x=cos(theta); y=sin(theta);
        
        % outer boundary as polyshape object
        xv1=1*x+0; yv1=1*y+0; P1=polyshape(xv1,yv1);
        
        % inner boundary (hole) as polyshape object
        xv2=0.15*x+0.5; yv2=0.15*y+0.5; P2=polyshape(xv2,yv2);
        
        % inner boundary (hole) as polyshape object
        xv3=0.15*x-0.5; yv3=0.15*y-0.5; P3=polyshape(xv3,yv3);
        
        % outer boundary of an inner polyshape object
        xv4=0.10*x-0.5; yv4=0.10*y-0.5; P4=polyshape(xv4,yv4);
        
        % outer boundary of an outer polyshape object
        xv5=0.50*x+2; yv5=0.50*y+2; P5=polyshape(xv5,yv5);
        
        % domain as difference/union of polyshape objects
        P=subtract(P1,P2); P=subtract(P,P3); P=union(P,P4); P=union(P,P5);
        
    case 0
        % fprintf('\n \t [POLYGON]: SHIFTED UNIT SQUARE [0,1]^2');
        
        % polygon defined via its vertices
        polygon_sides=[0 0; 1 0; 1 1; 0 1]-1/2;
        
        % polygon in polyshape form
        P=polyshape(polygon_sides(:,1),polygon_sides(:,2));
        
        
        
        
    case 1
        % fprintf('\n \t [POLYGON]: UNIT SQUARE [0,1]^2');
        
        % polygon defined via its vertices
        k=1;
        polygon_sides=k*[0 0; 1 0; 1 1; 0 1];
        
        % polygon in polyshape form
        P=polyshape(polygon_sides(:,1),polygon_sides(:,2));
        
    case 2
        % fprintf('\n \t [POLYGON]: CONVEX POLYGON');
        
        % polygon defined via its vertices
        polygon_sides=[0.1 0; 0.7 0.2; 1 0.5; 0.75 0.85; 0.5 1; 0 0.25];
        
        % polygon in polyshape form
        P=polyshape(polygon_sides(:,1),polygon_sides(:,2));
        
    case 3
        % fprintf('\n \t [POLYGON]: NON CONVEX POLYGON H1');
        
        % polygon defined via its vertices
        polygon_sides=(1/4)*[1 0; 3 2; 3 0; 4 2; 3 3; 3 0.85*4; 2 4; ...
            0 3; 1 2];
        
        % polygon in polyshape form
        P=polyshape(polygon_sides(:,1),polygon_sides(:,2));
        
    case 4
        % fprintf('\n \t [POLYGON]: NON CONVEX POLYGON H2');
        
        % polygon defined via its vertices
        polygon_sides=(1/4)*[1 0; 3 2; 3 0; 4 2; 3 3; 3 4; 2 4; 0 3; 1 2];
        
        % polygon in polyshape form
        P=polyshape(polygon_sides(:,1),polygon_sides(:,2));
        
    case 5
        % fprintf('\n \t [POLYGON]: NON SIMPLE POLYGON');
        
        % polygon defined via its vertices
        k=1/4;
        polygon_sides=k*[0 0; 1 2; 2 0; 3 2; 4 0; 4 4; 3 2; 2 4; ...
            1 2; 0 4];
        
        % polygon in polyshape form
        P=polyshape(polygon_sides(:,1),polygon_sides(:,2));
        
    case 6
        % fprintf('\n \t [POLYGON]: UNIT SQUARE [-1,1]^2');
        
        % polygon defined via its vertices
        polygon_sides=[-1 -1; 1 -1; 1 1; -1 1];
        
        % polygon in polyshape form
        P=polyshape(polygon_sides(:,1),polygon_sides(:,2));
        
    case 7
        % fprintf('\n \t [POLYGON]: UNIT TRIANGLE');
        
        % polygon defined via its vertices
        polygon_sides=[0 0; 1 0; 1 1];
        
        % polygon in polyshape form
        P=polyshape(polygon_sides(:,1),polygon_sides(:,2));
        
    case 8
        % fprintf('\n \t [POLYGON]: PSEUDO-CIRCLE');
        
        % polygon defined via its vertices
        N=20; a=0; b=2*pi;
        h=(b-a)/N; theta=(a:h:b-h)';
        x=cos(theta); y=sin(theta);
        x=0.5*x+0.5; y=0.5*y+0.5;
        polygon_sides=[x y];
        
        % polygon in polyshape form
        P=polyshape(polygon_sides(:,1),polygon_sides(:,2));
        
    otherwise
        
        % fprintf('\n \t [POLYGON]: PSEUDO-CIRCLE');
        
        % polygon defined via its vertices
        N=10; a=0; b=2*pi;
        h=(b-a)/N; theta=(a:h:b-h)';
        x=cos(theta); y=sin(theta);
        x=0.5*x+0.5; y=0.5*y+0.5;
        polygon_sides=[x y];
        
        % polygon in polyshape form
        P=polyshape(polygon_sides(:,1),polygon_sides(:,2));
        
        
        
end



function [f,fstr]=integrands(function_type)

%--------------------------------------------------------------------------
% OBJECT
%--------------------------------------------------------------------------
% Definition of the integrand.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% function_type: parameter that determines the integrand.
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% f: bivariate function
% fstr: string of the bivariate function for final statistics.
%--------------------------------------------------------------------------

switch function_type
    case 1
        fstr='f=@(x,y) franke(x,y);';
        f=@(x,y) franke(x,y);
    case 2
        fstr='f=@(x,y) ( (x-0.5).^2 +(y-0.5).^2 ).^(1/2);';
        f=@(x,y) ( (x-0.5).^2 +(y-0.5).^2 ).^(1/2);
    case 3
        fstr='f=@(x,y) (x+y).^19;';
        f=@(x,y) (x+y).^19;
    case 4
        fstr='f=@(x,y) exp(x-0.5).^2+(y-0.5).^2;';
        f=@(x,y) exp(x-0.5).^2+(y-0.5).^2;
    case 5
        fstr='f=@(x,y) exp(-100*((x-0.5).^2+(y-0.5).^2));';
        f=@(x,y) exp(-100*((x-0.5).^2+(y-0.5).^2));
    case 6
        fstr='f=@(x,y) cos(10*(x+y));';
        f=@(x,y) cos(10*(x+y));
    case 7
        fstr='f=@(x,y) ones(size(x))+zeros(size(y));';
        f=@(x,y) ones(size(x))+zeros(size(y));
    case 8
        fstr='f=@(x,y) exp(5*(x-y));';
        f=@(x,y) exp(5*(x-y));
    case 9
        fstr='f=@(x,y) exp(-x.^2-y.^2);';
        f=@(x,y) exp(-x.^2-y.^2);
    case 10
        fstr='f=@(x,y) 1./(1+25*(x.^2+y.^2));';
        f=@(x,y) 1./(1+25*(x.^2+y.^2));
    case 11
        fstr='f=@(x,y) 1./(1+25*(abs(x)+abs(y)));';
        f=@(x,y) 1./(1+25*(abs(x)+abs(y)));
    case 12
        fstr='f=@(x,y) cos(20*(x+y));';
        f=@(x,y) cos(20*(x+y));
    otherwise
        fstr=strcat('f=@(x,y) (x+y).^function_type')'
        f=@(x,y) (x+y).^function_type;
end



function centers=generate_centers_polyshape(P,N)

%--------------------------------------------------------------------------
% OBJECT
%--------------------------------------------------------------------------
% This routine defines "N" RBF centers, inside the polygonal domain
% represented by the polyshape object P.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% P       : polygon definition (polyshape form)
% N       ; number of centers.
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% centers: N x 2 matrix, describing the coordinates of the centers.
%--------------------------------------------------------------------------

% .......................... random points ................................
[xv,yv]=boundary(P);
xmin=min(xv); xmax=max(xv); ymin=min(yv); ymax=max(yv);

PP=rand_pointset_15000;
X0=xmin+(xmax-xmin)*PP(:,1);
Y0=ymin+(ymax-ymin)*PP(:,2);

% ........................ inpolygon vertices .............................

in=isinterior(P,X0,Y0);
X0=X0(find(in == 1));
Y0=Y0(find(in == 1));

centers=[X0(1:N) Y0(1:N)];


