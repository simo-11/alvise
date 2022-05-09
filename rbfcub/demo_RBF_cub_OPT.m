
function demo_RBF_cub_OPT

%--------------------------------------------------------------------------
% * Object:
%--------------------------------------------------------------------------
% In this demo we test routines on RBF cubature over polygons via
% triangulation ("polyshape" version).
%
% With respect to the basic version, an "almost" optimal shape parameter is
% determined.
%
% Note that the battery of tests is used in the reference paper mentioned
% below.
%--------------------------------------------------------------------------
% Reference papers:
% 1. R. Cavoretto, A. De Rossi, A. Sommariva, M. Vianello
% RBFCUB: a numerical package for near-optimal meshless cubature on general
% polygons
%
% 2. A. Sommariva and M. Vianello,
% RBF moment computation and meshless cubature on general polygonal regions.
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
%--------------------------------------------------------------------------
% Note: In the tests of the reference paper: "RBF_typeV=[4 5 6]".
%    In other words, we make experiments by means of Wendland 2, TPS and
%    r^3 polyharmonic spline.
%--------------------------------------------------------------------------

RBF_typeV=[1:4 10 11 13 14]; %[2 3 11 4]; %[4 5 6];

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
% Note:
%---------
% 1. In the tests of the reference paper: "polygon_typeV=[-6 4]".
% 2. To introduce new polygons, see "define_polygon_polyshape" below.
%--------------------------------------------------------------------------
% Example:
%---------
% By choosing "polygon_typeV=[-6 4];" we make all experiments using
%
% A. polygon with 2 polygonal holes (non simply connected polygonal reg.)
%    and a polygon inside a hole, plus a disconnected polygon.
%
% B. not convex polygon (9 sides),
%--------------------------------------------------------------------------

polygon_typeV=[-6 4];

%--------------------------------------------------------------------------
% * function_typeV
%--------------------------------------------------------------------------
% 1: f=@(x,y) franke(x,y);
% 2: f=@(x,y) ( (x-0.5).^2 +(y-0.5).^2 ).^(1/2);
% 3: f=@(x,y) (x+y).^19;
% 4: f=@(x,y) exp(x-0.5).^2+(y-0.5).^2;
% 5: f=@(x,y) exp(-100*((x-0.5).^2+(y-0.5).^2));
% 6: f=@(x,y) cos(10*(x+y));
% 7: f=@(x,y) ones(size(x))+zeros(size(y));
% 8: f=@(x,y) exp(5*(x-y));
% 9: f=@(x,y) exp(-x.^2-y.^2);
% 10: f=@(x,y) 1./(1+25*(x.^2+y.^2));
% 11: f=@(x,y) 1./(1+25*(abs(x)+abs(y)));
% 12: f=@(x,y) cos(20*(x+y));
% 13: f=@(x,y) exp(x-y);
% 14: f=@(x,y) ( (x-0.3).^2 +(y-0.3).^2 ).^(1/2);
%--------------------------------------------------------------------------
% Note:
% 1. In the tests of the reference paper: "function_typeV=[13 8 14]".
% 2. If you need to add some new integrands see subroutine "integrands"
%    below.
%--------------------------------------------------------------------------

function_typeV=1; %14; % choose one function alone

%--------------------------------------------------------------------------
% Note: scattered data are loaded from the file "rand_pointset_15000.m".
%--------------------------------------------------------------------------
NcentersV=[200 400 800];            % number of centers

fixed_scale=2;          % RBF scale: 0: random
                        %            1: fixed
                        %            2: optimal by a variant of LOOCV

doplot=0;               % plot domain & centers: 0: no 1: yes

ade=100; % algebraic degree of precision of the reference rule




% ...................... main code below ..................................

% .................. perform battery of experiments .......................

aeV=[];
reV=[];
cpus2V=[];
VcondV=[];
RBFV=[];
function_typeVV=[];
NcentersVV=[];
polygonV=[];



% ... define range for searching a (near) optimal RBF shape parameter,
%     e.g. epsilon in [minep,maxep]=[0.5,15]
minep=.5; maxep=15;

for jj=1:length(polygon_typeV)
    
    % .... define polygonal region as polyshape object
    polygon_type=polygon_typeV(jj);
    P=define_polygon_polyshape(polygon_type);
    
    % .... reference cubature rule of algebraic type ....
    fprintf('\n \t ... computing reference high order cubature rule ...');
    [xyw,xvc,yvc,pgon,tri]=polygauss_2018(ade,P);
    
    % .... RBF cubature ....
    fprintf('\n \t ... tests on RBF cubature rule ...');
    for ii=1:length(RBF_typeV)
        RBF_type=RBF_typeV(ii);
        for ss=1:length(NcentersV)
            
            % .... Define RBF centers ....
            Ncenters=NcentersV(ss);
            centers=generate_centers_polyshape(P,Ncenters);
            
            % .... define function ....
            mm=1;
            function_type=function_typeV(mm);
            [f,fstr]=integrands(function_type);
            fcub=feval(f,centers(:,1),centers(:,2));
            
            % .... RBF scale ....
            % In this part we determine the scale of the RBF. If
            % 1. fixed_scale == 0 it is random,
            % 2. fixed_scale == 1 it is 1,
            % otherwise it is determined by a variant of LOOCV.
            
            if fixed_scale == 0
                RBF_scale=0.35+rand(Ncenters,1);
            else
                if fixed_scale == 1
                    RBF_scale=ones(Ncenters,1);
                else
                    if (RBF_type >= 5 & RBF_type <= 8)
                        RBF_scale=ones(Ncenters,1);
                    else
                        RBF_scale=RBF_optimal_scale(centers,fcub,...
                            RBF_type,minep,maxep);
                    end
                end
            end
            
            % .... RBF cubature weights ....
            [weights,cpus,Vcond]=RBF_cub_polygon_OPT(P,centers,RBF_type,...
                RBF_scale);
            
            % .... reference integral IR ....
            fxy=feval(f,xyw(:,1),xyw(:,2)); IR=(xyw(:,3))'*fxy;
            
            % .... integration by scattered data ....
            cubature_value=weights'*fcub;
            
            % .... statistics ....
            [ae,re]=statistics(RBF_type,RBF_scale,function_type,...
                polygon_type,cubature_value,IR,weights,cpus,...
                centers,P,fixed_scale,doplot,fstr,Ncenters,Vcond);
            
            % .... statistics storage for final summary ....
            aeV=[aeV; ae]; reV=[reV; re]; cpus2V=[cpus2V; cpus];
            VcondV=[VcondV; Vcond]; RBFV=[RBFV; RBF_type];
            function_typeVV=[function_typeVV; function_type];
            NcentersVV=[NcentersVV; Ncenters];
            polygonV=[polygonV; polygon_type];
        end
    end
end


% ........................ summary of experiments .........................

fprintf('\n \t ............. summary cpus ...................................');
fprintf('\n \t test | int.mat. | cmp.mom. | add.ter. | lin.sys. |  tot    |');
fprintf('\n \t .............................................................. ');
for k=1:size(cpus2V,1)
    fprintf('\n \t %3.0f  |  %1.1e |  %1.1e |  %1.1e |  %1.1e | %1.1e |',k,...
        cpus2V(k,1),cpus2V(k,2),cpus2V(k,3),cpus2V(k,4),sum(cpus2V(k,1:4)));
end

fprintf('\n \n');

fprintf('\n \t ............. summary ...................................................');
fprintf('\n \t test |abs.err.|rel.err.| cond i.m.| RBF | ftype | Ncents | Polyg. |');
fprintf('\n \t ......................................................................... ');
for k=1:length(aeV)
    aeL=aeV(k); reL=reV(k);
    VcondVL=VcondV(k);
    function_typeVVL=function_typeVV(k);
    NcentersL=NcentersVV(k);
    polygonVL=polygonV(k);
    fprintf('\n \t %3.0f  |  %1.0e |  %1.0e |  %1.1e | %2.0f  |  %2.0f   |  %2.0f   |  %2.0f   |',...
        k,aeL,reL, VcondVL,RBFV(k),function_typeVVL,NcentersL,polygonVL);
end
fprintf('\n \t ......................................................................... \n \n ');








function [ae,re,cpus,Vcond]=statistics(RBF_type,RBF_scale,function_type,...
    polygon_type,cubature_value,IR,weights,cpus,centers,P,...
    fixed_scale,doplot,fstr,Ncenters,Vcond)


%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% Statistics of our numerical experiments.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% MAIN PROGRAM STARTS HERE.
%--------------------------------------------------------------------------


% .......................... statistics ...................................

fprintf('\n \n \t *** STATISTICS ***');

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
    case 10, RBF_string=' WENDLAND, W4:(35*r^2+18*r+3)*max(0,(1-r))^6';
    case 11, RBF_string=' WENDLAND, W6:(32*r^3+25*r^2+8*r+1)*max(0,(1-r))^8';
    case 12, RBF_string=' MISSING WENDLAND:';
    case 13, RBF_string=' MATERN B1=(d+1)/2: exp(-r)';
    case 14, RBF_string=' MATERN B2=(d+3)/2: (1+r)*exp(-r)';
end
fprintf('\n \t RBF TYPE     : %2.0f',RBF_type);
fprintf('\n \t *** '); disp(RBF_string);

% ... RBF scale ...

switch fixed_scale
    
    case 0
        fprintf('\n \t RBF SCALE    : RANDOM, min: %1.1e max: %1.1e',...
            min(RBF_scale),max(RBF_scale));
    case 1
        fprintf('\n \t RBF SCALE    : % 1.5e',RBF_scale(1));
    otherwise
        fprintf('\n \t RBF SCALE    : OPTIMAL, %1.1e'RBF_scale));
end

% ... integrand ...
fprintf('\n \n \t FUNCTION TYPE: %2.0f',function_type);
fprintf('\n \t *** '); disp(fstr);

% ... experiment data ...
fprintf('\n \t POLYGON TYPE: %4.0f',polygon_type);
fprintf('\n \t NODES CARDIN.: %4.0f',Ncenters);

% ... cubature values ...

fprintf('\n \n \t * CUBATURE RES.: %1.15e',cubature_value);
fprintf('\n \t * EXACT VALUE  : %1.15e',IR);

% ... cubature errors ...

ae=abs(IR-cubature_value);
fprintf('\n \t ** ABSOLUTE ERR.: %1.3e',ae);
if abs(IR) > 0
    re=ae/abs(IR);
    fprintf('\n \t ** RELATIVE ERR.: %1.3e',re);
else
    re=NaN;
end

% ... cubature condition ...

fprintf('\n \t * norm(w,1) : %1.3e',norm(weights,1));
fprintf('\n \t * COND. VAND.MATRIX  : %1.3e',Vcond);
if min(weights) > 0
    fprintf('\n \t NOTE: POSITIVE WEIGHTS FORMULA');
else
    fprintf('\n \t NOTE: FORMULA WITH NEGATIVE WEIGHTS');
end

% ... cputimes ...

fprintf('\n \n \t CPU: CUBATURE MATRIX %1.1e',cpus(1));
fprintf('\n \t CPU: MOMENTS         %1.1e',cpus(2));
fprintf('\n \t CPU: POLYNO. MOMENTS %1.1e',cpus(3));
fprintf('\n \t CPU: WEIGHTS COMPUT. %1.1e',cpus(4));
fprintf('\n \t CPU: TOTAL TIME      %1.1e',sum(cpus));
fprintf('\n \n');

fprintf('\n \n');

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
        
        % square doughnut
        
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
        
    case -7
        
        % polygonal cat: we show here the interesting features of polyshape
        
        N=20;
        
        t=linspace(0,2*pi,N); t(end)=0;
        
        % building face
        P0=[cos(t') sin(t')]; P0F=polyshape(P0(:,1),P0(:,2));
        
        % building ears
        P1x=[0.95*cos(pi/3) 0.95*cos(pi/6) 1]';
        P1y=[0.95*sin(pi/3) 0.95*sin(pi/6) 1]';
        P1F=polyshape(P1x,P1y);
        
        % building right eye
        P2x=[0.95*cos(pi/3+pi/2) 0.95*cos(pi/6+pi/2) -1]';
        P2y=[0.95*sin(pi/3+pi/2) 0.95*sin(pi/6+pi/2) 1]';
        P2F=polyshape(P2x,P2y);
        
        % building left eye
        P3x=0.3+0.2*cos(t'); P3y=0.3+0.2*sin(t'); P3F=polyshape(P3x,P3y);
        P4x=-0.3+0.2*cos(t'); P4y=0.3+0.2*sin(t'); P4F=polyshape(P4x,P4y);
        
        % building pupils
        P5x=0.3+0.1*cos(t'); P5y=0.3+0.1*sin(t'); P5F=polyshape(P5x,P5y);
        P6x=-0.3+0.1*cos(t'); P6y=0.3+0.1*sin(t'); P6F=polyshape(P6x,P6y);
        
        % building mouth
        t0=linspace(3*pi/2-pi/3,3*pi/2+pi/3,N); t=t0'; tf=flipud(t);
        P7x=0.45*cos(t); P7x=[P7x; 0.4*cos(tf)];
        P7y=0.45*sin(t); P7y=[P7y; 0.4*sin(tf)];
        P7F=polyshape(P7x,P7y);
        
        % building cat as polyshape object
        P=union(P0F,P1F); P=union(P,P2F);
        P=subtract(P,P3F); P=subtract(P,P4F);
        P=union(P, P5F); P=union(P, P6F);
        P=subtract(P, P7F);
        
        % format long e;
        % [XX,YY]=boundary(P);
        % plot(P)
        % INPOLY = isinterior(P, 0.3, 0.3)
        
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
    case 13
        fstr='f=@(x,y) exp(x-y);';
        f=@(x,y) exp(x-y);
    case 14
        fstr='f=@(x,y) ( (x-0.3).^2 +(y-0.3).^2 ).^(1/2);';
        f=@(x,y) ( (x-0.3).^2 +(y-0.3).^2 ).^(1/2);
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


