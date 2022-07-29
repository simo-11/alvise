
function demo_polygauss(example)

%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% IN THIS DEMO WE COMPARE, IF POSSIBLE, THE RESULTS OBTAINED BY THE OLD
% POLYGAUSS (i.e. polygauss_2017) AND THE NEW ONE (i.e. polygauss_2018).
% THE OLD CODES COULD INTEGRATE ON QUADRILATERALS OR VIA TRIANGULATION ON
% SIMPLE POLYGONS.
% THESE RULES HAVE ALG. DEGREE OF EXACTNESS EQUAL TO ade.
% THE CODE DISPLAYS CARDINALITIES, RESULTS ON INTEGRATING A NON SYMMETRIC 
% POLYNOMIAL OF DEGREE ade.
% SOME PLOTS OF POINTSETS, DOMAINS, TRIANGULATIONS, ARE GIVEN.
%--------------------------------------------------------------------------
% IMPORTANT:
% * THE CODE REQUIRES THE TOOLBOX CONTAINING POLYSHAPE (AVAILABLE FROM
% MATLAB R2017b VERSIONS).
% * TEST EXAMPLES COMES FROM THE FILE mydomain.m
% * COMPRESSION COMES FROM comprexcub.m THAT MUST BE APPLIED SEPARATELY
% FROM polygauss_2018.
% * COMPARISONS WITH OLDER VERSIONS REQUIRE polygauss_2017.m
%--------------------------------------------------------------------------
% INPUT:
% example: 1: convex domain. 2: non convex domain.
%--------------------------------------------------------------------------
% ROUTINES:
% 1. polygauss_2017 (external)
% 2. polygauss_2018 (external) 
% 3. comprexcub (external)
%--------------------------------------------------------------------------
% EXAMPLE:
% >>demo_polygauss(1)
% 
%  *** DATA ***
%   EXAMPLE:       1
%   ADE    :      10
%  *** CARDINALITIES OF THE FORMULA ***
%   #FULL:      96
%   #COMP:      66
%   #GG  :     660
%  *** CUBATURE RESULTS ***
%   IFULL: 4.383008816374430e+03
%   ICOMP: 4.383008816374429e+03
%   IGG  : 4.383008816374429e+03
%  *** LEGEND ***
%   IFULL: NEW FULL FORMULA
%   ICOMP: NEW COMPRESSED FORMULA
%   IGG  : OLD GAUSS-GREEN FORMULA 
%  
% demo_polygauss(2)
% 
%  *** DATA ***
%   EXAMPLE:       2
%   ADE    :      10
%  *** CARDINALITIES OF THE FORMULA ***
%   #FULL:     168
%   #COMP:      66
%   #GG  :     870
%  *** CUBATURE RESULTS ***
%   IFULL: 5.710403316810376e+03
%   ICOMP: 5.710403316810357e+03
%   IGG  : 5.710403316810378e+03
%  *** LEGEND ***
%   IFULL: NEW FULL FORMULA
%   ICOMP: NEW COMPRESSED FORMULA
%   IGG  : OLD GAUSS-GREEN FORMULA 
%  
% >>
%--------------------------------------------------------------------------
%% Copyright (C) 2018- Brian J. Bauman, Alvise Sommariva, Marco Vianello.
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
%% Author:  Brian J. Bauman  <bauman3@llnl.gov>
%%          Alvise Sommariva <alvise@math.unipd.it>
%%          Marco Vianello   <marcov@math.unipd.it>
%%
%% Date: December 07, 2018
%--------------------------------------------------------------------------

if nargin < 1, example=1; end

ade=10; % algebraic degree of precision.
pos=2;  % 0: qr compression, 1: lsqnonneg compression 2: fast L-H compr.

%--------------------------
% defining domain
%--------------------------
switch example
    case 1
        polygon_sides=[0.1 0; 0.7 0.2; 1 0.5; 0.75 0.85; 0.5 1; 0 0.25; 0.1 0];
    case 2
        polygon_sides=(1/4)*[1 2; 1 0; 3 2; 3 0; 4 2; 3 3; 3 0.85*4; 2 4;
            0 3; 1 2];
end
xv=polygon_sides(:,1); yv=polygon_sides(:,2);
iv=length(xv); % This variable depends on the holes or not connected domain.
               % In these simple cases the domains are without holes and
               % domains are connected. 


%--------------------------
% defining function.
%--------------------------
f=@(x,y) exp(1)+(0.5*x+pi*y).^ade;


%--------------------------
% determining full rule
%--------------------------
[xyw,xvc,yvc,P,tri]=polygauss_2018(ade,xv,yv,iv);

% full rule result.
X=xyw(:,1); Y=xyw(:,2); W=xyw(:,3);
fXY=feval(f,X,Y);
IF=W'*fXY;


%------------------------------
% determining compressed rule.
%------------------------------
[pts,w,momerr] = comprexcub(ade,xyw(:,1:2),xyw(:,3),pos);

% compressed rule result
X=pts(:,1); Y=pts(:,2); W=w;
fXY=feval(f,X,Y);
IC=W'*fXY;


%------------------------------
% determining Gauss-Green rule
%------------------------------
xyw17=polygauss_2017(polygon_sides,ade,'quadrangulation');

% Gauss-Green rule result.
X=xyw17(:,1); Y=xyw17(:,2); W=xyw17(:,3);
fXY=feval(f,X,Y);
IGG=W'*fXY;


%------------------------------
% results
%------------------------------
fprintf('\n *** DATA ***');
fprintf('\n  EXAMPLE: %7.0f',example);
fprintf('\n  ADE    : %7.0f',ade);


fprintf('\n *** CARDINALITIES OF THE FORMULA ***');
fprintf('\n  #FULL: %7.0f',size(xyw,1));
fprintf('\n  #COMP: %7.0f',size(pts,1));
fprintf('\n  #GG  : %7.0f',size(xyw17,1));

fprintf('\n *** CUBATURE RESULTS ***');
fprintf('\n  IFULL: %1.15e',IF);
fprintf('\n  ICOMP: %1.15e',IC);
fprintf('\n  IGG  : %1.15e',IGG);

fprintf('\n *** LEGEND ***');
fprintf('\n  IFULL: NEW FULL FORMULA');
fprintf('\n  ICOMP: NEW COMPRESSED FORMULA');
fprintf('\n  IGG  : OLD GAUSS-GREEN FORMULA \n \n');



%--------------------------------------------------------------------------
% plots
%--------------------------------------------------------------------------






% plotting triangulation

clf(figure(1));
figure(1)
hold on;
axis off
axis equal
titlestr1='Domain triangulation.';
titlestr2=strcat('Triangles: ',num2str(size(tri,1)),'.');
title({[titlestr1];[ titlestr2]})
triplot(tri);
plot(P,'FaceColor',[0.5 0.5 0.5]);
hold off;


% plotting full cub. nodes

clf(figure(2));
figure(2)
hold on;
axis off;
axis equal;
titlestr1='Quadrature points of the new formula (before compression).';
titlestr2=strcat('Cardinality: ',num2str(size(xyw,1)),'. Ade: ',num2str(ade),'.');
title({[titlestr1];[ titlestr2]});
plot(P,'FaceColor',[0.5 0.5 0.5]);
plot(xyw(:,1),xyw(:,2),'ko','LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',2)
% plot(xywc(:,1),xywc(:,2),'ko','LineWidth',1,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','g',...
%     'MarkerSize',6)
hold off;

% plotting compressed cub. nodes

clf(figure(3));
figure(3)
hold on;
axis off
axis equal
titlestr1='Quadrature points of the new formula (after compression).';
titlestr2=strcat('Cardinality: ',num2str(size(pts,1)),'. Ade: ',num2str(ade),'.');
title({[titlestr1];[ titlestr2]})
plot(P,'FaceColor',[0.5 0.5 0.5]);
plot(pts(:,1),pts(:,2),'ko','LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','c',...
    'MarkerSize',6)
hold off;



% plotting Gauss-Green nodes

clf(figure(4));
figure(4)
hold on;
axis off;
axis square
titlestr1='Quadrature points of the old formula, via Gauss-Green';
titlestr2=strcat('Cardinality: ',num2str(size(pts,1)),'. Ade: ',num2str(ade),'.');
title({[titlestr1];[ titlestr2]})
plot([xv; xv(:,1)],[yv; yv(:,1)],'k-','LineWidth',4,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',2);
plot(xyw17(:,1),xyw17(:,2),'ko','LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',4)
hold off;








