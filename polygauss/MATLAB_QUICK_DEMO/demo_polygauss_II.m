
function demo_polygauss_II(example)

%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% IN THIS DEMO WE COMPARE, IF POSSIBLE, THE RESULTS OBTAINED BY
% polygauss_2018).
% THESE RULES HAVE ALG. DEGREE OF EXACTNESS EQUAL TO ade.
% THE CODE DISPLAYS CARDINALITIES, RESULTS ON INTEGRATING A NON SYMMETRIC
% POLYNOMIAL OF DEGREE ade.
% SOME PLOTS OF POINTSETS, DOMAINS, TRIANGULATIONS, ARE GIVEN.
%--------------------------------------------------------------------------
% IMPORTANT:
% * THE CODE REQUIRES THE TOOLBOX CONTAINING POLYSHAPE (AVAILABLE FROM
% MATLAB R2017b VERSIONS).
% * COMPRESSION COMES FROM comprexcub.m THAT MUST BE APPLIED SEPARATELY
% FROM polygauss_2018.
% * IN THIS DEMO WE SHOW HOW TO TREAT DOMAINS THAT MAY BE NOT SIMPLY
% CONNECTED AND/OR EVEN NOT CONNECTED, VIA polyshape OPERATIONS.
%--------------------------------------------------------------------------
% INPUT:
% example: 1: not simply connected and not connected domain. 
%          2: not simply connected domain.
%--------------------------------------------------------------------------
% ROUTINES:
% 1. polygauss_2017 (external)
% 2. polygauss_2018 (external)
% 3. comprexcub (external)
%--------------------------------------------------------------------------
% EXAMPLE:
% >> demo_polygauss_II(1)
% 
%  *** DATA ***
%   EXAMPLE:       1
%   ADE    :      10
%  *** CARDINALITIES OF THE FORMULA ***
%   #FULL:    7080
%   #COMP:      66
%  *** CUBATURE RESULTS ***
%   IFULL: 5.558449132098099e+09
%   ICOMP: 5.558449132098125e+09
%  *** LEGEND ***
%   IFULL: NEW FULL FORMULA
%   ICOMP: NEW COMPRESSED FORMULA 
%  
% >> demo_polygauss_II(2)
% 
%  *** DATA ***
%   EXAMPLE:       2
%   ADE    :      10
%  *** CARDINALITIES OF THE FORMULA ***
%   #FULL:    4848
%   #COMP:      66
%  *** CUBATURE RESULTS ***
%   IFULL: 9.318964220022317e+03
%   ICOMP: 9.318964220022319e+03
%  *** LEGEND ***
%   IFULL: NEW FULL FORMULA
%   ICOMP: NEW COMPRESSED FORMULA 
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
    case 1 % domain not connected and not simply connected.
        Nsides=100;
        th=linspace(0,2*pi,Nsides); th=(th(1:end-1))';
        polygon1=[cos(th) sin(th)]; P1=polyshape(polygon1); % first polygon.
        polygon2=2+[cos(th) sin(th)]; P2=polyshape(polygon2); % second polygon.
        polygon3=0.5*[cos(th) sin(th)]; P3=polyshape(polygon3); % first polygon.
        
        P=subtract(P1,P3);
        P=union(P,P2);
        
    case 2 % domain not simply connected (optics)
        Nsides=100;
        y=[0         0   -0.1184   -0.1184   -0.3761];
        r=[1.0000    0.6120    0.5663    1.0761    1.2810];
        th=linspace(0,2*pi,Nsides); th=(th(1:end-1))';
        C1=[0 y(1)]; P1v=C1+r(1)*[cos(th) sin(th)]; P1=polyshape(P1v);
        C2=[0 y(2)]; P2v=C2+r(2)*[cos(th) sin(th)]; P2=polyshape(P2v);
        C3=[0 y(3)]; P3v=C3+r(3)*[cos(th) sin(th)]; P3=polyshape(P3v);
        C4=[0 y(4)]; P4v=C4+r(4)*[cos(th) sin(th)]; P4=polyshape(P4v);
        C5=[0 y(5)]; P5v=C5+r(5)*[cos(th) sin(th)]; P5=polyshape(P5v);
        Pout=intersect(P1,P4);
        Pout=intersect(Pout,P5);
        Pin=union(P2,P3);
        P=subtract(Pout,Pin);
        
end



%--------------------------
% defining function.
%--------------------------
f=@(x,y) exp(1)+(0.5*x+pi*y).^ade;


%--------------------------
% determining full rule
%--------------------------
[xyw,xvc,yvc,pgon,tri]=polygauss_2018(ade,P);

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
% results
%------------------------------
fprintf('\n *** DATA ***');
fprintf('\n  EXAMPLE: %7.0f',example);
fprintf('\n  ADE    : %7.0f',ade);


fprintf('\n *** CARDINALITIES OF THE FORMULA ***');
fprintf('\n  #FULL: %7.0f',size(xyw,1));
fprintf('\n  #COMP: %7.0f',size(pts,1));

fprintf('\n *** CUBATURE RESULTS ***');
fprintf('\n  IFULL: %1.15e',IF);
fprintf('\n  ICOMP: %1.15e',IC);

fprintf('\n *** LEGEND ***');
fprintf('\n  IFULL: NEW FULL FORMULA');
fprintf('\n  ICOMP: NEW COMPRESSED FORMULA \n \n');


%--------------------------------------------------------------------------
% plots
%--------------------------------------------------------------------------

clf;
close;

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












