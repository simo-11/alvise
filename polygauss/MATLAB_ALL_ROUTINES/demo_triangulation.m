
function demo_triangulation(domain)

%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% IN THIS DEMO WE COMPARE, IF POSSIBLE, THE TRIANGULATIONS OF SEVERAL
% POLYGONS, USING polyshape CLASS.
%--------------------------------------------------------------------------
% INPUT.
%--------------------------------------------------------------------------
% domain: 1: disk, 2: cardioid, 3: Bernoulli Lemniscate
%--------------------------------------------------------------------------
% REQUIRED ROUTINES.
%--------------------------------------------------------------------------
% * triangulation_part (ATTACHED HERE).
% * define_domain (ATTACHED HERE).
%--------------------------------------------------------------------------
% IMPORTANT.
%--------------------------------------------------------------------------
% * THE CODE REQUIRES THE TOOLBOX CONTAINING POLYSHAPE (AVAILABLE FROM 
% MATLAB R2017b VERSIONS).
%--------------------------------------------------------------------------
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
%% Date: December 04, 2018
%--------------------------------------------------------------------------

warning on;

if nargin < 1
    domain=1;
end

spec_settingsV=[100 500 1000 2000 3000 4000 5000 10000 20000];

ntest=3; % number of tests


%--------------------------------------------------------------------------
%                        MAIN CODE STARTS HERE
%--------------------------------------------------------------------------

hist=[];

for ii=1:length(spec_settingsV)
    spec_settings=spec_settingsV(ii);
    fprintf('\n \t =====================================================');
    fprintf('\n \t **** spec_settings: %5.0f',spec_settings);
    
    [xv,yv,iv]=define_domain(domain,spec_settings);
    
    for jj=1:ntest
        [xvc,yvc,pgon,tri,N,tL]=triangulation_part(xv,yv,iv);
        tt(jj,:)=tL;
    end
    
    
    if ntest > 0
        fprintf('\n \t SIDES: %5.0f TRIANGLES: %5.0f TP: %2.2e TT: %2.2e %2.2e',...
            length(xv),N,mean(tt(:,1)),mean(tt(:,1)));
        hist=[hist; length(xv) N mean(tt(:,1)) mean(tt(:,2))];
    end
    
end

fprintf('\n \t =====================================================');
fprintf('\n \t                    LEGEND');
fprintf('\n \t *** TP: AVERAGE POLYSHAPE CPUTIME');
fprintf('\n \t *** TT: AVERAGE TRIANGULATION CPUTIME');
fprintf('\n \t =====================================================');


clf;
hold on;
axis square
triplot(tri);
hold off;


if ntest > 0
    fprintf('\n \t    L  |    N   |     TP    |    TT   |');
    for ii=1:size(hist,1)
        fprintf('\n \t %5.0f | % 5.0f  |  %1.1e  | %1.1e |',hist(ii,1),...
            hist(ii,2),hist(ii,3),hist(ii,4));
    end
end

fprintf('\n \t =====================================================');
fprintf('\n \t                    LEGEND');
fprintf('\n \t *** L : NUMBER OF SIDES');
fprintf('\n \t *** N : NUMBER OF TRIANGLES');
fprintf('\n \t *** TP: AVERAGE POLYSHAPE CPUTIME');
fprintf('\n \t *** TT: AVERAGE TRIANGULATION CPUTIME');
fprintf('\n \t =====================================================');


fprintf('\n \n');




function [xvc,yvc,pgon,tri,N,cput]=triangulation_part(xv,yv,iv)

%--------------------------------------------------------------------------
% Important: this polygauss version requires at least Matlab 9.3.0.713579.
%--------------------------------------------------------------------------
% Input:
%--------------------------------------------------------------------------
% ade: algebraic degree of precision of the rule
% xv: abscissae of the vertices.
% yv: ordinates of the vertices.
% iv: index of the polygon components. Example: if part of the boundary is
%     made by the first 5 coordinates and two other non connected
%     components by other 10 and 12 coordinates, then iv=[5,10,12].
% Important: differently from classical codes on polygons, the first and
% last vertex must NOT be equal.
%--------------------------------------------------------------------------
% Output:
%--------------------------------------------------------------------------
% xyw: N x 3 matrix, where (x,y) with x=xyw(:,1), y=xyw(:,2) are the nodes
%      and w=xyw(:,3) are the weights.
% xvc, yvc: cells of the vertices.
% pgon: polygon coded in polyshape form.
% tri: polygon triangulation in polyshape form.
% N: number of triangles.
% cput: 1 x 2 vector containing polyshape and triangulation time.
%--------------------------------------------------------------------------



% TROUBLESHOOTING.
if nargin < 2
    pts=xv;
    if size(pts,1) <= size(pts,2) % working with column vectors.
        pts=xv';
    end
    if pts(1,:) == pts(end,:)
        pts=pts(1:end-1,:);
    end
    xv=pts(:,1); yv=pts(:,2); iv=length(xv);
else
    if size(xv,1) <= size(xv,2) % working with column vectors.
        xv=xv';
    end
    if size(yv,1) <= size(yv,2) % working with column vectors.
        yv=yv';
    end
    %     if [xv(1) yv(1)] == [xv(end) yv(end)]
    %         xv=xv(1:end-1); yv=yv(1:end-1);
    %     end
    if nargin < 3
        iv=length(xv);
    end
end

% Here we convert the given information into cells of coordinates, where
% each component represent a geometrical component of the polygon, not
% connected or at most tangent to other geometrical components.

xvc={xv(1:iv(1))}; yvc={yv(1:iv(1))};
ii1=1; ii2=iv(1);
for ii=2:length(iv)
    ii1=ii2+1;
    ii2=ii2+iv(ii);
    xvl=xv(ii1:ii2);
    yvl=yv(ii1:ii2);
    
    xvc{end+1}=xvl; yvc{end+1}=yvl;
end

tic;
pgon = polyshape(xvc,yvc);
tp=toc;

tic;
tri = triangulation(pgon);
tt=toc;

% fprintf('\n \t \t t: %1.3e tp: %1.3e',tt,tp);

TCL=tri.ConnectivityList;
X = tri.Points(:,1);
Y = tri.Points(:,2);

N=size(TCL,1);

cput=[tp tt];








function [xv,yv,iv]=define_domain(domain,spec_settings)

%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% IN THIS DEMO WE COMPARE, IF POSSIBLE, THE TRIANGULATIONS OF SEVERAL
% POLYGONS, USING polyshape CLASS.
%--------------------------------------------------------------------------
% INPUT.
%--------------------------------------------------------------------------
% domain: 1: disk, 2: cardioid, 3: Bernoulli Lemniscate
% spec_settings: number of vertices 
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
%% Date: November 30, 2018
%--------------------------------------------------------------------------


xv=[];

switch domain
    case 1 % DISK. 
        if isempty(spec_settings) == 1
            spec_settings=32;
        end
        th=linspace(0,2*pi,spec_settings);
        %th=(th(1:end-1))';
        polygon_sides=[cos(th') sin(th')];
        polygon_sides=polygon_sides(1:end-1,:);
        
    case 2 % CARDIOID.
        if isempty(spec_settings) == 1
            spec_settings=32;
        end
        th=linspace(0,2*pi,spec_settings+1);
        %th=(th(1:end-1))';
        polygon_sides=[cos(th').*(1-cos(th')) sin(th').*(1-cos(th'))];
        polygon_sides=polygon_sides(1:end-1,:);
        
    case 3 % BERNOULLI LEMNISCATE.
        a=1;
        if isempty(spec_settings) == 1
            spec_settings=33;
        end
        th=linspace(0,2*pi,spec_settings); th=(th(1:end-1))';
        xv=a*sqrt(2)*cos(th)./(1+(sin(th)).^2);
        yv=(a*sqrt(2)*cos(th).*sin(th))./(1+(sin(th)).^2);
        iv=length(xv);
end


if isempty(xv)
    xv=polygon_sides(:,1); yv=polygon_sides(:,2); iv=length(xv);
end
