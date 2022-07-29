
function demo_polygauss_quatrefoil(example)

%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% IN THIS DEMO WE COMPARE, IF POSSIBLE, THE RESULTS OBTAINED BY 
% polygauss_2018 AND ITS COMPRESSED VERSION.
% THE CODE DISPLAYS MOMENTS ERRORS, CARDINALITIES AND RATIOS.
%--------------------------------------------------------------------------
% IMPORTANT:
% * THE CODE REQUIRES THE TOOLBOX CONTAINING POLYSHAPE (AVAILABLE FROM
% MATLAB R2017b VERSIONS).
% * TEST EXAMPLES COMES FROM THE FILE mydomain.m
% * COMPRESSION COMES FROM comprexcub.m
% * COMPARISONS WITH OLDER VERSIONS REQUIRE polygauss_2017.m
%--------------------------------------------------------------------------
% INPUT:
% example: 1: quatrefoil with M=129, 2: quatrefoil with M=512 (see paper).
%--------------------------------------------------------------------------
% ROUTINES:
% 1. define_domain (attached)
% 2. polygauss_2018 (external)
% 3. comprexcub (external)
%--------------------------------------------------------------------------
% EXAMPLE:
% >> demo_polygauss_quatrefoil(1)
% 
%        =====================================
%              Domain: Omega(QF1)   
%        =====================================
%         delta E(T,TC) M(new) M^(new) ratio
%        -------------------------------------
%  	  5   1e-15     840     21    40.0
%  	 10   2e-15    2880     66    43.6
%  	 15   3e-15    5520    136    40.6
%  	 20   4e-15    9360    231    40.5
%  	 25   7e-15   14040    351    40.0
%  	 30   4e-15   20520    496    41.4
%  	 35   4e-15   27360    666    41.1
%        ------------------------------------------------------
%        *** LEGEND ***
%        delta  : algebraic degree of precision
%        E(T,TC): moms err. (norm 2), triang. vs compr.
%        M(new) : cardinality polygauss via triangulation
%        M^(new): cardinality polygauss via triangulation+compr.
%        ratio  : M(new)/M^(new)
%        ------------------------------------------------------
%  	 ** NUMBER OF TRIANGLES:   120
%
% >> demo_polygauss_quatrefoil(2)
% 
%        =====================================
%              Domain: Omega(QF2)   
%        =====================================
%         delta E(T,TC) M(new) M^(new) ratio
%        -------------------------------------
%  	  5   1e-15    3528     21   168.0
%  	 10   2e-15   12096     66   183.3
%  	 15   3e-15   23184    136   170.5
%  	 20   5e-15   39312    231   170.2
%  	 25   6e-15   58968    351   168.0
%  	 30   1e-14   86184    496   173.8
%  	 35   1e-14  114912    666   172.5
%        ------------------------------------------------------
%        *** LEGEND ***
%        delta  : algebraic degree of precision
%        E(T,TC): moms err. (norm 2), triang. vs compr.
%        M(new) : cardinality polygauss via triangulation
%        M^(new): cardinality polygauss via triangulation+compr.
%        ratio  : M(new)/M^(new)
%        ------------------------------------------------------
%  	 ** NUMBER OF TRIANGLES:   504
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
%% Date: December 5, 2018
%--------------------------------------------------------------------------


adeV=5:5:35; % ALGEBRAIC DEGREE OF PRECISION.
pos=2;


%--------------------------------------------------------------------------
%                        MAIN CODE STARTS HERE
%--------------------------------------------------------------------------

warning off;
if nargin < 1
    example=1;
end

[xv,yv,iv]=define_domain(example);
vertices=[xv yv];

%--------------------------------------------------------------------------
% DEFINE RULES.
%--------------------------------------------------------------------------

cards=[];
momerrV=[];
ratio=[];

fprintf('\n       =====================================');
switch example
    case 1
        fprintf('\n             Domain: Omega(QF1)   ');
    case 2
        fprintf('\n             Domain: Omega(QF1)   ');
end
fprintf('\n       =====================================');

fprintf('\n        delta E(T,TC) M(new) M^(new) ratio');
fprintf('\n       -------------------------------------');
for ade=adeV
    
    % New polygauss and its compressed version.
    [xyw,xvc,yvc,pgon,tri]=polygauss_2018(ade,xv,yv,iv);
    [xyc,wc,momerrL] = comprexcub(ade,xyw(:,1:2),xyw(:,3),pos);
    
    momerrV=[momerrV; momerrL];
    
    % Cardinalities and ratios
    cards=[cards; size(xyw,1) size(xyc,1)];
    ratio=[ratio; size(xyw,1)/size(xyc,1)];
    
    fprintf('\n \t %2.0f   %1.0e  %6.0f %6.0f   %5.1f',ade,momerrV(end),...
        size(xyw,1),size(xyc,1),ratio(end));
    
end


fprintf('\n       ------------------------------------------------------');
fprintf('\n       *** LEGEND ***');
fprintf('\n       delta  : algebraic degree of precision');
fprintf('\n       E(T,TC): moms err. (norm 2), triang. vs compr.');
fprintf('\n       M(new) : cardinality polygauss via triangulation');
fprintf('\n       M^(new): cardinality polygauss via triangulation+compr.');
fprintf('\n       ratio  : M(new)/M^(new)');
fprintf('\n       ------------------------------------------------------');


fprintf('\n \t ** NUMBER OF TRIANGLES: %5.0f',size(tri,1))

fprintf('\n \n');


clf;
axis off;
plot(pgon,'FaceColor',[0.5 0.5 0.5]);




function [xv,yv,iv]=define_domain(domain)

%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% IN THIS DEMO WE COMPARE, IF POSSIBLE, THE TRIANGULATIONS OF SEVERAL
% POLYGONS, USING polyshape CLASS.
%--------------------------------------------------------------------------
% INPUT.
%--------------------------------------------------------------------------
% domain: 1: quatrefoil with M=129, 2: quatrefoil with M=512.
%
% Important note: choose M=13+4*k.
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
%% Date: December 05, 2018
%--------------------------------------------------------------------------

switch domain
    case 1
        % fprintf('\n \t [POLYGON]: QUATERFOIL LIKE');
        M=129;
        th=linspace(0,2*pi,M);
        %th=(th(1:end-1))';
        polygon_sides=[cos(th').*(sin(2*th')) sin(th').*(sin(2*th'))];
        polygon_sides=polygon_sides(1:end-1,:);
        
    case 2
        % fprintf('\n \t [POLYGON]: QUATERFOIL LIKE');
        spec_settings=513;
        
        th=linspace(0,2*pi,spec_settings);
        %th=(th(1:end-1))';
        polygon_sides=[cos(th').*(sin(2*th')) sin(th').*(sin(2*th'))];
        polygon_sides=polygon_sides(1:end-1,:);
end

% adjusting origin matching.
toll=10^(-14);
for ii=1:size(polygon_sides,1)
    if norm(polygon_sides(ii,:)) < toll
        polygon_sides(ii,:)=[0 0];
    end
end

xv=polygon_sides(:,1); yv=polygon_sides(:,2); iv=length(xv);

