
function demo_polygauss_cardinality(example)

%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% IN THIS DEMO WE COMPARE, IF POSSIBLE, THE RESULTS OBTAINED BY THE OLD
% POLYGAUSS (i.e. polygauss_2017) AND THE NEW ONE (i.e. polygauss_2018).
% THE OLD CODES COULD INTEGRATE ON QUADRILATERALS OR VIA TRIANGULATION ON
% SIMPLE POLYGONS, PROVIDING A P.I. RULE.
% THE CODE DISPLAYS RELATIVE ERRORS, CARDINALITY AND CPUTIMES IF THE CODES
% CAN DETERMINE THE INTEGRALS USING P.I. RULES.
% AS NUMERICAL TEST WE USE A NON SYMMETRIC POLYNOMIAL OF DEGREE ade.
%--------------------------------------------------------------------------
% IMPORTANT:
%--------------------------------------------------------------------------
% * THE CODE REQUIRES THE TOOLBOX CONTAINING POLYSHAPE (AVAILABLE FROM
% MATLAB R2017b VERSIONS).
% * TEST EXAMPLES COMES FROM THE FILE mydomain.m
% * COMPRESSION COMES FROM comprexcub.m
% * COMPARISONS WITH OLDER VERSIONS REQUIRE polygauss_2017.m
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% example: 1: convex, 6 sides. 
%          2: non convex, 9 sides.
%--------------------------------------------------------------------------
% ROUTINES:
%--------------------------------------------------------------------------
% 1. define_domain (attached)
% 2. polygauss_2017 (external)
% 3. polygauss_2018 (external)
% 4. comprexcub (external)
%--------------------------------------------------------------------------
% EXAMPLE:
%--------------------------------------------------------------------------
% >>demo_polygauss_cardinality(1)
% 
%        ==============================
%              Domain: Omega(C)   
%        ==============================
%  
%         CARDINALITIES OF THE RULES 
% 
%         delta  M(old) M(new) M^(new)
%        ------------------------------
%  	  5    180     28     21 
%  	 10    660     96     66 
%  	 15   1440    184    136 
%  	 20   2520    312    231 
%  	 25   3900    468    351 
%  	 30   5580    684    496 
%  	 35   7560    912    666 
%  	 40   9840   1180    861 
%        ------------------------------------------------------
%        *** LEGEND ***
%        delta: algebraic degree of precision
%        M(old) : cardinality polygauss via Gauss-Green
%        M(new) : cardinality polygauss via triangulation
%        M^(new): cardinality polygauss via triangulation+compr.
%        ------------------------------------------------------
%  
% 
%          MOMENTS ERRORS (NORM 2) 
% 
%         delta  E(T,TC)  E(T,GG)
%        ------------------------------
%  	  5     3e-16    4e-16 
%  	 10     2e-16    5e-15 
%  	 15     3e-16    2e-15 
%  	 20     6e-16    4e-15 
%  	 25     5e-16    4e-14 
%  	 30     1e-15    1e-14 
%  	 35     9e-16    4e-14 
%  	 40     1e-15    8e-14 
%        ------------------------------------------------------
%        *** LEGEND ***
%        delta: algebraic degree of precision
%        E(T,TC) : triang. vs compr.
%        E(T,GG) : triang. vs Gauss-Green
%        ------------------------------------------------------
%  
% >> demo_polygauss_cardinality(2)
% 
%        ==============================
%              Domain: Omega(NC)   
%        ==============================
%  
%         CARDINALITIES OF THE RULES 
% 
%         delta  M(old) M(new) M^(new)
%        ------------------------------
%  	  5    235     49     21 
%  	 10    870    168     66 
%  	 15   1905    322    136 
%  	 20   3340    546    231 
%  	 25   5175    819    351 *
%  	 30   7410   1197    496 
%  	 35  10045   1596    666 
%  	 40  13080   2065    861 
%        ------------------------------------------------------
%        *** LEGEND ***
%        delta: algebraic degree of precision
%        M(old) : cardinality polygauss via Gauss-Green
%        M(new) : cardinality polygauss via triangulation
%        M^(new): cardinality polygauss via triangulation+compr.
%        ------------------------------------------------------
%  
% 
%          MOMENTS ERRORS (NORM 2) 
% 
%         delta  E(T,TC)  E(T,GG)
%        ------------------------------
%  	  5     3e-16    1e-15 
%  	 10     6e-16    1e-15 
%  	 15     8e-16    8e-15 
%  	 20     8e-16    5e-03 
%  	 25     1e-15    1e-11 
%  	 30     1e-15    1e-05 
%  	 35     5e-15    6e-05 
%  	 40     4e-15    1e-03 
%        ------------------------------------------------------
%        *** LEGEND ***
%        delta: algebraic degree of precision
%        E(T,TC) : triang. vs compr.
%        E(T,GG) : triang. vs Gauss-Green
%        ------------------------------------------------------
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


adeV=5:5:40; % ALGEBRAIC DEGREE OF PRECISION.
pos=2; %0: QR, 1:lsqnonneg. 2: fast Lawson-Hanson.


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

fprintf('\n       ==============================');
switch example
    case 1
        fprintf('\n             Domain: Omega(C)   ');
    case 2
        fprintf('\n             Domain: Omega(NC)   ');
end
fprintf('\n       ==============================');

fprintf('\n \n        CARDINALITIES OF THE RULES \n');

fprintf('\n        delta  M(old) M(new) M^(new)');
fprintf('\n       ------------------------------');
for ade=adeV
    
    % New polygauss and its compressed version.
    [xyw,xvc,yvc,pgon,tri]=polygauss_2018(ade,xv,yv,iv);
    [xyc,wc,momerrL] = comprexcub(ade,xyw(:,1:2),xyw(:,3),pos);
    
    % Old polygauss and its compressed version.
    xyw17=polygauss_2017(vertices,ade,'quadrangulation');
    [xyc_17,wc_17,momerrL17] = comprexcub(ade,xyw17(:,1:2),xyw17(:,3),pos);
    
    
    momerrV=[momerrV; momerrL momerrL17];
    
    % Cardinalities
    cardsL=[size(xyw17,1) size(xyw,1) size(xyc,1)];
    cards=[cards; cardsL];
    
    fprintf('\n \t %2.0f %6.0f %6.0f %6.0f ',ade,cardsL(1),...
        cardsL(2),cardsL(3));
    
end

fprintf('\n       ------------------------------------------------------');
fprintf('\n       *** LEGEND ***');
fprintf('\n       delta: algebraic degree of precision');
fprintf('\n       M(old) : cardinality polygauss via Gauss-Green');
fprintf('\n       M(new) : cardinality polygauss via triangulation');
fprintf('\n       M^(new): cardinality polygauss via triangulation+compr.');
fprintf('\n       ------------------------------------------------------');

fprintf('\n \n');

fprintf('\n         MOMENTS ERRORS (NORM 2) \n');
fprintf('\n        delta  E(T,TC)  E(T,GG)');
fprintf('\n       ------------------------------');
for ii=1:length(adeV)
    ade=adeV(ii);
    
    fprintf('\n \t %2.0f     %1.0e    %1.0e ',ade,momerrV(ii,1),...
        momerrV(ii,2));
    
end

fprintf('\n       ------------------------------------------------------');
fprintf('\n       *** LEGEND ***');
fprintf('\n       delta: algebraic degree of precision');
fprintf('\n       E(T,TC) : triang. vs compr.');
fprintf('\n       E(T,GG) : triang. vs Gauss-Green');
fprintf('\n       ------------------------------------------------------');

fprintf('\n \n');



function [xv,yv,iv]=define_domain(domain)

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

switch domain
    case 1
        % fprintf('\n \t [POLYGON]: 6 SIDES CONVEX POLYGON');
        polygon_sides=[0.1 0; 0.7 0.2; 1 0.5; 0.75 0.85; 0.5 1; 0 0.25; 0.1 0];
        
    case 2
        % fprintf('\n \t [POLYGON]: 9 SIDES NON CONVEX POLYGON');
        polygon_sides=(1/4)*[1 2; 1 0; 3 2; 3 0; 4 2; 3 3; 3 0.85*4; 2 4;
            0 3; 1 2];
end

xv=polygon_sides(:,1); yv=polygon_sides(:,2); iv=length(xv);

