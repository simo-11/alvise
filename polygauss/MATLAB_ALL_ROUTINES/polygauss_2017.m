
function xyw=polygauss_2017(polygon_vertices,degree,method)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine computes a quadrature rule over a polygon with a prescribed
% degree of precision.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% polygon_vertices:
% The variable "polygon_vertices" is a matrix N x 2 of N points. 
% First and last row of polygon_vertices must be equal.
%
% degree:
% Algebraic degree of precision of the rule.
%
% method:
% It is a string. Alternatives are 
% 'triangulation': formula based on minimal triangulation of the domain.
% 'quadrangulation': formula based on a quadrangulation of the domain.
% 'quadrangulation_allin': decomposition of the domain in convex subdomains 
%                  where it uses the classical polygauss with quadrangulations.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% xyw: The cubature rule is stored in this variable.
%      Nodes are (xyw(:,1),xyw(:,2)). Weights are stored in xyw(:,3).
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% NOTE:
%--------------------------------------------------------------------------
% First and last row of "polygon_vertices" must be equal.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Copyright (C) 2007-2016 Alvise Sommariva, Marco Vianello.
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
%% Author:  Federica Basaglia
%%          Alvise Sommariva <alvise@euler.math.unipd.it>
%%          Marco Vianello   <marcov@euler.math.unipd.it>
%%
%% Date: August 09, 2016
%--------------------------------------------------------------------------

warning off;

if nargin < 3
    method = 'triangulation';
end

validStrings = {'triangulation','quadrangulation','quadrangulation_allin'};

if ( any(strcmpi(method, validStrings)) )
    
    if ( strcmpi(method, 'triangulation') )
        %fprintf('\n \t polygon_cubature_triangulation');
        [xyw,polygon_vertices]=polygon_cubature_triangulation(...
            polygon_vertices,degree);
    end
    
    if ( strcmpi(method, 'quadrangulation') )
        % fprintf('\n \t polygauss_2013');
        if norm(polygon_vertices(1,:)-polygon_vertices(end,:)) > 0
            polygon_vertices=[polygon_vertices; polygon_vertices(1,:)];
        end
        [nodes_x, nodes_y, weights]=polygauss_2013(degree,polygon_vertices);
        % P=[0 0.75]; Q=[1 0.5]; [nodes_x, nodes_y, weights]=polygauss_2013(degree,polygon_vertices,2,P,Q);
        xyw=[nodes_x nodes_y weights];
    end
    
    if ( strcmpi(method, 'quadrangulation_allin') )
        % fprintf('\n \t polygauss_better');
        [nodes_x, nodes_y, weights,vertices_sequence, ...
            vertices_sequence_pointer]=polygauss_quad_allin(degree,polygon_vertices);
        xyw=[nodes_x nodes_y weights];
    end
    
    
else
    warning('wrong string as method, choosen better polygauss');
    tw=trigauss(n,alpha,beta,'classic');
end















function [xyw,polygon_vertices]=polygon_cubature_triangulation(...
    polygon_vertices,degree)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine computes a quadrature rule over a polygon via triangulation,
% with a prescribed degree of precision.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
%
% polygon_vertices: Polygon vertices written in cartesian coordinates,
%                   counterclockwise.
%                   It is a N x 2 matrix, where N is the number of nodes.
%
% degree          : degree of precision of the cubature rule.
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
%
% xyw             : Cubature rule \sum_{k=1}^M w_k f(x_k,y_k) stored as
%                   M x 3 matrix, in which the k-th row is [x_k y_k w_k].
% polygon_vertices_closed: Polygon vertices written in cartesian
%                   coordinates, counterclockwise.
%                   It is a N x 2 matrix, where N is the number of nodes.
%                   First and last row are equal.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% ROUTINES CALLED:
%--------------------------------------------------------------------------
% 1. minimal_triangulation
% 2. stroud_conical_rules
% 3. polygauss_2014
%--------------------------------------------------------------------------

N=length(polygon_vertices(:,1)); % NUMBER OF VERTICES.

% COMPUTING TRIANGULATION
XV=polygon_vertices(:,1); YV=polygon_vertices(:,2);
[TRI,err]=minimal_triangulation(XV,YV);
N_subtriangles=size(TRI,1);

% DETERMINING CUBATURE RULE.
xyw=[];


for ii=1:N_subtriangles
    TRI_ii=TRI(ii,:); TRI_ii=TRI_ii';
    XV_ii=XV(TRI_ii); YV_ii=YV(TRI_ii);
    xyw_ii=stroud_conical_rules(degree,[XV_ii YV_ii]);
    
    % VERTICES MUST NOT BE REPEATED.
    xyw=[xyw; xyw_ii];
    
end

















function [tri,err]=minimal_triangulation(x,y)

%--------------------------------------------------------------------------
% OBJECT.
%----------
%
% COMPUTE MINIMAL TRIANGULATION (BASIC ALGORITHM) ON SIMPLE DOMAINS,
% I.E. WITHOUT SELF-INTERSECTIONS.
%
%--------------------------------------------------------------------------
% INPUT.
%----------
%
% x,y: VERTICES OF THE POLYGONS. THE FIRST VERTEX IS NOT REPEATED.
%      THEY MUST BE COLUMN VECTORS.
%
%--------------------------------------------------------------------------
% OUTPUT.
%----------
%
% tri: MATRIX IN WHICH EACH ROW IS FORMED BY 3 POSITIVE INTEGERS IN WHICH
%      EACH ONE REPRESENTS A POINTER TO A VERTEX. EXAMPLE: [1 2 3; 1 3 4]
%      DESCRIBE TWO TRIANGLES T1=[x(1) y(1); x(2) y(2); x(3) y(3)] AND
%      T2=[x(1) y(1); x(3) y(3); x(4) y(4)].
%
%--------------------------------------------------------------------------
% FUNCTIONS USED IN THIS ROUTINE.
%---------------------------------
%
% 1. find_ear (USED BY "minimal_triangulation")
% 3. polyarea (MATLAB BUILT_IN FUNCTION, USED BY "minimal_triangulation").
% 4. setdiff (MATLAB BUILT_IN FUNCTION, USED BY "minimal_triangulation").
% 6. verLessThan (MATLAB BUILT_IN FUNCTION, USED BY "find_ear").
% 7. inpolygon (MATLAB BUILT_IN FUNCTION, USED BY "find_ear").
%
% THE FUNCTION "find_ear" IS ATTACHED TO THIS FILE.
%
% THE FUNCTIONS "polyarea", "setdiff", "verLessThan", "inpolygon" ARE
% MATLAB BUILT-IN.
%
%--------------------------------------------------------------------------
% TEST.
%-------
%
% THIS ROUTINE HAS BEEN TESTED IN MATLAB 7.6.0.324 (R2008a).
%
%--------------------------------------------------------------------------
% REFERENCES.
%-------------
%
% [1] Federica Basaglia, Un nuovo metodo di cubatura su poligoni (in
%     italian)
%
% [2] Chazelle and J. Incerpi, Triangulation and shape complexity,
%     ACM Trans. on Graphics, vol. 3, pp. 135-152, 1984.
%
% [3] B. Chazelle, Triangulating a simple polygon in linear time,
%     Discrete Comput. Geom., vol. 6, pp. 485-524, 1991.
%
% [4] David Eberly, Triangulation by Ear Clipping,
%           http://www.geometrictools.com/
%
% [5] R. Seidel, A simple and fast incremental randomized algorithm for
%     computing trapezoidal decompositions and for triangulating polygons,
%     Computational Geometry: Theory and Applications, vol. 1, no.1,
%     pp. 51-64, 1991.
%
% [6] Subhash Suri, Polygon Triangulation
%
% [7] http://cgm.cs.mcgill.ca/~godfried/teaching/cg-projects/97/
%            Ian/cutting_ears.html
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Copyright (C) 2007-2009 Alvise Sommariva, Marco Vianello.
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
%% Author:  Federica Basaglia
%%          Alvise Sommariva <alvise@euler.math.unipd.it>
%%          Marco Vianello   <marcov@euler.math.unipd.it>
%%
%% Date: November 04, 2010
%--------------------------------------------------------------------------

err=0;


pts=[x y];
number_pts=length(x);

% if norm(pts(1,:)-pts(end,:)) == 0
%     fprintf('\n \t reps');
%     x=x(1:end-1,:); y=y(1:end-1,:);
%     pts=[x y];
%     number_pts=length(x);
% end


pts_remaining=pts;
pointer_pts_remaining=(1:number_pts);
tri=[];

while length(pts_remaining(:,1)) > 3
    
    [tri_loc,err]=find_ear(pts_remaining);
    
    if (err == 1)
        tri=[];
        fprintf('\n \t [WARNING]: THE FUNCTION IS NOT ABLE');
        fprintf(' TO COMPUTE THE TRIANGULATION');
        return;
    end
    
    tri_add=pointer_pts_remaining(tri_loc);
    
    tri=[tri; tri_add];
    pointer_pt_remove=tri_add(2);
    pointer_pts_remaining=setdiff(pointer_pts_remaining,pointer_pt_remove);
    pts_remaining=pts(pointer_pts_remaining,:);
    
end

tri=[tri; pointer_pts_remaining];

% CUTTING DEGENERATE TRIANGLES.
vv=[1 number_pts+1];
index=[];
for ii=1:size(tri,1)
    
    ints=intersect(tri(ii,:),vv);
    if length(ints) <= 1
        index=[index; ii];
    end
end

tri=tri(index,:);




%--------------------------------
% find ear.
%--------------------------------
function [tri,err]=find_ear(vertices_pts)

%--------------------------------------------------------------------------
% OBJECT.
%----------
%
% FIND ONE EAR IN A POLYGON DESCRIBED BY THE ORDER SEQUENCE vertices_pts,
% BY A BASIC ALGORITHM, TO COMPUTE A MINIMAL TRIANGULATION (SEE CHAZELLE -
% SEIDEL ALGORITHM FOR A FASTER IMPLEMENTATION). IF THE ROUTINE IS NOT
% COMPLETED, THEN err=1.
%
% PS. WE USED A BASIC TRIANGULATION ALGORITHM FOR THIS PURPOSE. A FASTER
% ONE HAS BEEN STUDIED BY CHAZELLE AND SEIDEL. OUR ALGORITHM WORKS IN
% POLYGONS THAT ARE "SIMPLE", I.E. WITHOUT SELF-INTERSECTIONS. SEE THE
% REFERENCES ABOVE FOR MORE DETAILS.
%
%--------------------------------------------------------------------------
% INPUTS.
%----------
%
% vertices_pts: vertices of the polygon. Two column matrix.
%
%--------------------------------------------------------------------------
% OUTPUT.
%----------
%
% tri: vertices of the triangle that is "ear".
%
% err: "0" means "no warning", "1" means "the ear has not been found".
%
%--------------------------------------------------------------------------

number_pts=size(vertices_pts,1);
vertices_index=[ (1:number_pts)'; 1];
pts_remaining=vertices_pts;
curr_starting_index=1;
ear_found=0;
err=0;

while (ear_found == 0) & (curr_starting_index < number_pts)
    
    trial_triangle_vertices_pointer=...
        vertices_index(curr_starting_index:curr_starting_index+2);
    
    % TRIANGLE TO TEST.
    pts_trial_triangle=vertices_pts(trial_triangle_vertices_pointer,:);
    
    pts_remaining_pointer=...
        setdiff(1:number_pts,curr_starting_index:curr_starting_index+2);
    
    pts_remaining=vertices_pts(pts_remaining_pointer,:);
    
    xv=[vertices_pts(trial_triangle_vertices_pointer,1); ...
        vertices_pts(curr_starting_index,1)];
    yv=[vertices_pts(trial_triangle_vertices_pointer,2); ...
        vertices_pts(curr_starting_index,2)];
    
    % CHECKING POSITION OF BARYCENTER (TO SEE IF IT IS AN EAR).
    barycenter_triangle=...
        [sum(vertices_pts(trial_triangle_vertices_pointer,1)) ...
        sum(vertices_pts(trial_triangle_vertices_pointer,2))]/3;
    
    xt=pts_remaining(:,1);
    yt=pts_remaining(:,2);
    
    in = inpolygon(xt,yt,xv,yv);
    
    barycenter_triangle=[sum(vertices_pts...
        (trial_triangle_vertices_pointer,1)) ...
        sum(vertices_pts(trial_triangle_vertices_pointer,2))]/3;
    
    
    % inpolygon has been modified by MATLAB. We try to patch this problem.
    
    % local_version=verLessThan('matlab', '7.6.1');
    local_version=0;
    if local_version == 0
        
        %******************************************************************
        %
        % IN = inpolygon(X,Y,xv,yv) returns a matrix IN the
        %      same size as X and Y.
        %      Each element of IN is assigned the value 1 or 0 depending on
        %      whether the point (X(p,q),Y(p,q)) is inside the polygonal
        %     region whose vertices are specified by the vectors xv and yv.
        %
        % In particular:
        %
        % IN(p,q) = 1: if (X(p,q),Y(p,q)) is inside the polygonal region
        %              or on the polygon boundary
        % IN(p,q) = 0: if (X(p,q),Y(p,q)) is outside the polygonal region
        %
        % [IN ON] = inpolygon(X,Y,xv,yv) returns a second matrix ON the
        %           same size as X and Y.
        %           Each element of ON is assigned the value 1 or 0
        %           depending on whether the point (X(p,q),Y(p,q)) is on
        %           the boundary of the polygonal region whose vertices are
        %           specified by the vectors xv and yv.
        %
        % In particular:
        %
        % ON(p,q) = 1: if (X(p,q),Y(p,q)) is on the polygon boundary
        % ON(p,q) = 0: if (X(p,q),Y(p,q)) is inside or outside the polygon
        %              boundary.
        %
        %******************************************************************
        
        % Working in Matlab 7.6.
        [bar_in,bar_on]= inpolygon(barycenter_triangle(1),...
            barycenter_triangle(2),vertices_pts(vertices_index,1),...
            vertices_pts(vertices_index,2));
        
        if (sum(in(1:end-1)) == 0) & (bar_in == 1) & (bar_on == 0)
            ear_found=1;
            tri=curr_starting_index:curr_starting_index+2;
            return;
        end
        
    else
        
        %******************************************************************
        %
        % IN = INPOLYGON(X, Y, XV, YV) returns a matrix IN the size of
        % X and Y.  IN(p,q) = 1 if the point (X(p,q), Y(p,q)) is
        % strictly inside the polygonal region whose vertices are
        % specified by the vectors XV and YV;  IN(p,q) is 0.5 if
        % the point is on the polygon; otherwise IN(p,q) = 0.
        %
        %******************************************************************
        
        % Working in Matlab 6.1.
        [bar_in]= inpolygon(barycenter_triangle(1),...
            barycenter_triangle(2),vertices_pts(vertices_index,1),...
            vertices_pts(vertices_index,2));
        
        if (sum(in(1:end-1)) == 0) & (bar_in > 0.5)
            ear_found=1;
            tri=curr_starting_index:curr_starting_index+2;
            return;
        end
        
    end
    
    curr_starting_index=curr_starting_index+1;
end

err=1;





function [nodes_x, nodes_y, weights,vertices_sequence, ...
    vertices_sequence_pointer]=polygauss_quad_allin(N,polygon_vertices)


%--------------------------------------------------------------------------
% INPUT.
%--------
% [N]      : DEGREE OF PRECISION OF THE ALGEBRAIC CUBATURE RULE.
%
% [polygon_sides]: IF THE POLYGON HAS "L" SIDES, "boundary.pts" IS A
%            VARIABLE CONTAINING ITS VERTICES, ORDERED COUNTERCLOCKWISE. AS
%            LAST ROW MUST HAVE THE COMPONENTS OF THE FIRST VERTEX.
%            IN OTHER WORDS, THE FIRST ROW AND LAST ROW ARE EQUAL.
%            "polygon_sides" IS A "L+1 x 2" MATRIX.
%
%--------------------------------------------------------------------------
% OUTPUT.
%--------
% [nodes_x, nodes_y]: THE CUBATURE RULE PRODUCES THE NODES
%            "(nodes_x,nodes_y)". IF THE CUBATURE RULE GENERATES "K" NODES,
%             "nodes_x" AND "nodes_y" ARE "K x 1" COLUMN VECTORS.
%
% [weights]: cubature weights.
%
%--------------------------------------------------------------------------
% SUBROUTINES.
%-------------
% 1. convex_decomposition (external).
% 2. polygauss_2013 (internal).
%
%--------------------------------------------------------------------------
% CODE FEATURES.
%----------------
% THIS CODE FIRST PERFORMS A DECOMPOSITION OF THE SIMPLE POLYGON HAVING
% VERTICES polygon_vertices INTO CONVEX SUBPOLYGONS "R_k" AND THEN APPLY
% OUR splinegauss ROUTINE INTO "R_k" OBTAINING A RULE WITH POSITIVE WEIGHTS
% AND DEGREE OF PRECISION "N".
% THE MAIN DIFFERENCE WITH RESPECT TO splinegauss IS THAT ALL THE NODES ARE
% INTERNAL TO THE DOMAIN DEFINED BY polygon_vertices.
%
%--------------------------------------------------------------------------
%% Copyright (C) 2011 Alvise Sommariva, Marco Vianello.
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
%% Date: April 01, 2011 - April 11, 2011.
%--------------------------------------------------------------------------

% Try first polygauss_2013. If some points are outside the domain, we
% choose an alternative version with more points, but all inside the
% domain.


% xyw=polygauss_2013(N,polygon_vertices);
% nodes_x=xyw(:,1); nodes_y=xyw(:,2); weights=xyw(:,3);
%
% Xv=polygon_vertices(:,1);
% Yv=polygon_vertices(:,2);
%
% in2013=inpolygon(nodes_x,nodes_y,Xv,Yv);
% out2013=length(nodes_x)-sum(in2013);

[convex_flag,nonconvex_vertex]=is_convex(polygon_vertices(1:end-1,:));

% fprintf('\n \t CONVEX: %1.0f',convex_flag)

if convex_flag == 0
    [vertices_sequence, vertices_sequence_pointer]=convex_decomposition...
        (polygon_vertices);
    [nodes_x,nodes_y,weights]=polygauss_2013_in(N,polygon_vertices);
    % [nodes_x, nodes_y, weights]=polygauss_2013(N,polygon_vertices);
else
    vertices_sequence=[]; vertices_sequence_pointer=[];
    
    [nodes_x, nodes_y, weights]=polygauss_2013(N,polygon_vertices);
end

if nargout == 1
    nodes_x=[nodes_x nodes_y weights];
end






function [nodes_x, nodes_y, weights]=polygauss_2013_in(N,polygon_vertices)


%--------------------------------------------------------------------------
% INPUT.
%--------
% [N]      : DEGREE OF PRECISION OF THE ALGEBRAIC CUBATURE RULE.
%
% [polygon_sides]: IF THE POLYGON HAS "L" SIDES, "boundary.pts" IS A
%            VARIABLE CONTAINING ITS VERTICES, ORDERED COUNTERCLOCKWISE. AS
%            LAST ROW MUST HAVE THE COMPONENTS OF THE FIRST VERTEX.
%            IN OTHER WORDS, THE FIRST ROW AND LAST ROW ARE EQUAL.
%            "polygon_sides" IS A "L+1 x 2" MATRIX.
%
%--------------------------------------------------------------------------
% OUTPUT.
%--------
% [nodes_x, nodes_y]: THE CUBATURE RULE PRODUCES THE NODES
%            "(nodes_x,nodes_y)". IF THE CUBATURE RULE GENERATES "K" NODES,
%             "nodes_x" AND "nodes_y" ARE "K x 1" COLUMN VECTORS.
%
% [weights]: cubature weights.
%
%--------------------------------------------------------------------------
% SUBROUTINES.
%-------------
% 1. convex_decomposition (external).
% 2. polygauss_2013 (internal).
%
%--------------------------------------------------------------------------
% CODE FEATURES.
%----------------
% THIS CODE FIRST PERFORMS A DECOMPOSITION OF THE SIMPLE POLYGON HAVING
% VERTICES polygon_vertices INTO CONVEX SUBPOLYGONS "R_k" AND THEN APPLY
% OUR splinegauss ROUTINE INTO "R_k" OBTAINING A RULE WITH POSITIVE WEIGHTS
% AND DEGREE OF PRECISION "N".
% THE MAIN DIFFERENCE WITH RESPECT TO splinegauss IS THAT ALL THE NODES ARE
% INTERNAL TO THE DOMAIN DEFINED BY polygon_vertices.
%
%--------------------------------------------------------------------------
%% Copyright (C) 2011 Alvise Sommariva, Marco Vianello.
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
%% Date: April 01, 2011 - April 11, 2011.
%--------------------------------------------------------------------------

[vertices_sequence, vertices_sequence_pointer]=convex_decomposition...
    (polygon_vertices);

nodes_x=[]; nodes_y=[]; weights=[];

i_end=0;

for ii=1:length(vertices_sequence_pointer)
    
    i_init=i_end+1;
    i_end=i_init+vertices_sequence_pointer(ii)-1;
    vertices_loc=[vertices_sequence(i_init:i_end,:); ...
        vertices_sequence(i_init,:)];
    [nodes_x_loc,nodes_y_loc,weights_loc]=polygauss_2013(N,vertices_loc);
    
    nodes_x=[nodes_x; nodes_x_loc];
    nodes_y=[nodes_y; nodes_y_loc];
    weights=[weights; weights_loc];
    
end


function [nodes_x, nodes_y, weights]=polygauss_2013(N,polygon_sides,rotation,P,Q)

%--------------------------------------------------------------------------
% REFERENCE PAPER:
% [1] A. SOMMARIVA and M. VIANELLO
% "Gauss-like and triangulation-free cubature over polygons".
%
% INPUT:
%
% N     : DEGREE OF THE 1 DIMENSIONAL GAUSS-LEGENDRE RULE.
%
% polygon_sides: IF THE POLYGON HAS "L" SIDES, "boundary.pts" IS A
%         VARIABLE CONTAINING ITS VERTICES, ORDERED COUNTERCLOCKWISE.
%         AS LAST ROW MUST HAVE THE COMPONENTS OF THE FIRST VERTEX.
%         IN OTHER WORDS, THE FIRST ROW AND LAST ROW ARE EQUAL.
%         "polygon_sides" IS A "L+1 x 2" MATRIX.
%
%            --------- NOT MANDATORY VARIABLES ---------
%
% rotation: 0: NO ROTATION.
%           1: AUTOMATIC.
%           2: PREFERRED DIRECTION ROTATION BY P, Q.
%
% P, Q: DIRECTION THAT FIXES THE ROTATION.
%
% OUTPUT:
%
% nodes_x, nodes_y, weights : THE GAUSS LIKE FORMULA PRODUCES THE NODES AND
%          THE WEIGHTS OFA CUBATURE RULE ON THE POLYGON.
%
%--------------------------------------------------------------------------
% EXAMPLE 1 (NO ROTATION.)
%---------------------------
%
% >> xyw=polygauss_2013(2,[0 0; 1 0; 1 1; 0 1; 0 0],0)
%
% xyw =
%
%     0.2113    0.2113    0.2500
%     0.2113    0.7887    0.2500
%     0.7887    0.2113    0.2500
%     0.7887    0.7887    0.2500
%
% >>
%
%--------------------------------------------------------------------------
% EXAMPLE 2 (AUTO ROTATION.)
%-----------------------------
%
% >> xyw=polygauss_2013(2,[0 0; 1 0; 1 1; 0 1; 0 0])
%
% xyw =
%
%     0.0683    0.0444    0.0078
%     0.3028    0.1972    0.0556
%     0.5374    0.3499    0.0616
%     0.6501    0.4626    0.0616
%     0.8028    0.6972    0.0556
%     0.9556    0.9317    0.0078
%     0.9317    0.9556    0.0078
%     0.6972    0.8028    0.0556
%     0.4626    0.6501    0.0616
%     0.3499    0.5374    0.0616
%     0.1972    0.3028    0.0556
%     0.0444    0.0683    0.0078
%     0.1008    0.0119    0.0078
%     0.4472    0.0528    0.0556
%     0.7935    0.0938    0.0616
%     0.9062    0.2065    0.0616
%     0.9472    0.5528    0.0556
%     0.9881    0.8992    0.0078
%     0.8992    0.9881    0.0078
%     0.5528    0.9472    0.0556
%     0.2065    0.9062    0.0616
%     0.0938    0.7935    0.0616
%     0.0528    0.4472    0.0556
%     0.0119    0.1008    0.0078
%
% >>

%--------------------------------------------------------------------------
%% Copyright (C) 2007-2013 Marco Vianello and Alvise Sommariva
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

%% Authors:
%% Marco Vianello    <marcov@euler.math.unipd.it>
%% Alvise Sommariva  <alvise@euler.math.unipd.it>
%% Date: April 30, 2013.
%--------------------------------------------------------------------------

%----------------------------------------------------------------------
% BOUNDARY PTS.
%----------------------------------------------------------------------
x_bd=polygon_sides(:,1);
y_bd=polygon_sides(:,2);

%----------------------------------------------------------------------
% "MINIMUM" RECTANGLE CONTAINING POLYGON.
%----------------------------------------------------------------------
x_min=min(x_bd); x_max=max(x_bd);
y_min=min(y_bd); y_max=max(y_bd);

%----------------------------------------------------------------------
% SOME AUTOMATIC SETTINGS.
%----------------------------------------------------------------------
if nargin < 3
    rotation=1;
end
cubature_type=4;

%--------------------------------------------------------------------------
% POLYGON ROTATION (IF NECESSARY).
%--------------------------------------------------------------------------
switch rotation
    case 0
        %         fprintf('\n \t [ROTATION]: NO.');
        rot_matrix=eye(2);
        axis_abscissa=[x_min y_max]-[x_min y_min];
    case 1
        %         fprintf('\n \t [ROTATION]: AUTOMATIC');
        [polygon_sides,rot_matrix,rot_angle,axis_abscissa,P,Q]=...
            auto_rotation(polygon_sides,[],[]);
        %         fprintf(' [ANGLE CLOCKWISE (RESPECT Y, IN DEGREES)]: %5.5f',...
        %             rot_angle*180/pi);
    case 2
        %         fprintf('\n \t [ROTATION]: PREFERRED DIRECTION');
        nrm_vect=norm(Q-P);
        if (nrm_vect > 0)
            direction_axis=(Q-P)/nrm_vect;
            [polygon_sides,rot_matrix,rot_angle,axis_abscissa,P,Q]=...
                auto_rotation(polygon_sides,P,Q);
        else
            %             fprintf('\n \t [WARNING]: THE DIRECTION VECTOR IS NULL. ')
            %             fprintf('USING AUTOMATIC ROTATION.');
            [polygon_sides,rot_matrix,rot_angle,axis_abscissa,P,Q]=...
                auto_rotation(polygon_sides,P,Q);
        end
        %         fprintf(' [ANGLE CLOCKWISE (RESPECT Y)]: %5.5f',rot_angle*180/pi);
end

%--------------------------------------------------------------------------
% COMPUTE NODES AND WEIGHTS OF 1D GAUSS-LEGENDRE RULE.
% TAKEN FROM TREFETHEN PAPER "Is ... Clenshaw-Curtis?".
%--------------------------------------------------------------------------

% DEGREE "N".
[s_N,w_N]=cubature_rules_1D((N-1),cubature_type);
N_length=length(s_N);

% DEGREE "M".
M=N+1;
[s_M,w_M]=cubature_rules_1D((M-1),cubature_type);

%----------------------------------------------------------------------
% L: NUMBER OF SIDES OF THE POLYGON.
% M: ORDER GAUSS INTEGRATION.
% N: ORDER GAUSS PRIMITIVE.
%----------------------------------------------------------------------
L=length(polygon_sides(:,1))-1;

%a=0.5;
a=axis_abscissa(1);

%----------------------------------------------------------------------
% COMPUTE 2D NODES (nodes_x,nodes_y) AND WEIGHTS "weights".
%----------------------------------------------------------------------

nodes_x=[];
nodes_y=[];
weights=[];

for index_side=1:L
    x1=polygon_sides(index_side,1); x2=polygon_sides(index_side+1,1);
    y1=polygon_sides(index_side,2); y2=polygon_sides(index_side+1,2);
    if ~(x1 == a & x2 == a)
        if (y2-y1) ~=0
            
            if (x2-x1) ~=0
                s_M_loc=s_M;
                w_M_loc=w_M;
            else
                s_M_loc=s_N;
                w_M_loc=w_N;
            end
            
            M_length=length(s_M_loc);
            
            half_pt_x=(x1+x2)/2; half_pt_y=(y1+y2)/2;
            half_length_x=(x2-x1)/2; half_length_y=(y2-y1)/2;
            
            
            % GAUSSIAN POINTS ON THE SIDE.
            x_gauss_side=half_pt_x+half_length_x*s_M_loc; %SIZE: (M_loc,1)
            y_gauss_side=half_pt_y+half_length_y*s_M_loc; %SIZE: (M_loc,1)
            
            scaling_fact_plus=(x_gauss_side+a)/2; %SIZE: (M_loc,1)
            scaling_fact_minus=(x_gauss_side-a)/2;%SIZE: (M_loc,1)
            
            local_weights=...
                (half_length_y*scaling_fact_minus).*w_M_loc;%SIZE:(M_loc,1)
            
            term_1=repmat(scaling_fact_plus,1,N_length); % SIZE: (M_loc,N)
            
            term_2=repmat(scaling_fact_minus,1,N_length); % SIZE: (M_loc,N)
            
            rep_s_N=repmat(s_N',M_length,1);
            
            % x, y ARE STORED IN MATRICES. A COUPLE WITH THE SAME INDEX
            % IS A POINT, i.e. "P_i=(x(k),y(k))" FOR SOME "k".
            x=term_1+term_2.*rep_s_N;
            y=repmat(y_gauss_side,1,N_length);
            
            number_rows=size(x,1);
            number_cols=size(x,2);
            
            x=x(:); x=x';
            y=y(:); y=y';
            
            rot_gauss_pts=rot_matrix'*[x;y]; % THE INVERSE OF A ROTATION
            % MATRIX IS ITS TRANSPOSE.
            
            x_rot=rot_gauss_pts(1,:); % GAUSS POINTS IN THE ORIGINAL SYSTEM.
            y_rot=rot_gauss_pts(2,:);
            
            x_rot=reshape(x_rot',number_rows,number_cols);
            y_rot=reshape(y_rot',number_rows,number_cols);
            
            nodes_x=[nodes_x; x_rot];
            nodes_y=[nodes_y; y_rot];
            weights=[weights; local_weights];
            
            
        end
    end
end

weights=weights*w_N';
weights=weights(:);

nodes_x=nodes_x(:);
nodes_y=nodes_y(:);
xyw=[nodes_x nodes_y weights];

% method_used=3;
%
% switch method_used
%     case 1
%         % 2007 METHOD: FASTER BUT COMPLICATED.
%         f = fcnchk(intfcn);
%         f_xy = feval(f,nodes_x,nodes_y, varargin{:});
%         cubature_val=(weights'*f_xy)*w_N; % COMPUTING CUBATURE.
%     case 2
%         % MESHGRID LIKE METHOD
%         f = fcnchk(intfcn);
%         f_xy = feval(f,nodes_x,nodes_y, varargin{:});
%         cubature_val=sum(sum(weights.*f_xy));
%     case 3
%         % CLASSICAL VECTOR DESCRIPTION.
%         f = fcnchk(intfcn);
%         nodes_x=nodes_x(:);
%         nodes_y=nodes_y(:);
%         f_xy=feval(f,nodes_x,nodes_y, varargin{:});
%         cubature_val=weights'*f_xy;
% end





%----------------------------------------------------------------------
% FUNCTIONS USED IN THE ALGORITHM.
%----------------------------------------------------------------------


%----------------------------------------------------------------------
% 1. "auto_rotation"
%----------------------------------------------------------------------
function [polygon_bd_rot,rot_matrix,rot_angle,axis_abscissa,vertex_1,vertex_2]=...
    auto_rotation(polygon_bd,vertex_1,vertex_2)


% AUTOMATIC ROTATION OF A CONVEX POLYGON SO THAT "GAUSSIAN POINTS",
% AS IN THE PAPER THEY ARE ALL CONTAINED IN THE CONVEX POLYGON.
% SEE THE PAPER FOR DETAILS.


% FIND DIRECTION AND ROTATION ANGLE.
if length(vertex_1) == 0
    % COMPUTING ALL THE DISTANCES BETWEEN POINTS.A LITTLE TIME CONSUMING
    % AS PROCEDURE.
    distances = points2distances(polygon_bd);
    [max_distances,max_col_comp]=max(distances,[],2);
    [max_distance,max_row_comp]=max(max_distances,[],1);
    vertex_1=polygon_bd(max_col_comp(max_row_comp),:);
    vertex_2=polygon_bd(max_row_comp,:);
    direction_axis=(vertex_2-vertex_1)/max_distance;
else
    direction_axis=(vertex_2-vertex_1)/norm(vertex_2-vertex_1);
end

rot_angle_x=acos(direction_axis(1));
rot_angle_y=acos(direction_axis(2));

if rot_angle_y <= pi/2
    if rot_angle_x <= pi/2
        rot_angle=-rot_angle_y;
    else
        rot_angle=rot_angle_y;
    end
else
    if rot_angle_x <= pi/2
        rot_angle=pi-rot_angle_y;
    else
        rot_angle=rot_angle_y;
    end
end


% CLOCKWISE ROTATION.
rot_matrix=[cos(rot_angle) sin(rot_angle);
    -sin(rot_angle) cos(rot_angle)];

number_sides=size(polygon_bd,1)-1;

polygon_bd_rot=(rot_matrix*polygon_bd')';

axis_abscissa=rot_matrix*vertex_1';



%----------------------------------------------------------------------
% 3. "cubature_rules_1D"
%----------------------------------------------------------------------

function [nodes,weights]=cubature_rules_1D(n,cubature_type)

% SEE WALDVOGEL PAPER. ADDED NODES

% Weights of the Fejer2, Clenshaw-Curtis and Fejer1 quadrature by DFTs
% n>1. Nodes: x_k = cos(k*pi/n)

N=[1:2:n-1]'; l=length(N); m=n-l; K=[0:m-1]';

switch cubature_type
    
    case 1 % FEJER 1.
        v0=[2*exp(i*pi*K/n)./(1-4*K.^2); zeros(l+1,1)];
        v1=v0(1:end-1)+conj(v0(end:-1:2));
        weights=ifft(v1);
        k=(1/2):(n-(1/2)); nodes=(cos(k*pi/n))';
        
    case 2 % FEJER 2.
        v0=[2./N./(N-2); 1/N(end); zeros(m,1)];
        v2=-v0(1:end-1)-v0(end:-1:2);
        wf2=ifft(v2); weights=[wf2;0];
        k=0:n; nodes=(cos(k*pi/n))';
        
    case 3 % CLENSHAW CURTIS.
        g0=-ones(n,1); g0(1+l)=g0(1+l)+n; g0(1+m)=g0(1+m)+n;
        g=g0/(n^2-1+mod(n,2));
        v0=[2./N./(N-2); 1/N(end); zeros(m,1)];
        v2=-v0(1:end-1)-v0(end:-1:2);
        wcc=ifft(v2+g); weights=[wcc;wcc(1,1)];
        k=0:n; nodes=(cos(k*pi/n))';
        
    case 4 % GAUSS LEGENDRE
        beta=0.5./sqrt(1-(2*(1:n)).^(-2));
        T=diag(beta,1)+diag(beta,-1);
        [V,D]=eig(T);
        x=diag(D); [x,index]=sort(x); x=x';
        w=2*V(1,index).^2;
        nodes=x';
        weights=w';
        
end








%----------------------------------------------------------------------
% 3. "points2distances"
%----------------------------------------------------------------------

function distances = points2distances(points)

% Create euclidean distance matrix from point matrix.

% Get dimensions.
[numpoints,dim]=size(points);

% All inner products between points.
distances=points*points';

% Vector of squares of norms of points.
lsq=diag(distances);

% Distance matrix.
distances=sqrt(repmat(lsq,1,numpoints)+repmat(lsq,1,numpoints)'-2*distances);




function [done_table,done_nsize]=convex_decomposition(pts)


%--------------------------------------------------------------------------
% INPUT.
%------------
% pts: POINTS DESCRIBING THE POLYGON (VERTICES ORDERED ANTI-CLOCKWISE).
%      THEY MUST BE A "N x 2" MATRIX. THE LAST VERTEX MUST NOT BE EQUAL TO
%      THE FIRST ONE.
%
%--------------------------------------------------------------------------
% OUTPUT.
%------------
% done_table: SEQUENCE OF POINTS OF THE DECOMPOSITION.
% done_nsize: SEQUENCE OF SUBPOLYGON DIMENSIONS.
%
%--------------------------------------------------------------------------
% ROUTINES.
%------------
% 1. is_convex (attached)
% 2. subdivide_polygon (attached)
%
%--------------------------------------------------------------------------
% EXAMPLE.
%----------
% >> pts=[1 1; 5 1; 6 3; 8 1; 5 7; 2 6; 1 1 ];
% >> [done_table,done_nsize]=convex_decomposition(pts)
% done_table =
%    6.0000    3.0000
%    8.0000    1.0000
%    6.5000    4.0000
%    5.0000    1.0000
%    6.5000    4.0000
%    5.0000    7.0000
%    2.0000    6.0000
%    1.0000    1.0000
% done_nsize =
%     3
%     5
%
% THE FIRST 3 COORDINATES OF done_table ARE THE VERTICES OF THE FIRST
% POLYGON, WHILE THE LAST 5 ARE THE VERTICES OF THE SECOND AND LAST
% POLYGON.
%
%--------------------------------------------------------------------------
% IMPORTANT.
%------------
% FIRST VERTEX MUST NOT BE REPEATED AS LAST VERTEX.
%--------------------------------------------------------------------------
%%--------------------------------------------------------------------------
%% Copyright (C) 2011 Alvise Sommariva, Marco Vianello.
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
%% Author:  Mariano Gentile
%%          Alvise Sommariva <alvise@euler.math.unipd.it>
%%          Marco Vianello   <marcov@euler.math.unipd.it>
%% Date: Apr. 02, 2011 - Apr. 11, 2011.
%--------------------------------------------------------------------------

if norm(pts(1,:)-pts(end,:)) == 0
    pts=pts(1:end-1,:);
end

% TABLE OF POLYGONS ALREADY FOUND AND NUMBER OF SIDES OF EACH POLYGON.
done_table=[];
done_nsize=[];

% TABLE OF POLYGONS TO BE FOUND AND NUMBER OF SIDES OF EACH POLYGON.
todo_table=pts; % vector.
todo_nsize=size(todo_table,1); % vector.

% tisc=[]; tsp=[];

% DECOMPOSING THE POLYGON INTO CONVEX SETS.
while length(todo_nsize) > 0
    
    % FIRST POLYGON OF THE LIST TO BE DECOMPOSED.
    subpolygon=todo_table(1:todo_nsize(1),:);
    
    
    % t(1)=cputime;
    if size(subpolygon,1) > 3
        % fprintf('\n \t %% CONVEX');
        [convex_flag,nonconvex_vertex]=is_convex(subpolygon);
    else
        % fprintf('\n \t %% NON CONVEX');
        convex_flag=1; nonconvex_vertex=[];
    end
    % t(2)=cputime; tisc=[tisc diff(t)];
    
    
    if convex_flag == 1 | length(nonconvex_vertex) == 0 % CONVEX POLYGON.
        
        
        % fprintf('\n \t %% CONVEX ACTION');
        done_table=[done_table; subpolygon];
        done_nsize=[done_nsize; size(subpolygon,1)];
        todo_table=todo_table(todo_nsize(1)+1:end,:);
        todo_nsize=[todo_nsize(2:end)];
        
    else % NON CONVEX POLYGON.
        
        % fprintf('\n \t %% NON CONVEX ACTION');
        %t(1)=cputime;
        [polyg1,polyg2]=subdivide_polygon(subpolygon,...
            nonconvex_vertex);
        %t(2)=cputime; tsp=[tsp diff(t)];
        
        L1=size(polyg1,1); L2=size(polyg2,1);
        todo_table_purged=todo_table(todo_nsize(1)+1:end,:);
        todo_table=[todo_table_purged; polyg1; polyg2];
        todo_nsize=[todo_nsize(2:end); L1; L2];
        
    end
    
end

%fprintf('\n \t is_convex: %2.2e',sum(tisc));
%fprintf('\n \t subdivide_polygon: %2.2e',sum(tsp));











%--------------------------------------------------------------------------
% ATTACHED ROUTINES.
%
% 1. is_convex
% 1a. isleft
%
% 2. subdivide_polygon
% 2a. compute_intersection
%--------------------------------------------------------------------------

function [convex_flag,nonconvex_vertex]=is_convex(subpolygon)

L=size(subpolygon,1);

subpolygon_extended=[subpolygon; subpolygon(1,:); subpolygon(2,:)];

convex_flag=1; nonconvex_vertex=[];

ii=1;
while convex_flag > 0 & ii <= L
    
    P1=subpolygon_extended(ii,:);
    P2=subpolygon_extended(ii+1,:);
    P3=subpolygon_extended(ii+2,:);
    
    convex_flag=isleft(P1,P2,P3);
    ii=ii+1;
    
end

if convex_flag == 0
    nonconvex_vertex=rem(ii,L);
end







function res=isleft(P1,P2,P3)

%--------------------------------------------------------------------------
% PURPOSE.
%----------
% P1, P2, P3 ARE DIFFERENT POINTS OF THE PLANE, STORED AS ROW VECTORS.
% THIS ROUTINE SAYS IF P3 IS ON THE LEFT SIDE W.R.T. THE STRAIGHT
% LINE "r" GOING FROM P1 TO P2 (WITH OBVIOUS ORIENTATION).
% IF res(K)=1 THEN P3(k,:) IS ON THE LEFT OF THE STRAIGHTLINE "r".
% IF res(K)=0 THEN P3(k,:) IS ON THE RIGHT  OF THE STRAIGHTLINE "r".
% IF res(K)=0.5 IS OVER THE LINE "r".
%--------------------------------------------------------------------------

% STAY ON THE LEFT IS EQUIVALENT THAT THE CROSS PRODUCT VECTOR STAYS
% "OVER" THE xy PLANE. WE MODIFY THE VECTORS INVOLVED SINCE THE BUILT-IN
% MATLAB FUNCTION REQUIRES A VECTOR AT LEAST OF DIMENSION 3.

v1=P2-P1; v2=P3-P2;

res=(sign(v1(1)*v2(2)-v1(2)*v2(1))+1)/2;





function  [poly1,poly2,P]=subdivide_polygon(pts,index)

% SUBDIVIDE POLYGON WITH VERTICES pts IN TWO NON OVERLAPPING POLYGONS.
% THE index VERTEX IS A "NOTCH". THE TWO SUBPOLYGONS HAVE A NUMBER OF
% ANGLES INFERIOR THAN pi BIGGER THEN THE ORIGINAL POLYGON.

% REWRITING INITIAL POLYGON WITH FIRST CONCAVE VERTEX AS SECOND INDEX.
L=size(pts,1);
ii=rem((1:L)'-index+2,L);
ss=sign(ii); tt=abs(sign(ss+0.5)-1)/2;
ii=ii+L*(1-abs(ss))+L*tt; [ii1 ii]=sort(ii);
pts_r=pts(ii,:);

v3=(3:L-1)'; v4=(4:L)';

% COMPUTING INTERSECTIONS.
pts_int=[]; t_int=[]; s_int=[];
P1=pts_r(1,:); P2=pts_r(2,:);

for jj=1:size(v3,1)
    P3=pts_r(v3(jj),:);  P4=pts_r(v4(jj),:);
    [pts_loc,t_loc,s_loc]=compute_intersection(P1,P2,P3,P4);
    pts_int=[pts_int; pts_loc];
    t_int=[t_int; t_loc];
    s_int=[s_int; s_loc];
end

[v3 v4 t_int s_int];

% COMPUTING INTERSECTIONS INSIDE THE SEGMENTS.
i_int_01=find(t_int >= 0 & t_int <= 1);
t_int_01=t_int(i_int_01);

pts_int_01=pts_int(i_int_01,:);


% COMPUTING INTERSECTION CLOSER TO THE CONCAVE VERTEX.
s_int_01=s_int(i_int_01);
[smin,imin]=min(abs(s_int_01));
is=i_int_01(imin); ibf=v3(is); iaf=v4(is);
tmin=t_int_01(imin);
P=pts_int(is,:);

% COMPUTING SUBDIVISION.
switch tmin
    case 1 % INTERSECTION IS A VERTEX.
        poly1=pts_r(2:iaf,:);
        vv=[(iaf:L)'; 1];
        poly2=pts_r(vv,:);
    otherwise % STEINER POINT.
        poly1=[pts_r(2:ibf,:); P];
        poly2=[pts_r(1,:); P; pts_r(iaf:L,:)];
        Pa1=poly2(1,:); Pa2=poly2(2,:); Pa3=poly2(3,:);
        delta1=Pa2-Pa1; delta2=Pa3-Pa1;
        % CUTTING POSSIBLE ALIGNMENTS.
        if delta1(1)*delta2(2) == delta1(2)*delta2(1)
            fprintf('\n \t Alignment correction.')
            poly2=[poly2(1,:); poly2(3:end,:)];
        end
end





function [P,t,s]=compute_intersection(P1,P2,P3,P4)

left_matrix=(P2-P1);
right_matrix=(P4-P3);

% PARALLEL LINES CASE.
if left_matrix(1)*right_matrix(2) == ...
        left_matrix(2)*right_matrix(1)
    
    P=[realmax realmax]; t=realmax; s=realmax;
    
else
    
    AA=[left_matrix' right_matrix'];
    bb=(P1-P3)';
    st=AA\bb;
    s=-st(1);
    t=st(2);
    P=P3+t*right_matrix;
    
end







% function flag=check_axis(vertices,Q1,Q2)
%
% % LAST VERTEX AND FIRST ONE ARE EQUAL!
%
% L=size(vertices,1)-1;
% u=Q2-Q1;
% flag=1;
% k=1;
% islefts=[];
%
% while k <= L
%     P1=vertices(k,:);
%     P2=vertices(k+1,:);
%     v=P2-P1;
%
%     if (P1 == Q1) | (P1 == Q2)
%         left_loc=res=isleft(P1,P2,P3);
%     else
%         if (P1 == Q1) | (P1 == Q2)
%             left_loc=res=isleft(P1,P2,P3);
%         else
%         end
%     end
%
%
% end









function [xw]=stroud_conical_rules(ade,vertices)

% INPUT:
% ade: ALGEBRAIC DEGREE OF EXACTNESS.
% vertices: 3 x 2 MATRIX OF VERTICES OF THE SIMPLEX.

% OUTPUT:
% xw: NODES AND WEIGHTS OF STROUD CONICAL RULE TYPE OF ADE ade ON THE SIMPLEX
%     WITH VERTICES vertices.


if nargin < 2 % SEE LYNESS, COOLS,
    % "A survey on numerical cubature over triangles", p.4.
    vertices=[0 0; 1 0; 1 1];
end

[xw]=stroud_conical_rules_ref(ade+1);
bar_coord=[1-xw(:,1) xw(:,1)-xw(:,2) xw(:,2)];
xy=bar_coord*vertices;

A = polyarea(vertices(:,1),vertices(:,2));

w=xw(:,3)*A*2;
xw=[xy w];



function [xw]=stroud_conical_rules_ref(ade)

% SEE LYNESS, COOLS, "A survey on numerical cubature over triangles", p.4.

N=ceil((ade+1)/2);

[t,T]=gauss_jacobi(N,0,0);
t=(t+1)/2; T=T/2;

[x,w]=gauss_jacobi(N,0,1);
x=(x+1)/2; X=w/4;

[wx,wy]=meshgrid(T,X);
ww=wx.*wy;

[yt,xx]=meshgrid(x,x);

[yt,yx]=meshgrid(t,x);
yy=yt.*yx;

xw=[xx(:) yy(:) ww(:)];






function [x,w]=gauss_jacobi(N,a,b,gl)

% GAUSS-JACOBI (LOBATTO) RULE ON (-1,1).
% N IS ...
% a,b ARE THE GAUSS-JACOBI EXPONENTS.
% gl: 0: GAUSS POINTS. 1: GAUSS-LOBATTO POINTS.
% x, w ARE COLUMN VECTORS OF NODES AND WEIGHTS.
%      THE LENGTH OF x AND w IS "N" IF gl=0, "N+2" IF "gl=1".

if nargin < 2
    a=0; b=0;
end

if nargin < 4
    gl = 0;
end

if gl == 0
    ab=r_jacobi(N,a,b);
    xw=gauss(N,ab);
else
    xw=lobatto_jacobi(N,a,b);
end

x=xw(:,1);
w=xw(:,2);




%--------------------------------------------------------------------------
% ADDITIONAL FUNCTIONS BY D.LAURIE AND W.GAUTSCHI.
%--------------------------------------------------------------------------

function ab=r_jacobi(N,a,b)

nu=(b-a)/(a+b+2);
mu=2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2);
if N==1
    ab=[nu mu]; return
end

N=N-1;
n=1:N;
nab=2*n+a+b;
nuadd=(b^2-a^2)*ones(1,N)./(nab.*(nab+2));
A=[nu nuadd];
n=2:N;
nab=nab(n);
B1=4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
B=4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
abadd=[mu; B1; B'];
ab=[A' abadd];


function xw=gauss(N,ab)
N0=size(ab,1); if N0<N, error('input array ab too short'), end
J=zeros(N);
for n=1:N, J(n,n)=ab(n,1); end
for n=2:N
    J(n,n-1)=sqrt(ab(n,2));
    J(n-1,n)=J(n,n-1);
end
[V,D]=eig(J);
[D,I]=sort(diag(D));
V=V(:,I);
xw=[D ab(1,2)*V(1,:)'.^2];


function xw=lobatto_jacobi(N,a,b)

if nargin<2, a=0; end;
if nargin<3, b=a; end
ab=r_jacobi(N+2,a,b);
ab(N+2,1)=(a-b)/(2*N+a+b+2);
ab(N+2,2)=4*(N+a+1)*(N+b+1)*(N+a+b+1)/((2*N+a+b+1)*(2*N+a+b+2)^2);
xw=gauss(N+2,ab);











