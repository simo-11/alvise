
function demo_adaptive_polygons

%--------------------------------------------------------------------------
% Object:
% Demo on algebraic polynomial cubature on bivariate polygonal domains.
% The determination of the X pointset varies with the choice of the
% parameter "pts_type".
%--------------------------------------------------------------------------
% Dates:
% Written on 29/10/2020: M. Vianello;
%
% Modified on:
% 15/07/2021: A. Sommariva.
%--------------------------------------------------------------------------

clear;

%--------------------------------------------------------------------------
% Function to study. The variable "function_example" can be:
%
% case 1, f=@(x,y) 1./(1+25*(x+2*y).^2);
% case 2, f=@(x,y) sin(x+2*y);
% case 3, f=@(x,y) exp(x+2*y).*cos(10*(x+2*y)).*tanh(4*(x+2*y));
% case 4, f=@(x,y) exp(x+2*y);
% case 5, f=@(x,y) exp(-(x+2*y));
% case 6, f=@(x,y) (1.2*(x+2*y).^5+3*(x+2*y)+0.02*(x+2*y).^4-2).^4;
% case 7, f=@(x,y) (x+2*y).^20+1;
% case 8, f=@(x,y) (1.2*(x+2*y).^3+0.02*(x+2*y).^4+3*(x+2*y)-2).^2;
% case 9, f=@(x,y) (1.2*(x+2*y).^2+0.02*(x+2*y).^2+3*(x+2*y)-2).^2
%--------------------------------------------------------------------------

function_example=2;


%--------------------------------------------------------------------------
% Domain to be considered.
%
% case 1, convex domain,
% case 2, non convex and non simply connected domain.
%--------------------------------------------------------------------------

domain_example=2;

%--------------------------------------------------------------------------
% tri_type : in case a triangle must be subdivided, we use
%            0: triangulation via midpoint of each side  (default)
%            1: triangulation via midpoint of longer side
%            2: triangulation via centroid
%--------------------------------------------------------------------------

tri_type=2;

% ........................ Main code below ................................

[g,gstr]=define_function(function_example);
vertices=define_domain(domain_example);

tol=10^(-12);

tic;
[IH,IL,flag,iters,PG,tri_type]=cub_polygon_adaptive(vertices,g,tol,...
    tri_type);
cpus=toc;

% ... statistics ...

fprintf('\n \t f=@(x,y)'); disp(gstr);
fprintf('\n \t IH  : %1.15e',IH);
fprintf('\n \t IL  : %1.15e',IL);
fprintf('\n \t AE  : %1.15e',abs(IL-IH));
if abs(IH) > 0
    fprintf('\n \t RE  : %1.15e',abs(IL-IH)/abs(IH));
end
fprintf('\n');

fprintf('\n \t Iter: %6.0f',iters);
fprintf('\n \t Tol : %1.1e',tol);
fprintf('\n \t Cpu : %1.1e',cpus);
fprintf('\n');

switch tri_type
    case 0
        fprintf('\n \t * Triangulation via side midpoints');
    case 1
        fprintf('\n \t * Triangulation via midpoint longer side');
    case 2
        fprintf('\n \t * Triangulation via centroid');
end
fprintf('\n');

switch flag
    case 0
        fprintf('\n \t * Correct termination');
    case 1
        fprintf('\n \t * Early termination: too many triangles to be analysed');
end
fprintf('\n \n');

% ... plot domain ...
plot(PG);










function [f,fs]=define_function(example)

switch example
    case 1
        f=@(x,y) 1./(1+25*(x+2*y).^2);
        fs='1./(1+25*(x+2*y).^2)';
    case 2
        f=@(x,y) sin(x+2*y);
        fs='sin(x+2*y)';
    case 3
        f=@(x,y) exp(x+2*y).*cos(10*(x+2*y)).*tanh(4*(x+2*y));
        fs='exp(x+2*y).*cos(10*(x+2*y)).*tanh(4*(x+2*y))';
    case 4
        f=@(x,y) exp(x+2*y);
        fs='exp(x+2*y)';
    case 5
        f=@(x,y) exp(-(x+2*y));
        fs='exp(-(x+2*y))';
    case 6
        f=@(x,y) (1.2*(x+2*y).^5+3*(x+2*y)+0.02*(x+2*y).^4-2).^4;
        fs='(1.2*(x+2*y).^5+3*(x+2*y)+0.02*(x+2*y).^4-2).^4';
    case 7
        f=@(x,y) (x+2*y).^20+1;
        fs='(x+2*y).^20+1';
    case 8
        f=@(x,y) (1.2*(x+2*y).^3+0.02*(x+2*y).^4+3*(x+2*y)-2).^2;
        fs='(1.2*(x+2*y).^3+0.02*(x+2*y).^4+3*(x+2*y)-2).^2';
    case 9
        f=@(x,y) (1.2*(x+2*y).^2+0.02*(x+2*y).^2+3*(x+2*y)-2).^2;
        fs='(1.2*(x+2*y).^2+0.02*(x+2*y).^2+3*(x+2*y)-2).^2';
end










function vertices=define_domain(example)


switch example
    
    case 1 % simply connected domain
        vertices=[0.1 0; 0.7 0.2; 1 0.5; 0.75 0.85; 0.5 1; 0 0.25; ...
            0.5 1; 0 0.25; 0.1 0];
        
    case 2 % domain with holes (using polyshape)
        Nsides=10;
        th=linspace(0,2*pi,Nsides); th=(th(1:end-1))';
        % first polygon.
        polygon1=[cos(th) sin(th)]; P1=polyshape(polygon1);
        % second polygon.
        polygon2=2+[cos(th) sin(th)]; P2=polyshape(polygon2);
        % third polygon.
        polygon3=0.5*[cos(th) sin(th)]; P3=polyshape(polygon3);
        
        P=subtract(P1,P3);
        P=union(P,P2);
        
        [XX, YY] = boundary(P);
        vertices=[XX YY];
        
end



