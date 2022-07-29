
function demo_optics_V(example)

%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% Cubature on a "polygonal" lens standard domain (PLS) with 5 disks D1, D2,
% D3, D4, D5.
% By "standard" we mean that an external boundary "Bext" is given by the
% intersection of the possibly overlapping "polygonal disks" (PD) D1, D4,
% D5, while the internal boundary "Bint" is determined by the boundary of
% the union of D2 and D3 (the so called "pupils").
%
% We suppose that
% 1. Bint is in the open domain spanned by Bext, i.e. of the domain
% intersection(D1,D4,D5).
% 2. Bint is connected, i.e. D2 intersects D3.
% 3. The "PDs" D1, D2, D3, D4, D5 have center on the straight line "x=0",
%    i.e. on the y-axes. Thus, to determine them, we need to assign the
%    coordinate of the center of the k-th disk, say "y(k)" and its radius
%    say r(k).
%
% Given the vectors "y", "r" of dimension 1 x 5, it computes first
% an algebraic cubature formula of degree "deg" on the standard domain
% defined by the P.D.s "Dk", with center (0,y(k)) and radius "r(k)".
%
% The cubature rule has nodes xyw(:,1:2) and weights xyw(:,3)>=0. Next it
% determines a compressed rule with (deg+1)(deg+2)/2 nodes based in xyw
% rule. Its nodes/weights are stored in the vector "xywc".
%
% We compare these results with those on the relative "circular lens
% standard domain". Such rules are exact, and depending on the case may
% have/have not internal nodes, positive weights.
%
% In each instance, if "I" is the computed integral, we define the
% parameter
%
%              sigma=sqrt(I/area_poly-(I/area_poly)^2)
%
% for polygonal lenses ("area_poly" is the area of the polygonal lens),
% while for circular lenses,
%
%              sigma=sqrt(I/area_circ-(I/area_circ)^2)
%
% for polygonal lenses ("area_circ" is the area of the circular lens).
%
% .........................................................................
% VERY IMPORTANT:
% .........................................................................
% 1. We compute sigma coeffs, mimicking Fourier expansions.
%    The coefficients are in decreasing order.
% 2. We use orthonormal Zernike polynomials:
%
%       n    m    Zernike function           Normalization
%       --------------------------------------------------
%       0    0    1                                 1
%       1    1    r * cos(theta)                    2
%       1   -1    r * sin(theta)                    2
%       2   -2    r^2 * cos(2*theta)             sqrt(6)
%       2    0    (2*r^2 - 1)                    sqrt(3)
%       2    2    r^2 * sin(2*theta)             sqrt(6)
%       3   -3    r^3 * cos(3*theta)             sqrt(8)
%       3   -1    (3*r^3 - 2*r) * cos(theta)     sqrt(8)
%       3    1    (3*r^3 - 2*r) * sin(theta)     sqrt(8)
%       3    3    r^3 * sin(3*theta)             sqrt(8)
%       4   -4    r^4 * cos(4*theta)             sqrt(10)
%       4   -2    (4*r^4 - 3*r^2) * sin(2*theta) sqrt(10)
%       4    0    6*r^4 - 6*r^2 + 1              sqrt(5)
%       4    2    (4*r^4 - 3*r^2) * cos(2*theta) sqrt(10)
%       4    4    r^4 * sin(4*theta)             sqrt(10)
%       --------------------------------------------------
%
% Single index numeration of these polynomials may be different is some
% contexts.
%
% 3. The polynomials we test are such that the coefficients relative to the
% first three Zernike polynomials are null.
%
%--------------------------------------------------------------------------
% INPUT.
%--------------------------------------------------------------------------
% example: DOMAIN OF INTEREST (SEE FUNCTION lens_data AT LINE 287).
%--------------------------------------------------------------------------
% ROUTINES.
%--------------------------------------------------------------------------
% * make_experiment (ATTACHED HERE)
% * lens_data      (ATTACHED HERE)
% * polygauss_2018 (EXTERNAL)
% * compresscub    (EXTERNAL)
% * gqlune         (EXTERNAL)
% * trigauss       (EXTERNAL)
% * rjacobi        (EXTERNAL)
% * gauss          (EXTERNAL)
% * vap            (EXTERNAL)
% * pupilscub_PI (EXTERNAL)
% * compute_integral (ATTACHED HERE)
% * do_plots       (ATTACHED HERE)
% * define_domain  (ATTACHED HERE)
%--------------------------------------------------------------------------
% IMPORTANT:
% * THE CODE REQUIRES THE TOOLBOX CONTAINING POLYSHAPE (AVAILABLE FROM
% MATLAB R2017b VERSIONS).
% * USED IN THE PAPER
%   B. BAUMAN , A. SOMMARIVA y, AND M. VIANELLO,
%   "COMPRESSED CUBATURE OVER TRIANGULATED DOMAINS",
%   TO DETERMINE TABLE 6.1.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Copyright (C) 2019- Brian J. Bauman, Alvise Sommariva, Marco Vianello.
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
%% Date: December 4, 2018
%% Date: May 6, 2019
%--------------------------------------------------------------------------


% ............................ settings ...................................

% Parameters explanations:
% * example: performs the relative example, relatively to what is set in
%     "define_domain" routine.
% * ade: degree of precision.
% * Ntest: number of tests that are performed for cputime purposes.
% * Nsides_factor: sides on each polygonal disk.
% * do_plots: for plotting domain and sets, use do_plots=1.

example=0;
ade=8;
Ntest=100;
Nsides_factorV=[100 1600];
do_plots=0; 

% ........................... main call ...................................
diary on
for Nsides_factor=Nsides_factorV
    make_experiment(example,ade,Ntest,Nsides_factor);
end
diary off;





%==========================================================================
% ATTACHED ROUTINES
%==========================================================================

%==========================================================================
% demo_optics_III_slave
%==========================================================================

function make_experiment(example,ade,Ntest,Nsides_factor,do_plots)

% INPUT:
% example: performs the relative example, relatively to what is set in
%     "define_domain" routine.
% ade: degree of precision.
% Ntest: number of tests that are performed for cputime purposes.
% Nsides_factor: sides on each polygonal disk.
% do_plots: for plotting domain and sets, use do_plots=1.
%
% IMPORTANT: degree of precision will be 2*ade since for computing "sigma"
% parameter we must compute the integral of a polynomial and of its square.

% .......................... troubleshooting ..............................

if nargin < 1
    example=0;
end

if nargin < 2
    ade=10;
end

if nargin < 3
    Ntest=100; % cputests are done, averaging on Ntests experiments.
end

if nargin < 4
    Nsides_factor=8000;
end

if nargin < 5
    do_plots=0;
end

% .......................... settings .....................................

pos=1; % 1: lsqnonneg, 2: NNLSlab routine.
draw_disks=1; % plot of all the disks determining the domain.
interstats=0; % write statistics of any test: 0: no, 1: yes.

%==========================================================================
%                       MAIN PROGRAM STARTS HERE.
%==========================================================================

fprintf('\n \t .................... EXAMPLE: %3.0f .......... \n',...
    example)

%==========================================================================
% GENERAL ALGEBRAIC QUADRATURE.
% FORMULA WITH NON NEGATIVE WEIGHTS.
%==========================================================================

% .......................... define domain ................................

tic;
[P,y,r,P1v,P2v,P3v,P4v,P5v]=define_domain(example,Nsides_factor);
tP=toc;

% .......................... cubature rule ................................

tic;
[xyw,xvc,yvc,pgon,tri]=polygauss_2018(ade,P);
t(1)=toc;

% ...................... compressed cubature rule .........................

tic;
[xyc,w,momerr1] = comprexcub(ade,xyw(:,1:2),xyw(:,3),pos);
t(2)=toc;
xywc=[xyc w];

% ................ reference rule circular domain (no PI type) ............

[xywvap,ptsvap,wvap,momerr]=vap(ade,y,r);

% ................ reference rule circular domain (PI type) ...............

[xyw_nsl,xywc_nsl]=pupilscub_PI(y,r,ade);

% ...................... stats initialisations ............................

aeV=[];
reV=[];
sigma_aeV=[];
sigmac_aeV=[];
sigma_reV=[];
sigmac_reV=[];



for kk=1:Ntest
    
    % .................... defining integrand coeffs ......................
    ade_poly=floor(ade/2);
    dim=(ade_poly+1)*(ade_poly+2)/2;
    
    % define coefficients of ORTHONORMAL Zernike polynomials:
    % we mimick coefficients of Fourier coeffs of an L2 function.
    coeffs=sort(rand(dim,1),'descend'); 
    coeffs(1:3)=zeros(3,1);
    
    
    
    % ............... computing integrals (several methods) ...............
    
    % Note that we compute the integral of a polynomial and of its square.
    
    % ... current method ...
    [Ifull,Ifull2]=compute_integrals(ade_poly,coeffs,xyw(:,1:2),xyw(:,3));
    % ... compressed method ...
    [Ic,Ic2]=compute_integrals(ade_poly,coeffs,xyc,w);
    % ... non PI method ...
    [Idisk,Idisk2]=compute_integrals(ade_poly,coeffs,xywvap(:,1:2),...
        xywvap(:,3));
    % ... non PI compressed method ...
    [Idiskc,Idiskc2]=compute_integrals(ade_poly,coeffs,ptsvap,wvap);
    % ... PI method ...
    [Idisk_nsl,Idisk_nsl2]=compute_integrals(ade_poly,coeffs,...
        xyw_nsl(:,1:2),xyw_nsl(:,3));
    % ... PI compressed method ...
    [Idiskc_nsl,Idiskc_nsl2]=compute_integrals(ade_poly,coeffs,...
        xywc_nsl(:,1:2),xywc_nsl(:,3));
    
    % ... polygonal/circular domain areas ...
    area_poly=sum(xyw(:,3));
    area_circ=sum(xywc_nsl(:,3));
    
    
    % ............... comparing integrals (several methods) ...............
    
    aef=abs(Ifull-Idisk);
    aec=abs(Ic-Idisk);
    ref=abs(Ifull-Idisk)/abs(Idisk);
    rec=abs(Ic-Idisk)/abs(Idisk);
    aeV=[aeV; aef];
    reV=[reV; ref];
    
    
    % ..................... computing sigmas ..............................
    
    sigma_full(kk)=sqrt(Ifull2/area_poly-(Ifull/area_poly)^2);
    sigma_c(kk)=sqrt(Ic2/area_poly-(Ic/area_poly)^2);
    sigma_disk(kk)=sqrt(Idisk2/area_circ-(Idisk/area_circ)^2);
    sigma_diskc(kk)=sqrt(Idiskc2/area_circ-(Idiskc/area_circ)^2);
    sigma_disk_nls(kk)=sqrt(Idisk_nsl2/area_circ-(Idisk_nsl/area_circ)^2);
    sigma_diskc_nls(kk)=sqrt(Idiskc_nsl2/area_circ-(Idiskc_nsl/area_circ)^2);
    
    % ..................... comparing sigmas ..............................
    
    sigma_aeV=[sigma_aeV; abs(sigma_full(kk)-sigma_disk_nls(kk))];
    sigmac_aeV=[sigmac_aeV; abs(sigma_c(kk)-sigma_disk_nls(kk))];
    sigma_reV=[sigma_aeV; abs(sigma_full(kk)-sigma_disk_nls(kk))...
        /abs(sigma_disk_nls(kk))];
    sigmac_reV=[sigma_aeV; abs(sigma_c(kk)-sigma_disk_nls(kk))...
        /abs(sigma_disk_nls(kk))];
    
    
    % .................... current example statistics .....................
    
    if interstats == 1
        
        fprintf('\n \n \t ----------------------------------------------');
        
        fprintf('\n \n \t ** Integrals: numerical values');
        fprintf('\n \t Ifull polyg.: %1.15e',Ifull);
        fprintf('\n \t Icomp polyg.: %1.15e',Ic);
        fprintf('\n \t Ifull vap   : %1.15e',Idisk);
        fprintf('\n \t Icomp vap   : %1.15e',Idiskc);
        fprintf('\n \t Ifull c.nsl.: %1.15e',Idisk_nsl);
        fprintf('\n \t Icomp c.nsl.: %1.15e',Idiskc_nsl);
        
        fprintf('\n \n \t ** Integrals: absolute errors');
        
        
        fprintf('\n \t Ifull polyg.: %1.15e',aef);
        fprintf('\n \t Icomp polyg.: %1.15e',aec);
        
        fprintf('\n \n \t ** Integrals: relative errors');
        
        fprintf('\n \t Ifull polyg.: %1.15e',ref);
        fprintf('\n \t Icomp polyg.: %1.15e \n \n',rec);
        
        fprintf('\n \n \t ** Integrals: RMS');
        
        fprintf('\n \t Ifull polyg. RMS: %1.15e',sigma_full(kk));
        fprintf('\n \t Icomp polyg. RMS: %1.15e',sigma_c(kk));
        fprintf('\n \t Ivap. circ.  RMS: %1.15e',sigma_disk(kk));
        fprintf('\n \t Icomp. vap   RMS: %1.15e',sigma_diskc(kk));
        fprintf('\n \t Inls  circ.  RMS: %1.15e',sigma_disk_nls(kk));
        fprintf('\n \t Icomp  nls   RMS: %1.15e',sigma_diskc_nls(kk));
        
    end
    
end


% .......................... final statistics .............................

fprintf('\n \n \t ********** Final summary **********');
% stats
fprintf('\n \n \t ** Main data');
fprintf('\n \t ade         : %3.0f',ade);
fprintf('\n \t sides x poly: %5.0f',Nsides_factor);
[xP,yP] = boundary(P);
distinct_polys=sum(isnan(xP)); % some nans say if there are disconn. polygs.
tot_poly_sides=length(xP)-distinct_polys;
fprintf('\n \t lens sides  : %7.0f',tot_poly_sides);

fprintf('\n \n \t ** Cardinalities');
fprintf('\n \t polyg. lens full: %7.0f',size(xyw(:,1),1));
fprintf('\n \t polyg. lens comp: %7.0f',size(xywc(:,1),1));
fprintf('\n \t vap lens full   : %7.0f',size(xywvap(:,1),1));
fprintf('\n \t vap lens comp   : %7.0f',size(ptsvap(:,1),1));
fprintf('\n \t nsl lens full: %7.0f',size(xyw_nsl(:,1),1));
fprintf('\n \t nsl lens comp: %7.0f',size(xywc_nsl(:,1),1));

fprintf('\n \n \t ** CPUtimes');
fprintf('\n \t TRIANGULATION : %1.15e',tP);
fprintf('\n \t FULL RULE     : %1.15e',t(1));
fprintf('\n \t COMP RULE     : %1.15e',t(2));

fprintf('\n \n \t ** Cubature Errors');
fprintf('\n \t BETTER ABS.ERR. : %1.15e',min(aeV));
fprintf('\n \t WORST ABS.ERR.  : %1.15e',max(aeV));
fprintf('\n \t BETTER REL.ERR. : %1.15e',min(reV));
fprintf('\n \t WORST REL.ERR.  : %1.15e',max(reV));

fprintf('\n \n \t ** Average RMS');
fprintf('\n \t Ifull polyg. RMS: %1.15e',sum(sigma_full)/Ntest);
fprintf('\n \t Icomp polyg. RMS: %1.15e',sum(sigma_c)/Ntest);
fprintf('\n \t Ivap. circ.  RMS: %1.15e',sum(sigma_disk)/Ntest);
fprintf('\n \t Icomp. vap   RMS: %1.15e',sum(sigma_diskc)/Ntest);
fprintf('\n \t Inls  circ.  RMS: %1.15e',sum(sigma_disk_nls)/Ntest);
fprintf('\n \t Icomp  nls   RMS: %1.15e',sum(sigma_diskc_nls)/Ntest);

fprintf('\n \n \t ** Sigma Errors (poly. vs circ., full version)');
fprintf('\n \t BETTER ABS.ERR. : %1.15e',min(sigma_aeV));
fprintf('\n \t WORST ABS.ERR.  : %1.15e',max(sigma_aeV));
fprintf('\n \t AVER. ABS.ERR.  : %1.15e',mean(sigma_aeV));
fprintf('\n \t BETTER REL.ERR. : %1.15e',min(sigma_reV));
fprintf('\n \t WORST REL.ERR.  : %1.15e',max(sigma_reV));
fprintf('\n \t AVER. REL.ERR.  : %1.15e',mean(sigma_reV));

fprintf('\n \n \t ** Sigma Errors (poly. vs circ., compressed version)');
fprintf('\n \t BETTER ABS.ERR. : %1.15e',min(sigmac_aeV));
fprintf('\n \t WORST ABS.ERR.  : %1.15e',max(sigmac_aeV));
fprintf('\n \t AVER. ABS.ERR.  : %1.15e',mean(sigmac_aeV));
fprintf('\n \t BETTER REL.ERR. : %1.15e',min(sigmac_reV));
fprintf('\n \t WORST REL.ERR.  : %1.15e',max(sigmac_reV));
fprintf('\n \t AVER. REL.ERR.  : %1.15e',mean(sigmac_reV));

fprintf('\n \n');


% ....................... plotting nodes/domain ...........................

if do_plots == 1
    clf(figure(1));
    figure(1)
    do_plots(P,xyw,xywc,draw_disks,P1v,P2v,P3v,P4v,P5v);
    
    clf(figure(2));
    figure(2)
    do_plots(P,xyw_nsl,xywc_nsl,draw_disks,P1v,P2v,P3v,P4v,P5v);
end





%--------------------------------------------------------------------------
% compute_integrals
%--------------------------------------------------------------------------

function [I,I2]=compute_integrals(deg,coeffs,xy,w)

V=vandermonde_zernike(xy,deg);
fxy=V*coeffs;
I=w'*fxy;
I2=w'*(fxy.^2);




%--------------------------------------------------------------------------
% define_domain
%--------------------------------------------------------------------------

function [P,y,r,P1v,P2v,P3v,P4v,P5v]=define_domain(example,Nsides_factor)

[y,r]=lens_data(example);

switch example
    case 0
        % NUMBER OF SIDES OF POLYGONS D1,D2,D3,D4,D5.
        NsidesV=Nsides_factor*ones(5,1);
        th=linspace(0,2*pi,NsidesV(1)); th=(th(1:end-1))';
        C1=[0 y(1)]; P1v=C1+r(1)*[cos(th) sin(th)]; P1=polyshape(P1v);
        
        th=linspace(0,2*pi,NsidesV(2)); th=(th(1:end-1))';
        C2=[0 y(2)]; P2v=C2+r(2)*[cos(th) sin(th)]; P2=polyshape(P2v);
        
        th=linspace(0,2*pi,NsidesV(3)); th=(th(1:end-1))';
        C3=[0 y(3)]; P3v=C3+r(3)*[cos(th) sin(th)]; P3=polyshape(P3v);
        
        th=linspace(0,2*pi,NsidesV(4)); th=(th(1:end-1))';
        C4=[0 y(4)]; P4v=C4+r(4)*[cos(th) sin(th)]; P4=polyshape(P4v);
        
        th=linspace(0,2*pi,NsidesV(5)); th=(th(1:end-1))';
        C5=[0 y(5)]; P5v=C5+r(5)*[cos(th) sin(th)]; P5=polyshape(P5v);
        
        Pout=intersect(P1,P4);
        Pout=intersect(Pout,P5);
        Pin=union(P2,P3);
        P=subtract(Pout,Pin);
    case 1
        % NUMBER OF SIDES OF POLYGONS D1,D2,D3,D4,D5.
        NsidesV=Nsides_factor*ones(5,1);
        th=linspace(0,2*pi,NsidesV(1)); th=(th(1:end-1))';
        C1=[0 y(1)]; P1v=C1+r(1)*[cos(th) sin(th)]; P1=polyshape(P1v);
        
        th=linspace(0,2*pi,NsidesV(2)); th=(th(1:end-1))';
        C2=[0 y(2)]; P2v=C2+r(2)*[cos(th) sin(th)]; P2=polyshape(P2v);
        
        th=linspace(0,2*pi,NsidesV(3)); th=(th(1:end-1))';
        C3=[0 y(3)]; P3v=C3+r(3)*[cos(th) sin(th)]; P3=polyshape(P3v);
        
        th=linspace(0,2*pi,NsidesV(4)); th=(th(1:end-1))';
        C4=[0 y(4)]; P4v=C4+r(4)*[cos(th) sin(th)]; P4=polyshape(P4v);
        
        th=linspace(0,2*pi,NsidesV(5)); th=(th(1:end-1))';
        C5=[0 y(5)]; P5v=C5+r(5)*[cos(th) sin(th)]; P5=polyshape(P5v);
        
        Pout=intersect(P1,P4);
        Pout=intersect(Pout,P5);
        Pin=union(P2,P3);
        P=subtract(Pout,Pin);
end




%--------------------------------------------------------------------------
% lens_data
%--------------------------------------------------------------------------

function [y,r]=lens_data(example)


% For each disk Dk we record the coordinates Dk=[y(k) r(k)] so that
% its center is (0,y(k)) and the radius is r(k).

switch example
    case 0
        alpha=1;
        theta=alpha*1.75*pi/180;
        y3=-3.876*tan(theta);
        y4=y3;
        y5=-12.310*tan(theta);
        y=[0 0 y3 y4 y5];
        r=[1 0.6120 0.5663 1.0761 1.2810];
    case 1
        y=[0 0 -0.5 -0.7 -1.8];r=[1 0.5 0.4 1.6 2.5];
    otherwise
        fprintf('\n \t WARNING: NO DEMO. \n \n'); return;
end





%--------------------------------------------------------------------------
% do_plots
%--------------------------------------------------------------------------

function do_plots(P,xyw,xywc,draw_disks,P1v,P2v,P3v,P4v,P5v)

% PLOT.
plot(P,'FaceColor',[0.5 0.5 0.5]);
hold on;
plot(xyw(:,1),xyw(:,2),'+','LineWidth',0.5,...
    'MarkerEdgeColor',[0.25 0.25 0.25],...
    'MarkerFaceColor',[0.25 0.25 0.25],...
    'MarkerSize',0.5)
plot(xywc(:,1),xywc(:,2),'go','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',6);
%title('Algebraic cubature over polygonal and vignetted pupils');
%legend('Domain','Full set','Compressed set');
v=1:size(P1v,1);
v=[v 1];

plot(P1v(v,1),P1v(v,2),'LineWidth',0.5,'Color',[0.5 0.5 0.5]);
plot(P2v(v,1),P2v(v,2),'LineWidth',0.5,'Color',[0.5 0.5 0.5]);
plot(P3v(v,1),P3v(v,2),'LineWidth',0.5,'Color',[0.5 0.5 0.5]);
plot(P4v(v,1),P4v(v,2),'LineWidth',0.5,'Color',[0.5 0.5 0.5]);
plot(P5v(v,1),P5v(v,2),'LineWidth',0.5,'Color',[0.5 0.5 0.5]);
axis off
axis equal
hold off;



