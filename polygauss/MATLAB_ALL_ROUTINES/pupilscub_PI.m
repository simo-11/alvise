
function [xyw,xywc,momerr,flag,D,xyw_nr,xyw_sr,xyw_midT,xyw_midL,xyw_midC]=...
    pupilscub_PI(y,r,deg)

%--------------------------------------------------------------------------
% OBJECT:
% Algebraic cubature of degree "deg" on standard lenses. 
% We consider cubature on a domain formed by 5 disks Dk whose centers are
% (0,y(k)) and radii r(k).
% The disks must have the following properties:
% 1. "D1, D4, D5 must not intersect the union of D2 with D3".
% 2. "The disks D1 and D2 must be always defined".
% 3. "The disks D2 and D3 can be disjoined or not".
% The external component of the boundary of the domain consists of the 
% intersection of D1, D4, D5.
% The internal component of the boundary consists on the boundary of the 
% union of the disks D2 and D3.
%--------------------------------------------------------------------------
% INPUT:
% y: n x 1 vector (with n <= 5) such that the disk Dk has center (0,y(k)).
% r: n x 1 vector (with n <= 5) such that the disk Dk has radius r(k).
% deg: algebraid degree of the quadrature rule.
%--------------------------------------------------------------------------
% OUTPUT:
% xyw: M x 3 vector that stores the quadrature rule with algebraic degree 
%     of precision "deg". The nodes are
%     (x(k),y(k)) where x(k)=xyw(k,1), y(k)=xyw(k,2), while the weights 
%     w(k)=xyw(k,3), i.e. they are stored in the third column of xyw.
%     We stress that the nodes are in the domain and that the weights are
%     nonnegative.
% xywc: N x 3 vector that stores the quadrature rule with algebraic degree 
%      of precision "deg". In particular N <= (deg+1)*(deg+2)/2. The nodes 
%      are(x(k),y(k)) where x(k)=xywc(k,1), y(k)=xywc(k,2), while the  
%      weights  w(k)=xywc(k,3), i.e. they are stored in the third column of 
%      xyw.
%     We stress that the nodes are in the domain and that the weights are
%     nonnegative. The rule is obtained by compression from xyw.
% momerr: it is a quantity that shows how much the compressed rule matches
%     the moments of a certain polynomial basis.
% flag: 1 if the cubature rule is stisfactory, 0 otherwise.
%--------------------------------------------------------------------------
% note: trigauss.m, lawsonhanson.m are the only external routine required 
%      to run this file.
%      The routine has been tested on Matlab 2017A on a Mac Book Pro.
%--------------------------------------------------------------------------
% authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15 2017
%--------------------------------------------------------------------------

flag=1;

xyw=[];

% INITIAL CHECKS.
[P1,P2,Q1,Q2,D1,D2,D3,D4,D5,flag]=troubleshooting(y,r);

% [P1,P2,Q1,Q2,D1,D2,D3,D4,D5]=compute_fundamental_pts(y,r);
[AA1,AA2,BB1,BB2,CC1,CC2,AAA1,AAA2,BBB1,BBB2,CCC1,CCC2]=...
    fundamental_pts_pupils(D2,D3);

%fprintf('\n \t northern_region');
ytop=BB1(2);
xyw_nr=northern_region(D1,D4,D5,P1,P2,Q1,Q2,ytop,deg);

%fprintf('\n \t southern_region');
ybottom=AAA1(2);
xyw_sr=southern_region(D1,D4,D5,P1,P2,Q1,Q2,ybottom,deg);

if norm(D2-D3) > 0 % 2 overlapping pupils.
    
    %fprintf('\n \t midupper_region');
    if size(CC1,1) > 0
        xyw_T=central_region(D2,D1,D4,D5,P1,P2,Q1,Q2,CC1,BB1,deg);
        xyw_L=central_region(D2,D1,D4,D5,P1,P2,Q1,Q2,AA1,CC1,deg);
        xyw_midT=[xyw_T; xyw_L];
    else
        xyw_midT=central_region(D2,D1,D4,D5,P1,P2,Q1,Q2,AA1,BB1,deg);
    end
    
    %fprintf('\n \t midlower_region');
    if size(CCC1,1) > 0
        xyw_T=central_region(D3,D1,D4,D5,P1,P2,Q1,Q2,CCC1,BBB1,deg);
        xyw_L=central_region(D3,D1,D4,D5,P1,P2,Q1,Q2,AAA1,CCC1,deg);
        xyw_midL=[xyw_T; xyw_L];
    else
        xyw_midL=central_region(D3,D1,D4,D5,P1,P2,Q1,Q2,AAA1,BBB1,deg);
    end
    
    if flag == 0 % NON OVERLAPPING PUPILS
        ytop=D2(2)-D2(3); ybottom=D3(2)+D3(3);
        % fprintf('\n \t cushion_region')
        xyw_midC=cushion_region(D1,D4,D5,P1,P2,Q1,Q2,ytop,ybottom,deg);
    else
        xyw_midC=[];
    end
    
    
else % pupil consisting of one disk (simple pupil).
    % when there is one pupil, D2, D3 are previously managed to be equal.
    AA1=[D2(1) D2(2)-D2(3)];
    BB1=[D2(1) D2(2)+D2(3)];
    CC1=[D2(1)+D2(3) D2(2)];
    
    xyw_T=central_region(D2,D1,D4,D5,P1,P2,Q1,Q2,CC1,BB1,deg);
    xyw_L=central_region(D2,D1,D4,D5,P1,P2,Q1,Q2,AA1,CC1,deg);
    xyw_midT=[xyw_T; xyw_L];
    xyw_midL=[];
    xyw_midC=[];
end

xyw=[xyw_nr; xyw_sr; xyw_midT; xyw_midL; xyw_midC];

X=xyw(:,1:2); omega=xyw(:,3); pos=1;
[pts,w,momerr] = comprexcub(deg,X,omega,pos);

xywc=[pts w];

D=[D1; D2; D3; D4; D5];








%--------------------------------------------------------------------------
% compute_fundamental_pts
%--------------------------------------------------------------------------

function [P1,P2,Q1,Q2,D1,D2,D3,D4,D5]=compute_fundamental_pts(y,r)

%--------------------------------------------------------------------------
% input: 
% y: circle center quotas
% r: circle radii
%--------------------------------------------------------------------------
% output:
% P1, P2: intersection of the boundary of D1 with the boundary of D4. In
%         particular P2 is on the left of P1.
% Q1, Q2: intersection of the boundary of D1 with the boundary of D4. In
%         particular Q2 is on the left of Q1.
% D1, D2, D3, D4, D5: disks determining the standard domain. In particular
%         Dk=[0 y(k) r(k)] and D2 is "above" D3.
%--------------------------------------------------------------------------
% by Alvise Sommariva and Marco Vianello
% University of Padova, November 2 2017
%--------------------------------------------------------------------------

D1=[0 y(1) r(1)]; D4=[0 y(4) r(4)]; D5=[0 y(5) r(5)];

if y(4)+r(4) < y(5)+r(5)
    [D4,D5]=deal(D5,D4);
end


% Determining relevant points.
if norm(D4-D1) > 0
    [xout14,yout14] = circcirc(0,y(1),r(1),0,y(4),r(4));
    
    C1=[0 y(1)];
    P14_1=[xout14(1) yout14(1)];
    P14_1s=P14_1-C1;  th14(1)=mycart2pol([P14_1s(1) P14_1s(2)]);
    P14_2=[xout14(2) yout14(2)];
    P14_2s=P14_2-C1;  th14(2)=mycart2pol([P14_2s(1) P14_2s(2)]);
    [min_th14]=min(th14(1),th14(2));
else
    min_th14=pi/2;
end

if norm(D5-D1) > 0
    [xout15,yout15] = circcirc(0,y(1),r(1),0,y(5),r(5));
    C1=[0 y(1)];
    P15_1=[xout15(1) yout15(1)];
    P15_1s=P15_1-C1;  th15(1)=mycart2pol([P15_1s(1) P15_1s(2)]);
    P15_2=[xout15(2) yout15(2)];
    P15_2s=P15_2-C1;  th15(2)=mycart2pol([P15_2s(1) P15_2s(2)]);
    [min_th15]=min(th15(1),th15(2));
else
    min_th15=pi/2;
end

D1=[0 y(1) r(1)];
if min_th14 < min_th15
    D4=[0 y(4) r(4)];
    D5=[0 y(5) r(5)];
    C4=[0 y(4)];
    C5=[0 y(5)];
    P1=[r(1)*cos(min_th14) y(1)+r(1)*sin(min_th14)];
    P2=[-P1(1) P1(2)];
else
    D4=[0 y(5) r(5)];
    D5=[0 y(4) r(4)];
    C4=[0 y(5)];
    C5=[0 y(4)];
    P1=[r(1)*cos(min_th15) y(1)+r(1)*sin(min_th15)];
    P2=[-P1(1) P1(2)];
end

if norm(D5-D4) > 0
    [xout45,yout45] = circcirc(C4(1),C4(2),D4(3),C5(1),C5(2),D5(3));
    P45_1=[xout45(1) yout45(1)];
    P45_2=[xout45(2) yout45(2)];
    
    if xout45(1) > xout45(2)
        Q1=P45_1; Q2=P45_2;
    else
        Q1=P45_2; Q2=P45_1;
    end
else
    Q1=P1; Q2=P2;
end

if norm(D4-D1) == 0 && norm(D5-D1) == 0
    NP=[0 y(1)+r(1)];
    P1=NP; P2=NP; Q1=NP; Q2=NP;
end

if y(2) > y(3)
    D2=[0 y(2) r(2)];
    D3=[0 y(3) r(3)];
else
    D2=[0 y(3) r(3)];
    D3=[0 y(2) r(2)];
end

P1=real(P1);
P2=real(P2);
Q1=real(Q1);
Q2=real(Q2);
D1=real(D1);
D2=real(D2);
D3=real(D3);
D4=real(D4);
D5=real(D5);








%--------------------------------------------------------------------------
% fundamental_pts_pupils
%--------------------------------------------------------------------------

function [AA1,AA2,BB1,BB2,CC1,CC2,AAA1,AAA2,BBB1,BBB2,CCC1,CCC2]=...
    fundamental_pts_pupils(D1,D2)

%--------------------------------------------------------------------------
% Object: 
% This routine determines some points useful for cubature on a standard
% lens.
% The points are on the inner pupil of the lens, i.e. on the boundary of
% the union D12 of the disk D1 with the disk D2.
%
% The boundary of D12, insisting on D1, counterclockwise, starting from the
% right, on the bottom, is AA1, CC1 (possibly, if the arc AA1,BB1 has 
% maximum y-coordinate in CC1, distinct from AA1 and BB1), BB1, BB2, CC2
% (possibly, if CC1 is defined), AA2, AA1.
%
% The boundary of D12, insisting on D2, counterclockwise, starting from the
% right, on the bottom, is AAA1, CCC1 (possibly, if the arc AAA1,BBB1 has 
% maximum y-coordinate in CCC1, distinct from AAA1 and BBB1), BBB1, BBB2, 
% CCC2 (possibly, if CCC1 is defined), AAA2, AAA1.
%
% Important: fundamental points are defined even if the disks are not
% connected.
%--------------------------------------------------------------------------
% Input: 
% D1, D2: 3 x 1 vector, so that the center of the disk D1 is D1(:,1:2) and
%  the radius is D1(3), while the center of the disk D2 is D2(:,1:2) and
%  the radius is D2(3).
%--------------------------------------------------------------------------
% Output:
% AA1,AA2,BB1,BB2,CC1,CC2,AAA1,AAA2,BBB1,BBB2,CCC1,CCC2: see the object for
% their description. They are 2 x 1 vectors.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

x1=D1(1); y1=D1(2);
x2=D2(1); y2=D2(2);
r1=D1(3); r2=D2(3);

d=norm(D1(1:2)-D2(1:2));
sumr=r1+r2;

if sumr >= d
    [xout,yout] = circcirc(x1,y1,r1,x2,y2,r2);
    [v,i]=max(xout);
    j=setdiff([1 2],1);
    S1=[xout(i) yout(i)];
    S2=[xout(j) yout(j)];
    
    % ordering from lower right, counterclockwise, top disk D1 over the
    % intersection with D2.
    if y1 > S1(2)
        AA1=S1; AA2=S2;
        CC1=[x1+r1 y1]; CC2=[x1-r1 y1];
        BB1=[x1 y1+r1]; BB2=BB1;
    else
        AA1=S1; AA2=S2;
        CC1=[]; CC2=[];
        BB1=[x1 y1+r1]; BB2=BB1;
    end
    % V1=[x1+r1 S1(2); x1+r1 y1+r1; x1-r1 y1+r1; x1-r1 S1(2)];
    
    
    % ordering from lower right, counterclockwise, lower disk D2 below the
    % intersection with D1.
    if y2 < S1(2)
        AAA1=[x2 y2-r2]; AAA2=AAA1;
        CCC1=[x2+r2 y2]; CCC2=[x2-r2 y2];
        BBB1=S1; BBB2=S1;
    else
        AAA1=[x2 y2-r2]; AAA2=AAA1;
        CCC1=[]; CCC2=[];
        BBB1=S1; BBB2=S1;
    end
    
else % if sumr <=d
    
    % fprintf('\n \t ESTABLISHED POINTS DISCONNECTED')
    
    AA1=[x1 y1-r1];
    AA2=[x1 y1-r1];
    CC1=[x1+r1 y1];
    CC2=[x1-r1 y1];
    BB1=[x1 y1+r1];
    BB2=[x1 y1+r1];
    
    AAA1=[x2 y2-r2];
    AAA2=[x2 y2-r2];
    CCC1=[x2+r2 y2];
    CCC2=[x2-r2 y2];
    BBB1=[x2 y2+r2];
    BBB2=[x2 y2+r2];
    
end








%--------------------------------------------------------------------------
% northern_region
%--------------------------------------------------------------------------

function xyw=northern_region(D1,D4,D5,P1,P2,Q1,Q2,ytop,deg)

%--------------------------------------------------------------------------
% Object:
% This routine determines a cubature rule in the northern region of the
% standard lens, determined by the straight line  y=ytop, where "ytop" is the
% highest ordinate obtained on the inner pupil (i.e., with the previous
% notations, the union of the disks D2 and D3) and the intersection of the
% disks D1,D4,D5, say D145. 
% Here P1, P2 are the intersection of D4 with D1, while Q1, Q2 are the
% intersection of D4 with D5. Depending on their location, the boundary is
% defined in part on D1, in part, possibly on D4, in part, possibly, on D5.
%--------------------------------------------------------------------------
% Input: 
% D1, D4, D5: 3 x 1 vector, so that the center of the disk D1 is D1(:,1:2) 
%  and the radius is D1(3), while the center of the disk D4 is D4(:,1:2) 
%  the radius is D4(3) and finally the center of the disk D5 is D5(:,1:2) 
%  and the radius is D5(3).
% P1, P2: are the intersection of D4 with D1, P1 on the right and P2 on the
% left half plane, wrt the y-axis.
% Q1, Q2: are the intersection of D4 with D5, Q1 on the right and Q2 on the
% left half plane, wrt the y-axis.
%--------------------------------------------------------------------------
% Output:
% xyw: M x 3 vector that stores the quadrature rule with algebraic degree 
%     of precision "deg". The nodes are
%     (x(k),y(k)) where x(k)=xyw(k,1), y(k)=xyw(k,2), while the weights 
%     w(k)=xyw(k,3), i.e. they are stored in the third column of xyw.
%     We stress that the nodes are in the domain and that the weights are
%     nonnegative.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

C1=D1(:,1:2); r1=D1(:,3);
C4=D4(:,1:2); r4=D4(:,3);
C5=D5(:,1:2); r5=D5(:,3);

if Q1(2) > ytop
    if P1(2) > ytop
        %fprintf('\n \t Northern region: if P1(2) > ytop');
        NPy=D1(2)+r1;
        if (Q1(2) == NPy) && (P1(2) == NPy) % defaults if 3 circles.
            %fprintf('\n \t Northern region: if P1(2) > ytop A');
            [H2,H1] = mylinecirc(0,ytop,C1(1),C1(2),r1);
            xyw=gq_circularsegment_2017(deg,C1,r1,H1,H2);
        else
            %fprintf('\n \t Northern region: if P1(2) > ytop B');
            NPy=D1(2)+D1(3);
            if Q1(2) < NPy
                xyw1=gq_circularsegment_2017(deg,C5,r5,Q1,Q2);
            else
                xyw1=[];
            end
            if norm(Q1-P1) > 0
                xyw2=gq_circularzone_2017(deg,C4,r4,P1,Q1);
            else
                xyw2=[];
            end
            [H2,H1] = mylinecirc(0,ytop,C1(1),C1(2),r1);
            xyw3=gq_circularzone_2017(deg,C1,r1,H1,P1);
            xyw=[xyw1; xyw2; xyw3];
        end
    else
        %fprintf('\n \t Northern region: if P1(2) > ytop else');
        xyw1=gq_circularsegment_2017(deg,C5,r5,Q1,Q2);
        [H2,H1] = mylinecirc(0,ytop,C4(1),C4(2),r4);
        xyw2=gq_circularzone_2017(deg,C4,r4,H1,Q1);
        xyw=[xyw1; xyw2];
    end
else % if Q1(1) > ytop
    % fprintf('\n \t Northern region: if P1(2) > ytop else else');
    [H2,H1] = mylinecirc(0,ytop,C5(1),C5(2),r5);
    xyw=gq_circularsegment_2017(deg,C5,r5,H1,H2);
end








%--------------------------------------------------------------------------
% southern_region
%--------------------------------------------------------------------------

function xyw=southern_region(D1,D4,D5,P1,P2,Q1,Q2,ybottom,deg)

%--------------------------------------------------------------------------
% Object:
% This routine determines a cubature rule in the southern region of the
% standard lens, determined by the straight line  y=ybottom, where 
% "ybottom" is the lowest ordinate obtained on the inner pupil (i.e., 
% with the previous notations, the union of the disks D2 and D3) and the 
% intersection of the disks D1,D4,D5, say D145. 
% Here P1, P2 are the intersection of D4 with D1, while Q1, Q2 are the
% intersection of D4 with D5. Depending on their location, the boundary is
% defined in part on D1, in part, possibly on D4, in part, possibly, on D5.
%--------------------------------------------------------------------------
% Input: 
% D1, D4, D5: 3 x 1 vector, so that the center of the disk D1 is D1(:,1:2) 
%  and the radius is D1(3), while the center of the disk D4 is D4(:,1:2) 
%  the radius is D4(3) and finally the center of the disk D5 is D5(:,1:2) 
%  and the radius is D5(3).
% P1, P2: are the intersection of D4 with D1, P1 on the right and P2 on the
% left half plane, wrt the y-axis.
% Q1, Q2: are the intersection of D4 with D5, Q1 on the right and Q2 on the
% left half plane, wrt the y-axis.
%--------------------------------------------------------------------------
% Output:
% xyw: M x 3 vector that stores the quadrature rule with algebraic degree 
%     of precision "deg". The nodes are
%     (x(k),y(k)) where x(k)=xyw(k,1), y(k)=xyw(k,2), while the weights 
%     w(k)=xyw(k,3), i.e. they are stored in the third column of xyw.
%     We stress that the nodes are in the domain and that the weights are
%     nonnegative.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

C1=D1(:,1:2); r1=D1(:,3);
C4=D4(:,1:2); r4=D4(:,3);
C5=D5(:,1:2); r5=D5(:,3);

if P1(2) < ybottom
    if Q1(2) < ybottom
        [H2,H1] = mylinecirc(0,ybottom,C5(1),C5(2),r5);
        xyw1=gq_circularzone_2017(deg,C5,r5,Q1,H1);
        if norm(Q1-P1) > 0
            xyw2=gq_circularzone_2017(deg,C4,r4,P1,Q1);
        else
            xyw2=[];
        end
        xyw3=gq_circularsegment_2017(deg,C1,r1,P2,P1);
        xyw=[xyw1; xyw2; xyw3];
    else
        [H2,H1] = mylinecirc(0,ybottom,C4(1),C4(2),r4);
        xyw1=gq_circularzone_2017(deg,C4,r4,P1,H1);
        xyw2=gq_circularsegment_2017(deg,C1,r1,P2,P1);
        xyw=[xyw1; xyw2];
    end
else % if P1(2) < ybottom
    [H2,H1] = mylinecirc(0,ybottom,C1(1),C1(2),r1);
    xyw=gq_circularsegment_2017(deg,C1,r1,H2,H1);
end








%--------------------------------------------------------------------------
% central_region
%--------------------------------------------------------------------------

function xyw=central_region(Dk,D1,D4,D5,P1,P2,Q1,Q2,AA1,BB1,deg)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

ytop=BB1(2); ybottom=AA1(2);

C1=D1(:,1:2); r1=D1(:,3);
C4=D4(:,1:2); r4=D4(:,3);
C5=D5(:,1:2); r5=D5(:,3);
Ck=Dk(:,1:2); rk=Dk(:,3);

AA2=[-AA1(1) AA1(2)];
BB2=[-BB1(1) BB1(2)];

% REGION BETWEEN ycenter AND ytop.

if (Q1(2) >= ytop) && (P1(2) >= ytop)
    [DD2,DD1] = mylinecirc(0,ytop,C1(1),C1(2),r1);
    [CC2,CC1] = mylinecirc(0,ybottom,C1(1),C1(2),r1);
    xywL=new_cubature_rect_clcl(BB1,AA1,CC1,DD1,Ck,rk,C1,r1,deg);
    xywR=new_cubature_rect_clcl(AA2,BB2,DD2,CC2,Ck,rk,C1,r1,deg);
    xyw=[xywL; xywR];
end


if (Q1(2) >= ytop) && (P1(2) < ytop) && (P1(2) >= ybottom)
    [DD2,DD1] = mylinecirc(0,ytop,C4(1),C4(2),r4);
    [CC2,CC1] = mylinecirc(0,ybottom,C1(1),C1(2),r1);
    [PP2,PP1] = mylinecirc(0,P1(2),Ck(1),Ck(2),rk);
    xywR1=new_cubature_rect_clcl(BB1,PP1,P1,DD1,Ck,rk,C4,r4,deg);
    xywL1=new_cubature_rect_clcl(PP2,BB2,DD2,P2,Ck,rk,C4,r4,deg);
    if P1(2) > ybottom
        xywR2=new_cubature_rect_clcl(PP1,AA1,CC1,P1,Ck,rk,C1,r1,deg);
        xywL2=new_cubature_rect_clcl(AA2,PP2,P2,CC2,Ck,rk,C1,r1,deg);
    else
        xywR2=[]; xywL2=[];
    end
    xyw=[xywL1; xywR1; xywL2; xywR2];
end


if (Q1(2) >= ytop) && (P1(2) < ybottom)
    [DD2,DD1] = mylinecirc(0,ytop,C4(1),C4(2),r4);
    [CC2,CC1] = mylinecirc(0,ybottom,C4(1),C4(2),r4);
    xywR=new_cubature_rect_clcl(BB1,AA1,CC1,DD1,Ck,rk,C4,r4,deg);
    xywL=new_cubature_rect_clcl(AA2,BB2,DD2,CC2,Ck,rk,C4,r4,deg);
    xyw=[xywL; xywR];
end


if (Q1(2) < ytop) && (Q1(2) >= ybottom) && (P1(2) < ytop) && (P1(2) >= ybottom)
    [DD2,DD1] = mylinecirc(0,ytop,C5(1),C5(2),r5);
    [CC2,CC1] = mylinecirc(0,ybottom,C1(1),C1(2),r1);
    [PP2,PP1] = mylinecirc(0,P1(2),Ck(1),Ck(2),rk);
    [QQ2,QQ1] = mylinecirc(0,Q1(2),Ck(1),Ck(2),rk);
    xywR1=new_cubature_rect_clcl(BB1,QQ1,Q1,DD1,Ck,rk,C5,r5,deg);
    xywL1=new_cubature_rect_clcl(QQ2,BB2,DD2,Q2,Ck,rk,C5,r5,deg);
    if norm(Q1-P1) > 0
        xywR2=new_cubature_rect_clcl(QQ1,PP1,P1,Q1,Ck,rk,C4,r4,deg);
        xywL2=new_cubature_rect_clcl(PP2,QQ2,Q2,P2,Ck,rk,C4,r4,deg);
    else
        xywR2=[]; xywL2=[];
    end
    if P1(2) > ybottom
        xywR3=new_cubature_rect_clcl(PP1,AA1,CC1,P1,Ck,rk,C1,r1,deg);
        xywL3=new_cubature_rect_clcl(AA2,PP2,P2,CC2,Ck,rk,C1,r1,deg);
    else
        xywR3=[]; xywL3=[];
    end
    xyw=[xywL1; xywR1; xywL2; xywR2; xywL3; xywR3];
end


if (Q1(2) < ytop) && (Q1(2) >= ybottom) && (P1(2) < ybottom)
    [DD2,DD1] = mylinecirc(0,ytop,C5(1),C5(2),r5);
    [CC2,CC1] = mylinecirc(0,ybottom,C4(1),C4(2),r4);
    [QQ2,QQ1] = mylinecirc(0,Q1(2),Ck(1),Ck(2),rk);
    xywR1=new_cubature_rect_clcl(BB1,QQ1,Q1,DD1,Ck,rk,C5,r5,deg);
    xywL1=new_cubature_rect_clcl(QQ2,BB2,DD2,Q2,Ck,rk,C5,r5,deg);
    if Q1(2) > ybottom
        xywR2=new_cubature_rect_clcl(QQ1,AA1,CC1,Q1,Ck,rk,C4,r4,deg);
        xywL2=new_cubature_rect_clcl(AA2,QQ2,Q2,CC2,Ck,rk,C4,r4,deg);
    else
        xywR2=[]; xywL2=[];
    end
    xyw=[xywL1; xywR1; xywL2; xywR2];
end



if (Q1(2) < ybottom) && (P1(2) < ybottom)
    [DD2,DD1] = mylinecirc(0,ytop,C5(1),C5(2),r5);
    [CC2,CC1] = mylinecirc(0,ybottom,C5(1),C5(2),r5);
    xywR=new_cubature_rect_clcl(BB1,AA1,CC1,DD1,Ck,rk,C5,r5,deg);
    xywL=new_cubature_rect_clcl(AA2,BB2,DD2,CC2,Ck,rk,C5,r5,deg);
    xyw=[xywL; xywR];
end








%--------------------------------------------------------------------------
% cushion_region
%--------------------------------------------------------------------------

function xyw=cushion_region(D1,D4,D5,P1,P2,Q1,Q2,ytop,ybottom,deg)

%--------------------------------------------------------------------------
% Object:
% If the internal pupil consist of two disks D2 and D3 that are disjoined
% (here we suppose that the disk D2 is above the disk D3) in
% determining a cubature rule on the standard, one must provide as well
% a cubature rule on a "cushion" region, defined by the part of the
% intersection D145 of D1 with D4 and D5, below the straight line y=ytop
% and above the straight line y=ybottom, being ybottom the higher ordinate
% of the disk D3, and ytop the lower ordinate of the disk D2.
%--------------------------------------------------------------------------
% Input: 
% D1, D4, D5: 3 x 1 vector, so that the center of the disk D1 is D1(:,1:2) 
%  and the radius is D1(3), while the center of the disk D4 is D4(:,1:2) 
%  the radius is D4(3) and finally the center of the disk D5 is D5(:,1:2) 
%  and the radius is D5(3).
% P1, P2: are the intersection of D4 with D1, P1 on the right and P2 on the
% left half plane, wrt the y-axis.
% Q1, Q2: are the intersection of D4 with D5, Q1 on the right and Q2 on the
% left half plane, wrt the y-axis.
% ybottom: the higher ordinate of the disk D3.
% ytop: the lower ordinate of the disk D2.
% deg: algebraic degree of precision of the cubature rule, having 
% nonnegative weights.
%--------------------------------------------------------------------------
% Output:
% xyw: M x 3 vector that stores the quadrature rule with algebraic degree 
%     of precision "deg". The nodes are
%     (x(k),y(k)) where x(k)=xyw(k,1), y(k)=xyw(k,2), while the weights 
%     w(k)=xyw(k,3), i.e. they are stored in the third column of xyw.
%     We stress that the nodes are in the domain and that the weights are
%     nonnegative.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

% 1.
if (Q1(2) >= ytop) && (P1(2) >= ytop)
    % fprintf('\n \t type 1');
    C1=D1(1:2); r1=D1(3);
    [DD2,DD1] = mylinecirc(0,ytop,C1(1),C1(2),r1);
    [CC2,CC1] = mylinecirc(0,ybottom,C1(1),C1(2),r1);
    
    xyw=gq_circularzone_2017(deg,C1,r1,CC1,DD1);
end

%2.
if (Q1(2) >= ytop) && (P1(2) <= ytop) && (P1(2) > ybottom)
    % fprintf('\n \t type 2');
    C1=D1(1:2); r1=D1(3);
    C4=D4(1:2); r4=D4(3);
    [DD2,DD1] = mylinecirc(0,ytop,C4(1),C4(2),r4);
    [CC2,CC1] = mylinecirc(0,ybottom,C1(1),C1(2),r1);
    
    if P1(2) < ytop
        xywT=gq_circularzone_2017(deg,C4,r4,P1,DD1);
    else
        xywT=[];
    end
    
    xywL=gq_circularzone_2017(deg,C1,r1,CC1,P1);
    
    xyw=[xywT; xywL];
end

% 3.
if (Q1(2) >= ytop) && (P1(2) <= ybottom)
    % fprintf('\n \t type 3');
    C4=D4(1:2); r4=D4(3);
    [DD2,DD1] = mylinecirc(0,ytop,C4(1),C4(2),r4);
    [CC2,CC1] = mylinecirc(0,ybottom,C4(1),C4(2),r4);
    
    xyw=gq_circularzone_2017(deg,C4,r4,CC1,DD1);
    
end

% 4.
if (Q1(2) <= ytop) && (Q1(2) > ybottom) && (P1(2) <= ytop) && (P1(2) > ybottom)
    % fprintf('\n \t type 4');
    C1=D1(1:2); r1=D1(3);
    C4=D4(1:2); r4=D4(3);
    C5=D5(1:2); r5=D5(3);
    
    [DD2,DD1] = mylinecirc(0,ytop,C5(1),C5(2),r5);
    [CC2,CC1] = mylinecirc(0,ybottom,C1(1),C1(2),r1);
    
    if Q1(2) < ytop
        xywT=gq_circularzone_2017(deg,C5,r5,Q1,DD1);
    else
        xywT=[];
    end
    
    if Q1(2) > P1(2)
        xywM=gq_circularzone_2017(deg,C4,r4,P1,Q1);
    else
        xywM=[];
    end
    
    if P1(2) > ybottom
        xywL=gq_circularzone_2017(deg,C1,r1,CC1,P1);
    else
        xywL=[];
    end
    
    xyw=[xywT; xywM; xywL];
end

% 5.
if (Q1(2) <= ytop) && (Q1(2) > ybottom) && (P1(2) <= ybottom)
    % fprintf('\n \t type 5');
    C5=D5(1:2); r5=D5(3);
    C4=D4(1:2); r4=D4(3);
    
    [DD2,DD1] = mylinecirc(0,ytop,C5(1),C5(2),r5);
    [CC2,CC1] = mylinecirc(0,ybottom,C4(1),C4(2),r4);
    
    if Q1(2) < ytop
        xywT=gq_circularzone_2017(deg,C5,r5,Q1,DD1);
    else
        xywT=[];
    end
    
    xywL=gq_circularzone_2017(deg,C4,r4,CC1,Q1);
    
    xyw=[xywT; xywL];
end

% 6.
if (Q1(2) <= ybottom)
    % fprintf('\n \t type 6');
    C5=D5(1:2); r5=D5(3);
    [DD2,DD1] = mylinecirc(0,ytop,C5(1),C5(2),r5);
    [CC2,CC1] = mylinecirc(0,ybottom,C5(1),C5(2),r5);
    
    xyw=gq_circularzone_2017(deg,C5,r5,CC1,DD1);
end








%--------------------------------------------------------------------------
% mylinecirc
%--------------------------------------------------------------------------

function [ptL,ptR]=mylinecirc(m,q,centerx,centery,radius)

%--------------------------------------------------------------------------
% Object:
% In this variant of Matlab built-in function linecirc.m, we produce two 
% points ptL and ptR, where the abscissa of ptL is inferior or equal to
% that of ptR.
% We recall that linecirc.m computes the intersection of the circle with
% center (centerx,centery) and radius "radius" with the straight line of
% equations y=mx+q.
%--------------------------------------------------------------------------
% Input: 
% m,q: scalars defining the straight line y=mx+q.
% centerx,centery: coordinates of the center of the circle.
% radius: radius of the circle.
%--------------------------------------------------------------------------
% Output:
% ptL,ptR: intersection of the line and of the circle defined by inputs. 
% The abscissa of ptL is inferior or equal to that of ptR.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

[xout,yout] = linecirc(m,q,centerx,centery,radius);

if xout(1) < xout(2)
    ptL=[xout(1) yout(1)];
    ptR=[xout(2) yout(2)];
else
    ptL=[xout(2) yout(2)];
    ptR=[xout(1) yout(1)];
end








%--------------------------------------------------------------------------
% check_pupils
%--------------------------------------------------------------------------

function [y,r,flag]=check_pupils(y,r)

%--------------------------------------------------------------------------
% Object:
% This routine checks if the internal pupils defined via y, r, i.e.  
% D1=[0 y(2) r(2)] and D2=[0 y(3) r(3)], are pathological, i.e. disjoined or
% one inside the other one. 
% In general the flag variable is 1, being 0 if the pupils are disjoined.
% In case one pupil is inside the other one, the procedure sets y(2)=y(3),
% r(2)=r(3), depending on the "larger" disk in the pupil.
%--------------------------------------------------------------------------
% Input: 
% y: n x 1 vector (with n <= 5) such that the disk Dk has center (0,y(k)).
% r: n x 1 vector (with n <= 5) such that the disk Dk has radius r(k).
%--------------------------------------------------------------------------
% Output:
% y: n x 1 vector (with n <= 5) such that the disk Dk has center (0,y(k)).
% r: n x 1 vector (with n <= 5) such that the disk Dk has radius r(k).
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

flag=1;

D1=[0 y(2) r(2)];
D2=[0 y(3) r(3)];

C1=D1(1:2); r1=D1(3);
C2=D2(1:2); r2=D2(3);

% check if one disk contains the other one.
if r1 > r2
    chk_pts=[C2+[r2 0]; C2+[0 r2]; C2+[-r2 0]; C2+[0 -r2]];
    flagL=point_in_disk(chk_pts,D1);
    
    if sum(flagL) == size(chk_pts,1)
        fprintf('\n \t Warning: Disk 3 is contained in disk 2. Normalizing.');
        y(3)=y(2); r(3)=r(2);
        return;
    end
    
    if sum(flagL) == 0
        flag=0;
        fprintf('\n \t WARNING: The disks are disconnected.');
        
        if D1(2) < D2(2)
            [y(3),y(2)]=deal(y(2),y(3));
            [r(3),r(2)]=deal(r(2),r(3));
        end
        return;
    end
    
else
    chk_pts=[C1+[r1 0]; C1+[0 r1]; C1+[-r1 0]; C1+[0 -r1]];
    flagL=point_in_disk(chk_pts,D2);
   
    if sum(flagL) == size(chk_pts,1)
        fprintf('\n \t Warning: Disk 2 is contained in disk 3. Normalizing.');
        y(2)=y(3); r(2)=r(3);
        return;
    end
    
    if sum(flagL) == 0
        flag=0;
        
        if D1(2) < D2(2)
            [y(3),y(2)]=deal(y(2),y(3));
            [r(3),r(2)]=deal(r(2),r(3));
        end
        fprintf('\n \t WARNING: The disks are disconnected.');
        return;
    end
    
end








%--------------------------------------------------------------------------
% point_in_disk
%--------------------------------------------------------------------------

function flag=point_in_disk(P,D)

%--------------------------------------------------------------------------
% Object:
% Given a point P in the plane and a disk D, the procedure checks if P is
% inside the closed disk D (i.e. also on the in the boundary).
%--------------------------------------------------------------------------
% Input: 
% P: 2 x 1 vector.
% D: 3 x 1 vector, so that the center of the disk D is D(:,1:2) 
%    and the radius is D(3).
%--------------------------------------------------------------------------
% Output:
% flag: 0 if the point is outside the closed disk D, 1 otherwise.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

C=D(1:2); r=D(3);
dv=sqrt( (P(:,1)-C(1)).^2 + (P(:,2)-C(2)).^2 );
flag= (dv <= r);








%--------------------------------------------------------------------------
% troubleshooting
%--------------------------------------------------------------------------

function [P1,P2,Q1,Q2,D1,D2,D3,D4,D5,flag]=troubleshooting(y,r)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

if size(y,1) < size(y,2)
    y=y';
end

if size(r,1) < size(r,2)
    r=r';
end

if length(y) == 2
    fprintf('\n \t Warning: Only 2 active disks. Normalizing.');
    y=[y; y(2); y(1); y(1)]; r=[r; r(2); r(1); r(1)];
end

if length(y) == 3
    fprintf('\n \t Warning: Only 3 active disks. Normalizing.');
    y=[y; y(1); y(1)]; r=[r; r(1); r(1)];
end

if length(y) == 4
    fprintf('\n \t Warning: Only 4 active disks. Normalizing.');
    y=[y; y(4)]; r=[r; r(4)];
end

% Checking radii.
if r(1) == 0
    fprintf('\n \t Fatal error: radius of first disk must be larger than 0');
    P1=[]; P2=[]; Q1=[]; Q2=[]; D1=[]; D2=[]; D3=[]; D4=[]; D5=[];
    return;
end

if r(2) == 0
    fprintf('\n \t Fatal error: radius of second disk must be larger than 0');
    P1=[]; P2=[]; Q1=[]; Q2=[]; D1=[]; D2=[]; D3=[]; D4=[]; D5=[];
    return;
end

if r(3) == 0
    fprintf('\n \t Warning: radius of third disk must be larger than 0, Using simple pupil.');
    y(3)=y(2); r(3)=r(2);
end

if r(4) == 0
    fprintf('\n \t  Warning: radius of fourth disk must be larger than 0, Using third fourth equal to first disk.');
    y(4)=y(1); r(4)=r(1);
end

if r(5) == 0
    fprintf('\n \t  Warning: radius of fifth disk must be larger than 0, Using fifth disk equal to fourth disk.');
    y(5)=y(4); r(5)=r(4);
end

[y,r,flag]=check_pupils(y,r);

% COMPUTE FUNDAMENTAL POINTS.
[P1,P2,Q1,Q2,D1,D2,D3,D4,D5]=compute_fundamental_pts(y,r);

flagL=point_in_disk(Q1,D1);
if flagL == 0
    if isnan(Q1(1)) ||  isnan(Q1(2))
        fprintf('\n \t WARNING 2: No intersection between D4 and D5. \n \n');
        D5=D4; Q1=P1; Q2=P2;
    else
        fprintf('\n \t WARNING 2: D5 intersects D4 outside D1 \n \n');
        D5=D4; Q1=P1; Q2=P2;
    end
    
end








%--------------------------------------------------------------------------
% cubature_trapezoid
%--------------------------------------------------------------------------

function xyw=cubature_trapezoid(ade,vertices)

%--------------------------------------------------------------------------
% Object:
% Quadrature on a trapezoid with vertices "vertices" and algebraic degree
% of precision "ade".
% Vertices are defined without the requirement that the first and the last 
% vertices are equal.
%--------------------------------------------------------------------------
% Input: 
% ade: algebraic degree of precision of the cubature rule.
% vertices: vertices of the trapezoid, defined counterclockwise, without
% the restriction requirement that the first and the last vertices are 
% equal.
%--------------------------------------------------------------------------
% Output:
% xyw: M x 3 vector that stores the quadrature rule with algebraic degree 
%     of precision "deg". The nodes are
%     (x(k),y(k)) where x(k)=xyw(k,1), y(k)=xyw(k,2), while the weights 
%     w(k)=xyw(k,3), i.e. they are stored in the third column of xyw.
%     We stress that the nodes are in the domain and that the weights are
%     nonnegative.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

T1=vertices(1:3,:);
xyw1=stroud_conical_rules(ade,T1);

T2=vertices([3 4 1],:);
xyw2=stroud_conical_rules(ade,T2);

xyw=[xyw1; xyw2];








%--------------------------------------------------------------------------
% cubature_trapezoid
%--------------------------------------------------------------------------

function [xyw]=stroud_conical_rules(ade,vertices)

%--------------------------------------------------------------------------
% Object:
%--------------------------------------------------------------------------
% Input: 
% ade: ALGEBRAIC DEGREE OF EXACTNESS.
% vertices: 3 x 2 MATRIX OF VERTICES OF THE SIMPLEX.
%--------------------------------------------------------------------------
% Output:
% xyw: NODES AND WEIGHTS OF STROUD CONICAL RULE TYPE OF ADE ade ON THE 
%      SIMPLEX WITH VERTICES vertices.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------


if nargin < 2 % SEE LYNESS, COOLS,
    % "A survey on numerical cubature over triangles", p.4.
    vertices=[0 0; 1 0; 1 1];
end

[xw]=stroud_conical_rules_ref(ade+1);
bar_coord=[1-xw(:,1) xw(:,1)-xw(:,2) xw(:,2)];
xy=bar_coord*vertices;

A = polyarea(vertices(:,1),vertices(:,2));

w=xw(:,3)*A*2;
xyw=[xy w];








%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

function [xw]=stroud_conical_rules_ref(ade)

% SEE LYNESS, COOLS, "A survey on numerical cubature over triangles", p.4.

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

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








%--------------------------------------------------------------------------
% gauss_jacobi
%--------------------------------------------------------------------------

function [x,w]=gauss_jacobi(N,a,b,gl)

%--------------------------------------------------------------------------
% Object:
% GAUSS-JACOBI (LOBATTO) RULE ON (-1,1).
%--------------------------------------------------------------------------
% Input: 
% N: determines the number of points of the rule, N if gl=0, N+2 if gl=1.
% a,b ARE THE GAUSS-JACOBI EXPONENTS.
% gl: 0: GAUSS POINTS. 1: GAUSS-LOBATTO POINTS.
%--------------------------------------------------------------------------
% Output:
% x, w ARE COLUMN VECTORS OF NODES AND WEIGHTS.
%      THE LENGTH OF x AND w IS "N" IF gl=0, "N+2" IF "gl=1".
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------


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
% r_jacobi
%--------------------------------------------------------------------------

function ab=r_jacobi(N,a,b)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Dirk Laurie and Walter Gautschi.
%--------------------------------------------------------------------------

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








%--------------------------------------------------------------------------
% gauss
%--------------------------------------------------------------------------

function xw=gauss(N,ab)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Dirk Laurie and Walter Gautschi.
%--------------------------------------------------------------------------

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








%--------------------------------------------------------------------------
% lobatto_jacobi
%--------------------------------------------------------------------------

function xw=lobatto_jacobi(N,a,b)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Dirk Laurie and Walter Gautschi.
%--------------------------------------------------------------------------

if nargin<2, a=0; end
if nargin<3, b=a; end
ab=r_jacobi(N+2,a,b);
ab(N+2,1)=(a-b)/(2*N+a+b+2);
ab(N+2,2)=4*(N+a+1)*(N+b+1)*(N+a+b+1)/((2*N+a+b+1)*(2*N+a+b+2)^2);
xw=gauss(N+2,ab);








%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

function xyw=new_cubature_rect_clcl(A,B,C,D,C1,r1,C2,r2,deg,method)

%--------------------------------------------------------------------------
% Object:
% A,B,C,D: vertices, where A,B are the inner ones. Counterclockwise.
%
%
% Important: the tangents on the circles passing for A,B with center C1 and
% radius r1 must have form y=mx+q with m of costant sign. Similarly those
% passing for C,D with center C2 and radius r2.

% First method: the domain can be subdivided as two sectors and a rectangle.
% xyw=first_method(A,B,C,D,C1,r1,C2,r2,deg);

% Second method: the domain can be subdivided as two sectors.
% xyw=second_method(A,B,C,D,C1,r1,C2,r2,deg);
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------


if nargin < 10
    method=3;
end

switch method
    case 1
        xyw=first_method(A,B,C,D,C1,r1,C2,r2,deg);
    case 2
        xyw=second_method(A,B,C,D,C1,r1,C2,r2,deg);
    otherwise
        xyw=third_method(A,B,C,D,C1,r1,C2,r2,deg);
end








%--------------------------------------------------------------------------
% first_methods
%--------------------------------------------------------------------------

function xyw=first_method(A,B,C,D,C1,r1,C2,r2,deg)

%--------------------------------------------------------------------------
% Object:
% Cubature rule over a region where AB and CD are arcs and BC, DA segments
% parallel to the x-axis.
% AB is the arc relative to (C1,r1) and CD is the arc relative to (C2,r2).
% The boundary of the region is obtained by running counterclockwise
% 1. the arc AB,
% 2. the segment BC,
% 3. the arc CD,
% 4. the segment DA.
% Let y=mx+q be the tangent in the generic point P of each arc. Here we
% suppose that m is of constant sign on each arc, varying P, possibly null.
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------


xyw=[];

% 1. right halfplane, upper halfplane.
if (A(2) > B(2)) && (D(2) > C(2)) && (A(1) < B(1)) && (D(1) < C(1))
    if B(1) < D(1)
        P=[B(1) A(2)];
        CC1=C1; rr1=r1; th1=B; th2=A;
        xyw1=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        AA=P; BB=B; CC=[D(1) C(2)]; DD=D;
        vertices=[AA; BB; CC; DD];
        xyw2=cubature_trapezoid(deg,vertices);
        
        P=CC;
        CC1=C2; rr1=r2; th1=C; th2=D;
        xyw3=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        xyw=[xyw1; xyw2; xyw3];
        
    end
end

% 2. right halfplane, lower halfplane.
if (A(2) > B(2)) && (A(1) > B(1)) && (D(1) > C(1)) && (D(2) > C(2))
    if A(1) < C(1)
        P=[A(1) B(2)];
        CC1=C1; rr1=r1; th1=B; th2=A;
        xyw1=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        AA=A; BB=P; CC=C; DD=[C(1) D(2)];
        vertices=[AA; BB; CC; DD];
        xyw2=cubature_trapezoid(deg,vertices);
        
        P=DD;
        CC1=C2; rr1=r2; th1=C; th2=D;
        xyw3=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        xyw=[xyw1; xyw2; xyw3];
        
    end
end


% 3. left halfplane, upper halfplane.
if (A(2) < B(2)) && (A(1) < B(1)) && (D(1) < C(1)) && (D(2) < C(2))
    if C(1) < A(1)
        P=[A(1) B(2)];
        CC1=C1; rr1=r1; th1=B; th2=A;
        xyw1=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        AA=A; BB=P; CC=C; DD=[C(1) D(2)];
        vertices=[AA; BB; CC; DD];
        xyw2=cubature_trapezoid(deg,vertices);
        
        P=DD;
        CC1=C2; rr1=r2; th1=C; th2=D;
        xyw3=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        xyw=[xyw1; xyw2; xyw3];
        
    end
end

% 4. left halfplane, lower halfplane.
if (A(2) < B(2)) && (A(1) > B(1)) && (D(1) > C(1)) && (D(2) < C(2))
    if B(1) > D(1)
        P=[B(1) A(2)];
        CC1=C1; rr1=r1; th1=B; th2=A;
        xyw1=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        AA=P; BB=B; CC=[D(1) C(2)]; DD=D;
        vertices=[AA; BB; CC; DD];
        xyw2=cubature_trapezoid(deg,vertices);
        
        P=CC;
        CC1=C2; rr1=r2; th1=C; th2=D;
        xyw3=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        xyw=[xyw1; xyw2; xyw3];
        
    end
end







% 5. other regions. (=)
if (A(2) > B(2)) && (A(1) > B(1)) && (D(1) < C(1)) && (D(2) > C(2))
    if A(1) < D(1)
        P=[A(1) B(2)];
        CC1=C1; rr1=r1; th1=A; th2=B;
        xyw1=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        AA=A; BB=P; CC=[D(1) C(2)]; DD=D;
        vertices=[AA; BB; CC; DD];
        xyw2=cubature_trapezoid(deg,vertices);
        
        P=CC;
        CC1=C2; rr1=r2; th1=C; th2=D;
        xyw3=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        xyw=[xyw1; xyw2; xyw3];
        
    end
end







% 6. other regions. (=)
if (A(1) < B(1)) && (A(2) > B(2)) && (C(1) < D(1)) && (C(2) < D(2))
    if B(1) < C(1)
        P=[B(1) A(2)];
        CC1=C1; rr1=r1; th1=A; th2=B;
        xyw1=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        AA=P; BB=B; CC=C; DD=[C(1) D(2)];
        vertices=[AA; BB; CC; DD];
        xyw2=cubature_trapezoid(deg,vertices);
        
        P=DD;
        CC1=C2; rr1=r2; th1=C; th2=D;
        xyw3=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        xyw=[xyw1; xyw2; xyw3];
        
    end
end





% 7. other regions. (=)
if (A(1) < B(1)) && (A(2) > B(2)) && (C(1) < D(1)) && (C(2) < D(2))
    if B(1) < C(1)
        P=[B(1) A(2)];
        CC1=C1; rr1=r1; th1=A; th2=B;
        xyw1=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        AA=P; BB=B; CC=C; DD=[C(1) D(2)];
        vertices=[AA; BB; CC; DD];
        xyw2=cubature_trapezoid(deg,vertices);
        
        P=DD;
        CC1=C2; rr1=r2; th1=C; th2=D;
        xyw3=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        xyw=[xyw1; xyw2; xyw3];
        
    end
end





% 8. other regions. (=)
if (A(1) > B(1)) && (A(2) > B(2)) && (C(1) > D(1)) && (C(2) < D(2))
    if A(1) < D(1)
        P=[A(1) B(2)];
        CC1=C1; rr1=r1; th1=A; th2=B;
        xyw1=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        AA=A; BB=P; CC=[D(1) C(2)]; DD=D;
        vertices=[AA; BB; CC; DD];
        xyw2=cubature_trapezoid(deg,vertices);
        
        P=CC;
        CC1=C2; rr1=r2; th1=C; th2=D;
        xyw3=gq_circularsector_2017(P,CC1,rr1,th1,th2,deg);
        
        xyw=[xyw1; xyw2; xyw3];
        
    end
end








%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

function xyw=second_method(A,B,C,D,C1,r1,C2,r2,deg)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

xyw=[];

% See if if the segment AC is completely inside the domain.
x0=A(1); y0=A(2); x1=C(1); y1=C(2);
m=(y1-y0)/(x1-x0); q=y0-m*x0;

% part1 of the test.
[xout,yout] = linecirc(m,q,C1(1),C1(2),r1);
P0=[xout(1) yout(1)]; P1=[xout(2) yout(2)];

if norm(P0-A) > norm(P1-A)
    Ptest=P0;
else
    Ptest=P1;
end
% flag0=new_pointindomain_clcl(Ptest,A,B,C,D,C1,r1,C2,r2); % 1: ko.
flag0=new_pointinarc(Ptest,C1,r1,B,A);

% part 2 of the test
[xout,yout] = linecirc(m,q,C2(1),C2(2),r2);
P0=[xout(1) yout(1)]; P1=[xout(2) yout(2)];
if norm(P0-C) > norm(P1-C)
    Ptest=P0;
else
    Ptest=P1;
end
% flag1=new_pointindomain_clcl(Ptest,A,B,C,D,C1,r1,C2,r2); % 1: ko.
flag1=new_pointinarc(Ptest,C2,r2,C,D);

% Determine rule.

if flag0+flag1 == 0
    % fprintf('\n \t CASE AC');
    xyw1=gq_circularsector_2017(C,C1,r1,B,A,deg);
    xyw2=gq_circularsector_2017(A,C2,r2,C,D,deg);
    xyw=[xyw1; xyw2];
    flag=new_pointindomain_clcl(xyw(:,1:2),A,B,C,D,C1,r1,C2,r2);
    return;
end



% In case the previous procedure does not work, seee if if the segment BD
% is completely inside the domain.
x0=B(1); y0=B(2); x1=D(1); y1=D(2);
m=(y1-y0)/(x1-x0); q=y0-m*x0;

% part1 of the test.
[xout,yout] = linecirc(m,q,C1(1),C1(2),r1);
P0=[xout(1) yout(1)]; P1=[xout(2) yout(2)];
if norm(P0-B) > norm(P1-B)
    Ptest=P0;
else
    Ptest=P1;
end
flag0=new_pointinarc(Ptest,C1,r1,B,A);
% flag0=new_pointindomain_clcl(Ptest,A,B,C,D,C1,r1,C2,r2); % 1: ko.

% part 2 of the test
[xout,yout] = linecirc(m,q,C2(1),C2(2),r2);
P0=[xout(1) yout(1)]; P1=[xout(2) yout(2)];
if norm(P0-D) > norm(P1-D)
    Ptest=P0;
else
    Ptest=P1;
end
% flag1=new_pointindomain_clcl(Ptest,A,B,C,D,C1,r1,C2,r2); % 1: ko.
flag1=new_pointinarc(Ptest,C2,r2,C,D);
% Determine rule.
if flag0+flag1 == 0
    % fprintf('\n \t CASE BD');
    xyw1=gq_circularsector_2017(D,C1,r1,B,A,deg);
    xyw2=gq_circularsector_2017(B,C2,r2,C,D,deg);
    xyw=[xyw1; xyw2];
    flag=new_pointindomain_clcl(xyw(:,1:2),A,B,C,D,C1,r1,C2,r2);
    return;
end








%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

function xyw=third_method(A,B,C,D,C1,r1,C2,r2,deg)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

subtype=2; % method 1 not always work.

switch subtype
    case 1
        xyw=first_method(A,B,C,D,C1,r1,C2,r2,deg);
    case 2
        xyw=second_method(A,B,C,D,C1,r1,C2,r2,deg);
end



if size(xyw,1) == 0
    % fprintf('\n \t THIRD METHOD ITERATIONS: ');
    subdivide_again=1; iter=0;
    
    while (subdivide_again == 1)
        [xywL,A,B,C,D]=third_method_iteration(A,B,C,D,...
            C1,r1,C2,r2,deg);
        xyw=[xyw; xywL];
        switch subtype
            case 1
                xywL=first_method(A,B,C,D,C1,r1,C2,r2,deg);
            case 2
                xywL=second_method(A,B,C,D,C1,r1,C2,r2,deg);
        end
        if size(xywL,1) > 0
            xyw=[xyw; xywL];
            subdivide_again=0;
        end
        iter=iter+1;
        fprintf(' %3.0f',iter);
    end
end








%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

function [xywL,A,B,C,D]=third_method_iteration(A,B,C,D,C1,r1,C2,r2,deg)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

% 1. right halfplane, upper halfplane.
if (A(2) > B(2)) && (D(2) > C(2)) && (A(1) < B(1)) && (D(1) < C(1))
    
    % fprintf('\n \t SUBCASE 1');
    
    AA(1)=B(1);
    cc=(AA(1)-C2(1))/r2; ss=sqrt(1-cc^2);
    AA(2)=C2(2)+r2*ss;
    
    CC(2)=AA(2);
    ss=(CC(2)-C1(2))/r1; cc=sqrt(1-ss^2);
    CC(1)=C1(1)+r1*cc;
    
    xywL1=gq_circularsector_2017(B,C2,r2,C,AA,deg);
    xywL2=gq_circularsector_2017(AA,C1,r1,B,CC,deg);
    xywL=[xywL1; xywL2];
    B=CC; C=AA;
    
end



% 2. right halfplane, lower halfplane.
if (A(2) > B(2)) && (A(1) > B(1)) && (D(1) > C(1))
    
    % fprintf('\n \t SUBCASE 2');
    
    AA(1)=A(1);
    cc=(AA(1)-C2(1))/r2; ss=-sqrt(1-cc^2);
    AA(2)=C2(2)+r2*ss;
    
    CC(2)=AA(2);
    ss=(CC(2)-C1(2))/r1; cc=sqrt(1-ss^2);
    CC(1)=C1(1)+r1*cc;
    
    xywL1=gq_circularsector_2017(A,C2,r2,AA,D,deg);
    xywL2=gq_circularsector_2017(AA,C1,r1,CC,A,deg);
    xywL=[xywL1; xywL2];
    A=CC; D=AA;
    
end


% 3. left halfplane, upper halfplane.
if (A(2) < B(2)) && (A(1) < B(1)) && (D(1) < C(1))
    
    % fprintf('\n \t SUBCASE 3');
    
    AA(1)=A(1);
    cc=(AA(1)-C2(1))/r2; ss=sqrt(1-cc^2);
    AA(2)=C2(2)+r2*ss;
    
    CC(2)=AA(2);
    ss=(CC(2)-C1(2))/r1; cc=-sqrt(1-ss^2);
    CC(1)=C1(1)+r1*cc;
    
    xywL1=gq_circularsector_2017(A,C2,r2,AA,D,deg);
    xywL2=gq_circularsector_2017(AA,C1,r1,CC,A,deg);
    xywL=[xywL1; xywL2];
    A=CC; D=AA;
    
end

% 4. left halfplane, lower halfplane.
if (A(2) < B(2)) && (A(1) > B(1)) && (D(1) > C(1))
    
    % fprintf('\n \t SUBCASE 4');
    
    AA(1)=B(1);
    cc=(AA(1)-C2(1))/r2; ss=-sqrt(1-cc^2);
    AA(2)=C2(2)+r2*ss;
    
    CC(2)=AA(2);
    ss=(CC(2)-C1(2))/r1; cc=-sqrt(1-ss^2);
    CC(1)=C1(1)+r1*cc;
    
    xywL1=gq_circularsector_2017(B,C2,r2,C,AA,deg);
    xywL2=gq_circularsector_2017(AA,C1,r1,B,CC,deg);
    xywL=[xywL1; xywL2];
    B=CC; C=AA;
end








%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

function xyw=gq_rectangular_sector_2017(deg,C,r,P1,P2,Q1,Q2)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

% Input:
% deg: degree of the cubature rule.
% P1,P2: vertices of the linear part (counterclockwise direction).
% Q1,Q2: vertices of the circular part  (counterclockwise direction).
% P1,P2,Q1,Q2 are 1 x 2 vectors.
%
% Output:
% xyw: cubature nodes (x,y) and weights "w" stored as [x y w].

polygon_vertices=[P1; P2; Q1; Q2; P1];
if polyarea(polygon_vertices(:,1),polygon_vertices(:,2)) > 0
    xywP=cubature_trapezoid(deg,polygon_vertices);
    % xywP=cubature_polygons(polygon_vertices,deg);
else
    xywP=[];
end
xywC=gq_circularsegment_2017(deg,C,r,Q1,Q2);
xyw=[xywP; xywC];








%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

function xyw=gq_circularsegment_2017(deg,CC,r,th1,th2)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

if length(th1) == 2
    P=th1-CC;
    th1=cart2pol(P(1),P(2));
end

if length(th2) == 2
    P=th2-CC; th2=cart2pol(P(1),P(2));
end

if th1 > th2
    th2=th2+2*pi;
end


omega=(th2-th1)/2;
xyw0 = gqcircsegm(deg,omega,r);
X0=xyw0(:,1); Y0=xyw0(:,2); W=xyw0(:,3);
[TH0,R]=cart2pol(X0,Y0);
TH=TH0+(th2+th1)/2;
X=CC(1)+R.*cos(TH); Y=CC(2)+R.*sin(TH);
xyw=[X Y W];








%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

function xyw = gqcircsegm(n,omega,r)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

% by Gaspare Da Fies and Marco Vianello, University of Padova
% 2 Dec 2011

% computes the nodes and weights of a product gaussian formula
% on a circular segment of a disk centered at the origin
% with angles in [-omega,omega]

% uses the routines:
%
% r_jacobi.m, gauss.m from
% www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
%
% trigauss.m
% http://www.math.unipd.it/~marcov/mysoft/trigauss.m

% input:
% n: algebraic degree of exactness
% omega: half-length of the angular interval, 0<omega<=pi
% r: radius of the disk

% output:
% xyw: (ceil((n+1)/2) x ceil((n+2)/2) x 3 array of (xnodes,ynodes,weights)


% trigonometric gaussian formula on the arc
tw=trigauss(n+2,-omega,omega);

% algebraic gaussian formula on [-1,1]
ab=r_jacobi(ceil((n+1)/2),0,0);
xw=gauss(ceil((n+1)/2),ab);

% creating the grid
Ltw1=size(tw,1); M1=floor(Ltw1/2);
[t,theta]=meshgrid(xw(:,1),tw(1:M1,1));
[w1,w2]=meshgrid(xw(:,2),tw(1:M1,2));

% [t,theta]=meshgrid(xw(:,1),tw(1:ceil((n+2)/2),1));
% [w1,w2]=meshgrid(xw(:,2),tw(1:ceil((n+2)/2),2));

% nodal cartesian coordinates and weights
s=sin(theta(:));
xyw(:,1)=r*cos(theta(:));
xyw(:,2)=r*t(:).*s;
xyw(:,3)=r^2*s.^2.*w1(:).*w2(:);








%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

function xyw=gq_circularsector_2017(P,C,r,th1,th2,deg)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

% Quadrature rule on the circular sector. P is the point external to the
% disk with center C and radius r.

% If th1, th2, phi1, phi2 are cartesian points, we convert them into polar
% coordinates relatively to C1 for th1, th2 and C2 for phi1,phi2.
[th1,th2]=initial_check(C,r,th1,th2);

xyw=gq_sect(P,C,r,th1,th2,deg);








%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

function [th1,th2]=initial_check(C,r,th1,th2)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

% working with row vectors.

if size(C,1) == 2
    C=C';
end

if length(th1) > 1
    P=th1-C;
    th1=cart2pol(P(1),P(2));
end

if length(th2) > 1
    P=th2-C;
    th2=cart2pol(P(1),P(2));
end

if th2 < th1
    th2=th2+2*pi;
end








%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

function xyw=gq_sect(C1,C2,r2,th1,th2,deg)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

% sector with center C1 and arc on the disk with center C2 and arc spanned
% by the angles th1 and th2 with th2 > th1.

AA1=[0 0]; BB1=[0 0]; CC1=C1; % point
AA2=[r2 0]; BB2=[0 r2]; CC2=C2; % arc
A=[AA1; AA2]; B=[BB1; BB2]; C=[CC1; CC2];
alpha=th1;
beta=th2;
xyw = gqellblend(deg,A,B,C,alpha,beta);








%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

function xyw=gq_circularzone_2017(deg,CC,r,th1,th2)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

if length(th1) == 2
    P=th1-CC; th1=cart2pol(P(1),P(2));
end

if length(th2) == 2
    P=th2-CC; th2=cart2pol(P(1),P(2));
end

if th1 > th2
    th2=th2+2*pi;
end

alpha=th1; beta=th2;

A1=[r 0]; B1=[0 r]; C1=[CC(1) CC(2)];
A2=[-r 0]; B2=[0 r]; C2=[CC(1) CC(2)];

A=[A1; A2]; B=[B1; B2]; C=[C1; C2];

xyw = gqellblend(deg,A,B,C,alpha,beta);








%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

function xyw = gqcirczone(n,alpha,beta,r)

%--------------------------------------------------------------------------
% Object:
%
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

% by Gaspare Da Fies and Marco Vianello, University of Padova
% 8 Nov 2011

% computes the nodes and weights of a product gaussian formula
% on a circular zone of a disk centered at the origin namely, the
% portion of the disk between the lines x=r*cos(alpha) and x=r*cos(beta)

% in the special case alpha=0 the zone is a circular segment

% uses the routines:
%
% r_jacobi.m, gauss.m from
% www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
%
% trigauss.m
% http://www.math.unipd.it/~marcov/mysoft/trigauss.m

% input:
% n: algebraic degree of exactness
% [alpha,beta]: angular interval, 0<=alpha<beta<=pi
% r: radius of the disk

% output:
% xyw: (ceil((n+1)/2) x (n+3)) x 3 array of (xnodes,ynodes,weights)


% trigonometric gaussian formula on the arc
tw=trigauss(n+2,alpha,beta);

% algebraic gaussian formula on [-1,1]
ab=r_jacobi(ceil((n+1)/2),0,0);
xw=gauss(ceil((n+1)/2),ab);

% creating the polar grid
[t,theta]=meshgrid(xw(:,1),tw(:,1));
[w1,w2]=meshgrid(xw(:,2),tw(:,2));

% nodal cartesian coordinates and weights
s=sin(theta(:));
xyw(:,1)=r*cos(theta(:));
xyw(:,2)=r*t(:).*s;
xyw(:,3)=r^2*s.^2.*w1(:).*w2(:);








%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

function flag=new_pointinarc(P,C,r,A,B)

%--------------------------------------------------------------------------
% Object:
% Let X be the circle with center "C" and radius "r".
%--------------------------------------------------------------------------
% Input: 
% 
%--------------------------------------------------------------------------
% Output:
%
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

% point P must be on the circle with center C and radius r, as well as A,
% B. We test if P is in the arc AB. Angle A is supposed to be inferior of
% angle B and they have the same sign.

PP=P-C; AA=A-C; BB=B-C;
[thPP,rPP]=cart2pol(PP(:,1),PP(:,2));
[thA,rA]=cart2pol(AA(:,1),AA(:,2));
[thB,rB]=cart2pol(BB(:,1),BB(:,2));

if thB < thA
    s=1;
    if sign(thPP) == sign(thB)
        thPP=thPP+2*pi;
    end
    thB=thB+2*pi;
end

% angle thB must be bigger than thA.
flag=(thPP >= thA) && (thPP <= thB);








%--------------------------------------------------------------------------
% new_pointindomain_clcl
%--------------------------------------------------------------------------

function flag=new_pointindomain_clcl(P,A,B,C,D,C1,r1,C2,r2)

%--------------------------------------------------------------------------
% Object:
% Let X be a domain of "clcl" type, i.e. whose closed boundary is defined, 
% counterclockwise and in this order, by an arc AB, a segment BC parallel 
% to the x axis, another arc CD and another segment DA parallel to the x 
% axis.
% 
% Suppose that P is a point in the plane. This routine checks if P is in
% the closed domain X, and in this case "flag=1", otherwise "flag=0".
%--------------------------------------------------------------------------
% Input: 
% P: n x 2 vector. The routine tests if P(:,k) is in the "clcl" domain X,
% and in this case flag(:,k)=1, otherwise flag(:,k)=0.
% A,B,C,D: points defining counterclockwise the "clcl" domain X.
% C1, r1: the arc AB is on the circle with center C1 and radius r1.
% C2, r2: the arc CD is on the circle with center C2 and radius r2.
%--------------------------------------------------------------------------
% Output:
% flag: the routine tests if the P(:,k) is in the "clcl" domain X,
%       and in this case flag(:,k)=1, otherwise flag(:,k)=0.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

XX=[A(1) B(1) C(1) D(1)];
YY=[A(2) B(2) C(2) D(2)];
a=min(XX); b=max(XX); c=min(YY); d=max(YY);

flag=[];
flag(:,1)=(P(:,1) >= a);
flag(:,2)=(P(:,1) <= b);
flag(:,3)=(P(:,2) >= c);
flag(:,4)=(P(:,2) <= d);

[th1,rr1]=mycart2pol(P-C1);

flag(:,5)=(rr1 >= r1);

[th2,rr2]=mycart2pol(P-C2);

flag(:,6)=(rr2 <= r2);
flag=floor(sum(flag,2)/size(flag,2));








%--------------------------------------------------------------------------
% mycart2pol
%--------------------------------------------------------------------------

function [th,r]=mycart2pol(PP)

%--------------------------------------------------------------------------
% Object:
% Variant of cart2pol, using PP as vector (and not two variables, one for
% coordinate). The angle is always in [0,2*pi).
%--------------------------------------------------------------------------
% Input: 
% PP: n x 2 vector.
%--------------------------------------------------------------------------
% Output:
% th,r: th(k),r(k) are, respectively, angles and radius of the generic
% point PP(:,k).
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

x=real(PP(:,1)); y=real(PP(:,2));
[th,r]=cart2pol(x,y);
if th < 0
    th=2*pi+th;
end








%--------------------------------------------------------------------------
% gqellblend
%--------------------------------------------------------------------------

function xyw = gqellblend(n,A,B,C,alpha,beta)

%--------------------------------------------------------------------------
% Object:
% The routine computes the nodes and weights of a product gaussian formula
% exact on total-degree bivariate polynomials of degree <=n
% on the planar region R obtained by linear blending (convex combination)
% of two trigonometric arcs with parametric equations
% P(theta)=A1*cos(theta)+B1*sin(theta)+C1
% Q(theta)=A2*cos(theta)+B2*sin(theta)+C2
% namely
% R = {(x,y)=t*P(theta)+(1-t)*Q(theta), t in [0,1], theta in [alpha,beta],
% 0<beta-alpha<=2*pi}
%--------------------------------------------------------------------------
% Input: 
% n: algebraic degree of exactness
% A,B,C: 2x2 matrices of the parametric arc coefficients:
% A1=A(1,:), B1=B(1,:), C1=C(1,:)
% A2=A(2,:), B2=B(2,:), C2=C(2,:)
% [alpha,beta]: angular interval, 0<beta-alpha<=2*pi
%--------------------------------------------------------------------------
% Output:
% xyw: 3 columns array of (xnodes,ynodes,weights)
%--------------------------------------------------------------------------
% External routines:
% r_jacobi.m, gauss.m from
% www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
%
% trigauss.m
% http://www.math.unipd.it/~marcov/mysoft/trigauss.m
%--------------------------------------------------------------------------
% Authors:
% Gaspare Da Fies, Alvise Sommariva and Marco Vianello
% University of Padova
% 2 June 2011
%--------------------------------------------------------------------------


% computing the algebraic and trigonometric degree increment
S1=abs((A(1,1)-A(2,1))*(B(1,2)-B(2,2))+(A(1,2)-A(2,2))*(B(2,1)-B(1,1)))>10*eps;
S2=abs((C(1,1)-C(2,1))*(B(1,2)-B(2,2))+(C(1,2)-C(2,2))*(B(2,1)-B(1,1)))>10*eps;
S3=abs((A(1,1)-A(2,1))*(C(1,2)-C(2,2))+(A(1,2)-A(2,2))*(C(2,1)-C(1,1)))>10*eps;

if (S1 || S2 || S3)
    h=1;
else
    h=0;
end

S4=abs(A(1,2)*A(2,1)-A(1,1)*A(2,2)-B(1,2)*B(2,1)+B(1,1)*B(2,2))>10*eps;
S5=abs(A(1,2)*B(2,1)-A(1,1)*B(2,2)+B(1,2)*A(2,1)-B(1,1)*A(2,2))>10*eps;
S6=abs(B(2,1)*(C(1,2)-C(2,2))-B(2,2)*(C(1,1)-C(2,1)))>10*eps;
S7=abs(A(2,1)*(C(1,2)-C(2,2))-A(2,2)*(C(1,1)-C(2,1)))>10*eps;
S8=abs((C(1,1)-C(2,1))*(B(1,2)-B(2,2))+(C(1,2)-C(2,2))*(B(2,1)-B(1,1)))>10*eps;
S9=abs((A(1,1)-A(2,1))*(C(1,2)-C(2,2))+(A(1,2)-A(2,2))*(C(2,1)-C(1,1)))>10*eps;

if (S4 || S5)
    k=2;
elseif (S6 || S7 || S8 || S9)
    k=1;
else
    k=0;
end

% trigonometric gaussian formula on the arc
tw=trigauss(n+k,alpha,beta);

% algebraic gaussian formula on [0,1]
ab=r_jacobi(ceil((n+h+1)/2),0,0);
xw=gauss(ceil((n+h+1)/2),ab);
xw(:,1)=xw(:,1)/2+1/2;
xw(:,2)=xw(:,2)/2;

% creating the grid
[t,theta]=meshgrid(xw(:,1),tw(:,1));
[w1,w2]=meshgrid(xw(:,2),tw(:,2));

% nodal cartesian coordinates and weights
s=sin(theta(:));
c=cos(theta(:));
p1=A(1,1)*c+B(1,1)*s+C(1,1);
p2=A(1,2)*c+B(1,2)*s+C(1,2);
q1=A(2,1)*c+B(2,1)*s+C(2,1);
q2=A(2,2)*c+B(2,2)*s+C(2,2);
dp1=-A(1,1)*s+B(1,1)*c;% plot(xyw(:,1)+trasl,xyw(:,2),'.','MarkerSize',6);
dp2=-A(1,2)*s+B(1,2)*c;
dq1=-A(2,1)*s+B(2,1)*c;
dq2=-A(2,2)*s+B(2,2)*c;

xyw(:,1)=p1.*t(:)+q1.*(1-t(:));
xyw(:,2)=p2.*t(:)+q2.*(1-t(:));
xyw(:,3)=abs((p1-q1).*(dp2.*t(:)+dq2.*(1-t(:))) - ...
    (p2-q2).*(dp1.*t(:)+dq1.*(1-t(:)))).* w1(:).*w2(:);








%--------------------------------------------------------------------------
% comprexcub
%--------------------------------------------------------------------------

function [pts,w,momerr] = comprexcub(deg,X,omega,pos)

%--------------------------------------------------------------------------
% Object:
% Compression of bivariate cubature formulas by Tchakaloff points
% or approximate Fekete points
% useful, for example, in node reduction of algebraic cubature formulas
% see the web page: http://www.math.unipd.it/~marcov/Tchakaloff.html
%--------------------------------------------------------------------------
% Input: 
% deg: polynomial exactness degree
% X: 2-column array of point coordinates
% omega: 1-column array of weights
% pos: NNLS for pos=1, QR with column pivoting for pos=0, Fast Lawson-Hanson
%      for pos=2.
%--------------------------------------------------------------------------
% Output:
% pts: 2-column array of extracted points
% w: 1-column array of corresponding weights (positive for pos=1)
% momerr: moment reconstruction error
%--------------------------------------------------------------------------
% Authors:
% Federico Piazzon, Alvise Sommariva and Marco Vianello
% University of Padova, May 2016
%--------------------------------------------------------------------------

% FUNCTION BODY
% total-degree Chebyshev-Vandermonde matrix at X


rect=[min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))];



V=chebvand(deg,X,rect);
% polynomial basis orthogonalization
[Q,R]=qr(V,0);
% tiny imaginary parts could appear
Q=real(Q);
% possible re-orthogonalization
% [Q,R]=qr(Q,0);

% moments of the orthogonal basis
orthmom=Q'*omega;
% weigths computation
switch pos 
    case 1
        % Tchakaloff points (positive weights)
        weights=lsqnonneg(Q',orthmom);
    case 2
        weights = lawsonhanson(Q',orthmom);
    otherwise
        % approximate Fekete points (possible negative weights)
        weights=Q'\orthmom;
end
% indexes of nonvanishing weights and compression
ind=find(abs(weights)>0);
pts=X(ind,:);
w=weights(ind);

% moment reconstruction error
% bivariate Chebyshev basis
mom=V'*omega;
momerr=norm(V(ind,:)'*w-mom);
% discrete OP basis
% momerr=norm(Q(ind,:)'*w-orthmom);








%--------------------------------------------------------------------------
% chebvand
%--------------------------------------------------------------------------

function V = chebvand(deg,X,rect)

%--------------------------------------------------------------------------
% Object:
% Vandermonde matrix at X points, relative to Chebyshev basis of degree deg
% on rectangle "rect". The basis is of tensorial type, but of total degree.
%--------------------------------------------------------------------------
% Input: 
% deg: polynomial degree.
% X: 2-column array of point coordinates.
% rect: 4-component vector such that the rectangle.
% [rect(1),rect(2)] x [rect(3),rect(4)] contains X.
%--------------------------------------------------------------------------
% Output:
% V: Chebyshev-Vandermonde matrix at X, graded lexic. order.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

% FUNCTION BODY
% rectangle containing the mesh
if isempty(rect)
    rect=[min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))];
end

% couples with length less or equal to deg
% graded lexicographical order
j=(0:1:deg);
[j1,j2]=meshgrid(j);
dim=(deg+1)*(deg+2)/2;
couples=zeros(dim,2);
for s=0:deg
    good=find(j1(:)+j2(:)==s);
    couples(1+s*(s+1)/2:(s+1)*(s+2)/2,:)=[j1(good) j2(good)];
end

% mapping the mesh in the square [-1,1]^2
a=rect(1);b=rect(2);c=rect(3);d=rect(4);
map=[(2*X(:,1)-b-a)/(b-a) (2*X(:,2)-d-c)/(d-c)];

% Chebyshev-Vandermonde matrix on the mesh
T1=chebpolys(deg,map(:,1));
T2=chebpolys(deg,map(:,2));
V=T1(:,couples(:,1)+1).*T2(:,couples(:,2)+1);








%--------------------------------------------------------------------------
% chebpolys
%--------------------------------------------------------------------------

function T=chebpolys(deg,x)

%--------------------------------------------------------------------------
% Object:
% This routine computes the Chebyshev-Vandermonde matrix on the real line
% by recurrence.
%--------------------------------------------------------------------------
% Input: 
% deg: maximum polynomial degree
% x: 1-column array of abscissas 
%--------------------------------------------------------------------------
% Output:
% T: Chebyshev-Vandermonde matrix at x, T(i,j+1)=T_j(x_i), j=0,...,deg.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

T=zeros(length(x),deg+1);
t0=ones(length(x),1);
T(:,1)=t0;
t1=x;
T(:,2)=t1;

for j=2:deg
    t2=2*x.*t1-t0;
    T(:,j+1)=t2;
    t0=t1;
    t1=t2;
end


