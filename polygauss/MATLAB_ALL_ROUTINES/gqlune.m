function xyw = gqlune(n,x1,y1,r1,x2,y2,r2)

% by Gaspare Da Fies and Marco Vianello, University of Padova
% March 1, 2013 

% computes nodes and weights of a product Gaussian formula of 
% polynomial exactness degree n on the circular lune obtained as 
% difference of the disk with center (x1,y1) and radius r1, and the disk 
% with center (x2,y2) and radius r2

% the formula is chosen to have the lowest possible cardinality among 
% the three formulas studied in: G. Da Fies and M. Vianello, Product 
% Gaussian quadrature on circular lunes  
% preprint online at: http://www.math.unipd.it/~marcov/pdf/lune.pdf

% uses the routines:
%
% r_jacobi.m, gauss.m from
% www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
%
% trigauss.m
% http://www.math.unipd.it/~marcov/mysoft/trigauss.m


% input:
% n: degree of polynomial exactness
% x1,y1: center coordinates of the first disk
% r1: radius of the first disk
% x2,y2: center coordinates of the second disk
% r2: radius of the second disk

% output:
% xyw: (n+1) x 3 array of (nodes,weights)


% distance between the centers
d=sqrt((x1-x2)^2+(y1-y2)^2);

% by translation and rotation we can always reduce the lune to the case
% disk(0,0,r1) setminus disk(-d,0,r2)

if d>=(r1+r2)
    
    % zero measure intersection
    % the difference is the first disk
    
    % trigonometric gaussian formula on [-pi,pi]
    tw=[pi*(-n-2:2:n+2)'/(n+3), repmat(2*pi/(n+3),n+3,1)];
    
    % algebraic gaussian formula on [-1,1]
    ab=r_jacobi(ceil((n+1)/2),0,0);
    xw=gauss(ceil((n+1)/2),ab);
    
    % creating the grid
    [t,theta]=meshgrid(xw(:,1),tw((1:ceil((n+2)/2)),1));
    [w1,w2]=meshgrid(xw(:,2),tw((1:ceil((n+2)/2)),2));
    
    % nodal cartesian coordinates and weights
    c=cos(theta(:));
    s=sin(theta(:));
    xyw(:,1)=r1*c;
    xyw(:,2)=r1*t(:).*s;
    xyw(:,3)=r1^2*s.^2.*w1(:).*w2(:);
    
elseif d<=abs(r1-r2)
    
    % one disk is included in the other
    
    if(r1<=r2)
        
        % zero measure difference
        
        % returns the first center with null weight
        xyw=[0 0 0];
        
    else
        
        % the difference is a generalized annulus
        
        % trigonometric gaussian formula on [-pi,pi]
        tw=[pi*(-n-1:2:n+1)'/(n+2), repmat(2*pi/(n+2),n+2,1)];
        
        % algebraic gaussian formula on [0,1] 
        ab=r_jacobi(ceil((n+2)/2),0,0);
        xw=gauss(ceil((n+2)/2),ab);
        xw(:,1)=xw(:,1)/2+1/2;
        xw(:,2)=xw(:,2)/2;
        
        % creating the grid
        [t,theta]=meshgrid(xw(:,1),tw(:,1));
        [w1,w2]=meshgrid(xw(:,2),tw(:,2));
        
        % nodal cartesian coordinates and weights
        c=cos(theta(:));
        s=sin(theta(:));
        xyw(:,1)=(r1*c).*t(:)+(r2*c-d).*(1-t(:));
        xyw(:,2)=(r1*s).*t(:)+(r2*s).*(1-t(:));
        xyw(:,3)=(r1*t(:)+r2*(1-t(:))).*(d*c-r2+r1).*w1(:).*w2(:);
        
    end
    
else
    
    % the difference is a proper lune
    
    % angles of the two arcs
    omega2=pi-acos((d^2+r1^2-r2^2)/(2*d*r1));
    omega1=acos((d^2+r2^2-r1^2)/(2*d*r2));
    
    if omega1<=atan(2*(1-cos(omega2))/sin(omega2))
        
        % trigonometric gaussian formula on the arc
        tw1=trigauss(n+2,-omega2,omega2);
        tw2=trigauss(n+2,-omega1,omega1);
        
        % creating the grid
        [theta1,theta2]=meshgrid(tw1((1:ceil((n+2)/2)),1),tw2(:,1));
        [w1,w2]=meshgrid(tw1((1:ceil((n+2)/2)),2),tw2(:,2));
        
        % nodal cartesian coordinates and weights
        c1=cos(theta1(:)); 
        s1=sin(theta1(:));
        c2=cos(theta2(:)); 
        s2=sin(theta2(:));
        
        xyw(:,1)=r1*(c1 + (1-c1)/(1-cos(omega2))*sin(omega2)/sin(omega1).*(c2-cos(omega1)));
        xyw(:,2)=r1/sin(omega1)*s1.*s2;
        xyw(:,3)=-(r1/sin(omega1))^2*sin(omega2)/(1-cos(omega2))...
        *( (-sin(omega1)/sin(omega2)*(1-cos(omega2))-cos(omega1))...
        *s1.^2.*c2 + s1.^2.*c2.^2 + c1.*s2.^2 - c1.^2.*s2.^2 ).*w1(:).*w2(:);
                 
    elseif (cos(omega2)+cos(omega2-omega1))^2<=4*cos(omega1)
        
        % trigonometric gaussian formula on the arcs
        tw1=trigauss(n+2,-omega2,omega2);
        tw2=trigauss(n+2,-omega1,omega1);
        
        % creating the grid
        [theta1,theta2]=meshgrid(tw1(:,1),tw2((1:ceil((n+2)/2)),1));
        [w1,w2]=meshgrid(tw1(:,2),tw2((1:ceil((n+2)/2)),2));
        
        % nodal cartesian coordinates and weights
        c1=cos(theta1(:)); 
        s1=sin(theta1(:));
        c2=cos(theta2(:)); 
        s2=sin(theta2(:));
        
        xyw(:,1)=r1*((c2-cos(omega1))*(cos(omega2)/(1-cos(omega1)) ...
            +sin(omega2)/sin(omega1)) + c1.*(1-c2)/(1-cos(omega1)));
        xyw(:,2)=r1/sin(omega1)*s1.*s2;
        xyw(:,3)=r1^2/sin(omega1)/(1-cos(omega1))*((-sin(omega2)/sin(omega1)...
            *(1-cos(omega1))-cos(omega2))*s2.^2.*c1 + s2.^2.*c1.^2 ...
            + c2.*s1.^2 - c2.^2.*s1.^2 ) .*w1(:).*w2(:);
        
    else
        
        % trigonometric gaussian formula on the arcs
        tw1=trigauss(n+2,omega1,omega2);
        tw2=trigauss(n+1,-omega1,omega1);
        
        
        % creating the grid     
        [theta1,theta2]=meshgrid(tw1(:,1),tw2(:,1));
        [w1,w2]=meshgrid(tw1(:,2),tw2(:,2));
        
        % nodal cartesian coordinates and weights
        c1=cos(theta1(:));
        s1=sin(theta1(:));
        c2=cos(theta2(:)); 
        s2=sin(theta2(:));
        xyw(:,1)=r1/sin(omega1)*(sin(omega1)*c1-cos(omega1)*s1+s1.*c2);
        xyw(:,2)=r1/sin(omega1)*(s1.*s2);
        xyw(:,3)=(r1/sin(omega1))^2*((cos(omega1)*c1+sin(omega1)*s1).*c2-c1).*s1.*w1(:).*w2(:);
               
    end
    
end

% rotating and translating the nodes to go back to the original lune
% disk(x1,y1,r1) setminus disk(x2,y2,r2)

% rotation
if(d~=0)
    
    R=[(x1-x2)/d (y1-y2)/d; (y1-y2)/d -(x1-x2)/d];
    xyw(:,(1:2))=xyw(:,(1:2))*R';
    
end

% translation
xyw(:,1)=x1+xyw(:,1);
xyw(:,2)=y1+xyw(:,2);

% test plot of the quadrature nodes
%figure;
%hold on;
%plot(x1+r1*cos(linspace(-pi,pi,500)),y1+r1*sin(linspace(-pi,pi,500)),'b');
%plot(x2+r2*cos(linspace(-pi,pi,500)),y2+r2*sin(linspace(-pi,pi,500)),'b');
%plot(xyw(:,1),xyw(:,2),'or','MarkerSize',4);
%axis equal;


end


