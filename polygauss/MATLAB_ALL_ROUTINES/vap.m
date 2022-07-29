
function [xyw,pts,w,momerr]=vap(n,y,r)

% quadrature by lune subtraction
% on vignetted annular pupils

% by Alvise Sommariva and Marco Vianello
% University of Padova, November 2 2017

% input
% n: degree of polynomial exactness
% y: circle center quotas
% r: circle radii

% there are five circles
% the circles are: 1 the main pupil, 2 & 3 the internal holes,
% 4 & 5 the upper circles

% WARNING: r(1),r(2),r(4)>0; put r(3)=0 or r(5)=0
% to eliminate circle 3 or circle 5

% WARNING: compatibility checks are still missing
% circles 2 & 3 are inside circle 1
% circles 4 & 5 intersect circle 1
% circles 4 & 5 do not intersect circles 2 & 3

% annulus, as a lune: diff(disk1,disk2)
xyw=gqlune(n,0,y(1),r(1),0,y(2),r(2));

% lune subtraction: diff(disk3,disk2)
if r(3)>0
    nw=gqlune(n,0,y(3),r(3),0,y(2),r(2));
    xyw=[xyw;[nw(:,1:2) -nw(:,3)]];
end

% upper circles: checking which has the lower
% intersections with the main circle
if r(5)>0
    [xout14,yout14] = circcirc(0,y(1),r(1),0,y(4),r(4));
    [xout15,yout15] = circcirc(0,y(1),r(1),0,y(5),r(5));
    if yout14(1)<yout15(1)
        a=4;b=5;
    else
        a=5;b=4;
    end
   
else
    a=4;
end

% lune subtraction with the lower: diff(disk1,diska)
nw=gqlune(n,0,y(1),r(1),0,y(a),r(a));
xyw=[xyw;[nw(:,1:2) -nw(:,3)]];

% checking whther another lunette has to be subtracted:
% if the intersection of the upper circles 4 & 5
% are inside the main circle the lunette diff(diska,diskb)
% is active, otherwise it is not
if r(5)>0
    [xout,yout] = circcirc(0,y(4),r(4),0,y(5),r(5));
    if (xout(1)^2+(yout(1)-y(1))^2<r(1)^2)
        nw=gqlune(n,0,y(a),r(a),0,y(b),r(b));
        xyw=[xyw;[nw(:,1:2) -nw(:,3)]];
    end
end

% compression of the cubature formula
[pts,w,momerr] = comprexcub(n,xyw(:,1:2),xyw(:,3),0);

end





