function temp = dist(x,y,c1,c2)
%'''Distance(s) x,y away from a point c1,c2 in 2D'''
    tx=abs(c1-x);
    ty=abs(c2-y);
    temp = sqrt((tx).^2 + (ty).^2);
end