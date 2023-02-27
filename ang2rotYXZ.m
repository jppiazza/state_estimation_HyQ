function [R] = ang2rotYXZ(x, y, z)

g = deg2rad(x); 
b = deg2rad(y);
a = deg2rad(z);

cg = cos(g); sg = sin(g);
cb = cos(b); sb = sin(b);
ca = cos(a); sa = sin(a);

R = [sa*sb*sg+ca*cg sa*sb*cg-ca*sg sa*cb;
     cb*sg          cb*cg          -sb;
     ca*sb*sg-sa*cg          ca*sb*cg+sa*sg      ca*cb];

end

