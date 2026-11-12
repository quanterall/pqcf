#-
#include- forcer.h
L M1 = vx(-Q,p1,p2)*vx(Q,-p1,-p2)*p1.p1^-1*p2.p2^-1;
L M2 = vx(-p1,-p2,Q)*vx(p2,p1,-Q)*p1.p1^-1*p2.p2^-1;
#call Forcer
B ep;
P +s;
.end
