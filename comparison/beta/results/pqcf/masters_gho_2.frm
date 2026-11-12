#-
#include- forcer.h
L M1 = vx(-p1,Q,-p2)*vx(-Q,p1,-p3)*vx(p2,-p4,-p5)*vx(p3,p4,p5)*p1.p1^-1*p4.p4^-1*p5.p5^-1;
L M2 = vx(-p1,Q,-p2)*vx(-p4,-p5,p2)*vx(-Q,p1,-p3)*vx(p5,p4,p3)*p1.p1^-1*p4.p4^-1*p5.p5^-1;
L M3 = vx(-p2,Q,-p1)*vx(-p4,-p3,p1)*vx(-Q,p4,p5)*vx(p3,p2,-p5)*p1.p1^-1*p3.p3^-1*p5.p5^-1;
L M4 = vx(-p2,Q,-p1)*vx(-p4,-p3,p1)*vx(-Q,p4,p5)*vx(p3,p2,-p5)*p1.p1^-1*p2.p2^-1*p4.p4^-1*p5.p5^-1;
L M5 = vx(-p2,Q,-p1)*vx(-p4,p2,-p5)*vx(-Q,-p3,p1)*vx(p3,p4,p5)*p1.p1^-1*p4.p4^-1*p5.p5^-1;
L M6 = vx(-p2,Q,-p1)*vx(-p5,p2,p3)*vx(-Q,p5,p4)*vx(p1,-p3,-p4)*p1.p1^-1*p3.p3^-1*p5.p5^-1;
L M7 = vx(-p2,Q,-p1)*vx(-p5,p2,p3)*vx(-Q,p5,p4)*vx(p1,-p3,-p4)*p1.p1^-1*p2.p2^-1*p4.p4^-1*p5.p5^-1;
#call Forcer
B ep;
P +s;
.end
