#-
#include- forcer.h
L M1 = vx(-Q,-p3,-p4)*vx(p1,p2,p3,p4)*vx(Q,-p1,-p2)*p1.p1^-1*p2.p2^-1*p3.p3^-1*p4.p4^-1;
L M2 = vx(-Q,p1,-p3)*vx(p2,-p4,-p5)*vx(p3,p4,p5)*vx(Q,-p1,-p2)*p1.p1^-1*p4.p4^-1*p5.p5^-1;
L M3 = vx(-Q,p1,p2,p3)*vx(Q,-p1,-p2,-p3)*p1.p1^-1*p2.p2^-1*p3.p3^-1;
L M4 = vx(-Q,p2,p3,p4)*vx(p1,-p3,-p4)*vx(Q,-p1,-p2)*p2.p2^-1*p3.p3^-1*p4.p4^-1;
L M5 = vx(-Q,p4,p5)*vx(p1,-p3,-p4)*vx(p2,p3,-p5)*vx(Q,-p1,-p2)*p1.p1^-1*p3.p3^-1*p5.p5^-1;
L M6 = vx(-Q,p4,p5)*vx(p1,-p3,-p4)*vx(p2,p3,-p5)*vx(Q,-p1,-p2)*p1.p1^-1*p2.p2^-1*p4.p4^-1*p5.p5^-1;
L M7 = vx(-p1,-p2,-p3)*vx(-Q,p1,-p4)*vx(Q,p2,p3,p4)*p2.p2^-1*p3.p3^-1*p4.p4^-1;
L M8 = vx(-p1,-p2,Q)*vx(-p3,p1,-Q)*vx(p2,-p4,-p5)*vx(p4,p3,p5)*p1.p1^-1*p4.p4^-1*p5.p5^-1;
L M9 = vx(-p2,-p1,Q)*vx(-Q,p4,p5)*vx(p1,-p3,-p4)*vx(p3,p2,-p5)*p1.p1^-1*p3.p3^-1*p5.p5^-1;
L M10 = vx(-p2,-p1,Q)*vx(-Q,p4,p5)*vx(p1,-p3,-p4)*vx(p3,p2,-p5)*p1.p1^-1*p2.p2^-1*p4.p4^-1*p5.p5^-1;
L M11 = vx(-p2,-p1,Q)*vx(-p4,p2,-p5)*vx(p1,-p3,-Q)*vx(p3,p4,p5)*p1.p1^-1*p4.p4^-1*p5.p5^-1;
L M12 = vx(-p2,-p1,Q)*vx(-p5,p2,p3)*vx(p1,-p4,-p3)*vx(p4,p5,-Q)*p1.p1^-1*p3.p3^-1*p5.p5^-1;
L M13 = vx(-p2,-p1,Q)*vx(-p5,p2,p3)*vx(p1,-p4,-p3)*vx(p4,p5,-Q)*p1.p1^-1*p2.p2^-1*p4.p4^-1*p5.p5^-1;
L M14 = vx(-p4,-p3,p1)*vx(p3,-p5,p2)*vx(p5,p4,-Q)*vx(Q,-p1,-p2)*p1.p1^-1*p3.p3^-1*p5.p5^-1;
L M15 = vx(-p4,-p3,p1)*vx(p3,-p5,p2)*vx(p5,p4,-Q)*vx(Q,-p1,-p2)*p1.p1^-1*p2.p2^-1*p4.p4^-1*p5.p5^-1;
L M16 = vx(-p4,-p5,p2)*vx(-Q,p1,-p3)*vx(p5,p4,p3)*vx(Q,-p1,-p2)*p1.p1^-1*p4.p4^-1*p5.p5^-1;
#call Forcer
B ep;
P +s;
.end
