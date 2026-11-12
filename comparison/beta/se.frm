#: WorkSpace 200M
#: MaxTermSize 1M
#-
* Gluon (FIELD=glu) or ghost (FIELD=gho) self-energy at LOOPS loops, QCD, Landau gauge, Minkowski Feynman rules
* (Peskin-Schroeder), Wick-rotated to Euclidean integrals, then Forcer (MS-bar).
Off Statistics;
#include color.h
#include- forcer.h
Symbols g,a,D,nf,xi;
Dimension D;
Indices mu,nu,mu1,...,mu20,nu1,...,nu20;
Indices m1,...,m4;
Dimension NA;
Indices aa=NA,bb=NA,a1,...,a20=NA,b1,...,b20=NA,gc1,...,gc20=NA;
Indices cq1,...,cq20=NR;
Symbols sq1,...,sq20,sQ1,...,sQ20;
Vectors Qout,k1,...,k20;
Set kk:k1,...,k20;
CFunctions node,G,QQ,QB,GH,GB,DG,SQ,DGH,gam,vxtmp,dotE,dotP,sp,ff,fg;
Symbols n,cc,ss,ss1,...,ss4;
Indices ix1,...,ix4,c1,c2;
Indices e1,...,e20=NA;
Set ee:e1,...,e20;
CFunction cnt;
Symbol nn;
Vectors w1,...,w4,wk;
Model QCD;
  Particle qua,QUA,-2;
  Particle gho,GHO,-1;
  Particle glu,+3;
  Vertex QUA,qua,glu:g;
  Vertex GHO,gho,glu:g;
  Vertex glu,glu,glu:g;
  Vertex glu,glu,glu,glu:g^2;
EndModel;
* FORM's NOSNAIL_ option also removes non-zero diagrams at three loops (e.g. a four-gluon vertex with an external
* leg and two lines to one other vertex); only the scaleless ones are removed here, explicitly: a line from a vertex
* to itself, or both external legs on one vertex (no momentum through the rest)
L F = diagrams_(QCD,{`FIELD'},{`FIELD'},{Q,Qout},kk,`LOOPS',`ONEPI_'+`NOTADPOLE_');
.sort
#$n = termsin_(F);
#message `LOOPS' loops, `FIELD': `$n' diagrams
#ifdef `ONLY'
  if ( match(node_(?x,`ONLY'(?y),?z)) == 0 ) discard;
#endif
#ifdef `WITHOUT'
  if ( match(node_(?x,`WITHOUT'(?y),?z)) ) discard;
#endif
* external legs dropped; the vertex momenta kept for Forcer's topology (vx)
id node_(n?,1,?x) = 1;
id topo_(?x) = 1;
id node_(n?,cc?,?x) = node_(n,cc,?x)*vxtmp(?x);
repeat id vxtmp(?x,glu(?m),?y) = vxtmp(?x,?m,?y);
repeat id vxtmp(?x,qua(?m),?y) = vxtmp(?x,?m,?y);
repeat id vxtmp(?x,QUA(?m),?y) = vxtmp(?x,?m,?y);
repeat id vxtmp(?x,gho(?m),?y) = vxtmp(?x,?m,?y);
repeat id vxtmp(?x,GHO(?m),?y) = vxtmp(?x,?m,?y);
* the scaleless diagrams (they vanish): a line from a vertex to itself, or both external legs on one vertex
#do i = 1,20
  if ( match(vxtmp(?x,k`i',?y,-k`i',?z)) || match(vxtmp(?x,-k`i',?y,k`i',?z)) ) discard;
#enddo
if ( match(vxtmp(?x,Q,?y,Qout,?z)) || match(vxtmp(?x,Qout,?y,Q,?z)) ) discard;
.sort:nonzero;
#$nz = termsin_(F);
#message `LOOPS' loops, `FIELD': `$nz' diagrams not scaleless
#ifdef `COUNTONLY'
  .end
#endif
* the external fields
id node_(?x,glu(Q),?y) = node_(?x,G(1,Q,mu,aa),?y);
id node_(?x,glu(Qout),?y) = node_(?x,G(1,Qout,nu,bb),?y);
id node_(?x,gho(Q),?y) = node_(?x,GH(1,Q,aa),?y);
id node_(?x,GHO(Qout),?y) = node_(?x,GB(1,Qout,bb),?y);
id node_(?x,gho(Qout),?y) = node_(?x,GH(1,Qout,bb),?y);
id node_(?x,GHO(Q),?y) = node_(?x,GB(1,Q,aa),?y);
* the internal lines: indices per line end, one propagator per line
#do i = 1,20
  repeat id node_(?x,glu(k`i'),?y) = node_(?x,G(1,k`i',mu`i',a`i'),?y)*DG(k`i',mu`i',nu`i')*d_(a`i',b`i');
  repeat id node_(?x,glu(-k`i'),?y) = node_(?x,G(-1,k`i',nu`i',b`i'),?y);
  repeat id node_(?x,qua(k`i'),?y) = node_(?x,QQ(sq`i',cq`i'),?y)*SQ(sq`i',sQ`i',1,k`i');
  repeat id node_(?x,qua(-k`i'),?y) = node_(?x,QQ(sq`i',cq`i'),?y)*SQ(sq`i',sQ`i',-1,k`i');
  repeat id node_(?x,QUA(k`i'),?y) = node_(?x,QB(sQ`i',cq`i'),?y);
  repeat id node_(?x,QUA(-k`i'),?y) = node_(?x,QB(sQ`i',cq`i'),?y);
  repeat id node_(?x,gho(k`i'),?y) = node_(?x,GH(1,k`i',gc`i'),?y)*DGH(k`i');
  repeat id node_(?x,gho(-k`i'),?y) = node_(?x,GH(-1,k`i',gc`i'),?y)*DGH(k`i');
  repeat id node_(?x,GHO(k`i'),?y) = node_(?x,GB(1,k`i',gc`i'),?y);
  repeat id node_(?x,GHO(-k`i'),?y) = node_(?x,GB(-1,k`i',gc`i'),?y);
#enddo
* vertices (all momenta incoming; s the sign of the line momentum)
id node_(n?,cc?,QB(ss1?,c1?),QQ(ss2?,c2?),G(ss3?,w3?,m3?,ix3?)) = cc*i_*gam(ss1,ss2,m3)*T(c1,c2,ix3);
id node_(n?,cc?,GB(ss1?,w1?,c1?),GH(ss2?,w2?,c2?),G(ss3?,w3?,m3?,ix3?)) = cc*(-1)*f(c1,ix3,c2)*(-ss1)*w1(m3);
id node_(n?,cc?,G(ss1?,w1?,m1?,ix1?),G(ss2?,w2?,m2?,ix2?),G(ss3?,w3?,m3?,ix3?)) = cc*f(ix1,ix2,ix3)*(
     d_(m1,m2)*(ss1*w1(m3)-ss2*w2(m3)) + d_(m2,m3)*(ss2*w2(m1)-ss3*w3(m1)) + d_(m3,m1)*(ss3*w3(m2)-ss1*w1(m2)));
* the four-gluon vertex: one summed colour index per vertex, from a counter (two vertices in one term need two)
multiply cnt(1);
repeat id,once node_(n?,cc?,G(ss1?,w1?,m1?,ix1?),G(ss2?,w2?,m2?,ix2?),G(ss3?,w3?,m3?,ix3?),G(ss4?,w4?,m4?,ix4?))*cnt(nn?) = cc*(-i_)*(
     f(ix1,ix2,ee[nn])*f(ix3,ix4,ee[nn])*(d_(m1,m3)*d_(m2,m4)-d_(m1,m4)*d_(m2,m3))
   + f(ix1,ix3,ee[nn])*f(ix2,ix4,ee[nn])*(d_(m1,m2)*d_(m3,m4)-d_(m1,m4)*d_(m2,m3))
   + f(ix1,ix4,ee[nn])*f(ix2,ix3,ee[nn])*(d_(m1,m2)*d_(m3,m4)-d_(m1,m3)*d_(m2,m4)))*cnt(nn+1);
id cnt(nn?) = 1;
* propagators (Landau gauge)
id DG(wk?,m1?,m2?) = -i_*(d_(m1,m2)-(1-xi)*wk(m1)*wk(m2)/wk.wk)/wk.wk;
#ifdef `LANDAU'
  id xi = 0;
#endif
id SQ(ss1?,ss2?,ss?,wk?) = i_*ss*gam(ss1,ss2,wk)/wk.wk;
id DGH(wk?) = i_/wk.wk;
* the Dirac strings closed into traces
repeat id gam(ss1?,ss2?,?x)*gam(ss2?,ss3?,?y) = gam(ss1,ss3,?x,?y);
* each closed quark loop: one flavour sum, nf
id,once gam(ss1?,ss1?,?x) = g_(1,?x)*nf;
id,once gam(ss1?,ss1?,?x) = g_(2,?x)*nf;
id,once gam(ss1?,ss1?,?x) = g_(3,?x)*nf;
.sort:rules;
tracen,3;
tracen,2;
tracen,1;
.sort:traces;
* projection: gluon d^ab g_mu,nu / (NA (D-1) Q.Q); ghost d^ab / (NA Q.Q)
#if "`FIELD'" == "glu"
  #ifdef `LONG'
    multiply d_(aa,bb)*Q(mu)*Q(nu)/NA;
  #else
    multiply d_(aa,bb)*d_(mu,nu)/NA;
  #endif
#else
  multiply d_(aa,bb)/NA;
#endif
multiply replace_(Qout,-Q);
.sort:projected;
#call docolor
.sort:color;
* Wick rotation: every Minkowski scalar product -> minus the Euclidean one, i per loop measure
repeat id,once w1?.w2?^ss? = sign_(ss)*dotP(w1,w2,ss);
id dotP(w1?,w2?,ss?) = w1.w2^ss;
* Gamma/Q^2 = 1 + P: gluon P (D-1) = i X g_mu,nu / Q^2, ghost P = -i X / Q^2 (Q^2 Minkowski = -Q.Q Euclidean)
multiply i_^`LOOPS'/Q.Q*(-1);
#if "`FIELD'" == "glu"
  multiply i_;
#else
  multiply -i_;
#endif
id D = 4-2*ep;
id g = 1;
.sort:wick;
* the vertex functions for Forcer
id vxtmp(?x) = vx(?x);
* Forcer's momentum names
#do i = 1,11
  multiply replace_(k`i',p`i');
#enddo
.sort:forcerinput;
#ifdef `DUMP'
* the integrand for PQCF: scalar products as sp(a,b), bracketed by everything but the coefficient
  .sort:dump-0;
  repeat id,once w1?.w2?^ss? = sp(w1,w2,ss);
  Bracket vx,sp;
  .sort:dump-1;
  Format nospaces;
  Format 250;
  #write <`DUMP'> "%E", F
  .end
#endif
#call Forcer(msbarexpand=3)
B ep,NA,cA,cR,I2R,nf,xi;
Print +s;
.end
