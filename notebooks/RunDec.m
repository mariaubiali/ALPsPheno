(* ::Package:: *)

(*
 RunDec: a Mathematica package for running and decoupling of the
 strong coupling and quark masses

 Copyright: K.G. Chetyrkin, J.H. K"uhn and M. Steinhauser (Jan. 2000)
*)

(*
18Jul06: mOS2mSI[]: cut[api,loops] -> cut[apiM,loops]
23May07: mMS2mRGImod[] introduced
06Mar08: Rationalize[] introduced in AlphasExact[]
28Oct11: AsmMSrunexact[] included
28Oct11: mMS2mRGI[]: bug for nl=0 fixed
*)

(* ************************************************************ *)

(*
 The examples shown in the paper can be found in the the file
 RunDec_ex.m
*)

(* ************************************************************ *)

Print["RunDec: a Mathematica package for running and decoupling of the"];
Print["        strong coupling and quark masses"];
Print["by K.G. Chetyrkin, J.H. K\\\"uhn and M. Steinhauser (January 2000)"];

(* ************************************************************ *)
(* ************************************************************ *)

(*
 default values of the numerical constants
*)

NumDef = {
    asMz -> 0.118,
    Mz -> 91.18,
    Mt -> 172.5,
    Mb -> 4.75,
    Mc -> 1.44,
    muc -> 1.2,
    mub -> 3.97,
    Mtau -> 1.777
    };

(* ************************************************************ *)
(* ************************************************************ *)

BeginPackage["RunDec`"];

 LamExpl::usage =
    "LamExpl[als, mu, nf, l] computes \\Lambda^(nf) with l-loop accuracy
using the explicite formulae \\ln(\\mu/\\Lambda) = f(als,mu,nf).";

 LamImpl::usage =
    "LamImpl[als, mu, nf, l] computes \\Lambda^(nf) with l-loop accuracy
using the implicite formulae \\alpha_s = f(mu,\\Lambda,nf).";

 AlphasLam::usage = 
    "AlphasLam[lam ,mu ,nf ,l] computes \\alpha_s^(nf)(mu) to l-loop accuracy
using the formulae \\alpha_s = f(mu,\\Lambda,nf).";

 AlphasExact::usage = 
    "AlphasExact[als0 ,mu0 ,mu ,nf ,l] computes \\alpha_s^(nf)(mu) integrating
the \beta function numerically. The initial condition \\alpha_s(mu0)=als0
is used.";

 AsmMSrunexact::usage = 
    "AsmMSrunexact[mmu0, asmu0, mu0, mu, nf, l] computes \\alpha_s^(nf)(mu) 
and m_q^(nf)(mu) solving the coupled system of differential equations 
simultaneously. The initial conditions are \\alpha_s(mu0)=als0 
and m_q(mu)=mmu0.";

 mOS2mMS::usage =
    "mOS2mMS[mOS, {mq}, asmu, mu, nf, l] computes the MS-bar mass at the
scale mu. {mq} represents a set of light quark masses in the on-shell 
scheme leading to corrections of O(as^2 mq/mOS).
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.
mOS2mMS[mOS, nf, l] computes the MS-bar mass at the the scale mOS. 
\\alpha_s^(nf)(M) is computed from \\alpha_s^(5)(Mz) as defined in NumDef.";

 mMS2mOS::usage =
    "mMS2mOS[mMS, {mq}, asmu, mu, nf, l] computes the on-shell mass.
{mq} represents a set of light quark masses in the on-shell scheme
leading to corrections of O(as^2 mq/mMS).
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.
mMS2mOS[mum, nf, l] computes the on-shell mass where mum is 
defined through mum=m^(nf)(mum).
\\alpha_s^(nf)(mum) is computed from \\alpha_s^(5)(Mz) as defined in NumDef.";

 mOS2mMSrun::usage =
    "mOS2mMSrun[mOS, {mq}, asmu, mu, nf, l] computes the MS-bar mass at the
scale mu.
In contrast to mOS2mMS[..] in a first step mMS(mMS) is computed
and only then the conversion to mOS is performed.
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mMS2mOSrun::usage =
    "mMS2mOSrun[mMS, {mq}, asmu, mu, nf, l] computes the on-shell mass.
In contrast to mMS2mOS[..] in a first step mMS(mMS) is computed
and afterwards mMS(mu) is evaluated.
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mOS2mMSit::usage =
    "mOS2mMSit[mOS, {mq}, asmu, mu, nf, l] computes the MS-bar mass at the
scale mu. However, in contrast to mOS2mMS[..], the relation
M_OS = m_MS * (1 + ...) is solved iteratively. This has the advantage
that the light quark masses can be evaluated in the MS-bar scheme.";

 mOS2mSI::usage =
    "mOS2mSI[mOS, {mq}, asM, nf, l] computes the mass
\\mu_m = m_MS(\\mu_m). {mq} represents a set of light quark masses in 
the on-shell scheme leading to corrections of O(as^2 mq/mOS).
asM = \\alpha_s^(nf)(M_OS), nf is the number
of active flavours and l represents the number of loops.";

 mMS2mMS::usage =
    "mMS2mMS[mmu0, asmu0, asmu, nf, l] computes m(mu) from the knowledge
of m(mu0). asmu = \\alpha_s(\\mu), asmu0 = \\alpha_s(\\mu0), nf is the number 
of active flavours and l represents the number of loops.";

 mMS2mSI::usage =
    "mMS2mSI[mMS, asmu, mu, nf, l] computes the scale
invarant mass mMS(mMS). mMS is the MS-bar mass at the scale mu,
asmu = \\alphaA_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mMS2mRI::usage =
    "mMS2mRI[mMS, asmu, nf, l] computes the regularization invariant
mass. mMS is the MS-bar mass,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mRI2mMS::usage =
    "mRI2mMS[mRI, asmu, nf, l] computes the MS-bar quark mass.
mRI is the regularization invariant mass,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mMS2mRGI::usage =
    "mMS2mRGI[mMS, asmu, nf, l] computes the renormalization group
invariant mass. mMS is the MS-bar mass,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 mMS2mRGImod::usage =
    "mMS2mRGImod[mMS, asmu, nf, l] computes the renormalization group
invariant mass. mMS is the MS-bar mass,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.
(slightly modified version of mMS2mRGI[], factor '2 beta_0')";

 mRGI2mMS::usage =
    "mRGI2mMS[mRGI, asmu, nf, l] computes the MS-bar quark mass.
mRGI is the renormalization group invariant mass,
asmu = \\alpha_s^(nf)(\\mu), nf is the number
of active flavours and l represents the number of loops.";

 DecAsUpOS::usage = 
    "DecAsUpOS[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl+1)(muth)
from the knowledge of als=\\alpha_s^(nl)(muth). massth is the on-shell
mass value of the heavy quark. l specifies the number of loops.
It coincides with the parameter used for the running, i.e. for l=1
one has \\alpha_s^(nl+1)(muth)=\\alpha_s^(nl)(muth), for l=2 the
one-loop decoupling formula is used, ...";

 DecAsDownOS::usage = 
    "DecAsDownOS[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl)(muth)
from the knowledge of als=\\alpha_s^(nl+1)(muth). For the other parameters see
DecAsUpOS[als ,massth ,muth ,nl ,l ]";

 DecAsUpMS::usage = 
    "DecAsUpMS[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl+1)(muth)
from the knowledge of als=\\alpha_s^(nl)(muth). massth is the MS-bar
mass value of the heavy quark. l specifies the number of loops.
It coincides with the parameter used for the running, i.e. for l=1
one has \\alpha_s^(nl+1)(muth)=\\alpha_s^(nl)(muth), for l=2 the
one-loop decoupling formula is used, ...";

 DecAsDownMS::usage = 
    "DecAsDownMS[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl)(muth)
from the knowledge of als=\\alpha_s^(nl+1)(muth). For the other parameters see
DecAsUpMS[als ,massth ,muth ,nl ,l ]";

 DecAsUpSI::usage = 
    "DecAsUpSI[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl+1)(muth)
from the knowledge of als=\\alpha_s^(nl)(muth). massth is the
scale invariant
mass value of the heavy quark. l specifies the number of loops.
It coincides with the parameter used for the running, i.e. for l=1
one has \\alpha_s^(nl+1)(muth)=\\alpha_s^(nl)(muth), for l=2 the
one-loop decoupling formula is used, ...";

 DecAsDownSI::usage = 
    "DecAsDownSI[als ,massth ,muth ,nl ,l ] computes \\alpha_s^(nl)(muth)
from the knowledge of als=\\alpha_s^(nl+1)(muth). For the other parameters see
DecAsUpSI[als ,massth ,muth ,nl ,l ]";

 DecMqUpOS::usage = 
    "DecMqUpOS[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl+1)(muth)
from the knowledge of mq=m_q^(nl)(muth). als is \\alpha_s^(nl)(muth).
massth is the on-shell mass value
of the heavy quark. l specifies the number of loops. It coincides with the
parameter used for the running, i.e. for l=1 one has 
m_q^(nl+1)(muth)=m_q^(nl)(muth), for l=2 the one-loop 
decoupling formula is used, ...";

 DecMqDownOS::usage = 
    "DecMqDownOS[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl)(muth)
from the knowledge of mq=m_q^(nl+1)(muth). als is \\alpha_s^(nl+1)(muth).
For the other parameters see DecMqUpOS[mq, als ,massth ,muth ,nl ,l ]";

 DecMqUpMS::usage = 
    "DecMqUpMS[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl+1)(muth)
from the knowledge of mq=m_q^(nl)(muth). als is \\alpha_s^(nl)(muth).
massth is the MS-bar mass value
of the heavy quark. l specifies the number of loops. It coincides with the
parameter used for the running, i.e. for l=1 one has 
m_q^(nl+1)(muth)=m_q^(nl)(muth), for l=2 the one-loop 
decoupling formula is used, ...";

 DecMqDownMS::usage = 
    "DecMqDownMS[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl)(muth)
from the knowledge of mq=m_q^(nl+1)(muth). als is \\alpha_s^(nl+1)(muth).
For the other parameters see DecMqUpMS[mq, als ,massth ,muth ,nl ,l ]";

 DecMqUpSI::usage = 
    "DecMqUpSI[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl+1)(muth)
from the knowledge of mq=m_q^(nl)(muth). als is \\alpha_s^(nl)(muth).
massth is the scale invariant mass value
of the heavy quark. l specifies the number of loops. It coincides with the
parameter used for the running, i.e. for l=1 one has 
m_q^(nl+1)(muth)=m_q^(nl)(muth), for l=2 the one-loop 
decoupling formula is used, ...";

 DecMqDownSI::usage = 
    "DecMqDownSI[mq, als ,massth ,muth ,nl ,l ] computes m_q^(nl)(muth)
from the knowledge of mq=m_q^(nl+1)(muth). als is \\alpha_s^(nl+1)(muth).
For the other parameters see DecMqUpSI[mq, als ,massth ,muth ,nl ,l ]";

 DecLambdaUp::usage = 
    "DecLambdaUp[lam, massth, nl ,l ] computes \\Lambda^(nl+1) 
from the knowledge of lam=\\Lambda^(nl).
massth is the scale invariant mass value
of the heavy quark.
l specifies the number of loops. It coincides with the
parameter used for the running.";

 DecLambdaDown::usage = 
    "DecLambdaDown[lam, massth, nl ,l ] computes \\Lambda^(nl) 
from the knowledge of lam=\\Lambda^(nl+1).
massth is the scale invariant mass value
of the heavy quark.
l specifies the number of loops. It coincides with the
parameter used for the running.";

 AlL2AlH::usage = 
    "AlL2Alh[als, mu1, decpar, mu2, l] computes alphas(mu2) with h active
flavours from the knowledge of als=alphas(mu1) with l active flavours.
h > l must be fulfilled. decpar is a set which 
may contain several triples  indicating the
number of flavours, the heavy (on-shell)
quark mass and the scale at which the
decoupling is performed. The running is performed with
AlphasExact[..] and the decoupling with DecAsUpOS[..].
    ";

 AlH2AlL::usage = 
    "AlH2AlL[als, mu1, decpar, mu2, l] computes alphas(mu2) with l active
flavours from the knowledge of als=alphas(mu1) with h active flavours.
h > l must be fulfilled. decpar is a set which 
may contain several triples  indicating the
number of flavours, the heavy (on-shell)
quark mass and the scale at which the
decoupling is performed. The running is performed with
AlphasExact[..] and the decoupling with DecAsDownOS[..].
    ";

 mL2mH::usage = 
    "mL2mH[mql, asl, mu1, decpar, mu2, l] computes mq(mu2) with h active
flavours from the knowledge of mql=mq(mu1) with l active flavours.
asl=\\alpha_s(mu1) with l active flavours.
h > l must be fulfilled. decpar is a set which 
may contain several triples  indicating the
number of flavours, the heavy (on-shell)
quark mass and the scale at which the
decoupling is performed. The running of the coupling constant and 
the quark mass is performed with
AlphasExact[..] and mMS2mMS[..], respectively, and 
the decoupling with DecMqUpOS[..] and DecAsUpOS[..].
    ";

 mH2mL::usage = 
    "mH2L[mqh ,ash ,mu1 , decpar, mu2, l] computes mq(mu2) with l active
flavours from the knowledge of mqh=mq(mu1) with h active flavours.
ash=\\alpha_s(mu1) with h active flavours.
h > l must be fulfilled. decpar is a set which 
may contain several triples  indicating the
number of flavours, the heavy (on-shell)
quark mass and the scale at which the
decoupling is performed. The running of the coupling constant and 
the quark mass is performed with
AlphasExact[..] and mMS2mMS[..], respectively, and 
the decoupling with DecMqDownOS[..] and DecAsDownOS[..].
    ";

 AsRunDec::usage = 
    "AsRunDec[als ,mu0 ,mu , l] computes \\alpha_s(mu) from the knowledge of
als=\\alpha_s(mu0) The number of active flavours is determined automatically
from the choices of mu and mu0, respectively. The decoupling is performed at 
the pole mass of the respective heavy quark.";

 Mc5Mzfrommuc4::usage = 
    "";

 AsMbfromAsMz::usage = 
    "";

 As4McfromAs5Mz::usage = 
    "";

(* ************************************************************ *)

Begin["`Modules`"];

cut[x_, n_Integer] := (x^p_Integer /; p > n -> 0);

(*
 argument used in N[] 
*)

$NumPrec=20;
numprec=$NumPrec;

(*
 minimal allowed ratio of \mu/\Lambda for which no warning is printed
*)
rmulam = 1.5;

(*
 numerical expression of some symbols ... 
*)

num1 = {
    z2 -> Zeta[2],
    z3 -> Zeta[3],
    z4 -> Zeta[4],
    z5 -> Zeta[5],
    log2 -> Log[2],
    B4->16*PolyLog[4,1/2]+2/3*Log[2]^4-2/3*Pi^2*Log[2]^2-13/180*Pi^4
};

(*
 numerical expression for colour factors (QCD)
*)

cf2num = {
    cf -> 4/3,
    ca -> 3,
    tr -> 1/2
};

(*
 coefficients of \beta function 
*)

setbeta = {
    b0 -> 11/4 - nf/6, 
    b1 -> 51/8 - (19*nf)/24,
    b2 -> 2857/128 - (5033*nf)/1152 + (325*nf^2)/3456,
    b3 -> 149753/1536 - (1078361*nf)/41472 + (50065*nf^2)/41472 +
	(1093*nf^3)/186624 + (891*Zeta[3])/64 - (1627*nf*Zeta[3])/1728 +
	    (809*nf^2*Zeta[3])/2592
	    };

(*
 coefficients of \gamma_m function 
*)

setgamma = {
    g0 -> 1,
    g1 -> (101/2 - (5*nf)/3)/12,
    g2 -> (3747/4 - (nf*(1108/9 + (70*nf)/27 + 80*Zeta[3]))/2)/48,
    g3 -> 4603055/41472 - (91723*nf)/6912 + (2621*nf^2)/31104 - 
	(83*nf^3)/15552 +
     (11*nf*Pi^4)/288 - (nf^2*Pi^4)/432 + (530*Zeta[3])/27 -
	 (2137*nf*Zeta[3])/144 + (25*nf^2*Zeta[3])/72 + (nf^3*Zeta[3])/108 -
	     (275*Zeta[5])/8 + (575*nf*Zeta[5])/72
};

(*
 \alpha_s as function of \mu and \Lambda expanded in 1/\Lambda\beta_0
*)

setasL = {
    as0 -> 1/L/b0,
    as1 -> 1/L/b0 - 1/(b0*L)^2*b1/b0*Log[L],
    as2 -> ( 1/L/b0 - 1/(b0*L)^2*b1/b0*Log[L] +
	    1/(b0*L)^3*(b1^2/b0^2*(Log[L]^2-Log[L]-1)+b2/b0)
	    ),
    as3 -> (
	    1/(b0*L) - (b1*Log[L])/(b0^3*L^2) + 
	    (b2/b0^4 + (b1^2*(-1 - Log[L] + Log[L]^2))/b0^5)/L^3 + 
	    (b3/(2*b0^5) - (3*b1*b2*Log[L])/b0^6 + 
	     (b1^3*(-1/2 + 2*Log[L] + (5*Log[L]^2)/2 - Log[L]^3))/b0^7)/L^4
	    )
    };


(*
 \zeta_g, 1/\zeta_g, ...
 as6to5:  api^(nf) = api^(nf-1)*(1 + O(api^(nf-1)) + ...)
*)

setdec = {

    as6to5os ->
  1 + (7*api^2)/24 + (58933*api^3)/124416 + (api*lmm)/6 + (19*api^2*lmm)/24 +
   (8941*api^3*lmm)/1728 + (api^2*lmm^2)/36 + (511*api^3*lmm^2)/576 +
   (api^3*lmm^3)/216 - (2479*api^3*nl)/31104 - (409*api^3*lmm*nl)/1728 +
   (2*api^3*z2)/3 - (api^3*nl*z2)/9 + (2*api^3*z2*Log[2])/9 +
   (80507*api^3*Zeta[3])/27648,

    as6to5ms ->
  1 + (api*lmm)/6 + api^2*(-11/72 + (11*lmm)/24 + lmm^2/36) + 
   api^3*(-564731/124416 + (2645*lmm)/1728 + (167*lmm^2)/576 + lmm^3/216 + 
   (2633/31104 - (67*lmm)/576 + lmm^2/36)*nl + (82043*Zeta[3])/27648),

    as6to5si ->
  1 + (api*lmmu)/6 + api^2*(-11/72 + (19*lmmu)/24 + lmmu^2/36) + 
   api^3*(-564731/124416 + (2191*lmmu)/576 + (511*lmmu^2)/576 + lmmu^3/216 + 
   (2633/31104 - (281*lmmu)/1728)*nl + (82043*Zeta[3])/27648),

    as5to6os ->
  1 - (api*lmm)/6 + api^2*(-7/24 - (19*lmm)/24 + lmm^2/36) +
   api^3*(-58933/124416 - (8521*lmm)/1728 - (131*lmm^2)/576 - lmm^3/216 +
      nl*(2479/31104 + (409*lmm)/1728 + z2/9) - (2*z2)/3 - (2*z2*Log[2])/9 -
      (80507*Zeta[3])/27648),

    as5to6ms ->
  1 - (api*lmm)/6 + api^2*(11/72 - (11*lmm)/24 + lmm^2/36) +
   api^3*(564731/124416 - (955*lmm)/576 + (53*lmm^2)/576 - lmm^3/216 +
	(-2633/31104 + (67*lmm)/576 - lmm^2/36)*nl - (82043*Zeta[3])/27648),

    as5to6si ->
  1 - (api*lmmu)/6 + api^2*(11/72 - (19*lmmu)/24 + lmmu^2/36) + 
   api^3*(564731/124416 - (6793*lmmu)/1728 - (131*lmmu^2)/576 - lmmu^3/216 + 
   (-2633/31104 + (281*lmmu)/1728)*nl - (82043*Zeta[3])/27648),

    mq5to6os ->
  1 + api^2*(89/432 - (5*lmm)/36 + lmm^2/12) +
   api^3*(1871/2916 - B4/36 + (121*lmm)/2592 + (319*lmm^2)/432 +
      (29*lmm^3)/216 + nl*(1327/11664 - (53*lmm)/432 - lmm^3/108 -
         (2*z3)/27) - (407*z3)/864 - (5*lmm*z3)/6 + (5*z4)/4),

    mq5to6ms ->
  1 + api^2*(89/432 - (5*lmm)/36 + lmm^2/12) + 
   api^3*(2951/2916 - B4/36 + (175*lmm^2)/432 + (29*lmm^3)/216 + 
   lmm*(-311/2592 - (5*z3)/6) + nl*(1327/11664 - (53*lmm)/432 - lmm^3/108 - 
     (2*z3)/27) - (407*z3)/864 + (5*z4)/4),

    mq5to6si ->
  1 + api^2*(89/432 - (5*lmmu)/36 + lmmu^2/12) + 
   api^3*(2951/2916 - B4/36 - (1031*lmmu)/2592 + (319*lmmu^2)/432 + 
   (29*lmmu^3)/216 + nl*(1327/11664 - (53*lmmu)/432 - lmmu^3/108 - 
     (2*z3)/27) - (407*z3)/864 - (5*lmmu*z3)/6 + (5*z4)/4),

    mq6to5os ->
  1 + api^2*(-89/432 + (5*lmm)/36 - lmm^2/12) +
   api^3*(-1871/2916 + B4/36 - (299*lmm)/2592 - (299*lmm^2)/432 -
      (35*lmm^3)/216 - (1327*nl)/11664 + (53*lmm*nl)/432 + (lmm^3*nl)/108 +
      (407*z3)/864 + (5*lmm*z3)/6 + (2*nl*z3)/27 - (5*z4)/4),

    mq6to5ms ->
  1 + api^2*(-89/432 + (5*lmm)/36 - lmm^2/12) + 
   api^3*(-2951/2916 + B4/36 - (155*lmm^2)/432 - (35*lmm^3)/216 + 
   nl*(-1327/11664 + (53*lmm)/432 + lmm^3/108 + (2*z3)/27) + 
   lmm*(133/2592 + (5*z3)/6) + (407*z3)/864 - (5*z4)/4),

    mq6to5si ->
  1 + api^2*(-89/432 + (5*lmmu)/36 - lmmu^2/12) + 
   api^3*(-2951/2916 + B4/36 + (853*lmmu)/2592 - (299*lmmu^2)/432 - 
   (35*lmmu^3)/216 + nl*(-1327/11664 + (53*lmmu)/432 + lmmu^3/108 + 
     (2*z3)/27) + (407*z3)/864 + (5*lmmu*z3)/6 - (5*z4)/4)

   };


setMS2OS = {

  ms2msc ->  1 + ( g1/b0 - b1*g0/b0^2 ) * x
			  + 1/2 * ( (g1/b0-b1*g0/b0^2)^2 + g2/b0 
				   + b1^2*g0/b0^3 - b1*g1/b0^2 
				   - b2*g0/b0^2 ) * x^2
		       + ( 1/6*(g1/b0-b1*g0/b0^2)^3
			  +1/2*(g1/b0-b1*g0/b0^2)
			  *(g2/b0+b1^2*g0/b0^3
			    -b1*g1/b0^2-b2*g0/b0^2)
			  +1/3*(g3/b0-b1^3*g0/b0^4
				+2*b1*b2*g0/b0^3
				-b3*g0/b0^2+b1^2*g1/b0^3
				-b2*g1/b0^2
				-b1*g2/b0^2) ) * x^3 ,
    msfromos -> 1 + 
	api*(-cf - (3*cf*lmM)/4) + api^2*((-1111*ca*cf)/384 + (7*cf^2)/128 -
   (185*ca*cf*lmM)/96 + (21*cf^2*lmM)/32 - (11*ca*cf*lmM^2)/32 +
   (9*cf^2*lmM^2)/32 + (143*cf*tr)/96 + (13*cf*lmM*tr)/24 + (cf*lmM^2*tr)/8 +
   (71*cf*nl*tr)/96 + (13*cf*lmM*nl*tr)/24 + (cf*lmM^2*nl*tr)/8 +
   (ca*cf*z2)/2 - (15*cf^2*z2)/8 - (3*ca*cf*log2*z2)/2 + 3*cf^2*log2*z2 -
   cf*tr*z2 + (cf*nl*tr*z2)/2 + (3*ca*cf*z3)/8 - (3*cf^2*z3)/4) +
 api^3*((lmM^2*(-2341*ca^2*cf + 1962*ca*cf^2 - 243*cf^3 + 1492*ca*cf*tr -
      468*cf^2*tr + 1492*ca*cf*nl*tr - 468*cf^2*nl*tr - 208*cf*tr^2 -
      416*cf*nl*tr^2 - 208*cf*nl^2*tr^2))/1152 +
   (lmM^3*(-242*ca^2*cf + 297*ca*cf^2 - 81*cf^3 + 176*ca*cf*tr -
      108*cf^2*tr + 176*ca*cf*nl*tr - 108*cf^2*nl*tr - 32*cf*tr^2 -
      64*cf*nl*tr^2 - 32*cf*nl^2*tr^2))/1152 +
   (lmM*(-105944*ca^2*cf + 52317*ca*cf^2 - 13203*cf^3 + 74624*ca*cf*tr -
      5436*cf^2*tr + 55616*ca*cf*nl*tr + 2340*cf^2*nl*tr - 12608*cf*tr^2 -
      18304*cf*nl*tr^2 - 5696*cf*nl^2*tr^2 + 12672*ca^2*cf*z2 -
      52704*ca*cf^2*z2 + 19440*cf^3*z2 - 38016*ca^2*cf*log2*z2 +
      91584*ca*cf^2*log2*z2 - 31104*cf^3*log2*z2 - 29952*ca*cf*tr*z2 +
      27648*cf^2*tr*z2 + 13824*ca*cf*log2*tr*z2 - 27648*cf^2*log2*tr*z2 +
      8064*ca*cf*nl*tr*z2 + 12096*cf^2*nl*tr*z2 + 13824*ca*cf*log2*nl*tr*z2 -
      27648*cf^2*log2*nl*tr*z2 + 9216*cf*tr^2*z2 + 4608*cf*nl*tr^2*z2 -
      4608*cf*nl^2*tr^2*z2 + 9504*ca^2*cf*z3 - 22896*ca*cf^2*z3 +
      7776*cf^3*z3 + 6912*ca*cf*tr*z3 - 3456*cf^2*tr*z3 +
      6912*ca*cf*nl*tr*z3 - 3456*cf^2*nl*tr*z3))/13824
	),

    msfromos3set[0] -> -202,
    msfromos3set[1] -> -176,
    msfromos3set[2] -> -150,
    msfromos3set[3] -> -126,
    msfromos3set[4] -> -103,
    msfromos3set[5] -> -82,

    msfromos3 -> -9478333/93312 - (644201*Pi^2)/38880 + (695*Pi^4)/7776 + 
 (587*Pi^2*Log[2])/162 + (22*Pi^2*Log[2]^2)/81 + (55*Log[2]^4)/162 + 
 (220*PolyLog[4, 1/2])/27 + nl^2*(-2353/23328 - (13*Pi^2)/324 - 
   (7*Zeta[3])/54) - (61*Zeta[3])/27 + (1439*Pi^2*Zeta[3])/432 + 
 nl*(246643/23328 + (967*Pi^2)/648 - (61*Pi^4)/1944 + (11*Pi^2*Log[2])/81 - 
   (2*Pi^2*Log[2]^2)/81 - Log[2]^4/81 - (8*PolyLog[4, 1/2])/27 + 
   (241*Zeta[3])/72) - (1975*Zeta[5])/216,


    osfromms -> 1 + 
	api*(cf + (3*cf*lmm)/4) + api^2*((1111*ca*cf)/384 - (71*cf^2)/128 -
   (143*cf*tr)/96 - (71*cf*nl*tr)/96 +
   lmm*((185*ca*cf)/96 - (9*cf^2)/32 - (13*cf*tr)/24 - (13*cf*nl*tr)/24) +
   lmm^2*((11*ca*cf)/32 + (9*cf^2)/32 - (cf*tr)/8 - (cf*nl*tr)/8) -
   (ca*cf*z2)/2 + (15*cf^2*z2)/8 + (3*ca*cf*log2*z2)/2 - 3*cf^2*log2*z2 +
   cf*tr*z2 - (cf*nl*tr*z2)/2 - (3*ca*cf*z3)/8 + (3*cf^2*z3)/4) +
   api^3*(lmm^3*((121*ca^2*cf)/576 + (33*ca*cf^2)/128 +
     (9*cf^3)/128 - (11*ca*cf*tr)/72 - (3*cf^2*tr)/32 - (11*ca*cf*nl*tr)/72 -
     (3*cf^2*nl*tr)/32 + (cf*tr^2)/36 + (cf*nl*tr^2)/18 +
     (cf*nl^2*tr^2)/36) + lmm^2*((2341*ca^2*cf)/1152 + (21*ca*cf^2)/64 -
     (63*cf^3)/128 - (373*ca*cf*tr)/288 - (3*cf^2*tr)/32 -
     (373*ca*cf*nl*tr)/288 - (3*cf^2*nl*tr)/32 + (13*cf*tr^2)/72 +
     (13*cf*nl*tr^2)/36 + (13*cf*nl^2*tr^2)/72) +
   lmm*((13243*ca^2*cf)/1728 - (4219*ca*cf^2)/1536 + (495*cf^3)/512 -
     (583*ca*cf*tr)/108 - (307*cf^2*tr)/384 - (869*ca*cf*nl*tr)/216 -
     (91*cf^2*nl*tr)/384 + (197*cf*tr^2)/216 + (143*cf*nl*tr^2)/108 +
     (89*cf*nl^2*tr^2)/216 - (11*ca^2*cf*z2)/12 + (49*ca*cf^2*z2)/16 +
     (45*cf^3*z2)/32 + (11*ca^2*cf*log2*z2)/4 - (35*ca*cf^2*log2*z2)/8 -
     (9*cf^3*log2*z2)/4 + (13*ca*cf*tr*z2)/6 - (cf^2*tr*z2)/2 -
     ca*cf*log2*tr*z2 + 2*cf^2*log2*tr*z2 - (7*ca*cf*nl*tr*z2)/12 -
     (13*cf^2*nl*tr*z2)/8 - ca*cf*log2*nl*tr*z2 + 2*cf^2*log2*nl*tr*z2 -
     (2*cf*tr^2*z2)/3 - (cf*nl*tr^2*z2)/3 + (cf*nl^2*tr^2*z2)/3 -
     (11*ca^2*cf*z3)/16 + (35*ca*cf^2*z3)/32 + (9*cf^3*z3)/16 -
     (ca*cf*tr*z3)/2 + (cf^2*tr*z3)/4 - (ca*cf*nl*tr*z3)/2 +
     (cf^2*nl*tr*z3)/4)),

    osfrommsset -> 1 + 
	api*(cf + (3*cf*lmm)/4) + api^2*((1111*ca*cf)/384 - (71*cf^2)/128 -
   (143*cf*tr)/96 - (71*cf*nl*tr)/96 +
   lmm*((185*ca*cf)/96 - (9*cf^2)/32 - (13*cf*tr)/24 - (13*cf*nl*tr)/24) +
   lmm^2*((11*ca*cf)/32 + (9*cf^2)/32 - (cf*tr)/8 - (cf*nl*tr)/8) -
   (ca*cf*z2)/2 + (15*cf^2*z2)/8 + (3*ca*cf*log2*z2)/2 - 3*cf^2*log2*z2 +
   cf*tr*z2 - (cf*nl*tr*z2)/2 - (3*ca*cf*z3)/8 + (3*cf^2*z3)/4) +
   api^3*(lmm^3*((121*ca^2*cf)/576 + (33*ca*cf^2)/128 +
     (9*cf^3)/128 - (11*ca*cf*tr)/72 - (3*cf^2*tr)/32 - (11*ca*cf*nl*tr)/72 -
     (3*cf^2*nl*tr)/32 + (cf*tr^2)/36 + (cf*nl*tr^2)/18 +
     (cf*nl^2*tr^2)/36) + lmm^2*((2341*ca^2*cf)/1152 + (21*ca*cf^2)/64 -
     (63*cf^3)/128 - (373*ca*cf*tr)/288 - (3*cf^2*tr)/32 -
     (373*ca*cf*nl*tr)/288 - (3*cf^2*nl*tr)/32 + (13*cf*tr^2)/72 +
     (13*cf*nl*tr^2)/36 + (13*cf*nl^2*tr^2)/72) - (ca*cf^2*z2)/4 +
   (15*cf^3*z2)/16 + (3*ca*cf^2*log2*z2)/4 - (3*cf^3*log2*z2)/2 +
   (cf^2*tr*z2)/2 - (cf^2*nl*tr*z2)/4 - (3*ca*cf^2*z3)/16 + (3*cf^3*z3)/8 +
   lmm*((13243*ca^2*cf)/1728 - (4219*ca*cf^2)/1536 + (495*cf^3)/512 -
     (583*ca*cf*tr)/108 - (307*cf^2*tr)/384 - (869*ca*cf*nl*tr)/216 -
     (91*cf^2*nl*tr)/384 + (197*cf*tr^2)/216 + (143*cf*nl*tr^2)/108 +
     (89*cf*nl^2*tr^2)/216 - (11*ca^2*cf*z2)/12 + (49*ca*cf^2*z2)/16 +
     (45*cf^3*z2)/32 + (11*ca^2*cf*log2*z2)/4 - (35*ca*cf^2*log2*z2)/8 -
     (9*cf^3*log2*z2)/4 + (13*ca*cf*tr*z2)/6 - (cf^2*tr*z2)/2 -
     ca*cf*log2*tr*z2 + 2*cf^2*log2*tr*z2 - (7*ca*cf*nl*tr*z2)/12 -
     (13*cf^2*nl*tr*z2)/8 - ca*cf*log2*nl*tr*z2 + 2*cf^2*log2*nl*tr*z2 -
     (2*cf*tr^2*z2)/3 - (cf*nl*tr^2*z2)/3 + (cf*nl^2*tr^2*z2)/3 -
     (11*ca^2*cf*z3)/16 + (35*ca*cf^2*z3)/32 + (9*cf^3*z3)/16 -
     (ca*cf*tr*z3)/2 + (cf^2*tr*z3)/4 - (ca*cf*nl*tr*z3)/2 +
     (cf^2*nl*tr*z3)/4)),

    osfromms3set[0] -> 194,
    osfromms3set[1] -> 168,
    osfromms3set[2] -> 143,
    osfromms3set[3] -> 119,
    osfromms3set[4] -> 96,
    osfromms3set[5] -> 75,

    osfromms3 -> 8481925/93312 + 
 (137*nl)/216 + (652841*Pi^2)/38880 - (nl*Pi^2)/27 - 
 (695*Pi^4)/7776 - (575*Pi^2*Log[2])/162 - (22*Pi^2*Log[2]^2)/81 - 
 (55*Log[2]^4)/162 - (220*PolyLog[4, 1/2])/27 - 
 nl^2*(-2353/23328 - (13*Pi^2)/324 - (7*Zeta[3])/54) + (58*Zeta[3])/27 - 
 (1439*Pi^2*Zeta[3])/432 - nl*(246643/23328 + (967*Pi^2)/648 - 
   (61*Pi^4)/1944 + (11*Pi^2*Log[2])/81 - (2*Pi^2*Log[2]^2)/81 - 
   Log[2]^4/81 - (8*PolyLog[4, 1/2])/27 + (241*Zeta[3])/72) + 
 (1975*Zeta[5])/216,

    mumfromos -> 1 
   - apiM*cf + apiM^2*((-1111*ca*cf)/384 + (199*cf^2)/128 + (143*cf*tr)/96 +
   (71*cf*nl*tr)/96 + (ca*cf*z2)/2 - (15*cf^2*z2)/8 - (3*ca*cf*log2*z2)/2 +
   3*cf^2*log2*z2 - cf*tr*z2 + (cf*nl*tr*z2)/2 + (3*ca*cf*z3)/8 -
   (3*cf^2*z3)/4),

    mumfromos3set[0] -> -170,
    mumfromos3set[1] -> -146,
    mumfromos3set[2] -> -123,
    mumfromos3set[3] -> -101,
    mumfromos3set[4] -> -81,
    mumfromos3set[5] -> -62,

    mumfromos3 -> -7172965/93312 - 
	(293*nl)/216 - (618281*Pi^2)/38880 - (nl*Pi^2)/9 + 
 (695*Pi^4)/7776 + (623*Pi^2*Log[2])/162 + (22*Pi^2*Log[2]^2)/81 + 
 (55*Log[2]^4)/162 + (220*PolyLog[4, 1/2])/27 + 
 nl^2*(-2353/23328 - (13*Pi^2)/324 - (7*Zeta[3])/54) - (70*Zeta[3])/27 + 
 (1439*Pi^2*Zeta[3])/432 + nl*(246643/23328 + (967*Pi^2)/648 - 
   (61*Pi^4)/1944 + (11*Pi^2*Log[2])/81 - (2*Pi^2*Log[2]^2)/81 - 
   Log[2]^4/81 - (8*PolyLog[4, 1/2])/27 + (241*Zeta[3])/72) - 
 (1975*Zeta[5])/216,

    msfromri -> 1 - (4*api)/3 + api^2*(-995/72 + (89*nf)/144 + (19*z3)/6) + 
 api^3*(-6663911/41472 + (118325*nf)/7776 - (4459*nf^2)/23328 + 
   (5*nf*z4)/12 - (185*z5)/36 + (408007*z3)/6912 - 
   (617*nf*z3)/216 - (nf^2*z3)/54),

    rifromms -> 1 + (4*api)/3 + api^2*(1123/72 - (89*nf)/144 - (19*z3)/6) + 
 api^3*(6663911/41472 - (118325*nf)/7776 + (4459*nf^2)/23328 + 
   (4*(1123/72 - (89*nf)/144 - (19*z3)/6))/3 - (408007*z3)/6912 + 
   (617*nf*z3)/216 + (nf^2*z3)/54 - (4*(-995/72 + (89*nf)/144 + (19*z3)/6))/
    3 - (5*nf*z4)/12 + (185*z5)/36)

    };

(* ************************************************************ *)

(* Compute \\Lambda^(nf) using 3 methods:
   1. Explicit formulae for \\Lambda: lamexpl[alphas_,mu_,nn_,loops_]
   2. Implicit formulae for \\Lambda: lamimpl[alphas_,mu_,nn_,loops_]
   3. Solving the integral exactly:   lamexact[alphas_,mu_,nn_]
*)

(* Determine \Lamba for given \alpha_s				*)

LamExpl[alphas_,mu_,nn_,loops_] := Module[
    {logmu2lam2,nnn},
    If[(loops<1)||(loops>4), Print["Invalid # loops."];
       Return[];];
    (* \log (\mu^2/\Lambda^2) =	
     *)
    logmu2lam2 = 1/(as*b0) + (as*(-b1^2 + b0*b2))/b0^3 +
	(as^2*(b1^3 - 2*b0*b1*b2 + b0^2*b3))/(2*b0^4) + (b1*Log[as])/b0^2 +
	    (b1*Log[b0])/b0^2;
    (* \Lambda/\mu = Exp[ lamomu[nn_] ]; nn: # active flavours
     *)
    lamomu[nnn_]:=Collect[Expand[ 
	-1/2*logmu2lam2/.setbeta/.nf->nnn
	],as]/.as:>ToExpression["api"<>ToString[nnn]];
    Return[N[mu*Exp[
	Expand[
	    (Expand[ lamomu[nn] * ToExpression["api"<>ToString[nn]]^2 ]
	     /.cut[ToExpression["api"<>ToString[nn]],loops]
	     )
		*1/ToExpression["api"<>ToString[nn]]^2
		]
	/.ToExpression["api"<>ToString[nn]]->alphas/Pi] 
	     ,numprec]];
];

(* ************************************************************ *)

LamImpl[alphas_,mu_,nn_,loops_] := Module[
    {lam,astmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==1, astmp=as0];
    If[loops==2, astmp=as1];
    If[loops==3, astmp=as2];
    If[loops==4, astmp=as3];
    Return[lam/.FindRoot[
	(Expand[
	    astmp/.setasL
		/.setbeta/.nf->nn]/.L->Log[mu^2/lam^2])==alphas/Pi,
	    {lam,LamExpl[alphas,mu,nn,loops]} 
    ]
       ]
];

(* ************************************************************ *)

(* Compute \alpha_s(\mu);
   Input: \Lambda, \mu, # active flavours and # loops
*)

AlphasLam[lambda_,mu_,nn_,loops_] := Module[
    {},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[ N[mu/lambda] < rmulam ,
	Print["WARNING: the ratio \\mu/\\Lambda = ",
	      N[mu/lambda]," is very small!"];
    ];
    If[loops==1, Return[ N[Pi*as0/.setasL/.setbeta/.nf->nn
			       /.L->Log[mu^2/lambda^2],numprec]] ];
    If[loops==2, Return[ N[Pi*as1/.setasL/.setbeta/.nf->nn
			       /.L->Log[mu^2/lambda^2],numprec]] ];
    If[loops==3, Return[ N[Pi*as2/.setasL/.setbeta/.nf->nn
			       /.L->Log[mu^2/lambda^2],numprec]] ];
    If[loops==4, Return[ N[Pi*as3/.setasL/.setbeta/.nf->nn
			       /.L->Log[mu^2/lambda^2],numprec]] ];
];

(* ************************************************************ *)

(* Compute \alpha_s(\mu);
   Solve differential equation numerically
   Input: \alpha_s(\mu0), \mu0, \mu, # active flavours and # loops
*)

AlphasExact[alphasmu0_,mu0_,muend_,nnf_,loops_] := Module[
    {rhs,solx,a,x,lambda,alphasmu0tmp},
    alphasmu0tmp = Rationalize[alphasmu0,10^-numprec];
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    lambda = LamExpl[alphasmu0,mu0,nnf,loops];
    If[ N[muend/lambda] < rmulam ,
	Print["WARNING: the ratio \\mu/\\Lambda = ",
	      N[muend/lambda]," is very small!"];
    ];
    If[ mu0 == muend , Return[alphasmu0] ];
    If[ loops==1, rhs = -a[x]^2*(b0) ];
    If[ loops==2, rhs = -a[x]^2*(b0 + b1*a[x]) ];
    If[ loops==3, rhs = -a[x]^2*(b0 + b1*a[x] + b2*a[x]^2) ];
    If[ loops==4, rhs = -a[x]^2*(b0 + b1*a[x] + b2*a[x]^2 + b3*a[x]^3) ];
    sol = N[(Pi*a[x]/.
	     NDSolve[{a[mu0^2]==alphasmu0tmp/Pi,
		      x*a'[x] == rhs
         /.setbeta/.nf->nnf}, 
         a[x], 
         {x,muend^2,mu0^2},        WorkingPrecision->26,
                                   AccuracyGoal->16,
                                   PrecisionGoal->16,
                                   MaxSteps->1500]
/.x->muend^2)[[1]],numprec];
Return[sol];
];

(* ************************************************************ *)

(*--#[ AsmMSrunexact: *)

(* Solve simultaneously RGE for \alpha_s and m_q. *)
AsmMSrunexact[mmu0_, asmu0_, mu0_, mu_, nf_, l_] := Module[
    {res,b0,b1,g0,g1,beta,gammam,x,xl,api},
    If[l>4,Print["AsmMSrunexact: l>4 not (yet) implemented."];Abort[];];
    If[l<=0,Return[{mmu0,asmu0}]];
    beta   = Expand[-api[x]^2*(xl*b0+xl^2*api[x]*b1+
			       xl^3*api[x]^2*b2+xl^4*api[x]^3*b3)
		    ]/.cut[xl,l]/.{xl->1};
    gammam = Expand[-api[x]  *(xl*g0+xl^2*api[x]*g1+
			       xl^3*api[x]^2*g2+xl^4*api[x]^3*g3)
		    ]/.cut[xl,l]/.{xl->1};
    b0 = 11/4 - nf/6;
    b1 = 51/8 - (19*nf)/24;
    b2 = 2857/128 - (5033*nf)/1152 + (325*nf^2)/3456;
    (b3 = 149753/1536 - (1078361*nf)/41472 + (50065*nf^2)/41472 +
     (1093*nf^3)/186624 + (891*Zeta[3])/64 - (1627*nf*Zeta[3])/1728 +
     (809*nf^2*Zeta[3])/2592);
    g0 = 1;
    g1 = (101/2 - (5*nf)/3)/12;
    g2 = (3747/4 - (nf*(1108/9 + (70*nf)/27 + 80*Zeta[3]))/2)/48;
    (g3 = 4603055/41472 - (91723*nf)/6912 + (2621*nf^2)/31104 - 
     (83*nf^3)/15552 +
     (11*nf*Pi^4)/288 - (nf^2*Pi^4)/432 + (530*Zeta[3])/27 -
     (2137*nf*Zeta[3])/144 + (25*nf^2*Zeta[3])/72 + (nf^3*Zeta[3])/108 -
     (275*Zeta[5])/8 + (575*nf*Zeta[5])/72);
    
    res = NDSolve[ {x*api'[x]      == beta,
		    x/mq[x]*mq'[x] == gammam,
		    api[mu0^2]     == asmu0/Pi, 
		    mq[mu0^2]      == mmu0},    {api[x],mq[x]}, {x,mu^2,mu0^2}
(*,
			WorkingPrecision->26,
			AccuracyGoal->16,
			PrecisionGoal->16,
			MaxSteps->5000
*)
		    ];
    Return[{mq[x],Pi*api[x]}/.res[[1]]/.{x->mu^2}];
];

(*--#] *)

(* ************************************************************ *)

(*
 Light-mass corrections at order \alpha_s^2
*)

delta[mOS_,mq_] := Module[
    {i,del,delta,r},
    del=0;
    delta[r_] := If[(r<=1) && (r>=0),
		    Pi^2/8*r-0.597*r^2+0.230*r^3,
		    Print["\\Delta(",N[r],
                          ") IS CALLED; THE FUNCTION IS NOT IMPLEMENTED 
FOR ARGUMENTS OUTSIDE THE INTERVAL [0,1]."];
		    Abort[]
		];
    For[i=1,i<=Length[mq],i++,
	del=del+delta[mq[[i]]/mOS];
    ];
    Return[del];
];

(* ************************************************************ *)

(*
 m_MS = M_OS * ( 1 + ... )
*)

mOS2mMS[mOS_,mq_,asmu_,mu_,nnf_,loops_] := Module[
    {extmp},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,extmp=1,
       extmp = 
	   ( msfromos 
	    + api^2*(-4/3 * delta[mOS,mq])
	    (* + api^3*msfromos3set[nnf-1] /.setMS2OS *)
	    + api^3*msfromos3 /.setMS2OS
	    )/.cut[api,loops] /.{lmM->Log[mu^2/mOS^2],nl->nnf-1
				 } /.num1 /.cf2num;
   ];
    Return[ N[ mOS * (extmp/.api->asmu/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 m_OS = M_MS * ( 1 + ... )
*)

mMS2mOS[mMS_,mq_,asmu_,mu_,nnf_,loops_] := Module[
    {extmp},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,extmp=1,
       extmp = 
	   ( osfromms 
	    + api^2*(4/3 * delta[mMS,mq])
	    (* + api^3*osfromms3set[nnf-1] /.setMS2OS *)
	    + api^3*osfromms3 /.setMS2OS
          )/.cut[api,loops] /.{lmm->Log[mu^2/mMS^2],nl->nnf-1} /.num1 /.cf2num;
   ];
    Return[ N[ mMS * (extmp/.api->asmu/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 Compute in a first step m_MS(m_MS) and afterwards m_MS(mu).
*)

mOS2mMSrun[mOS_,mq_,asmu_,mu_,nnf_,loops_] := Module[
    {mum,asM,asmum},
    asM = AlphasExact[asmu,mu,mOS,nnf,loops];
    mum = mOS2mSI[mOS,mq,asM,nnf,loops];
    asmum = AlphasExact[asmu,mu,mum,nnf,loops];
    Return[ N[ mMS2mMS[mum,asmum,asmu,nnf,loops] , numprec ] ];
];

(* ************************************************************ *)

(*
 Compute in a first step m_MS(m_MS). Then use m_OS = m_MS * ( 1 + ... )
 for mu=m_MS.
*)

mMS2mOSrun[mMS_,mq_,asmu_,mu_,nnf_,loops_] := Module[
    {asmMS,mMS2mSI2,xx,res1},
    asmMS = AlphasLam[LamImpl[asmu,mu,nnf,loops],xx,nnf,loops];
    mMS2mSI2 = mMS2mMS[mMS,asmu,asmMS,nnf,loops];
    res1 = FindRoot[ xx == mMS2mSI2 , {xx,mMS} ];
    (* print mMS2mSI2 and \alpha_s(mMS2mSI2)
     Print[ mMS2mSI2/.res1 ];
     Print[ asmMS/.res1 ];
     *)
    Return[ N[ mMS2mOS[mMS2mSI2/.res1,mq,asmMS/.res1,mMS2mSI2/.res1,nnf,loops] 
	       , numprec ] ];
];

(* ************************************************************ *)

mOS2mMSit[mOS_,mq_,asmu_,mu_,nnf_,loops_] := Module[
    {extmp},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,extmp=1,
       extmp = 
	   ( osfromms 
	    + api^2*(4/3 * delta[mOS,mq])
	    (* + api^3*osfromms3set[nnf-1] /.setMS2OS *)
	    + api^3*osfromms3 /.setMS2OS
          )/.cut[api,loops] /.{lmm->Log[mu^2/mMS^2],nl->nnf-1} /.num1 /.cf2num;
   ];
    Return[ N[ mMS/.FindRoot[ mOS == mMS * (extmp/.api->asmu/Pi) , {mMS,mOS}]
	       , numprec ] ];
];

(* ************************************************************ *)

(*
 \mu_m = M_OS * ( 1 + ... )
*)

mOS2mSI[mOS_,mq_,asM_,nnf_,loops_] := Module[
    {extmp},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,extmp=1,
       extmp = 
	   ( mumfromos 
	    + apiM^2*(-4/3 * delta[mOS,mq])
	    (* + apiM^3*mumfromos3set[nnf-1] /.setMS2OS *)
	    + apiM^3*mumfromos3 /.setMS2OS
          )/.cut[apiM,loops] /.{nl->nnf-1} /.num1 /.cf2num;
   ];
    Return[ N[ mOS * (extmp/.apiM->asM/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 Running MS-bar mass
*)

mMS2mMS[mmu0_,asmu0_,asmu_,nnf_,loops_] := Module[
    {ccc},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,
       Return[ N[mmu0, numprec] ]
   ];
    If[loops==1, 
       ccc = x^(g0/b0) /.setbeta /.setgamma /.nf->nnf;
       ,
       ccc = (x^(g0/b0)*
	      (Expand[ ms2msc /.setMS2OS ]/.cut[x,loops-1])
	      ) /.setbeta /.setgamma /.nf->nnf;
   ];
    Return[ N[ mmu0 * (ccc/.x->asmu/Pi)/(ccc/.x->asmu0/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 Compute in a first step m_MS(m_MS). Then use m_OS = m_MS * ( 1 + ... )
 for mu=m_MS.
*)

mMS2mSI[mMS_,asmu_,mu_,nnf_,loops_] := Module[
    {asmMS,mMS2mSI2,xx,res1},
    asmMS = AlphasLam[LamImpl[asmu,mu,nnf,loops],xx,nnf,loops];
    mMS2mSI2 = mMS2mMS[mMS,asmu,asmMS,nnf,loops];
    res1 = FindRoot[ xx == mMS2mSI2 , {xx,mMS} ];

    Return[ N[ mMS2mSI2/.res1 , numprec ] ];
];

(* ************************************************************ *)

(*
 m_RI = m_MS * ( 1 + ... )
*)

mMS2mRI[mMS_,asmu_,nnf_,loops_] := Module[
    {extmp},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,extmp=1,
       extmp = 
	   ( rifromms 
	    )/.setMS2OS/.cut[api,loops] /.{nf->nnf} /.num1 /.cf2num;
   ];
    Return[ N[ mMS * (extmp/.api->asmu/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 m_MS = m_RI * ( 1 + ... )
*)

mRI2mMS[mRI_,asmu_,nnf_,loops_] := Module[
    {extmp},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,extmp=1,
       extmp = 
	   ( msfromri 
	    )/.setMS2OS/.cut[api,loops] /.{nf->nnf} /.num1 /.cf2num;
   ];
    Return[ N[ mRI * (extmp/.api->asmu/Pi) , numprec ] ];
];

(* ************************************************************ *)
(*
 m_RGI = m_MS * ( 1 + ... )
*)
(*
mMS2mRGI[1,0.1,10,5,1]
*)
mMS2mRGI[mMS_,asmu_,nnf_,loops_] := Module[
    {ccc},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,
       Return[ N[mMS, numprec] ]
   ];
    If[loops==1, 
       ccc = x^(g0/b0) /.setbeta /.setgamma /.nf->nnf;
       ,
       ccc = (x^(g0/b0)*
	      (Expand[ ms2msc /.setMS2OS ]/.cut[x,loops-1])
	      ) /.setbeta /.setgamma /.nf->nnf;
   ];
    Return[ N[ mMS / (ccc/.x->asmu/Pi) , numprec ] ];
];

mMS2mRGImod[mMS_,asmu_,nnf_,loops_] := Module[
    {ccc,bet0},
    bet0 = 11/4 - nnf/6;
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,
       Return[ N[mmu0, numprec] ]
   ];
    If[loops==1, 
       ccc = x^(g0/b0) /.setbeta /.setgamma /.nf->nnf;
       ,
       ccc = (x^(g0/b0)*
	      (Expand[ ms2msc /.setMS2OS ]/.cut[x,loops-1])
	      ) /.setbeta /.setgamma /.nf->nnf;
   ];
    Return[ N[ mMS / (ccc/.x-> 2*bet0*asmu/Pi ) , numprec ] ];
(*                             ^^^^^^ *)
];


(* ************************************************************ *)

(*
 m_MS = m_RGI * ( 1 + ... )
*)

mRGI2mMS[mRGI_,asmu_,nnf_,loops_] := Module[
    {ccc},
    If[(loops<0)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    If[loops==0,
       Return[ N[mmu0, numprec] ]
   ];
    If[loops==1, 
       ccc = x^(g0/b0) /.setbeta /.setgamma /.nf->nnf;
       ,
       ccc = (x^(g0/b0)*
	      (Expand[ ms2msc /.setMS2OS ]/.cut[x,loops-1])
	      ) /.setbeta /.setgamma /.nf->nnf;
   ];
    Return[ N[ mRGI * (ccc/.x->asmu/Pi) , numprec ] ];
];

(* ************************************************************ *)

(*
 decoupling of \alpha_s: on-shell scheme (for heavy mass)
 compute: \alpha_s^(nf+1)
 input: \alpha_s^(nf), ...
 the decoupling is performed at "loops-1" order
*)

DecAsUpOS[alsp_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as6to5os/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as6to5os/.setdec/.api->0 ];
    Return[N[
	alsp*(dectmp/.num1/.nl->nnl/.api->alsp/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(* 
 decoupling of \alpha_s: MS-bar scheme (for heavy mass)
*)

DecAsUpMS[alsp_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as6to5ms/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as6to5ms/.setdec/.api->0 ];
    Return[N[
	alsp*(dectmp/.num1/.nl->nnl/.api->alsp/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(* 
 decoupling of \alpha_s: scale invariant mass (for heavy mass)
*)

DecAsUpSI[alsp_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as6to5si/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as6to5si/.setdec/.api->0 ];
    Return[N[
	alsp*(dectmp/.num1/.nl->nnl/.api->alsp/Pi/.lmmu->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of \alpha_s: on-shell scheme (for heavy mass)
 compute: \alpha_s^(nf-1)
 input: \alpha_s^(nf), ...
 the decoupling is performed at "loops-1" order
*)

DecAsDownOS[als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as5to6os/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as5to6os/.setdec/.api->0 ];
    Return[N[
	als*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(* 
 decoupling of \alpha_s: MS-bar scheme (for heavy mass)
*)

DecAsDownMS[als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as5to6ms/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as5to6ms/.setdec/.api->0 ];
    Return[N[
	als*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(* 
 decoupling of \alpha_s: scale invariant mass (for heavy mass)
*)

DecAsDownSI[als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ as5to6si/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=as5to6si/.setdec/.api->0 ];
    Return[N[
	als*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmmu->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of quark mass: on-shell scheme (for heavy mass)
 compute: m_q^(nf)
 input: m_q^(nf-1), ...
 the decoupling is performed at "loops-1" order
*)

DecMqUpOS[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq6to5os/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq6to5os/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of quark mass: MS-bar scheme (for heavy mass)
 compute: m_q^(nf)
 input: m_q^(nf-1), ...
*)

DecMqUpMS[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq6to5ms/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq6to5ms/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of quark mass: scale invariant mass (for heavy mass)
 compute: m_q^(nf)
 input: m_q^(nf-1), ...
*)

DecMqUpSI[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq6to5si/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq6to5si/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmmu->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decouplingof quark mass: on-shell scheme (for heavy mass)
 compute: m_q^(nf-1)
 input: m_q^(nf), ...
*)

DecMqDownOS[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq5to6os/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq5to6os/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];


(* ************************************************************ *)

(*
 decoupling of quark mass: MS-bar scheme (for heavy mass)
 compute: m_q^(nf-1)
 input: m_q^(nf), ...
*)

DecMqDownMS[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq5to6ms/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq5to6ms/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmm->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of quark mass: scale invariant mass (for heavy mass)
 compute: m_q^(nf-1)
 input: m_q^(nf), ...
*)

DecMqDownSI[mq_,als_,massth_,muth_,nnl_,loops_] := Module[
    {dectmp},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    dectmp = Expand[ mq5to6si/.setdec ]/.cut[api,loops-1];
    If[ loops==1 , dectmp=mq5to6si/.setdec/.api->0 ];
    Return[N[
	mq*(dectmp/.num1/.nl->nnl/.api->als/Pi/.lmmu->Log[muth^2/massth^2]),
	numprec]];
];

(* ************************************************************ *)

(*
 decoupling of \Lambda
 compute: \Lambda^(nl+1)
 input: \Lambda^(nl), ...
*)

DecLambdaUp[lam_,massth_,nnl_,loops_] := Module[
    {exra,exra2,iLb0m,ex1,ex2,ex3,ex4},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    exra = Expand[(L*(1/2 - b0[-1 + nf]/b0times2[nf]) +
 ((-b1p[-1 + nf] + b1p[nf])*Log[L] +
   (-b1p[-1 + nf]^2 + b1p[nf]^2 + b2p[-1 + nf] - b2p[nf] + c2[nf, -1 + nf] +
     (-b1p[-1 + nf]^2 + b1p[-1 + nf]*b1p[nf])*Log[L])/Lb0m +
   (-b1p[-1 + nf]^3/2 - b1p[nf]^3/2 + b3p[-1 + nf]/2 - b3p[nf]/2 +
     b1p[nf]*(b1p[-1 + nf]^2 - b2p[-1 + nf] + b2p[nf] - c2[nf, -1 + nf]) +
     c3[nf, -1 + nf] + (b1p[-1 + nf]^2*b1p[nf] - b1p[-1 + nf]*b1p[nf]^2 +
       b1p[-1 + nf]*(-b2p[-1 + nf] + b2p[nf] - c2[nf, -1 + nf]))*Log[L] +
     (b1p[-1 + nf]^3/2 - (b1p[-1 + nf]^2*b1p[nf])/2)*Log[L]^2)/Lb0m^2 +
  b1p[nf]*Log[b0[-1 + nf]/b0[nf]])/b0times2[nf]
		   )/.{b0times2[nn_]->2*b0[nn]}
		  /.{Lb0m->(L*b0[nf-1])}
	      ];
    ex1 = Coefficient[exra,L,1];
    ex2 = Coefficient[exra,L,0];
    ex3 = Coefficient[exra,L,-1];
    ex4 = Coefficient[exra,L,-2];

    If[ loops==4 , exra = ex1*L+ex2+ex3/L+ex4/L^2 ];
    If[ loops==3 , exra = ex1*L+ex2+ex3/L ];
    If[ loops==2 , exra = ex1*L+ex2 ];
    If[ loops==1 , exra = ex1*L ];

    exra2 = Expand[ exra  
		    /.{nf->nnl+1}
		    /.{b0[nn_] -> (b0/.{setbeta/.nf->nn}),
		       b1p[nn_]-> ((b1/b0)/.{setbeta/.nf->nn}),
		       b2p[nn_]-> ((b2/b0)/.{setbeta/.nf->nn}),
		       b3p[nn_]-> ((b3/b0)/.{setbeta/.nf->nn})}
		    /.{c2[nn__]->11/72}
		    /.{c3[mm_,nn_]->-564731/124416+82043/27648*Zeta[3]
		       +2633/31104*nn}
		][[1]];
    Return[ lam*Exp[exra2 /. L->Log[massth^2/lam^2]] ] ;
];

(* ************************************************************ *)

(*
 decoupling of \Lambda
 compute: \Lambda^(nl)
 input: \Lambda^(nl+1), ...
*)

DecLambdaDown[lam_,massth_,nnl_,loops_] := Module[
    {exra,exra2,iLb0},
    If[(loops<1)||(loops>4), Print["PROCEDURE IS NOT IMPLEMENTED FOR ",
				   loops," LOOPS."];
       Return[];];
    exra = Expand[(L*(1/2 - b0[nf]/b0times2[-1 + nf]) +
 ((b1p[-1 + nf] - b1p[nf])*Log[L] + (b1p[-1 + nf]^2 - b1p[nf]^2 -
     b2p[-1 + nf] + b2p[nf] + c2[-1 + nf, nf] +
     (b1p[-1 + nf]*b1p[nf] - b1p[nf]^2)*Log[L])/Lb0 +
   (-b1p[-1 + nf]^3/2 - b1p[nf]^3/2 - b3p[-1 + nf]/2 + b3p[nf]/2 +
     b1p[-1 + nf]*(b1p[nf]^2 + b2p[-1 + nf] - b2p[nf] - c2[-1 + nf, nf]) +
     c3[-1 + nf, nf] + (-(b1p[-1 + nf]^2*b1p[nf]) + b1p[-1 + nf]*b1p[nf]^2 +
       b1p[nf]*(b2p[-1 + nf] - b2p[nf] - c2[-1 + nf, nf]))*Log[L] +
     (-(b1p[-1 + nf]*b1p[nf]^2)/2 + b1p[nf]^3/2)*Log[L]^2)/Lb0^2 +
   b1p[-1 + nf]*Log[b0[nf]/b0[-1 + nf]])/b0times2[-1 + nf]
		   )/.{b0times2[nn_]->2*b0[nn]}
		  /.{Lb0->(L*b0[nf])}
	      ];
    
    ex1 = Coefficient[exra,L,1];
    ex2 = Coefficient[exra,L,0];
    ex3 = Coefficient[exra,L,-1];
    ex4 = Coefficient[exra,L,-2];

    If[ loops==4 , exra = ex1*L+ex2+ex3/L+ex4/L^2 ];
    If[ loops==3 , exra = ex1*L+ex2+ex3/L ];
    If[ loops==2 , exra = ex1*L+ex2 ];
    If[ loops==1 , exra = ex1*L ];
    
    exra2 = Expand[ exra  /.{iLb0->1/(L*b0[nf])}
		    /.{nf->nnl+1}
		    /.{b0[nn_] -> (b0/.{setbeta/.nf->nn}),
		       b1p[nn_]-> ((b1/b0)/.{setbeta/.nf->nn}),
		       b2p[nn_]-> ((b2/b0)/.{setbeta/.nf->nn}),
		       b3p[nn_]-> ((b3/b0)/.{setbeta/.nf->nn})}
		    /.{c2[nn__]->11/72}
		    /.{c3[mm_,nn_]->-564731/124416+82043/27648*Zeta[3]
		       +2633/31104*nn}
		][[1]];
    
    Return[lam*Exp[exra2 /. L->Log[massth^2/lam^2]] ] ;
];

(* ************************************************************ *)

(* 
 running-decoupling-running-decoupling-running-...
 input: \alpha_s^(l)(mu1), ...
 output: \alpha_s^(h)(mu2)
*)

AlL2AlH[asl_,mu1_,decpar_,mu2_,loops_] := Module[
    {i,asini,muini,res1,res2,res3,decpar2},
    asini=asl;
    muini=mu1;
    decpar2 = Sort[decpar];
    For[i=2,i<=Length[decpar2],i++,
	If[(decpar2[[i]][[1]]-decpar2[[i-1]][[1]]=!=1) ||
	   (decpar2[[i]][[1]]-decpar2[[i-1]][[1]]==0) ,
	   Print["WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT."];
	   Abort[];
       ];
    ];
    For[i=1,i<=Length[decpar2],i++,
	res1 = AlphasExact[asini,muini,decpar2[[i]][[3]],
			   decpar2[[i]][[1]]-1,loops];
	res2 = DecAsUpOS[res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			 decpar2[[i]][[1]]-1,loops];
	asini = res2;
	muini = decpar2[[i]][[3]];
    ];
    res3 = AlphasExact[asini,muini,mu2,decpar2[[i-1]][[1]],loops];
    Return[res3];
];

(* ************************************************************ *)

(* 
 running-decoupling-running-decoupling-running-...
 input: \alpha_s^(h)(mu1), ...
 output: \alpha_s^(l)(mu2)
*)

AlH2AlL[ash_,mu1_,decpar_,mu2_,loops_] := Module[
    {i,asini,muini,res1,res2,res3,decpar2},
    asini=ash;
    muini=mu1;
    decpar2 = Reverse[Sort[decpar]];
    For[i=2,i<=Length[decpar2],i++,
	If[(decpar2[[i]][[1]]-decpar2[[i-1]][[1]]=!=-1) ||
	   (decpar2[[i]][[1]]-decpar2[[i-1]][[1]]==0) ,
	   Print["WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT."];
	   Abort[];
       ];
    ];
    For[i=1,i<=Length[decpar2],i++,
	res1 = AlphasExact[asini,muini,decpar2[[i]][[3]],
			   decpar2[[i]][[1]],loops];
	res2 = DecAsDownOS[res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			   decpar2[[i]][[1]]-1,loops];
	asini = res2;
	muini = decpar2[[i]][[3]];
    ];
    res3 = AlphasExact[asini,muini,mu2,decpar2[[i-1]][[1]]-1,loops];
    Return[res3];
];

(* ************************************************************ *)

(* 
 example: running-decoupling-running-decoupling-running-...
 input: m_q^(l)(mu1), ...
 output: m_q^(h)(mu2)
*)

mL2mH[mql_,asl_,mu1_,decpar_,mu2_,loops_] := Module[
    {i,asini,muini,res1,res2,res3,res4,res5,res6,decpar2},
    asini=asl;
    muini=mu1;
    mqini=mql;
    decpar2 = Sort[decpar];
    For[i=2,i<=Length[decpar2],i++,
	If[(decpar2[[i]][[1]]-decpar2[[i-1]][[1]]=!=1) ||
	   (decpar2[[i]][[1]]-decpar2[[i-1]][[1]]==0) ,
	   Print["WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT."];
	   Abort[];
       ];
    ];
    For[i=1,i<=Length[decpar2],i++,
	(* res1 and res2: alpha_s *)
	(* res3 and res4: masses  *)
	res1 = AlphasExact[asini,muini,decpar2[[i]][[3]],
			   decpar2[[i]][[1]]-1,loops];
	res3 = mMS2mMS[mqini,asini,res1,decpar2[[i]][[1]]-1,loops];
	res2 = DecAsUpOS[res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			 decpar2[[i]][[1]]-1,loops];
	res4 = DecMqUpOS[res3,res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			 decpar2[[i]][[1]]-1,loops];
	asini = res2;
	mqini = res4;
	muini = decpar2[[i]][[3]];
    ];
    res5 = AlphasExact[asini,muini,mu2,decpar2[[i-1]][[1]],loops];
    res6 = mMS2mMS[mqini,asini,res5,decpar2[[i-1]][[1]],loops];
    Return[res6];
];

(* ************************************************************ *)

(* 
 example: running-decoupling-running-decoupling-running-...
 input: m_q^(h)(mu1), ...
 output: m_q^(l)(mu2)
*)

mH2mL[mqh_,ash_,mu1_,decpar_,mu2_,loops_] := Module[
    {i,asini,muini,res1,res2,res3,res4,res5,res6,decpar2},
    asini=ash;
    muini=mu1;
    mqini=mqh;
    decpar2 = Reverse[Sort[decpar]];
    For[i=2,i<=Length[decpar2],i++,
	If[(decpar2[[i]][[1]]-decpar2[[i-1]][[1]]=!=-1) ||
	   (decpar2[[i]][[1]]-decpar2[[i-1]][[1]]==0) ,
	   Print["WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT."];
	   Abort[];
       ];
    ];
    For[i=1,i<=Length[decpar2],i++,
	(* res1 and res2: alpha_s *)
	(* res3 and res4: masses  *)
	res1 = AlphasExact[asini,muini,decpar2[[i]][[3]],
			   decpar2[[i]][[1]],loops];
	res3 = mMS2mMS[mqini,asini,res1,decpar2[[i]][[1]],loops];
	res2 = DecAsDownOS[res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			   decpar2[[i]][[1]]-1,loops];
	res4 = DecMqDownOS[res3,res1,decpar2[[i]][[2]],decpar2[[i]][[3]],
			   decpar2[[i]][[1]]-1,loops];
	asini = res2;
	mqini = res4;
	muini = decpar2[[i]][[3]];
    ];
    res5 = AlphasExact[asini,muini,mu2,decpar2[[i-1]][[1]]-1,loops];
    res6 = mMS2mMS[mqini,asini,res5,decpar2[[i-1]][[1]]-1,loops];
    Return[res6];
];

(* ************************************************************ *)

(*
 AsRunDec[als, mu0, mu, l] computes \alpha_s^(m)(mu) from the knowledge
 of \alpha_s^(n)(mu0)  
*)

AsRunDec[als_,mu0_,mu_,loops_] := Module[
    {m,n,kk,decpar},
    If[mu>=Global`Mt/.Global`NumDef,m=6,
       If[mu>=Global`Mb/.Global`NumDef,m=5,
	  If[mu>=Global`Mc/.Global`NumDef,m=4,
	     If[mu<Global`Mc/.Global`NumDef,m=3
		]]]];
    If[mu0>=Global`Mt/.Global`NumDef,n=6,
       If[mu0>=Global`Mb/.Global`NumDef,n=5,
	  If[mu0>=Global`Mc/.Global`NumDef,n=4,
	     If[mu0<Global`Mc/.Global`NumDef,n=3
		]]]];
    If[m==n,
       (* Return[{AlphasExact[als,mu0,mu,n,loops],m,n}]; *)
       Return[AlphasExact[als,mu0,mu,n,loops]];
   ];
    decpar = {};
    If[m>n,
       For[kk=n+1,kk<=m,kk++,
	   If[kk==4,decpar=Join[decpar,{{4,Global`Mc/.Global`NumDef,
					 Global`Mc/.Global`NumDef}}]];
	   If[kk==5,decpar=Join[decpar,{{5,Global`Mb/.Global`NumDef,
					 Global`Mb/.Global`NumDef}}]];
	   If[kk==6,decpar=Join[decpar,{{6,Global`Mt/.Global`NumDef,
					 Global`Mt/.Global`NumDef}}]];
       ];
       (* Return[{AlL2AlH[als,mu0,decpar,mu,loops],m,n}]; *)
       Return[AlL2AlH[als,mu0,decpar,mu,loops]];
   ];
    If[m<n,
       For[kk=n,kk>=m+1,kk--,
	   If[kk==4,decpar=Join[decpar,{{4,Global`Mc/.Global`NumDef,
					 Global`Mc/.Global`NumDef}}]];
	   If[kk==5,decpar=Join[decpar,{{5,Global`Mb/.Global`NumDef,
					 Global`Mb/.Global`NumDef}}]];
	   If[kk==6,decpar=Join[decpar,{{6,Global`Mt/.Global`NumDef,
					 Global`Mt/.Global`NumDef}}]];
       ];
       (* Print[decpar]; *)
       (* Return[{AlH2AlL[als,mu0,decpar,mu,loops],m,n}]; *)
       Return[AlH2AlL[als,mu0,decpar,mu,loops]];
   ];
];

(* ************************************************************ *)

(*
 mOS2mMS[mOS, nf, l] computes MS-bar mass corresponding to mOS
*)

mOS2mMS[mOS_,nnf_,loops_] := Module[
    {as},
    as=AsRunDec[Global`asMz/.Global`NumDef,Global`Mz/.Global`NumDef,mOS,loops];
    Return[mOS2mMS[mOS,{},as,mOS,nnf,loops]];
];

(* ************************************************************ *)

(*
 mMS2mOS[mOS, nf, l] computes the on-shell mass corresponding to mMS
*)

mMS2mOS[mMS_,nnf_,loops_] := Module[
    {as},
    as=AsRunDec[Global`asMz/.Global`NumDef,Global`Mz/.Global`NumDef,mMS,loops];
    Return[mMS2mOS[mMS,{},as,mMS,nnf,loops]];
];

(* ************************************************************ *)

(* 
 input: \mu_c^(4)
 output: m_c(MZ)^(5)
*)

Mc5Mzfrommuc4[asMz_,muc4_,Mb_,mub_,Mz_,loops_] := Module[
    {alsmuth,alsmuthp,alsmuc,mcthp,mcth,mcMZ},

    alsmuth = AlphasExact[asMz,Mz,mub,5,loops];
    alsmuthp = DecAsDownOS[alsmuth,Mb,mub,4,loops];
    alsmuc = AlphasExact[alsmuthp,mub,muc4,4,loops];

    mcthp = mMS2mMS[muc4,alsmuc,alsmuthp,4,loops];
    mcth = DecMqUpOS[mcthp,alsmuthp,Mb,mub,4,loops];
    mcMZ = mMS2mMS[mcth,alsmuth,asMz,5,loops];
    
    Return[mcMZ];
];

(* ************************************************************ *)

(* 
 compare running of "exact" vs "\lambda" 
 input: \alpha_s^(5)(Mz)
 output: \alpha_s^(5)(Mb)
 rem: if no value is given for Mb the default one is chosen
*)

AsMbfromAsMz[asMz_,Mb_,loops_] := Module[
    {res1,res2},
    res1 = AlphasExact[asMz,Global`Mz/.Global`NumDef,Mb,5,loops];
    res2 = AlphasLam[LamImpl[asMz,Global`Mz/.Global`NumDef,5,loops],Mb,5,loops];
    Return[ {res1, res2}];
];

AsMbfromAsMz[asMz_,loops_] := AsMbfromAsMz[asMz,Global`Mb/.Global`NumDef,
					   loops];

(* ************************************************************ *)

(* 
 running-decoupling-running
 input: \alpha_s^(5)(Mz)
 output: \alpha_s^(4)(Mc)
*)

As4McfromAs5Mz[asMz_,Mb_,mub_,Mc_,loops_] := Module[
    {res1,res2,res3},
    res1 = AlphasExact[asMz,Global`Mz/.Global`NumDef,mub,5,loops];
    res2 = DecAsDownOS[res1,Mb,mub,4,loops];
    res3 = AlphasExact[res2,mub,Mc,4,loops];
    Return[res3];
];

(* ************************************************************ *)

End[];

(* ************************************************************ *)

EndPackage[];

(* ************************************************************ *)
(* ************************************************************ *)

