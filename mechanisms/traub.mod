COMMENT
All the channels are taken from same good old classic articles.
The arrengment was done after:
Kang, S., Kitano, K., and Fukai, T. (2004). 
  Self-organized two-state membrane potential 
  transitions in a network of realistically modeled 
  cortical neurons. Neural Netw 17, 307-312.

Whenever available I used the same parameters they used,
except in n gate:
  n' = phi*(ninf-n)/ntau

Kang used phi = 12
I used phi = 1

Written by Albert Gidon & Leora Menhaim (2004).
ENDCOMMENT

NEURON {
  SUFFIX traub
  NONSPECIFIC_CURRENT i
  RANGE iL, iNa, iK
  RANGE eL, eNa, eK
  RANGE gLbar, gNabar, gKbar
  RANGE v_shft
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
}

PARAMETER {
    gNabar = .03 (S/cm2)    :Traub et. al. 1991
    gKbar = .015 (S/cm2)    :Traub et. al. 1991
    gLbar = 0.00014 (S/cm2) :Siu Kang - by email.
    eL = -62.0 (mV) :Siu Kang - by email.
    eK = -80 (mV)   :Siu Kang - by email.
    eNa = 90 (mV)   :Leora
    v_shft = 49.2 (mV) : shift to apply to all curves
    Q10 = 3 (1)  : temperature sensitivity
}

STATE {
    m
    h
    n
}

ASSIGNED {
    v (mV)
    i (mA/cm2)
    iL (mA/cm2)
    iNa (mA/cm2)
    iK (mA/cm2)
    gNa (S/cm2)
    gK (S/cm2)
    minf (1)
    hinf (1)
    ninf (1) 
    mtau (ms)
    htau (ms)
    ntau (ms) 
}


BREAKPOINT {
    SOLVE states METHOD cnexp 
    :-------------------------
    :Traub et. al. 1991
    gNa = gNabar*h*m*m
    iNa = gNa*(v - eNa)
    gK = gKbar*n : - Traub et. al. 1991
    iK = gK*(v - eK)
    :-------------------------
    iL = gLbar*(v - eL)
    :-------------------------
    i = iL + iK + iNa
}

INITIAL {
    rates(v)
    m = minf
    h = hinf
    n = ninf
}

DERIVATIVE states {  
    rates(v)
    :Traub Spiking channels
    m' = (minf-m)/mtau
    h' = (hinf-h)/htau
    n' = 2*(ninf-n)/ntau :phi=12 from Kang et. al. 2004
}

PROCEDURE rates(v(mV)) {  
  :Computes rate and other constants at current v.
  :Call once from HOC to initialize inf at resting v.
  LOCAL  alpha, beta, sum, vt, qt
  : see Resources/The unreliable Q10.htm for details
  : remember that not only Q10 is temprature dependent 
  : and just astimated here, but also the calculation of
  : Q is itself acurate only in about 10% in this range of
  : temperatures. the transformation formulation is:
  : Q = Q10^(( new(degC) - from_original_experiment(degC) )/ 10)
  
  :--------------------------------------------------------
  
  : This part was taken **directly** from:
  : Traub, R. D., Wong, R. K., Miles, R., and Michelson, H. (1991). 
  : A model of a CA3 hippocampal pyramidal neuron incorporating 
  : voltage-clamp data on intrinsic conductances. 
  : J Neurophysiol 66, 635-650.
  : Experiments were done in >=32degC for m,h
  : Traub et al uses their -60mV as 0mV thus here is the shift
  vt = v + v_shft :49.2
  qt = Q10^((35 - 32)/ 10)
  :"m" sodium activation system
  if(vt == 13.1(mV)){alpha = 0.32*4(/ms)}
  else{alpha = 0.32(/ms/mV)*(13.1(mV) - vt)/(exp((13.1(mV) - vt)/4(mV)) - 1)}
  if(vt == 40.1){beta = 0.28*5}
  else{beta = 0.28(/ms/mV)*(vt - 40.1(mV))/(exp((vt - 40.1(mV))/5(mV))-1)}
  sum = alpha + beta
  mtau = 1/sum
  mtau = mtau/qt
  minf = alpha/sum

  :"h" sodium inactivation system
  alpha = 0.128(/ms)*exp((17(mV) - vt)/18(mV))
  beta = 4(/ms)/(1 + exp((40(mV) - vt)/5(mV)))
  sum = alpha + beta
  htau = 1/sum
  htau = htau/qt
  hinf = alpha/sum

  :"n" potassium activation system
  if(vt == 35.1(mV)){ alpha = 0.016*5 }
  else{alpha = 0.016(/mV/ms)*(35.1(mV) - vt)/(exp((35.1(mV) - vt)/5(mV)) - 1)}
  beta = 0.25(/ms)*exp((20(mV) - vt)/40(mV))
  sum = alpha + beta
  ntau = 1/sum
  ntau = ntau/qt
  ninf = alpha/sum
}