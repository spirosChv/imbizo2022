TITLE Sodium and Potassium channels

COMMENT
All the channels are taken from same good old classic articles.
The arrengment was done after:
Kang, S., Kitano, K., and Fukai, T. (2004). Self-organized two-state membrane potential transitions in a network of realistically modeled cortical neurons. Neural Netw 17, 307-312. doi: https://doi.org/10.1016/j.neunet.2003.11.010

Whenever available I used the same parameters they used, except in n gate:
n' = phi*(ninf-n)/ntau

Kang used phi = 12, I used phi = 1

Written by Albert Gidon & Leora Menhaim (2004).
ENDCOMMENT

NEURON {
  SUFFIX traub
  NONSPECIFIC_CURRENT i
  RANGE iL,iNa,iK
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
}

STATE {
    m h n a b
}

ASSIGNED {
    v (mV)
    i (mA/cm2)
    iL (mA/cm2)
    iNa (mA/cm2)
    iK (mA/cm2)
    gNa (S/cm2)
    gK (S/cm2)
    minf
    hinf
    ninf 
    mtau (ms)
    htau (ms)
    ntau (ms)
}


BREAKPOINT {
    SOLVE states METHOD cnexp 
    :-------------------------
    :Traub et. al. 1991
    gNa = gNabar*h*pow(m, 2)
    iNa = gNa*(v - eNa)
    gK = gKbar*n
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


FUNCTION vtrap(x (mV), y (mV)) (1) {
    :Traps for 0 in denominator of rate eqns. Taylor expansion is used.
    if (fabs(x/y) < 1e-6) {
        vtrap = 1(/mV)*y*(1 - x/y/2)
    } else {  
        vtrap = 1(/mV)*x/(exp(x/y) - 1)
    }
}

PROCEDURE rates(v (mV)) {
    :Computes rate and other constants at current v.
    :Call once from HOC to initialize inf at resting v.
    
    LOCAL Q10, alpha, beta, vt, qt
    : see Resources/The unreliable Q10.htm for details
    : remember that not only Q10 is temprature dependent 
    : and just astimated here, but also the calculation of
    : Q is itself acurate only in about 10% in this range of
    : temperatures. the transformation formulation is:
    : Q = Q10^((new(degC) - from_original_experiment(degC))/10)

    :--------------------------------------------------------

    : This part was taken **directly** from:
    : Traub, R. D., Wong, R. K., Miles, R., and Michelson, H. (1991).
    : A model of a CA3 hippocampal pyramidal neuron incorporating voltage-clamp data on intrinsic conductances.
    : J Neurophysiol 66, 635-650.


    : Experiments were done in >=32degC for m, and h state variables.
    : Traub et al uses their -60mV as 0mV thus here is the shift.
    vt = v + v_shft :49.2
    qt = Q10^((35(degC) - 32(degC))/10(degC))
    
    : "m" sodium activation system
    alpha = 0.32(/ms)*vtrap(-(vt - 13.1(mV)), 4(mV))
    beta = 0.28(/ms)*vtrap(vt - 40.1(mV), 5(mV))

    mtau = 1/(alpha + beta)
    mtau = mtau/qt
    minf = alpha/(alpha + beta)

    : "h" sodium inactivation system
    alpha = 0.128(/ms)*exp(-(vt - 17(mV))/18(mV))
    beta = 4(/ms)/(1 + exp(-(vt - 40(mV))/5(mV)))
    
    htau = 1/(alpha + beta)
    htau = htau/qt
    hinf = alpha/(alpha + beta)

    : "n" potassium activation system
    alpha = 0.016(/ms)*vtrap(-(vt - 35.1(mV)), 5(mV))
    beta = 0.25(/ms)*exp(-(vt - 20(mV))/40(mV))
    
    ntau = 1/(alpha + beta)
    ntau = ntau/qt
    ninf = alpha/(alpha + beta)
}