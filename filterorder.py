import numpy as np
from scipy import signal as sig

## Functions ##
def get_lp_filter_order(Ad, As, fg, fs):

    ws = fs/fg;
    wg = fg/fg;

    n_bw = butterworth_order(Ad, As, wg, ws);
    n_ch = tschebyscheff_order(Ad, As, wg, ws);
    n_be = bessel_order(3, As, wg, ws);
    n_ca, Wn = cauer_order(Ad, As, wg, ws);

    return n_bw, n_ch, n_be, n_ca

def get_bp_filter_order(Ad, As, dfd, dfs):

    ws = dfs;
    wg = dfd;

    n_bw = butterworth_order(Ad, As, wg, ws);
    n_bw = 2*n_bw;
    n_ch = tschebyscheff_order(Ad, As, wg, ws);
    n_ch = 2*n_ch;
    n_be = bessel_order(Ad, As, wg, ws);
    n_be = 2*n_be;
    n_ca, Wn = cauer_order(Ad, As, wg, ws);
    n_ca = 2*n_ca;

    return n_bw, n_ch, n_be, n_ca


def get_hp_filter_order(Ad, As, fg, fs):
    
    fg2 = fg**2
    ws = fg2/fs;
    wg = fg2/fg;
    # ws_hp = 1/ws;
    # wg_hp = 1/wg;

    n_bw = butterworth_order(Ad, As, wg, ws);
    n_ch = tschebyscheff_order(Ad, As, wg, ws);
    n_be = bessel_order(Ad, As, wg, ws);
    n_ca, Wn = cauer_order(Ad, As, wg, ws);

    return n_bw, n_ch, n_be, n_ca

def butterworth_order(Ad, As, wg, ws):
    ks = 10**(0.1*As)-1;
    kd = 10**(0.1*Ad)-1;

    num = np.log10(ks/kd);
    den = 2*np.log10(ws/wg);

    n_bw = np.ceil(num/den);
    
    return int(n_bw)

def tschebyscheff_order(Ad, As, wg, ws):

    ks = 10**(0.1*As)-1;
    kd = 10**(0.1*Ad)-1;

    num = np.arccosh(np.sqrt(ks/kd));
    den = np.arccosh(ws/wg);

    n_ch = np.ceil(num/den);

    return int(n_ch)

def bessel_order(Ad, As, wg, ws):

    NMAX = 25;
    A0 = np.arange(1, NMAX+2,dtype=np.float64);
    print
    i = NMAX
    while i >= 0:
        if i%2 == 0:
            A0[i] = 1
        else:
            A0[i] = Ad
        i=i-1

    B0 = bessel_polynom(0,wg);
    BS = bessel_polynom(ws,wg);
    H  = A0*(B0/BS);
    Habs = np.absolute(H);
    A = -20*np.log10(Habs);

    if A[NMAX] < As:
        x = 0; # --> Error
    else:
        i = NMAX
        while A[i] > As:
            x = i;
            i = i-1;
    n_be = x;

    return n_be
 
def cauer_order(Ad, As, wg, ws):
    n_ca = sig.ellipord(wg,ws,Ad,As, True)

    return n_ca

def bessel_polynom(w, wg):
    NMAX = 25;
    N = np.arange(1, NMAX+2);
    B = np.arange(1, NMAX+2,dtype=np.float64);

    # Berechne Gruppenlaufzeit
    wgtg0 = np.sqrt( ((2*N) -1)*np.log(2) );
    tg0   = wgtg0/wg;

    S = w*tg0;

    B[0]= 1;
    B[1]= 1+S[0];
    B[2]= 3+3*S[1]+(S[1]**2);

    i = 3;
    while i <= NMAX:
        B[i] = (2*i-1)*B[i-1]+(S[i]**2)*B[i-2];
        i=i+1;
    
    return B
