import numpy as np
from scipy import signal as sig

## Funktionen ##
# Errechne Filterordnungen von TP-Filtern
def get_lp_filter_order(Ad, As, fg, fs):
    # Normierung der Frequenzen
    ws = fs/fg;
    wg = fg/fg;
    
    # Errechne Ordnung fuer verschiedene Approximationen
    n_bw = butterworth_order(Ad, As, wg, ws);
    n_ch = tschebyscheff_order(Ad, As, wg, ws);
    n_be = bessel_order(3, As, wg, ws);
    n_ca, Wn = cauer_order(Ad, As, wg, ws);

    return n_bw, n_ch, n_be, n_ca

# Errechne Filterordnungen von BP-Filtern
def get_bp_filter_order(Ad, As, dfd, dfs):
    # Frequenzwandlung Bandpass -> TP-Prototyp
    ws = dfs;
    wg = dfd;

    # Errechne Ordnung fuer verschiedene Approximationen (Doppelte Ordnung des TP-Prototypen)
    n_bw = butterworth_order(Ad, As, wg, ws);
    n_bw = 2*n_bw; 
    n_ch = tschebyscheff_order(Ad, As, wg, ws);
    n_ch = 2*n_ch;
    n_be = bessel_order(Ad, As, wg, ws);
    n_be = 2*n_be;
    n_ca, Wn = cauer_order(Ad, As, wg, ws);
    n_ca = 2*n_ca;

    return n_bw, n_ch, n_be, n_ca

# Errechne Filterordnungen von HP-Filtern
def get_hp_filter_order(Ad, As, fg, fs):
    # Frequenzwandlung Hochpass -> TP-Prototyp
    fg2 = fg**2
    ws = fg2/fs;
    wg = fg2/fg;
    # ws_hp = 1/ws;
    # wg_hp = 1/wg;

    # Errechne Ordnung fuer verschiedene Approximationen (Doppelte Ordnung des TP-Prototypen)
    n_bw = butterworth_order(Ad, As, wg, ws);
    n_ch = tschebyscheff_order(Ad, As, wg, ws);
    n_be = bessel_order(Ad, As, wg, ws);
    n_ca, Wn = cauer_order(Ad, As, wg, ws);

    return n_bw, n_ch, n_be, n_ca

# Errechne Filterordnungen von Butterworth-TP
def butterworth_order(Ad, As, wg, ws):
    ks = 10**(0.1*As)-1;
    kd = 10**(0.1*Ad)-1;

    num = np.log10(ks/kd);
    den = 2*np.log10(ws/wg);

    n_bw = np.ceil(num/den);
    
    return int(n_bw)

# Errechne Filterordnungen von Tschebyscheff-TP
def tschebyscheff_order(Ad, As, wg, ws):

    ks = 10**(0.1*As)-1;
    kd = 10**(0.1*Ad)-1;

    num = np.arccosh(np.sqrt(ks/kd));
    den = np.arccosh(ws/wg);

    n_ch = np.ceil(num/den);

    return int(n_ch)

# Errechne Filterordnungen von Bessel-TP
def bessel_order(Ad, As, wg, ws):

    NMAX = 25; # maximal errechenbare Ordnung
    A0 = np.arange(1, NMAX+2,dtype=np.float64);
    i = NMAX
    while i >= 0:
        if i%2 == 0:
            A0[i] = 1 #A0 fuer gerade Filterordnungen -> 1
        else:
            A0[i] = Ad #A0 fuer ungerade Filterordnungen -> Ad
        i=i-1
    # Errechne Besselpolynome fuer w=0
    B0 = bessel_polynom(0,wg);
    # Errechne Besselpolynome fuer w=ws
    BS = bessel_polynom(ws,wg);
    # Errechne Uebertragungsfunktionen fuer ws und wg fuer alle Filterordnungen
    H  = A0*(B0/BS);
    # Errechne Betrag der Uebertragungsfunktion
    Habs = np.absolute(H);
    # Betraf der Uebertragungsfunktion in dB
    A = -20*np.log10(Habs);

    if A[NMAX] < As:
    # Falls As größer als Dämpfung der maximalen Filterordnung bei fs
        x = 0; # --> Error
    else:
        i = NMAX
        while A[i] > As:
        # Falls As größer als Dämpfung der aktuellen Filterordnung bei fs -> Filterordnungen = i
            x = i;
            i = i-1; # Dekrementiere i
    n_be = x;

    return n_be
 
# Errechne Filterordnungen von Cauer-TP
def cauer_order(Ad, As, wg, ws):
    n_ca = sig.ellipord(wg,ws,Ad,As, True)

    return n_ca

# Berechne Besselpolynom
def bessel_polynom(w, wg):
    NMAX = 25; # maximal errechenbare Ordnung
    N = np.arange(1, NMAX+2);                   # Array fuer Filterordnungen
    B = np.arange(1, NMAX+2,dtype=np.float64);  # Array fuer Besselpolynome

    # Berechne Gruppenlaufzeit(tg0)*Grenzfrequenz(wg)
    wgtg0 = np.sqrt( ((2*N) -1)*np.log(2) );
    # Berechne Gruppenlaufzeit(tg0)
    tg0   = wgtg0/wg;
    
    # Multipliziere Frequenz fuer die Besselpolynom berechnet werden soll mit tg0
    S = w*tg0;
    
    # Berechne Besselpolynome bis Grad 2
    B[0]= 1;
    B[1]= 1+S[0];
    B[2]= 3+3*S[1]+(S[1]**2);

    # Berechne Besselpolynome von Grad 3 bis NMAX rekursiv
    i = 3;
    while i <= NMAX:
        B[i] = (2*i-1)*B[i-1]+(S[i]**2)*B[i-2];
        i=i+1;
    
    return B
