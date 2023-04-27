import numpy as np

## Berechne Widerstandswerte in kOhm
# Eingabe: 
#   N               Filterordnung
#   Filtertyp       Bandpass oder Tiefpass
#   f0              
#   Q
#   A               Gesamtverstärkung des Filters
# Ausgabe:
#   R1,R2,R3,R4     Arrays mit Widerstandswerten [kOhm]
#   H0              Array mit Verstärkung je Sektion
#   fc              Array mit Potential, das mit FC verbunden wird 
def berechne_widerstand(N, Filtertyp, f0, Q, A):
    ## Variablen Deklarieren ##
    faktor  = 2*(10**9);# konstanter Faktor fuer Widerstandsberechnung
    r5k     = 5000;     # 5000 kOhm
    r4M     = 4000000;  # 4 MOhm
    RxRyG   = 0.2;      # Rx/Ry, wenn FC mit GND verbunden
    RxRyP   = 4;        # Rx/Ry, wenn FC mit V+ verbunden
    RxRyM   = 0.04;     # Rx/Ry, wenn FC mit V- verbunden
    numSek  = N;        # Anzahl an Filtersektionen (entspricht N)

    # Erstelle Arrays fuer Widerstandswerte und FC
    R1 = np.arange(0,numSek);
    R2 = np.arange(0,numSek);
    R3 = np.arange(0,numSek);
    R4 = np.arange(0,numSek);
    fc = np.empty(( numSek,1),dtype=object);
    
    ## Berechnungen ##
    # Berechne Verstärkung fuer jede Sektion
    H0 = np.tile(A**(1/N), numSek); # H0 = A^(1/N)
    # Berechne Widerstaende
    R2 = faktor/f0; # R2 = (2*10^9) / f0
    R4 = R2 - r5k;  # R2 = R2 - 5kOhm

    R3_raw = (Q*faktor)/f0; # R3 = (Q*2*10^9)/ f0 * Rx/Ry
    R3G = R3_raw * RxRyG;   # R3, wenn FC mit GND verbunden
    R3P = R3_raw * RxRyP;   # R3, wenn FC mit V+ verbunden
    R3M = R3_raw * RxRyM;   # R3, wenn FC mit V- verbunden

    if Filtertyp == 'Tiefpass':
        #Bei Tiefpass
        R1_raw = faktor/(f0*H0);    # R1 = (2*10^9)/(f0*H0) * Rx/Ry
        R1G = R1_raw * RxRyG;       # R1, wenn FC mit GND verbunden
        R1P = R1_raw * RxRyP;       # R1, wenn FC mit V+ verbunden
        R1M = R1_raw * RxRyM;       # R1, wenn FC mit V- verbunden
    else:
        #Bei Bandpass
        R1G = R3G / H0;             # R1 = R3 / H0
        R1P = R3P / H0;
        R1M = R3M / H0;
    
    # Ermittle korrektes Potential, mit dem FC verbunden wird und ordne korrekten Widerstandswert zu
    i = 0
    while i < numSek:
        if R3G[i] < r5k:
            if (R3G[i] > r4M) or (R3P[i] > r4M):
                R3[i] = R3M[i];
                R1[i] = R1M[i];
                fc[i] = 'V-';
            else:
                R3[i] = R3P[i]
                R1[i] = R1P[i]
                fc[i] = 'V+'
        else:
            R3[i] = R3G[i];
            R1[i] = R1G[i];
            fc[i] = 'GND';
        i = i+1

    # Rückgabe in kOhm
    return R1/1000, R2/1000, R3/1000, R4/1000, H0, fc

## Berechne Polstellen aus Widerstandswerten
# Eingabe:
#   R1k,R2k,R3k,R4k Arrays mit Widerstandswerten [kOhm]
#   fc              Array mit Potential, das mit FC verbunden wird 
#   Filtertyp       Bandpass oder Tiefpass
# Ausgabe:
#   f0 
#   Q 
#   p               Polstellen
#   z               Nullstellen
def berechne_pole(R1k,R2k,R3k,R4k,fc, Filtertyp):
    ## Variablen Deklarieren ##
    R1      = R1k*1000  # R1 in Ohm 
    R2      = R2k*1000  # R2 in Ohm
    R3      = R3k*1000  # R3 in Ohm
    R4      = R4k*1000  # R4 in Ohm
    numSek  = np.size(R1)   # Anzahl an Filtersektionen
    numSek2 = 2*numSek      # 2* Anzahl an Filtersektionen
    p       = np.zeros(numSek2, dtype='complex_') # Array für Polstellen
    z       = np.zeros(numSek, dtype='complex_')  # Array für Nullstellen
    faktor  = 2*(10**9) # konstanter Faktor 2*10^9
    r5k     = 5000      # 5kOhm
    r4M     = 4000000   # 4MOhm
    root    = np.sqrt(1/(R2*(R4+r5k))) # sqrt( 1/(R2*(R4+5kOhm)) )
    RxRy    = 65/13     # Rx/Ry
    
    ## Berechnungen ##
    # Bestimme korrekten Wert für Rx/Ry je nach Potential an FC
    i=0
    while i < numSek:
        if fc[i] == 'GND':
            RxRy=65/13
        elif fc[i] == 'V+':
            RxRy=13/52
        elif fc[i] == 'V-':
            RxRy=325/13
        i = i+1

    f0 = root * faktor      # f0 = (2*10^9)* sqrt( 1/(R2*(R4+5kOhm)) )
    Q  = root * R3*RxRy     # Q  = R3* Rx/Ry * sqrt( 1/(R2*(R4+5kOhm)) )
    real = -f0 / (2*Q)      # Berechne Realteil der Polstellen: re = -f0 / (2*Q) 
    imag = np.sqrt( (f0**2)-(real**2) ) # Berechne Imaginärteil der Polstellen = sqrt(f0^2-re^2)

    k = 0
    while k < numSek:
        # Erstelle Array mit Polstellen p = re + im
        p[k]                = real[k] + 1j*imag[k] 
        p[numSek2-k-1]      = real[k] - 1j*imag[k]
        if Filtertyp == 'Bandpass':
            # Wenn Bandpass -> Nullstellen bei 0
            z[k]= 0
        k = k+1

    if Filtertyp == 'Tiefpass':
        # Wenn Tiefpass keine Nullstellen
        z = np.empty(0)
    # Rückgabe
    return f0, Q, p, z


def f0Q_to_PZ(f0, Q, Typ):
    P = np.zeros(2,dtype='complex_')
    real = -f0 / (2*Q)
    imag = np.sqrt( (f0**2)-(real**2) )
    P[0] = real + 1j*imag
    P[1] = real - 1j*imag
    if Typ == 'Bandpass':
        Z = np.zeros(2,dtype='complex_')
        Z[0] = 0
        Z[1] = 0
    else:
        Z = np.empty(0)
    return P, Z

