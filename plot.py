# import libraries
import numpy as np
from scipy import signal as sig
import matplotlib.pyplot as plt
import matplotlib as mpl

## Plotte Uebertragungsfunktion
# Eingabe:
#   pol     Polstellen
#   null    Nullstellen
#   title   Titel für Plot
def plot(pol,null,title):
    # errechne k so, dass das maximum der 
    # Uebertragungsfunktion = 1 ist
    w_1, h_1 = sig.freqs_zpk(null, pol, 1);
    k = 1/(max(abs(h_1)));
    # Erstelle Array für Frequenzachse
    f=np.linspace(1,100000,num=100000)

    # Berechne die Uebertragungsfunktion
    # w, h = sig.freqs_zpk(null, pol, k, f)
    b, a = sig.zpk2tf(null, pol, k)
    w, h = sig.freqs(b, a, f)
    # Berechne den Betrag der Uebertragungsfunktion
    hdB = 20* np.log10(abs(h));
    # Berechne die Phase der Uebertragungsfunktion
    phase_rad = np.angle(h);
    phase_deg = np.rad2deg(np.angle(h));
    # Berechne die Phase der Uebertragungsfunktion
    # tg0 = -np.diff(phase_rad, prepend=phase_rad[0]); #alt: Fehler bei Phasendrehung von 180° auf -180°
    # Skalierung auf ms
    # tg0[0] = tg0[1]
    wg,tg0 = sig.group_delay((b,a),100000 ,fs=20000, whole=True) # neue Berechnung der Gruppenlaufzeit
    
    ## PLOT ##
    fig, ax = plt.subplots(num=title);
    # Erzeuge zwei zusaetzliche Ordinaten
    fig.subplots_adjust(right=0.75);
    twin1 = ax.twinx() # Phase
    twin2 = ax.twinx() # Gruppenlaufzeit
    twin2.spines.right.set_position(("axes", 1.2))
    # Erzeuge Grid
    ax.grid(True)
    plt.title(title)

    # plotte den Betrag der Uebertragungsfunktion
    p1,= ax.semilogx(w, hdB, color='blue', label=r'$|H(j\omega)|$')
    ax.set_xlabel('Frequenz in Hz')         # Beschriftung Abzisse
    ax.set_ylabel(r'$|H(j\omega)|$ in dB')  # Beschriftung Ordinate 1

    p2,= twin1.semilogx(w, phase_deg, color='orange', label=r'$arg(H(j\omega))$')
    twin1.set_ylabel(r'$arg[ H(j\omega) ]$ in $\degree$') # Beschriftung Ordinate 2

    # plotte die Gruppenlaufzeit der Uebertragungsfunktion
    # p3,= twin2.semilogx(wg, tg0, color='green', label='Gruppenlaufzeit in s')
    p3,= twin2.semilogx(w, tg0, color='green', label='Gruppenlaufzeit in ms') 
    twin2.set_ylabel('Gruppenlaufzeit in ms') # Beschriftung Ordinate 3
    
    # Farbe der Label stimmt mit Graph ueberein 
    ax.yaxis.label.set_color(p1.get_color())
    twin1.yaxis.label.set_color(p2.get_color())
    twin2.yaxis.label.set_color(p3.get_color())
    
    # Farbe der Ticks stimmt mit Graph ueberein
    ax.tick_params(axis='y', colors=p1.get_color())
    twin1.tick_params(axis='y', colors=p2.get_color())
    twin2.tick_params(axis='y', colors=p3.get_color())
    ax.tick_params(axis='x')
    # Legende
    ax.legend(handles=[p1, p2, p3])
    # Plot blockt keine neuen plots
    # plt.ion()
    # zeige Plot an
    plt.show()
