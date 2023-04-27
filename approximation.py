import numpy as np
from scipy import signal as sig
import filterorder as ord

def lp_approximations(Ad,As,fg,fs,N,app):
    # Normierung von fg und fs
    wg = 1; # fg/fg = 1
    ws = fs/fg;

    # Überprüfe, welche Approximation genutzt wird
    if app == 'Butterworth':
            p, z = butterworth_pol_null(Ad,wg,N);
    elif app == 'Tschebyscheff':
            p, z = tschebyscheff_pol_null(Ad,wg,N);
    elif app == 'Bessel':
            p, z = bessel_pol_null(wg,N);
    elif app == 'Cauer':
            p, z = cauer_pol_null(Ad,As,wg,N);
    else:
        # Fehlermeldung
        print('invalid approximation type, choose bw, ch, be, ca')
    # Entnormierung von Pol- und Nullstellen
    p = p*fg;
    z = z*fg;
    
    # Ermittle f0 (p_abs)   f0 = |p|
    p_abs = np.absolute(p); 
    # Ermittle Q            Q  = -|p|/(2*Re{p})
    Q = -p_abs / ( 2*np.real(p) );
    return p, z, p_abs, Q

def hp_approximations(Ad,As,fg,fs,N,app):
    # Frequenzwandlung HP->TP-Prototyp
    fg2 = fg**2;
    # Normierung von fg und fs
    wg = fg2/fg;
    ws = fg2/fs;

    # Überprüfe, welche Approximation genutzt wird
    if app == 'Butterworth':
            p, z = butterworth_pol_null(Ad,wg,N);
    elif app == 'Tschebyscheff':
            p, z = tschebyscheff_pol_null(Ad,wg,N);
    elif app == 'Bessel':
            p, z = bessel_pol_null(wg,N);
    elif app == 'Cauer':
            p, z = cauer_pol_null(Ad,As,wg,N);
    else:
        # Fehlermeldung
        print('invalid approximation type, choose bw, ch, be, ca')
    # Ermittle Q            Q  = -|p|/(2*Re{p})
    Q = -abs(p) / ( 2*np.real(p) ) ;
    # Frequenzwandlung
    p = fg2/p;

    if app == 'ca':
        z = fg2/z;
    else:
        z = np.zeros(N);

    # Ermittle f0 (p_abs)   f0 = |p|
    p_abs = abs(p);
    return p, z, p_abs, Q

def bp_approximations(Ad,As,dfd,dfs,fc,N,app):
    # Frequenzwandlung BP->TP-Prototyp
    wM2 = fc**2; 
    ws = dfs;
    wg = dfd;
     
    NTP = int(N/2); # Ordnung des TP-Prototypen = N/2

    p = np.arange(N, dtype='complex_')
    z = np.arange(NTP, dtype='complex_')

    # Überprüfe, welche Approximation genutzt wird
    if app == 'Butterworth':
            p1, z1 = butterworth_pol_null(Ad,wg,NTP);
    elif app == 'Tschebyscheff':
            p1, z1 = tschebyscheff_pol_null(Ad,wg,NTP);
    elif app == 'Bessel':
            p1, z1 = bessel_pol_null(wg,NTP);
    elif app == 'Cauer':
            p1, z1 = cauer_pol_null(Ad,As,wg,NTP);
    else:
        # Fehlermeldung
        print('invalid approximation type, choose bw, ch, be, ca')
    # Frequenzwandlung der Polstellen
    P_1 = 0.5*p1
    P_2 = np.sqrt( np.square(P_1)-wM2 )
    if app == 'ca':
        # Wenn Cauer Filter gewählt wurde -> Frequenzwandlung der Nullstellen
        z = np.zeros(int(NTP), dtype='complex_')
        Z_1 = 0.5*z1
        Z_2 = np.sqrt(np.square(Z_1)-wM2 )

    # Sortiere Reihenfolge der Polstellen
    i = 0;
    while i < NTP: 
        p[i]      = P_1[i] + P_2[i]
        p[NTP+i]  = P_1[i] - P_2[i]
        z[i]    = 0
        if app == 'ca':
            # Sortiere Reihenfolge der Nullstellen, falls Cauer Filter gewählt wurde
            z[i] = Z_1[i] - Z_2[i]
        i = i+1
    # Sortiere Polstellen für einheitliche Reihenfolge zwischen BP, HP, TP
    p = np.sort(p)
    # Ermittle f0 (p_abs)   f0 = |p|
    p_abs = np.absolute(p);
    # Ermittle Q            Q  = -|p|/(2*Re{p})
    Q = -p_abs / ( 2*np.real(p) ) ;
    return p, z, p_abs, Q

def get_second_first_order(f0,Q,N):   
    numSek2 = int(np.floor(N/2)) # Ermittle Anzahl der Sektionen 2ten Grades
    numSek1 = int(N%2)           # Ermittle ob Sektion ersten Grades existiert
    f02 = np.zeros(numSek2)      # Erstelle Array für f0 der Sektionen 2ten Grades
    f01 = np.empty(0)            # Erstelle Array für f0 der Sektion 1ten Grades
    Q2  = np.zeros(numSek2)      # Erstelle Array für Q der Sektionen 2ten Grades
    Q1  = np.empty(0)            # Erstelle Array für Q der Sektion 1ten Grades
    
    # Ermittle f0 und Q für Sektionen 2ten Grades
    i=0
    while i < numSek2:
        if (i%2 == 0):
            f02[i] = f0[i]
            Q2[i]  = Q[i]
        else:
            f02[i] = f0[N-i-1]
            Q2[i]  = Q[N-i-1]
        i = i+1
    # Ermittle f0 und Q für Sektionen 1ten Grades
    if numSek1 > 0:
        f01 = np.zeros(numSek1)
        Q1  = np.zeros(numSek1)
        f01[0] = f0[numSek2]
        Q1[0]  = Q[numSek2]

    return f01, Q1, f02, Q2

def butterworth_pol_null(Ad, wg, N):
    # Erzeuge Vektor k[1:N]
    k = np.arange(1,N+1)
    kd = 10**(0.1*Ad)-1;
    epsilon = np.sqrt(kd);

    pol_radius = wg*epsilon**(-1/N);
    x = (np.pi*( k*2 + N-1 ) ) / (N*2);

    real = np.cos(x);
    imag = np.sin(x);

    p=pol_radius*(real+1j*imag);

    # z=np.tile(np.inf, N)
    z=np.empty(0)
    return p, z

def tschebyscheff_pol_null(Ad, wg, N):
    # Erzeuge Vektor k[1:N]
    k = np.arange(1,N+1)
    kd = 10**(0.1*Ad)-1;
    epsilon = np.sqrt(kd);

    pol_radius = wg*epsilon**(-1/N);
    x = (np.pi*( k*2 -1 ) ) / (N*2);

    a = np.sinh( (1/N) * np.arcsinh(1/epsilon) );
    b = np.cosh( (1/N) * np.arcsinh(1/epsilon) );

    real = -a*np.sin(x);
    imag = -b*np.cos(x);

    p=pol_radius*(real+1j*imag);

    # Nullstellen 
    # z=np.tile(np.inf, N);
    z=np.empty(0)
    return p, z

def bessel_pol_null(wg,N):
    # Look-up Tabelle für Polstellen nach Filtergrad
    if N == 1 :
        p =  np.array([-1.+0.j])
    elif N == 2 :
        p =  np.array([-1.10160133+0.63600982j, -1.10160133-0.63600982j])
    elif N == 3 :
        p =  np.array([-1.04740916+0.99926444j, -1.3226758 +0.j        , -1.04740916-0.99926444j])
    elif N == 4 :
        p =  np.array([-0.99520876+1.25710574j, -1.37006783+0.41024972j, -1.37006783-0.41024972j \
                      ,-0.99520876-1.25710574j])
    elif N == 5 :
        p =  np.array([-0.95767655+1.47112432j, -1.38087733+0.71790959j, -1.50231627+0.j         \
                      ,-1.38087733-0.71790959j, -0.95767655-1.47112432j])
    elif N == 6 :
        p =  np.array([-0.93065652+1.66186327j, -1.3818581 +0.97147189j, -1.5714904 +0.32089637j \
                      ,-1.5714904 -0.32089637j, -1.3818581 -0.97147189j, -0.93065652-1.66186327j])
    elif N == 7 :
        p =  np.array([-0.90986778+1.83645135j, -1.37890322+1.19156678j, -1.61203877+0.58924451j \
                      ,-1.68436818+0.j        , -1.61203877-0.58924451j, -1.37890322-1.19156678j \
                      ,-0.90986778-1.83645135j])
    elif N == 8 :
        p =  np.array([-0.89286972+1.99832584j, -1.37384122+1.38835658j, -1.63693942+0.82279563j \
                      ,-1.7574084 +0.27286758j, -1.7574084 -0.27286758j, -1.63693942-0.82279563j \
                      ,-1.37384122-1.38835658j, -0.89286972-1.99832584j])
    elif N == 9 :
        p =  np.array([-0.87839928+2.14980052j, -1.36758831+1.56773371j, -1.65239648+1.03138957j \
                      ,-1.80717053+0.51238373j, -1.8566005 +0.j        , -1.80717053-0.51238373j \
                      ,-1.65239648-1.03138957j, -1.36758831-1.56773371j, -0.87839928-2.14980052j])
    elif N == 10 :
        p =  np.array([-0.8657569 +2.29260483j, -1.36069228+1.73350574j, -1.66181024+1.22110022j \
                      ,-1.84219624+0.7272576j , -1.92761969+0.24162347j, -1.92761969-0.24162347j \
                      ,-1.84219624-0.7272576j , -1.66181024-1.22110022j, -1.36069228-1.73350574j \
                      ,-0.8657569 -2.29260483j])
    elif N == 11 :
        p =  np.array([-0.85451258+2.42805947j, -1.35348668+1.88829684j, -1.66719364+1.3959629j  \
                      ,-1.86736124+0.92311558j, -1.98016065+0.45959874j, -2.01670147+0.j         \
                      ,-1.98016065-0.45959874j, -1.86736124-0.92311558j, -1.66719364-1.3959629j  \
                      ,-1.35348668-1.88829684j, -0.85451258-2.42805947j])
    elif N == 12 :
        p =  np.array([-0.84437887+2.55718897j, -1.34617468+2.03399851j, -1.66980359+1.5588027j  \
                      ,-1.88564962+1.10381488j, -2.01994593+0.6589965j , -2.08464451+0.21916154j \
                      ,-2.08464451-0.21916154j, -2.01994593-0.6589965j , -1.88564962-1.10381488j \
                      ,-1.66980359-1.5588027j , -1.34617468-2.03399851j, -0.84437887-2.55718897j])
    elif N == 13 :
        p =  np.array([-0.83515201+2.6808028j , -1.33888032+2.17202295j, -1.67045856+1.71167829j \
                      ,-1.89898612+1.27211944j, -2.05058088+0.84338311j, -2.13764829+0.42041631j \
                      ,-2.16608271+0.j        , -2.13764829-0.42041631j, -2.05058088-0.84338311j \
                      ,-1.89898612-1.27211944j, -1.67045856-1.71167829j, -1.33888032-2.17202295j \
                      ,-0.83515201-2.6808028j ])
    elif N == 14 :
        p =  np.array([-0.82668133+2.79955221j, -1.3316792 +2.30345534j, -1.66971016+1.85613874j \
                      ,-1.90866458+1.43007973j, -2.07445158+1.01536709j, -2.17970952+0.60702983j \
                      ,-2.23093074+0.20200027j, -2.23093074-0.20200027j, -2.17970952-0.60702983j \
                      ,-2.07445158-1.01536709j, -1.90866458-1.43007973j, -1.66971016-1.85613874j \
                      ,-1.3316792 -2.30345534j, -0.82668133-2.79955221j])
    elif N == 15 :
        p =  np.array([-0.81885183+2.91396993j, -1.32461676+2.42914977j, -1.6679408 +1.99338081j \
                      ,-1.91558435+1.57926036j, -2.09319972+1.17691618j, -2.21352749+0.78142998j \
                      ,-2.28342656+0.38982894j, -2.30637006+0.j        , -2.28342656-0.38982894j \
                      ,-2.21352749-0.78142998j, -2.09319972-1.17691618j, -1.91558435-1.57926036j \
                      ,-1.6679408 -1.99338081j, -1.32461676-2.42914977j, -0.81885183-2.91396993j])
    elif N == 16 :
        p =  np.array([-0.81157347+3.02449808j, -1.31771944+2.54979207j, -1.66542193+2.12434959j \
                      ,-1.92038809+1.72088333j, -2.10799083+1.32955251j, -2.24099322+0.94547405j \
                      ,-2.32647906+0.56573669j, -2.36834668+0.18833296j, -2.36834668-0.18833296j \
                      ,-2.32647906-0.56573669j, -2.24099322-0.94547405j, -2.10799083-1.32955251j \
                      ,-1.92038809-1.72088333j, -1.66542193-2.12434959j, -1.31771944-2.54979207j \
                      ,-0.81157347-3.02449808j])
    elif N == 17 :
        p =  np.array([-0.80477429+3.13150834j, -1.31100147+2.66594249j, -1.66235006+2.24980556j \
                      ,-1.92354583+1.85592213j, -2.11967427+1.47447847j, -2.26347141+1.10061605j \
                      ,-2.3621559 +0.73146474j, -2.41990692+0.36508673j, -2.43892722+0.j         \
                      ,-2.41990692-0.36508673j, -2.3621559 -0.73146474j, -2.26347141-1.10061605j \
                      ,-2.11967427-1.47447847j, -1.92354583-1.85592213j, -1.66235006-2.24980556j \
                      ,-1.31100147-2.66594249j, -0.80477429-3.13150834j])
    elif N == 18 :
        p =  np.array([-0.7983958 +3.23531668j, -1.30446926+2.77806544j, -1.65886985+2.37037051j \
                      ,-1.92540814+1.98516539j, -2.12888269+1.61266094j, -2.28197142+1.24801645j \
                      ,-2.39196895+0.88838996j, -2.46324195+0.53189096j, -2.49830137+0.17711353j \
                      ,-2.49830137-0.17711353j, -2.46324195-0.53189096j, -2.39196895-0.88838996j \
                      ,-2.28197142-1.24801645j, -2.12888269-1.61266094j, -1.92540814-1.98516539j \
                      ,-1.65886985-2.37037051j, -1.30446926-2.77806544j, -0.7983958 -3.23531668j])
    elif N == 19 :
        p =  np.array([-0.79238974+3.33619438j, -1.2981241 +2.88655089j, -1.65508928+2.48655994j \
                      ,-1.926241  +2.10926147j, -2.13609626+1.74488941j, -2.29725552+1.38861665j \
                      ,-2.41704998+1.03762101j, -2.49997144+0.69004482j, -2.54874515+0.34453452j \
                      ,-2.56484699+0.j        , -2.54874515-0.34453452j, -2.49997144-0.69004482j \
                      ,-2.41704998-1.03762101j, -2.29725552-1.38861665j, -2.13609626-1.74488941j \
                      ,-1.926241  -2.10926147j, -1.65508928-2.48655994j, -1.2981241 -2.88655089j \
                      ,-0.79238974-3.33619438j])
    elif N == 20 :
        p =  np.array([-0.78671579+3.4343764j , -1.29196404+2.99172993j, -1.65108992+2.59880653j \
                      ,-1.92624918+2.22875022j, -2.14168539+1.87181719j, -2.30990971+1.52319125j \
                      ,-2.43826313+1.18006442j, -2.53131995+0.84060226j, -2.59195625+0.50348811j \
                      ,-2.62187627+0.16768821j, -2.62187627-0.16768821j, -2.59195625-0.50348811j \
                      ,-2.53131995-0.84060226j, -2.43826313-1.18006442j, -2.30990971-1.52319125j \
                      ,-2.14168539-1.87181719j, -1.92624918-2.22875022j, -1.65108992-2.59880653j \
                      ,-1.29196404-2.99172993j, -0.78671579-3.4343764j ])
    elif N == 21 :
        p =  np.array([-0.78133993+3.53006786j, -1.28598502+3.09388646j, -1.64693399+2.70747729j \
                      ,-1.92559238+2.34408624j, -2.14593988+1.99399131j, -2.32039149+1.65238549j \
                      ,-2.4562796 +1.31647125j, -2.55823111+0.98443092j, -2.62922902+0.65497428j \
                      ,-2.6711425 +0.32710771j, -2.68500354+0.j        , -2.6711425 -0.32710771j \
                      ,-2.62922902-0.65497428j, -2.55823111-0.98443092j, -2.4562796 -1.31647125j \
                      ,-2.32039149-1.65238549j, -2.14593988-1.99399131j, -1.92559238-2.34408624j \
                      ,-1.64693399-2.70747729j, -1.28598502-3.09388646j, -0.78133993-3.53006786j])
    elif N == 22 :
        p =  np.array([-0.77623311+3.62344908j, -1.28018174+3.19326609j, -1.64266938+2.81288661j \
                      ,-1.92439656+2.45565624j, -2.14908932+2.11187457j, -2.32906289+1.7767427j  \
                      ,-2.47162839+1.44747102j, -2.58144409+1.12225408j, -2.66156868+0.79982098j \
                      ,-2.71399201+0.47920723j, -2.73991884+0.15962437j, -2.73991884-0.15962437j \
                      ,-2.71399201-0.47920723j, -2.66156868-0.79982098j, -2.58144409-1.12225408j \
                      ,-2.47162839-1.44747102j, -2.32906289-1.7767427j , -2.14908932-2.11187457j \
                      ,-1.92439656-2.45565624j, -1.64266938-2.81288661j, -1.28018174-3.19326609j \
                      ,-0.77623311-3.62344908j])
    elif N == 23 :
        p =  np.array([-0.77137039+3.7146796j , -1.27454813+3.29008293j, -1.63833315+2.91530602j \
                      ,-1.92276207+2.56379208j, -2.15131761+2.22586219j, -2.33621388+1.89672489j \
                      ,-2.48473196+1.57359672j, -2.60154614+1.25468134j, -2.68976834+0.9387221j  \
                      ,-2.75147737+0.62477604j, -2.78800663+0.31208598j, -2.80010283+0.j         \
                      ,-2.78800663-0.31208598j, -2.75147737-0.62477604j, -2.68976834-0.9387221j  \
                      ,-2.60154614-1.25468134j, -2.48473196-1.57359672j, -2.33621388-1.89672489j \
                      ,-2.15131761-2.22586219j, -1.92276207-2.56379208j, -1.63833315-2.91530602j \
                      ,-1.27454813-3.29008293j, -0.77137039-3.7146796j ])
    elif N == 24 :
        p =  np.array([-0.76673015+3.80390135j, -1.26907773+3.38452499j, -1.63395419+3.0149719j  \
                      ,-1.92076955+2.66878091j, -2.15277351+2.33629451j, -2.34207917+2.01272829j \
                      ,-2.49593149+1.69530376j, -2.61900972+1.38223142j, -2.7144624 +1.07226532j \
                      ,-2.78443426+0.76447858j, -2.83035558+0.45813734j, -2.85310631+0.15262274j \
                      ,-2.85310631-0.15262274j, -2.83035558-0.45813734j, -2.78443426-0.76447858j \
                      ,-2.7144624 -1.07226532j, -2.61900972-1.38223142j, -2.49593149-1.69530376j \
                      ,-2.34207917-2.01272829j, -2.15277351-2.33629451j, -1.92076955-2.66878091j \
                      ,-1.63395419-3.0149719j , -1.26907773-3.38452499j, -0.76673015-3.80390135j])
    else:
        print('ERROR: Filterordnung > 24 oder <1')
        p =np.empty(0)
    # Skalierung der Polstellen mit wg
    p = p*wg
    z=np.empty(0)
    return p, z

def cauer_pol_null(Ad, As, wg, N):
    # Nutze integrierte scipy funktion zum ermitelln der Pole und Nullstellen der Cauer Approximation
    z, p, k = sig.ellipap(N, Ad, As)
    z = z*wg
    p = p*wg
    return p, z
