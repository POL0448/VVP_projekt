"""
Modul pro definice funkce liftline, vykresleni grafu a vypis vypocitanych hodnot.
"""

import numpy as np
import matplotlib.pyplot as plt
import json 

class liftline_trida:    
    """
    Trida reprezentujici metodu liftline.

    Obsahuje funkce liftline, vykresleni grafu rozlozeni vztlaku a vraceni vystupnich hodnot cy_kr a Y.
    """

    def nacteni_dat(data_kridla):
        with open(data_kridla) as data_kridla:
            konst = json.load(data_kridla)
        return konst

    #nacteni vstupnich dat
    def nacteni_csv(jmeno1, jmeno2):
        data_koren = np.loadtxt(jmeno1, delimiter=',')
        data_konec = np.loadtxt(jmeno2, delimiter=',')
        return data_koren, data_konec

    def vypocet_cy(data_koren, data_konec, konst):
        alfainf_koren = konst['alfainf_koren']
        alfainf_konec = konst['alfainf_konec']
        indexy_koren = np.nonzero(np.isin(data_koren[:,0], alfainf_koren))
        cy_koren = data_koren[indexy_koren, 1]
        indexy_konec = np.nonzero(np.isin(data_konec[:,0], alfainf_konec))
        cy_konec = data_koren[indexy_konec, 1]
        return cy_koren, cy_konec

    def liftline(konst, cy_koren, cy_konec):
        m = konst['m']
        b_koren = konst['b_koren']
        b_konec = konst['b_konec']
        alfainf_koren = ['alfainf_koren']
        alfainf_konec = ['alfainf_konec']
        alfa0_koren = ['alfa0_koren']
        alfa0_konec = ['alfa0_konec']
        l = ['l']
        cy_koren_1 = cy_koren
        del cy_koren[-1] #vymaze posledni polozku vektoru
        del cy_koren_1[0] #vymaze prvni polozku seznamu
        cy_konec_1 = cy_konec
        del cy_konec[-1]
        del cy_konec_1[0]
        #vektor uhlu teta (rovnomerne rozdeleni 90/m-90 stupnu)
        teta = np.linspace(np.pi/(2*m), np.pi/2, m)
        #vzdalenost od korene, kde budou jednotlive rezy kridla (b)
        a = l/2*np.cos(teta)
        #vektor hloubky kridla, pouze pro polorozpeti, zacatek vektoru na konci kridla
        b_elipsa = np.sqrt(1-4*a**2/l**2)
        b = (b_koren-b_konec)*b_elipsa+b_konec
        #max pokles nabezne hrany (rozdil uhlu nabehu u korene a na konci kridla)
        e_konec = (alfainf_koren-alfainf_konec)*np.pi/180
        #25% korenove tetivy
        b4 = b_koren/4
        #pokles nabezne hrany na polorozpeti
        delta_h_konec = np.sqrt(2*b4**2-2*b4**2*np.cos(e_konec))
        delta_h = delta_h_konec*2*a/l
        #vypocet uhlu zkrouceni v konkretnich rezech
        e = np.arccos((2*b4**2 - delta_h**2)/(2*b4**2))
        #vypocet aerodynamickych uhlu nabehu na polorozpeti
        alfainf = alfainf_koren*np.pi/180 - e
        #prubeh rozlozeni uhlu nuloveho vztlaku
        alfa0 = ((alfa0_koren-alfa0_konec)*2*a/l+alfa0_konec)*np.pi/180
        #vypocet aerodynamickeho uhlu nabehu
        alfaa = alfainf - alfa0
        #numericke derivace soucinitele vztlaku podle uhlu nabehu
        cyalfainf_koren = 10*(cy_koren_1-cy_koren)/(np.pi/180)
        cyalfainf_konec = 10*(cy_konec_1-cy_konec)/(np.pi/180)
        #rozlozeni derivaci podel rozpeti (polorozpeti),zacatek vektoru je na konci kridla
        cyalfainf = (cyalfainf_koren-cyalfainf_konec)*2*a/l+cyalfainf_konec
        #vypocet vektoru konstant mi
        mi = b*cyalfainf/(4*l)
        #vektor lichych cisel
        n = np.linspace(1, 2*m, m)
        #matice ze soustavy rovnic
        matice = np.zeros((m, m))
        for i in range(m):
            for j in range(m):
                matice[i,j] = (np.sin(teta[i])+n[j]*mi[i])*np.sin(n[j]*teta[i])
        #vektor prave strany (transpozice na sloupcovy)
        vektor = (mi*alfaa*np.sin(teta)).reshape(-1, 1)
        #vypocet koeficientu Fourierovy rady
        A = np.linalg.solve(matice, vektor)
        print(A)

    def vykresleni(konst, A):
        m = konst['m']
        vinf = konst['vinf']
        nazev = konst['nazev']
        l = konst['l']
        alfainf_koren = konst['alfainf_koren']
        teta = np.linspace(np.pi/(2*m), np.pi/2, m)
        n = np.linspace(1, 2*m, m)
        #vypocet cirkulace pro profil
        soucet = np.zeros((m, 1))
        gamma = np.zeros((m, 1))
        for k in range(m):
            soucet[k] = np.sin(n*teta[k])@A
            gamma[k] = 2*l*vinf*soucet[k]
        #vypocet vztlaku pro profil (bezrozmerny rez)
        Y_gamma = 1.225*vinf*gamma
        #zrcadleni vektoru Y_gamma
        druhapulka = Y_gamma[::-1]
        #vykresleni rozlozeni vztlaku na kridle
        plt.figure()
        plt.plot(np.concatenate((Y_gamma, druhapulka)))
        plt.title(f"Rozlozeni vztlaku podel rozpeti kridla {nazev} pro alfainf_koren = {alfainf_koren}°")
        plt.xlabel('Rozpeti vyjadreno pomoci deleni m')
        plt.ylabel('Hodnota vztlaku')

    def hodnoty(konst, A):
        vinf = konst['vinf']
        l = konst['l']
        SAT = konst['SAT']
        #vypocet hodnoty vztlaku
        Y = 0.5*1.225*vinf**2*l**2*A[0]*np.pi
        #vypocet sounicetele vztlaku pro kridlo
        cy_kr = np.pi*l/SAT*A[0]
        #vrati hodnoty cy_kr a Y
        return cy_kr, Y