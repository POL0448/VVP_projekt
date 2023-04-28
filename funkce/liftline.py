"""
Modul pro definice funkce liftline, vykresleni grafu a vypis vypocitanych hodnot.
"""

import numpy as np
import matplotlib.pyplot as plt
import json 

class Parametry:    
    """
    Trida reprezentujici metodu liftline.

    Obsahuje funkce liftline, vykresleni grafu rozlozeni vztlaku a vraceni vystupnich hodnot cy_kr a Y.
    """

    def __init__(self, data_kridla: str):
        with open(data_kridla) as self.data_kridla:
            konst = json.load(self.data_kridla)
        
        self.alfainf_koren = konst['alfainf_koren']
        self.alfainf_konec = konst['alfainf_konec']
        self.nazev = konst['nazev']
        self.alfa0_koren = konst['alfa0_koren']
        self.alfa0_konec = konst['alfa0_konec']
        self.b_koren = konst['b_koren']
        self.b_konec = konst['b_koren']
        self.m = konst['m']
        self.l = konst['l']
        self.vinf = konst['vinf']
        self.SAT = konst['SAT']

        self.j = len(self.alfainf_koren) - 1

        self.cy_kr = []
        self.Y = []

class Liftline(Parametry):
    #nacteni vstupnich dat
    def __init__(self, parametry):
        self.parametry = parametry

    def nacteni_dat(self, jmeno1, jmeno2):
        self.data_koren = np.loadtxt(jmeno1, delimiter=',')
        self.data_konec = np.loadtxt(jmeno2, delimiter=',')

    def vypocet_cy(self):
        indexy_koren = np.nonzero(np.isin(self.data_koren[:,0], self.parametry.alfainf_koren))
        self.cy_koren = self.data_koren[indexy_koren, 1]
        indexy_konec = np.nonzero(np.isin(self.data_konec[:,0], self.parametry.alfainf_konec))
        self.cy_konec = self.data_konec[indexy_konec, 1]

    def funkce_liftline(self, m, b_koren, b_konec, cy_koren, cy_koren_1, cy_konec, cy_konec_1,\
                         alfainf_koren, alfainf_konec, alfa0_koren, alfa0_konec, l):
        #vektor uhlu teta (rovnomerne rozdeleni 90/m-90 stupnu)
        self.teta = np.linspace(np.pi/(2*m), np.pi/2, m)
        #vzdalenost od korene, kde budou jednotlive rezy kridla (b)
        a = l/2*np.cos(self.teta)
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
        self.n = np.linspace(1, 2*m, m)
        #matice ze soustavy rovnic
        matice = np.zeros((m, m))
        for i in range(m):
            for j in range(m):
                matice[i,j] = (np.sin(self.teta[i])+self.n[j]*mi[i])*np.sin(self.n[j]*self.teta[i])
        #vektor prave strany (transpozice na sloupcovy)
        vektor = (mi*alfaa*np.sin(self.teta)).reshape(-1, 1)
        #vypocet koeficientu Fourierovy rady
        self.A = np.linalg.solve(matice, vektor)

    def vypocet_liftline(self):
        self.parametry.cy_kr = []
        self.parametry.Y = []
        self.vypocet_cy()
        for i in range(self.parametry.j):
            self.funkce_liftline(self.parametry.m, self.parametry.b_koren, self.parametry.b_konec,\
                self.cy_koren[0][i], self.cy_koren[0][i+1], self.cy_konec[0][i], self.cy_konec[0][i+1],\
                self.parametry.alfainf_koren[i], self.parametry.alfainf_konec[i], \
                self.parametry.alfa0_koren, self.parametry.alfa0_konec, self.parametry.l)
            self.cy_kr_uloz = np.pi*self.parametry.l/self.parametry.SAT*self.A[0]
            self.Y_uloz = 0.5*1.225*self.parametry.vinf**2*self.parametry.l**2*self.A[0]*np.pi
            self.parametry.cy_kr.append(self.cy_kr_uloz)
            self.parametry.Y.append(self.Y_uloz)

    def Y_gamma_fce(self):
        #vypocet cirkulace pro profil
        soucet = np.zeros((self.parametry.m, 1))
        gamma = np.zeros((self.parametry.m, 1))
        for k in range(self.parametry.m):
            soucet[k] = np.sin(self.n*self.teta[k])@self.A
            gamma[k] = 2*self.parametry.l*self.parametry.vinf*soucet[k]
        #vypocet vztlaku pro profil (bezrozmerny rez)
        self.Y_gamma = 1.225*self.parametry.vinf*gamma
        #zrcadleni vektoru Y_gamma
        self.druhapulka = self.Y_gamma[::-1]

    def vykresleni(self):
        #vykresleni rozlozeni vztlaku na kridle
        self.vypocet_cy()
        for i in range(self.parametry.j):
            self.vypocet_liftline()
            self.Y_gamma_fce()
            plt.figure()
            plt.plot(np.concatenate((self.Y_gamma, self.druhapulka)))
            plt.title(f"Rozlozeni vztlaku podel rozpeti kridla {self.parametry.nazev} pro alfainf_koren = {self.parametry.alfainf_koren[i]}°")
            plt.xlabel('Rozpeti vyjadreno pomoci deleni m')
            plt.ylabel('Hodnota vztlaku')

    def vypsani_hodnot(self):
        #vypise hodnoty cy_kr a Y
        self.vypocet_liftline()
        print(f"{self.parametry.nazev}")
        for i in range(self.parametry.j):
            print(f"Součinitel vztlaku: {self.parametry.cy_kr[i]}")
            print(f"Vztlak: {self.parametry.Y[i]}")
        
    def ulozeni_grafu(self, path):
        #ulozeni grafu do pdf souboru
        for i in range(self.parametry.j):
            plt.savefig(path + f'Rozlozeni vztlaku, {self.parametry.nazev}, alfainf_koren = {self.parametry.alfainf_koren[i]}°.pdf')

    #ulozeni do csv
    def ulozeni_dat(self, cesta):
        self.vypocet_liftline()
        np.savetxt(cesta + f'cy_kr_{self.parametry.nazev}.csv', self.parametry.cy_kr)
        np.savetxt(cesta + f'Y_{self.parametry.nazev}.csv', self.parametry.Y)