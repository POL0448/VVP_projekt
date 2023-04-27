"""
Modul pro konkretni vypocty.
"""

import numpy as np
import matplotlib.pyplot as plt
from .liftline import liftline_trida

class vypocty_trida(liftline_trida):  
    """
    Trida reprezentujici pouziti funkce liftline na konkretnich kridlech.
    """    
    def grafy(konst, A):
        nazev = konst['nazev']
        alfainf_koren = konst['alfainf_koren']
        x = len(alfainf_koren) - 1
        #vykresleni a ulozeni grafu do pdf souboru
        for i in range(x):
            liftline_trida.vykresleni(konst, A)
            plt.savefig(f'Rozlozeni vztlaku, {nazev}, alfainf_koren = {alfainf_koren}°.pdf')
        

    def vypsani_cy_kr_Y(konst, A):
        l = konst['l']
        vinf = konst['vinf']
        SAT = konst['SAT']
        vystup = liftline_trida.hodnoty(A, vinf, SAT, l)
        print(vystup)

    #ulozeni do csv
    def ulozeni_dat(konst, vystup):
        nazev = konst['nazev']
        np.savetxt(f'vvp_{nazev}.csv', vystup, delimiter=',')