"""
Balicek funkce obsahuje moduly liftline a vypocty, stejne jsou pouzity pro pouzivane tridy.
"""

from .vypocty import vypocty_trida
from .liftline import liftline_trida
__all__ = ["vypocty_trida", "liftline_trida"]