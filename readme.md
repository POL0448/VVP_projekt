# Použití Glauertova řešení Prandtlovy rovnice pro výpočet rozložení vztlaku na křídle

## Textový popis
Tématem tohoto projektu je Glauertovo řešení Prandtlovy rovnice. Pomocí tohoto algoritmu můžeme vykreslit rozložení vztlaku na křídle letadla, vypočítat celkový vztlak generovaný křídlem a také součinitel vztlaku pro celé křídlo. Takto napsaný kód je použitelný pro obdélníkový a lichoběžníkový půdorysný tvar křídla a také pro zkroucená křídla.

Vstupem jsou geometrické charakteristiky křídla, počet řezů křídla ve směru podélné osy (směr "čumák - ocas"), vztlakové charakteristiky krajních řezů (které lze nahrát z 𝑐𝑠𝑣 souboru), rychlost nabíhajícího proudu a úhel, pod kterým je křídlo obtékáno.
Z těchto vstupů lze po úpravách dostat soustavu 𝑚 rovnic. Jejím řešením jsou liché koeficienty Fourierovy řady. Pomocí těchto koeficientů dostaneme právě výše uvedené výstupy.

## Funkcionality

Načtení vztlakových charakteristik z 𝑐𝑠𝑣 souboru.
- Implementace funkcí pro výpočet
- délek jednotlivých řezů křídla,
- zkroucení křídla,
- absolutních úhlů náběhu z úhlu náběhu a úhlu, při kterém křídlo negeneruje žádný
vztlak a
- vztlakových charakteristik pro jednotlivé řezy ze zadaných krajních řezů.
- Sestavení a vyřešení soustavy rovnic, odkud dostaneme koeficienty Fourierovy řady.
- Užití Fourierovy řady pro výpočet požadovaných výstupů.
- Uložení grafu s rozložením vztlaku do samostatného souboru.