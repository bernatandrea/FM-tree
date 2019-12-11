# FM-tree
Uni project from collegium Bioinformatic

#pokretanje
prvo se pozicionirajte u folder libdivsufsort-lite i onda od tamo u cmd-u napisite:make, izgenerirat ce vam .o datoteke i suftest.exe
nakon toga za pokretanje upisite: sufsort ime_datoteke(npr test-short.txt)

#upute
dokumetacija se nalazi u folderu dokumetacija. tamo je detaljno opisan obradjen primjer 
algoritma fm-tree zajedno sa primjerom.

u folderu: libdivsufsort-lite nalaze se datoteke algoritma za stvaranje sufixsalnog polja na sto optimalniji nacin.
vazne datoteke:
 divsufsort.c
	-tu su implementirani algoritmi za izradu bwt i sa
	-takodjer tu sam implemetirala par metodi koje ce nam biti potrebne za nas fm-tree algoritam(one se nalaze na kraju datoteke)

 divsufsort.h
	-header zajedno s deklaracijom metoda iz divsufsort.c kako bismo te metode mogli ugraditi sa include u nas program

 suftest.c
	-file za testiranje dosta stvari main-u sam zakomenitrala da nemam nepotreban ispis, vi slobodno odkomenitajte da vidite kako
	rade metode i sto ispisuju

 test.txt 
	-primjer od milijun acgt nizova

 test-short.txt
	-primjer za kraci niz iz dokumentacije (za provjere)
	
#todo
	-implementirati od gotovih stvari count i locate metodu

#napomena!!!
	prije svakog rada kod sebe na racunalu obavezno: git pull!! pa tek onda radite ili si napravite posebnu granu: git checkout -b ime-grane
	pa tamo commitajte
	
	
