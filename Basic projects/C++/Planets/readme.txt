program do planet jest napisany obiektowo
każdy obiekt klasy reprezentuje niezafiksowany obiekt w przestrzeni kosmicznej
zawiera on w polach wskaźnik na tablice wskaźników wszytkich ciał niezafiksowanych(wspólne pole w c++ nazywa się to static)
zawiera także wektor z pozycjami i prędkościami w osi x i osi y (przechowuje cały dotychczasowy tor ciała , potrzebne do formuły 3 punktowej oraz energi i momentu pędu)
każdy obiekt ma też swoją ścieżke do pliku w którym wypisuje dane
konstruktor tworzy nowe ciało o masie i początkowych wartościach wektora  2D połóżenia i prędkości
Calculate_q1 liczy pierwsze przejście układu (założenie liniowej siły)
Calculate_Step liczy nastepne przejście układu jeśli ma przynajmniej policzony jeden pkt i pkt  z warunkór początkowych
Write data wypisuje wszystkie położenia i prędkości (do pliku o nazwie "obiekt n.data") (n argument metody)
Calc_and_Write_E_ANGM bierze wektor z wartościami czasowymi symulacji i wypisuje go łącznie z energią i momentem pędu układu w każdej chwili

aby program w mian zadziałał poprawnie nalzęy najpierw ustalić ilość symulowanych punktów i przedział czasowy a także wielkość tablicy ciał niebieskich(może być większa niż planowana , ale nie zalecam i trzeba się upewnić że nie używane elementy tablicy są NULL)
następnie należy stworzyć potrzebne niezafiksowane ciała niebieskie w tablicy (funkcja new)
następnie należy przypisać wskaźnik (tablice)ciał niebieskich do wspólnego pola obiektów lista_ob (nazwatablicy[istniejący indeks]->lista_ob=nazwatablicy)
następnie należy użyć gotowych metod do symulacji

q1 po wszytkich obiektach

pętla:
zwykły step po wszystkich obiektach

współczyniiki ciała zafikoswanego (pozycja masa) oraz wartośc stałej grawitacyjnej znajdują się w pliku Stellar_Body.cpp


!poniższy komentarz odnosi się raczej do wszytkich kompliatorów, jeśli nie chcemy wrzucać plików do projektu to w głównym pliku trzeba odkomentować:	!
//#include "Stellar_Body.cpp" //zawrzeć tą linijke jeśli plik nie jest w projekcie z resztą							!

jeśli użytkownik korzysta z code::blocks , prawdopodobne jest że program się nie skompiluje dopóki nie zepnie sie wszytkich plików w projekt code::blocks