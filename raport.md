# Sprawozdanie: Rekurencyjne mnozenie macierzy

## Wstep

Mnozenie macierzy jest jednym z podstawowych zadan obliczen numerycznych.
Dla dwoch macierzy kwadratowych `A` i `B` o wymiarze `n x n` klasyczny
algorytm mnozenia wykonuje rzad `n^3` mnozen i dodawan, a jego zlozonosc
czasowa wynosi `O(n^3)`. W wersji blokowej macierze dzieli sie na cztery
podmacierze i rekurencyjnie laczy wyniki osmiu mnozen blokowych.

Algorytm Strassena zmniejsza liczbe mnozen blokowych z osmiu do siedmiu.
Odbywa sie to kosztem wiekszej liczby dodawan i odejmowan, ale asymptotycznie
daje lepsza zlozonosc `O(n^(log2 7)) ~= O(n^2.81)`. Dla odpowiednio duzych
macierzy moze to prowadzic do przyspieszenia obliczen.

Celem projektu jest implementacja hybrydowego algorytmu rekurencyjnego:

- dla macierzy o rozmiarze `n <= 2^l` stosowany jest algorytm Strassena,
- dla wiekszych macierzy stosowany jest klasyczny algorytm blokowy (Binet),
- dodatkowo liczona jest liczba operacji zmiennoprzecinkowych oraz czas
  wykonania dla roznych wartosci parametru `l`.

## Opis algorytmu

Niech `multiply(A, B, l)` oznacza funkcje obliczajaca iloczyn dwoch macierzy
kwadratowych o rozmiarze bedacym potega dwojki.

1. Jesli `n = 1`, wynik jest rowny `[[A[0][0] * B[0][0]]]`.
2. Jesli `n <= 2^l`, wykonywany jest krok Strassena:
   - macierze sa dzielone na bloki `A11`, `A12`, `A21`, `A22`
     oraz `B11`, `B12`, `B21`, `B22`,
   - liczonych jest siedem iloczynow pomocniczych:

```text
P1 = (A11 + A22)(B11 + B22)
P2 = (A21 + A22)B11
P3 = A11(B12 - B22)
P4 = A22(B21 - B11)
P5 = (A11 + A12)B22
P6 = (A21 - A11)(B11 + B12)
P7 = (A12 - A22)(B21 + B22)
```

   - z iloczynow skladane sa bloki wyniku:

```text
C11 = P1 + P4 - P5 + P7
C12 = P3 + P5
C21 = P2 + P4
C22 = P1 - P2 + P3 + P6
```

3. Jesli `n > 2^l`, wykonywany jest krok klasyczny:

```text
M1 = A11B11   M2 = A12B21
M3 = A11B12   M4 = A12B22
M5 = A21B11   M6 = A22B21
M7 = A21B12   M8 = A22B22

C11 = M1 + M2
C12 = M3 + M4
C21 = M5 + M6
C22 = M7 + M8
```

4. Kazde dodawanie, odejmowanie i mnozenie skalarne zwieksza licznik operacji
   o `1`.

## Implementacja

Implementacja znajduje sie w pliku [binet_strassen.py](/Users/manoftheprincess/Helloworld/Uni/agh/sem1/rachunek macierzowy/program1/binet_strassen.py).
Najwazniejsze elementy programu to:

- klasa `MatrixMultiplier`, ktora przechowuje prog `threshold` oraz licznik
  operacji,
- metody `_add` i `_sub`, realizujace dodawanie i odejmowanie macierzy
  wraz z liczeniem operacji,
- metody `_binet` i `_strassen`, implementujace dwa warianty mnozenia
  rekurencyjnego,
- funkcja `run_experiment`, ktora generuje losowe macierze, uruchamia
  pomiary i zwraca czasy wykonania oraz liczbe operacji,
- funkcja `self_test`, ktora porownuje wyniki algorytmu z referencyjnym
  mnozeniem klasycznym.

Program zaklada, ze macierze sa kwadratowe i maja rozmiar bedacy potega
dwojki. Taka walidacja zostala dodana na wejsciu do funkcji `multiply`,
zeby uniknac blednych danych i ulatwic testowanie.

## Metodologia eksperymentu

Dla wybranych wartosci:

- `k` okreslajacych rozmiar macierzy `n = 2^k`,
- `l` okreslajacych prog `threshold = 2^l`,

generowane sa dwie losowe macierze `A` i `B`, po czym mierzony jest:

- czas wykonania algorytmu,
- laczny licznik operacji zmiennoprzecinkowych.

W notatnikach [code_notebook.ipynb](/Users/manoftheprincess/Helloworld/Uni/agh/sem1/rachunek macierzowy/program1/code_notebook.ipynb)
oraz [report_and_code.ipynb](/Users/manoftheprincess/Helloworld/Uni/agh/sem1/rachunek macierzowy/program1/report_and_code.ipynb)
znajduje sie kod do uruchomienia eksperymentow i narysowania wykresow.

## Oczekiwane obserwacje

Z teorii wynikaja nastepujace wnioski:

- dla malych macierzy przewaga Strassena moze byc niewielka lub zadna,
  poniewaz dodatkowe dodawania i narzut rekurencji sa kosztowne,
- wraz ze wzrostem rozmiaru macierzy liczba mnozen w Strassenie zaczyna
  dawac korzysc asymptotyczna,
- najlepszy prog przelaczania `l` jest zwykle kompromisem miedzy liczba
  operacji a rzeczywistym czasem wykonania.

W praktyce czysty model operacyjny i rzeczywisty czas moga dawac nieco inne
wnioski, poniewaz czas zalezy tez od narzutu interpretera Pythona,
alokacji pamieci i kopiowania podmacierzy.

## Wyniki praktyczne

Na podstawie uruchomionego eksperymentu otrzymano nastepujace wyniki
dla `k = 2, 3, 4, 5, 6` oraz `l = 1, 2, 3`.

### Czas wykonania

| `n = 2^k` | `l = 1` | `l = 2` | `l = 3` |
| --- | ---: | ---: | ---: |
| `4`  | `0.000079 s` | `0.000057 s` | `0.000063 s` |
| `8`  | `0.000527 s` | `0.000401 s` | `0.000373 s` |
| `16` | `0.005330 s` | `0.003192 s` | `0.003085 s` |
| `32` | `0.048467 s` | `0.025833 s` | `0.023831 s` |
| `64` | `0.211794 s` | `0.209029 s` | `0.191168 s` |

### Liczba operacji

| `n = 2^k` | `l = 1` | `l = 2` | `l = 3` |
| --- | ---: | ---: | ---: |
| `4`  | `216` | `247` | `247` |
| `8`  | `1792` | `2040` | `2017` |
| `16` | `14592` | `16576` | `16392` |
| `32` | `117760` | `133632` | `132160` |
| `64` | `946176` | `1073152` | `1061376` |

### Analiza wynikow

Wyniki praktyczne pokazuja, ze dla wszystkich badanych rozmiarow macierzy
najszybszy okazal sie wariant z `l = 3`, czyli z progiem `threshold = 8`.
Roznice sa male dla najmniejszych danych, ale rosna wraz z rozmiarem
problemu. Dla `n = 64` czas spadl z `0.211794 s` dla `l = 1` do
`0.191168 s` dla `l = 3`, czyli o okolo `9.7%`.

Jednoczesnie sama liczba operacji nie wskazuje na ten sam wariant jako
najlepszy. Najmniej operacji odnotowano dla `l = 1`, a mimo to ten wariant
nie byl najszybszy. Oznacza to, ze w implementacji w Pythonie narzut
rekurencji, tworzenia list i kopiowania podmacierzy ma istotny wplyw na
rzeczywisty czas wykonania. Innymi slowy, mniejsza liczba operacji
arytmetycznych nie musi od razu oznaczac krotszego czasu dzialania programu.

W badanym zakresie rozmiarow najlepszy kompromis miedzy liczba operacji a
czasem wykonania daje `l = 3`. Wariant `l = 1` minimalizuje liczbe operacji,
ale przegrywa czasowo. Wariant `l = 2` daje wynik posredni, a `l = 3`
wykorzystuje wieksze podproblemy Strassena na tyle skutecznie, ze przewaza to
nad dodatkowym kosztem organizacji obliczen.

Mozna wiec sformulowac praktyczny wniosek, ze dla tej konkretnej implementacji
i dla badanego zakresu `n` najlepszym wyborem jest prog `2^3 = 8`. Jest to
dobry przyklad tego, ze optymalny parametr hybrydowego algorytmu powinien byc
dobierany eksperymentalnie, a nie tylko na podstawie analizy asymptotycznej.

## Wnioski

Projekt realizuje wszystkie wymagane elementy:

- implementacje rekurencyjnego mnozenia blokowego,
- implementacje algorytmu Strassena,
- hybrydowe przelaczanie z parametrem `l`,
- zliczanie operacji zmiennoprzecinkowych,
- pomiar czasu wykonania,
- analiza wynikow praktycznych dla kilku progow przelaczania,
- kod do eksperymentow i wizualizacji wynikow.

Uzyskany program nadaje sie do porownania zachowania obu algorytmow dla
roznych rozmiarow danych. Dodatkowa walidacja danych oraz test zgodnosci z
algorytmem referencyjnym pomagaja upewnic sie, ze implementacja jest poprawna.
