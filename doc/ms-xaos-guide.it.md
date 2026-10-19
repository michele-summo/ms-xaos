# MS XaoS

> Traduzione della [guida in inglese](ms-xaos-guide.md) della versione 2.10.
> Dove le due differiscono, fa fede l'originale inglese.

Un fork di [XaoS](https://github.com/xaos-project/XaoS) 4.3.3. Tutto ciò che fa
l'originale lo fa ancora; questa guida descrive ciò che è stato aggiunto o
cambiato, e perché. Versione 2.4.

## Due eseguibili

`XaoS.exe` calcola con 64 bit di mantissa, `XaoS-quad.exe` con 113. Entrambi
sono compilati dagli stessi sorgenti — la precisione è un tipo scelto in
compilazione, non un'impostazione — e quello quad si chiama **MS XaoS Quad**,
così i due si distinguono mentre sono aperti.

Quad zooma circa quindici ordini di grandezza più in profondità. Renderizzando
la stessa posizione a 1e-18, l'eseguibile normale conserva il 21% del
dettaglio, quello quad il 92%. Non è il predefinito perché `__float128` è
emulato in software: circa otto volte più lento sulle formule incluse e
quindici su quelle definite dall'utente.

Una posizione salvata dall'eseguibile quad contiene `(precision 113)`. Aprirla
con una build meno precisa lo segnala e la disegna comunque; aprirla con una
versione che non conosce il comando viene rifiutato, ed è voluto — un'immagine
calcolata a 113 bit non si può riprodurre a 64, e non aprirla è più onesto che
disegnare qualcos'altro.

## I numeri nei dialoghi

Un campo che chiede un numero ne accetta uno scritto in qualunque modo
leggibile, 1e-18 come 0.5. I dialoghi a campo singolo usavano una casella a
scorrimento, che tiene due decimali se non le si dice altro, così 0.01
sopravviveva e 0.001 diventava zero — una convergenza di Newton di un
milionesimo non si poteva proprio inserire.

Ciò che viene mostrato è il testo più breve che si rilegge come il numero
esatto contenuto, così un milionesimo di milionesimo appare come 1e-18 e non
come 9.99999999999999999978e-19, che è ciò che stampare ogni cifra della build
fa di un decimale che un float binario non può rappresentare esattamente. Una
coordinata che ha bisogno di tutte e ventuno le cifre le riceve comunque:
averne bisogno equivale a non rileggersi invariata senza di esse.

La forma esponenziale si usa da mille miliardi in su e sotto i cento
miliardesimi, e la decimale semplice in mezzo — così venti è venti, cento
miliardi si scrivono per esteso, e 1e-8 no. Lasciata a sé, la stampa
sceglierebbe confrontando l'esponente con il numero di cifre richieste, il che
non serve quando le cifre sono scelte per brevità: basta una cifra per
rileggersi come venti, e venti usciva come 2e+01.

## Posizioni salvate

La vista viene scritta alla piena precisione della build — 21 cifre
significative normalmente, 36 per quad. Prima veniva scritta con un numero di
cifre ricavato dalla profondità dello zoom, che ne riportava troppo poche e si
fermava a 17, così una posizione profonda si riapriva a una frazione di pixel
da dove era stata salvata. Con molte iterazioni è un'immagine visibilmente
diversa.

Tutto ciò che una posizione può contenere ora sopravvive invariato a un
salvataggio e un caricamento: tutti i 169 esempi inclusi e 141 posizioni
generate, che coprono ogni formula, piano, modalità di colorazione e filtro,
fanno il giro identici pixel per pixel.

## Annulla

Modifica → **Annulla** (`U`, o `Ctrl` `Z`) torna indietro di un passo e
**Ripeti** (`Ctrl` `Y` su Windows) va di nuovo avanti. Un passo è qualunque
cosa cambi l'immagine: un comando — una formula, una tavolozza, una modalità di
colorazione, un dialogo, uno zoom fisso — una pressione del mouse con lo zoom,
lo spostamento o la rotazione che avvia, uno zoom a selezione, e le iterazioni
impostate con le frecce. Un comando che non cambia nulla, come Ricalcola, non è
un passo, e annulla non torna più indietro dell'immagine con cui XaoS si è
aperto, dopo la configurazione e la riga di comando. Ogni passo è conservato
come posizione salvata, così annulla rimette tutto ciò che una posizione
contiene, esattamente come riaprirla; si conservano gli ultimi 256.

Sbagliava in due modi. Il mouse conservava il suo passo al rilascio, mentre lo
zoom stava ancora rallentando, così annulla riportava un'immagine a metà; ora
il passo si conserva dove la vista si ferma. E ogni comando conservava un
passo, che cambiasse qualcosa o no, mentre annulla tornava indietro di due
presupponendo di aver conservato l'ultimo da sé — dopo un Ricalcola riportava
l'immagine già sullo schermo. Ora un'immagine si conserva solo dove differisce
dall'ultima conservata, e annulla e ripeti si muovono di un passo lungo di
esse. Un passo nuovo dopo un annulla taglia ciò che si sarebbe potuto ripetere,
come in un editor.

Ciò che si è mosso da sé dall'ultimo passo — il pilota automatico, una
rotazione continua — si conserva come passo a sé quando si preme Annulla, così
che Ripeti vi ritorni.

Annulla ha trovato un errore nelle posizioni. Caricarne una con una tavolozza
inclusa dopo una tavolozza personalizzata lasciava la tavolozza personalizzata
segnata come in uso, così la posizione veniva salvata di nuovo come tavolozza
personalizzata: non ciò che era stato caricato, e per annulla non il passo a
cui era tornato, il che costava il ripeti. Ora una tavolozza inclusa toglie il
segno, come aveva sempre fatto una generata.

## Formule utente

**Riferimento.** Aiuto → Riferimento delle formule utente elenca ogni funzione,
variabile e notazione che il linguaggio delle formule accetta, in una finestra
che si può lasciare aperta accanto alla formula che si sta scrivendo. È
raggruppato secondo lo scopo delle funzioni, con le funzioni speciali — `erf`,
`gamma`, `lambertw` — in una sezione a sé, perché nessuna è elementare e non le
si usa per caso. Una quarta scheda elenca i numeri che compaiono come argomenti
e hanno un significato particolare — rispetto a quale forma misura una
trappola, con quale specchio ripiega il caleidoscopio — invece di seppellirli
nella descrizione della funzione che li prende. Una quinta mostra un'immagine
di ogni tassellatura che `randsctile` accetta, perché quarantacinque numeri con
una riga di parole ciascuno dicono poco delle forme.

Una sesta, **Utilità**, scrive formule invece di descriverle. Per ora contiene
un `poly` casuale, da quattro campi:

| | |
| --- | --- |
| lunghezza minima e massima | m, il numero dei coefficienti, estratto tra le due; il grado è m − 1. Servono entrambe, come numeri interi da 1 a 254 — fin dove arriva `poly`, dato che il parser accetta 255 argomenti per chiamata e `z` è uno di essi |
| modulo minimo e massimo | la grandezza di ogni coefficiente, estratta uniformemente tra i due. 0 e 1 se lasciati vuoti; al massimo 1000000 |

Ogni coefficiente è complesso, `{re,im}`, con il suo argomento estratto a caso
sull'intero giro, e ogni parte è scritta con quattro cifre significative del
suo modulo, così che un intervallo da 0.001 a 0.002 non venga arrotondato via.
Nella formula non c'è `c`. Una formula utente fa partire `z` dal pixel sia in
modalità mandelbrot sia in modalità julia, quindi ciò che disegna è l'insieme
di Julia pieno del polinomio stesso, uguale in entrambe; `c` va dove serve,
come in `poly(...) + c`.

Ogni modifica a un campo genera una formula nuova, come fa **Generane un
altro**, e **Copia** la mette negli appunti per Frattale → Formula utente.
Finché una lunghezza è vuota, o un campo non è un numero, è fuori dal suo
intervallo o ha il minimo sopra il massimo, non c'è formula, solo una riga che
dice perché.

Ogni riga dell'elenco delle funzioni dice cosa prende la chiamata e in quale
ordine, con nomi che dicono a cosa serve ogni posto — "base, exponent", "z,
low, high", "before, after, n" — invece che con lettere, che non dicevano quale
posto fosse quale; gli stessi nomi sono quelli che mostra la riga sotto la
barra della formula. Una convenzione li attraversa tutti: `z` è il valore su
cui la funzione lavora, qualunque cosa vi sia scritta davvero. L'elenco e quei
conteggi sono entrambi controllati da un test rispetto alla tabella del parser
stesso, così una funzione che guadagna o perde un argomento non può continuare
a essere descritta con quelli che aveva. Un argomento che si può omettere è
mostrato tra parentesi con il valore che prende quando lo si omette:
`z, [power=1], [angle_factor=1]`.

**Gli argomenti sono separati da una virgola, e da nient'altro.** Un tempo
serviva anche il punto e virgola: il parser lo trasformava in virgola prima di
leggere qualunque cosa, così `LOGN(5;Z^2)+C` e `LOGN(5,Z^2)+C` erano la stessa
formula, e lo stesso valeva per le due parti di un numero complesso. Non è più
così. Un punto e virgola ora è rifiutato dove compare — *Operatore non valido:
;* — e la virgola è il separatore ovunque, come lo è in tutto il resto della
sintassi.

### La sonda della tavolozza

**Vista → Sonda della tavolozza (clic su un punto)** trasforma il puntatore in
una sonda: fai clic in qualunque punto del frattale e la riga dei messaggi dice
**da quale punto della tavolozza è stato disegnato quel punto**, come numero e
non come colore — *Posizione nella tavolozza 12.375*, e nient'altro.

L'editor della tavolozza dispone **31 colori**, e la tavolozza si costruisce
passando da ciascuno al successivo **otto celle alla volta**. Quindi il numero
è la cella divisa per otto: `12.000` è esattamente il dodicesimo colore
dell'editor, `12.375` è a tre ottavi del cammino dal dodicesimo al tredicesimo.
Una cella è un ottavo di colore, più fine di quanto chiunque legga dall'editor.

* Con la **tavolozza personalizzata** in uso il numero viene riportato tra
  **`0` e `31`**, che è la numerazione dell'editor stesso.
* Con qualunque altra tavolozza viene **lasciato dove cade e va ben oltre 31**
  — una tavolozza di 65534 celle arriva a 8191. Il numero fuori intervallo è
  ciò che dice che la tavolozza non è quella che disegna l'editor.
* Un punto che ha preso il **colore interno** legge *Posizione nella tavolozza:
  nessuna*: quel colore è un unico tono piatto che nessuna cella della
  tavolozza rappresenta.

Serve meno a sapere "che colore è questo" che **quale parte della tavolozza
l'immagine usa davvero**, e la risposta spesso fa riflettere. Con la semplice
colorazione `iter` a 60 iterazioni le celle vanno da 1 a 59 — da `0.125` a
`7.375` — quindi l'immagine non tocca mai i colori dall'8 al 30 dell'editor e
modificarli non cambia nulla. Una colorazione che si allarga, per esempio
`iter + real`, percorre invece tutto l'anello, più volte.

Legge il punto e non lo schermo, quindi è esatta: lo stesso colore può stare in
due celle, e un'immagine a colori reali mescola due celle vicine, quindi nessun
colore indica una cella sua. Come lo zoom a selezione, la sonda restituisce il
puntatore alla modalità Julia veloce finché questa è attiva.

**La barra.** Frattale → Formula utente apre un dialogo largo circa metà di uno
schermo 1280×720 — circa seicentocinquanta pixel, meno su uno schermo più
piccolo — la cui barra va a capo quando la formula non ci sta su una riga. Va a
capo **tra i nomi, mai dentro uno di essi**: un nome o un numero (`randsctile`,
`z_p3`, `0.6`) resta intero, e una riga può finire in qualunque altro punto —
dopo una parentesi, una virgola, un operatore. Il ritorno a capo di Qt non
serviva a questo: spezza dove lo consentono le regole della prosa, che in una
formula significa dopo qualche virgola e poco altro, e poi taglia un nome dove
deve. Solo un nome più largo dell'intera barra viene mai tagliato. La barra
cresce di una riga alla volta con la formula, e quando il dialogo non starebbe
più in uno schermo 720p smette di crescere e scorre. Per tutto il resto resta
una formula su una riga: Invio è OK e non un a capo, e in una formula incollata
gli a capo e le tabulazioni diventano spazi, che il parser legge come niente.

**Cronologia...** sotto la barra apre le formule usate in precedenza — le dieci
più recenti, dalla più nuova, ciascuna una volta — che la vecchia barra su una
riga teneva nel suo elenco a discesa. Una formula troppo lunga per la finestra
è elencata con il suo inizio, tre puntini e la sua lunghezza,
`poly(z, {0.01097,0.00699}, {-0.05... (157 caratteri)`, così che due che
iniziano allo stesso modo si possano ancora distinguere; l'inizio finisce tra i
nomi, come una riga della barra, e se ne vede di più quanto più la finestra è
larga. Fai doppio clic su una, oppure selezionala e premi Usa, e finisce intera
nella barra. Una formula viene aggiunta quando il dialogo è accettato con essa,
e solo allora: caricare una posizione che contiene una formula utente,
riprodurne una, o renderizzarne una con `-render` esegue lo stesso comando
`usrform` ma non la mette più nella cronologia, dove ciascuna di queste
prendeva il posto di una che avevi scritto tu.

**Freccia su e freccia giù** scorrono le stesse formule senza lasciare la
barra. Freccia su, dalla prima riga della barra, porta la formula usata prima
di quella mostrata; freccia giù, dall'ultima riga, torna verso la più recente e
poi a ciò che la barra conteneva prima che cominciassi. Le due piccole frecce
accanto alla barra fanno lo stesso, ciascuna grigia quando da quella parte non
c'è altro. Una formula che è già nella barra viene saltata, e ciò che scrivi
diventa il punto da cui si riparte.

L'ultima formula usata è mostrata in grigio chiaro nella barra: intera finché
la barra è vuota, e il resto dopo ciò che hai scritto finché ciò che scrivi ne
è l'inizio. Tab la mette intera nella barra. Appena ciò che scrivi prende
un'altra strada il testo in grigio sparisce, e Tab passa ai pulsanti come
prima. La barra si apre vuota, con quella formula in grigio, e OK su una barra
vuota usa quella; solo la primissima volta, quando non è ancora stata usata
nessuna formula, si apre con la formula in uso, `z^2+c`.

**Vista → Stato** mostra la formula allo stesso modo: su una riga se ci sta
nella finestra, su due se no — la seconda sotto la prima, andando a capo tra i
nomi — e quando nemmeno due bastano, la seconda finisce con tre puntini e la
sua lunghezza, `... (918 caratteri)`. Le righe sotto si spostano in giù per
fare spazio, e tutto viene ridisposto quando la finestra cambia larghezza.
Prima era una sola riga che usciva dal bordo dell'immagine.

**Sotto la barra c'è una riga che dice cosa prende la chiamata in cui ti
trovi**, con l'argomento che stai scrivendo in grassetto. Metti il cursore
dentro `randsc(` e si legge *randsc(**seed**, [size=1+i],
[degradation=0.5+0.5i], [kaleidoscope=1], [kaled_mirror=0], [skew=0])*; supera
una virgola e il grassetto si sposta con te. Non le serve altro della formula
che una parentesi aperta con un nome davanti, quindi c'è mentre la formula è
scritta a metà, che è quando serve. Metti il cursore fuori da ogni chiamata e
la riga è vuota.

Sei delle posizioni incluse erano scritte alla vecchia maniera e sono state
riscritte con le virgole: `circle`, `heart`, `helloween`, `pentafrac`,
`warriormask` e `burnship`, tutte sotto **File → Apri → examples/Malczak**.
**Una tua posizione che ha un punto e virgola nella formula non si caricherà**
finché la formula non sarà scritta con le virgole — il messaggio indica il
carattere, e la soluzione è aprire la formula e sostituirli.

**Lasciare un posto vuoto.** Un argomento mostrato tra parentesi si può
omettere in mezzo a una chiamata come alla fine, non scrivendo nulla tra i due
separatori: `julian(z, ,3)` dà il primo e il terzo e lascia dire alla funzione
qual è il secondo. Gli spazi non fanno differenza — `f(z, ,5)` e `f(z,,5)` sono
la stessa chiamata — e un argomento che non è tra parentesi va comunque
scritto, così `poly( ,1,2)` e `sin(,)` sono rifiutati invece di essere
indovinati.

Cosa significhi un posto vuoto lo decide la funzione. La maggior parte prende
il valore predefinito che dichiara; un coefficiente di `poly` lasciato vuoto è
un termine che non c'è, così `poly(z,1, ,1)` è `z^2 + 1`; e un ramo di `ifiter`
o `ifiterl` lasciato vuoto ripete quello prima, che è il modo per dare a un
ramo più di un passo nel ciclo. `ifiter(f(z), , ,g(z), , , , )` esegue `f` per
tre passi e `g` per cinque, e non costa più da valutare che scriverli per
esteso: dove porta ogni scelta si stabilisce una volta, quando la formula viene
analizzata.

Una formula utente non ancora scritta dice `z^2+c`. Prima diceva il burning
ship, che è un frattale a sé e un enigma da trovarsi come punto di partenza.

**`ifiterf(rest, last)`** — valuta `rest` in ogni passo tranne l'ultimo e
`last` in quello. L'ultimo passo è l'ultimo consentito dal limite delle
iterazioni: una formula non ha modo di sapere quale passo sarà quello che esce,
perché dipende dal valore che non ha ancora prodotto. Viene valutato solo
quello scelto.

**`ifiterr(before, after, n)`** — valuta `before` finché il numero del passo è
sotto `n` e `after` da `n` in poi. Viene valutato solo quello scelto, come per
gli altri due, anche se prima la soglia andava insegnata al parser: ora un
argomento può essere segnato come letto dal selettore anziché scelto da esso, e
allora viene valutato prima che la scelta sia fatta.

**Le variabili.** `z` è il valore corrente; `c` è il punto, che è il pixel in
modalità mandelbrot e la costante in modalità julia; `x` è la coordinata
semplice, che è il pixel in entrambe le modalità e non cambia da un passo
all'altro; `n` è il numero dell'iterazione. Nella *Inizializzazione utente* `z`
è dove `z` sarebbe partito se non ci fosse stata inizializzazione, che è lo
stesso valore di `x`, quindi un'inizializzazione fatta solo di `z` dice
esattamente ciò che dice non dire nulla.

Quest'ultima parte è nuova, e lo è anche poterla dire: `z` non era registrato
in nessuno dei due parser per cui passa l'inizializzazione, e `x` in quello
sbagliato, così un'inizializzazione che nominava uno dei due veniva rifiutata
dal thread che doveva calcolarla -- mentre il dialogo, il cui parser ha
entrambi, l'accettava e non dava alcun segno che poi venisse buttata via. `n`
restava con ciò che vi aveva lasciato il pixel precedente, ossia
un'inizializzazione che lo leggeva leggeva l'ordine in cui i pixel venivano
calcolati.

**`p1`, `p2` ... `p9999`** — il valore che `z` aveva in un passo precedente:
`p1` nel passo prima di questo, `p2` in quello prima ancora, `p9999`
novemilanovecentonovantanove passi indietro. `p` è un altro nome per `p1`.
Prima che ci siano stati così tanti passi valgono ciò da cui parte la
cronologia, che decide *Frattale → Imposta i valori p alla prima iterazione*:
il punto che si itera quando è attivo, zero quando non lo è.

Prima si fermavano a `p6` e venivano conservati spostando di un posto un array
di sei posti a ogni passo — una copia per posto per passo, pagata che la
formula ne nominasse qualcuno o no. Ecco perché si fermavano a sei: a novemila
sarebbero stati 320 KB di copie per passo, misurati in ventiquattro volte
l'intero costo di iterare `z^2+c`. Ora la cronologia è un anello, scritto una
volta per passo, e se ne leggono solo i posti che una formula nomina davvero,
così `p9999` costa quanto `p1` e una formula che non ne nomina nessuno non
costa nulla. La finestra del riferimento li elenca tra le variabili.

**Suffissi su una variabile.** `z`, `c`, `x`, `p` e da `p1` a `p9999` accettano
suffissi che stanno per le chiamate fatte più spesso su di essi:

| suffisso | sta per |
| --- | --- |
| `_b` | `bship(…)` |
| `_bi` | `bshipi(…)` |
| `_br` | `bshipr(…)` |
| `_pM` | `parchment(…, M)` |
| `_paM` | `parchmenta(…, M)` |

Si leggono da sinistra a destra, ciascuno avvolgendo ciò che hanno prodotto i
precedenti: `z_p3` è `parchment(z,3)`, `c_b_p2` è `parchment(bship(c),2)` e
`p12_p2_p3` è `parchment(parchment(p12,2),3)`. *M* fa parte del nome, quindi si
scrive in cifre — un numero intero, da uno in su; `z_p0` è rifiutato con un
messaggio che lo dice. Non su `n`, che è un conteggio e non ha una seconda
componente su cui uno di essi possa lavorare.

Vengono espansi nelle chiamate prima che la formula sia letta, così una
variabile con suffisso disegna esattamente ciò che disegnano le chiamate
scritte per esteso e costa quanto costano; la formula viene salvata come è
stata scritta. Tutto ciò che non è uno dei cinque viene lasciato com'è e
rifiutato come il nome sconosciuto che è.

**Rimossi.** `powi`, `powdc` e `logcn` erano secondi nomi di `pow` e `logn`;
`rad`, `deg` e `sign` erano elencati senza alcuna implementazione dietro, e un
secondo `trunc` era oscurato da quello funzionante. Una formula che usa un
alias rimosso deve usare il nome che resta — la funzione è identica.

**`erf(z)`** — la funzione degli errori sul piano complesso. Precisa a circa
tre ulp entro un bailout di due, dove un frattale itera davvero.

**`poly(z, k1, k2, ..., km)`** — un polinomio in `z`:

    k1*z^(m-1) + k2*z^(m-2) + ... + k(m-1)*z + km

Il primo coefficiente scritto moltiplica la potenza più alta e l'ultimo sta da
solo, così la chiamata si legge nell'ordine in cui si pronuncia il polinomio.
Calcolato con la regola di Horner, che fa m−1 moltiplicazioni invece delle
m(m−1)/2 che richiederebbe elevare ogni potenza separatamente, ed è per giunta
la più precisa delle due.

    poly(z,1,0,0)+c            il Mandelbrot, scritto per esteso
    poly(z,1,0,0,{0.7,0.2})    z^3 + 0.7+0.2i

**`randsc(seed, size, degradation, kaleidoscope, kaled_mirror, skew, skew_mode, selfsim)`**
— rumore coerente sul punto, che dà macchie invece di neve pixel per pixel.
`size` (predefinito `1+i`) è la larghezza media di una macchia lungo l'asse
reale e la sua altezza lungo quello immaginario. È **rumore a gradiente**,
un'ottava del moto che sommano le [modalità fbm](#moto-browniano-frazionario) e
`fbm()`. Era rumore di valore, e per la stessa ragione delle loro disegnava il
suo reticolo: le macchie diventavano piatte lungo ogni sua linea e uscivano
squadrate. Le macchie sono larghe come prima — il reticolo sta a tre quinti
della dimensione, e il campo smette di somigliare a sé stesso della metà a
undici sedicesimi della dimensione, prima e dopo — e i valori si distribuiscono
come prima, 0.2145 dal centro, con la curva di contrasto dell'fbm che riporta
la distribuzione più stretta del rumore a gradiente su quella del rumore di
valore. Ma ogni valore è diverso, e costa 167 ns a chiamata contro 147. I
cinque mosaici sono il loro reticolo e sono invariati fino all'ultimo bit.
`degradation` (predefinito `0.5+0.5i`, che le dimezza a ogni passo) le
rimpicciolisce man mano che l'iterazione procede: la dimensione viene
moltiplicata per esso a ogni passo, componente per componente, così `0.5+0.2i`
su `1+i` dà `1+i` al primo passo, poi `0.5+0.2i`, poi `0.25+0.04i`. Uno zero in
una qualunque componente di uno dei due argomenti restituisce zero invece di
dividere per zero. Poi vengono i due argomenti del caleidoscopio, `skew`
(predefinito `0`), `skew_mode` e `selfsim`, ciascuno descritto più avanti. Solo
il seme è obbligatorio.

### L'inclinazione, e perché questi campi sono difficili da colorare

Il motore colora con le due componenti dell'orbita: `real` e `zmag` ne leggono
una, `imag`, `angle` e `real / imag` l'altra. Questi campi restituiscono un
solo numero per un'intera cella, sull'asse reale, quindi ogni modalità disegna
**un tono per cella** — `zmag` e `iter + real` sembrano il valore troncato, e
le altre tre non hanno nulla da leggere. Misurato su novecento pixel: da 14 a
126 valori per le prime due, ed esattamente **uno** per le altre.

Variare il valore dentro una cella non si può fare con un solo numero reale. Lo
stesso numero decide se il punto esce, quindi un valore che si muove attraverso
una cella fa uscire metà della cella e lascia dentro l'altra metà, e la cella
esce **tagliata in due da una linea retta**. È stato provato due volte e ha
tagliato i mosaici entrambe le volte.

`skew` lo fa invece con la seconda componente, e **ruota** il valore anziché
scalarlo. Di quanto ruota segue **quanto il punto sta fuori dal centro della
sua cella**: con `out` che va da zero al centro a uno sul bordo, il valore
viene ruotato di due volte l'arcotangente di `skew_re*out + skew_im`.

`out` è misurato **nella geometria della cella stessa**, quindi ciò che il
colore disegna dentro una cella è la forma in cui il campo è diviso, non un
insieme di linee stese sopra — zero al centro di una cella, uno lungo tutto il
suo bordo:

| campo | cosa disegna il colore dentro una cella |
| --- | --- |
| `randsc` | le sue stesse macchie — la rotazione segue il livello, quindi i contorni sono quelli del campo |
| `randscq` | quadrati, attorno al centro del quadrato |
| `randsch` | esagoni, attorno al centro dell'esagono |
| `randsct` | triangoli, attorno al centro del triangolo |
| `randscp` | il suo poligono, rimpicciolito passo dopo passo — non ci sono due celle della stessa forma |

La **parte reale** dell'inclinazione è quel gradiente attraverso una cella. La
**parte immaginaria** è una rotazione uniforme condivisa da tutta la cella:
sposta il colore di una cella senza disegnarci dentro nulla, così `0.4+0.2i` dà
insieme i contorni della forma e uno spostamento per cella.

### Cosa disegna l'inclinazione: `skew_mode`

Il settimo argomento dice **cosa varia dentro la cella**. I numeri **si
sommano** anziché escludersi a vicenda, così una chiamata può chiederne due e
averli entrambi.

| | |
| --- | --- |
| **1** (e `0`) | **la forma** — il contorno stesso della cella, rimpicciolito passo dopo passo. Ciò che riceve una chiamata che non dà `skew_mode`, e ciò che si usa quando non si nomina né 1 né 2 |
| **2** | **la rosetta** — l'angolo attorno al centro della cella, ripiegato in tanti giri quanti spicchi ha il caleidoscopio, così ogni cella porta la simmetria dell'immagine stessa. Semplici raggi quando nulla è ripiegato |
| **4** | **attraverso lo spicchio** — una rotazione che segue dove sta il punto nel suo spicchio del caleidoscopio, uguale in ogni spicchio |
| **8** | **radiale** — si muove anche il modulo oltre all'angolo |

`3` è la forma e la rosetta insieme, che fa una spirale. `11` è queste due con
in più il modulo che si muove.

**`4` segue la geometria del caleidoscopio.** Ruota il valore di nulla sui
bordi di ogni spicchio e di tutta la rotazione che la parte immaginaria dà da
sola — due volte la sua arcotangente, 127 gradi per `2` — lungo la linea per il
suo centro, e in mezzo della rotazione della parte immaginaria per quanto si è
attraversato. Ogni spicchio riceve la stessa cosa, così un caleidoscopio a sei
è ancora un caleidoscopio a sei, specchi compresi, e la rotazione è continua
attraverso i bordi e i centri, che sono dove la piega rispecchia. Non fa nulla
quando nulla è ripiegato, o quando la parte immaginaria dell'inclinazione è
zero. Niente trigonometria, e non costa nulla a 64 bit di mantissa e un decimo
di chiamata a 113.

All'inizio era **una rotazione diversa in ogni spicchio** — la parte
immaginaria per il numero dello spicchio più uno — pensata per far smettere le
copie di essere identiche. È l'unica cosa che le copie di un caleidoscopio non
possono essere: un'immagine a simmetria sei usciva come sei fette scollegate
con un taglio dritto a ogni giunzione. Il multiplo andava anche a sbattere
contro il mezzo giro, così quattro fette su sei stavano quasi sull'asse reale,
e lì `real / imag` le disegnava come neve comunque `selfsim` facesse la media
del rumore.

Di ciò resta una proprietà di `real / imag`, che divide per la parte
immaginaria di *z*: dovunque una rotazione la porti attraverso lo zero la
colorazione salta da un capo all'altro della tavolozza, e gli ultimi frammenti
dei passi fini decidono il lato pixel per pixel. Una rotazione che spazza lo
spicchio lo attraversa lì per la maggior parte delle impostazioni. Su un
`randsct` a simmetria sei colorato con `real / imag`, con `selfsim` a uno, una
parte immaginaria di `0.25` non lasciava alcuna grana e `2` ne lasciava bande
accanto ai bordi; una colorazione che non divide per una componente non ne ha.

**`8` è quella che risponde a `zmag`.** Una rotazione lascia il modulo dov'è, e
`zmag` e il bailout leggono il modulo e nient'altro — quindi nessun'altra
modalità può raggiungerli. Il prezzo è quello che l'inclinazione era stata
progettata per evitare: l'uscita si muove con il modulo, quindi la figura viene
tagliata dovunque il bailout cada dentro una cella. Abbinala a un bailout
circolare e a una colorazione che legge le componenti separatamente, o accetta
il taglio di proposito.

Un'inclinazione nulla lascia ogni modalità senza effetto, quindi `skew_mode` da
solo non cambia nessuna immagine.

Prima era stata scritta una rampa dritta — `skew_re*du + skew_im*dv` attraverso
la cella — ed è ciò che questa sostituisce. Colorava, ma disegnava le stesse
bande diagonali, verticali od orizzontali in ogni cella allo stesso modo,
comunque fosse diviso il campo, e le bande rovinavano l'immagine.

A **`0`**, che è ciò che riceve una chiamata che non lo nomina, la rotazione è
nulla, la parte immaginaria resta a zero, e ogni valore è il numero che è
sempre stato **fino all'ultimo bit** — le firme di riferimento di prima che
l'argomento esistesse passano invariate, e una posizione salvata si renderizza
identica. La posizione `talc`, tre chiamate di rumore e una `c` su 345600
valori, dà la stessa firma della build prima dell'argomento.

**Quanto colore dà.** L'indice del motore è
`(iter + quantity) * speed + shift`, quindi una variazione di uno nella
grandezza vale una banda a una velocità di uno. La variazione di `imag` dentro
una cella:

| inclinazione | 0.02 | 0.05 | 0.1 | 0.2 | 0.4 | 0.6 | 1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `randsc`, `randsct` | 0.08 | 0.21 | 0.43 | 0.83 | 1.5 | 1.9 | 2.1 |
| i tre mosaici | 0.06 | 0.14 | 0.28 | 0.53 | 0.96 | 1.2 | 1.4 |

Quindi **da `0.3` a `0.6` è dove compare una banda a una velocità di uno**, che
è ciò che mancava a questi campi. Contare i valori distinti direbbe che `0.01`
basta — novecento valori — ma mille valori dentro un centesimo di banda sono
tutti un solo colore.

**Quanto costa.** Scalare spostava il modulo del valore, e un bailout rotondo
guarda esattamente quello: una cella il cui livello stava vicino a dove cade il
bailout usciva tagliata in due, e i mosaici perdevano le loro forme. Una
rotazione lascia il modulo dov'è, quindi **con un bailout circolare
l'inclinazione è gratuita a qualunque intensità** — su 90000 pixel, a 0.05, 0.3
e 1, e su tutti e cinque i campi: nemmeno un pixel esce diversamente.

Un poligono di bailout è un'altra cosa, perché non legge il modulo ma le
componenti, e una rotazione porta il valore lungo un cerchio che può
attraversare un lato. Di quanto può attraversarlo è quanto gli angoli della
forma sporgono oltre i suoi lati, e i numeri seguono questo — la frazione
dell'immagine la cui uscita cambia, sui cinque campi:

| forma del bailout | angoli su apotema | inclinazione 0.01 | inclinazione 0.05 | inclinazione 0.3 | inclinazione 1 |
| --- | --- | --- | --- | --- | --- |
| cerchio | 1.00 | **0%** | **0%** | **0%** | **0%** |
| esagono | 1.15 | **0%** | 0.0–0.8% | 0.7–2.2% | 0.4–2.4% |
| quadrato | 1.41 | **0%** | 0.0–1.0% | 0.8–2.8% | 1.2–12.8% |
| triangolo | 2.00 | 0.05–1.3% | 0.7–5.0% | 2.6–6.4% | 7.4–21.6% |

Anche il **numero** a cui è impostato il bailout lo sposta, e per la stessa
ragione: ciò che può cambiare lato è un valore il cui modulo cade tra l'apotema
del poligono e i suoi angoli, e dove quell'anello si trova tra le celle è ciò
che il numero decide.

Quindi con un poligono i due desideri tirano in direzioni opposte, e non c'è
modo di sistemarlo diversamente: l'uscita legge le due componenti del valore e
così fa la colorazione, quindi qualunque cosa dia al colore qualcosa da leggere
è qualcosa che l'uscita può vedere. **Abbina l'inclinazione a un bailout
circolare** e non costa nulla a qualunque intensità. Con un poligono, o accetti
quei pochi punti percentuali o tieni l'inclinazione a `0.05` e alzi invece la
velocità di colorazione.

Ciò che una rotazione non può toccare è **`zmag`**, che legge il modulo. Colora
con qualcosa che legga le due componenti separatamente: `iter + real`,
`iter + imag`, `angle`, `real / imag` all'esterno, e `real` o `real / imag`
nella colorazione interna.

Ne seguono due cose. Il valore è complesso quando l'inclinazione non è zero,
così `randsc(7,,,,,0.4)*z` ruota `z` oltre a scalarlo. E il caleidoscopio
ripiega in entrambi i casi: misurato su cinque campi, due, tre, cinque e sei
spicchi ed entrambi gli specchi, la rotazione di uno spicchio lascia ogni
valore dov'era.

**`randscq(...)`** — lo stesso campo senza l'interpolazione: un mosaico di
celle quadrate piatte invece di macchie. Stessi argomenti, stesso significato.

**`randscp(...)`** — di nuovo lo stesso campo, tolte le curve ma non
l'irregolarità: un seme è sparso dentro ogni cella e ogni posizione prende il
valore del seme più vicino, il che disegna un diagramma di Voronoi — poligoni
convessi piatti a lati dritti, non due della stessa forma.

**`randsch(...)`** — lo stesso ancora su esagoni: un alveare di celle piatte.
Dei poligoni regolari che tassellano il piano è quello senza una direzione
privilegiata — ogni cella ha gli stessi sei vicini agli stessi sei angoli —
quindi si legge come un materiale e non come una griglia.

**`randsct(...)`** — triangoli equilateri, di orientamento alterno. Un mosaico
triangolare ha una direzione privilegiata, ed è ciò che si chiede scegliendolo.

**`randsctile(tiling, seed, ...)`** — lo stesso campo su una qualunque di
**quarantacinque tassellature**, che sceglie il primo argomento; tutto ciò che
segue è ciò che prende la famiglia, spostato di un posto, quindi `selfsim` è il
nono. Ogni tassellatura è scalata a una tessera di area unitaria in media,
quindi il primo argomento cambia la forma delle celle e non la loro scala. Un
numero fuori da 1–45 disegna zero.

| | |
| --- | --- |
| 1–3 | le tassellature regolari: quadrati, triangoli, esagoni |
| 4–11 | le otto archimedee, poligoni regolari con tutti i vertici uguali: 4.8.8 (ottagoni e quadrati), 3.6.3.6, 3.4.6.4, 3.12.12, 4.6.12, 3.3.3.4.4, 3.3.4.3.4 (quadrato camuso), 3.3.3.3.6 (esagonale camusa) |
| 12–19 | le loro duali, che non sono regolari: tetrakis quadrata, rombille (i cubi), deltoidale triesagonale (aquiloni), triakis triangolare, kisrombille, e i pentagoni prismatici, del Cairo e a fiore |
| 20–32 | mattoni, opera fiamminga, spina di pesce, intreccio a cesto, pitagorica (quadrati di due dimensioni), chevron, quadrati e rombi, casette (un pentagono, dritto e capovolto a turno), file di quadrati e triangoli in due ritmi, esagoni tra triangoli, croci greche, tetramini a T |
| 33–37 | stelle islamiche: a otto punte con croci, a sei punte con esagoni, a otto punte da 4.8.8, a dodici punte da 3.12.12 e da 4.6.12 |
| 38 | celle di Voronoi, irregolari — quelle di `randscp`, con un valore proprio |
| 39–43 | tassellature che non si ripetono mai, dalle multigriglie di de Bruijn: i rombi di Penrose, gli aquiloni e le frecce di Penrose, i quadrati e rombi di Ammann–Beenker, e rombi in dodici e in sette direzioni |
| 44–45 | per sostituzione: la girandola, i cui triangoli sono rivolti in ogni direzione possibile, e la sedia, L tagliate in L |

La **scheda Tassellature** di Aiuto → Riferimento delle formule utente ha
un'immagine di ciascuna, con la didascalia come la elenca la scheda Valori, e
un doppio clic su una copia negli appunti l'inizio della sua chiamata —
`randsctile(18, ` — pronto per il seme. Le immagini sono disegnate da
`randsctile` stessa, otto unità per lato attorno all'origine, da
`tools/randsctile-thumbnails.cpp`, e conservate nell'eseguibile invece di
essere disegnate ogni volta che si apre la scheda. Le immagini conservate
possono diventare obsolete, quindi ognuna è salvata con un'impronta di ciò che
la funzione disegnava, e un test ridisegna gli stessi punti e fallisce, dicendo
quale tassellatura, finché
`cmake --build <build directory> --target randsctile-thumbnails` non l'ha
ridisegnata.

Quelle periodiche sono **tabelle**, scritte da `tools/randsctile-tables.py`,
che costruisce ogni tassellatura dalla sua geometria e la verifica — ventimila
punti a caso, ognuno dei quali deve cadere in esattamente una tessera — prima
di scriverla. Le stelle islamiche sono ricavate dalle tassellature archimedee:
ogni poligono di otto lati o più diventa una stella con le punte sui vertici
del poligono, e il sottile triangolo tra ogni rientranza della stella e il lato
va alla tessera al di là di quel lato, così due stelle che condividono un lato
formano un rombo tra loro e un piccolo poligono prende una rientranza su ogni
lato che condivide con una stella. Le tessere possono essere concave — stelle,
croci, la T — quindi un punto viene collocato nella sua tessera con il numero
di attraversamenti anziché lato per lato. Il punto viene diviso in periodi
interi nella precisione di lavoro, come `randscq` lo divide in celle, quindi
una tassellatura periodica è esatta comunque lontano.

Quelle che non si ripetono mai **non hanno tabella** e vengono calcolate da
capo in ogni punto. I rombi di Penrose, Ammann–Beenker e i rombi a dodici e a
sette direzioni vengono dalla **multigriglia di de Bruijn**: famiglie di linee
parallele in un secondo piano, una tessera dovunque due si incrociano, nominata
dalle due linee — quattro numeri interi, esatti per chiunque li chieda. Il
punto viene collocato provando gli incroci vicini a dove cade in quel piano,
prima le famiglie le cui linee passano più vicine, il che lo trova in
quattro-otto tentativi in media. Gli aquiloni e le frecce vengono dai triangoli
di Robinson tagliati livello per livello da una ruota di dieci attorno
all'origine, la girandola da un triangolo rettangolo tagliato in cinque, la
sedia da una L tagliata in quattro; a ogni livello una linea o due decidono in
quale figlio sta il punto.

**Quanto costa**, misurato per chiamata contro i 220 ns di `randscq` su questa
macchina: le tassellature periodiche da 255 a 380 ns, quelle che non si
ripetono mai da 410 a 680. Sono calcolate in `double` in entrambe le build,
compreso il posto entro un periodo — in long double quelle che non si ripetono
mai impiegavano il triplo — e così le due build disegnavano tutte le
quarantacinque tassellature uguali pixel per pixel, quarantamila punti
ciascuna. Oltre 2³² celle dall'origine, dove `double` comincerebbe a perdere le
loro tessere, le tassellature che non si ripetono mai diventano piatte come fa
la famiglia oltre la sua griglia; quelle periodiche resistono fin dove resiste
il resto della famiglia.

L'inclinazione misura una tessera come `randscp` misura il suo poligono,
rispetto al raggio del cerchio più grande che la tessera contiene, quindi i
contorni sono la tessera rimpicciolita, comprese quelle concave; la rosetta
ruota attorno al centro della tessera.

    z+randsctile(18,13,{0.3,0.3},{1,1})*0.25        pentagoni del Cairo
    z+randsctile(36,13,{0.3,0.3},{1,1})*0.25        stelle a dodici punte
    z+randsctile(39,13,{0.3,0.3},{1,1},5)*0.25      Penrose, ripiegato cinque volte

**`fbm(value, seed, [intensity=4], [frequency=8], [octaves=4], [roughness=0.5])`**
— un moto browniano frazionario: rumore a gradiente sommato in ottave, ognuna a
frequenza doppia della precedente e con `roughness` della sua altezza — lo
stesso moto che disegnano le modalità di colorazione con quel nome, vedi [Moto
browniano frazionario](#moto-browniano-frazionario). Ciò che disegna è
**usura** più che un motivo — macchie e ammaccature — perché nessuna ottava è
abbastanza grande da vedersi da sola e nessuna abbastanza piccola da sparire.

Dove `randsc` legge la posizione e può leggere solo la posizione, **questa
legge qualunque cosa le scrivi davanti**. È tutto il senso di averla come
funzione:

    fbm(z,7)          si muove con l'orbita
    fbm(x,7)          resta ferma sul piano
    fbm(z*3+c,7)      qualunque cosa sia
    z^2+c+fbm(z,7,0.1)*i    l'insieme con un bordo ruvido

Il valore e il seme sono obbligatori; gli altri hanno valori predefiniti. Il
moto va **da zero a `intensity` e mai sotto**, quindi sommarlo a qualcosa non
può portarlo sotto zero. Il seme si legge come lo legge `randsc`, quindi lo
stesso numero indica lo stesso campo in entrambe. `octaves` è tenuto tra uno e
ventiquattro, e conviene alzarlo solo quando si alza `roughness`: con la metà
semplice l'ottava ottava porta un duecentocinquantesimo del totale.

**Non ha un caleidoscopio suo**, a differenza della famiglia `randsc`, perché
ce n'è già uno con cui comporla: `parchment(z, sectors)` ripiega l'angolo di
`z` in `|sectors|` settori e `parchmenta` rispecchia le metà alterne, così
`fbm(parchmenta(z,6),7)` lo dice, e lo dice dove chi legge la formula può
vederlo.

Le modalità di colorazione sotto **Frattale → Colorazione esterna → Moto
browniano frazionario...** applicano lo stesso moto a un colore. Questa lo
mette dove una formula può usarlo.

| | |
| --- | --- |
| `randsc` | morbido, curvo |
| `randscq` | netto, regolare, quadrato |
| `randscp` | netto, irregolare |
| `randsch` | netto, regolare, esagonale |
| `randsct` | netto, regolare, triangolare |
| `randsctile` | netto, una qualunque di quarantacinque tassellature |

Tutte e sei dispongono una cella su ogni quadrato unitario della dimensione in
vigore, quindi `size` significa la stessa cosa ovunque: cambiare una lettera
cambia la forma delle celle e non la scala dell'immagine. Un esagono di
circoraggio uno copre 2.6 quadrati unitari e un triangolo di lato uno ne copre
0.43, quindi quelle due griglie sono scalate per corrispondere alle altre
anziché lasciate come vengono.

Ogni passo ha un campo suo. L'iterazione viene mescolata nell'hash insieme alla
posizione, quindi una formula che chiama una di queste una volta per iterazione
ottiene ogni volta un valore nuovo — anche con una degradazione di uno, che
lascia la dimensione com'è. La degradazione stabilisce quanto sono grandi le
macchie in un dato passo; non è ciò che rende diverso il passo.

C'è un limite inferiore alla dimensione. Una cella si trova dividendo la
posizione per la dimensione e il risultato deve cadere in un intero, quindi la
dimensione non può utilmente scendere sotto circa la posizione divisa per ciò
che contiene un intero — un miliardesimo di miliardesimo di essa. Una
degradazione di un mezzo ci arriva in una sessantina di passi. Oltre non resta
struttura di celle da risolvere e il campo è un unico valore piatto su tutto il
piano; cambia comunque a ogni passo e differisce comunque tra le sei funzioni,
quindi una formula che ne sottrae una da un'altra non si ferma su zero iterando
fino al limite per niente. Per una dissolvenza che resti un'immagine fino in
fondo, usa una degradazione vicina a uno.

Tutte fanno l'hash della posizione del punto e dell'iterazione, mai di `z`, e
mai di uno stato globale. Quindi la stessa immagine torna a ogni ridisegno, con
qualunque numero di thread, e in modalità Mandelbrot o Julia allo stesso modo.
`randsc` è continuo, quindi le due precisioni concordano fino a circa 1e-19; i
mosaici sono funzioni a gradini e possono discordare su un filo lungo i bordi
delle celle, il che è inevitabile quando si chiedono bordi netti.

Una formula che chiama una qualunque di esse disattiva il tracciamento del
bordo. Quell'ottimizzazione riempie una regione che ha trovato circondata da un
solo colore senza calcolarla — vero per un frattale, falso per un campo di
rumore.

Il vecchio `rand` è invariato e dipende ancora da quante volte è stato
chiamato, quindi la stessa posizione non si ridisegna allo stesso modo.
Preferisci `randsc`.

### Caleidoscopi

I due argomenti dopo la degradazione ripiegano il piano prima che il campo vi
venga campionato. Il primo è in quanti spicchi ripiegarlo, con `1` — il
predefinito — che lo lascia com'è; il secondo, `kaled_mirror`, è quale specchio
fa il ripiegamento:

| kaled_mirror | |
| --- | --- |
| `0` (e qualunque altro valore) | la metà lontana di ogni spicchio riflette quella vicina, quindi ogni spicchio è simmetrico rispetto alla propria bisettrice |
| `1` | lo stesso al contrario, con la metà vicina che riflette quella lontana |

    randsc(13,{0.6,0.6},{1,1},6,0)     sei spicchi, ciascuno specchio di sé stesso
    randsc(13,{0.6,0.6},{1,1},3,1)     tre spicchi, ripiegati nell'altro verso

Entrambi i ripiegamenti sono continui dove gli spicchi si incontrano, quindi il
rumore resta coerente e le due precisioni continuano a concordare; un
ripiegamento che incontrasse sé stesso in modo irregolare si vedrebbe come una
cucitura.

È l'unica parte della famiglia che costa trigonometria, e solo un valore di due
o più la raggiunge. Una chiamata che non ne dice nulla paga un confronto.

### Ottenere un'immagine dal rumore

`randsc` restituisce un valore in `[0, 1)`. Usato da solo, l'iterazione non
lascia mai un bailout di 4, quindi nulla esce e l'immagine è piatta.
Moltiplicalo, o lascia che perturbi un'iterazione che esce davvero:

    z+randsc(13,{0.35,0.35},{1,1})*0.25       macchie
    randsc(13,{0.15,0.15},{0.5,0.5})*2.5      moto browniano
    randscq(13,{0.15,0.15},{0.5,0.5})*2.5     quadrati sparsi, sempre più piccoli
    z+randscp(13,{0.3,0.3},{1,1})*0.25        poligoni irregolari
    z+randsch(13,{0.3,0.3},{1,1})*0.25        un favo
    z+randsct(13,{0.3,0.3},{1,1})*0.25        triangoli
    z^2+c+randsc(13,{0.25,0.25},{1,1})*1.2    l'insieme stesso deformato

Una degradazione di 0.5 dimezza le macchie a ogni passata, quindi dopo venti
iterazioni sono un milionesimo della loro dimensione e sotto il pixel. Per una
dissolvenza lenta su una lunga corsa, usa qualcosa vicino a 1 — `0.97^60` è
ancora 0.16.

Un seme intero è esatto in entrambe le compilazioni; uno frazionario viene
quantizzato ai suoi bit iniziali, che concordano tranne che per circa un seme
su cinquanta milioni.

### Tutte le passate finora: `selfsim`

Ogni passata di queste funzioni è un campo a sé — l'iterazione entra nell'hash
— con celle pari alla degradazione per quelle della passata precedente. Una
formula che ne chiama una a ogni passata legge quindi ogni volta un campo nuovo
e più fine, e una modalità di colorazione che legge l'ultima passata legge solo
l'ultimo campo. A una degradazione di 0.75 le celle di una vista larga cinque
sono sotto il pixel entro la sedicesima passata, e oltre l'interno è neve:
cento passate di `randsct(672,,{0.75,0.75},,,{2,2},6)+c` in modalità Julia
mettono il 70% dei pixel interni vicini a un sedicesimo della tavolozza o più
l'uno dall'altro.

L'ottavo argomento restituisce invece **le passate finora, mediate**, che
stanno l'una rispetto all'altra come le ottave di un moto browniano
frazionario: la passata *n* pesa d^(nH), dove *d* è la degradazione — la media
geometrica delle sue due componenti prese senza segno, quindi `0.75+0.75i` vale
0.75 — e *H* è l'argomento. I pesi sono divisi per il loro totale, quindi la
risposta resta tra i valori mediati: in [0, 1] senza inclinazione.

| selfsim | |
| --- | --- |
| omesso | spento — l'unica passata, fino all'ultimo bit. Ciò che riceve una chiamata che non nomina nulla, e un posto lasciato vuoto |
| `0` | la media semplice: ogni passata pesa uguale, ed è ciò che dà ogni *H* vicino a zero — non spento |
| `1` | il moto semplice: ogni passata pesa *d* volte la precedente |
| `0.5` | più ruvido: le passate fini conservano più del loro peso, quindi si vede più dettaglio |
| più grande | più liscio: le prime passate si prendono quasi tutto |
| negativo | le passate fini pesano più di quelle grossolane, ed ecco di nuovo la neve |
| complesso, `{1,2}` | la parte reale pesa come sopra, la parte immaginaria ruota le passate — vedi sotto |

    randsct(672,,{0.75,0.75},,,{2,2},6,1)+c     lo stesso, autosimile

Su quella formula la neve sparisce: i pixel interni vicini escono in media a
1.8 celle di tavolozza di distanza invece che a 35, e lo 0.6% di essi a più di
sedici invece del 70%. Ciò che resta sono i triangoli a ogni scala insieme, con
i grandi che portano i piccoli.

* **`selfsim` può essere complesso.** Allora anche d^(nH) è complesso: d^(n·Hr)
  in grandezza e una rotazione di n·Hi·ln d, quindi ogni passata è ruotata di
  tanto in più della precedente e poi mediata come la media un *H* reale — una
  spirale di ottave. A 0.75 una parte immaginaria di uno fa circa sedici gradi
  a passata. I pesi sono divisi per **la somma delle loro grandezze**, non per
  la loro somma complessa: così la risposta resta entro il più grande dei
  valori mediati. La somma complessa si avvicinerebbe a zero a qualche passata
  e farebbe esplodere ogni valore a quella passata per tutti i pixel insieme,
  dato che i pesi non hanno nulla a che fare con il pixel. Le rotazioni danno
  al valore una parte immaginaria anche senza inclinazione, quindi `imag` e
  `angle` hanno qualcosa da leggere. Una degradazione di uno ha logaritmo zero
  e non ruota nulla; `{1,0}` è `1` fino all'ultimo bit; e `{0,2}` è ogni
  passata pesata uguale, e ruotata.
* **Zero non è spento.** d^(n·0) vale uno per ogni passata, quindi `0` è la
  media semplice, e lo è ogni *H* vicino — da entrambi i lati, e anche
  complesso. All'inizio lo zero era l'interruttore, il che ne faceva l'unico
  valore per cui la curva non passa: `0.000001` era la media semplice e `0` la
  sola ultima passata, a un quarto dell'intervallo di distanza. La sola ultima
  passata è l'estremo lontano della curva, con *H* che corre a meno infinito.
  Per non avere media, ometti l'argomento, o lascia vuoto il suo posto.
* **È il valore inclinato a essere mediato**, quindi l'inclinazione e ogni bit
  di `skew_mode` continuano a fare ciò che facevano, e la risposta è complessa
  quando l'inclinazione la rende tale.
* **Una media avvicina i valori.** A 0.75 e *H* uguale a uno lo scarto è circa
  due quinti di quello di una passata, quindi alza la velocità di colorazione
  per riavere il contrasto.
* **Una degradazione di uno** pesa ogni passata uguale e dà la loro media
  semplice, che si assesta verso un centro piatto man mano che le passate si
  sommano. È ciò che dice la definizione, ed è fatto come detto.
* **Le passate che una chiamata non ha visto vengono calcolate quando viene
  interrogata.** Una chiamata su un ramo che la formula prende solo dalla
  settima passata riceve comunque le prime sette, e la risposta a una passata è
  la stessa qualunque strada vi sia arrivata. Un pixel nuovo si riconosce dalla
  posizione oltre che dalla passata: la sola passata, che è ciò su cui si basa
  `trap`, permetterebbe a una tale chiamata di continuare la media lasciata dal
  pixel precedente.
* **Nulla trabocca.** La media è tenuta come una sola, con ogni passata che ne
  prende una quota che segue dalla quota della passata precedente, invece che
  come somma e totale — quindi un *H* negativo, o una degradazione sopra uno, i
  cui pesi supererebbero ciò che il tipo contiene, non ha dove farlo.
* **Costa** circa un decimo in più a chiamata a 64 bit di mantissa e un settimo
  a 113. Spento, niente che si possa misurare.

### Figure invece di rumore

`sierpinskyt`, `sierpinskyc` e `snowflake` stanno dove sta il rumore, e non
hanno nulla di casuale: la figura è dove sembra essere, e la modalità
Mandelbrot e la modalità Julia disegnano la stessa immagine, perché una forma
nel piano sta nello stesso posto in entrambi i casi.

**Sono frattali a sé, non campi da moltiplicare dentro uno.** Scritta da sola —
`sierpinskyt()` e nient'altro — ognuna disegna la sua figura, come Frattale →
Altre formule disegna il suo **Sierpinski**, il **Sierpinski Carpet** e il
**Koch Snowflake**.

Tutte e tre lo fanno **portando il punto al suo genitore**. Ogni parte di una
di queste figure ne ha uno: un buco della guarnizione di Sierpinski sta dentro
un buco più grande un livello sopra, un buco del tappeto dentro la cella che è
stata tagliata allo stesso modo, un triangolo del fiocco di neve su un lato
dell'esagono al centro o sul bordo libero di un altro triangolo. Il passo sul
genitore è un movimento del piano — un raddoppio allontanandosi dall'angolo più
vicino per la guarnizione, l'ingrandimento di una cella per il numero di celle
per il tappeto. La parte più alta non ha genitore dentro la figura, quindi il
suo passo la porta fuori dal bailout, e una parte *n* livelli più giù impiega
*n* passi per arrivarci: **la passata in cui un punto esce è il livello a cui
si trova**, e il conteggio delle iterazioni è l'immagine.

Ciò che torna indietro nella passata in cui un punto esce conta quanto la
passata stessa: **ogni modalità di colorazione esterna tranne il conteggio
delle iterazioni lo legge**. Quindi la parte più alta di ogni figura viene
lanciata fuori invece di essere semplicemente dichiarata uscita — la
guarnizione raddoppia il suo buco centrale allontanandolo da un angolo e lo fa
atterrare fuori dal triangolo, il tappeto scala la sua cella centrale rispetto
alla cella accanto, il fiocco di neve lancia fuori il suo esagono centrale
nella direzione verso cui guarda quel settore. Ognuno porta con sé da dove
veniva, quindi `real`, `imag`, `angle`, la colorazione graduale e il resto
hanno tutti qualcosa su cui lavorare. Solo a un punto che non è mai stato nella
figura viene dato un numero che sta per «uscito».

Un fiocco di neve si legge **dall'esagono verso l'esterno**, ed è la cosa da
sapere su di esso. Si scompone esattamente in un esagono regolare al centro,
sei triangoli del lato stesso di quell'esagono poggiati sui suoi sei lati,
dodici di un terzo di quello sui loro bordi liberi, quarantotto di un nono sui
loro, e così via per sempre. Tre dei sei sono gli angoli del triangolo da cui
la figura è cresciuta e tre sono le prime protuberanze messe sui suoi lati, e
**nulla li distingue**: dall'esagono verso l'esterno un fiocco di neve ha
simmetria sei.

Letto invece dal primo triangolo, cinque ottavi della figura sono il livello
uno — un unico triangolo piatto che riempie l'immagine, con il fiocco di neve
solo nella frangia intorno, e nulla da dire per il colore esterno. Letto
dall'esagono, i primi due livelli sono cinque dodicesimi ciascuno e gli altri
scendono in modo uniforme. Un triangolo poggiato sull'esagono passa sul
triangolo più grande che l'esagono contiene — angolo su angolo, tre volte
l'area — quindi viene ingrandito della radice quadrata di tre e ruotato di un
dodicesimo di giro; uno poggiato sul bordo libero di un altro triangolo viene
ingrandito di tre e ruotato verso la direzione in cui guarda quel bordo.

**Ogni passo espande**, ed è una regola più che un caso. La guarnizione
raddoppia, il tappeto ingrandisce una cella per il numero di celle, il fiocco
di neve per tre o per la radice di tre. Un triangolo poggiato sull'esagono
potrebbe invece essere ribaltato sul lato su cui poggia, il che lo fa atterrare
su uno dei sei di cui è fatto l'esagono ed è più ordinato; ma un ribaltamento è
un'isometria, e un'isometria non ha sensibilità per dove il punto è partito.
Scritta dentro una formula più grande, la figura darebbe allora a quella
formula un'immagine rigida di sé stessa su cinque dodicesimi della sua area, e
la formula vi disegnerebbe una regione piatta. Un test verifica le tre
espansioni.

Il punto viaggia davvero, ed è ciò che le rende utilizzabili e non solo
corrette: `z` lungo la strada è un punto del piano come un altro, quindi le
modalità di colorazione che lo leggono, la colorazione graduale e lo scrivere
una figura dentro una formula più grande hanno tutti un significato.

Una guarnizione e un tappeto riempiono la loro forma, quindi ogni punto in uno
di essi ha un livello ed esce a quel livello. Un fiocco di neve non riempie il
suo esagono: lascia sei angoli di **sfondo**, e quello sfondo non fa parte
della figura e non ha livello. Viene **ruotato di mezzo giro**, il che porta
sfondo su sfondo, quindi non esce mai e risulta nel **colore interno** — lo
spazio vuoto resta vuoto invece di prendersi una banda del colore esterno
togliendola alla figura.

Ruotato, e non restituito dov'era, perché **un punto fisso è un'orbita ferma e
non un'orbita limitata**. Restando fermo, il ciclo eseguiva ogni passata fino
al limite senza che nulla cambiasse; le modalità di colorazione interna che
confrontano una passata con la precedente non avevano nulla da confrontare; e
`snowflake()` scritta dentro una formula più grande restituiva a quella formula
il punto che già teneva, mentre `sierpinskyt()` e `sierpinskyc()` muovono ogni
punto che ricevono. Il mezzo giro è preferito al sesto di giro, che la
simmetria della figura permetterebbe ugualmente: entrambi portano sfondo su
sfondo, ma il mezzo giro tiene il punto dentro un quadrato e un cerchio oltre
che dentro l'esagono. Lascia invariata la distanza dal centro, quindi le
modalità di colorazione interna che la leggono non ne sono toccate, e
l'immagine che disegna `snowflake()` è invariata al pixel.

È l'unico punto in cui queste figure costano qualcosa. A un pixel di sfondo
viene posta la domanda di appartenenza una volta per passata per tutto il tempo
in cui lo si guarda, mentre le altre due figure e il corpo del fiocco di neve
escono alla passata che il loro livello indica. Disegnato direttamente nel
ciclo di iterazione, un fiocco di neve intero costa circa 2.4 volte quanto
costava quando allo sfondo era data una banda propria. Il tracciamento del
bordo resta attivo per queste formule — le loro bande sono davvero piene, che è
la condizione per usarlo — e lo sfondo è un'unica regione piena, quindi ciò che
arriva sullo schermo è meno di quanto suggerisca quel numero.

Scritte dentro una formula più grande — `snowflake()^-2+snowflake()` e simili —
le tre **non** si comportano allo stesso modo, e vale la pena dirne chiaramente
la ragione perché nessuna quantità di lavoro sulla figura la cambierà.

Una guarnizione è una **polvere**: non ha area, i suoi punti non escono mai, e
ognuno di essi viene raddoppiato a ogni passata. Una formula costruita su di
essa ha sotto qualcosa di caotico ovunque, a ogni profondità. Un fiocco di neve
è **pieno**: dieci dodicesimi sono l'esagono e i sei triangoli, che escono alla
prima e alla seconda passata, quindi una formula costruita su di esso esaurisce
la figura quasi subito e disegna ciò che disegna da sola. Il dettaglio sta
lungo il bordo di Koch, dove i livelli vanno in profondità, e le regioni piatte
sono gli interni, dove non ci vanno.

**Ciò che una passata sa è ciò che quella passata ha pagato.** Il fiocco di
neve decide se un punto è sfondo scendendo lungo la curva di Koch, e a quella
discesa è concesso **un livello per passata**. Prima eseguiva tutti i suoi
ventiquattro livelli in una sola valutazione, quindi una formula scritta
intorno alla figura riceveva l'intero bordo alla prima passata e lo disegnava
per intero per quante poche passate le venissero date — a `Iterazioni: 2` si
otteneva ancora ogni protuberanza della curva di Koch, cosa che nessun'altra
formula di XaoS fa e che non è lo scopo di un numero di iterazioni. Esaurire i
livelli non equivale a essere fuori: un punto fuori dalla cornice o sopra
l'involucro è sfondo e si sa già che lo è, mentre uno non deciso viene portato
avanti come parte della figura e interrogato di nuovo alla passata successiva
con un livello in più da spendere. L'immagine che disegna `snowflake()` da sola
è invariata — un punto nella figura non ha mai avuto bisogno della discesa
profonda, solo della certezza di non essere sfondo.

Ciò che resta è la differenza onesta tra le due figure, e nessun lavoro
sull'una o sull'altra la colmerà. Una guarnizione è una **polvere**: nessuna
area, i suoi punti non escono mai, ognuno raddoppiato a ogni passata, quindi
una formula costruita su di essa ha sotto qualcosa di caotico a ogni
profondità. Un fiocco di neve è **pieno**, e dieci dodicesimi sono l'esagono e
i sei triangoli, che escono alla prima e alla seconda passata. Misurato a un
solo zoom — la frazione di pixel che stanno sul bordo di una banda —
`sierpinskyt()^-2+sierpinskyt()` arriva al 17 per cento entro dodici passate e
`snowflake()^-2+snowflake()` al 2.5; entrambi crescono con il numero di
passate, che è il punto, e dove si fermano è ciò che le figure sono.

`radius` significa ciò che significa `bailout` e si legge allo stesso modo:
**come quadrato della distanza**. Ciò che disegna è la forma che disegna un
bailout di quel numero, **inscritta in essa esattamente** — imposta il bailout
sulla forma corrispondente e sullo stesso numero, e la figura e il bailout si
sovrappongono:

| | |
| --- | --- |
| `sierpinskyt` | il triangolo del bailout **triangolo −90°**, angolo per angolo |
| `sierpinskyc` | il quadrato del bailout **quadrato**, lato per lato |
| `snowflake` | l'esagono del bailout **esagono 0°**, con le sue sei punte sui sei angoli dell'esagono |

La cosa da sapere è che un poligono di bailout pone i suoi **lati** alla radice
quadrata del bailout dal centro — è il suo apotema, non il raggio circoscritto
— quindi i suoi angoli sono più lontani di quanto dica il numero: il doppio per
un triangolo, e di un settimo per un esagono. `sierpinskyt(4)` arriva quindi a
4 con gli angoli mentre `sierpinskyc(4)` arriva a 2 con i lati, e sono giusti
entrambi.

Imposta il bailout più grande del raggio e la figura non riempie più la forma:
c'è un **margine** tra le due, e impostalo più piccolo e la figura viene
tagliata. Entrambe le cose sono ciò che è stato chiesto e nessuna viene
corretta. Ma il margine non fa parte della figura e non ha livello, quindi non
gliene viene dato uno: viene **ruotato secondo la simmetria propria della
figura** — un terzo di giro per una guarnizione, un quarto per un tappeto,
mezzo giro per un fiocco di neve — il che lo tiene fuori dalla figura e dentro
il bailout, quindi non esce mai ed è disegnato nel **colore interno** con le
modalità di colorazione interna che lavorano su di esso. Prima usciva alla
prima passata e prendeva la banda più bassa del colore esterno, il che
disegnava l'intero margine in un unico tono piatto e toglieva una banda alla
figura.

Una rotazione e non uno specchio, anche se uno specchio servirebbe altrettanto
a tenerlo dov'è: ciò che una figura restituisce deve ruotare con il punto, così
che ruotare l'immagine ruoti ciò che vi è disegnato, e uno specchio lo ruota
nel verso opposto. Per la stessa ragione il centro di un tappeto — il blocco di
celle che l'anello del bordo racchiude, che non ha genitore e quindi deve
essere lanciato fuori — viene lanciato nella direzione verso cui guarda il suo
**quarto** del quadrato, i quarti tagliati dalle diagonali. Lanciato sempre
nella stessa direzione, com'era, disegnava quel centro come una rampa da un
lato all'altro con nulla di traverso; il quadrato ha simmetria quattro e così
dovrebbe avere ciò che vi è disegnato. Un test verifica che ogni figura
risponda a una rotazione del proprio ordine con la stessa rotazione.

**`sierpinskyt([radius=4], [kaleidoscope=1], [kaled_mirror=0])`** — la
guarnizione di Sierpinski, in quel triangolo.

**`sierpinskyc([radius=4], [squares=3], [kaleidoscope=1], [kaled_mirror=0])`**
— il tappeto di Sierpinski, in quel quadrato. Il quadrato è tagliato in
`squares` per `squares`, l'anello di celle lungo il bordo viene tenuto, tutto
ciò che quell'anello racchiude viene scartato, e lo stesso viene fatto a ogni
cella tenuta. Quindi l'immagine è **un quadrato al centro e 4×`squares`−4
intorno**: otto a tre, che è il tappeto come lo si disegna di solito, dodici a
quattro, sedici a cinque. Due è l'unico numero senza un vero anello, e lì se ne
va invece l'angolo lontano, il che è di nuovo una guarnizione — un quadrato
tagliato in quattro con un angolo tolto è proprio ciò che è una guarnizione.

**`snowflake([radius=4], [kaleidoscope=1], [kaled_mirror=0])`** — il fiocco di
neve di Koch, con le punte sugli angoli di quell'esagono, a bande per
generazione dal centro verso l'esterno: l'esagono alla prima passata, i sei
triangoli sui suoi lati alla seconda, i dodici sui loro bordi liberi alla
terza, i quarantotto alla quarta. I sei angoli di sfondo che la figura non
copre non hanno alcuna banda — non escono mai, e sono disegnati nel colore
interno.

Ogni argomento ha un valore predefinito, quindi `snowflake()` è una chiamata, e
lo è anche `sierpinskyc( ,5)` — un tappeto più traforato alla dimensione
predefinita.

**Il caleidoscopio è quello della famiglia del rumore**, scritto per ultimo in
ognuna delle tre e con lo stesso significato: in quanti spicchi viene tagliato
il piano, e quale specchio li ripiega. `1` — ciò che riceve una chiamata che
non dice nulla — non ripiega nulla e lascia la figura il numero che è sempre
stata, fino all'ultimo bit.

Ciò che ripiega è però diverso, e deve esserlo. Il rumore ripiega la
**posizione**, che resta ferma per tutta l'orbita; una figura non sa che
farsene della posizione e legge **z**, il punto che sta portando giù verso la
parte in cui si trova, quindi è z a essere ripiegato. L'immagine risulta allo
stesso modo — la figura disegnata in uno spicchio e ripetuta intorno
all'origine — e regge fino in fondo, perché ogni passo è fatto sul punto
ripiegato e ripiegato di nuovo alla passata successiva. Dentro il mezzo
spicchio che il ripiegamento lascia com'è, la figura è la figura.

    sierpinskyt(4,6,0)      la guarnizione sei volte intorno all'origine
    sierpinskyc(4,3,5,1)    il tappeto in cinque spicchi, lo specchio vicino
    snowflake(4,5)          cinque spicchi di fiocco di neve

Un fiocco di neve ha già la simmetria sei con lo specchio lungo ogni
bisettrice, quindi ripiegarlo in sei spicchi in quel modo è l'identità
sull'immagine: `snowflake(4,6,0)` esce **pixel per pixel** uguale a
`snowflake(4)` su 90000 pixel. Chiedi cinque spicchi, o quattro, per vederlo
ripiegare.

Una cosa da sapere sul tappeto: ciò che taglia via lo lancia in una di quattro
direzioni, i quarti tagliati dalle diagonali del quadrato, e sulla diagonale
stessa la scelta è un pareggio. Un ripiegamento può portare un intero raggio
del piano su quella diagonale, e lì i due lati del pareggio vengono lanciati in
direzioni opposte — una linea sottile lungo la diagonale, che è il gradino che
il tappeto ha sempre avuto e non qualcosa che fa il ripiegamento. Misurato su
diecimila punti a sei spicchi: ogni disaccordo stava su una diagonale, e i
punti ripiegati stessi concordavano fino a una parte su 10^16.

    sierpinskyt()        con Frattale -> Forma della fuga -> triangolo -90
    sierpinskyc()        con forma del bailout quadrato
    sierpinskyc( ,5)     un tappeto più traforato
    snowflake()          con forma del bailout esagono 0

Nessuna costa più del rumore accanto, passata per passata: misurate rispetto a
`randsc`, la guarnizione 0.44, il tappeto 0.54 e il fiocco di neve 0.50 di esso
a 64 bit di mantissa, e 0.50, 0.87 e 0.65 a 113. Leggere un fiocco di neve dal
centro lo ha reso più economico oltre che più bello da vedere: sei settimi
della figura sono l'esagono o un triangolo poggiato su di esso, e a quelli si
risponde con quattro moltiplicazioni e un confronto senza percorrere affatto la
curva. Ciò che un fiocco di neve costa su un'immagine intera è il paragrafo
sopra, e riguarda quante passate prende lo sfondo più che quanto costa una
passata. Non era gratis, e a 113 bit era quasi perso — una radice quadrata e
una divisione sono entrambe software lì, e ciascuna costa circa quanto una
figura intera. Quindi il reciproco della radice quadrata del raggio viene
tenuto sul punto di chiamata invece di essere ricalcolato a ogni passata, e la
radice stessa torna moltiplicando il raggio per quel reciproco invece che
dividendo; la guarnizione confronta i suoi tre pesi baricentrici scalati per un
numero positivo, il che non cambia né quale sia il più grande né quale sia
negativo; il tappeto conta in celle invece che nel piano, dove la
moltiplicazione per `squares` e la divisione per esso si annullano; e il fiocco
di neve prende il figlio su cui passare dal lato sotto la cui protuberanza è
stato trovato, che doveva comunque trovare. Percorre inoltre la curva di Koch
in `double`, poiché la figura ha una dimensione fissa e 52 bit di una frazione
di un lato sono più di quanti ne possano usare i 24 livelli che disegna. Un
test verifica i rapporti, quindi una riscrittura che li perde lo dice.

## Tavolozze

Frattale → Tavolozza chiede un numero di algoritmo, e dopo il numero dice come
si chiama ciascuno. Da 1 a 3 sono quelli di XaoS; da 4 a 7 sono stati
ricostruiti dopo aver misurato di cosa sono fatti i primi tre, e l'8 dopo aver
misurato dei dipinti.

| | |
| --- | --- |
| 1 | **Scuro e colore** — un quasi nero con una tinta propria e un colore, alternati |
| 2 | **Nero, colore, bianco** |
| 3 | **Caldi su verde scuro** — viola, rossi e arancioni accesi su verdi scuri e verde acqua, nero o bianco ogni terzo |
| 4 | **Smog** — scura e malinconica: grigi, ocra, ardesia e viola polveroso, con sangue di bue e scarlatto in mezzo |
| 5 | **Caldi sulla notte** — rosa, rosso, bordeaux, marrone, arancione, ocra e giallo su verde acqua, petrolio, blu e indaco, nello scheletro della 3 |
| 6 | **Coppie preferite** — le coppie di colori che i primi tre mettono più spesso fianco a fianco, poste intere tra nero e bianco |
| 7 | **Colori casuali** — ogni tappa un colore davvero casuale, e nient'altro |
| 8 | **Kandinsky** — i gradienti dei suoi dipinti, dall'opaco al vivido e ogni tanto quasi al grigio |

### Di cosa sono fatti i primi tre

Ventimila tavolozze di ciascuno, generate come le genera il programma, ogni
colore chiamato col nome del più vicino fra tre dozzine di colori di
riferimento in CIELAB, contando i colori e le coppie di colori:

* **Ancore.** Nero e bianco sono da un terzo a metà di ogni tappa, e i colori
  stanno tra loro con un periodo di due o tre tappe, così lo schema si vede
  nelle prime bande — che, in una tavolozza lunga circa tremila tappe, sono
  tutto ciò che la maggior parte delle immagini usa.
* **L'1 e il 2 contengono 256 colori e non di più.** Prendono il byte basso del
  generatore, e il byte basso di quel generatore è da solo un generatore di
  periodo 256: il rosso decide verde e blu, e ogni colore è sempre seguito
  dallo stesso. Su ventimila tavolozze del 2, nessun colore è stato seguito
  anche una sola volta da uno diverso. Da qui vengono le loro coppie ricorrenti
  — rosa accanto a giallo ventiquattro volte più spesso del caso, bordeaux
  accanto a oliva sei, rosso accanto a petrolio otto.
* **Il 3 attinge da due finestre**: colori accesi dal viola attraverso il rosso
  fino all'arancione, e scuri dal verde fino al verde acqua. Le sue coppie più
  comuni sono un colore caldo acceso su un verde scuro.

### I quattro dopo di loro

Ognuno è stato misurato rispetto ai primi tre prima di restare: la divergenza
di Jensen-Shannon tra i colori che due algoritmi mostrano, per nome, e tra le
coppie di vicini che mostrano, su cinquemila tavolozze — 0 per la stessa
tavolozza, 1 per nulla in comune. Con questa misura l'1 e il 2 distano 0.11,
l'1 e il 3 0.16.

* **Il 4, Smog**, è scuro e spento dove l'1 è vivido sul nero. Le tappe
  profonde sono un colore appena sopra il nero, quelle chiare un grigio coperto
  con una tinta e mai bianco, e ogni quarta tappa un grigio semplice. Le sue
  tinte sono percorse da un anello di 256 come l'1 e il 2 percorrono le loro,
  quindi i vicini ricorrono come ricorrono i loro, ma solo sulle tinte dello
  smog — ocra e oliva, e acciaio, ardesia e un viola polveroso — e, per un
  quarto dell'anello, dal cremisi allo scarlatto: sangue di bue quando è
  profondo e uno scarlatto smorzato quando è chiaro, tenuto saturo, perché un
  rosso scuro che non lo è si legge come marrone. Una tappa su cinque si legge
  come rosso, contro una su trenta nell'1. Dista da 0.43 a 0.56 dai primi tre.
  All'inizio era il meccanismo dell'1 su un altro anello di 256 colori, scelto
  il più lontano possibile dall'anello classico, e misurava 0.05 dall'1 — più
  vicino del 2: la somiglianza era il meccanismo, nero metà delle volte e un
  colore da qualunque parte l'altra metà, non i colori.
* **Il 5, Caldi sulla notte**, è lo scheletro del 3 con la temperatura
  invertita: colori caldi di ogni luminosità, dal rosa attraverso rosso,
  bordeaux, marrone e arancione fino a ocra e giallo, su colori notturni tenuti
  scuri, dal verde acqua attraverso petrolio e blu fino all'indaco. Entrambe le
  finestre sono larghe 125 gradi, vicino ai 135 del 3; all'inizio erano un
  terzo di tanto e si leggevano come due colori. Il blu marina è il suo colore
  più comune, e blu marina con oliva, bordeaux o marrone le sue coppie più
  comuni.
* **Il 6, Coppie preferite**, prende ogni coppia di colori che uno dei primi
  tre mette uno dopo l'altro almeno tre volte più spesso del caso —
  sessantasei, con i due colori come sono usciti davvero — e pone ogni coppia
  intera, i due colori a contatto, tra una tappa nera e una bianca.
* **Il 7, Colori casuali**, è ogni tappa un colore estratto a caso dalla parte
  alta del generatore, dove l'1 e il 2 estraggono dalla bassa: settantanovemila
  colori diversi in ventimila tavolozze dove il 2 ne ha 256, e nessuna coppia
  più comune del caso di oltre un quinto. All'inizio aveva nero e bianco ogni
  terza tappa, e misurava 0.00 dal 2: il ritmo di molto scuro e molto chiaro
  era tutto ciò che si vedeva. Senza, ha meno escursione dallo scuro al chiaro
  di quelli ancorati, che è ciò che vuol dire casuale.

I quattro che sostituiscono erano uno spettro, un bicromo, una triade e una
coppia complementare, ragionati a partire dal cerchio dei colori, e sullo
schermo ciascuno era un solo colore. Distribuivano le loro tinte lungo tutta la
tavolozza, e la tavolozza nel programma è lunga circa tremila tappe — da
quattro a nove nel test che doveva accorgersene, ed ecco perché non se n'è
accorto — quindi un'immagine vedeva solo la prima tinta. E tenevano ogni tappa
un colore, profondo o pallido, senza nero e senza bianco, quindi nulla in esse
aveva un bordo. Il test ora genera le sue tavolozze come fa il programma, e
chiede a tutti e sette che la parte che un'immagine vede vada dallo scuro al
chiaro e contenga più di un colore; all'8, venuto dopo, chiede un po' meno, e
il perché è più sotto.

Una posizione salvata con i numeri da 4 a 7 prima di questo cambiamento torna
con colori diversi: il numero e il seme sono tutto ciò che una posizione
registra della sua tavolozza, e ciò che il numero significa è cambiato. Una
salvata con i numeri da 1 a 3 non ne è toccata, fino all'ultimo bit, e una
salvata con i numeri da 4 a 7 nomina un algoritmo che lo XaoS originale non ha
e rifiuterà.

### 8, Kandinsky

Fatta da due serie di immagini nel suo stile: quattro vignette vivide — nero,
vermiglione, giallo, ceruleo, viola e rosa su carta crema — e un'immagine
opaca, arancione bruciato, ocra, marrone, salvia e petrolio sulla stessa carta.
Le vignette danno diciotto colori, con k-means in CIELAB, e ognuno ha una
**forma opaca**: il colore dell'immagine opaca più vicino in luminosità e
tinta, croma esclusa, perché la croma è ciò in cui le due differiscono. Il
vermiglione diventa arancione bruciato, il giallo ocra, il bordeaux marrone, il
nero un marrone scuro, il ceruleo salvia, il blu petrolio. L'immagine opaca non
ha grigio, quindi il grigio resta grigio. Il verde acqua che sembra avere è un
salvia, 128,146,121: si legge come verde acqua solo accanto all'arancione.

Ogni tavolozza estrae un **umore**. Quattro su cinque cadono in modo uniforme
tra le forme opache e quelle vivide; la quinta va oltre le forme opache verso
il grigio, fino a quindici centesimi della loro croma, e una tavolozza su venti
circa è quasi grigia.

**Le ancore arrivano per caso.** Una tavolozza su dieci può trasformare una
tappa qua e là in nero puro o bianco puro: ogni tappa dopo la prima con un
proprio lancio, a una frequenza che la tavolozza estrae tra una tappa su venti
e una su dieci, nero o bianco secondo come cade il lancio — nessun ritmo e
nessuna alternanza, a differenza delle ancore da 1 a 3. Le altre nove tavolozze
su dieci sono esattamente ciò che erano prima che arrivassero le ancore: il
lancio era già fatto per ogni tappa, e lasciato inutilizzato.

Ciò che segue un colore è dove lo portano **i suoi gradienti**. Dovunque il
colore delle vignette cambia dolcemente — su otto pixel, senza un passo più
grande di un terzo del cambiamento — sono stati contati i colori ai due
estremi, e ogni tappa estrae la successiva da quei conteggi: la carta prosegue
in una velatura di giallo e poi nel giallo, il nero nel vermiglione (tre volte
su dieci) e il vermiglione nel bordeaux, il blu nel malva e il malva nel rosa,
il ceruleo in un verde acqua grigio. Anche i bordi netti sono stati contati, e
lasciati fuori: la maggior parte sono carta accanto a carta di un'altra
sfumatura, e le tavolozze fatte con essi venivano carta e grigio. Una tavolozza
comincia sulla carta o sul nero.

Su ventimila tavolozze, quelle del decimo più alto di umori distano 0.06 dalle
vignette, con la stessa divergenza di sopra, e quelle del decimo sopra lo zero
0.07 dall'immagine opaca; dalle altre sette tavolozze distano da 0.39 a 0.67.

Poiché un gradiente resta per un po' in una parte del cerchio dei colori, e l'8
non ha nero e bianco ogni poche tappe su cui ripiegare, i test gli chiedono
meno che agli altri. Dove la tavolozza è stesa lunga, quaranta voci per tappa e
più, un'immagine vede una decina di tappe, e una tavolozza su cento è allora un
unico campo — ocra con dentro salvia, petrolio con vermiglione — che non arriva
a un terzo della strada dal nero al bianco né a un secondo colore. E tre tappe,
che è ciò che spesso è una tavolozza corta di 256 voci, sono un unico
gradiente.

## Osservare l'orbita

`trap` e `stripe` tengono un numero sull'intera orbita invece che sul punto,
che è il tipo di grandezza che le modalità di colorazione del motore non
possono raggiungere: lì un ciclo di calcolo è compilato per restituire un
colore, e solo l'ultimo `z` e quello prima sopravvivono per essere usati nella
colorazione. Fatto nella formula, al motore non costa nulla.

Entrambe restituiscono il valore ricevuto fino all'ultima passata che il limite
di iterazioni permette, e in quella ciò che hanno raccolto. Quindi una formula
intera è

    trap(z^2+c, 3)

e nient'altro: il frattale itera come farebbe, e all'ultima passata il valore
diventa la trappola, che poi disegnano le modalità di colorazione interna. Un
punto che esce non arriva mai a quella passata e conserva il suo normale colore
esterno, quindi ciò che queste disegnano è l'interno.

**`trap(orbit, shape, centre, size)`** — il punto più vicino a cui l'orbita sia
mai arrivata a una forma. `shape` è predefinito a 0, `centre` all'origine,
`size` a 1.

| shape | |
| --- | --- |
| `0` | il centro stesso, un punto |
| `1` | una retta orizzontale che lo attraversa |
| `2` | una verticale |
| `3` | entrambe, una croce |
| `4` | un anello di raggio `size` |
| `5` | un quadrato di semilato `size` |
| `6` | un rombo di semidiagonale `size` |

**`stripe(orbit, density)`** — la media di `(sin(density * arg orbit) + 1) / 2`
lungo l'orbita. `density` è predefinito a 4 ed è quante strisce fanno un giro;
un numero intero, altrimenti non si incontrano dove il giro si chiude. Una
media sull'orbita cambia dolcemente con il punto anche dove il conteggio delle
iterazioni salta, ed è ciò che disegna le fibre per cui il metodo è noto.

    trap(z^2+c,0)                    quanto vicino all'origine è passata l'orbita
    trap(z^2+c,3)                    una croce, che disegna raggi
    trap(z^2+c,4,{0,0},{0.5,0})      anelli dentro la cardioide
    stripe(z^2+c,6)                  sei strisce per giro

La grandezza accumulata vive sul punto di chiamata, quindi due trappole in una
formula tengono ciascuna la propria e un thread non può disturbarne un altro.
Entrambe disattivano il tracciamento del bordo, per la stessa ragione del
rumore: due vicini che impiegano lo stesso numero di passate possono comunque
aver visto orbite molto diverse, e una regione riempita senza essere calcolata
sarebbe sbagliata.

Queste non sono modalità di colorazione, e non possono esserlo: le modalità
vivono dentro il ciclo compilato di ogni formula, che restituisce direttamente
un colore, e ce ne sono trentuno, ciascuno con la propria tabella di simmetria
per posizione. Ciò che c'è qui funziona solo per una formula utente.

## Forma della fuga

Calcolo → Raggio di fuga → Forma della fuga. L'iterazione si ferma quando `z`
lascia una regione, che è sempre stata un cerchio. La forma si vede nelle bande
fuori dall'insieme, quindi è uno strumento di disegno quanto uno numerico.

| forma | esce quando |
| --- | --- |
| Cerchio (classico) | \|re\|² + \|im\|² supera il bailout |
| Quadrato | lo fa una delle due componenti da sola |
| Rombo | lo fa \|re\| + \|im\| |
| Asse reale | la sola componente reale |
| Asse immaginario | la sola componente immaginaria |
| Entrambi gli assi | entrambe le componenti insieme |
| Triangolo 0° / 90° / 180° / −90° | il punto lascia un triangolo con quell'orientamento |
| Esagono 0° / 90° | il punto lascia un esagono con quell'orientamento |
| Ottagono | il punto lascia un ottagono |

I poligoni sono misurati dal loro **apotema** — la distanza dal centro al punto
medio di un lato — quindi un poligono e il cerchio dello stesso bailout si
toccano lungo i lati e non negli angoli. L'orientamento fa parte della forma:
un triangolo ruotato di 120° è lo stesso triangolo, quindi 0°, 90°, 180° e −90°
sono quattro triangoli diversi, mentre un esagono si ripete ogni 60° e quindi
ne ha due, e un ottagono ogni 45°, il che non gli lascia nulla da scegliere. Il
triangolo e l'esagono hanno quindi un sottomenu proprio per l'orientamento. Le
loro normali sono calcolate una volta per fotogramma, quindi il ciclo interno
costa da tre a otto moltiplicazioni-addizioni invece di un seno e un coseno per
lato.

Rispetto al cerchio, il rombo cambia circa il 39% dei pixel e l'asse reale il
44%. La sagoma dell'insieme non si sposta — la decide ciò che non esce mai —
quindi guarda le bande.

Viene scritta in una posizione salvata solo quando non è il cerchio, quindi una
posizione che non chiede una forma resta caricabile da qualunque versione
precedente.

## Moto browniano frazionario

Un dialogo su ogni lato — **Frattale → Colorazione interna → Moto browniano
frazionario...** e lo stesso sotto Colorazione esterna — che contiene le tre
modalità che offre e i cinque numeri da cui sono disegnate. Era un sottomenu,
con le tre modalità come voci a scelta esclusiva e i numeri dietro un dialogo
Impostazioni sotto di esse; ora è un dialogo unico, spuntato quando una delle
tre modalità è in uso, e OK imposta la modalità e i numeri insieme come un solo
passo dell'annullamento. Una posizione salva ancora la modalità come
`(incoloring N)` o `(outcoloring N)` e i numeri a parte, come `(infbmset ...)`
o `(outfbmset ...)`, quindi un file di una versione precedente si apre
invariato.

Ciò che fa è rendere **usata** una colorazione pulita: le bande mantengono la
loro forma e il loro ordine e smettono di essere perfette, come se la cosa
fosse stata lasciata alle intemperie. Non è casualità pixel per pixel, che si
legge come grana ed è stata provata e scartata; è un campo di rumore con una
scala, quindi ciò che disegna sono segni.

**All'esterno** le tre modalità sono `fbm + smooth`, `fbm + iter` e `fbm` da
solo — il conteggio graduale consumato, il conteggio a bande consumato, e il
rumore da solo, che disegna nuvole intorno all'insieme senza alcun conteggio
dentro.

**All'interno** sono `fbm`, `fbm + zmag` e `fbm + decomposition`. La prima è
quella interessante: l'interno è per impostazione predefinita un unico tono
piatto, e questa gli dà una superficie.

### Rumore a gradiente, e nessuna griglia

Il moto è **rumore a gradiente**: ogni angolo di un reticolo riceve una
direzione, e un punto prende da ogni angolo la pendenza che quella direzione dà
alla sua distanza. Prima era rumore di valore — un numero a ogni angolo,
raccordato da una smoothstep — e quello disegnava una griglia. Una smoothstep
non ha pendenza a nessuna delle due estremità, quindi il campo diventava piatto
lungo ogni linea del reticolo, e poiché il reticolo di ogni ottava conteneva le
linee della prima, tutte diventavano piatte lì insieme. Attraverso una linea la
pendenza era zero; l'immagine usciva a quadrati, otto pixel di lato in
`FBM_ERROR.xpf`, come se fosse stata calcolata a una risoluzione più bassa di
quella a cui era mostrata. Il rumore a gradiente è ripido sulle linee quanto
tra di esse — misurato, 0.99 alle impostazioni predefinite e 1.08 a quelle di
quella posizione — e il reticolo di ogni ottava è spostato rispetto al
precedente, di quantità che nessun raddoppio rende mai intere, quindi due
ottave non condividono mai una linea né un angolo. `fbm()` nelle formule
disegna lo stesso moto, e `randsc` ne è un'ottava.

Ciò che è stato mantenuto: il valore va ancora da zero all'intensità e mai
fuori; i **segni hanno la dimensione che avevano**, con il reticolo a tre
quinti della frequenza perché il rumore a gradiente li disegni grandi quanto li
disegnava il rumore di valore (misurato da quanto si sposta il campo prima di
somigliare a sé stesso solo la metà: 0.56 di cella alle impostazioni
predefinite, prima e dopo); il **contrasto**, di cui il rumore a gradiente da
solo ha solo sette decimi — ogni ottava passa per una curva che porta i suoi
valori sulla distribuzione che aveva il rumore di valore, monotona perché non
sposti alcun contorno, e misurato sul moto alle impostazioni predefinite
oscilla di 0.133 dal suo centro contro 0.132 di prima; e il costo, più o meno —
le modalità un po' più veloci di prima, `fbm()` circa un ottavo più lenta. Ciò
che non lo è stato: una posizione salvata con una modalità fbm conserva i suoi
numeri, la sua scala e la sua forza e **disegna segni diversi**.

### I cinque numeri

Lo stesso dialogo della modalità, sotto di essa. Sono salvati con la posizione.

| | |
| --- | --- |
| **intensità** | di quante bande di colore il rumore sposta il valore. `4` è il predefinito; sotto uno è un sussurro, sopra dieci le bande smettono di essere bande |
| **frequenza** | la dimensione dei segni, grandi quanto li disegnerebbe un reticolo di altrettante celle per unità del piano. **È quella da alzare man mano che zoomi**: alla prima vista `8` si legge bene, con un'ampiezza di `0.02` ne servono circa `200` |
| **ottave** | quante se ne sommano, ognuna al doppio della frequenza della precedente |
| **ruvidità** | quanto ogni ottava conserva dell'altezza della precedente. `0.5` è il moto semplice; più alto è più granuloso |
| **seme** | la stessa immagine ogni volta. Cambialo per averne un'altra dello stesso carattere |

**Ottave e ruvidità si sostengono a vicenda.** Con una ruvidità di `0.5`
l'ottava ottava porta un duecentocinquantesimo del totale, quindi le ottave
oltre la quarta o la quinta non cambiano nulla che si possa vedere — misurato,
quattro e diciotto davano la stessa immagine. Alza la ruvidità e le ottave fini
hanno qualcosa da spendere, e allora il numero comincia a contare.

### Dimensionare i segni con il conteggio

All'esterno, due campi in più in fondo al dialogo. **Dimensionamento del
rumore** dice che cosa decide la grandezza dei segni: **fisso**, dove parte, la
sola frequenza, come sopra; **iterazione**, la frequenza alzata con il
conteggio delle iterazioni; **iterazione smooth**, la frequenza alzata con il
conteggio graduale. **Fattore** dice di quanto: la frequenza raddoppia quel
numero di volte ogni dieci iterazioni, quindi a `1` i segni sono grandi la metà
dove l'orbita ha impiegato dieci iterazioni in più a uscire, e un quarto dove
ne ha impiegate venti. A `0` il conteggio non cambia nulla. Dove il conteggio
sale, verso il bordo dell'insieme, i segni si fanno più fini e seguono il
dettaglio.

Una frequenza che cambia da un pixel all'altro non si può semplicemente
stendere sul piano: il campo si stirerebbe, tanto più quanto più ci si
allontana dall'origine. Così il moto è letto alle due potenze di due a cavallo
del conteggio, ciascuna un campo a sé a frequenza fissa, e i due sono mescolati
come una texture tra le sue mipmap, riportando la miscela al contrasto di un
campo solo. Con **iterazione smooth** i segni cambiano grandezza senza
cuciture; con **iterazione** cambiano ai bordi delle bande, ogni banda con una
grana sua. All'interno il conteggio è lo stesso ovunque, quindi l'interno non
ha nessuno dei due campi.

Tra due raddoppi si leggono due campi, e la colorazione costa circa il doppio
che con **fisso**: 248 ns per pixel contro 127, a quattro ottave. Dove il
conteggio cade esattamente su un raddoppio, come fa ogni conteggio intero a un
fattore di `10`, si legge un campo solo e il costo è quello di **fisso**. Con
**iterazione smooth**, `fbm + iter` e `fbm` sono disegnati dal ciclo smooth,
che fa qualche iterazione oltre il bailout per calcolare il conteggio graduale.
Un fattore troppo alto porta i segni sotto il pixel, dove diventano grana.

Con **fisso** disegna come prima, bit per bit. Una posizione salva i due come
`(outfbmscale 'smooth 1.5)`, e solo quando sono stati spostati, quindi una
posizione che non li tocca si apre ancora nelle versioni precedenti.

### Quanto costa, e l'unica cosa che non fa

Il rumore è calcolato **una volta per pixel**, nella colorazione, non una volta
per iterazione. Su un'immagine di 90000 punti con un'ampiezza di `0.02` e 400
iterazioni, `fbm + smooth` impiegava un quinto in più del semplice `smooth` a
quattro ottave e un quarto in più a diciotto; a 4000 iterazioni, dove l'orbita
fa più lavoro, le stesse quattro ottave costavano un decimo.

I segni sono ancorati **al piano**, non alla vista: è sporco sulla cosa, non
sull'obiettivo. Quindi zoomare li ingrandisce, e se vai abbastanza in fondo
finisci dentro un solo segno e l'immagine torna pulita. Alzare la frequenza è
la risposta, ed è per questo che la frequenza è un'impostazione e non una
costante. Non la si può invece legare allo zoom: il motore riusa i pixel
calcolati alla scala precedente, e una colorazione che si muovesse con la vista
lascerebbe loro i segni sbagliati.

Le modalità sono numerate **dopo** quelle a colori reali, quindi ogni numero di
modalità che una posizione salvata contiene già significa ancora ciò che
significava.

## Tassellazione

**Frattale → Colorazione interna → Tassellazione...** colora l'interno con le
tessere di una delle quarantacinque tassellature di `randsctile`, posate sulla
`z` finale, quaranta per unità come i quadrati di `squares`. Il suo dialogo ha
quattro campi.

| campo | |
| --- | --- |
| **Tassellazione** | quale tassellatura, con il numero con cui la conosce `randsctile` e il suo nome: **3: esagoni**, dove parte. La scheda Tassellature di Aiuto → Riferimento delle formule utente ha un'immagine di ciascuna |
| **Colorazione** | da cosa viene il colore, dappertutto: **angolo di z**, dove parte; **angolo attorno al centro della tessera**, quello della tessera in cui cade il punto, che fa di ogni tessera una piccola ruota di colore; **zmag**, il valore della modalità con quel nome; **nessuna**, solo le tessere |
| **Tessere** | come ogni tessera si stacca dalle vicine: **ognuna spostata di un tanto proprio**, lungo la tavolozza, dove parte; **sfumate dal centro al bordo**, lo stesso, e spostata di un ottavo dei colori di un giro in più dal centro al bordo, come una smussatura; **con il bordo marcato**, lo stesso, e la parte esterna della tessera disegnata nell'altra metà della tavolozza; **un colore pieno ciascuna**, un colore suo, da qualunque cosa venga il colore |
| **Intensità** | quanto sono distanti lungo la tavolozza le tessere: moltiplica ogni quantità propria di una tessera, il suo spostamento, la smussatura di **sfumate dal centro al bordo** e il bordo di **con il bordo marcato**. A `1`, dove parte, le tessere sono come descritte qui; a `0` resta solo ciò da cui viene il colore, così **angolo di z** è una decomposizione semplice e **zmag** lo zmag semplice; in mezzo le tessere si vedono appena, e sopra `1` di più. Misurato con **Vista → Sonda della tavolozza (clic su un punto)**, a `1` una tessera di **angolo di z** e **ognuna spostata di un tanto proprio** si sposta fino a una decina di posti, così `0.1` la sposta fino a uno e `0.2` fino a due; con **un colore pieno ciascuna**, **zmag** o **nessuna** si sposta fino all'intera tavolozza, 8191 posti in quella predefinita, così lì un posto vale circa `0.0001` |

La modalità sta tra le modalità che si ripetono e le decomposizioni, e le
impostazioni di entrambe la raggiungono. **Frequenza (modalità periodiche)** e
**Funzione (modalità periodiche)** raggiungono le tessere come raggiungono ogni
argomento periodico dell'interno: a frequenza 2 le tessere sono grandi la metà.
Le [impostazioni delle decomposizioni](#impostazioni-delle-decomposizioni)
dell'interno raggiungono i due angoli, che sono tagliati e percorsi come quelli
delle decomposizioni. **Settori**, **Rotazione**, **Ripetizioni**, **Verso** e
**Velocità** raggiungono entrambi; **Angolo di** raggiunge solo l'angolo di
`z`, perché quello attorno al centro della tessera non sa che farsene; **Anelli
per unità** fa anelli del modulo per l'angolo di `z`, attraverso la frequenza e
la funzione come per le decomposizioni, e per l'angolo attorno al centro
altrettanti anelli dal centro al bordo, che seguono la forma della tessera.
**zmag**, **nessuna** e una tessera piena non ne leggono nessuna.

La quantità propria di una tessera viene dall'hash che le dà il nome. Per un
angolo copre tanta tavolozza quanta ne copre un giro dell'angolo, così
un'impostazione della decomposizione spostata lascia i colori delle tessere
dove erano; per **zmag**, **nessuna** e una tessera piena copre tutta la
tavolozza, così **nessuna** con **ognuna spostata di un tanto proprio** è lo
stesso di **un colore pieno ciascuna**. **angolo attorno al centro della
tessera** mette un giro intero di colori in ogni tessera, che alla grandezza
delle tessere sono moltissimi colori in pochi pixel: una **Velocità** sotto
uno, o una frequenza sotto uno, rende leggibili le ruote.

Con **Colorazione** a **angolo di z**, **Tessere** a **ognuna spostata di un
tanto proprio** e le impostazioni delle decomposizioni dove partono, la
modalità disegna quello che disegnava prima di avere una scelta, bit per bit.
Le tassellature che non si ripetono mai diventano piatte, senza una quantità
propria, oltre quattro miliardi di tessere dall'origine. Questo vale con
**Intensità** a `1`. A `0` le tessere si vedono ancora con **angolo attorno al
centro della tessera**, il cui angolo è misurato dal centro di ogni tessera.

Le tassellature periodiche costano circa quanto le altre modalità interne, da
120 a 210 ns per pixel; quelle che non si ripetono mai, che vengono localizzate
da capo in ogni punto, da 250 a 490. **angolo attorno al centro della
tessera**, **sfumate dal centro al bordo** e **con il bordo marcato** misurano
dove sta il punto nella sua tessera, il che costa fino a due terzi in più sulle
tassellature periodiche e poco di più sulle altre. Una posizione salva la
modalità come `(incoloring 30)`, la tassellazione come `(intessellation 39)`,
solo quando non è quella degli esagoni, e la colorazione come
`(intessellstyle 'tile-angle 'edged)`, solo quando è stata spostata.
L'intensità si salva come `(intessellintensity 0.5)`, solo quando non è uno.

### All'esterno

**Frattale → Colorazione esterna → Tassellazione...** ha tre modalità, come
Moto browniano frazionario: `tessellation + smooth`, le tessere sopra il
conteggio smooth; `tessellation + iter`, sopra il conteggio intero;
`tessellation`, le tessere da sole. Le tessere sono posate sul posto del pixel
nel piano, non sul punto di fuga, che salta da una fascia all'altra e darebbe
solo grana, come accadde quando ci fu provato il moto.

| campo | |
| --- | --- |
| **Modalità** | una delle tre |
| **Tassellazione**, **Colorazione**, **Tessere**, **Intensità** | come all'interno, con impostazioni proprie. **angolo di z** è l'angolo del punto di fuga e **zmag** il logaritmo del suo modulo, come lo aggiunge `iter+log(mag)`; tranne che in `tessellation + iter` con dimensionamento **fisso** o **iterazione**, si leggono entrambi qualche iterazione oltre il bailout, dove si calcola il conteggio smooth, e sfumano nello stesso punto un'iterazione dopo quanto il conteggio smooth è avanzato verso l'iterazione seguente, così proseguono attraverso i bordi delle iterazioni invece di spezzarsi lì |
| **Frequenza** | tessere per unità del piano, `8` dove parte: la grandezza delle tessere, che va alzata quando si zooma, come la frequenza del moto |
| **Dimensionamento delle tessere** e **Fattore** | cosa dà la grandezza alle tessere, come [Dimensionare i segni con il conteggio](#dimensionare-i-segni-con-il-conteggio) la dà ai segni del moto: **fisso**, la sola frequenza; **iterazione** o **iterazione smooth**, la frequenza raddoppiata tante volte quante dice il fattore ogni dieci iterazioni di quel conteggio |

Tra due raddoppi le tessere sono lette a entrambe le grandezze e i due colori
si fondono, una dissolvenza da una grandezza alla successiva: le tessere non
hanno un valore di mezzo verso cui fondersi, come ha il moto. Una tavolozza di
256 colori, che non ne ha tra i suoi, prende la grandezza più vicina. Tra due
raddoppi la tassellatura è letta due volte, il che ne raddoppia il costo. Un
fattore di `1` raddoppia la frequenza dieci volte ogni cento iterazioni, il che
porta presto le tessere sotto il pixel: qualche decimo è un inizio migliore. Le
[impostazioni delle decomposizioni](#impostazioni-delle-decomposizioni)
dell'esterno raggiungono i due angoli: **Iterazioni per fascia** fa fasce del
conteggio per l'angolo di z e anelli dal centro al bordo per l'angolo attorno
al centro.

Con **Intensità** a `0` e **Colorazione** a **nessuna**, `tessellation + iter`
disegna come `iter` e `tessellation` un solo colore. Una posizione salva la
modalità come da `(outcoloring 21)` a `(outcoloring 23)`, e le impostazioni in
un solo comando, `(outtessellation 5 'tile-angle 'shaded 0.6 40 'smooth 0.4)`,
solo quando una di esse è stata spostata.

## Altre modalità di colorazione

Non è il calcolo a decidere il colore di un pixel: il colore viene da una
formula nella `z` finale, nel parametro `c` e nel conteggio delle iterazioni, e
c'è una formula del genere per ogni modalità. Ventidue in più, in un sottomenu
proprio sotto Frattale → Colorazione interna e Frattale → Colorazione esterna,
sopra i sottomenu a colori reali.

**All'interno**, dove l'orbita non è mai uscita e solo la `z` finale e `c`
dicono qualcosa:

| modalità | |
| --- | --- |
| `real/mag` | il coseno dell'angolo: un'ombreggiatura che si richiude su sé stessa dove `atan2` ha una cucitura |
| `max(\|real\|,\|imag\|)` | la norma quadrata, quindi le bande sono quadrati dove quelle di `zmag` sono cerchi |
| `\|real\|+\|imag\|` | le stesse bande ruotate di un ottavo di giro |
| `min(\|real\|,\|imag\|)` | vicino a zero lungo l'uno o l'altro asse, quindi raggi |
| `\|z-c\|` | quanto lontano dal parametro si è assestata l'orbita, il che distingue i bulbi |
| `\|z*c\|` | i due moduli l'uno contro l'altro |
| `angle(z)-angle(c)` | la direzione di assestamento misurata rispetto al parametro |
| `real*imag` | una sella, quattro lobi intorno agli assi |
| `sin(real)*sin(imag)` | una griglia posata su dove l'orbita si è assestata |
| `sign(imag)` | due toni piatti: una decomposizione dell'interno, come la decomposizione binaria fa con l'esterno |
| `frac(mag)` | curve di livello a ogni ottavo del modulo |
| `log(mag)` | `zmag` con il contrasto spostato sui valori piccoli, dove passa il suo tempo un'orbita attratta |
| `cos(\|z\|)` | anelli attorno a dove si è fermata l'orbita, otto per unità del modulo |
| `hexagons` | tre onde a un terzo di giro l'una dall'altra, che sommate danno un reticolo di esagoni |
| `spiral` | cinque bracci attorno a dove si è fermata l'orbita, avvolti una volta a ogni fattore e del modulo |
| `cos(\|z-c\|)` | anelli attorno al parametro, che distinguono i bulbi come fa `\|z-c\|` |

Undici delle modalità interne si ripetono — `cos(mag)`, `mag*cos(real^2)`,
`sin(real^2-imag^2)`, `squares`, `sin(real)*sin(imag)`, `frac(mag)`, le quattro
in fondo alla tabella e `tessellation` — e **Frequenza (modalità periodiche)**,
il primo campo del dialogo Velocità della colorazione interno, dice quanto
spesso. Moltiplica l'argomento di ciò che si ripete, quindi a 2 i seni girano
il doppio più in fretta, `squares` disegna quadrati grandi la metà e
`frac(mag)` sedici contorni per unità del modulo dove ne disegnava otto. In
`cos(mag)` restringe anche le bande che scelgono tra il seno e il coseno. A
uno, dove parte, disegnano come hanno sempre fatto, bit per bit, e le altre
modalità non lo leggono mai. Una posizione lo salva come `(incolorfrequency 2)`
solo quando non è uno, quindi una posizione che non lo tocca si apre ancora
nelle versioni precedenti.

**Funzione (modalità periodiche)**, sotto di essa, piega l'argomento con una
delle [funzioni di velocità della
colorazione](#funzioni-di-velocità-della-colorazione) prima che la frequenza lo
moltiplichi, così che le bande girino più in fretta dove la funzione sale
ripida e più piano dove si appiattisce: `square` e `cube` le addensano verso il
bordo e le diradano attorno allo zero, `log` e `atan` al contrario, e `sin`,
`cos` e `triangle` ripiegano l'argomento su se stesso, così che le bande
tornino indietro. Raggiunge ogni funzione periodica delle undici, il cui
argomento `a` entra come `k·f(a)`, dove `k` è la frequenza: i seni e i coseni,
il modulo di `frac(mag)`, la parte reale e quella immaginaria le cui strisce
fanno i quadrati di `squares` e le tessere di `tessellation`, e in `cos(mag)`
anche le bande che scelgono tra il seno e il coseno. In `spiral` i cinque
bracci restano cinque, perché un numero di giri non intero lascerebbe una
cucitura dove l'angolo si chiude: la frequenza e la funzione agiscono
sull'avvolgimento. Una funzione periodica dentro l'argomento di un'altra
avrebbe piegato anche il proprio. A `ident`, dove parte, non cambia nulla. Una
posizione la salva come `(incolorperiodicfun 6)`, con il posto della funzione
nell'elenco contando da zero, solo quando non è ident.

**All'esterno**, dove ciò che c'è sono il conteggio e il punto di uscita:

| modalità | |
| --- | --- |
| `iter+angle` | la direzione da cui è uscita, fusa nel conteggio: le bande acquistano una torsione |
| `iter+log(mag)` | di quanto ha oltrepassato il bailout, una sfumatura più economica del logaritmo di un logaritmo |
| `iter+real*imag` | biomorfi resi continui invece che a soglia |
| `max(\|real\|,\|imag\|)` | da quale lato è uscita, senza alcun conteggio |
| `iter banded` | il conteggio ripiegato in otto bande, che mostra la forma degli insiemi di livello invece del loro numero |
| `\|real\|-\|imag\|` | quanto è sbilanciato il punto di uscita |

I numeri delle modalità più vecchie non si spostano, quindi una posizione
salvata continua a caricarsi come prima.

La velocità e lo spostamento della colorazione ora raggiungono smooth e smooth
log, cosa che non avevano mai fatto: quelle due restituivano un pixel proprio
invece di passare per il passo con cui finiscono le altre modalità, quindi
entrambi i controlli non facevano nulla mentre una delle due era scelta. Ai
valori con cui vengono — velocità uno, spostamento nullo — il passo è
l'identità, quindi nessuna immagine salvata prima si sposta.

Smooth e smooth log seguono la forma della fuga. Interpolano tra una passata e
la successiva chiedendo dove l'orbita ha attraversato il bailout, e prima lo
chiedevano al modulo qualunque forma fosse stata scelta — così uscire su un
quadrato e interpolare su un cerchio metteva più di un quinto dei punti
uscenti, e la metà con un rombo, fuori dalla banda a cui appartengono. Ogni
forma ora passa la grandezza che verifica davvero e la soglia con cui la
confronta. Per il cerchio le due cose coincidono, quindi nulla di ciò che era
disegnato prima si sposta.

Funzionano anche in modalità Newton, dove un'orbita converge invece di uscire:
la grandezza è allora di quanto si è mosso l'ultimo passo e la soglia il limite
di convergenza, e l'interpolazione è la stessa domanda posta a un numero che
scende invece di salire. I due frattali di Newton, che non avevano alcuna
variante graduale, ora ne hanno una.

**Convergenza di Newton**, nel menu Calcolo, viene finalmente letta. Il limite
era scritto nell'iterazione come un milionesimo, quindi l'impostazione veniva
offerta, salvata e ricaricata mentre nulla la consultava mai.

**Smooth e smooth log** si sono spostate in quel sottomenu nel menu esterno.
Interpolano tra un'iterazione e la successiva usando di quanto l'orbita ha
oltrepassato il bailout, il che richiede una formula che esca sul bailout —
venti delle trentuno lo fanno, e le altre undici o convergono (le Newton) o si
fermano su un test proprio (le Sierpinski, Koch, il Clock). Per quelle la
modalità non fa nulla, e non è qualcosa che il primo elenco dovrebbe offrire.
La formula utente, che avrebbe sempre potuto averla e non l'aveva, ora ce l'ha.

Ciò che queste modalità non possono fare è qualunque cosa abbia bisogno
dell'orbita invece che della sua fine: una trappola d'orbita, una media a
strisce, un esponente di Lyapunov, la stima della distanza. Queste hanno
bisogno di un valore portato da un'iterazione alla successiva, il che in questo
motore significa una seconda copia compilata di ogni ciclo di calcolo. `trap` e
`stripe` in una formula utente lo fanno invece.

### Impostazioni delle decomposizioni

**Frattale → Colorazione interna → Decomposizione**, e lo stesso sotto
Colorazione esterna, imposta come le decomposizioni tagliano il giro: binary
decomposition e color decomposition all'esterno, `sign(imag)`,
`Decomposition-like` e `fbm + decomposition` all'interno, e gli angoli delle
tassellazioni, interna ed esterna. Con ogni campo dove parte — zero, e uno,
destra e uno gli ultimi tre — disegnano come hanno sempre fatto, bit per bit.

| campo | |
| --- | --- |
| **Settori** | in quanti settori si taglia ogni ripetizione del giro. Zero è quello della modalità: i due semipiani di quelle binarie, il giro intero delle altre. Quelle binarie colorano settori alterni, quindi un numero dispari passa al pari successivo — due settori dello stesso colore si toccherebbero. Quelle a colori diventano altrettanti gradini piatti |
| **Rotazione** | di quanti gradi i settori sono ruotati, in senso antiorario |
| **Angolo di** | solo all'interno: l'angolo di `z`, di `z-c` o di `z*c`. `z-c` centra i settori sul parametro, `z*c` li ruota dell'angolo del parametro |
| **Anelli per unità** (interno) | anelli del modulo di ciò di cui si prende l'angolo, distinti uno sì e uno no: quelle binarie diventano una scacchiera polare. Passano per `k·f(a)` come ogni argomento periodico dell'interno, quindi li raggiungono la frequenza e la funzione delle modalità periodiche |
| **Iterazioni per fascia** (esterno) | fasce di tante iterazioni attraverso i settori: con la binary decomposition, la griglia di raggi esterni ed equipotenziali |
| **Ripetizioni** | quante volte si ripete il giro: i settori, e per quelle a colori la tavolozza, ricominciano in ognuna. Un numero qualsiasi, ma una parte di ripetizione lascia una cucitura dove l'angolo si chiude |
| **Verso** | destra, come sono sempre andate; sinistra, al contrario; ambidestro, avanti e poi indietro in ogni ripetizione — uno specchio, senza cucitura quando le ripetizioni sono intere |
| **Velocità** | per quelle a colori, quanta tavolozza attraversa una ripetizione: sopra uno, più colori in ciascuna, non più ripetizioni. Per quelle binarie, il distacco fra i colori dei settori alterni: a uno come sono sempre stati, a zero un colore solo, a due il doppio — la binary decomposition sposta i settori dispari dal conteggio verso il conteggio rovesciato, `sign(imag)` la seconda tinta dalla prima |

Una posizione le salva come `(indecomposition 8 22.5 'z-c 6 2.5 'both 1.5)` e
`(outdecomposition 6 30 2 3 'left 2)` solo quando sono cambiate, quindi una
posizione che non le tocca si apre ancora nella 2.5. Cambiate, una
decomposizione non è più simmetrica rispetto all'asse reale, e XaoS smette di
disegnarne una metà come specchio dell'altra.

## Funzioni di velocità della colorazione

**Frattale → Colorazione interna → Velocità della colorazione**, e lo stesso
sotto Colorazione esterna, è un unico dialogo per tutto ciò che porta il valore
consegnato da una modalità di colorazione a un posto nella tavolozza. La
**funzione**, scelta da un elenco, piega il valore prima che si applichino la
**velocità** e lo **spostamento** della colorazione: la tavolozza viene letta a
`f(v) × speed + shift`. In una colonna propria, accanto, c'è il [cambio di
velocità dopo ITER](#cambio-di-velocità-dopo-iter), che agisce sul valore prima
che la funzione lo riceva. Quello interno si apre con la frequenza e la
funzione delle modalità che si ripetono, che agiscono sulla modalità stessa;
vedi [Altre modalità di colorazione](#altre-modalità-di-colorazione). OK
imposta tutto come un solo passo dell'annullamento. Per le modalità che contano
— iter, smooth e le altre — `v` è il conteggio delle iterazioni stesso, quindi
una funzione si può pensare in iterazioni. Ognuna è dispari, `f(−v) = −f(v)`,
quindi una grandezza che cambia segno lo conserva.

Erano quattro voci del menu: la velocità, lo spostamento, un sottomenu delle
funzioni e il dialogo proprio del cambio di velocità. Una posizione viene
scritta e letta esattamente come prima, e i comandi che quelle voci eseguivano
ci sono ancora per uno script.

`sin`, `cos` e `triangle` girano: la velocità della colorazione, applicata dopo
di esse, dice quante voci della tavolozza copre un giro, e lascia com'era
quanto spesso il giro ritorna. Questo lo fa la **velocità iniziale** del cambio
di velocità, applicata prima: con entrambi i gruppi spenti e l'esponente
iniziale a uno la funzione riceve `s0·v`, quindi 2 avvicina i giri il doppio e
0,5 li allontana il doppio.

| | |
| --- | --- |
| `ident`, `log`, `atan`, `sin`, `cos`, `square`, `cube`, `square root`, `cube root` | ciò che c'è sempre stato, invariato fino all'ultimo bit |
| `log log` | `log(1 + log(1 + v))`, per gli zoom profondi dove il solo log stringe ancora troppo le bande: a un milione di iterazioni log dà 13.8 bande e questa 2.7 |
| `exp` | `e^v − 1`, al contrario del log: apre i valori piccoli. Per le grandezze che restano piccole — le componenti che leggono i colori reali e le decomposizioni; a un conteggio di novanta supera già 10³⁸ |
| `triangle` | il seno con le curve raddrizzate: lo stesso periodo e la stessa altezza, la tavolozza percorsa avanti e indietro a passo costante invece di indugiare nelle svolte |

Le nuove vengono dopo le vecchie, quindi ogni numero che una posizione salvata
contiene nomina ancora ciò che nominava.

### Cambio di velocità dopo ITER

La parte inferiore dello stesso dialogo non è una funzione ma una modifica a
qualunque funzione sia scelta: fino a un'iterazione che indichi, la funzione
continua esattamente com'era; oltre, il conteggio che riceve corre più veloce o
più lento, quindi i colori accelerano o rallentano da lì e non prima. Con
`square root` scelta e ITER a 40, i colori sono quelli della radice quadrata
fino a 40 iterazioni e la radice quadrata di un conteggio che prende velocità
dopo. Contiene un esponente e una velocità iniziali e poi due gruppi, ognuno
separato da una linea:

| | |
| --- | --- |
| **Esponente iniziale** | l'esponente dall'inizio al primo ITER di un gruppo attivo: il conteggio diventa la semplice potenza `v^p0`, la stessa che i gruppi siano attivi o no, e con entrambi spenti vale per ogni conteggio. A uno, dove parte, il conteggio è intatto fino al primo ITER |
| **Velocità iniziale** | ciò per cui si moltiplica il passo di quel tratto: `s0·v^p0`. A uno, dove parte, non cambia nulla |
| **ITER** | dove comincia il tratto del gruppo, nelle unità che la funzione riceve. A zero o sotto comincia all'inizio — il che rende `ident` una potenza con qualunque esponente |
| **Esponente** | l'esponente del tratto, dal suo ITER al successivo o oltre: come si piega il tratto. Parte dal conteggio a cui era arrivato quello prima, `f(ITER)`, e aggiunge `s·ITER·((v/ITER)^p − 1)`, il cui passo a ITER è `s·p`: sopra uno i colori accelerano da lì, sotto uno rallentano |
| **Velocità** | ciò per cui si moltiplica il passo del tratto — lineare dove l'esponente è una potenza: a un esponente di uno il tratto è la retta `f(ITER) + s·(v − ITER)`. A uno, dove parte, non cambia nulla |
| **Modalità** | **off**, **jump** o **smooth**: se il gruppo è usato, e come il suo tratto si raccorda con quello precedente |
| **Raggio** | solo in smooth, e mostrato solo allora: quante iterazioni prende il raccordo da ciascun lato di ITER. Qualunque numero da zero, che è il jump |

Con entrambi i gruppi attivi e il primo ITER sotto il secondo, il conteggio è
`v^p0` fino a `i1`, prosegue da lì al passo di `p1`, e da `i2` al passo di `p2`
— **ogni tratto con il proprio esponente e la propria velocità**, mai `p0·p1` o
`p0·p1·p2`. Con `p0` = 2, `p1` = 1 a 50 e `p2` = 2 a 100, tutti in jump e ogni
velocità a uno, il conteggio su `ident` è `v²` fino a 50; `v + 2450` da 50 a
100, un esponente di uno che prosegue a passo uno dai 2500 a cui era arrivato
il quadrato; e da 100 `2550 + 100·((v/100)² − 1)`, di nuovo al quadrato da un
passo di 2 lì. Una velocità di 3 sul primo gruppo rende il tratto di mezzo
`2500 + 3·(v − 50)`. L'ITER più basso va per primo qualunque gruppo lo
contenga, e **allo stesso ITER il secondo gruppo sostituisce il primo**,
raggio, velocità e tutto; un gruppo con ITER a zero o sotto mette il proprio
esponente e la propria velocità sul primo tratto.

**Un gruppo spento viene lasciato fuori**, insieme al suo ITER:

| primo gruppo | secondo gruppo | il conteggio corre a |
| --- | --- | --- |
| attivo | attivo | `p0` fino all'ITER più basso, poi l'esponente di quel gruppo fino all'altro ITER, poi quello dell'altro |
| spento | attivo | `p0` fino all'ITER del secondo, poi `p2` |
| attivo | spento | `p0` fino all'ITER del primo, poi `p1` |
| spento | spento | `p0` per ogni conteggio, la semplice potenza `v^p0` |

Entrambi partono spenti, e l'esponente iniziale e ogni velocità a uno, dove al
valore non viene fatto proprio nulla, quindi ogni funzione dà fino all'ultimo
bit ciò che dava prima che esistesse una cosa del genere.

I due modi di raccordare, per come cambiano il conteggio passato alla funzione
(scritti per `ident`, dove il conteggio è ciò che vede la tavolozza, con
`f(ITER)` il conteggio arrivando a ITER — `ITER^p0` al primo — `p` l'esponente
del gruppo, `s` la sua velocità e `r` il suo raggio):

| modalità | da ITER − r a ITER + r | oltre |
| --- | --- | --- |
| **jump** | — | `f(ITER) + s·ITER·((v/ITER)^p − 1)`: il passo svolta a ITER di colpo, il che si vede come un anello dove le bande cambiano passo |
| **smooth** | il passo va in linea retta da quello che il conteggio ha a `ITER − r` a quello che ha a `ITER + r`, a metà strada a ITER, con una parabola intorno a ITER aggiunta che restituisce ciò che la retta ha tagliato dall'angolo del jump | il conteggio del jump, fino all'ultima cifra |

**Smooth** è simmetrico rispetto a ITER e cambia i colori dentro il raggio e da
nessun'altra parte: prima di `ITER − r` e oltre `ITER + r` l'immagine è quella
del jump. Dove entrambi i passi sono costanti la parabola è nulla e l'angolo è
tagliato come lo taglia una parabola — con il conteggio intatto fino a un ITER
di 50 e un esponente di zero oltre, un raggio di 10 lo porta a 47.5 a ITER, un
quarto del raggio in meno, e di nuovo a 50 da 60 in poi; da un passo di uno a
uno di tre è due a ITER. Dove un passo si piega, la parabola restituisce ciò
che la piega aggiunge: dopo `v²` fino a 50 e un esponente di uno, un raggio di
20 porta il passo da 60 a 30 a 1 a 70, 45.5 a ITER, e prosegue da 70 come
`v + 2450`, come fa il jump.

Prima fondeva i passi dei due tratti come ciascuno sarebbe proseguito: il passo
del quadrato, ancora in crescita, continuava a crescere ben dentro
l'arrotondamento, quindi calava a malapena prima di ITER ed era a metà discesa
solo mezzo raggio dopo, e ciò che guadagnava spostava ogni colore oltre
l'arrotondamento — di 133 bande con un raggio di 20.

La parabola non porta mai il passo sotto il più lento dei due estremi. Dove lo
farebbe — un passo lento trasformato in uno molto veloce su un raggio corto,
come una centesima potenza su dieci iterazioni — viene tenuta lì, e il
conteggio oltre l'arrotondamento viene spostato di ciò che questo trattiene. Un
raggio più largo del primo ITER si ferma prima di zero, e due che si
sovrapporrebbero vengono ridotti allo stesso modo perché ognuno resti tra i
suoi vicini — raggi di 40 e 40 tra ITER a 30 di distanza diventano 15 ciascuno.

Un esponente di 100 è lecito e inutile oltre ITER: poche iterazioni dopo, le
bande sono più fini dei pixel e il colore è rumore. Con il secondo gruppo a
ITER 50 ed esponente 2 e il primo a 100 e 1.5, entrambi in jump, e un esponente
iniziale di uno, il conteggio su `ident` è 30 a 30, `50·1.5²` = 112.5 a 75, e
`200 + 100·(1.5^1.5 − 1)` = 284 a 150: il tratto da 100 parte dai 200 a cui era
arrivato quello prima, a un passo di 1.5.

Prima c'erano quattro modalità, e le posizioni salvate prima del raggio le
nominano — `outcolorfunafter` della 1.8.0, e `outcolorfunchange` con due gruppi
e senza raggio. Vengono lette nelle due che ci sono, con ogni velocità a uno:

| salvato come | letto come |
| --- | --- |
| off, jump | lo stesso |
| smooth, che arrotondava l'angolo nel logaritmo del conteggio, un po' da ciascun lato di ITER per quanto lontano | smooth con un raggio di tre quarti di ITER: entro quattro centesimi di ITER da ciò che disegnava fino a un esponente di uno e mezzo, un decimo fino a due |
| continuous, che proseguiva da ITER al passo con cui il conteggio arrivava | smooth a tre quarti di ITER, con l'esponente `0.66 + p/4 + 0.09·p²` — uno a uno, dove anche continuous non cambiava nulla — entro dodici centesimi di ITER da esso fino a un mezzo |

Trovato mettendo le curve una accanto all'altra su quattro ITER: un'immagine
del genere esce simile, non uguale al pixel. Una posizione il cui cambio di
velocità non è più di un gruppo in jump, con l'esponente iniziale e ogni
velocità a uno e il resto come parte, viene salvata come la salvava la 1.8.0,
`(outcolorfunafter 40 1.5 'jump)`, e si apre nella 1.8.0 e lì disegna uguale.
Qualunque altra cosa viene salvata come la salvava la 2.0, nell'ordine del
dialogo senza le velocità,
`(outcolorspeedchange 1 40 1.5 'smooth 10 50 2 'off 10)`; e le velocità, solo
dove una di esse non è uno, a parte come `(outcolorspeedlinear 1 3 1)` — la
velocità iniziale e quelle dei due gruppi — così che una posizione che non le
usa si apra ancora nella 2.0. Il solo cambio di velocità, come lo aveva il suo
dialogo, è il comando `outcolorspeedgroups`, con i suoi dodici campi nel loro
ordine.

### Oltre la fine di un int

La tavolozza viene raggiunta in 256esimi di voce attraverso un int, e un valore
oltre ciò che un int contiene — 8.4 milioni di bande — veniva convertito
comunque, il che è indefinito: ogni pixel del genere prendeva una sola e
medesima voce. `cube` ci arriva a 204 iterazioni, quindi un'immagine colorata
con il cubo oltre quel punto era una macchia piatta; `exp` e un forte cambio di
velocità ci arrivano subito. Ora finisce dove dice il suo valore, avvolto
intorno alla tavolozza come qualunque valore più piccolo; uno troppo grande
perché i suoi 256esimi stiano in un double finisce sulla prima voce, come ha
sempre fatto un infinito. Dentro un int nulla si sposta. Fuori, l'unico
cambiamento a un'immagine salvata prima è nella colorazione graduale di Magnet
e Magnet2 lontano dall'insieme, che era un unico colore piatto e ora mostra le
bande che il conteggio descrive lì.

Con essa sono arrivate due riparazioni più piccole. Scegliere una funzione di
velocità spostava il segno di spunta nel menu della colorazione interna invece
che nel proprio. E **initstate**, con cui comincia ogni posizione, non
rimetteva ai valori predefiniti i cinque numeri dell'fbm — una posizione li
salva solo dove differiscono da quelli, quindi una caricata dopo un'altra che
li aveva cambiati disegnava con quelli dell'altra; ora li reimposta, e con essi
spegne il cambio di velocità dopo ITER.

### Il conteggio, intero

Il valore che una modalità di colorazione consegna perdeva cifre sulla strada
verso la tavolozza in due punti. La colorazione graduale rendeva il suo
conteggio un int di 256esimi di iterazione prima che lo ricevessero la velocità
della colorazione, la funzione o il cambio di velocità; e ciò che la funzione
restituiva era un float, che oltre 65 536 voci conserva meno di un 256esimo di
una e oltre sedici milioni nemmeno una intera. Nessuno dei due si vede a
velocità uno, dove un 256esimo di iterazione è un 256esimo di voce. Entrambi si
vedono appena qualcosa fa correre i colori: con un cambio di velocità che porta
il conteggio da zero a diecimila voci per iterazione, una vista a trenta
iterazioni, dove il passo era 2500, usciva a strisce di un solo colore larghe
otto pixel, ciascuna un 256esimo di iterazione e dieci voci dalla successiva,
dove la tavolozza avrebbe dovuto scorrere ogni voce in mezzo.

Ora il conteggio prosegue intero, e la funzione, la velocità della colorazione
e il cambio di velocità sono calcolati nel tipo proprio della compilazione,
come le formule — long double nella compilazione normale, `__float128` in
quella quad — e restituiscono ciò che calcolano. Le due compilazioni colorano
simili quanto prima: sulle immagini provate differiscono dove hanno sempre
differito, in uno zoom profondo le cui orbite quad calcola più a fondo, e da
nessun'altra parte.

Due cose sono rimaste com'erano.

* **Un infinito** — il centro di un potenziale andato all'infinito, lo sfondo
  intorno a una formula che esce subito — finisce sulla prima voce, come un
  valore che non è un numero. Reso in 256esimi in un float, ogni valore del
  genere usciva infinito e finiva lì; quale colore debba prendere un infinito
  non lo sa nessuno, e così prende quello che ha sempre preso. Un valore finito
  oltre ciò che un float contiene finisce dove dice.
* **Un valore sul bordo di un 256esimo** ci resta: il valore viene arrotondato
  a un 256esimo di un 256esimo prima di essere tagliato a un 256esimo, dato che
  la colorazione non dice nulla di più fine. La parte reale di un'orbita di
  Newton si ferma su una radice, a uno o a meno un mezzo, il che mette
  `iter+real` esattamente su un bordo del genere, e le ultime cifre della
  convergenza metterebbero ogni pixel da una parte o dall'altra — un 256esimo
  di voce, che nessuno vede, finché un filtro dei bordi disegna ogni pixel del
  genere come bordo e le curve di livello escono dentellate come ingranaggi. Il
  float ne arrotondava via la maggior parte e lasciava tacche; ora non ce n'è
  nessuna.

Renderizzati rispetto alla compilazione precedente a 640×480 — prima del
cambiamento della sezione successiva — dei 169 esempi 83 disegnavano al pixel
come prima e altri 71 differivano di un solo livello in al massimo 113 pixel.
Degli altri quindici, quelli a colori reali e smooth log — `perb4`, `spiral4`,
`truecol1`, `truecol3` — differiscono fino a cinque livelli in fino a un
settimo dei pixel, dato che anche le loro componenti tornavano come float;
`edge` ha perso le sue tacche; e il resto differisce in una manciata di pixel
dove un valore stava sul bordo di una voce e un filtro o una tavolozza
disegnata senza sfumatura mostra il gradino. Magnet e Magnet2 in colorazione
graduale cambiano di più: lontano dall'insieme il loro conteggio è negativo,
cosa che l'int, letto come senza segno, spingeva avanti di quattro miliardi di
256esimi, e che il float poi disegnava a gradini di una voce o due — macchie
piatte con bordi netti dove ora c'è il gradiente che il conteggio descrive.

Su un'immagine che esce uguale costa dal due al cinque per cento nella
compilazione normale, e dal sette al venti in quella quad, dove i logaritmi del
conteggio e le funzioni sono calcolati in software. Dove le strisce diventano
un gradiente anche il rendering richiede più tempo: il riempimento a stima
riempie un pixel solo da vicini di un unico colore, e l'immagine scritta è più
grande — la vista a trenta iterazioni, a 1200×1200, impiegava 1.07 secondi dove
le sue strisce ne impiegavano 0.30.

### Il passo del conteggio graduale

La colorazione graduale aggiunge al conteggio delle iterazioni intere quanto
dentro l'ultima il punto è uscito. Lo calcolava dove il bailout cade tra la
grandezza prima della passata di uscita e dopo, linearmente nel suo logaritmo.
In quel modo il conteggio prosegue senza interruzioni da un'iterazione alla
successiva, ma il suo passo no: da una passata alla successiva la grandezza
viene quasi elevata al quadrato, quindi il passo cambia lungo ogni iterazione e
svolta di colpo dove comincia la successiva. In una vista larga 2.6×10⁻⁷, dove
il conteggio va da 59.94 a 60.02, svoltava di 4.26 sul bordo tra le due. A una
velocità di colorazione di uno ciò sposta un colore di qualche centesimo di
voce; con un cambio di velocità che corre a diecimila voci per iterazione, sono
le bande che cambiano densità di colpo lungo ogni bordo di un'iterazione.

Ora l'orbita viene portata oltre il bailout finché il suo modulo al quadrato
supera un milione — tre o quattro passate in più per un quadrato — e la
frazione viene presa nel logaritmo del logaritmo, dove un'elevazione al
quadrato è un passo di lunghezza uno:

    n + ln(ln T / ln s) / ln(ln m / ln s)

con `s` e `m` il modulo al quadrato da una parte e dall'altra di un milione,
`T` il bailout e `n` le passate fatte. Le passate oltre il bailout si tolgono
da sé, e il grado è quello dell'orbita stessa, `ln m / ln s`, quindi un cubo
viene contato come un cubo. Il conteggio resta entro un decimo di iterazione da
ciò che era, e il suo passo prosegue: la svolta di 4.26 in quella vista esce a
1.0000, e lungo un raggio da fuori dell'insieme le svolte da 1.3 a 1.8 a ogni
bordo escono a 1.0025.

| | il conteggio |
| --- | --- |
| le Mandelbrot, Triceratops, Mandelbar, Lambda, Manowar, Spider, Phoenix e Octo, sul cerchio | portato oltre il bailout |
| una formula utente che esce sul cerchio | portato oltre finché l'orbita cresce almeno come una potenza di uno e mezzo, e contato al bailout dove non lo fa |
| le tre Barnsley | al bailout: moltiplicano l'orbita per `c`, così che il suo logaritmo cresce di un passo invece di raddoppiare, e il conteggio al bailout è già quello giusto per loro |
| Newton, Magnet, Catseye | al bailout: convergono, o verificano l'uscita a modo loro |
| ogni forma della fuga tranne il cerchio | al bailout |

Le forme restano com'erano perché ciascuna misura l'orbita con qualcosa che
ruota con essa, mentre ogni passata ne raddoppia l'angolo: portato avanti con
la sua misura, il grado usciva abbastanza instabile da peggiorare le svolte —
quelle dell'esagono da 1.28 a tra 1.67 e 3.57 lungo quel raggio — e quelle del
quadrato persino tornavano indietro.

Dei 169 esempi, i 154 che non colorano in modo graduale disegnano al pixel come
prima di questo. I quindici che lo fanno cambiano fino a tredici livelli in
fino a nove decimi dei pixel — `eye` fino a 39, dove le voci vicine della sua
tavolozza sono molto distanti, e `spiral4`, il cui smooth log su una tavolozza
a colori reali ingrandisce qualunque cambiamento nel conteggio, nel mezzo di un
piccolo dettaglio. Rispetto alla 2.1, con i cambiamenti della sezione
precedente, 84 disegnano al pixel come prima e altri 59 differiscono al massimo
di tre livelli in al massimo 500 pixel.

Costa quelle passate e due logaritmi in più per ogni pixel che esce: sulle
immagini graduali provate, dall'11 al 28 per cento nella compilazione normale e
dal 18 al 20 in quella quad, dove un logaritmo in `__float128` richiede 946 ns
contro 39 in long double; un'immagine che non colora in modo graduale è
invariata. Calcolati in double, i quattro logaritmi sarebbero costati dall'8 al
19 per cento nella compilazione normale e avrebbero risparmiato l'8 in quella
quad — un double porta la frazione a sedici cifre, dove una velocità di
diecimila ne richiede nove — ma la colorazione è calcolata per intero nel tipo
proprio della compilazione, come le formule.

## Zoom

**Zoom a selezione** (Vista → Zoom) sostituisce lo zoom continuo con un
rettangolo: trascina con il tasto sinistro e rilascia per andarci, come fa
Ultra Fractal. Solo in ingrandimento. Il rettangolo viene forzato alle
proporzioni della finestra allungandone il lato corto, mai tagliando quello
lungo, e funziona correttamente con il frattale ruotato. La modalità Julia
veloce prende il trascinamento mentre è attiva, e lo restituisce all'uscita.

**Passi fissi**, nello stesso sottomenu:

| | |
| --- | --- |
| Ingrandisci 2× | `Ctrl` `+` |
| Riduci 2× | `Ctrl` `-` |
| Ingrandisci 10× | `Ctrl` `Alt` `+` |
| Riduci 10× | `Ctrl` `Alt` `-` |

Scalano rispetto al centro della vista, quindi si annullano esattamente a
vicenda.

## Ruotare e centrare la vista

**Frattale → Vista** imposta centro, raggio e angolo con dei numeri. Sotto
l'angolo ha due strumenti per quegli stessi campi. Riempiono i campi e non
fanno altro: OK applica ciò che mostrano come un solo passo di annullamento,
Annulla niente, e una posizione salvata non ne sa nulla.

**Fette** divide un giro completo in altrettanti giri uguali, e **−** e **+**
ruotano l'angolo indietro o avanti di uno: 60 gradi per sei, 360 per uno. Con
nessuna, i due sono disattivati. I giri si contano dall'angolo da cui sono
partiti, quindi sei sesti sono un giro intero fino all'ultima cifra; scrivere
un angolo, o cambiare le fette, fa ripartire il conteggio da lì. Il numero di
fette viene conservato finché XaoS è aperto, per la prossima volta che si apre
il dialogo.

**Centra sull'asse** sposta il centro dritto attraverso l'immagine (**R**) o
dritto verso l'alto (**I**) finché l'origine è in mezzo in quel senso,
all'angolo del dialogo. Senza rotazione, è la parte reale o immaginaria che
diventa zero e l'altra che resta. Con la rotazione, è lo stesso lungo gli assi
propri dell'immagine: R a 90 gradi azzera la parte immaginaria, e R e poi I
portano l'origine in mezzo a qualunque angolo. Un centro che è già lì viene
lasciato com'è, e una parte che la rotazione lascia come arrotondamento viene
presa per zero, così il campo dice 0 e non 1e-17.

La rotazione stessa è calcolata nel tipo proprio della compilazione, pi greco
compreso. Era calcolata da un pi greco double, il che nella compilazione quad
voleva dire anche un seno e un coseno double: una coppia che è una rotazione
solo fino a sedici cifre, quindi il centro di una vista ruotata veniva
disegnato a 5e-17 della sua distanza dall'origine lontano da dove si trovava.
Sotto 1e-16, dove serve quad, una vista ruotata di 30 gradi mostrava tutt'altro
posto. Un quarto di giro ora è esatto. Nella compilazione normale un'immagine
ruotata è la stessa salvo pixel al limite della sua precisione: 52 dei 64
esempi con un angolo sono invariati, gli altri spostano pixel sparsi dove il
dettaglio è più fine dell'aritmetica.

## La stessa immagine con qualunque numero di thread

Un'immagine usciva ogni tanto diversa di qualche pixel quando era disegnata con
più thread: `examples/Hubicka/julia3.xpf` renderizzato con quattro thread su
una macchina carica differiva in 26 render su 3000, fino a 37 pixel lungo i
suoi filamenti, e mai con un thread.

I pixel venivano dal riempimento a stima. Riempie un pixel dalle righe e dalle
colonne già calcolate intorno quando concordano, e in modalità Mandelbrot usa
gli stessi vicini per decidere se vale la pena controllare la periodicità di un
pixel. La stima procede per passate — prima una riga ogni sei, poi le righe in
mezzo — e le passate giravano in tutti i thread insieme senza nulla tra di
loro: un thread passava alla passata successiva mentre gli altri erano ancora
all'ultima, e una riga veniva segnata come fatta prima di essere calcolata.
Quindi se un vicino contasse come calcolato dipendeva da quanto avanti fossero
arrivati gli altri thread.

Ora ogni passata è in due fasi, prima le sue righe e poi le sue colonne, e ogni
thread finisce una fase prima che qualcuno cominci la successiva; ciò che una
fase calcola non viene usato come vicino al suo interno. I vicini sono allora
gli stessi qualunque sia il numero di thread, e i 3000 render escono tutti
uguali. Con un thread l'ordine è quello di sempre.

Anche i filtri escono uguali. **3D** disegna una fascia di colonne per thread,
e l'ombra di una colonna è portata dalla colonna due a destra nella riga sotto,
quindi le ultime due colonne di una fascia leggono le prime due della fascia
alla sua destra; prima ogni fascia leggeva lì copie della propria ultima
colonna, e l'immagine cambiava lungo i bordi delle fasce con il numero di
thread. Ora una fascia aspetta, riga per riga, ciò che la fascia alla sua
destra ha disegnato lì, e le fasce girano ancora fianco a fianco. **Campo
stellare** disegnava le sue stelle da un unico generatore casuale che ogni
thread faceva avanzare insieme; ora ogni riga parte dal punto della sequenza
che avrebbe raggiunto un solo thread che disegna l'intera immagine, quindi le
stelle sono le stesse con qualunque numero di thread. Cambiano ancora da un
avvio di XaoS al successivo, come hanno sempre fatto: il generatore prende il
seme dall'orologio.

**L'antialiasing** non è mai stato la causa, anche se tre esempi disegnati con
esso — `circle`, `helloween` e `warriormask` — uscivano diversi da un render
all'altro. Ciò che differiva era il messaggio con il nome della formula:
spariva dopo cinque secondi dell'orologio a muro, e un render rallentato da una
macchina carica lo disegnava in meno fotogrammi, o in nessuno. In un render un
messaggio ora dura cinque secondi del tempo proprio dell'animazione, qualunque
sia la macchina; nella finestra non cambia nulla.

Dei 169 esempi renderizzati con un thread, 167 disegnano al pixel come li
disegnava la 1.8.0: `galaxy` prende le sue stelle dall'orologio, e `heart` ora
conserva un messaggio che la 1.8.0, cronometrandolo con l'orologio a muro,
aveva tolto prima che l'immagine fosse disegnata. Renderizzati otto alla volta
su quattro core, una volta con un thread e tre volte con quattro, disegnano
tutti uguale tranne `galaxy`.

## Proporzioni

**Vista → Proporzioni** ridimensiona la finestra una volta, perché l'immagine
prenda la forma scelta: 1:1, 4:3, 3:2, 16:10, 16:9, 2:1, 21:9, o 9:16 in
verticale. La finestra non viene vincolata: dopo la si può trascinare a
qualunque forma, e il segno di spunta nel menu sta sulla forma che l'immagine
ha in quel momento, o su nessuna.

L'immagine conserva il numero di pixel che aveva e cambia solo la forma —
oppure, quando non starebbe sullo schermo, diventa la più grande di quella
forma che ci sta. Da schermo intero o da una finestra massimizzata torna prima
a una finestra normale, sempre con i pixel che aveva, quindi lì 16:9 è la
finestra 16:9 più larga che lo schermo contiene. La forma è rispettata al
pixel, anche su uno schermo in scala: al 125% un numero di pixel su cinque non
è ottenibile, e la dimensione viene scelta tra quelli che lo sono.

Ciò a cui serve è il render. Una posizione viene disegnata in modo che l'intera
vista che contiene stia nell'immagine, e un'immagine di un'altra forma mostra
più piano sul suo lato lungo. Quindi una vista inquadrata nella normale
finestra 4:3 e renderizzata a 1920×1080 ha un terzo in più ai lati di quanto ci
sia mai stato sullo schermo; inquadrata in una finestra 16:9, il render è ciò
che si è visto.

## Antialiasing in luce lineare

Filtri → Antialiasing in luce lineare, o `-linearaa`. L'antialiasing normale fa
la media dei colori dei campioni come sono memorizzati, il che scurisce i bordi
perché i valori memorizzati non sono proporzionali alla luce. Questo converte
in luce lineare, fa la media lì, e riconverte, così un bordo conserva la
luminosità di ciò tra cui sta.

## Lingue

Vista → Lingua passa XaoS a un'altra lingua, dopo un riavvio. Sei sono
complete: tedesco, spagnolo, francese, italiano, portoghese e russo — i menu, i
dialoghi e i messaggi, il riferimento delle formule utente, e i dialoghi OK,
Annulla, dei file e dei colori che vengono da Qt stesso. Questi ultimi ora sono
inclusi, mentre prima venivano cercati dove Qt è installato e nella lingua del
sistema qualunque fosse quella scelta. Le altre lingue che XaoS ha ereditato
arrivano fin dove le hanno portate i loro traduttori, e oltre sono in inglese.

Nel riferimento delle formule utente le descrizioni sono tradotte e i nomi e
gli elenchi degli argomenti no: sono la sintassi in cui si scrive una formula.

Questa guida è tradotta nelle stesse sei lingue, accanto alla guida inglese
`ms-xaos-guide.md` come `ms-xaos-guide.de.md`, `.es.md`, `.fr.md`, `.it.md`,
`.pt.md` e `.ru.md`, e Aiuto apre quella della lingua in cui è XaoS. Ognuna
dice da quale versione della guida inglese è stata fatta; dove differiscono, fa
fede quella inglese.

## Compilazione

    cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build

Entrambi gli eseguibili finiscono in `bin/`. Opzioni:

| | |
| --- | --- |
| `-DQUAD=OFF` | salta il secondo eseguibile |
| `-DDEEPZOOM=ON` | rende quad l'eseguibile principale, e non compila il secondo |
| `-DDIAG=ON` | compila anche `XaoS-diag`, che riporta in `xaos-diag.txt` l'ordine di sovrapposizione e di calcolo |

`cmake --build build --target deploy` riempie `bin/` con il runtime di Qt e del
compilatore, così che giri su una macchina senza Qt.

## Test

    ctest --test-dir build

351 test: il parser delle formule in entrambe le precisioni, l'accuratezza di
`gamma`, `lambertw`, `erf` e `randsc` rispetto a riferimenti esatti dove
esistono, i cicli di iterazione rispetto a checksum registrate, le funzioni di
velocità della colorazione rispetto a ciò che davano prima e a ciò a cui
servono, le posizioni salvate che vanno e tornano e ogni render che termina in
modo pulito, il salvataggio e ripristino della sovrapposizione, la cronologia
di annullamento, la rotazione e la centratura del dialogo Vista, il riferimento
dell'aiuto rispetto alla tabella del parser stesso, e il `poly` casuale della
sua scheda Utilità rispetto al polinomio scritto termine per termine.

Uno di essi misura la velocità invece delle risposte. Misurare in nanosecondi
asserirebbe su quale macchina gira la suite, quindi asserisce invece un
rapporto: quanto costa una degradazione alle funzioni di rumore, contato in
moltiplicazioni. Dovrebbe essere due, una per componente, perché la dimensione
viene moltiplicata per la degradazione una volta per passata; elevarla invece
alla potenza della passata dà sei o sette e cresce con il limite di iterazioni.
Il test fallisce sulle due implementazioni che questa ha sostituito.
