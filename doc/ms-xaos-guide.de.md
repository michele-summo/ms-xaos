# MS XaoS

> Übersetzung der [englischen Anleitung](ms-xaos-guide.md) der Version 2.7. Wo
> sie voneinander abweichen, gilt das englische Original.

Ein Fork von [XaoS](https://github.com/xaos-project/XaoS) 4.3.3. Alles, was das
Original kann, kann er weiterhin; diese Anleitung beschreibt, was hinzugekommen
oder geändert ist, und warum. Version 2.4.

## Zwei Programme

`XaoS.exe` rechnet mit 64 Bit Mantisse, `XaoS-quad.exe` mit 113. Beide werden
aus denselben Quellen gebaut — die Genauigkeit ist ein beim Kompilieren
gewählter Typ, keine Einstellung — und das Quad-Programm nennt sich **MS XaoS
Quad**, damit man die beiden im Betrieb unterscheiden kann.

Quad zoomt etwa fünfzehn Größenordnungen tiefer. Beim Berechnen derselben
Position bei 1e-18 behält das normale Programm 21% seiner Details, das
Quad-Programm 92%. Es ist nicht die Voreinstellung, weil `__float128` in
Software emuliert wird: etwa achtmal langsamer bei den eingebauten Formeln und
fünfzehnmal bei benutzerdefinierten.

Eine vom Quad-Programm gespeicherte Position enthält `(precision 113)`. Öffnet
man sie mit einem weniger genauen Build, wird das gemeldet und sie trotzdem
gezeichnet; eine Version, die den Befehl nicht kennt, verweigert das Öffnen,
und das ist Absicht — ein mit 113 Bit berechnetes Bild lässt sich mit 64 nicht
wiedergeben, und nicht zu öffnen ist ehrlicher, als etwas anderes zu zeichnen.

## Zahlen in den Dialogen

Ein Feld, das nach einer Zahl fragt, nimmt sie in jeder lesbaren Schreibweise,
1e-18 ebenso wie 0.5. Die Dialoge mit einem einzigen Feld benutzten ein
Drehfeld, das zwei Nachkommastellen behält, wenn man ihm nichts anderes sagt,
also überlebte 0.01 und 0.001 wurde zu nichts — eine Newton-Konvergenz von
einem Millionstel ließ sich überhaupt nicht eingeben.

Angezeigt wird der kürzeste Text, der sich als genau die gespeicherte Zahl
zurücklesen lässt, also sieht ein Millionstel eines Millionstels wie 1e-18 aus
und nicht wie 9.99999999999999999978e-19 — das macht das Ausgeben aller Stellen
des Builds aus einer Dezimalzahl, die eine binäre Gleitkommazahl nicht exakt
darstellen kann. Eine Koordinate, die alle einundzwanzig Stellen braucht,
bekommt sie weiterhin: sie zu brauchen heißt dasselbe, wie sich ohne sie nicht
unverändert zurücklesen zu lassen.

Die Exponentialform wird ab einer Billion aufwärts und unterhalb von hundert
Milliardsteln verwendet, dazwischen die einfache Dezimalform — also ist zwanzig
zwanzig, hundert Milliarden werden ausgeschrieben, und 1e-8 nicht. Sich selbst
überlassen, würde die Ausgabe wählen, indem sie den Exponenten mit der
verlangten Stellenzahl vergleicht, was nichts nützt, wenn die Stellen nach
Kürze gewählt werden: eine Stelle genügt, um sich als zwanzig zurückzulesen,
und zwanzig kam als 2e+01 heraus.

## Gespeicherte Positionen

Der Ausschnitt wird mit der vollen Genauigkeit des Builds geschrieben —
normalerweise 21 signifikante Stellen, 36 bei Quad. Früher wurde er mit einer
aus der Zoomtiefe abgeleiteten Stellenzahl geschrieben, die zu niedrig ausfiel
und bei 17 endete, sodass sich eine tiefe Position einen Bruchteil eines Pixels
neben der Stelle öffnete, an der sie gespeichert worden war. Bei hoher
Iterationszahl ist das ein sichtbar anderes Bild.

Alles, was eine Position enthalten kann, übersteht nun Speichern und Laden
unverändert: alle 169 mitgelieferten Beispiele und 141 erzeugte Positionen, die
jede Formel, Ebene, jeden Färbungsmodus und Filter abdecken, kommen pixelgenau
zurück.

## Rückgängig

Bearbeiten → **Rückgängig** (`U` oder `Strg` `Z`) geht einen Schritt zurück und
**Wiederholen** (`Strg` `Y` unter Windows) wieder vor. Ein Schritt ist alles,
was das Bild ändert: ein Befehl — eine Formel, eine Palette, ein Färbungsmodus,
ein Dialog, ein fester Zoom —, ein Mausdruck mit dem Zoom, Verschieben oder
Drehen, das er in Gang setzt, ein Auswahlzoom und die mit den Pfeiltasten
gesetzten Iterationen. Ein Befehl, der nichts ändert, wie Neu berechnen, ist
kein Schritt, und Rückgängig geht nicht weiter zurück als zu dem Bild, mit dem
XaoS nach Konfiguration und Kommandozeile geöffnet wurde. Jeder Schritt wird
als gespeicherte Position aufbewahrt, also stellt Rückgängig alles wieder her,
was eine Position enthält, so genau wie ihr erneutes Öffnen; die letzten 256
werden aufbewahrt.

Es ging auf zwei Arten schief. Die Maus bewahrte ihren Schritt beim Loslassen
auf, während der Zoom noch abbremste, also brachte Rückgängig ein Bild auf
halbem Weg zurück; jetzt wird der Schritt dort aufbewahrt, wo der Ausschnitt
zur Ruhe kommt. Und jeder Befehl bewahrte einen Schritt auf, ob er etwas
änderte oder nicht, während Rückgängig zwei zurückging, in der Annahme, den
letzten selbst aufbewahrt zu haben — nach einem Neu berechnen brachte es das
Bild zurück, das schon auf dem Bildschirm war. Jetzt wird ein Bild nur
aufbewahrt, wo es sich vom zuletzt aufbewahrten unterscheidet, und Rückgängig
und Wiederholen bewegen sich einen Schritt entlang dieser Bilder. Ein neuer
Schritt nach einem Rückgängig schneidet ab, was sich hätte wiederholen lassen,
wie in einem Editor.

Was sich seit dem letzten Schritt von selbst bewegt hat — der Autopilot, eine
fortlaufende Rotation —, wird beim Drücken von Rückgängig als eigener Schritt
aufbewahrt, damit Wiederholen dorthin zurückkommt.

Rückgängig fand einen Fehler in den Positionen. Lud man eine mit eingebauter
Palette nach einer Benutzerpalette, blieb die Benutzerpalette als verwendet
markiert, also wurde die Position wieder als Benutzerpalette gespeichert: nicht
das, was geladen worden war, und für Rückgängig nicht der Schritt, zu dem es
zurückgegangen war, was das Wiederholen kostete. Jetzt löscht eine eingebaute
Palette die Markierung, wie es eine erzeugte schon immer tat.

## Benutzerformeln

**Referenz.** Hilfe → Referenz der Benutzerformeln listet jede Funktion,
Variable und Schreibweise auf, die die Formelsprache annimmt, in einem Fenster,
das man neben der gerade geschriebenen Formel offen lassen kann. Sie ist danach
gruppiert, wozu die Funktionen dienen, mit den speziellen Funktionen — `erf`,
`gamma`, `lambertw` — in einem eigenen Abschnitt, da keine davon elementar ist
und man nicht zufällig zu ihnen greift. Ein vierter Reiter listet die Zahlen
auf, die als Argumente vorkommen und etwas Bestimmtes bedeuten — an welcher
Form eine Falle misst, mit welchem Spiegel das Kaleidoskop faltet —, statt sie
in der Beschreibung der Funktion zu vergraben, die sie nimmt. Ein fünfter zeigt
ein Bild jeder Kachelung, die `randsctile` nimmt, denn fünfundvierzig Zahlen
mit je einer Zeile Worte sagen wenig über die Formen.

Ein sechster, **Werkzeuge**, schreibt Formeln, statt sie zu beschreiben.
Vorerst enthält er ein zufälliges `poly` aus vier Feldern:

| | |
| --- | --- |
| minimale und maximale Länge | m, die Anzahl der Koeffizienten, zwischen beiden gezogen; der Grad ist m − 1. Beide werden gebraucht, als ganze Zahlen von 1 bis 254 — so weit `poly` reicht, da der Parser 255 Argumente pro Aufruf nimmt und `z` eines davon ist |
| minimaler und maximaler Betrag | die Größe jedes Koeffizienten, gleichmäßig zwischen beiden gezogen. 0 und 1, wenn leer gelassen; höchstens 1000000 |

Jeder Koeffizient ist komplex, `{re,im}`, mit einem über die ganze Umdrehung
zufällig gezogenen Argument, und jeder Teil wird mit vier signifikanten Stellen
seines Betrags geschrieben, damit ein Bereich von 0.001 bis 0.002 nicht
weggerundet wird. In der Formel steht kein `c`. Eine Benutzerformel lässt `z`
im Mandelbrot- wie im Juliamodus beim Pixel beginnen, also zeichnet sie die
gefüllte Julia-Menge des Polynoms selbst, in beiden gleich; `c` kommt dahin, wo
man es haben will, wie in `poly(...) + c`.

Jede Änderung an einem Feld erzeugt eine neue Formel, ebenso **Ein weiteres
erzeugen**, und **Kopieren** legt sie in die Zwischenablage für Fraktal →
Benutzerformel. Solange eine Länge leer ist oder ein Feld keine Zahl ist,
außerhalb seines Bereichs liegt oder sein Minimum über dem Maximum hat, gibt es
keine Formel, nur eine Zeile, die sagt, warum.

Jede Zeile der Funktionsliste sagt, was der Aufruf nimmt und in welcher
Reihenfolge, mit Namen, die sagen, wofür jede Stelle da ist — "base, exponent",
"z, low, high", "before, after, n" — statt mit Buchstaben, die nichts darüber
sagten, welche Stelle welche ist; dieselben Namen zeigt die Zeile unter der
Formelleiste. Eine Konvention zieht sich durch alle: `z` ist der Wert, an dem
eine Funktion arbeitet, was auch immer tatsächlich dort steht. Die Liste und
diese Anzahlen werden beide von einem Test gegen die eigene Tabelle des Parsers
geprüft, also kann eine Funktion, die ein Argument hinzugewinnt oder verliert,
nicht weiter mit den alten beschrieben werden. Ein Argument, das weggelassen
werden darf, steht in Klammern mit dem Wert, den es dann annimmt:
`z, [power=1], [angle_factor=1]`.

**Argumente werden durch ein Komma getrennt, und durch nichts anderes.** Früher
ging auch ein Semikolon: der Parser machte es zu einem Komma, bevor er
irgendetwas las, also waren `LOGN(5;Z^2)+C` und `LOGN(5,Z^2)+C` eine Formel,
und dasselbe galt für die beiden Teile einer komplexen Zahl. Das ist vorbei.
Ein Semikolon wird jetzt dort abgelehnt, wo es steht — *Falsche Operation: ;* —
und das Komma ist überall das Trennzeichen, wie in der ganzen übrigen Syntax.

### Die Palettensonde

**Ausschnitt → Palettensonde (auf einen Punkt klicken)** macht den Mauszeiger
zu einer Sonde: klicken Sie irgendwo auf das Fraktal, und die Nachrichtenzeile
sagt, **von welcher Stelle der Palette dieser Punkt gezeichnet wurde**, als
Zahl statt als Farbe — *Palettenstelle 12.375*, und sonst nichts.

Der Paletteneditor legt **31 Farben** fest, und die Palette entsteht, indem man
von jeder zur nächsten **acht Zellen auf einmal** geht. Die Zahl ist also die
Zelle geteilt durch acht: `12.000` ist genau die zwölfte Farbe des Editors,
`12.375` liegt drei Achtel des Wegs von der zwölften zur dreizehnten. Eine
Zelle ist ein Achtel einer Farbe, feiner, als irgendjemand am Editor abliest.

* Mit der **Benutzerpalette** wird die Zahl auf **`0` bis `31`** gebracht, die
  Nummerierung des Editors selbst.
* Mit jeder anderen Palette wird sie **gelassen, wo sie hinfällt, und reicht
  weit über 31 hinaus** — eine Palette mit 65534 Zellen reicht bis 8191. Die
  Zahl außerhalb des Bereichs sagt, dass die Palette nicht die ist, die der
  Editor zeichnet.
* Ein Punkt, der die **innere Farbe** annahm, ergibt *Palettenstelle: keine*:
  diese Farbe ist ein einziger flacher Ton, für den keine Zelle der Palette
  steht.

Wozu sie dient, ist weniger „welche Farbe ist das“ als **welchen Teil der
Palette das Bild tatsächlich benutzt**, und die Antwort ist oft ernüchternd.
Mit der einfachen `iter`-Färbung bei 60 Iterationen laufen die Zellen von 1 bis
59 — `0.125` bis `7.375` —, also berührt das Bild die Farben 8 bis 30 des
Editors nie, und sie zu bearbeiten ändert nichts. Eine Färbung, die sich
ausbreitet, etwa `iter + real`, durchläuft dagegen den ganzen Ring, mehrmals.

Sie liest den Punkt und nicht den Bildschirm, ist also genau: dieselbe Farbe
kann in zwei Zellen stehen, und ein Echtfarbenbild mischt zwei benachbarte
Zellen, also benennt keine Farbe eine eigene Zelle. Wie der Auswahlzoom gibt
die Sonde den Mauszeiger an den schnellen Juliamodus zurück, solange dieser
läuft.

**Die Leiste.** Fraktal → Benutzerformel öffnet einen Dialog, der etwa halb so
breit ist wie ein 1280×720-Bildschirm — rund sechshundertfünfzig Pixel, weniger
auf einem kleineren Bildschirm —, dessen Leiste in eine neue Zeile umbricht,
wenn die Formel nicht in eine passt. Sie bricht **zwischen Namen um, nie
innerhalb eines Namens**: ein Name oder eine Zahl (`randsctile`, `z_p3`, `0.6`)
bleibt ganz, und eine Zeile darf überall sonst enden — nach einer Klammer,
einem Komma, einem Operator. Der eigene Zeilenumbruch von Qt taugte dafür
nicht: er bricht dort um, wo die Regeln für Fließtext es erlauben, was in einer
Formel nach ein paar Kommas und wenig sonst ist, und schneidet dann einen
Namen, wo es sein muss. Nur ein Name, der breiter ist als die ganze Leiste,
wird je geschnitten. Die Leiste wächst mit der Formel um je eine Zeile, und
sobald der Dialog nicht mehr auf einen 720p-Bildschirm passen würde, hört sie
auf zu wachsen und scrollt. Für alles andere bleibt es eine Formel in einer
Zeile: Eingabe ist OK und keine neue Zeile, und in einer eingefügten Formel
werden Zeilenumbrüche und Tabulatoren zu Leerzeichen, die der Parser als nichts
liest.

**Verlauf...** unter der Leiste öffnet die früher verwendeten Formeln — die
zehn letzten, die neueste zuerst, jede einmal —, die die alte einzeilige Leiste
in ihrer Aufklappliste hielt. Eine für das Fenster zu lange Formel wird mit
ihrem Anfang, drei Punkten und ihrer Länge aufgeführt,
`poly(z, {0.01097,0.00699}, {-0.05... (157 Zeichen)`, damit sich zwei, die
gleich beginnen, noch unterscheiden lassen; der Anfang endet zwischen Namen,
wie eine Zeile der Leiste, und es ist umso mehr davon zu sehen, je breiter das
Fenster ist. Doppelklicken Sie auf eine, oder wählen Sie sie aus und drücken
Sie Verwenden, und sie kommt ganz in die Leiste. Eine Formel wird hinzugefügt,
wenn der Dialog mit ihr bestätigt wird, und nur dann: eine Position mit
Benutzerformel zu laden, eine abzuspielen oder eine mit `-render` zu berechnen
führt denselben Befehl `usrform` aus, schiebt sie aber nicht mehr in den
Verlauf, wo jede davon früher den Platz einer selbst getippten einnahm.

**Pfeil nach oben und Pfeil nach unten** gehen dieselben Formeln durch, ohne
die Leiste zu verlassen. Pfeil nach oben bringt aus der ersten Zeile der Leiste
die Formel, die vor der gezeigten verwendet wurde; Pfeil nach unten geht aus
der letzten Zeile zurück zur neuesten und dann zu dem, was in der Leiste stand,
bevor Sie begannen. Die zwei kleinen Pfeile neben der Leiste tun dasselbe,
jeder ausgegraut, wenn es in seiner Richtung nicht weitergeht. Eine Formel, die
schon in der Leiste steht, wird übersprungen, und was Sie tippen, wird zum
neuen Ausgangspunkt.

Die zuletzt verwendete Formel steht blass in der Leiste: ganz, solange die
Leiste leer ist, und ihr Rest hinter dem Getippten, solange das Getippte ihr
Anfang ist. Tab setzt sie ganz in die Leiste. Sobald das Getippte anders
weitergeht, verschwindet der blasse Text, und Tab geht wie bisher zu den
Schaltflächen weiter.

**Ausschnitt → Status** zeigt die Formel genauso: in einer Zeile, wenn sie ins
Fenster passt, sonst in zwei — die zweite unter der ersten, zwischen Namen
umgebrochen —, und wenn auch zwei sie nicht fassen, endet die zweite mit drei
Punkten und ihrer Länge, `... (918 Zeichen)`. Die Zeilen darunter rücken nach
unten, um Platz zu machen, und alles wird neu angeordnet, wenn das Fenster
seine Breite ändert. Früher war es eine einzige Zeile, die über den Bildrand
hinauslief.

**Unter der Leiste steht eine Zeile, die sagt, was der Aufruf nimmt, in dem man
sich befindet**, mit dem gerade geschriebenen Argument fett. Setzen Sie den
Cursor in `randsc(`, und sie lautet *randsc(**seed**, [size=1+i],
[degradation=0.5+0.5i], [kaleidoscope=1], [kaled_mirror=0], [skew=0])*; gehen
Sie über ein Komma, und das Fette wandert mit. Sie braucht von der Formel nicht
mehr als eine offene Klammer mit einem Namen davor, ist also da, während die
Formel halb geschrieben ist, und genau dann wird sie gebraucht. Setzen Sie den
Cursor außerhalb jedes Aufrufs, und die Zeile ist leer.

Sechs der mitgelieferten Positionen waren auf die alte Weise geschrieben und
wurden mit Kommas neu geschrieben: `circle`, `heart`, `helloween`, `pentafrac`,
`warriormask` und `burnship`, alle unter **Datei → Öffnen → examples/Malczak**.
**Eine eigene Position mit einem Semikolon in der Formel lädt nicht**, bis die
Formel mit Kommas geschrieben ist — die Meldung nennt das Zeichen, und die
Abhilfe ist, die Formel zu öffnen und sie zu ersetzen.

**Eine Stelle leer lassen.** Ein in Klammern gezeigtes Argument darf mitten in
einem Aufruf ebenso wie am Ende weggelassen werden, indem man zwischen den
beiden Trennzeichen nichts schreibt: `julian(z, ,3)` gibt das erste und das
dritte und lässt die Funktion sagen, was das zweite ist. Leerzeichen spielen
keine Rolle — `f(z, ,5)` und `f(z,,5)` sind derselbe Aufruf —, und ein
Argument, das nicht in Klammern steht, muss trotzdem geschrieben werden, also
werden `poly( ,1,2)` und `sin(,)` abgelehnt statt erraten.

Was eine leere Stelle bedeutet, ist Sache der Funktion. Die meisten nehmen den
Standardwert, den sie angeben; ein leer gelassener Koeffizient von `poly` ist
ein Term, der nicht da ist, also ist `poly(z,1, ,1)` gleich `z^2 + 1`; und ein
leer gelassener Zweig von `ifiter` oder `ifiterl` wiederholt den davor, so
bekommt ein Zweig mehr als einen Durchlauf im Zyklus.
`ifiter(f(z), , ,g(z), , , , )` führt `f` drei Durchläufe lang und `g` fünf
lang aus und kostet bei der Auswertung nicht mehr, als sie auszuschreiben:
wohin jede Wahl führt, wird einmal festgelegt, wenn die Formel geparst wird.

Eine noch nicht geschriebene Benutzerformel lautet `z^2+c`. Früher lautete sie
Burning Ship, ein eigenes Fraktal und ein Rätsel als Ausgangspunkt.

**`ifiterf(rest, last)`** — wertet `rest` in jedem Durchlauf außer dem letzten
aus und `last` in diesem. Der letzte Durchlauf ist der letzte, den das
Iterationslimit erlaubt: eine Formel kann nicht wissen, welcher Durchlauf der
sein wird, in dem sie austritt, da das vom Wert abhängt, den sie noch nicht
erzeugt hat. Nur der gewählte wird ausgewertet.

**`ifiterr(before, after, n)`** — wertet `before` aus, solange die
Durchlaufnummer unter `n` liegt, und `after` ab `n`. Nur der gewählte wird
ausgewertet, wie bei den anderen beiden, obwohl dem Parser die Schwelle erst
beigebracht werden musste: ein Argument kann jetzt als vom Selektor gelesen
statt von ihm gewählt markiert werden und wird dann ausgewertet, bevor die Wahl
getroffen wird.

**Die Variablen.** `z` ist der laufende Wert; `c` ist der Punkt, also das Pixel
im Mandelbrotmodus und die Konstante im Juliamodus; `x` ist die einfache
Koordinate, also das Pixel in beiden Modi, und ändert sich von Durchlauf zu
Durchlauf nicht; `n` ist die Iterationsnummer. In der *Benutzerinitialisierung*
ist `z` dort, wo `z` ohne Initialisierung begonnen hätte, was derselbe Wert ist
wie `x`, also sagt eine Initialisierung nur aus `z` genau das, was nichts zu
sagen sagt.

Dieser letzte Teil ist neu, und ebenso, dass man es sagen kann: `z` war bei
keinem der beiden Parser registriert, durch die die Initialisierung läuft, und
`x` beim falschen, also wurde eine Initialisierung, die eines davon nannte, von
dem Thread abgelehnt, der sie berechnen musste -- während der Dialog, dessen
Parser beide kennt, sie annahm und kein Zeichen gab, dass sie danach verworfen
wurde. `n` behielt, was das vorige Pixel dort gelassen hatte, das heißt, eine
Initialisierung, die es las, las die Reihenfolge, in der die Pixel zufällig
berechnet wurden.

**`p1`, `p2` ... `p9999`** — der Wert, den `z` in einem früheren Durchlauf
hatte: `p1` im Durchlauf vor diesem, `p2` in dem davor, `p9999`
neuntausendneunhundertneunundneunzig Durchläufe zurück. `p` ist ein anderer
Name für `p1`. Bevor es so viele Durchläufe gab, stehen sie auf dem, womit der
Verlauf beginnt, und das entscheidet *Fraktal → p-Werte bei der ersten
Iteration setzen*: der iterierte Punkt, wenn es eingeschaltet ist, null, wenn
es ausgeschaltet ist.

Früher hörten sie bei `p6` auf und wurden aufbewahrt, indem ein sechsstelliges
Array in jedem Durchlauf um eins weitergeschoben wurde — eine Kopie pro Stelle
pro Durchlauf, bezahlt, ob die Formel eine davon nannte oder nicht. Deshalb
hörten sie bei sechs auf: bei neuntausend wären es 320 KB Kopieren pro
Durchlauf gewesen, gemessen das Vierundzwanzigfache der gesamten Kosten,
`z^2+c` zu iterieren. Der Verlauf ist jetzt ein Ring, einmal pro Durchlauf
geschrieben, und nur die Stellen, die eine Formel tatsächlich nennt, werden
daraus gelesen, also kostet `p9999`, was `p1` kostet, und eine Formel, die
keine nennt, kostet gar nichts. Das Referenzfenster führt sie unter den
Variablen auf.

**Suffixe an einer Variablen.** `z`, `c`, `x`, `p` und `p1` bis `p9999` nehmen
Suffixe, die für die am häufigsten darauf angewandten Aufrufe stehen:

| Suffix | steht für |
| --- | --- |
| `_b` | `bship(…)` |
| `_bi` | `bshipi(…)` |
| `_br` | `bshipr(…)` |
| `_pM` | `parchment(…, M)` |
| `_paM` | `parchmenta(…, M)` |

Sie werden von links nach rechts gelesen, jedes umschließt, was die vorigen
ergeben haben: `z_p3` ist `parchment(z,3)`, `c_b_p2` ist
`parchment(bship(c),2)` und `p12_p2_p3` ist `parchment(parchment(p12,2),3)`.
*M* ist Teil des Namens, wird also in Ziffern geschrieben — eine ganze Zahl,
eins oder mehr; `z_p0` wird mit einer Meldung abgelehnt, die das sagt. Nicht an
`n`, das eine Anzahl ist und keine zweite Komponente hat, an der eines davon
arbeiten könnte.

Sie werden zu den Aufrufen ausgeschrieben, bevor die Formel gelesen wird, also
zeichnet eine Variable mit Suffix genau, was die ausgeschriebenen Aufrufe
zeichnen, und kostet, was sie kosten; die Formel wird so gespeichert, wie sie
geschrieben wurde. Alles, was keines der fünf ist, bleibt, wie es ist, und wird
als der unbekannte Name abgelehnt, der es ist.

**Entfernt.** `powi`, `powdc` und `logcn` waren zweite Namen für `pow` und
`logn`; `rad`, `deg` und `sign` waren ohne jede Implementierung dahinter
aufgeführt, und ein zweites `trunc` wurde vom funktionierenden verdeckt. Eine
Formel, die einen entfernten Alias benutzt, sollte den verbliebenen Namen
verwenden — die Funktion ist identisch.

**`erf(z)`** — die Fehlerfunktion in der komplexen Ebene. Genau auf etwa drei
ulp innerhalb eines Fluchtradius von zwei, wo ein Fraktal tatsächlich iteriert.

**`poly(z, k1, k2, ..., km)`** — ein Polynom in `z`:

    k1*z^(m-1) + k2*z^(m-2) + ... + k(m-1)*z + km

Der erste geschriebene Koeffizient multipliziert die höchste Potenz, und der
letzte steht allein, also liest sich der Aufruf in der Reihenfolge, in der man
das Polynom ausspricht. Berechnet nach dem Horner-Schema, also m−1
Multiplikationen statt der m(m−1)/2, die das getrennte Potenzieren jeder Potenz
bräuchte, und obendrein das genauere der beiden.

    poly(z,1,0,0)+c            die Mandelbrot-Menge, ausgeschrieben
    poly(z,1,0,0,{0.7,0.2})    z^3 + 0.7+0.2i

**`randsc(seed, size, degradation, kaleidoscope, kaled_mirror, skew, skew_mode, selfsim)`**
— kohärentes Rauschen über dem Punkt, das Flecken statt Schnee pro Pixel
ergibt. `size` (Standard `1+i`) ist die mittlere Breite eines Flecks entlang
der reellen Achse und seine Höhe entlang der imaginären. Es ist
**Gradientenrauschen**, eine Oktave der Bewegung, die die
[fbm-Modi](#fraktionale-brownsche-bewegung) und `fbm()` summieren. Es war
Wertrauschen und zeichnete aus demselben Grund wie ihres sein Gitter: die
Flecken wurden entlang jeder seiner Linien flach und kamen eckig heraus. Die
Flecken sind so breit wie vorher — das Gitter liegt bei drei Fünfteln der
Größe, und das Feld hört bei elf Sechzehnteln der Größe auf, sich zur Hälfte
ähnlich zu sein, vorher wie nachher —, und die Werte streuen wie vorher, 0.2145
von der Mitte, wobei die Kontrastkurve des fbm die engere Streuung des
Gradientenrauschens auf die des Wertrauschens zurückbringt. Aber jeder Wert ist
ein anderer, und ein Aufruf kostet 167 ns gegen 147. Die fünf Mosaike sind ihr
Gitter und sind bis aufs Bit unverändert. `degradation` (Standard `0.5+0.5i`,
halbiert sie in jedem Durchlauf) verkleinert sie im Lauf der Iteration: die
Größe wird in jedem Durchlauf komponentenweise damit multipliziert, also ergibt
`0.5+0.2i` über `1+i` im ersten Durchlauf `1+i`, dann `0.5+0.2i`, dann
`0.25+0.04i`. Eine Null in einer Komponente eines der beiden Argumente ergibt
null, statt durch null zu teilen. Dann kommen die beiden Kaleidoskop-Argumente,
`skew` (Standard `0`), `skew_mode` und `selfsim`, jedes unten beschrieben. Nur
der Anfangswert ist Pflicht.

### Die Schräge, und warum diese Felder schwer zu färben sind

Die Engine färbt mit den beiden Komponenten des Orbits: `real` und `zmag` lesen
die eine, `imag`, `angle` und `real / imag` die andere. Diese Felder geben für
eine ganze Zelle eine einzige Zahl zurück, auf der reellen Achse, also zeichnet
jeder Modus **einen Ton pro Zelle** — `zmag` und `iter + real` sehen aus wie
der abgeschnittene Wert, und die anderen drei haben überhaupt nichts zu lesen.
Gemessen über neunhundert Pixel: 14 bis 126 Werte für die ersten beiden und
genau **einer** für den Rest.

Den Wert innerhalb einer Zelle zu verändern, geht mit einer einzigen reellen
Zahl nicht. Dieselbe Zahl entscheidet, ob der Punkt austritt, also lässt ein
Wert, der sich über eine Zelle bewegt, die halbe Zelle austreten und die andere
Hälfte drin, und die Zelle kommt **von einer geraden Linie in zwei geteilt**
heraus. Das wurde zweimal versucht und hat beide Male die Mosaike zerschnitten.

`skew` tut es stattdessen mit der zweiten Komponente, und es **dreht** den
Wert, statt ihn zu skalieren. Wie weit es dreht, folgt dem, **wie weit der
Punkt außerhalb der Mitte seiner Zelle steht**: mit `out`, das von null in der
Mitte bis eins am Rand läuft, wird der Wert um das Doppelte des Arkustangens
von `skew_re*out + skew_im` gedreht.

`out` wird **in der eigenen Geometrie der Zelle** gemessen, also ist das, was
die Farbe in einer Zelle zeichnet, die Form, in die das Feld geteilt ist, und
kein darübergelegter Satz Linien — null in der Mitte einer Zelle, eins entlang
ihres ganzen Randes:

| Feld | was die Farbe in einer Zelle zeichnet |
| --- | --- |
| `randsc` | seine eigenen Flecken — die Drehung folgt dem Pegel, also sind die Konturen die des Feldes |
| `randscq` | Quadrate, um die Mitte des Quadrats |
| `randsch` | Sechsecke, um die Mitte des Sechsecks |
| `randsct` | Dreiecke, um die Mitte des Dreiecks |
| `randscp` | sein eigenes Vieleck, Schritt für Schritt verkleinert — keine zwei Zellen haben dieselbe Form |

Der **Realteil** der Schräge ist dieser Gradient über eine Zelle. Der
**Imaginärteil** ist eine gleichmäßige Drehung, die die ganze Zelle teilt: er
verschiebt die Farbe einer Zelle, ohne darin etwas zu zeichnen, also ergibt
`0.4+0.2i` die Konturen der Form und eine Verschiebung pro Zelle zugleich.

### Was die Schräge zeichnet: `skew_mode`

Das siebte Argument sagt, **was sich innerhalb der Zelle ändert**. Die Zahlen
**addieren sich**, statt einander auszuschließen, also kann ein Aufruf zwei
davon verlangen und beide bekommen.

| | |
| --- | --- |
| **1** (und `0`) | **die Form** — der eigene Umriss der Zelle, Schritt für Schritt verkleinert. Was ein Aufruf ohne `skew_mode` bekommt, und was verwendet wird, wenn weder 1 noch 2 genannt ist |
| **2** | **die Rosette** — der Winkel um die Mitte der Zelle, in so viele Umdrehungen gefaltet, wie das Kaleidoskop Keile hat, sodass jede Zelle die Symmetrie des Bildes selbst trägt. Einfache Speichen, wenn nichts gefaltet ist |
| **4** | **quer über den Keil** — eine Drehung, die dem folgt, wo der Punkt in seinem Keil des Kaleidoskops steht, in jedem Keil gleich |
| **8** | **radial** — der Betrag bewegt sich ebenso wie der Winkel |

`3` ist Form und Rosette zugleich, was eine Spirale ergibt. `11` ist diese
beiden mit sich ebenfalls bewegendem Betrag.

**`4` folgt der eigenen Geometrie des Kaleidoskops.** Es dreht den Wert an den
Rändern jedes Keils um nichts und entlang der Linie durch seine Mitte um die
ganze Drehung, die der Imaginärteil allein ergibt — das Doppelte seines
Arkustangens, 127 Grad für `2` —, und dazwischen um die Drehung des
Imaginärteils mal den durchquerten Anteil. Jeder Keil bekommt dasselbe, also
bleibt ein Kaleidoskop aus sechs ein Kaleidoskop aus sechs, samt Spiegeln, und
die Drehung ist über die Ränder und die Mitten hinweg stetig, wo die Faltung
spiegelt. Es tut nichts, wenn nichts gefaltet ist oder der Imaginärteil der
Schräge null ist. Keine Trigonometrie; es kostet bei 64 Bit Mantisse nichts und
bei 113 ein Zehntel eines Aufrufs.

Zuerst war es **eine andere Drehung in jedem Keil** — der Imaginärteil mal die
Nummer des Keils plus eins —, gedacht, damit die Kopien aufhören, gleich zu
sein. Das ist das Einzige, was die Kopien eines Kaleidoskops nicht sein können:
ein sechszähliges Bild kam als sechs unzusammenhängende Stücke mit einem
geraden Schnitt an jeder Naht heraus. Das Vielfache stieß außerdem an die halbe
Umdrehung, also lagen vier von sechs Stücken fast auf der reellen Achse, und
dort zeichnete `real / imag` sie als Schnee, wie auch immer `selfsim` das
Rauschen mittelte.

Was davon bleibt, ist eine Eigenschaft von `real / imag`, das durch den
Imaginärteil von *z* teilt: wo auch immer eine Drehung ihn über null trägt,
springt die Färbung von einem Ende der Palette zum anderen, und die letzten
Reste der feinen Durchläufe entscheiden die Seite Pixel für Pixel. Eine
Drehung, die über den Keil streicht, kreuzt dort bei den meisten Einstellungen.
Bei einem sechszähligen `randsct`, gefärbt mit `real / imag` und `selfsim` bei
eins, ließ ein Imaginärteil von `0.25` gar kein Korn zurück und `2` Bänder
davon neben den Rändern; eine Färbung, die nicht durch eine Komponente teilt,
hat keines.

**`8` ist der, der `zmag` antwortet.** Eine Drehung lässt den Betrag, wo er
ist, und `zmag` und der Fluchtradius lesen den Betrag und sonst nichts — also
kann kein anderer Modus sie erreichen. Der Preis ist der, den die Schräge
vermeiden sollte: das Austreten bewegt sich mit dem Betrag, also wird die Figur
überall dort geschnitten, wo der Fluchtradius in eine Zelle fällt. Kombinieren
Sie es mit einem kreisförmigen Fluchtradius und einer Färbung, die die
Komponenten getrennt liest, oder nehmen Sie das Schneiden absichtlich in Kauf.

Eine Schräge von null lässt jeden Modus wirkungslos, also ändert `skew_mode`
allein kein Bild.

Zuerst wurde eine gerade Rampe geschrieben — `skew_re*du + skew_im*dv` über die
Zelle —, und die ersetzt dies. Sie färbte, zeichnete aber dieselben diagonalen,
senkrechten oder waagerechten Bänder in jede Zelle gleich, wie auch immer das
Feld geteilt war, und die Bänder verdarben das Bild.

Bei **`0`**, was ein Aufruf bekommt, der es nicht nennt, ist die Drehung null,
der Imaginärteil bleibt null, und jeder Wert ist die Zahl, die er immer war,
**bis aufs Bit** — die Referenz-Prüfsummen von vor dem Argument laufen
unverändert durch, und eine gespeicherte Position wird identisch berechnet. Die
Position `talc`, drei Rauschaufrufe und ein `c` über 345600 Werte, ergibt
dieselbe Signatur wie der Build vor dem Argument.

**Wie viel Farbe sie gibt.** Der Index der Engine ist
`(iter + quantity) * speed + shift`, also ist eine Streuung von eins in der
Größe ein Band wert bei Geschwindigkeit eins. Die Streuung von `imag` innerhalb
einer Zelle:

| Schräge | 0.02 | 0.05 | 0.1 | 0.2 | 0.4 | 0.6 | 1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `randsc`, `randsct` | 0.08 | 0.21 | 0.43 | 0.83 | 1.5 | 1.9 | 2.1 |
| die drei Mosaike | 0.06 | 0.14 | 0.28 | 0.53 | 0.96 | 1.2 | 1.4 |

Also erscheint **bei `0.3` bis `0.6` ein Band bei Geschwindigkeit eins**, und
das fehlte diesen Feldern. Verschiedene Werte zu zählen, würde sagen, `0.01`
genügt — neunhundert davon —, aber tausend Werte innerhalb eines Hundertstels
eines Bandes sind alle eine Farbe.

**Was sie kostet.** Skalieren bewegte den Betrag des Wertes, und ein runder
Fluchtradius schaut genau darauf: eine Zelle, deren Pegel nahe dem lag, wo der
Fluchtradius fällt, kam in zwei geteilt heraus, und die Mosaike verloren ihre
Formen. Eine Drehung lässt den Betrag, wo er ist, also ist **die Schräge bei
kreisförmigem Fluchtradius in jeder Stärke kostenlos** — über 90000 Pixel, bei
0.05, 0.3 und 1 und über alle fünf Felder: kein einziges Pixel tritt anders
aus.

Ein Fluchtradius-Vieleck ist eine andere Sache, denn es liest nicht den Betrag,
sondern die Komponenten, und eine Drehung führt den Wert auf einem Kreis herum,
der eine Seite kreuzen kann. Wie weit er sie kreuzen kann, ist, wie weit die
Ecken der Form über ihre Seiten hinausragen, und die Zahlen folgen dem — der
Anteil des Bildes, dessen Austreten sich ändert, über die fünf Felder:

| Fluchtform | Ecken zu Inkreisradius | Schräge 0.01 | Schräge 0.05 | Schräge 0.3 | Schräge 1 |
| --- | --- | --- | --- | --- | --- |
| Kreis | 1.00 | **0%** | **0%** | **0%** | **0%** |
| Sechseck | 1.15 | **0%** | 0.0–0.8% | 0.7–2.2% | 0.4–2.4% |
| Quadrat | 1.41 | **0%** | 0.0–1.0% | 0.8–2.8% | 1.2–12.8% |
| Dreieck | 2.00 | 0.05–1.3% | 0.7–5.0% | 2.6–6.4% | 7.4–21.6% |

Die **Zahl**, auf die der Fluchtradius gesetzt ist, verschiebt es ebenfalls,
und aus demselben Grund: die Seite wechseln kann ein Wert, dessen Betrag
zwischen dem Inkreisradius des Vielecks und seinen Ecken liegt, und wo dieser
Ring zwischen den Zellen liegt, entscheidet die Zahl.

Also ziehen die beiden Wünsche bei einem Vieleck gegeneinander, und es lässt
sich nicht anders einrichten: das Austreten liest die beiden Komponenten des
Wertes und die Färbung auch, also ist alles, was der Farbe etwas zu lesen gibt,
etwas, das das Austreten sehen kann. **Kombinieren Sie die Schräge mit einem
kreisförmigen Fluchtradius**, und sie kostet in jeder Stärke nichts. Bei einem
Vieleck nehmen Sie entweder die paar Prozent hin oder halten die Schräge bei
`0.05` und erhöhen stattdessen die Färbungsgeschwindigkeit.

Was eine Drehung nicht berühren kann, ist **`zmag`**, das den Betrag liest.
Färben Sie mit etwas, das die beiden Komponenten getrennt liest: `iter + real`,
`iter + imag`, `angle`, `real / imag` außen und `real` oder `real / imag` in
der inneren Färbung.

Daraus folgen zwei Dinge. Der Wert ist komplex, solange die Schräge nicht null
ist, also dreht `randsc(7,,,,,0.4)*z` `z` ebenso, wie es es skaliert. Und das
Kaleidoskop faltet in beiden Fällen: gemessen über fünf Felder, zwei, drei,
fünf und sechs Keile und beide Spiegel, lässt eine Drehung um einen Keil jeden
Wert, wo er war.

**`randscq(...)`** — dasselbe Feld ohne die Interpolation: ein Mosaik flacher
quadratischer Zellen statt Flecken. Gleiche Argumente, gleiche Bedeutung.

**`randscp(...)`** — wieder dasselbe Feld, die Kurven entfernt, aber nicht die
Unregelmäßigkeit: in jeder Zelle wird ein Keim verstreut, und jede Position
nimmt den Wert des nächsten Keims, was ein Voronoi-Diagramm zeichnet — flache
konvexe Vielecke mit geraden Kanten, keine zwei von derselben Form.

**`randsch(...)`** — dasselbe noch einmal auf Sechsecken: eine Wabe flacher
Zellen. Von den regelmäßigen Vielecken, die die Ebene kacheln, ist dies das
ohne Vorzugsrichtung — jede Zelle hat dieselben sechs Nachbarn unter denselben
sechs Winkeln —, also liest es sich als Material und nicht als Gitter.

**`randsct(...)`** — gleichseitige Dreiecke mit abwechselnder Ausrichtung. Ein
Dreiecksmosaik hat eine Vorzugsrichtung, und genau die verlangt man, wenn man
es wählt.

**`randsctile(tiling, seed, ...)`** — dasselbe Feld über einer von
**fünfundvierzig Kachelungen**, die das erste Argument wählt; alles danach ist,
was die Familie nimmt, eine Stelle weiter, also ist `selfsim` das neunte. Jede
Kachelung ist im Mittel auf eine Kachel der Fläche eins skaliert, also ändert
das erste Argument die Form der Zellen und nicht ihre Größe. Eine Zahl
außerhalb von 1 bis 45 zeichnet null.

| | |
| --- | --- |
| 1–3 | die regelmäßigen Kachelungen: Quadrate, Dreiecke, Sechsecke |
| 4–11 | die acht archimedischen, regelmäßige Vielecke mit lauter gleichen Ecken: 4.8.8 (Achtecke und Quadrate), 3.6.3.6, 3.4.6.4, 3.12.12, 4.6.12, 3.3.3.4.4, 3.3.4.3.4 (abgeschrägte quadratische), 3.3.3.3.6 (abgeschrägte hexagonale) |
| 12–19 | ihre dualen, die nicht regelmäßig sind: Tetrakis-Quadrat, Rhombille (die Würfel), deltoidale trihexagonale (Drachenvierecke), Triakis-Dreieck, Kisrhombille und die prismatischen, die Kairo- und die Blütenfünfecke |
| 20–32 | Ziegel, flämischer Verband, Fischgrät, Korbgeflecht, pythagoreisch (Quadrate zweier Größen), Winkel, Quadrate und Rauten, Häuschen (ein Fünfeck, abwechselnd stehend und hängend), Reihen aus Quadraten und Dreiecken in zwei Rhythmen, Sechsecke zwischen Dreiecken, griechische Kreuze, T-Tetrominos |
| 33–37 | islamische Sterne: achtzackig mit Kreuzen, sechszackig mit Sechsecken, achtzackig aus 4.8.8, zwölfzackig aus 3.12.12 und aus 4.6.12 |
| 38 | Voronoi-Zellen, unregelmäßig — die von `randscp`, mit einem eigenen Wert |
| 39–43 | Kachelungen, die sich nie wiederholen, aus de Bruijns Multigittern: Penrose-Rauten, Penrose-Drachen und -Pfeile, die Quadrate und Rauten von Ammann–Beenker und Rauten in zwölf und in sieben Richtungen |
| 44–45 | durch Substitution: das Windrad, dessen Dreiecke in jede mögliche Richtung zeigen, und der Stuhl, Ls in Ls geteilt |

Der **Reiter Kachelungen** von Hilfe → Referenz der Benutzerformeln hat ein
Bild jeder, beschriftet so, wie der Reiter Werte sie aufführt, und ein
Doppelklick auf eines kopiert den Anfang seines Aufrufs — `randsctile(18, ` —
in die Zwischenablage, bereit für den Anfangswert. Die Bilder zeichnet
`randsctile` selbst, acht Einheiten im Quadrat um den Ursprung, mit
`tools/randsctile-thumbnails.cpp`, und sie werden im Programm aufbewahrt, statt
bei jedem Öffnen des Reiters gezeichnet zu werden. Aufbewahrte Bilder können
veralten, also wird jedes mit einem Fingerabdruck dessen gespeichert, was die
Funktion zeichnete, und ein Test zeichnet dieselben Punkte neu und schlägt
fehl, mit Angabe der Kachelung, bis
`cmake --build <build directory> --target randsctile-thumbnails` sie neu
gezeichnet hat.

Die periodischen sind **Tabellen**, geschrieben von
`tools/randsctile-tables.py`, das jede Kachelung aus ihrer Geometrie baut und
prüft — zwanzigtausend zufällige Punkte, von denen jeder in genau eine Kachel
fallen muss —, bevor es sie schreibt. Die islamischen Sterne entstehen aus den
archimedischen Kachelungen: jedes Vieleck mit acht oder mehr Seiten wird zu
einem Stern mit seinen Spitzen auf den Ecken des Vielecks, und das schmale
Dreieck zwischen jeder Kerbe des Sterns und der Seite geht an die Kachel
jenseits dieser Seite, also bilden zwei Sterne mit gemeinsamer Seite eine Raute
zwischen sich, und ein kleines Vieleck bekommt eine Kerbe an jeder Seite, die
es mit einem Stern teilt. Kacheln dürfen konkav sein — Sterne, Kreuze, das T —,
also wird ein Punkt über die Kreuzungszahl in seine Kachel gelegt und nicht
Seite für Seite. Der Punkt wird in der Arbeitsgenauigkeit in ganze Perioden
zerlegt, wie `randscq` ihn in Zellen zerlegt, also ist eine periodische
Kachelung beliebig weit draußen exakt.

Die, die sich nie wiederholen, haben **keine Tabelle** und werden an jedem
Punkt neu berechnet. Die Penrose-Rauten, Ammann–Beenker und die zwölf- und
siebenzähligen Rauten kommen aus **de Bruijns Multigitter**: Scharen paralleler
Linien in einer zweiten Ebene, eine Kachel überall, wo sich zwei kreuzen,
benannt nach den beiden Linien — vier ganze Zahlen, exakt, wer auch immer
fragt. Der Punkt wird gefunden, indem man die Kreuzungen neben der Stelle
ausprobiert, an der er in jener Ebene landet, zuerst die Scharen, deren Linien
am nächsten vorbeigehen, was ihn im Mittel in vier bis acht Versuchen findet.
Die Drachen und Pfeile kommen aus Robinsons Dreiecken, Ebene für Ebene aus
einem Rad von zehn um den Ursprung geschnitten, das Windrad aus einem
rechtwinkligen Dreieck, in fünf geteilt, der Stuhl aus einem L, in vier
geteilt; auf jeder Ebene entscheiden eine oder zwei Linien, in welchem Kind der
Punkt liegt.

**Was es kostet**, gemessen pro Aufruf gegen die 220 ns von `randscq` auf
dieser Maschine: die periodischen Kachelungen 255 bis 380 ns, die, die sich nie
wiederholen, 410 bis 680. Sie werden in beiden Builds in `double` berechnet,
die Stelle innerhalb einer Periode eingeschlossen — in long double brauchten
die nie wiederholenden dreimal so lange —, und dabei zeichneten beide Builds
alle fünfundvierzig Kachelungen pixelgleich, je vierzigtausend Punkte. Jenseits
von 2³² Zellen vom Ursprung, wo `double` anfinge, ihre Kacheln zu verlieren,
werden die nie wiederholenden flach, wie die Familie jenseits ihres Gitters;
die periodischen halten so weit durch wie der Rest der Familie.

Die Schräge misst eine Kachel so, wie `randscp` sein Vieleck misst, am Radius
des größten Kreises, den die Kachel enthält, also sind die Konturen die
verkleinerte Kachel, konkave eingeschlossen; die Rosette dreht um die Mitte der
Kachel.

    z+randsctile(18,13,{0.3,0.3},{1,1})*0.25        Kairo-Fünfecke
    z+randsctile(36,13,{0.3,0.3},{1,1})*0.25        zwölfzackige Sterne
    z+randsctile(39,13,{0.3,0.3},{1,1},5)*0.25      Penrose, fünffach gefaltet

**`fbm(value, seed, [intensity=4], [frequency=8], [octaves=4], [roughness=0.5])`**
— eine fraktionale Brownsche Bewegung: Gradientenrauschen, in Oktaven summiert,
jede mit der doppelten Frequenz der vorigen und mit `roughness` ihrer Höhe —
dieselbe Bewegung, die die gleichnamigen Färbungsmodi zeichnen, siehe
[Fraktionale Brownsche Bewegung](#fraktionale-brownsche-bewegung). Was sie
zeichnet, ist eher **Abnutzung** als ein Muster — Flecken und Dellen —, weil
keine Oktave groß genug ist, um für sich gesehen zu werden, und keine klein
genug, um zu verschwinden.

Wo `randsc` die Position liest und nur die Position lesen kann, **liest dies,
was auch immer Sie davor schreiben**. Das ist der ganze Sinn, es als Funktion
zu haben:

    fbm(z,7)          bewegt sich mit dem Orbit
    fbm(x,7)          steht still auf der Ebene
    fbm(z*3+c,7)      was auch immer das ist
    z^2+c+fbm(z,7,0.1)*i    die Menge mit einem rauen Rand

Wert und Anfangswert sind Pflicht; die übrigen haben Standardwerte. Die
Bewegung läuft **von null bis `intensity` und nie darunter**, also kann das
Hinzufügen zu etwas dieses nicht unter null ziehen. Der Anfangswert wird
gelesen wie bei `randsc`, also bedeutet dieselbe Zahl in beiden dasselbe Feld.
`octaves` wird zwischen eins und vierundzwanzig gehalten, und es lohnt sich nur
zu erhöhen, wenn `roughness` erhöht wird: bei der einfachen Hälfte trägt die
achte Oktave ein Zweihundertfünfzigstel des Ganzen.

**Es hat kein eigenes Kaleidoskop**, anders als die `randsc`-Familie, weil
schon eines da ist, mit dem man es zusammensetzen kann: `parchment(z, sectors)`
faltet den Winkel von `z` in `|sectors|` Sektoren, und `parchmenta` spiegelt
die abwechselnden Hälften, also sagt `fbm(parchmenta(z,6),7)` es, und zwar
dort, wo der Leser der Formel es sehen kann.

Die Färbungsmodi unter **Fraktal → Äußere Färbung → Fraktionale Brownsche
Bewegung...** wenden dieselbe Bewegung auf eine Farbe an. Dies hier legt sie
dahin, wo eine Formel sie benutzen kann.

| | |
| --- | --- |
| `randsc` | weich, rund |
| `randscq` | hart, regelmäßig, quadratisch |
| `randscp` | hart, unregelmäßig |
| `randsch` | hart, regelmäßig, sechseckig |
| `randsct` | hart, regelmäßig, dreieckig |
| `randsctile` | hart, jede von fünfundvierzig Kachelungen |

Alle sechs legen eine Zelle über jedes Einheitsquadrat der geltenden Größe,
also bedeutet `size` überall dasselbe: einen Buchstaben zu ändern, ändert die
Form der Zellen und nicht den Maßstab des Bildes. Ein Sechseck mit
Umkreisradius eins bedeckt 2.6 Einheitsquadrate und ein Dreieck mit Seitenlänge
eins 0.43, also werden diese beiden Gitter passend zu den übrigen skaliert,
statt sie zu lassen, wie sie kommen.

Jeder Durchlauf bekommt ein eigenes Feld. Die Iteration geht zusammen mit der
Position in den Hash ein, also bekommt eine Formel, die eine davon einmal pro
Iteration aufruft, jedes Mal einen neuen Wert — auch bei einer Degradation von
eins, die die Größe unverändert lässt. Die Degradation legt fest, wie groß die
Flecken in einem gegebenen Durchlauf sind; sie ist nicht das, was den Durchlauf
anders macht.

Unter der Größe gibt es eine Grenze. Eine Zelle findet man, indem man die
Position durch die Größe teilt, und das Ergebnis muss in eine ganze Zahl
passen, also kann die Größe sinnvoll nicht unter etwa die Position geteilt
durch das, was eine ganze Zahl fasst, gehen — ein Milliardstel eines
Milliardstels davon. Eine Degradation von einhalb erreicht das in etwa sechzig
Durchläufen. Danach bleibt keine Zellstruktur mehr aufzulösen, und das Feld ist
ein einziger flacher Wert über die ganze Ebene; es ändert sich weiterhin in
jedem Durchlauf und unterscheidet sich weiterhin zwischen den sechs Funktionen,
also bleibt eine Formel, die eine von einer anderen abzieht, nicht bei null
stehen und iteriert umsonst bis zum Limit. Für ein Ausblenden, das bis zum Ende
ein Bild bleibt, nehmen Sie eine Degradation nahe eins.

Alle hashen die Position des Punktes und die Iteration, nie `z` und nie einen
globalen Zustand. Also kommt bei jedem Neuzeichnen dasselbe Bild zurück, bei
jeder Thread-Anzahl und im Mandelbrot- wie im Juliamodus. `randsc` ist stetig,
also stimmen die beiden Genauigkeiten bis etwa 1e-19 überein; die Mosaike sind
Treppenfunktionen und können sich auf einer Haarlinie entlang der Zellkanten
unterscheiden, was untrennbar dazugehört, wenn man harte Kanten verlangt.

Eine Formel, die eine davon aufruft, schaltet das Boundary-Trace aus. Diese
Optimierung füllt eine Region, die sie von einer einzigen Farbe umgeben fand,
ohne sie zu berechnen — wahr für ein Fraktal, falsch für ein Rauschfeld.

Das alte `rand` ist unverändert und hängt weiterhin davon ab, wie oft es
aufgerufen wurde, also zeichnet sich dieselbe Position nicht gleich neu.
Bevorzugen Sie `randsc`.

### Kaleidoskope

Die beiden Argumente nach der Degradation falten die Ebene, bevor das Feld auf
ihr abgetastet wird. Das erste gibt an, in wie viele Keile sie gefaltet wird,
wobei `1` — der Standard — sie unverändert lässt; das zweite, `kaled_mirror`,
welcher Spiegel die Faltung vornimmt:

| kaled_mirror | |
| --- | --- |
| `0` (und jeder andere Wert) | die ferne Hälfte jedes Keils spiegelt die nahe, also ist jeder Keil symmetrisch zu seiner eigenen Winkelhalbierenden |
| `1` | dasselbe andersherum, die nahe Hälfte spiegelt die ferne |

    randsc(13,{0.6,0.6},{1,1},6,0)     sechs Keile, jeder ein Spiegel seiner selbst
    randsc(13,{0.6,0.6},{1,1},3,1)     drei Keile, andersherum gefaltet

Beide Faltungen sind stetig, wo die Keile aufeinandertreffen, also bleibt das
Rauschen zusammenhängend und die beiden Genauigkeiten stimmen weiter überein;
eine Faltung, die sich ungleichmäßig träfe, wäre als Naht zu sehen.

Das ist der einzige Teil der Familie, der Trigonometrie kostet, und nur ein
Wert von zwei oder mehr erreicht ihn. Ein Aufruf, der nichts dazu sagt, zahlt
einen Vergleich.

### Ein Bild aus dem Rauschen holen

`randsc` gibt einen Wert in `[0, 1)` zurück. Allein verwendet, verlässt die
Iteration nie einen Fluchtradius von 4, also tritt nichts aus und das Bild ist
flach. Multiplizieren Sie es, oder lassen Sie es eine Iteration stören, die
tatsächlich austritt:

    z+randsc(13,{0.35,0.35},{1,1})*0.25       Flecken
    randsc(13,{0.15,0.15},{0.5,0.5})*2.5      Brownsche Bewegung
    randscq(13,{0.15,0.15},{0.5,0.5})*2.5     verstreute Quadrate, immer kleiner
    z+randscp(13,{0.3,0.3},{1,1})*0.25        unregelmäßige Vielecke
    z+randsch(13,{0.3,0.3},{1,1})*0.25        eine Wabe
    z+randsct(13,{0.3,0.3},{1,1})*0.25        Dreiecke
    z^2+c+randsc(13,{0.25,0.25},{1,1})*1.2    die Menge selbst verformt

Eine Degradation von 0.5 halbiert die Flecken in jedem Durchlauf, also sind sie
nach zwanzig Iterationen ein Millionstel ihrer Größe und kleiner als ein Pixel.
Für ein langsames Ausblenden über viele Iterationen nehmen Sie etwas nahe 1 —
`0.97^60` ist immer noch 0.16.

Ein ganzzahliger Anfangswert ist in beiden Builds exakt; ein gebrochener wird
auf seine führenden Bits quantisiert, die bis auf etwa einen Anfangswert von
fünfzig Millionen übereinstimmen.

### Alle Durchläufe bisher: `selfsim`

Jeder Durchlauf dieser Funktionen ist ein eigenes Feld — die Iteration geht in
den Hash ein — mit Zellen, die die Degradation mal die des vorigen Durchlaufs
sind. Eine Formel, die eine davon in jedem Durchlauf aufruft, liest daher jedes
Mal ein neues, feineres Feld, und ein Färbungsmodus, der den letzten Durchlauf
liest, liest nur das letzte Feld. Bei einer Degradation von 0.75 sind die
Zellen einer fünf Einheiten breiten Ansicht beim sechzehnten Durchlauf kleiner
als ein Pixel, und danach ist das Innere Schnee: hundert Durchläufe von
`randsct(672,,{0.75,0.75},,,{2,2},6)+c` im Juliamodus setzen 70 % der
benachbarten inneren Pixel ein Sechzehntel der Palette oder mehr auseinander.

Das achte Argument gibt stattdessen **die bisherigen Durchläufe gemittelt**
zurück, die zueinander stehen wie die Oktaven einer fraktionalen Brownschen
Bewegung: Durchlauf *n* wiegt d^(nH), wobei *d* die Degradation ist — das
geometrische Mittel ihrer beiden Komponenten ohne Vorzeichen, also ist
`0.75+0.75i` gleich 0.75 — und *H* das Argument. Die Gewichte werden durch ihre
Summe geteilt, also bleibt das Ergebnis zwischen den gemittelten Werten: in [0,
1] ohne Schräge.

| selfsim | |
| --- | --- |
| weggelassen | aus — der eine Durchlauf, bis aufs Bit. Was ein Aufruf bekommt, der nichts nennt, und eine leer gelassene Stelle |
| `0` | der einfache Mittelwert: jeder Durchlauf wiegt gleich, und das ergibt jedes *H* nahe null — nicht aus |
| `1` | die einfache Bewegung: jeder Durchlauf wiegt *d*-mal den vorigen |
| `0.5` | rauer: die feinen Durchläufe behalten mehr Gewicht, also zeigt sich mehr Detail |
| größer | glatter: die ersten paar Durchläufe nehmen fast alles |
| negativ | die feinen Durchläufe überwiegen die groben, und der Schnee ist wieder da |
| komplex, `{1,2}` | der Realteil gewichtet wie oben, der Imaginärteil dreht die Durchläufe — siehe unten |

    randsct(672,,{0.75,0.75},,,{2,2},6,1)+c     dasselbe, selbstähnlich

Bei dieser Formel verschwindet der Schnee: benachbarte innere Pixel liegen im
Mittel 1.8 Palettenzellen auseinander statt 35, und 0.6 % von ihnen mehr als
sechzehn statt 70 %. Was bleibt, sind die Dreiecke auf allen Maßstäben
zugleich, die großen tragen die kleinen.

* **`selfsim` darf komplex sein.** Dann ist auch d^(nH) komplex: d^(n·Hr) im
  Betrag und eine Drehung um n·Hi·ln d, also wird jeder Durchlauf um so viel
  weiter gedreht als der vorige und dann gemittelt, wie ein reelles *H* mittelt
  — eine Spirale von Oktaven. Bei 0.75 sind das bei einem Imaginärteil von eins
  etwa sechzehn Grad pro Durchlauf. Die Gewichte werden durch **die Summe ihrer
  Beträge** geteilt, nicht durch ihre komplexe Summe: so bleibt das Ergebnis
  innerhalb des größten der gemittelten Werte. Die komplexe Summe käme bei
  irgendeinem Durchlauf nahe null und würde bei diesem Durchlauf jeden Wert für
  alle Pixel zugleich aufblähen, da die Gewichte nichts mit dem Pixel zu tun
  haben. Die Drehungen geben dem Wert auch ohne Schräge einen Imaginärteil,
  also haben `imag` und `angle` etwas zu lesen. Eine Degradation von eins hat
  den Logarithmus null und dreht nichts; `{1,0}` ist `1` bis aufs Bit; und
  `{0,2}` ist jeder Durchlauf gleich gewichtet, und gedreht.
* **Null ist nicht aus.** d^(n·0) ist für jeden Durchlauf eins, also ist `0`
  der einfache Mittelwert, und ebenso jedes *H* in seiner Nähe — von beiden
  Seiten und auch komplex. Zuerst war die Null der Schalter, was sie zum
  einzigen Wert machte, durch den die Kurve nicht geht: `0.000001` war der
  einfache Mittelwert und `0` nur der letzte Durchlauf, ein Viertel des
  Bereichs voneinander entfernt. Nur der letzte Durchlauf ist das ferne Ende
  der Kurve, wenn *H* gegen minus unendlich läuft. Um keinen Mittelwert zu
  haben, lassen Sie das Argument weg oder seine Stelle leer.
* **Gemittelt wird der geschrägte Wert**, also tun die Schräge und jedes Bit
  von `skew_mode` weiter, was sie taten, und das Ergebnis ist komplex, wenn die
  Schräge es dazu macht.
* **Ein Mittelwert zieht die Werte zusammen.** Bei 0.75 und *H* gleich eins ist
  die Streuung etwa zwei Fünftel der eines einzelnen Durchlaufs, also erhöhen
  Sie die Färbungsgeschwindigkeit, um den Kontrast zurückzuholen.
* **Eine Degradation von eins** gewichtet jeden Durchlauf gleich und ergibt
  ihren einfachen Mittelwert, der sich mit zunehmender Zahl der Durchläufe zu
  einer flachen Mitte hin setzt. Das sagt die Definition, und so wird es getan.
* **Durchläufe, die ein Aufruf nicht gesehen hat, werden berechnet, wenn er
  gefragt wird.** Ein Aufruf in einem Zweig, den die Formel erst ab dem siebten
  Durchlauf nimmt, bekommt trotzdem die ersten sieben, und das Ergebnis bei
  einem Durchlauf ist dasselbe, auf welchem Weg man auch dorthin kam. Ein neuer
  Pixel wird an seiner Position wie an seinem Durchlauf erkannt: der Durchlauf
  allein, nach dem `trap` geht, ließe einen solchen Aufruf den Mittelwert
  fortsetzen, den der vorige Pixel hinterließ.
* **Nichts läuft über.** Der Mittelwert wird als einer geführt, wobei jeder
  Durchlauf einen Anteil daran nimmt, der aus dem Anteil des vorigen folgt,
  statt als Summe und Gesamtgewicht — also hat ein negatives *H* oder eine
  Degradation über eins, deren Gewichte über das hinauslaufen würden, was der
  Typ fasst, keine Gelegenheit dazu.
* **Es kostet** etwa ein Zehntel mehr pro Aufruf bei 64 Bit Mantisse und ein
  Siebtel bei 113. Ausgeschaltet nichts Messbares.

### Figuren statt Rauschen

`sierpinskyt`, `sierpinskyc` und `snowflake` stehen dort, wo das Rauschen
steht, und nichts an ihnen ist zufällig: die Figur ist dort, wo sie zu sein
scheint, und Mandelbrotmodus und Juliamodus zeichnen dasselbe Bild, da eine
Form in der Ebene so oder so am selben Ort ist.

**Sie sind eigene Fraktale, keine Felder, die man in eines
hineinmultipliziert.** Allein geschrieben — `sierpinskyt()` und sonst nichts —
zeichnet jede ihre Figur, so wie Fraktal → Mehr Formeln sein **Sierpinski**,
**Sierpinski Carpet** und **Koch Snowflake** zeichnet.

Alle drei tun es, indem sie **den Punkt zu seinem Elternteil tragen**. Jeder
Teil einer dieser Figuren hat einen: ein Loch im Sierpinski-Dreieck liegt in
einem größeren Loch eine Ebene höher, ein Loch im Teppich in der Zelle, die auf
dieselbe Weise zerschnitten wurde, ein Dreieck der Schneeflocke auf einer Seite
des Sechsecks in der Mitte oder auf der freien Kante eines anderen Dreiecks.
Der Schritt auf den Elternteil ist eine Bewegung der Ebene — Verdoppeln weg von
der nächsten Ecke beim Dreieck, Vergrößern einer Zelle um die Zahl der Zellen
beim Teppich. Der oberste Teil hat innerhalb der Figur keinen Elternteil, also
trägt ihn sein Schritt aus dem Fluchtradius hinaus, und ein Teil *n* Ebenen
tiefer braucht *n* Schritte, um dorthin zu gelangen: **der Durchlauf, in dem
ein Punkt austritt, ist die Ebene, auf der er steht**, und die Iterationszahl
ist das Bild.

Was in dem Durchlauf zurückkommt, in dem ein Punkt austritt, zählt so viel wie
der Durchlauf selbst: **jeder äußere Färbungsmodus außer der Iterationszahl
liest es**. Also wird der oberste Teil jeder Figur hinausgeworfen, statt
einfach für ausgetreten erklärt zu werden — das Dreieck verdoppelt sein
mittleres Loch weg von einer Ecke, und es landet außerhalb des Dreiecks, der
Teppich skaliert seine mittlere Zelle um die Nachbarzelle, die Schneeflocke
wirft ihr mittleres Sechseck in die Richtung hinaus, in die dieser Sektor
zeigt. Jeder nimmt mit, woher er kam, also haben `real`, `imag`, `angle`, die
weiche Färbung und der Rest alle etwas, womit sie arbeiten können. Nur ein
Punkt, der nie in der Figur war, bekommt eine Zahl, die für „ausgetreten“
steht.

Eine Schneeflocke wird **von einem Sechseck aus nach außen** gelesen, und das
muss man über sie wissen. Sie zerfällt exakt in ein regelmäßiges Sechseck in
der Mitte, sechs Dreiecke mit der Seitenlänge dieses Sechsecks auf seinen sechs
Seiten, zwölf mit einem Drittel davon auf deren freien Kanten, achtundvierzig
mit einem Neuntel auf deren, und so fort ohne Ende. Drei der sechs sind die
Ecken des Dreiecks, aus dem die Figur wuchs, und drei sind die ersten Zacken
auf seinen Kanten, und **nichts unterscheidet sie**: vom Sechseck aus nach
außen hat eine Schneeflocke sechszählige Symmetrie.

Liest man stattdessen vom ersten Dreieck aus, sind fünf Achtel der Figur Ebene
eins — ein einziges flaches Dreieck, das das Bild füllt, die Schneeflocke nur
im Saum darum, und nichts, was die Außenfarbe sagen könnte. Vom Sechseck aus
gelesen sind die ersten beiden Ebenen je fünf Zwölftel, und die übrigen nehmen
gleichmäßig ab. Ein Dreieck auf dem Sechseck geht auf das größte Dreieck über,
das das Sechseck enthält — Ecke auf Ecke, die dreifache Fläche —, wird also um
die Wurzel aus drei vergrößert und um ein Zwölftel einer Umdrehung gedreht;
eines auf der freien Kante eines anderen Dreiecks wird um drei vergrößert und
in die Richtung gedreht, in die diese Kante zeigt.

**Jeder Schritt dehnt**, und das ist eine Regel und kein Zufall. Das Dreieck
verdoppelt, der Teppich vergrößert eine Zelle um die Zahl der Zellen, die
Schneeflocke um drei oder um die Wurzel aus drei. Ein Dreieck auf dem Sechseck
könnte stattdessen über die Seite, auf der es steht, flach geklappt werden, was
es auf eines der sechs legt, aus denen das Sechseck besteht, und ordentlicher
ist; aber ein Klappen ist eine Isometrie, und eine Isometrie ist unempfindlich
dafür, wo der Punkt anfing. In eine größere Formel geschrieben, gäbe die Figur
dieser Formel dann auf fünf Zwölfteln ihrer Fläche ein starres Abbild ihrer
selbst, und die Formel zeichnete dort eine flache Region. Ein Test prüft die
drei Dehnungen.

Der Punkt reist wirklich, und das macht sie brauchbar statt nur korrekt: `z`
ist unterwegs ein Punkt der Ebene wie jeder andere, also haben die
Färbungsmodi, die ihn lesen, die weiche Färbung und das Schreiben einer Figur
in eine größere Formel alle einen Sinn.

Ein Sierpinski-Dreieck und ein Teppich füllen ihre Form, also hat jeder Punkt
darin eine Ebene und tritt auf ihr aus. Eine Schneeflocke füllt ihr Sechseck
nicht: sie lässt sechs Ecken **Grund** übrig, und dieser Grund ist kein Teil
der Figur und hat keine Ebene. Er wird **um eine halbe Umdrehung gedreht**, was
Grund auf Grund legt, also tritt er nie aus und erscheint in der **Innenfarbe**
— der leere Raum bleibt leer, statt der Figur ein Band der Außenfarbe
wegzunehmen.

Gedreht, und nicht dort zurückgegeben, wo er stand, denn **ein Fixpunkt ist ein
angehaltener Orbit und kein beschränkter**. Stillstehend lief die Schleife
jeden Durchlauf bis zum Limit, ohne dass sich etwas änderte; die inneren
Färbungsmodi, die einen Durchlauf mit dem vorigen vergleichen, hatten nichts zu
vergleichen; und `snowflake()` in einer größeren Formel gab dieser Formel den
Punkt zurück, den sie schon hielt, während `sierpinskyt()` und `sierpinskyc()`
jeden Punkt bewegen, den sie bekommen. Die halbe Umdrehung wird dem Sechstel
vorgezogen, das die Symmetrie der Figur ebenfalls erlauben würde: beide legen
Grund auf Grund, aber die halbe Umdrehung hält den Punkt innerhalb eines
Quadrats und eines Kreises ebenso wie innerhalb des Sechsecks. Sie lässt den
Abstand zur Mitte unverändert, also sind die inneren Färbungsmodi, die ihn
lesen, nicht berührt, und das Bild, das `snowflake()` zeichnet, ist pixelgenau
unverändert.

Das ist die einzige Stelle, an der diese Figuren etwas kosten. Einem Pixel des
Grundes wird die Zugehörigkeitsfrage einmal pro Durchlauf gestellt, solange man
hinsieht, während die beiden anderen Figuren und der Körper der Schneeflocke in
dem Durchlauf austreten, den ihre Ebene nennt. Direkt durch die
Iterationsschleife gezeichnet, kostet eine ganze Schneeflocke etwa das
2.4-Fache dessen, was sie kostete, als der Grund ein eigenes Band bekam. Das
Boundary-Trace bleibt für diese Formeln eingeschaltet — ihre Bänder sind
wirklich massiv, und das ist die Bedingung dafür —, und der Grund ist eine
einzige massive Region, also ist, was auf den Bildschirm kommt, weniger, als
diese Zahl nahelegt.

In eine größere Formel geschrieben — `snowflake()^-2+snowflake()` und
dergleichen — verhalten sich die drei **nicht** gleich, und der Grund ist es
wert, deutlich gesagt zu werden, weil keine Arbeit an der Figur ihn ändern
wird.

Ein Sierpinski-Dreieck ist ein **Staub**: es hat keine Fläche, seine Punkte
treten nie aus, und jeder von ihnen wird in jedem Durchlauf verdoppelt. Eine
darauf gebaute Formel hat überall und in jeder Tiefe etwas Chaotisches unter
sich. Eine Schneeflocke ist **massiv**: zehn Zwölftel von ihr sind das Sechseck
und die sechs Dreiecke, die im ersten und zweiten Durchlauf austreten, also
geht einer darauf gebauten Formel die Figur fast sofort aus, und sie zeichnet,
was sie von sich aus zeichnet. Das Detail liegt entlang des Koch-Randes, wo die
Ebenen tief reichen, und die flachen Regionen sind die Innenflächen, wo sie es
nicht tun.

**Was ein Durchlauf weiß, ist, wofür dieser Durchlauf bezahlt hat.** Die
Schneeflocke entscheidet, ob ein Punkt Grund ist, indem sie die Koch-Kurve
hinabsteigt, und dieser Abstieg darf **eine Ebene pro Durchlauf** gehen. Früher
lief er alle vierundzwanzig Ebenen in einer einzigen Auswertung, also bekam
eine um die Figur geschriebene Formel im ersten Durchlauf den ganzen Rand und
zeichnete ihn vollständig, wie wenige Durchläufe man ihr auch gab — bei
`Iterationen: 2` bekam man immer noch jede Zacke der Koch-Kurve, was keine
andere Formel in XaoS tut und wofür eine Iterationszahl nicht da ist. Dass die
Ebenen ausgehen, ist nicht dasselbe wie draußen zu sein: ein Punkt außerhalb
des Rahmens oder über der Hülle ist Grund und jetzt schon als Grund bekannt,
während ein unentschiedener als Teil der Figur weitergetragen und im nächsten
Durchlauf mit einer Ebene mehr erneut gefragt wird. Das Bild, das `snowflake()`
selbst zeichnet, ist unverändert — ein Punkt in der Figur brauchte nie den
tiefen Abstieg, nur die Gewissheit, dass er nicht Grund ist.

Was bleibt, ist der ehrliche Unterschied zwischen den beiden Figuren, und keine
Arbeit an einer von ihnen wird ihn schließen. Ein Sierpinski-Dreieck ist ein
**Staub**: keine Fläche, seine Punkte treten nie aus, jeder in jedem Durchlauf
verdoppelt, also hat eine darauf gebaute Formel in jeder Tiefe etwas
Chaotisches unter sich. Eine Schneeflocke ist **massiv**, und zehn Zwölftel von
ihr sind das Sechseck und die sechs Dreiecke, die im ersten und zweiten
Durchlauf austreten. Gemessen bei einem einzigen Zoom — dem Anteil der Pixel,
die auf einer Bandkante stehen — erreicht `sierpinskyt()^-2+sierpinskyt()` nach
zwölf Durchläufen 17 Prozent und `snowflake()^-2+snowflake()` 2.5; beide
wachsen mit der Zahl der Durchläufe, worum es geht, und wo sie stehen bleiben,
ist das, was die Figuren sind.

`radius` bedeutet, was `bailout` bedeutet, und wird genauso gelesen: **als
Quadrat des Abstands**. Was es zeichnet, ist die Form, die ein Fluchtradius
dieser Zahl zeichnet, **ihr exakt einbeschrieben** — stellen Sie den
Fluchtradius auf die passende Form und dieselbe Zahl, und Figur und
Fluchtradius liegen übereinander:

| | |
| --- | --- |
| `sierpinskyt` | das Dreieck des Fluchtradius **Dreieck −90°**, Ecke auf Ecke |
| `sierpinskyc` | das Quadrat des Fluchtradius **Quadrat**, Seite auf Seite |
| `snowflake` | das Sechseck des Fluchtradius **Sechseck 0°**, seine sechs Spitzen auf den sechs Ecken des Sechsecks |

Wissen muss man, dass ein Fluchtradius-Vieleck seine **Seiten** in der
Quadratwurzel des Fluchtradius von der Mitte stellt — das ist sein
Inkreisradius, nicht sein Umkreisradius —, also liegen seine Ecken weiter
draußen, als die Zahl sagt: doppelt so weit bei einem Dreieck und um ein
Siebtel weiter bei einem Sechseck. `sierpinskyt(4)` reicht daher mit seinen
Ecken bis 4, während `sierpinskyc(4)` mit seinen Seiten bis 2 reicht, und beide
sind richtig.

Stellen Sie den Fluchtradius größer als den Radius, und die Figur füllt die
Form nicht mehr: zwischen beiden liegt ein **Zwischenraum**, und stellen Sie
ihn kleiner, wird die Figur abgeschnitten. Beides ist, was verlangt wurde, und
keines wird korrigiert. Aber der Zwischenraum ist kein Teil der Figur und hat
keine Ebene, also bekommt er keine: er wird **nach der eigenen Symmetrie der
Figur gedreht** — um ein Drittel für ein Dreieck, um ein Viertel für einen
Teppich, um eine halbe Umdrehung für eine Schneeflocke —, was ihn außerhalb der
Figur und innerhalb des Fluchtradius hält, also tritt er nie aus und wird in
der **Innenfarbe** gezeichnet, wobei die inneren Färbungsmodi auf ihm arbeiten.
Früher trat er im ersten Durchlauf aus und nahm das unterste Band der
Außenfarbe, was den ganzen Zwischenraum in einem einzigen flachen Ton zeichnete
und der Figur ein Band wegnahm.

Eine Drehung und kein Spiegel, obwohl ein Spiegel ebenso dazu taugen würde, ihn
dort zu halten: was eine Figur zurückgibt, muss sich mit dem Punkt drehen,
damit das Drehen des Bildes dreht, was darauf gezeichnet ist, und ein Spiegel
dreht es in die andere Richtung. Aus demselben Grund wird die Mitte eines
Teppichs — der Block von Zellen, den der Randring umschließt, der keinen
Elternteil hat und daher hinausgeworfen werden muss — in die Richtung geworfen,
in die sein eigenes **Viertel** des Quadrats zeigt, die Viertel, die von den
Diagonalen abgeschnitten werden. Immer in dieselbe Richtung geworfen, wie
früher, zeichnete sie diese Mitte als Rampe von einer Seite zur anderen ohne
etwas quer dazu; das Quadrat hat vierzählige Symmetrie, und das sollte auch
haben, was darauf gezeichnet wird. Ein Test prüft, dass jede Figur eine Drehung
ihrer eigenen Ordnung mit derselben Drehung beantwortet.

**`sierpinskyt([radius=4], [kaleidoscope=1], [kaled_mirror=0])`** — das
Sierpinski-Dreieck in diesem Dreieck.

**`sierpinskyc([radius=4], [squares=3], [kaleidoscope=1], [kaled_mirror=0])`**
— der Sierpinski-Teppich in diesem Quadrat. Das Quadrat wird in `squares` mal
`squares` zerschnitten, der Ring der Zellen entlang des Randes bleibt, alles,
was dieser Ring umschließt, wird verworfen, und dasselbe geschieht mit jeder
behaltenen Zelle. Also ist das Bild **ein Quadrat in der Mitte und
4×`squares`−4 darum herum**: acht bei drei, was der Teppich ist, wie er
gewöhnlich gezeichnet wird, zwölf bei vier, sechzehn bei fünf. Zwei ist die
einzige Zahl ohne einen nennenswerten Ring, und dort fällt stattdessen die
ferne Ecke weg, was wieder ein Sierpinski-Dreieck ist — ein in vier geteiltes
Quadrat mit einer weggenommenen Ecke ist genau das, was ein Sierpinski-Dreieck
ist.

**`snowflake([radius=4], [kaleidoscope=1], [kaled_mirror=0])`** — die
Koch-Schneeflocke, ihre Spitzen auf den Ecken dieses Sechsecks, nach Generation
von der Mitte aus in Bänder geteilt: das Sechseck im ersten Durchlauf, die
sechs Dreiecke auf seinen Seiten im zweiten, die zwölf auf deren freien Kanten
im dritten, die achtundvierzig im vierten. Die sechs Ecken Grund, die die Figur
nicht bedeckt, bekommen überhaupt kein Band — sie treten nie aus und werden in
der Innenfarbe gezeichnet.

Jedes Argument hat einen Standardwert, also ist `snowflake()` ein Aufruf, und
ebenso `sierpinskyc( ,5)` — ein filigranerer Teppich in Standardgröße.

**Das Kaleidoskop ist das der Rauschfamilie**, bei allen dreien zuletzt
geschrieben und mit derselben Bedeutung: in wie viele Keile die Ebene
geschnitten wird und welcher Spiegel sie faltet. `1` — was ein Aufruf bekommt,
der nichts sagt — faltet nichts und lässt die Figur die Zahl, die sie immer
war, bis aufs Bit.

Was es faltet, ist jedoch etwas anderes, und muss es sein. Das Rauschen faltet
die **Position**, die für den ganzen Orbit stillsteht; eine Figur kann mit der
Position nichts anfangen und liest **z**, den Punkt, den sie hinab zu dem Teil
trägt, in dem er steht, also wird z gefaltet. Das Bild kommt genauso heraus —
die Figur in einem Keil gezeichnet und um den Ursprung wiederholt —, und es
hält bis ganz hinunter, da jeder Schritt am gefalteten Punkt getan und im
folgenden Durchlauf erneut gefaltet wird. Innerhalb des halben Keils, den die
Faltung unverändert lässt, ist die Figur die Figur.

    sierpinskyt(4,6,0)      das Dreieck sechsmal um den Ursprung
    sierpinskyc(4,3,5,1)    der Teppich in fünf Keilen, der nahe Spiegel
    snowflake(4,5)          fünf Keile Schneeflocke

Eine Schneeflocke hat bereits die sechszählige Symmetrie mit dem Spiegel
entlang jeder Winkelhalbierenden, also ist es die Identität auf dem Bild, sie
so in sechs Keile zu falten: `snowflake(4,6,0)` kommt über 90000 Pixel
**pixelgenau** gleich `snowflake(4)` heraus. Verlangen Sie fünf Keile oder
vier, um sie falten zu sehen.

Eines muss man über den Teppich wissen: was er wegschneidet, wirft er in eine
von vier Richtungen, die von den Diagonalen des Quadrats abgeschnittenen
Viertel, und auf einer Diagonalen selbst ist die Wahl ein Gleichstand. Eine
Faltung kann einen ganzen Strahl der Ebene auf diese Diagonale legen, und dort
werden die beiden Seiten des Gleichstands in entgegengesetzte Richtungen
geworfen — eine Haarlinie entlang der Diagonalen, die die Stufe ist, die der
Teppich schon immer hatte, und nichts, was die Faltung tut. Gemessen über
zehntausend Punkte bei sechs Keilen: jede Abweichung lag auf einer Diagonalen,
und die gefalteten Punkte selbst stimmten bis auf einen Teil in 10^16 überein.

    sierpinskyt()        mit Fraktal -> Form des Fluchtradius -> Dreieck -90
    sierpinskyc()        mit Fluchtradius-Form Quadrat
    sierpinskyc( ,5)     ein filigranerer Teppich
    snowflake()          mit Fluchtradius-Form Sechseck 0

Keine kostet mehr als das Rauschen daneben, Durchlauf für Durchlauf: gemessen
an `randsc` das Dreieck 0.44, der Teppich 0.54 und die Schneeflocke 0.50 davon
bei 64 Bit Mantisse, und 0.50, 0.87 und 0.65 bei 113. Eine Schneeflocke von der
Mitte aus zu lesen, hat sie billiger und zugleich schöner gemacht: sechs
Siebtel der Figur sind das Sechseck oder ein Dreieck darauf, und die werden mit
vier Multiplikationen und einem Vergleich beantwortet, ohne die Kurve überhaupt
abzugehen. Was eine Schneeflocke über ein ganzes Bild kostet, ist der Absatz
oben, und hängt davon ab, wie viele Durchläufe der Grund braucht, nicht davon,
was ein Durchlauf kostet. Das war nicht umsonst, und bei 113 Bit war es fast
verloren — eine Quadratwurzel und eine Division sind dort beide Software, und
jede kostet etwa so viel wie eine ganze Figur. Also wird der Kehrwert der
Quadratwurzel des Radius an der Aufrufstelle gehalten, statt in jedem Durchlauf
neu genommen zu werden, und die Wurzel selbst kommt durch Multiplizieren des
Radius mit diesem Kehrwert zurück statt durch Dividieren; das Dreieck
vergleicht seine drei baryzentrischen Gewichte skaliert mit einer positiven
Zahl, was weder ändert, welches am größten, noch welches negativ ist; der
Teppich zählt in Zellen statt in der Ebene, wo sich die Multiplikation mit
`squares` und die Division dadurch aufheben; und die Schneeflocke nimmt das
Kind, auf das sie übergeht, von der Kante, unter deren Zacke der Punkt gefunden
wurde, die sie ohnehin finden musste. Außerdem geht sie die Koch-Kurve in
`double` ab, da die Figur eine feste Größe hat und 52 Bit eines Bruchteils
einer Kante mehr sind, als die 24 Ebenen, die sie zeichnet, nutzen können. Ein
Test prüft die Verhältnisse, also meldet eine Neufassung, die sie verliert,
das.

## Paletten

Fraktal → Palette fragt nach einer Algorithmusnummer und sagt nach der Nummer,
wie jeder heißt. 1 bis 3 sind die von XaoS; 4 bis 7 wurden neu gebaut, nachdem
gemessen worden war, woraus die ersten drei bestehen, und 8 nach der Vermessung
von Gemälden.

| | |
| --- | --- |
| 1 | **Dunkel und Farbe** — ein Fast-Schwarz mit eigener Tönung und eine Farbe, abwechselnd |
| 2 | **Schwarz, Farbe, Weiß** |
| 3 | **Warm über Dunkelgrün** — leuchtende Violetts, Rots und Oranges über dunklen Grüns und Blaugrüns, jedes dritte Mal Schwarz oder Weiß |
| 4 | **Smog** — dunkel und melancholisch: Graus, Ocker, Schiefer und staubiges Violett, mit Ochsenblut und Scharlach dazwischen |
| 5 | **Warm über Nacht** — Rosé, Rot, Burgunder, Braun, Orange, Ocker und Gelb über Blaugrün, Petrol, Blau und Indigo, im Gerüst von 3 |
| 6 | **Bevorzugte Paare** — die Farbpaare, die die ersten drei am häufigsten nebeneinandersetzen, ganz zwischen Schwarz und Weiß gesetzt |
| 7 | **Zufallsfarben** — jede Stützstelle eine wirklich zufällige Farbe, und sonst nichts |
| 8 | **Kandinsky** — die Verläufe seiner Gemälde, von matt bis leuchtend und hin und wieder fast bis Grau |

### Woraus die ersten drei bestehen

Zwanzigtausend Paletten von jedem, erzeugt, wie das Programm sie erzeugt, jede
Farbe benannt nach der nächsten von drei Dutzend Referenzfarben in CIELAB, die
Farben und die Farbpaare gezählt:

* **Anker.** Schwarz und Weiß sind ein Drittel bis die Hälfte aller
  Stützstellen, und die Farben sitzen mit einer Periode von zwei oder drei
  Stützstellen dazwischen, also zeigt sich das Muster in den ersten paar
  Bändern — und die sind, in einer Palette von etwa dreitausend Stützstellen,
  alles, was die meisten Bilder je benutzen.
* **1 und 2 enthalten 256 Farben und nicht mehr.** Sie nehmen das unterste Byte
  des Generators, und das unterste Byte dieses Generators ist für sich ein
  Generator der Periode 256: Rot bestimmt Grün und Blau, und auf jede Farbe
  folgt immer dieselbe. In zwanzigtausend Paletten von 2 folgte auf keine Farbe
  auch nur einmal eine andere. Daher kommen ihre wiederkehrenden Paare — Rosa
  neben Gelb vierundzwanzigmal so oft wie zufällig, Burgunder neben Oliv
  sechsmal, Rot neben Petrol achtmal.
* **3 schöpft aus zwei Fenstern**: leuchtende Farben von Violett über Rot bis
  Orange und dunkle von Grün bis Blaugrün. Seine häufigsten Paare sind eine
  leuchtende warme Farbe auf einem dunklen Grün.

### Die vier danach

Jede wurde gegen die ersten drei gemessen, bevor sie blieb: die
Jensen-Shannon-Divergenz zwischen den Farben, die zwei Algorithmen zeigen, nach
Namen, und zwischen den Nachbarpaaren, die sie zeigen, über fünftausend
Paletten — 0 für dieselbe Palette, 1 für nichts gemeinsam. 1 und 2 liegen nach
diesem Maß 0.11 auseinander, 1 und 3 0.16.

* **4, Smog**, ist dunkel und verbraucht, wo 1 leuchtend auf Schwarz ist. Tiefe
  Stützstellen sind eine Farbe knapp über Schwarz, helle ein bedecktes Grau mit
  einer Tönung und nie Weiß, und jede vierte ein schlichtes Grau. Ihre Farbtöne
  werden auf einem Ring von 256 abgeschritten, wie 1 und 2 ihre abschreiten,
  also wiederholen sich Nachbarn wie bei ihnen, aber nur über die Farbtöne des
  Smogs — Ocker und Oliv, sowie Stahl, Schiefer und ein staubiges Violett — und
  auf einem Viertel des Rings Karmesin bis Scharlach: Ochsenblut in der Tiefe
  und ein gedämpftes Scharlach im Hellen, gesättigt gehalten, da ein dunkles
  Rot, das es nicht ist, als Braun gelesen wird. Jede fünfte Stützstelle liest
  sich als Rot, gegen eine von dreißig bei 1. Sie liegt 0.43 bis 0.56 von den
  ersten dreien entfernt. Zuerst war sie der Mechanismus von 1 auf einem
  anderen Ring von 256 Farben, so weit vom klassischen Ring gewählt, wie ein
  Ring sein kann, und lag 0.05 von 1 entfernt — näher als 2: die Ähnlichkeit
  war der Mechanismus, die Hälfte der Zeit Schwarz und die andere Hälfte eine
  Farbe von irgendwoher, nicht die Farben.
* **5, Warm über Nacht**, ist das Gerüst von 3 mit umgedrehter Temperatur:
  warme Farben jeder Helligkeit, von Rosé über Rot, Burgunder, Braun und Orange
  bis Ocker und Gelb, über dunkel gehaltenen Nachtfarben, von Blaugrün über
  Petrol und Blau bis Indigo. Beide Fenster sind 125 Grad breit, nahe den 135
  von 3; zuerst waren sie ein Drittel davon und lasen sich als zwei Farben.
  Marineblau ist ihre häufigste Farbe, und Marineblau mit Oliv, Burgunder oder
  Braun ihre häufigsten Paare.
* **6, Bevorzugte Paare**, nimmt jedes Farbpaar, das eine der ersten drei
  mindestens dreimal so oft wie zufällig hintereinandersetzt — sechsundsechzig
  davon, mit den beiden Farben, wie sie tatsächlich herauskamen — und setzt
  jedes Paar ganz, die beiden Farben einander berührend, zwischen eine schwarze
  und eine weiße Stützstelle.
* **7, Zufallsfarben**, ist jede Stützstelle eine zufällig gezogene Farbe aus
  dem oberen Teil des Generators, wo 1 und 2 aus dem unteren ziehen:
  neunundsiebzigtausend verschiedene Farben in zwanzigtausend Paletten, wo 2
  256 hat, und kein Paar mehr als ein Fünftel häufiger als zufällig. Zuerst
  hatte sie jede dritte Stützstelle Schwarz und Weiß und lag 0.00 von 2
  entfernt: der Rhythmus von sehr dunkel und sehr hell war alles, was man sehen
  konnte. Ohne ihn hat sie weniger Spanne von dunkel nach hell als die
  verankerten, und genau das heißt zufällig.

Die vier, die sie ersetzen, waren ein Spektrum, ein Duoton, eine Triade und ein
Komplementärpaar, aus dem Farbkreis hergeleitet, und auf dem Bildschirm war
jede eine einzige Farbe. Sie verteilten ihre Farbtöne über die ganze Palette,
und die Palette ist im Programm etwa dreitausend Stützstellen lang — vier bis
neun in dem Test, der es hätte merken sollen, weshalb er es nicht merkte —,
also sah ein Bild nur den ersten Farbton. Und sie hielten jede Stützstelle als
Farbe, tief oder blass, ohne Schwarz und ohne Weiß, also hatte nichts an ihnen
eine Kante. Der Test erzeugt seine Paletten jetzt, wie das Programm es tut, und
verlangt von allen sieben, dass der Teil, den ein Bild sieht, von dunkel nach
hell geht und mehr als eine Farbe enthält; von 8, die danach kam, wird etwas
weniger verlangt, und warum, steht weiter unten.

Eine Position, die vor dieser Änderung mit 4 bis 7 gespeichert wurde, kommt in
anderen Farben zurück: Nummer und Anfangswert sind alles, was eine Position von
ihrer Palette festhält, und was die Nummer bedeutet, hat sich geändert. Eine
mit 1 bis 3 gespeicherte ist nicht betroffen, bis aufs Bit, und eine mit 4 bis
7 gespeicherte nennt einen Algorithmus, den das ursprüngliche XaoS nicht hat
und ablehnen wird.

### 8, Kandinsky

Gemacht aus zwei Bildreihen in seiner Art: vier leuchtende Vignetten — Schwarz,
Zinnoberrot, Gelb, Coelinblau, Violett und Rosa auf cremefarbenem Papier — und
ein mattes Bild, gebranntes Orange, Ocker, Braun, Salbei und Petrol auf
demselben Papier. Die Vignetten ergeben achtzehn Farben, per k-Means in CIELAB,
und jede hat eine **matte Form**: die Farbe des matten Bildes, die ihr in
Helligkeit und Farbton am nächsten liegt, die Buntheit außen vor, da die beiden
sich gerade in der Buntheit unterscheiden. Zinnober wird gebranntes Orange,
Gelb Ocker, Burgunder Braun, Schwarz ein dunkles Braun, Coelinblau Salbei, Blau
Petrol. Das matte Bild hat kein Grau, also bleibt Grau Grau. Das Aquamarin, das
es zu haben scheint, ist ein Salbei, 128,146,121: als Aquamarin liest es sich
nur neben dem Orange.

Jede Palette zieht eine **Stimmung**. Vier von fünf fallen gleichmäßig zwischen
die matten Formen und die leuchtenden; die fünfte geht über die matten Formen
hinaus Richtung Grau, bis auf fünfzehn Hundertstel ihrer Buntheit, und etwa
eine Palette von zwanzig ist fast grau.

**Anker kommen zufällig.** Eine Palette von zehn darf hier und da eine
Stützstelle zu reinem Schwarz oder reinem Weiß machen: jede Stützstelle nach
der ersten mit eigenem Wurf, mit einer Rate, die die Palette zwischen einer
Stützstelle von zwanzig und einer von zehn zieht, Schwarz oder Weiß, wie der
Wurf fällt — kein Rhythmus und kein Abwechseln, anders als die Anker von 1
bis 3. Die anderen neun von zehn Paletten sind genau, was sie waren, bevor die
Anker hinzukamen: der Wurf wurde schon für jede Stützstelle getan und blieb
ungenutzt.

Was auf eine Farbe folgt, ist, wohin **seine Verläufe** sie führen. Überall, wo
sich die Farbe der Vignetten sanft ändert — über acht Pixel, ohne einen Schritt
größer als ein Drittel der Änderung —, wurden die Farben an den beiden Enden
gezählt, und jede Stützstelle zieht die nächste aus diesen Zählungen: Papier
geht über in eine Lasur von Gelb und dann Gelb, Schwarz in Zinnober (drei von
zehn Malen) und Zinnober in Burgunder, Blau in Mauve und Mauve in Rosa,
Coelinblau in ein graues Aquamarin. Die harten Kanten wurden ebenfalls gezählt
und weggelassen: die meisten sind Papier neben Papier einer anderen
Schattierung, und daraus gemachte Paletten kamen als Papier und Grau heraus.
Eine Palette beginnt auf dem Papier oder auf dem Schwarz.

Über zwanzigtausend Paletten liegen die des obersten Zehntels der Stimmungen
0.06 von den Vignetten entfernt, nach derselben Divergenz wie oben, und die des
Zehntels über null 0.07 vom matten Bild; von den anderen sieben Paletten liegen
sie 0.39 bis 0.67 entfernt.

Weil ein Verlauf eine Weile in einem Teil des Farbkreises bleibt und 8 kein
Schwarz und Weiß alle paar Stützstellen hat, auf das sie zurückfallen kann,
verlangen die Tests weniger von ihr als von den anderen. Wo die Palette lang
ausgelegt ist, vierzig Einträge pro Stützstelle und mehr, sieht ein Bild etwa
zehn Stützstellen, und eine Palette von hundert ist dann ein einziges Feld —
Ocker mit Salbei darin, Petrol mit Zinnober —, das weder ein Drittel des Weges
von Schwarz nach Weiß noch eine zweite Farbe erreicht. Und drei Stützstellen,
die eine kurze Palette von 256 Einträgen oft ist, sind ein einziger Verlauf.

## Den Orbit beobachten

`trap` und `stripe` halten eine Zahl über einen ganzen Orbit fest statt über
den Punkt, und das ist die Art Größe, die die Färbungsmodi der Engine nicht
erreichen können: eine Berechnungsschleife ist dort so kompiliert, dass sie
eine Farbe zurückgibt, und nur das letzte `z` und das davor überleben, um damit
gefärbt zu werden. In der Formel getan, kostet es die Engine nichts.

Beide geben den erhaltenen Wert zurück bis zum letzten Durchlauf, den das
Iterationslimit erlaubt, und in diesem, was sie gesammelt haben. Also ist eine
ganze Formel

    trap(z^2+c, 3)

und sonst nichts: das Fraktal iteriert, wie es würde, und im letzten Durchlauf
wird der Wert zur Falle, die dann die inneren Färbungsmodi zeichnen. Ein
austretender Punkt erreicht diesen Durchlauf nie und behält seine gewöhnliche
Außenfarbe, also zeichnen diese das Innere.

**`trap(orbit, shape, centre, size)`** — das Nächste, was der Orbit je an eine
Form herankam. `shape` ist standardmäßig 0, `centre` der Ursprung, `size` 1.

| shape | |
| --- | --- |
| `0` | der Mittelpunkt selbst, ein Punkt |
| `1` | eine waagerechte Gerade durch ihn |
| `2` | eine senkrechte |
| `3` | beide, ein Kreuz |
| `4` | ein Ring mit Radius `size` |
| `5` | ein Quadrat mit halber Seitenlänge `size` |
| `6` | eine Raute mit halber Diagonale `size` |

**`stripe(orbit, density)`** — der Mittelwert von
`(sin(density * arg orbit) + 1) / 2` entlang des Orbits. `density` ist
standardmäßig 4 und gibt an, wie viele Streifen eine Umdrehung ausmachen; eine
ganze Zahl, sonst treffen sie sich nicht, wo sich die Umdrehung schließt. Ein
Mittelwert über den Orbit ändert sich mit dem Punkt glatt, auch wo die
Iterationszahl springt, und das zeichnet die Fasern, für die das Verfahren
bekannt ist.

    trap(z^2+c,0)                    wie nah der Orbit am Ursprung vorbeikam
    trap(z^2+c,3)                    ein Kreuz, das Strahlen zeichnet
    trap(z^2+c,4,{0,0},{0.5,0})      Ringe in der Kardioide
    stripe(z^2+c,6)                  sechs Streifen pro Umdrehung

Die laufende Größe lebt an der Aufrufstelle, also behalten zwei Fallen in einer
Formel jede ihre eigene, und ein Thread kann einen anderen nicht stören. Beide
schalten das Boundary-Trace aus, aus demselben Grund wie das Rauschen: zwei
Nachbarn, die gleich viele Durchläufe brauchen, können dennoch ganz
verschiedene Orbits gesehen haben, und eine Region, die gefüllt wird, ohne
berechnet zu werden, wäre falsch.

Dies sind keine Färbungsmodi und können es nicht sein: die Modi leben in der
kompilierten Schleife jeder Formel, die direkt eine Farbe zurückgibt, und davon
gibt es einunddreißig mit je einer eigenen Symmetrietabelle nach Position. Was
hier steht, funktioniert nur für eine Benutzerformel.

## Form des Fluchtradius

Berechnung → Fluchtradius → Form des Fluchtradius. Die Iteration endet, wenn
`z` eine Region verlässt, und die war immer ein Kreis. Die Form ist in den
Bändern außerhalb der Menge zu sehen, also ist sie ebenso ein Zeichenwerkzeug
wie ein numerisches.

| Form | tritt aus, wenn |
| --- | --- |
| Kreis (klassisch) | \|re\|² + \|im\|² den Fluchtradius übersteigt |
| Quadrat | eine der beiden Komponenten allein es tut |
| Raute | \|re\| + \|im\| es tut |
| Reelle Achse | die reelle Komponente allein |
| Imaginäre Achse | die imaginäre Komponente allein |
| Beide Achsen | beide Komponenten zugleich |
| Dreieck 0° / 90° / 180° / −90° | der Punkt ein Dreieck dieser Ausrichtung verlässt |
| Sechseck 0° / 90° | der Punkt ein Sechseck dieser Ausrichtung verlässt |
| Achteck | der Punkt ein Achteck verlässt |

Die Vielecke werden an ihrem **Inkreisradius** gemessen — dem Abstand von der
Mitte zur Mitte einer Seite —, also berühren sich ein Vieleck und der Kreis
desselben Fluchtradius entlang der Seiten statt an den Ecken. Die Ausrichtung
gehört zur Form: ein um 120° gedrehtes Dreieck ist dasselbe Dreieck, also sind
0°, 90°, 180° und −90° vier verschiedene, während sich ein Sechseck alle 60°
wiederholt und daher zwei hat, und ein Achteck alle 45°, was ihm nichts zu
wählen lässt. Dreieck und Sechseck haben deshalb ein eigenes Untermenü für die
Ausrichtung. Ihre Normalen werden einmal pro Bild berechnet, also kostet die
innere Schleife drei bis acht Multiplikations-Additionen statt eines Sinus und
eines Kosinus pro Seite.

Gegenüber dem Kreis ändert die Raute etwa 39 % der Pixel und die reelle Achse
44 %. Die Silhouette der Menge bewegt sich nicht — die entscheidet, was nie
austritt —, also schauen Sie auf die Bänder.

In eine gespeicherte Position wird sie nur geschrieben, wenn sie nicht der
Kreis ist, also bleibt eine Position, die keine Form verlangt, mit jeder
früheren Version ladbar.

## Fraktionale Brownsche Bewegung

Ein Dialog auf jeder Seite — **Fraktal → Innere Färbung → Fraktionale Brownsche
Bewegung...** und dasselbe unter Äußere Färbung —, mit den drei Modi, die er
bietet, und den fünf Zahlen, aus denen sie gezeichnet werden. Es war ein
Untermenü, die drei Modi als Auswahleinträge und die Zahlen hinter einem
Einstellungen-Dialog darunter; jetzt ist es ein einziger Dialog, angehakt, wann
immer einer der drei Modi in Gebrauch ist, und OK setzt den Modus und die
Zahlen zusammen als einen Rückgängig-Schritt. Eine Position speichert den Modus
weiterhin als `(incoloring N)` oder `(outcoloring N)` und die Zahlen getrennt,
als `(infbmset ...)` oder `(outfbmset ...)`, sodass eine Datei aus einer
früheren Version unverändert öffnet.

Sie lässt eine saubere Färbung **gebraucht** aussehen: die Bänder behalten Form
und Reihenfolge und hören auf, perfekt zu sein, als hätte man das Ding draußen
im Wetter stehen lassen. Es ist kein Zufall Pixel für Pixel, der sich als Korn
liest und versucht und verworfen wurde; es ist ein Rauschfeld mit einem
Maßstab, also zeichnet es Spuren.

**Außen** sind die drei Modi `fbm + smooth`, `fbm + iter` und `fbm` allein —
die weiche Zählung abgenutzt, die gebänderte Zählung abgenutzt, und das
Rauschen für sich, das Wolken um die Menge zeichnet, ganz ohne Zählung darin.

**Innen** sind es `fbm`, `fbm + zmag` und `fbm + decomposition`. Der erste ist
der interessante: das Innere ist standardmäßig ein einziger flacher Ton, und
dies gibt ihm eine Oberfläche.

### Gradientenrauschen, und kein Gitter

Die Bewegung ist **Gradientenrauschen**: jede Ecke eines Gitters bekommt eine
Richtung, und ein Punkt nimmt von jeder Ecke die Steigung, die diese Richtung
in seinem Abstand ergibt. Früher war es Wertrauschen — eine Zahl an jeder Ecke,
mit einem Smoothstep überblendet —, und das zeichnete ein Gitter. Ein
Smoothstep hat an keinem Ende eine Steigung, also wurde das Feld entlang jeder
Gitterlinie flach, und da das Gitter jeder Oktave die Linien der ersten
enthielt, wurden sie dort alle zugleich flach. Quer über eine Linie war die
Steigung null; das Bild kam in Quadraten heraus, acht Pixel Seitenlänge in
`FBM_ERROR.xpf`, als wäre es in geringerer Auflösung berechnet worden, als es
gezeigt wurde. Gradientenrauschen ist auf den Linien so steil wie zwischen
ihnen — gemessen 0.99 davon bei den Standardeinstellungen und 1.08 bei denen
dieser Position —, und das Gitter jeder Oktave ist gegen das vorige verschoben,
um Beträge, die keine Verdoppelung je ganz macht, also teilen keine zwei eine
Linie oder eine Ecke. `fbm()` in Formeln zeichnet dieselbe Bewegung, und
`randsc` ist eine Oktave davon.

Was blieb: der Wert läuft immer noch von null bis zur Intensität und nie
darüber hinaus; die **Spuren sind so groß wie vorher**, wobei das Gitter mit
drei Fünfteln der Frequenz läuft, damit Gradientenrauschen sie so groß zeichnet
wie Wertrauschen (gemessen daran, wie weit sich das Feld bewegt, bevor es sich
selbst nur noch halb so sehr ähnelt: 0.56 einer Zelle bei den
Standardeinstellungen, vorher und nachher); der **Kontrast**, von dem
Gradientenrauschen allein nur sieben Zehntel hat — jede Oktave geht durch eine
Kurve, die ihre Werte auf die Verteilung bringt, die Wertrauschen hatte,
monoton, damit sie keinen Umriss verschiebt, und gemessen über die Bewegung bei
den Standardeinstellungen schwingt sie 0.133 um ihre Mitte gegen 0.132 vorher;
und die Kosten, ungefähr — die Modi etwas schneller als vorher, `fbm()` etwa
ein Achtel langsamer. Was nicht blieb: eine mit einem fbm-Modus gespeicherte
Position behält ihre Zahlen, ihren Maßstab und ihre Stärke und **zeichnet
andere Spuren**.

### Die fünf Zahlen

Derselbe Dialog wie der Modus, darunter. Sie werden mit der Position
gespeichert.

| | |
| --- | --- |
| **Intensität** | um wie viele Farbbänder das Rauschen den Wert verschiebt. `4` ist der Standard; unter eins ist es ein Flüstern, über zehn hören die Bänder auf, Bänder zu sein |
| **Frequenz** | die Größe der Spuren, so groß, wie ein Gitter mit so vielen Zellen pro Einheit der Ebene sie zeichnen würde. **Das ist die, die man beim Hineinzoomen erhöht**: in der ersten Ansicht liest sich `8` gut, bei einer Spanne von `0.02` braucht es etwa `200` |
| **Oktaven** | wie viele summiert werden, jede mit der doppelten Frequenz der vorigen |
| **Rauheit** | was jede Oktave von der Höhe der vorigen behält. `0.5` ist die einfache Bewegung; höher ist körniger |
| **Anfangswert** | jedes Mal dasselbe Bild. Ändern Sie ihn für ein anderes vom selben Charakter |

**Oktaven und Rauheit stützen einander.** Bei einer Rauheit von `0.5` trägt die
achte Oktave ein Zweihundertfünfzigstel des Ganzen, also ändern Oktaven über
vier oder fünf hinaus nichts, was man sehen kann — gemessen ergaben vier und
achtzehn dasselbe Bild. Erhöhen Sie die Rauheit, und die feinen Oktaven haben
etwas auszugeben, und dann fängt die Anzahl an zu zählen.

### Die Größe der Spuren nach der Zählung

Außen zwei Felder mehr am Fuß des Dialogs. **Skalierung des Rauschens** sagt,
was die Größe der Spuren bestimmt: **fest**, wo es beginnt, die Frequenz
allein, wie oben; **Iteration**, die Frequenz, mit der Iterationszahl erhöht;
**geglättete Iteration**, die Frequenz, mit der weichen Zählung erhöht.
**Faktor** sagt, um wie viel: die Frequenz verdoppelt sich so oft je zehn
Iterationen, also sind die Spuren bei `1` halb so groß, wo der Orbit zehn
Iterationen länger zum Austreten brauchte, und ein Viertel, wo er zwanzig
länger brauchte. Bei `0` ändert die Zählung nichts. Wo die Zählung steigt, zum
Rand der Menge hin, werden die Spuren feiner und folgen dem Detail.

Eine Frequenz, die sich von einem Pixel zum nächsten ändert, lässt sich nicht
einfach über die Ebene legen: das Feld würde sich dehnen, umso mehr, je weiter
vom Ursprung. Also wird die Bewegung bei den beiden Zweierpotenzen zu beiden
Seiten der Zählung gelesen, jede ein Feld für sich mit fester Frequenz, und die
beiden werden überblendet, wie eine Textur zwischen ihren Mipmaps, und die
Mischung auf den Kontrast eines einzigen Feldes zurückgebracht. Mit
**geglättete Iteration** ändern die Spuren ihre Größe ohne Naht; mit
**Iteration** ändern sie sie an den Rändern der Bänder, jedes Band mit einer
eigenen Körnung. Innen ist die Zählung überall dieselbe, also hat die
Innenseite keines der beiden Felder.

Zwischen zwei Verdopplungen werden zwei Felder gelesen, und die Färbung kostet
etwa doppelt so viel wie mit **fest**: 248 ns pro Pixel gegen 127, bei vier
Oktaven. Wo die Zählung genau auf eine Verdopplung fällt, wie es jede ganze
Zählung bei einem Faktor von `10` tut, wird ein einziges Feld gelesen, und die
Kosten sind die von **fest**. Mit **geglättete Iteration** werden `fbm + iter`
und `fbm` von der Smooth-Schleife gezeichnet, die einige Iterationen über den
Fluchtradius hinaus rechnet, um die weiche Zählung zu bestimmen. Ein zu hoher
Faktor bringt die Spuren unter ein Pixel, wo sie zu Körnung werden.

Mit **fest** zeichnet es wie zuvor, Bit für Bit. Eine Position speichert die
beiden als `(outfbmscale 'smooth 1.5)`, und nur, wenn sie verändert wurden,
also öffnet sich eine Position, die sie nicht anrührt, weiterhin in früheren
Versionen.

### Was es kostet, und das Einzige, was es nicht tut

Das Rauschen wird **einmal pro Pixel** berechnet, in der Färbung, nicht einmal
pro Iteration. Über ein Bild von 90000 Punkten bei einer Spanne von `0.02` und
400 Iterationen brauchte `fbm + smooth` mit vier Oktaven ein Fünftel länger als
schlichtes `smooth` und mit achtzehn ein Viertel länger; bei 4000 Iterationen,
wo der Orbit mehr von der Arbeit tut, kosteten dieselben vier Oktaven ein
Zehntel.

Die Spuren sind **an der Ebene** verankert, nicht an der Ansicht: es ist
Schmutz auf dem Ding, nicht auf der Linse. Also vergrößert das Zoomen sie, und
wenn man weit genug geht, landet man in einer einzigen Spur, und das Bild ist
wieder sauber. Die Frequenz zu erhöhen, ist die Antwort, und deshalb ist die
Frequenz eine Einstellung und keine Konstante. Stattdessen an den Zoom binden
lässt sie sich nicht: die Engine verwendet Pixel wieder, die im vorigen Maßstab
berechnet wurden, und eine Färbung, die sich mit der Ansicht bewegte, ließe sie
die falschen Spuren tragen.

Die Modi sind **nach** den Echtfarbmodi nummeriert, also bedeutet jede
Modusnummer, die eine gespeicherte Position schon enthält, weiterhin, was sie
bedeutete.

## Kachelung

**Fraktal → Innere Färbung → Kachelung...** färbt das Innere nach den Kacheln
einer der fünfundvierzig Kachelungen von `randsctile`, die über dem letzten `z`
liegen, vierzig pro Einheit wie die Quadrate von `squares`. Ihr Dialog hat drei
Felder.

| Feld | |
| --- | --- |
| **Kachelung** | welche Kachelung, mit der Zahl, unter der `randsctile` sie kennt, und ihrem Namen: **3: Sechsecke**, wo es beginnt. Der Reiter Kachelungen von Hilfe → Referenz der Benutzerformeln hat ein Bild jeder |
| **Färbung** | woher die Farbe kommt, überall: **Winkel von z**, wo es beginnt; **Winkel um die Kachelmitte**, der um die Mitte der Kachel, in der der Punkt liegt, der jede Kachel zu einem kleinen Farbrad macht; **zmag**, der Wert des Modus dieses Namens; **keine**, nur die Kacheln |
| **Kacheln** | wie sich jede Kachel von ihren Nachbarn abhebt: **jede um einen eigenen Betrag verschoben**, entlang der Palette, wo es beginnt; **von der Mitte zum Rand schattiert**, dasselbe, und von der Mitte zum Rand um ein Achtel der Farben einer Drehung weiter verschoben, wie eine Fase; **mit betontem Rand**, dasselbe, und der äußere Teil der Kachel in der anderen Hälfte der Palette gezeichnet; **jede in einer eigenen Farbe**, eine Farbe für sich, woher die Farbe auch kommt |

Der Modus steht zwischen den Modi, die sich wiederholen, und den Zerlegungen,
und die Einstellungen beider erreichen ihn. **Frequenz (periodische Modi)** und
**Funktion (periodische Modi)** erreichen die Kacheln, wie sie jedes
periodische Argument des Inneren erreichen: bei Frequenz 2 sind die Kacheln
halb so groß. Die [Einstellungen der
Zerlegungen](#einstellungen-der-zerlegungen) des Inneren erreichen die beiden
Winkel, die geschnitten und durchlaufen werden wie die der Zerlegungen.
**Sektoren**, **Drehung**, **Wiederholungen**, **Richtung** und
**Geschwindigkeit** erreichen beide; **Winkel von** erreicht nur den Winkel von
`z`, da der um die Kachelmitte nichts damit anfangen kann; **Ringe pro
Einheit** macht beim Winkel von `z` Ringe des Betrags, durch Frequenz und
Funktion wie bei den Zerlegungen, und beim Winkel um die Mitte ebenso viele
Ringe von der Mitte zum Rand, die der Form der Kachel folgen. **zmag**,
**keine** und eine Kachel in eigener Farbe lesen keine davon.

Der eigene Betrag einer Kachel kommt aus dem Hash, der sie benennt. Bei einem
Winkel umfasst er so viel der Palette wie eine Drehung des Winkels, sodass eine
verschobene Einstellung der Zerlegung die Farben der Kacheln lässt, wo sie
waren; bei **zmag**, **keine** und einer Kachel in eigener Farbe umfasst er die
ganze Palette, sodass **keine** mit **jede um einen eigenen Betrag verschoben**
dasselbe ist wie **jede in einer eigenen Farbe**. **Winkel um die Kachelmitte**
legt eine ganze Drehung der Farben in jede Kachel, was bei der Größe der
Kacheln sehr viele Farben auf wenige Pixel sind: eine **Geschwindigkeit** unter
eins oder eine Frequenz unter eins macht die Räder lesbar.

Mit **Färbung** auf **Winkel von z**, **Kacheln** auf **jede um einen eigenen
Betrag verschoben** und den Einstellungen der Zerlegungen, wo sie beginnen,
zeichnet der Modus, was er zeichnete, bevor er eine Wahl hatte, bitgenau. Die
Kachelungen, die sich nie wiederholen, werden jenseits von viertausend
Millionen Kacheln vom Ursprung flach, ohne eigenen Betrag.

Die periodischen Kachelungen kosten etwa so viel wie die anderen inneren Modi,
120 bis 210 ns pro Pixel; die sich nie wiederholenden, die an jedem Punkt neu
gesucht werden, 250 bis 490. **Winkel um die Kachelmitte**, **von der Mitte zum
Rand schattiert** und **mit betontem Rand** messen, wo der Punkt in seiner
Kachel steht, was bei den periodischen Kachelungen bis zu zwei Drittel mehr
kostet und bei den anderen wenig mehr. Eine Position speichert den Modus als
`(incoloring 30)`, die Kachelung als `(intessellation 39)`, nur wenn sie nicht
die der Sechsecke ist, und die Färbung als
`(intessellstyle 'tile-angle 'edged)`, nur wenn sie verschoben wurde.

## Weitere Färbungsmodi

Nicht die Berechnung entscheidet über die Farbe eines Pixels: die Farbe kommt
aus einer Formel im letzten `z`, im Parameter `c` und in der Iterationszahl,
und es gibt eine solche Formel pro Modus. Zweiundzwanzig weitere, in einem
eigenen Untermenü unter Fraktal → Innere Färbung und Fraktal → Äußere Färbung,
über den Echtfarb-Untermenüs.

**Innen**, wo der Orbit nie austrat und nur das letzte `z` und `c` etwas sagen:

| Modus | |
| --- | --- |
| `real/mag` | der Kosinus des Winkels: eine Schattierung, die sich dort in sich schließt, wo `atan2` eine Naht hat |
| `max(\|real\|,\|imag\|)` | die Quadratnorm, also sind die Bänder Quadrate, wo die von `zmag` Kreise sind |
| `\|real\|+\|imag\|` | dieselben Bänder, um ein Achtel einer Umdrehung gedreht |
| `min(\|real\|,\|imag\|)` | nahe null entlang beider Achsen, also Strahlen |
| `\|z-c\|` | wie weit vom Parameter sich der Orbit gesetzt hat, was die Knospen unterscheidet |
| `\|z*c\|` | die beiden Beträge gegeneinander |
| `angle(z)-angle(c)` | die Richtung, in die er sich setzte, gemessen am Parameter |
| `real*imag` | ein Sattel, vier Lappen um die Achsen |
| `sin(real)*sin(imag)` | ein Gitter über der Stelle, an der sich der Orbit setzte |
| `sign(imag)` | zwei flache Töne: eine Zerlegung des Inneren, wie die binäre Zerlegung sie dem Äußeren antut |
| `frac(mag)` | Höhenlinien bei jedem Achtel des Betrags |
| `log(mag)` | `zmag` mit dem Kontrast auf die kleinen Werte verlegt, wo ein angezogener Orbit seine Zeit verbringt |
| `cos(\|z\|)` | Ringe um den Punkt, an dem sich der Orbit niedergelassen hat, acht pro Einheit des Betrags |
| `hexagons` | drei Wellen, eine Dritteldrehung voneinander, die zusammen ein Sechseckgitter ergeben |
| `spiral` | fünf Arme um den Punkt, an dem sich der Orbit niedergelassen hat, einmal gewunden pro Faktor e im Betrag |
| `cos(\|z-c\|)` | Ringe um den Parameter, die die Knospen unterscheiden, wie es `\|z-c\|` tut |

Elf der inneren Modi wiederholen sich — `cos(mag)`, `mag*cos(real^2)`,
`sin(real^2-imag^2)`, `squares`, `sin(real)*sin(imag)`, `frac(mag)`, die vier
am Fuß der Tabelle und `tessellation` —, und **Frequenz (periodische Modi)**,
das erste Feld des inneren Dialogs Färbungsgeschwindigkeit, sagt, wie oft. Sie
multipliziert das Argument dessen, was sich wiederholt, also drehen sich bei 2
die Sinusse doppelt so schnell, `squares` zeichnet halb so große Quadrate und
`frac(mag)` sechzehn Konturen pro Einheit des Betrags, wo es acht zeichnete. In
`cos(mag)` verengt sie auch die Bänder, die zwischen Sinus und Kosinus wählen.
Bei eins, wo sie anfängt, zeichnen sie wie immer, bitgenau, und die anderen
Modi lesen sie nie. Eine Position speichert sie als `(incolorfrequency 2)` nur,
wenn sie nicht eins ist, also öffnet sich eine Position, die sie nicht anrührt,
weiterhin in früheren Versionen.

**Funktion (periodische Modi)**, darunter, biegt das Argument mit einer der
[Geschwindigkeitsfunktionen der
Färbung](#geschwindigkeitsfunktionen-der-färbung), bevor die Frequenz es
multipliziert, sodass die Bänder schneller umlaufen, wo die Funktion steil
ansteigt, und langsamer, wo sie flach wird: `square` und `cube` drängen sie zum
Rand hin und ziehen sie um null auseinander, `log` und `atan` umgekehrt, und
`sin`, `cos` und `triangle` falten das Argument auf sich selbst zurück, sodass
die Bänder umkehren und zurückkommen. Sie erreicht jede periodische Funktion
der elf, deren Argument `a` als `k·f(a)` eingeht, wobei `k` die Frequenz ist:
die Sinus und Kosinus, den Betrag von `frac(mag)`, Real- und Imaginärteil,
deren Streifen die Quadrate von `squares` und die Kacheln von `tessellation`
bilden, und in `cos(mag)` auch die Bänder, die zwischen Sinus und Kosinus
wählen. In `spiral` bleiben die fünf Arme fünf, denn eine nicht ganze Zahl von
Umläufen hinterließe eine Naht, wo sich der Winkel schließt: Frequenz und
Funktion wirken auf die Windung. Eine periodische Funktion im Argument einer
anderen bekäme ihr eigenes Argument ebenso gebogen. Bei `ident`, wo sie
beginnt, ändert sich nichts. Eine Position speichert sie als
`(incolorperiodicfun 6)`, mit dem Platz der Funktion in der Liste von null an
gezählt, nur wenn sie nicht ident ist.

**Außen**, wo die Zählung und der Austrittspunkt das sind, was es gibt:

| Modus | |
| --- | --- |
| `iter+angle` | die Richtung, in der er austrat, in die Zählung gemischt: die Bänder bekommen einen Drall |
| `iter+log(mag)` | wie weit er über den Fluchtradius hinausging, eine billigere Glättung als der Logarithmus eines Logarithmus |
| `iter+real*imag` | Biomorphe, stetig gemacht statt mit Schwelle |
| `max(\|real\|,\|imag\|)` | auf welcher Seite er austrat, ganz ohne Zählung |
| `iter banded` | die Zählung in acht Bänder gefaltet, die die Form der Niveaumengen zeigt statt ihrer Nummer |
| `\|real\|-\|imag\|` | wie schief der Austrittspunkt liegt |

Die Nummern der älteren Modi verschieben sich nicht, also lädt eine
gespeicherte Position weiter wie bisher.

Färbungsgeschwindigkeit und Färbungsverschiebung erreichen jetzt smooth und
smooth log, was sie nie taten: diese beiden gaben ein eigenes Pixel zurück,
statt durch den Schritt zu gehen, mit dem die anderen Modi enden, also taten
beide Regler gar nichts, solange einer der beiden gewählt war. Bei den Werten,
mit denen sie kommen — Geschwindigkeit eins, Verschiebung null —, ist der
Schritt die Identität, also verschiebt sich kein zuvor gespeichertes Bild.

Smooth und smooth log folgen der Form des Fluchtradius. Sie interpolieren
zwischen einem Durchlauf und dem nächsten, indem sie fragen, wo der Orbit den
Fluchtradius kreuzte, und früher fragten sie das beim Betrag, welche Form auch
gewählt war — also brachte Austreten auf einem Quadrat und Interpolieren auf
einem Kreis mehr als ein Fünftel der austretenden Punkte, bei einer Raute die
Hälfte, außerhalb des Bandes, zu dem sie gehören. Jede Form übergibt jetzt die
Größe, die sie tatsächlich prüft, und die Schwelle, gegen die sie sie prüft.
Für den Kreis ist beides dasselbe, also verschiebt sich nichts, was zuvor
gezeichnet wurde.

Sie funktionieren auch im Newton-Modus, wo ein Orbit konvergiert statt
auszutreten: die Größe ist dann, wie weit sich der letzte Schritt bewegte, und
die Schwelle die Konvergenzgrenze, und die Interpolation ist dieselbe Frage, an
eine fallende statt an eine steigende Zahl gestellt. Die beiden
Newton-Fraktale, die gar keine weiche Variante hatten, haben eine.

**Newton-Konvergenz** im Menü Berechnung wird endlich gelesen. Die Grenze war
als ein Millionstel in die Iteration geschrieben, also wurde die Einstellung
angeboten, gespeichert und neu geladen, während nichts sie je befragte.

**Smooth und smooth log** sind in dieses Untermenü im äußeren Menü gezogen. Sie
interpolieren zwischen einer Iteration und der nächsten anhand dessen, wie weit
der Orbit über den Fluchtradius hinausging, was eine Formel braucht, die am
Fluchtradius austritt — zwanzig der einunddreißig tun das, und die anderen elf
konvergieren entweder (die Newtons) oder halten bei einem eigenen Test an (die
Sierpinskis, Koch, die Clock). Für die tut der Modus gar nichts, und das sollte
die erste Liste nicht anbieten. Die Benutzerformel, die ihn immer hätte haben
können und nicht hatte, hat ihn jetzt.

Was diese Modi nicht können, ist alles, was den Orbit statt seines Endes
braucht: eine Orbitfalle, ein Streifenmittel, ein Ljapunow-Exponent, die
Abstandsschätzung. Die brauchen einen Wert, der von einer Iteration zur
nächsten mitgetragen wird, und das heißt in dieser Engine eine zweite
kompilierte Kopie jeder Berechnungsschleife. `trap` und `stripe` in einer
Benutzerformel tun es stattdessen.

### Einstellungen der Zerlegungen

**Fraktal → Innere Färbung → Zerlegung**, und dasselbe unter Äußere Färbung,
legt fest, wie die Zerlegungen die Drehung schneiden: binary decomposition und
color decomposition außen, `sign(imag)`, `Decomposition-like` und
`fbm + decomposition` innen, und die Winkel von `tessellation`. Mit jedem Feld,
wo es beginnt — null, und eins, rechts und eins bei den letzten drei —,
zeichnen sie wie immer, bitgenau.

| Feld | |
| --- | --- |
| **Sektoren** | in wie viele jede Wiederholung der Drehung geschnitten wird. Null ist die des Modus: die zwei Halbebenen der binären, die ungeteilte Drehung der übrigen. Die binären färben abwechselnde Sektoren, also wird eine ungerade Zahl zur nächsten geraden aufgerundet — zwei Sektoren einer Farbe würden sich berühren. Die farbigen werden zu ebenso vielen flachen Stufen |
| **Drehung** | um wie viele Grad die Sektoren gedreht werden, gegen den Uhrzeigersinn |
| **Winkel von** | nur innen: der Winkel von `z`, von `z-c` oder von `z*c`. `z-c` zentriert die Sektoren auf den Parameter, `z*c` dreht sie um den Winkel des Parameters |
| **Ringe pro Einheit** (innen) | Ringe des Betrags dessen, dessen Winkel genommen wird, abwechselnd unterschieden: die binären werden zu einem polaren Schachbrett. Sie gehen durch `k·f(a)` wie jedes periodische Argument des Inneren, also erreichen sie Frequenz und Funktion der periodischen Modi |
| **Iterationen pro Band** (außen) | Bänder aus so vielen Iterationen quer über die Sektoren: mit der binary decomposition das Gitter aus äußeren Strahlen und Äquipotentiallinien |
| **Wiederholungen** | wie oft sich die Drehung wiederholt: die Sektoren, und bei den farbigen die Palette, beginnen in jeder von vorn. Jede Zahl, doch ein Bruchteil hinterlässt eine Naht, wo sich der Winkel schließt |
| **Richtung** | rechts, wie sie immer liefen; links, andersherum; beidseitig, in jeder Wiederholung vor und wieder zurück — ein Spiegel, ohne Naht, wenn die Wiederholungen ganz sind |
| **Geschwindigkeit** | bei den farbigen, wie viel der Palette eine Wiederholung überspannt: über eins mehr Farben in jeder, nicht mehr Wiederholungen. Bei den binären der Abstand zwischen den Farben abwechselnder Sektoren: bei eins wie immer, bei null eine einzige Farbe, bei zwei doppelt so weit — binary decomposition rückt die ungeraden Sektoren von der Zählung zur rückwärts gezählten, `sign(imag)` seinen zweiten Ton vom ersten |

Eine Position speichert sie als `(indecomposition 8 22.5 'z-c 6 2.5 'both 1.5)`
und `(outdecomposition 6 30 2 3 'left 2)` nur, wenn sie verändert sind, also
öffnet sich eine Position, die sie nicht anrührt, weiterhin in 2.5. Verändert
ist eine Zerlegung nicht mehr symmetrisch zur reellen Achse, und XaoS hört auf,
eine Hälfte als Spiegelbild der anderen zu zeichnen.

## Geschwindigkeitsfunktionen der Färbung

**Fraktal → Innere Färbung → Färbungsgeschwindigkeit**, und dasselbe unter
Äußere Färbung, ist ein einziger Dialog für alles, was den Wert, den ein
Färbungsmodus übergibt, an einen Platz in der Palette bringt. Die **Funktion**,
aus einer Liste gewählt, biegt den Wert, bevor **Färbungsgeschwindigkeit** und
**-verschiebung** angewandt werden: die Palette wird bei `f(v) × speed + shift`
gelesen. In einer eigenen Spalte daneben steht der [Geschwindigkeitswechsel
nach ITER](#geschwindigkeitswechsel-nach-iter), der auf den Wert wirkt, bevor
die Funktion ihn bekommt. Der innere beginnt mit der Frequenz und der Funktion
der Modi, die sich wiederholen, die auf den Modus selbst wirken; siehe [Weitere
Färbungsmodi](#weitere-färbungsmodi). OK setzt alles als einen einzigen
Rückgängig-Schritt. Für die Modi, die zählen — iter, smooth und die übrigen —,
ist `v` die Iterationszahl selbst, also kann man eine Funktion in Iterationen
denken. Jede ist ungerade, `f(−v) = −f(v)`, also behält eine Größe, die ihr
Vorzeichen wechselt, es.

Das waren vier Einträge des Menüs: die Geschwindigkeit, die Verschiebung, ein
Untermenü der Funktionen und der eigene Dialog des Geschwindigkeitswechsels.
Eine Position wird genau wie bisher geschrieben und gelesen, und die Befehle,
die diese Einträge ausführten, gibt es für ein Skript weiterhin.

`sin`, `cos` und `triangle` laufen im Kreis: die Färbungsgeschwindigkeit, nach
ihnen angewandt, sagt, wie viele Paletteneinträge eine Runde umfasst, und
lässt, wie oft sie wiederkehrt, wie es war. Das macht die
**Anfangsgeschwindigkeit** des Geschwindigkeitswechsels, vor ihnen angewandt:
mit beiden Gruppen aus und dem Anfangsexponenten bei eins bekommt die Funktion
`s0·v`, also rückt 2 die Runden doppelt so dicht zusammen und 0,5 doppelt so
weit auseinander.

| | |
| --- | --- |
| `ident`, `log`, `atan`, `sin`, `cos`, `square`, `cube`, `square root`, `cube root` | was es immer gab, unverändert bis aufs Bit |
| `log log` | `log(1 + log(1 + v))`, für die tiefen Zooms, wo log allein die Bänder noch zu eng packt: bei einer Million Iterationen gibt log 13.8 Bänder und diese 2.7 |
| `exp` | `e^v − 1`, die andere Richtung als log: sie öffnet die kleinen Werte. Für die Größen, die klein bleiben — die Komponenten, die Echtfarbe und die Zerlegungen lesen; bei einer Zählung von neunzig ist sie schon über 10³⁸ |
| `triangle` | der Sinus mit begradigten Kurven: dieselbe Periode und Höhe, die Palette in gleichmäßigem Tempo vor und zurück abgeschritten, statt an den Wenden zu verweilen |

Die neuen kommen nach den alten, also nennt jede Nummer, die eine gespeicherte
Position enthält, weiter, was sie nannte.

### Geschwindigkeitswechsel nach ITER

Der untere Teil desselben Dialogs ist keine Funktion, sondern eine Änderung an
der gewählten Funktion: bis zu einer Iteration, die Sie nennen, läuft die
Funktion genau weiter wie bisher; danach läuft die Zählung, die sie bekommt,
schneller oder langsamer, also werden die Farben von dort an schneller oder
langsamer und nicht vorher. Mit `square root` gewählt und ITER bei 40 sind die
Farben bis 40 Iterationen die der Quadratwurzel und danach die Quadratwurzel
einer Zählung, die an Tempo gewinnt. Er enthält einen Anfangsexponenten und
eine Anfangsgeschwindigkeit und dann zwei Gruppen, jede durch eine Linie
abgesetzt:

| | |
| --- | --- |
| **Anfangsexponent** | der Exponent vom Anfang bis zum ersten ITER einer eingeschalteten Gruppe: die Zählung wird zur einfachen Potenz `v^p0`, gleich, ob die Gruppen ein- oder ausgeschaltet sind, und mit beiden aus gilt sie für jede Zählung. Bei eins, wo er beginnt, bleibt die Zählung bis zum ersten ITER unberührt |
| **Anfangsgeschwindigkeit** | womit das Tempo dieses Abschnitts multipliziert wird: `s0·v^p0`. Bei eins, wo sie beginnt, ändert sie nichts |
| **ITER** | wo der Abschnitt der Gruppe beginnt, in den Einheiten, die die Funktion bekommt. Bei null oder darunter beginnt er am Anfang — was `ident` zu einer Potenz mit beliebigem Exponenten macht |
| **Exponent** | der Exponent des Abschnitts, von seinem ITER bis zum nächsten oder weiter: wie sich der Abschnitt biegt. Er beginnt bei der Zählung, bis zu der der vorige kam, `f(ITER)`, und addiert `s·ITER·((v/ITER)^p − 1)`, dessen Tempo bei ITER `s·p` ist: über eins werden die Farben von dort an schneller, unter eins langsamer |
| **Geschwindigkeit** | womit das Tempo des Abschnitts multipliziert wird — linear, wo der Exponent eine Potenz ist: bei einem Exponenten von eins ist der Abschnitt die Gerade `f(ITER) + s·(v − ITER)`. Bei eins, wo sie beginnt, ändert sie nichts |
| **Modus** | **off**, **jump** oder **smooth**: ob die Gruppe benutzt wird und wie ihr Abschnitt an den vorigen anschließt |
| **Radius** | nur bei smooth, und nur dann gezeigt: wie viele Iterationen der Anschluss auf jeder Seite von ITER einnimmt. Jede Zahl ab null, was der jump ist |

Mit beiden Gruppen ein und dem ersten ITER unter dem zweiten ist die Zählung
`v^p0` bis `i1`, läuft von dort im Tempo von `p1` weiter und ab `i2` im Tempo
von `p2` — **jeder Abschnitt mit eigenem Exponenten und eigener
Geschwindigkeit**, nie `p0·p1` oder `p0·p1·p2`. Mit `p0` = 2, `p1` = 1 bei 50
und `p2` = 2 bei 100, alle in jump und jede Geschwindigkeit eins, ist die
Zählung auf `ident` bis 50 `v²`; von 50 bis 100 `v + 2450`, ein Exponent von
eins, der im Tempo eins von den 2500 weiterläuft, bis zu denen das Quadrat kam;
und ab 100 `2550 + 100·((v/100)² − 1)`, wieder quadriert ab einem Tempo von 2
dort. Eine Geschwindigkeit von 3 in der ersten Gruppe macht den mittleren
Abschnitt zu `2500 + 3·(v − 50)`. Das niedrigere ITER geht zuerst, welche
Gruppe es auch hält, und **beim selben ITER ersetzt die zweite Gruppe die
erste**, Radius, Geschwindigkeit und alles; eine Gruppe mit einem ITER von null
oder darunter setzt ihren eigenen Exponenten und ihre eigene Geschwindigkeit
auf den ersten Abschnitt.

**Eine ausgeschaltete Gruppe fällt weg**, ihr ITER mit ihr:

| erste Gruppe | zweite Gruppe | die Zählung läuft mit |
| --- | --- | --- |
| ein | ein | `p0` bis zum niedrigeren ITER, dann dem Exponenten dieser Gruppe bis zum anderen ITER, dann dem der anderen |
| aus | ein | `p0` bis zum ITER der zweiten, dann `p2` |
| ein | aus | `p0` bis zum ITER der ersten, dann `p1` |
| aus | aus | `p0` für jede Zählung, der einfachen Potenz `v^p0` |

Beide beginnen ausgeschaltet, und der Anfangsexponent und jede Geschwindigkeit
bei eins, wo mit dem Wert gar nichts geschieht, also gibt jede Funktion bis
aufs Bit, was sie gab, bevor es so etwas gab.

Die beiden Arten des Anschlusses, wie sie die Zählung ändern, die der Funktion
übergeben wird (geschrieben für `ident`, wo die Zählung das ist, was die
Palette sieht, mit `f(ITER)` der Zählung beim Erreichen von ITER — `ITER^p0`
beim ersten —, `p` dem Exponenten der Gruppe, `s` ihrer Geschwindigkeit und `r`
ihrem Radius):

| Modus | von ITER − r bis ITER + r | danach |
| --- | --- | --- |
| **jump** | — | `f(ITER) + s·ITER·((v/ITER)^p − 1)`: das Tempo knickt bei ITER auf einmal ab, was sich als Ring zeigt, wo die Bänder ihr Tempo ändern |
| **smooth** | das Tempo läuft gerade von dem, das die Zählung bei `ITER − r` hat, zu dem, das sie bei `ITER + r` hat, auf halbem Weg bei ITER, mit einer hinzugefügten Parabel um ITER, die zurückgibt, was die Gerade von der Ecke des jump abschnitt | die Zählung des jump, bis auf die Ziffer |

**Smooth** ist symmetrisch um ITER und ändert die Farben innerhalb des Radius
und nirgends sonst: vor `ITER − r` und nach `ITER + r` ist das Bild das des
jump. Wo beide Tempi gleichmäßig sind, ist die Parabel null, und die Ecke wird
abgeschnitten, wie eine Parabel sie abschneidet — mit der Zählung unberührt bis
zu einem ITER von 50 und einem Exponenten von null danach macht ein Radius von
10 sie bei ITER zu 47.5, ein Viertel des Radius weniger, und ab 60 wieder zu
50; von einem Tempo von eins zu einem von drei ist es bei ITER zwei. Wo sich
ein Tempo biegt, gibt die Parabel zurück, was die Biegung hinzufügt: nach `v²`
bis 50 und einem Exponenten von eins führt ein Radius von 20 das Tempo von 60
bei 30 zu 1 bei 70, 45.5 bei ITER, und läuft ab 70 als `v + 2450` weiter, wie
jump es tut.

Früher mischte es die Tempi der beiden Abschnitte, wie jeder weitergelaufen
wäre: das Tempo des Quadrats, noch wachsend, wuchs weit in die Rundung hinein,
also fiel es vor ITER kaum und war erst einen halben Radius danach auf halbem
Weg nach unten, und was es gewann, verschob jede Farbe hinter der Rundung — um
133 Bänder bei einem Radius von 20.

Die Parabel bringt das Tempo nie unter das langsamere der beiden Enden. Wo sie
es täte — ein langsames Tempo, über einen kurzen Radius sehr schnell gemacht,
wie eine hundertste Potenz über zehn Iterationen —, wird es dort gehalten, und
die Zählung hinter der Rundung wird um das verschoben, was das zurückhält. Ein
Radius, der breiter ist als das erste ITER, hört vor null auf, und zwei, die
sich überlappen würden, werden gleich gekürzt, sodass jeder zwischen seinen
Nachbarn bleibt — Radien von 40 und 40 zwischen ITERs im Abstand 30 werden
je 15.

Ein Exponent von 100 ist erlaubt und nach ITER nutzlos: ein paar Iterationen
weiter sind die Bänder feiner als die Pixel, und die Farbe ist Rauschen. Mit
der zweiten Gruppe bei ITER 50 und einem Exponenten von 2 und der ersten bei
100 und 1.5, beide in jump, und einem Anfangsexponenten von eins ist die
Zählung auf `ident` 30 bei 30, `50·1.5²` = 112.5 bei 75 und
`200 + 100·(1.5^1.5 − 1)` = 284 bei 150: der Abschnitt ab 100 beginnt bei den
200, bis zu denen der vorige kam, mit einem Tempo von 1.5.

Früher gab es vier Modi, und vor dem Radius gespeicherte Positionen nennen sie
— `outcolorfunafter` aus 1.8.0 und `outcolorfunchange` mit zwei Gruppen und
ohne Radius. Sie werden in die zwei gelesen, die es gibt, mit jeder
Geschwindigkeit bei eins:

| gespeichert als | gelesen als |
| --- | --- |
| off, jump | dasselbe |
| smooth, das die Ecke im Logarithmus der Zählung rundete, ein wenig auf jeder Seite von ITER, wie weit auch immer | smooth mit einem Radius von drei Vierteln von ITER: innerhalb von vier Hundertsteln von ITER von dem, was es zeichnete, bis zu einem Exponenten von anderthalb, ein Zehntel bis zwei |
| continuous, das ab ITER mit dem Tempo weiterlief, mit dem die Zählung kam | smooth bei drei Vierteln von ITER, mit dem Exponenten `0.66 + p/4 + 0.09·p²` — eins bei eins, wo auch continuous nichts änderte —, innerhalb von zwölf Hundertsteln von ITER davon bis zu einhalb |

Gefunden, indem die Kurven über vier ITERs nebeneinandergelegt wurden: so ein
Bild kommt ähnlich heraus, nicht pixelgleich. Eine Position, deren
Geschwindigkeitswechsel nicht mehr als eine Gruppe in jump ist, mit dem
Anfangsexponenten und jeder Geschwindigkeit bei eins und dem Rest wie am
Anfang, wird gespeichert, wie 1.8.0 sie speicherte,
`(outcolorfunafter 40 1.5 'jump)`, und öffnet sich in 1.8.0 und zeichnet dort
dasselbe. Alles andere wird gespeichert, wie 2.0 es speicherte, in der
Reihenfolge des Dialogs ohne die Geschwindigkeiten,
`(outcolorspeedchange 1 40 1.5 'smooth 10 50 2 'off 10)`; und die
Geschwindigkeiten, nur wo eine davon nicht eins ist, getrennt als
`(outcolorspeedlinear 1 3 1)` — die Anfangsgeschwindigkeit und die der beiden
Gruppen —, sodass eine Position, die sie nicht benutzt, sich weiterhin in 2.0
öffnet. Der Geschwindigkeitswechsel allein, wie ihn sein eigener Dialog hatte,
ist der Befehl `outcolorspeedgroups` mit seinen zwölf Feldern in ihrer
Reihenfolge.

### Jenseits des Endes eines int

Die Palette wird in 256steln eines Eintrags über einen int erreicht, und ein
Wert über das hinaus, was ein int fasst — 8.4 Millionen Bänder —, wurde
trotzdem umgewandelt, was undefiniert ist: jedes solche Pixel bekam ein und
denselben Eintrag. `cube` erreicht das bei 204 Iterationen, also war ein mit
dem Kubus gefärbtes Bild darüber hinaus ein flacher Fleck; `exp` und ein großer
Geschwindigkeitswechsel erreichen es sofort. Jetzt landet er, wo sein Wert es
sagt, um die Palette gewickelt wie jeder kleinere Wert; einer, der zu groß ist,
als dass seine 256stel in einen double passten, landet auf dem ersten Eintrag,
wie eine Unendlichkeit es immer tat. Innerhalb eines int verschiebt sich
nichts. Außerhalb ist die einzige Änderung an einem zuvor gespeicherten Bild
die weiche Färbung von Magnet und Magnet2 fern der Menge, die eine einzige
flache Farbe war und jetzt die Bänder zeigt, die die Zählung dort beschreibt.

Zwei kleinere Reparaturen kamen mit. Eine Geschwindigkeitsfunktion zu wählen,
verschob das Häkchen im Menü der inneren Färbung statt in ihrem eigenen. Und
**initstate**, womit jede Position beginnt, setzte die fünf Zahlen des fbm
nicht auf ihre Standardwerte zurück — eine Position speichert sie nur, wo sie
davon abweichen, also zeichnete eine, die nach einer anderen geladen wurde, die
sie verändert hatte, mit denen der anderen; jetzt setzt es sie zurück und
schaltet mit ihnen den Geschwindigkeitswechsel nach ITER aus.

### Die Zählung, ganz

Der Wert, den ein Färbungsmodus übergibt, verlor auf seinem Weg zur Palette an
zwei Stellen Ziffern. Die weiche Färbung machte ihre Zählung zu einem int aus
256steln einer Iteration, bevor Färbungsgeschwindigkeit, Funktion oder
Geschwindigkeitswechsel sie bekamen; und was die Funktion zurückgab, war ein
float, der jenseits von 65 536 Einträgen weniger als ein 256stel eines Eintrags
behält und jenseits von sechzehn Millionen nicht einmal einen ganzen. Keines
von beiden zeigt sich bei Geschwindigkeit eins, wo ein 256stel einer Iteration
ein 256stel eines Eintrags ist. Beide zeigen sich, sobald etwas die Farben
schnell laufen lässt: mit einem Geschwindigkeitswechsel, der die Zählung von
null auf zehntausend Einträge pro Iteration bringt, kam eine Ansicht bei
dreißig Iterationen, wo das Tempo 2500 war, in einfarbigen Streifen von acht
Pixel Breite heraus, jeder ein 256stel einer Iteration und zehn Einträge vom
nächsten entfernt, wo die Palette jeden Eintrag dazwischen hätte durchlaufen
sollen.

Jetzt läuft die Zählung ganz weiter, und Funktion, Färbungsgeschwindigkeit und
Geschwindigkeitswechsel werden im eigenen Typ des Builds berechnet, wie die
Formeln — long double im gewöhnlichen Build, `__float128` im Quad — und geben
zurück, was sie berechnen. Die beiden Builds färben so ähnlich wie zuvor: auf
den getesteten Bildern unterscheiden sie sich dort, wo sie es immer taten, in
einem tiefen Zoom, dessen Orbits Quad weiter berechnet, und nirgends sonst.

Zwei Dinge bleiben, wie sie waren.

* **Eine Unendlichkeit** — die Mitte eines ins Unendliche gegangenen
  Potentials, der Grund um eine Formel, die sofort austritt — landet auf dem
  ersten Eintrag, wie ein Wert, der keine Zahl ist. Als 256stel in einem float
  kam jeder solche Wert früher unendlich heraus und landete dort; welche Farbe
  eine Unendlichkeit annehmen sollte, weiß niemand, und so nimmt sie die, die
  sie immer nahm. Ein endlicher Wert jenseits dessen, was ein float fasst,
  landet, wo er es sagt.
* **Ein Wert auf der Kante eines 256stels** bleibt darauf: der Wert wird auf
  ein 256stel eines 256stels gerundet, bevor er auf ein 256stel geschnitten
  wird, da die Färbung nichts Feineres sagt. Der Realteil eines Newton-Orbits
  kommt auf einer Wurzel zur Ruhe, bei eins oder bei minus einhalb, was
  `iter+real` genau auf eine solche Kante setzt, und die letzten Ziffern der
  Konvergenz würden jedes Pixel auf die eine oder andere Seite setzen — ein
  256stel eines Eintrags, das niemand sieht, bis ein Kantenfilter jedes solche
  Pixel als Kante zeichnet und die Höhenlinien gezahnt wie Zahnräder
  herauskommen. Der float rundete die meisten davon weg und ließ Kerben; jetzt
  gibt es keine.

Gegen den vorigen Build bei 640×480 gerendert — vor der Änderung im nächsten
Abschnitt — zeichneten von den 169 Beispielen 83 pixelgenau wie zuvor, und 71
weitere unterschieden sich um eine einzige Stufe in höchstens 113 Pixeln. Von
den anderen fünfzehn unterscheiden sich die Echtfarb- und smooth-log-Beispiele
— `perb4`, `spiral4`, `truecol1`, `truecol3` — um bis zu fünf Stufen in bis zu
einem Siebtel der Pixel, da auch ihre Komponenten als floats zurückkamen;
`edge` hat seine Kerben verloren; und der Rest unterscheidet sich in einer
Handvoll Pixel, wo ein Wert auf der Kante eines Eintrags lag und ein Filter
oder eine ohne Überblendung gezeichnete Palette die Stufe zeigt. Magnet und
Magnet2 in weicher Färbung ändern sich am meisten: fern der Menge ist ihre
Zählung negativ, was der int, vorzeichenlos gelesen, um viertausend Millionen
256stel weiterschob und der float dann in Stufen von ein oder zwei Einträgen
zeichnete — flache Flecken mit harten Kanten, wo jetzt der Verlauf ist, den die
Zählung beschreibt.

Auf einem Bild, das gleich herauskommt, kostet es zwei bis fünf Prozent im
gewöhnlichen Build und sieben bis zwanzig im Quad, wo die Logarithmen der
Zählung und die Funktionen in Software berechnet werden. Wo die Streifen zu
einem Verlauf werden, dauert das Rendern auch länger: das Solid Guessing füllt
ein Pixel nur aus Nachbarn einer Farbe, und das geschriebene Bild ist größer —
die Ansicht bei dreißig Iterationen brauchte bei 1200×1200 1.07 Sekunden, wo
ihre Streifen 0.30 brauchten.

### Das Tempo der weichen Zählung

Die weiche Färbung addiert zur Zählung der ganzen Iterationen, wie weit in die
letzte hinein der Punkt austrat. Das berechnete sie dort, wo der Fluchtradius
zwischen die Größe vor dem austretenden Durchlauf und danach fällt, linear in
ihrem Logarithmus. So läuft die Zählung ungebrochen von einer Iteration zur
nächsten, ihr Tempo aber nicht: von einem Durchlauf zum nächsten wird die Größe
fast quadriert, also ändert sich das Tempo über jede Iteration und knickt auf
einmal ab, wo die nächste beginnt. In einer Ansicht von 2.6×10⁻⁷ Breite, wo die
Zählung von 59.94 bis 60.02 geht, knickte es an der Kante zwischen beiden um
4.26 ab. Bei einer Färbungsgeschwindigkeit von eins verschiebt das eine Farbe
um ein paar Hundertstel eines Eintrags; mit einem Geschwindigkeitswechsel, der
zehntausend Einträge pro Iteration läuft, sind es die Bänder, die entlang jeder
Kante einer Iteration schlagartig ihre Dichte ändern.

Jetzt wird der Orbit über den Fluchtradius hinaus weitergeführt, bis sein
Betragsquadrat eine Million übersteigt — drei oder vier Durchläufe mehr bei
einem Quadrat —, und der Bruchteil wird im Logarithmus des Logarithmus
genommen, wo ein Quadrieren ein Schritt der Länge eins ist:

    n + ln(ln T / ln s) / ln(ln m / ln s)

mit `s` und `m` dem Betragsquadrat zu beiden Seiten einer Million, `T` dem
Fluchtradius und `n` den gemachten Durchläufen. Die Durchläufe über den
Fluchtradius hinaus ziehen sich selbst wieder ab, und der Grad ist der des
Orbits selbst, `ln m / ln s`, also wird ein Kubus als Kubus gezählt. Die
Zählung bleibt innerhalb eines Zehntels einer Iteration von dem, was sie war,
und ihr Tempo läuft weiter: der Knick von 4.26 in dieser Ansicht kommt bei
1.0000 heraus, und entlang eines Strahls von außerhalb der Menge die Knicke von
1.3 bis 1.8 an jeder Kante bei 1.0025.

| | die Zählung |
| --- | --- |
| die Mandelbrots, Triceratops, Mandelbar, Lambda, Manowar, Spider, Phoenix und Octo, auf dem Kreis | über den Fluchtradius hinaus weitergeführt |
| eine Benutzerformel, die auf dem Kreis austritt | weitergeführt, solange der Orbit mindestens wie eine Potenz von anderthalb wächst, und am Fluchtradius gezählt, wo nicht |
| die drei Barnsleys | am Fluchtradius: sie multiplizieren den Orbit mit `c`, sodass sein Logarithmus um einen Schritt wächst statt sich zu verdoppeln, und die Zählung am Fluchtradius ist für sie schon die richtige |
| Newton, Magnet, Catseye | am Fluchtradius: sie konvergieren oder prüfen ihr Austreten auf eigene Weise |
| jede Form des Fluchtradius außer dem Kreis | am Fluchtradius |

Die Formen bleiben, wie sie waren, weil jede den Orbit mit etwas misst, das
sich mit ihm dreht, während jeder Durchlauf seinen Winkel verdoppelt: mit
eigenem Maß weitergeführt, kam der Grad so unstet heraus, dass die Knicke
schlimmer wurden — die des Sechsecks von 1.28 auf zwischen 1.67 und 3.57
entlang dieses Strahls —, und die des Quadrats kehrten sogar um.

Von den 169 Beispielen zeichnen die 154, die nicht weich färben, pixelgenau wie
vor dieser Änderung. Die fünfzehn, die es tun, ändern sich um bis zu dreizehn
Stufen in bis zu neun Zehnteln der Pixel — `eye` um bis zu 39, wo die
benachbarten Einträge seiner Palette weit auseinanderliegen, und `spiral4`,
dessen smooth log auf einer Echtfarbpalette jede Änderung der Zählung
vergrößert, mitten in einem kleinen Detail. Gegen 2.1, mit den Änderungen des
vorigen Abschnitts, zeichnen 84 pixelgenau wie zuvor, und 59 weitere
unterscheiden sich um höchstens drei Stufen in höchstens 500 Pixeln.

Es kostet diese Durchläufe und zwei Logarithmen mehr für jedes Pixel, das
austritt: auf den getesteten weichen Bildern 11 bis 28 Prozent im gewöhnlichen
Build und 18 bis 20 im Quad, wo ein Logarithmus in `__float128` 946 ns braucht
gegen 39 in long double; ein Bild, das nicht weich färbt, ist unverändert. In
double berechnet, hätten die vier Logarithmen im gewöhnlichen Build 8 bis 19
Prozent gekostet und im Quad 8 gespart — ein double trägt den Bruchteil auf
sechzehn Stellen, wo eine Geschwindigkeit von zehntausend neun braucht —, aber
die Färbung wird durchweg im eigenen Typ des Builds berechnet, wie die Formeln.

## Zoomen

**Auswahlzoom** (Ausschnitt → Zoom) ersetzt den stufenlosen Zoom durch ein
Rechteck: mit der linken Taste aufziehen und loslassen, um dorthin zu gehen,
wie bei Ultra Fractal. Nur hinein. Das Rechteck wird auf das Seitenverhältnis
des Fensters gebracht, indem seine kurze Seite verlängert wird, nie indem die
lange beschnitten wird, und es funktioniert auch bei gedrehtem Fraktal richtig.
Der schnelle Juliamodus übernimmt das Ziehen, solange er läuft, und gibt es
beim Verlassen zurück.

**Feste Schritte**, im selben Untermenü:

| | |
| --- | --- |
| 2× vergrößern | `Strg` `+` |
| 2× verkleinern | `Strg` `-` |
| 10× vergrößern | `Strg` `Alt` `+` |
| 10× verkleinern | `Strg` `Alt` `-` |

Sie skalieren um die Mitte der Ansicht, also heben sie einander exakt auf.

## Die Ansicht drehen und zentrieren

**Fraktal → Ausschnitt** setzt Zentrum, Radius und Winkel als Zahlen. Unter dem
Winkel hat er zwei Werkzeuge für ebendiese Felder. Sie füllen die Felder aus
und tun sonst nichts: OK wendet an, was sie zeigen, als einen einzigen
Rückgängig-Schritt, Abbrechen nichts davon, und eine gespeicherte Position
erfährt nie von ihnen.

**Stücke** teilt eine volle Umdrehung in so viele gleiche, und **−** und **+**
drehen den Winkel um eines zurück oder weiter: 60 Grad bei sechs, 360 bei
einem. Bei keinem sind die beiden deaktiviert. Die Drehungen werden vom Winkel
aus gezählt, bei dem sie begannen, also sind sechs Sechstel eine ganze
Umdrehung bis auf die letzte Ziffer; einen Winkel eintippen oder die Stücke
ändern, beginnt die Zählung von dort neu. Die Zahl der Stücke wird behalten,
solange XaoS läuft, für das nächste Öffnen des Dialogs.

**An der Achse zentrieren** verschiebt das Zentrum geradewegs quer über das
Bild (**R**) oder geradewegs daran hinauf (**I**), bis der Ursprung in dieser
Richtung in der Mitte liegt, beim Winkel im Dialog. Ungedreht heißt das, dass
der Real- oder der Imaginärteil null wird und der andere bleibt. Gedreht ist es
dasselbe entlang der eigenen Achsen des Bildes: R bei 90 Grad löscht den
Imaginärteil, und R und dann I bringen den Ursprung bei jedem Winkel in die
Mitte. Ein Zentrum, das schon dort ist, bleibt unberührt, und ein Teil, den die
Drehung als Rundungsrest lässt, wird als null genommen, also zeigt das Feld 0
und nicht 1e-17.

Die Drehung selbst wird im eigenen Typ des Builds berechnet, samt Pi. Sie wurde
aus einem double-Pi berechnet, was im Quad-Build auch einen double-Sinus und
-Kosinus bedeutete: ein Paar, das nur bis auf sechzehn Stellen eine Drehung
ist, also wurde das Zentrum einer gedrehten Ansicht um 5e-17 seines Abstands
vom Ursprung neben seinem Ort gezeichnet. Unter 1e-16, wofür Quad da ist,
zeigte eine um 30 Grad gedrehte Ansicht einen ganz anderen Ort. Eine
Vierteldrehung ist jetzt exakt. Im gewöhnlichen Build ist ein gedrehtes Bild
dasselbe bis auf Pixel an der Grenze seiner Genauigkeit: 52 der 64 Beispiele
mit einem Winkel sind unverändert, die übrigen verschieben verstreute Pixel, wo
das Detail feiner ist als die Arithmetik.

## Dasselbe Bild mit jeder Anzahl von Threads

Ein Bild kam hin und wieder mit ein paar anderen Pixeln heraus, wenn es mit
mehreren Threads gezeichnet wurde: `examples/Hubicka/julia3.xpf`, mit vier
Threads auf einer ausgelasteten Maschine gerendert, unterschied sich in 26 von
3000 Renderings, um bis zu 37 Pixel entlang seiner Filamente, und nie mit einem
Thread.

Die Pixel kamen vom Solid Guessing. Es füllt ein Pixel aus den bereits
berechneten Zeilen und Spalten um es herum, wenn diese übereinstimmen, und
nutzt im Mandelbrotmodus dieselben Nachbarn, um zu entscheiden, ob es sich
lohnt, ein Pixel auf Periodizität zu prüfen. Das Raten läuft in Durchgängen —
erst jede sechste Zeile, dann die Zeilen dazwischen —, und die Durchgänge
liefen in allen Threads zugleich ohne etwas dazwischen: ein Thread ging zum
nächsten Durchgang über, während die anderen noch beim letzten waren, und eine
Zeile wurde als erledigt markiert, bevor sie berechnet war. Ob ein Nachbar als
berechnet galt, hing also davon ab, wie weit die anderen Threads gekommen
waren.

Jeder Durchgang hat jetzt zwei Phasen, erst seine Zeilen, dann seine Spalten,
und jeder Thread beendet eine, bevor irgendeiner die nächste beginnt; was eine
Phase berechnet, wird innerhalb dieser nicht als Nachbar benutzt. Die Nachbarn
sind dann dieselben, wie viele Threads es auch sind, und die 3000 Renderings
kommen alle gleich heraus. Mit einem Thread ist die Reihenfolge die, die sie
immer war.

Auch die Filter kommen gleich heraus. **3D** zeichnet ein Band von Spalten pro
Thread, und die Schattierung einer Spalte wird von der Spalte zwei rechts davon
in der Zeile darunter übernommen, also lesen die letzten beiden Spalten eines
Bandes die ersten beiden des Bandes rechts davon; früher las jedes Band dort
Kopien seiner eigenen letzten Spalte, und das Bild änderte sich entlang der
Bandkanten mit der Zahl der Threads. Jetzt wartet ein Band Zeile für Zeile auf
das, was das Band rechts davon dort gezeichnet hat, und die Bänder laufen
weiterhin nebeneinander. **Sternenfeld** zog seine Sterne aus einem einzigen
Zufallsgenerator, den jeder Thread zugleich weiterschaltete; jetzt beginnt jede
Zeile an der Stelle der Folge, die ein einziger Thread, der das ganze Bild
zeichnet, erreicht hätte, also sind die Sterne mit jeder Anzahl von Threads
dieselben. Sie ändern sich weiterhin von einem Start von XaoS zum nächsten, wie
sie es immer taten: der Generator bekommt seinen Anfangswert von der Uhr.

**Antialiasing** war nie die Ursache, obwohl drei damit gezeichnete Beispiele —
`circle`, `helloween` und `warriormask` — von einem Rendering zum nächsten
anders herauskamen. Was sich unterschied, war die Meldung mit dem Namen der
Formel: sie verschwand nach fünf Sekunden der Wanduhr, und ein von einer
ausgelasteten Maschine gebremstes Rendering zeichnete sie in weniger Bildern
oder in keinem. In einem Rendering dauert eine Meldung jetzt fünf Sekunden der
eigenen Zeit der Animation, welche Maschine es auch ist; im Fenster ändert sich
nichts.

Von den 169 mit einem Thread gerenderten Beispielen zeichnen 167 pixelgenau so,
wie 1.8.0 sie zeichnete: `galaxy` hat seine Sterne von der Uhr, und `heart`
behält jetzt eine Meldung, die 1.8.0, nach der Wanduhr getaktet, entfernt
hatte, bevor das Bild gezeichnet war. Zu acht auf vier Kernen gerendert, einmal
mit einem Thread und dreimal mit vier, zeichnen alle gleich bis auf `galaxy`.

## Seitenverhältnis

**Ausschnitt → Seitenverhältnis** ändert die Fenstergröße einmal, damit das
Bild die gewählte Form annimmt: 1:1, 4:3, 3:2, 16:10, 16:9, 2:1, 21:9 oder 9:16
hochkant. Das Fenster wird nicht darauf festgehalten: es kann danach auf jede
Form gezogen werden, und das Häkchen im Menü steht auf der Form, die das Bild
gerade hat, oder auf keiner.

Das Bild behält die Zahl seiner Pixel und ändert nur seine Form — oder wird,
wenn es so nicht auf den Bildschirm passen würde, das größte Bild dieser Form,
das passt. Aus dem Vollbild oder einem maximierten Fenster kommt es zuerst in
ein gewöhnliches Fenster zurück, noch mit den Pixeln, die es hatte, also ist
16:9 dort das breiteste 16:9-Fenster, das der Bildschirm fasst. Die Form wird
pixelgenau getroffen, auch auf einem skalierten Bildschirm: bei 125 % ist eine
von fünf Pixelzahlen nicht zu haben, und die Größe wird unter denen gewählt,
die es sind.

Wofür es da ist, ist das Rendering. Eine Position wird so gezeichnet, dass die
ganze Ansicht, die sie hält, im Bild liegt, und ein Bild anderer Form zeigt an
seiner längeren Seite mehr von der Ebene. Also hat eine Ansicht, die im
gewöhnlichen 4:3-Fenster eingerahmt und bei 1920×1080 gerendert wird, an den
Seiten ein Drittel mehr, als je auf dem Bildschirm war; in einem 16:9-Fenster
eingerahmt, ist das Rendering das, was man sah.

## Antialiasing in linearem Licht

Filter → Antialiasing in linearem Licht, oder `-linearaa`. Gewöhnliches
Antialiasing mittelt die Farben der Abtastpunkte, wie sie gespeichert sind, was
Kanten abdunkelt, weil die gespeicherten Werte nicht proportional zum Licht
sind. Dies wandelt in lineares Licht um, mittelt dort und wandelt zurück, also
behält eine Kante die Helligkeit dessen, wozwischen sie liegt.

## Sprachen

Ausschnitt → Sprache einstellen schaltet XaoS nach einem Neustart auf eine
andere Sprache um. Sechs sind vollständig: Deutsch, Spanisch, Französisch,
Italienisch, Portugiesisch und Russisch — die Menüs, die Dialoge und die
Meldungen, die Referenz der Benutzerformeln und die Dialoge OK, Abbrechen,
Datei und Farbe, die von Qt selbst kommen. Letztere sind jetzt eingebaut, wo
sie früher dort gesucht wurden, wo Qt installiert ist, und in der Sprache des
Systems, welche auch gewählt war. Die anderen Sprachen, die XaoS geerbt hat,
reichen so weit, wie ihre Übersetzer sie gebracht haben, und sind darüber
hinaus Englisch.

In der Referenz der Benutzerformeln sind die Beschreibungen übersetzt und die
Namen und Argumentlisten nicht: sie sind die Syntax, in der eine Formel
geschrieben wird.

Diese Anleitung ist in dieselben sechs Sprachen übersetzt, neben der englischen
`ms-xaos-guide.md` als `ms-xaos-guide.de.md`, `.es.md`, `.fr.md`, `.it.md`,
`.pt.md` und `.ru.md`, und Hilfe öffnet die für die Sprache, in der XaoS ist.
Jede nennt die Version der englischen Anleitung, aus der sie gemacht wurde; wo
sie sich unterscheiden, ist die englische maßgeblich.

## Bauen

    cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build

Beide Programme landen in `bin/`. Optionen:

| | |
| --- | --- |
| `-DQUAD=OFF` | das zweite Programm überspringen |
| `-DDEEPZOOM=ON` | das Hauptprogramm zu Quad machen und kein zweites bauen |
| `-DDIAG=ON` | auch `XaoS-diag` bauen, das Überlagerungs- und Berechnungsreihenfolge nach `xaos-diag.txt` meldet |

`cmake --build build --target deploy` füllt `bin/` mit der Laufzeit von Qt und
dem Compiler, damit es auf einer Maschine ohne Qt läuft.

## Tests

    ctest --test-dir build

351 Tests: der Formelparser in beiden Genauigkeiten, die Genauigkeit von
`gamma`, `lambertw`, `erf` und `randsc` gegen exakte Referenzen, wo es welche
gibt, die Iterationsschleifen gegen aufgezeichnete Prüfsummen, die
Geschwindigkeitsfunktionen der Färbung gegen das, was sie vorher ergaben, und
das, wofür sie da sind, gespeicherte Positionen hin und zurück und jedes
Rendering, das sauber endet, das Speichern und Wiederherstellen der
Überlagerung, der Rückgängig-Verlauf, das Drehen und Zentrieren des
Ausschnitt-Dialogs, die Hilfereferenz gegen die eigene Tabelle des Parsers und
das zufällige `poly` ihres Reiters Werkzeuge gegen das Glied für Glied
ausgeschriebene Polynom.

Einer davon misst Geschwindigkeit statt Antworten. Zeitmessung in Nanosekunden
würde behaupten, auf welcher Maschine die Suite läuft, also behauptet er
stattdessen ein Verhältnis: was eine Degradation die Rauschfunktionen kostet,
gezählt in Multiplikationen. Es sollten zwei sein, eine pro Komponente, weil
die Größe einmal pro Durchlauf mit der Degradation multipliziert wird; sie
stattdessen in die Potenz des Durchlaufs zu erheben, ergibt sechs oder sieben
und steigt mit dem Iterationslimit. Der Test schlägt bei den beiden
Implementierungen fehl, die diese ersetzt hat.
