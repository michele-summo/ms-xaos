# MS XaoS

> Traducción de la [guía en inglés](ms-xaos-guide.md) de la versión 2.10. Donde
> difieran, la referencia es el original en inglés.

Un fork de [XaoS](https://github.com/xaos-project/XaoS) 4.3.3. Todo lo que hace
el original lo sigue haciendo; esta guía describe lo que se ha añadido o
cambiado, y por qué. Versión 2.4.

## Dos ejecutables

`XaoS.exe` calcula con 64 bits de mantisa, `XaoS-quad.exe` con 113. Ambos se
compilan a partir de las mismas fuentes — la precisión es un tipo elegido al
compilar, no un ajuste — y el quad se llama **MS XaoS Quad**, para que los dos
se distingan mientras están abiertos.

Quad hace zoom unos quince órdenes de magnitud más hondo. Al renderizar la
misma posición a 1e-18, el ejecutable normal conserva el 21% del detalle y el
quad el 92%. No es el predeterminado porque `__float128` se emula por software:
unas ocho veces más lento en las fórmulas incluidas y quince en las definidas
por el usuario.

Una posición guardada por el ejecutable quad lleva `(precision 113)`. Abrirla
con una compilación de menos precisión lo avisa y la dibuja igualmente; abrirla
con una versión que no conoce el comando se rechaza, y es lo que se pretende —
una imagen calculada a 113 bits no se puede reproducir a 64, y no abrirla es
más honesto que dibujar otra cosa.

## Los números en los diálogos

Un campo que pide un número acepta uno escrito de cualquier forma legible,
1e-18 igual que 0.5. Los diálogos de un solo campo usaban una casilla
giratoria, que guarda dos decimales si no se le dice otra cosa, así que 0.01
sobrevivía y 0.001 se convertía en nada — una convergencia de Newton de una
millonésima no se podía introducir en absoluto.

Lo que se muestra es el texto más corto que se vuelve a leer como el número
exacto que se guarda, así que una millonésima de millonésima aparece como 1e-18
y no como 9.99999999999999999978e-19, que es lo que imprimir todas las cifras
de la compilación hace con un decimal que un float binario no puede representar
exactamente. Una coordenada que necesita las veintiuna cifras las recibe
igualmente: necesitarlas es lo mismo que no volver a leerse sin cambios sin
ellas.

La forma exponencial se usa a partir de un billón y por debajo de cien
milmillonésimas, y la decimal simple entre medias — así veinte es veinte, cien
mil millones se escribe entero, y 1e-8 no. Si se la dejara sola, la impresión
elegiría comparando el exponente con el número de cifras pedidas, lo que no
sirve cuando las cifras se eligen por brevedad: basta una cifra para volver a
leerse como veinte, y veinte salía como 2e+01.

## Posiciones guardadas

La vista se escribe con la precisión completa de la compilación — 21 cifras
significativas normalmente, 36 en quad. Antes se escribía con un número de
cifras derivado de la profundidad del zoom, que se quedaba corto y llegaba como
mucho a 17, así que una posición profunda se volvía a abrir a una fracción de
píxel de donde se había guardado. Con muchas iteraciones eso es una imagen
visiblemente distinta.

Todo lo que una posición puede contener sobrevive ahora sin cambios a guardar y
cargar: los 169 ejemplos incluidos y 141 posiciones generadas que cubren cada
fórmula, plano, modo de coloreado y filtro hacen el viaje de ida y vuelta píxel
a píxel.

## Deshacer

Editar → **Deshacer** (`U`, o `Ctrl` `Z`) retrocede un paso y **Rehacer**
(`Ctrl` `Y` en Windows) avanza de nuevo. Un paso es cualquier cosa que cambie
la imagen: un comando — una fórmula, una paleta, un modo de coloreado, un
diálogo, un zoom fijo —, una pulsación del ratón con el zoom, el desplazamiento
o la rotación que pone en marcha, un zoom por selección, y las iteraciones
fijadas con las flechas. Un comando que no cambia nada, como Recalcular, no es
un paso, y deshacer no retrocede más allá de la imagen con la que se abrió
XaoS, después de la configuración y la línea de comandos. Cada paso se guarda
como una posición guardada, así que deshacer devuelve todo lo que una posición
contiene, tan exactamente como volver a abrirla; se guardan los últimos 256.

Fallaba de dos maneras. El ratón guardaba su paso al soltar, mientras el zoom
aún estaba frenando, así que deshacer devolvía una imagen a medio camino; ahora
el paso se guarda donde la vista se detiene. Y cada comando guardaba un paso,
cambiara algo o no, mientras deshacer retrocedía dos dando por hecho que había
guardado el último por sí mismo — tras un Recalcular devolvía la imagen que ya
estaba en pantalla. Ahora una imagen solo se guarda donde difiere de la última
guardada, y deshacer y rehacer se mueven un paso a lo largo de ellas. Un paso
nuevo tras deshacer corta lo que se podría haber rehecho, como en un editor.

Lo que se ha movido por sí solo desde el último paso — el piloto automático,
una rotación continua — se guarda como un paso propio al pulsar Deshacer, para
que Rehacer vuelva a ello.

Deshacer encontró un fallo en las posiciones. Cargar una con una paleta
incluida tras una paleta personalizada dejaba la paleta personalizada marcada
como en uso, así que la posición se volvía a guardar como paleta personalizada:
no lo que se había cargado, y para deshacer no el paso al que había vuelto, lo
que costaba el rehacer. Ahora una paleta incluida quita la marca, como siempre
hizo una generada.

## Fórmulas de usuario

**Referencia.** Ayuda → Referencia de fórmulas de usuario enumera cada función,
variable y notación que acepta el lenguaje de fórmulas, en una ventana que se
puede dejar abierta junto a la fórmula que se está escribiendo. Está agrupada
según para qué sirven las funciones, con las funciones especiales — `erf`,
`gamma`, `lambertw` — en una sección propia, ya que ninguna es elemental y no
se recurre a ellas por accidente. Una cuarta pestaña enumera los números que
aparecen como argumentos y significan algo concreto — contra qué forma mide una
trampa, con qué espejo pliega el caleidoscopio — en lugar de enterrarlos en la
descripción de la función que los recibe. Una quinta muestra una imagen de cada
teselado que admite `randsctile`, ya que cuarenta y cinco números con una línea
de palabras cada uno dicen poco de las formas.

Una sexta, **Utilidades**, escribe fórmulas en lugar de describirlas. Por ahora
contiene un `poly` aleatorio, a partir de cuatro campos:

| | |
| --- | --- |
| longitud mínima y máxima | m, el número de coeficientes, sorteado entre las dos; el grado es m − 1. Hacen falta ambas, como números enteros de 1 a 254 — hasta donde llega `poly`, ya que el parser admite 255 argumentos por llamada y `z` es uno de ellos |
| módulo mínimo y máximo | el tamaño de cada coeficiente, sorteado uniformemente entre los dos. 0 y 1 si se dejan vacíos; como máximo 1000000 |

Cada coeficiente es complejo, `{re,im}`, con su argumento sorteado al azar en
toda la vuelta, y cada parte se escribe con cuatro cifras significativas de su
módulo, para que un intervalo de 0.001 a 0.002 no se pierda al redondear. En la
fórmula no hay `c`. Una fórmula de usuario hace empezar `z` en el píxel tanto
en modo mandelbrot como en modo julia, así que lo que dibuja es el conjunto de
Julia relleno del propio polinomio, igual en ambos; `c` va donde se quiera,
como en `poly(...) + c`.

Cada cambio en un campo genera una fórmula nueva, igual que **Generar otro**, y
**Copiar** la pone en el portapapeles para Fractal → Fórmula de usuario.
Mientras una longitud esté vacía, o un campo no sea un número, esté fuera de su
intervalo o tenga el mínimo por encima del máximo, no hay fórmula, solo una
línea que dice por qué.

Cada fila de la lista de funciones dice qué recibe la llamada y en qué orden,
con nombres que dicen para qué es cada lugar — "base, exponent", "z, low,
high", "before, after, n" — en lugar de letras, que no decían qué lugar era
cuál; los mismos nombres son los que muestra la línea bajo la barra de
fórmulas. Una convención las recorre todas: `z` es el valor sobre el que
trabaja una función, sea lo que sea lo que esté escrito ahí. La lista y esos
recuentos se comprueban con un test frente a la propia tabla del parser, así
que una función que gana o pierde un argumento no puede seguir describiéndose
con los que tenía. Un argumento que se puede omitir se muestra entre corchetes
con el valor que toma cuando se omite: `z, [power=1], [angle_factor=1]`.

**Los argumentos se separan con una coma, y con nada más.** Antes también
servía un punto y coma: el parser lo convertía en coma antes de leer nada, así
que `LOGN(5;Z^2)+C` y `LOGN(5,Z^2)+C` eran una misma fórmula, y lo mismo pasaba
con las dos partes de un número complejo. Eso se acabó. Ahora un punto y coma
se rechaza donde aparece — *Operador no válido: ;* — y la coma es el separador
en todas partes, como lo es en el resto de la sintaxis.

### La sonda de paleta

**Vista → Sonda de paleta (haz clic en un punto)** convierte el puntero en una
sonda: haz clic en cualquier punto del fractal y la línea de mensajes dice **de
qué lugar de la paleta se dibujó ese punto**, como un número y no como un color
— *Posición en la paleta 12.375*, y nada más.

El editor de paleta coloca **31 colores**, y la paleta se construye pasando de
cada uno al siguiente **ocho celdas cada vez**. Así que el número es la celda
entre ocho: `12.000` es exactamente el duodécimo color del editor, `12.375`
está a tres octavos del camino del duodécimo al decimotercero. Una celda es un
octavo de color, más fino de lo que nadie lee en el editor.

* Con la **paleta personalizada** en uso el número se lleva al rango **`0` a
  `31`**, que es la numeración del propio editor.
* Con cualquier otra paleta se **deja donde cae y va mucho más allá de 31** —
  una paleta de 65534 celdas llega a 8191. El número fuera de rango es lo que
  dice que la paleta no es la que dibuja el editor.
* Un punto que tomó el **color interior** dice *Posición en la paleta:
  ninguna*: ese color es un único tono plano que ninguna celda de la paleta
  representa.

Sirve menos para saber "qué color es este" que **qué parte de la paleta usa de
verdad la imagen**, y la respuesta a menudo da que pensar. Con el coloreado
`iter` simple a 60 iteraciones las celdas van de 1 a 59 — de `0.125` a `7.375`
— así que la imagen nunca toca los colores 8 a 30 del editor y editarlos no
cambia nada. Un coloreado que se extiende, `iter + real` por ejemplo, recorre
en cambio todo el anillo, varias veces.

Lee el punto y no la pantalla, así que es exacta: el mismo color puede estar en
dos celdas, y una imagen en color real mezcla dos celdas vecinas, así que
ningún color nombra una celda propia. Como el zoom por selección, la sonda
devuelve el puntero al modo Julia rápido mientras este está activo.

**La barra.** Fractal → Fórmula de usuario abre un diálogo de aproximadamente
la mitad del ancho de una pantalla de 1280×720 — unos seiscientos cincuenta
píxeles, menos en una pantalla más pequeña — cuya barra pasa a una línea nueva
cuando la fórmula no cabe en una. Salta **entre nombres, nunca dentro de uno**:
un nombre o un número (`randsctile`, `z_p3`, `0.6`) se mantiene entero, y una
línea puede acabar en cualquier otro sitio — tras un paréntesis, una coma, un
operador. El ajuste de línea de Qt no servía para esto: corta donde lo permiten
las reglas de la prosa, que en una fórmula es tras unas pocas comas y poco más,
y luego corta un nombre donde haga falta. Solo se corta un nombre más ancho que
toda la barra. La barra crece una línea cada vez con la fórmula, y cuando el
diálogo ya no cabría en una pantalla de 720p deja de crecer y se desplaza. Para
todo lo demás sigue siendo una fórmula en una línea: Intro es Aceptar y no una
línea nueva, y en una fórmula pegada los saltos de línea y tabulaciones se
convierten en espacios, que el parser lee como nada.

**Historial...** bajo la barra abre las fórmulas usadas antes — las diez más
recientes, la más nueva primero, cada una una vez — que la antigua barra de una
línea guardaba en su lista desplegable. Una fórmula demasiado larga para la
ventana aparece con su comienzo, tres puntos y su longitud,
`poly(z, {0.01097,0.00699}, {-0.05... (157 caracteres)`, para que dos que
empiezan igual se puedan seguir distinguiendo; el comienzo acaba entre nombres,
como una línea de la barra, y se ve más cuanto más ancha es la ventana. Haz
doble clic en una, o selecciónala y pulsa Usar, y entra entera en la barra. Una
fórmula se añade cuando se acepta el diálogo con ella, y solo entonces: cargar
una posición que contiene una fórmula de usuario, reproducir una o renderizar
una con `-render` ejecuta el mismo comando `usrform` pero ya no la mete en el
historial, donde cada una de ellas ocupaba el lugar de una que habías escrito
tú.

**Flecha arriba y flecha abajo** recorren las mismas fórmulas sin salir de la
barra. Flecha arriba, desde la primera línea de la barra, trae la fórmula usada
antes de la que se ve; flecha abajo, desde la última línea, vuelve hacia la más
reciente y luego a lo que tenía la barra antes de empezar. Las dos flechas
pequeñas junto a la barra hacen lo mismo, cada una en gris cuando no hay nada
más en esa dirección. Una fórmula que ya está en la barra se salta, y lo que
escribes pasa a ser el punto desde el que se vuelve a empezar.

La última fórmula usada se muestra tenue en la barra: entera mientras la barra
está vacía, y el resto detrás de lo que has escrito mientras lo que escribes
sea su comienzo. Tab la pone entera en la barra. En cuanto lo que escribes va
por otro camino el texto tenue desaparece, y Tab pasa a los botones como antes.
La barra se abre vacía, con esa fórmula tenue dentro, y Aceptar con la barra
vacía la usa; solo la primerísima vez, cuando aún no se ha usado ninguna
fórmula, se abre con la fórmula en uso, `z^2+c`.

**Vista → Estado** muestra la fórmula de la misma manera: en una línea si cabe
en la ventana, en dos si no — la segunda bajo la primera, saltando entre
nombres — y cuando ni dos bastan, la segunda acaba en tres puntos y su
longitud, `... (918 caracteres)`. Las líneas de debajo bajan para dejar sitio,
y se vuelve a distribuir cuando la ventana cambia de ancho. Antes era una sola
línea que se salía del borde de la imagen.

**Bajo la barra hay una línea que dice qué recibe la llamada en la que estás**,
con el argumento que estás escribiendo en negrita. Pon el cursor dentro de
`randsc(` y dice *randsc(**seed**, [size=1+i], [degradation=0.5+0.5i],
[kaleidoscope=1], [kaled_mirror=0], [skew=0])*; pasa una coma y la negrita se
mueve contigo. No necesita de la fórmula más que un paréntesis abierto con un
nombre delante, así que está ahí mientras la fórmula está a medio escribir, que
es cuando hace falta. Pon el cursor fuera de toda llamada y la línea queda
vacía.

Seis de las posiciones incluidas estaban escritas a la antigua y se han
reescrito con comas: `circle`, `heart`, `helloween`, `pentafrac`, `warriormask`
y `burnship`, todas en **Archivo → Abrir → examples/Malczak**. **Una posición
tuya que tenga un punto y coma en la fórmula no se cargará** hasta que la
fórmula se escriba con comas — el mensaje nombra el carácter, y la solución es
abrir la fórmula y sustituirlos.

**Dejar un lugar vacío.** Un argumento mostrado entre corchetes se puede omitir
en medio de una llamada además de al final, no escribiendo nada entre los dos
separadores: `julian(z, ,3)` da el primero y el tercero y deja que la función
diga cuál es el segundo. Los espacios no importan — `f(z, ,5)` y `f(z,,5)` son
la misma llamada — y un argumento que no está entre corchetes hay que
escribirlo, así que `poly( ,1,2)` y `sin(,)` se rechazan en lugar de
adivinarse.

Lo que significa un lugar vacío es cosa de cada función. La mayoría toma el
valor predeterminado que declara; un coeficiente de `poly` dejado vacío es un
término que no está, así que `poly(z,1, ,1)` es `z^2 + 1`; y una rama de
`ifiter` o `ifiterl` dejada vacía repite la anterior, que es la forma de dar a
una rama más de una pasada en el ciclo. `ifiter(f(z), , ,g(z), , , , )` ejecuta
`f` durante tres pasadas y `g` durante cinco, y no cuesta más evaluarlo que
escribirlas enteras: adónde lleva cada elección se decide una vez, cuando se
analiza la fórmula.

Una fórmula de usuario que aún no se ha escrito dice `z^2+c`. Antes decía el
burning ship, que es un fractal propio y un enigma como punto de partida.

**`ifiterf(rest, last)`** — evalúa `rest` en cada pasada salvo la última y
`last` en esa. La última pasada es la última que permite el límite de
iteraciones: una fórmula no tiene forma de saber qué pasada será la que escape,
ya que depende del valor que aún no ha producido. Solo se evalúa la elegida.

**`ifiterr(before, after, n)`** — evalúa `before` mientras el número de pasada
está por debajo de `n` y `after` de `n` en adelante. Solo se evalúa la elegida,
como con las otras dos, aunque primero hubo que enseñar el umbral al parser:
ahora un argumento se puede marcar como leído por el selector en lugar de
elegido por él, y entonces se evalúa antes de hacer la elección.

**Las variables.** `z` es el valor en curso; `c` es el punto, que es el píxel
en modo mandelbrot y la constante en modo julia; `x` es la coordenada simple,
que es el píxel en ambos modos y no cambia de una pasada a otra; `n` es el
número de iteración. En la *Inicialización de usuario* `z` es donde habría
empezado `z` si no hubiera inicialización, que es el mismo valor que `x`, así
que una inicialización hecha solo de `z` dice exactamente lo mismo que no decir
nada.

Esa última parte es nueva, y también lo es poder decirla: `z` no estaba
registrado en ninguno de los dos parsers por los que pasa la inicialización, y
`x` en el equivocado, así que una inicialización que nombrara cualquiera de los
dos era rechazada por el hilo que tenía que calcularla -- mientras que el
diálogo, cuyo parser tiene ambos, la aceptaba sin dar ninguna señal de que
luego se descartaba. `n` se quedaba con lo que hubiera dejado el píxel
anterior, es decir, una inicialización que lo leía leía el orden en que se
calculaban los píxeles.

**`p1`, `p2` ... `p9999`** — el valor que tenía `z` en una pasada anterior:
`p1` en la pasada antes de esta, `p2` en la de antes, `p9999` nueve mil
novecientas noventa y nueve pasadas atrás. `p` es otro nombre para `p1`. Antes
de que haya habido tantas pasadas valen aquello de lo que parte el historial,
que decide *Fractal → Fijar los valores p en la primera iteración*: el punto
que se itera cuando está activado, cero cuando no.

Antes se detenían en `p6` y se guardaban desplazando un array de seis
posiciones un lugar en cada pasada — una copia por posición y pasada, pagada
tanto si la fórmula nombraba alguna como si no. Por eso se detenían en seis:
con nueve mil habrían sido 320 KB de copias por pasada, medidos en veinticuatro
veces el coste total de iterar `z^2+c`. Ahora el historial es un anillo,
escrito una vez por pasada, y solo se leen de él las posiciones que una fórmula
nombra de verdad, así que `p9999` cuesta lo que `p1` y una fórmula que no
nombra ninguna no cuesta nada. La ventana de referencia las enumera entre las
variables.

**Sufijos en una variable.** `z`, `c`, `x`, `p` y de `p1` a `p9999` admiten
sufijos que equivalen a las llamadas que más se hacen sobre ellas:

| sufijo | equivale a |
| --- | --- |
| `_b` | `bship(…)` |
| `_bi` | `bshipi(…)` |
| `_br` | `bshipr(…)` |
| `_pM` | `parchment(…, M)` |
| `_paM` | `parchmenta(…, M)` |

Se leen de izquierda a derecha, y cada uno envuelve lo que hicieron los
anteriores: `z_p3` es `parchment(z,3)`, `c_b_p2` es `parchment(bship(c),2)` y
`p12_p2_p3` es `parchment(parchment(p12,2),3)`. *M* es parte del nombre, así
que se escribe en cifras — un número entero, de uno o más; `z_p0` se rechaza
con un mensaje que lo dice. No en `n`, que es una cuenta y no tiene una segunda
componente sobre la que trabajar.

Se expanden en las llamadas antes de leer la fórmula, así que una variable con
sufijo dibuja exactamente lo que dibujan las llamadas escritas enteras y cuesta
lo que cuestan; la fórmula se guarda tal como se escribió. Todo lo que no sea
uno de los cinco se deja como está y se rechaza como el nombre desconocido que
es.

**Eliminados.** `powi`, `powdc` y `logcn` eran segundos nombres de `pow` y
`logn`; `rad`, `deg` y `sign` estaban en la lista sin ninguna implementación
detrás, y un segundo `trunc` quedaba oculto por el que funcionaba. Una fórmula
que use un alias eliminado debe usar el nombre que queda — la función es
idéntica.

**`erf(z)`** — la función error en el plano complejo. Precisa a unos tres ulp
dentro de un bailout de dos, donde de verdad itera un fractal.

**`poly(z, k1, k2, ..., km)`** — un polinomio en `z`:

    k1*z^(m-1) + k2*z^(m-2) + ... + k(m-1)*z + km

El primer coeficiente escrito multiplica la potencia más alta y el último va
solo, así que la llamada se lee en el orden en que se dice el polinomio. Se
calcula con la regla de Horner, que son m−1 multiplicaciones en lugar de las
m(m−1)/2 que costaría elevar cada potencia por separado, y además es la más
precisa de las dos.

    poly(z,1,0,0)+c            el Mandelbrot, escrito entero
    poly(z,1,0,0,{0.7,0.2})    z^3 + 0.7+0.2i

**`randsc(seed, size, degradation, kaleidoscope, kaled_mirror, skew, skew_mode, selfsim)`**
— ruido coherente sobre el punto, que da manchas en lugar de nieve píxel a
píxel. `size` (predeterminado `1+i`) es el ancho medio de una mancha a lo largo
del eje real y su alto a lo largo del imaginario. Es **ruido de gradiente**,
una octava del movimiento que suman los [modos
fbm](#movimiento-browniano-fraccionario) y `fbm()`. Era ruido de valor, y por
la misma razón que el de ellos dibujaba su retícula: las manchas se aplanaban a
lo largo de cada una de sus líneas y salían cuadradas. Las manchas son tan
anchas como antes — la retícula está a tres quintos del tamaño, y el campo deja
de parecerse a sí mismo la mitad a once dieciseisavos del tamaño, antes y
después — y los valores se reparten como antes, 0.2145 desde el centro, con la
curva de contraste del fbm llevando el reparto más estrecho del ruido de
gradiente de vuelta al del ruido de valor. Pero cada valor es distinto, y
cuesta 167 ns por llamada frente a 147. Los cinco mosaicos son su retícula y no
cambian ni un bit. `degradation` (predeterminado `0.5+0.5i`, que las reduce a
la mitad en cada pasada) las encoge a medida que avanza la iteración: el tamaño
se multiplica por él en cada pasada, componente a componente, así que
`0.5+0.2i` sobre `1+i` da `1+i` en la primera pasada, luego `0.5+0.2i`, luego
`0.25+0.04i`. Un cero en cualquier componente de cualquiera de los dos
argumentos devuelve cero en lugar de dividir por cero. Luego vienen los dos
argumentos del caleidoscopio, `skew` (predeterminado `0`), `skew_mode` y
`selfsim`, cada uno descrito más abajo. Solo la semilla es obligatoria.

### La inclinación, y por qué estos campos son difíciles de colorear

El motor colorea con las dos componentes de la órbita: `real` y `zmag` leen
una, `imag`, `angle` y `real / imag` la otra. Estos campos devuelven un solo
número para toda una celda, sobre el eje real, así que cada modo dibuja **un
tono por celda** — `zmag` e `iter + real` parecen el valor truncado, y los
otros tres no tienen nada que leer. Medido sobre novecientos píxeles: de 14 a
126 valores para los dos primeros, y exactamente **uno** para el resto.

Variar el valor dentro de una celda no se puede hacer con un solo número real.
El mismo número decide si el punto escapa, así que un valor que se mueve a
través de una celda saca media celda y deja dentro la otra mitad, y la celda
sale **cortada en dos por una línea recta**. Se intentó dos veces y cortó los
mosaicos las dos veces.

`skew` lo hace en cambio con la segunda componente, y **gira** el valor en
lugar de escalarlo. Cuánto gira sigue a **lo lejos del centro de su celda que
está el punto**: con `out` yendo de cero en el centro a uno en el borde, el
valor se gira el doble del arcotangente de `skew_re*out + skew_im`.

`out` se mide **en la geometría propia de la celda**, así que lo que el color
dibuja dentro de una celda es la forma en que está dividido el campo, no un
conjunto de líneas tendidas encima — cero en el centro de una celda, uno a lo
largo de todo su borde:

| campo | lo que dibuja el color dentro de una celda |
| --- | --- |
| `randsc` | sus propias manchas — el giro sigue al nivel, así que los contornos son los del campo |
| `randscq` | cuadrados, alrededor del centro del cuadrado |
| `randsch` | hexágonos, alrededor del centro del hexágono |
| `randsct` | triángulos, alrededor del centro del triángulo |
| `randscp` | su propio polígono, encogido paso a paso — no hay dos celdas con la misma forma |

La **parte real** de la inclinación es ese gradiente a través de una celda. La
**parte imaginaria** es un giro plano que comparte toda la celda: desplaza el
color de una celda sin dibujar nada dentro, así que `0.4+0.2i` da a la vez los
contornos de la forma y un desplazamiento por celda.

### Lo que dibuja la inclinación: `skew_mode`

El séptimo argumento dice **qué varía dentro de la celda**. Los números **se
suman** en lugar de excluirse entre sí, así que una llamada puede pedir dos y
obtener ambos.

| | |
| --- | --- |
| **1** (y `0`) | **la forma** — el propio contorno de la celda, encogido paso a paso. Lo que recibe una llamada que no da `skew_mode`, y lo que se usa cuando no se nombra ni 1 ni 2 |
| **2** | **la roseta** — el ángulo alrededor del centro de la celda, plegado en tantas vueltas como cuñas tiene el caleidoscopio, así que cada celda lleva la simetría de la propia imagen. Simples radios cuando no hay nada plegado |
| **4** | **a lo ancho de la cuña** — un giro que sigue dónde está el punto en su cuña del caleidoscopio, igual en cada cuña |
| **8** | **radial** — el módulo se mueve además del ángulo |

`3` es la forma y la roseta a la vez, que forma una espiral. `11` es esas dos
con el módulo moviéndose también.

**`4` sigue la geometría propia del caleidoscopio.** Gira el valor nada en los
bordes de cada cuña y todo el giro que da por sí sola la parte imaginaria — el
doble de su arcotangente, 127 grados para `2` — a lo largo de la línea por su
centro, y entre medias el giro de la parte imaginaria por lo que se ha
recorrido. Cada cuña recibe lo mismo, así que un caleidoscopio de seis sigue
siendo un caleidoscopio de seis, espejos incluidos, y el giro es continuo a
través de los bordes y los centros, que es donde refleja el pliegue. No hace
nada cuando no hay nada plegado, o cuando la parte imaginaria de la inclinación
es cero. Sin trigonometría, no cuesta nada con 64 bits de mantisa y una décima
de llamada con 113.

Al principio era **un giro distinto en cada cuña** — la parte imaginaria por el
número de la cuña más uno — pensado para que las copias dejaran de ser
idénticas. Es lo único que no pueden ser las copias de un caleidoscopio: una
imagen de simetría seis salía como seis porciones sin relación con un corte
recto en cada unión. El múltiplo además chocaba con la media vuelta, así que
cuatro porciones de seis quedaban casi sobre el eje real, y allí `real / imag`
las dibujaba como nieve por mucho que `selfsim` promediara el ruido.

Lo que queda de aquello es una propiedad de `real / imag`, que divide por la
parte imaginaria de *z*: dondequiera que un giro la lleve a través de cero el
coloreado salta de un extremo de la paleta al otro, y los últimos restos de las
pasadas finas deciden el lado píxel a píxel. Un giro que barre la cuña cruza
ahí en la mayoría de los ajustes. En un `randsct` de simetría seis coloreado
con `real / imag`, con `selfsim` a uno, una parte imaginaria de `0.25` no dejó
grano alguno y `2` dejó bandas junto a los bordes; un coloreado que no divide
por una componente no tiene ninguno.

**`8` es el que responde a `zmag`.** Un giro deja el módulo donde está, y
`zmag` y el bailout leen el módulo y nada más — así que ningún otro modo puede
alcanzarlos. El precio es el que la inclinación se diseñó para evitar: el
escape se mueve con el módulo, así que la figura se corta dondequiera que el
bailout caiga dentro de una celda. Combínalo con un bailout circular y un
coloreado que lea las componentes por separado, o acepta el corte a propósito.

Una inclinación nula deja a cada modo sin efecto, así que `skew_mode` por sí
solo no cambia ninguna imagen.

Primero se escribió una rampa recta — `skew_re*du + skew_im*dv` a través de la
celda — y es lo que esto sustituye. Coloreaba, pero dibujaba las mismas bandas
diagonales, verticales u horizontales en cada celda por igual, fuera cual fuera
la división del campo, y las bandas estropeaban la imagen.

En **`0`**, que es lo que recibe una llamada que no lo nombra, el giro es nulo,
la parte imaginaria se queda en cero, y cada valor es el número que siempre fue
**hasta el último bit** — las sumas de comprobación de antes de que existiera
el argumento pasan sin cambios, y una posición guardada se renderiza idéntica.
La posición `talc`, tres llamadas de ruido y una `c` sobre 345600 valores, da
la misma firma que la compilación anterior al argumento.

**Cuánto color da.** El índice del motor es
`(iter + quantity) * speed + shift`, así que una variación de uno en la
magnitud vale una banda a una velocidad de uno. La variación de `imag` dentro
de una celda:

| inclinación | 0.02 | 0.05 | 0.1 | 0.2 | 0.4 | 0.6 | 1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `randsc`, `randsct` | 0.08 | 0.21 | 0.43 | 0.83 | 1.5 | 1.9 | 2.1 |
| los tres mosaicos | 0.06 | 0.14 | 0.28 | 0.53 | 0.96 | 1.2 | 1.4 |

Así que **de `0.3` a `0.6` es donde aparece una banda a una velocidad de uno**,
que es lo que les faltaba a estos campos. Contar valores distintos diría que
`0.01` basta — novecientos — pero mil valores dentro de una centésima de banda
son todos un mismo color.

**Lo que cuesta.** Escalar movía el módulo del valor, y un bailout redondo mira
exactamente eso: una celda cuyo nivel estaba cerca de donde cae el bailout
salía cortada en dos, y los mosaicos perdían sus formas. Un giro deja el módulo
donde está, así que **con un bailout circular la inclinación es gratis a
cualquier intensidad** — sobre 90000 píxeles, a 0.05, 0.3 y 1, y en los cinco
campos: ni un píxel escapa de forma distinta.

Un polígono de bailout es otra cosa, porque no lee el módulo sino las
componentes, y un giro lleva el valor por un círculo que puede cruzar un lado.
Cuánto puede cruzarlo es cuánto sobresalen las esquinas de la forma más allá de
sus lados, y las cifras lo siguen — la fracción de la imagen cuyo escape
cambia, en los cinco campos:

| forma del bailout | esquinas entre apotema | inclinación 0.01 | inclinación 0.05 | inclinación 0.3 | inclinación 1 |
| --- | --- | --- | --- | --- | --- |
| círculo | 1.00 | **0%** | **0%** | **0%** | **0%** |
| hexágono | 1.15 | **0%** | 0.0–0.8% | 0.7–2.2% | 0.4–2.4% |
| cuadrado | 1.41 | **0%** | 0.0–1.0% | 0.8–2.8% | 1.2–12.8% |
| triángulo | 2.00 | 0.05–1.3% | 0.7–5.0% | 2.6–6.4% | 7.4–21.6% |

El **número** al que está fijado el bailout también lo mueve, y por la misma
razón: lo que puede cambiar de lado es un valor cuyo módulo cae entre la
apotema del polígono y sus esquinas, y dónde queda ese anillo entre las celdas
es lo que decide el número.

Así que con un polígono los dos deseos tiran en sentidos opuestos, y no hay
forma de arreglarlo de otro modo: el escape lee las dos componentes del valor y
también lo hace el coloreado, así que cualquier cosa que dé al color algo que
leer es algo que el escape puede ver. **Combina la inclinación con un bailout
circular** y no cuesta nada a ninguna intensidad. Con un polígono, o aceptas
esos pocos por ciento o mantienes la inclinación en `0.05` y subes en su lugar
la velocidad de coloreado.

Lo que un giro no puede tocar es **`zmag`**, que lee el módulo. Colorea con
algo que lea las dos componentes por separado: `iter + real`, `iter + imag`,
`angle`, `real / imag` en el exterior, y `real` o `real / imag` en el coloreado
interior.

De ahí se siguen dos cosas. El valor es complejo mientras la inclinación no sea
cero, así que `randsc(7,,,,,0.4)*z` gira `z` además de escalarlo. Y el
caleidoscopio pliega igualmente: medido sobre cinco campos, dos, tres, cinco y
seis cuñas y ambos espejos, un giro de una cuña deja cada valor donde estaba.

**`randscq(...)`** — el mismo campo sin la interpolación: un mosaico de celdas
cuadradas planas en lugar de manchas. Mismos argumentos, mismo significado.

**`randscp(...)`** — de nuevo el mismo campo, quitadas las curvas pero no la
irregularidad: se esparce una semilla dentro de cada celda y cada posición toma
el valor de la semilla más cercana, lo que dibuja un diagrama de Voronoi —
polígonos convexos planos de lados rectos, no hay dos con la misma forma.

**`randsch(...)`** — otra vez lo mismo sobre hexágonos: un panal de celdas
planas. De los polígonos regulares que teselan el plano es el que no tiene una
dirección preferente — cada celda tiene los mismos seis vecinos en los mismos
seis ángulos — así que se lee como un material y no como una rejilla.

**`randsct(...)`** — triángulos equiláteros, de orientación alterna. Un mosaico
triangular sí tiene una dirección preferente, y eso es lo que se pide al
elegirlo.

**`randsctile(tiling, seed, ...)`** — el mismo campo sobre cualquiera de
**cuarenta y cinco teselados**, que elige el primer argumento; todo lo que
viene después es lo que recibe la familia, un lugar más allá, así que `selfsim`
es el noveno. Cada teselado está escalado a una tesela de área unidad de media,
así que el primer argumento cambia la forma de las celdas y no su escala. Un
número fuera de 1 a 45 dibuja cero.

| | |
| --- | --- |
| 1–3 | los teselados regulares: cuadrados, triángulos, hexágonos |
| 4–11 | los ocho arquimedianos, polígonos regulares con todos los vértices iguales: 4.8.8 (octógonos y cuadrados), 3.6.3.6, 3.4.6.4, 3.12.12, 4.6.12, 3.3.3.4.4, 3.3.4.3.4 (cuadrado romo), 3.3.3.3.6 (hexagonal romo) |
| 12–19 | sus duales, que no son regulares: tetrakis cuadrado, rombille (los cubos), deltoidal trihexagonal (cometas), triakis triangular, kisrombille, y los pentágonos prismáticos, de El Cairo y en flor |
| 20–32 | ladrillos, aparejo flamenco, espiga, tejido de cesta, pitagórico (cuadrados de dos tamaños), chevrones, cuadrados y rombos, casitas (un pentágono, derecho y colgando por turnos), filas de cuadrados y triángulos en dos ritmos, hexágonos entre triángulos, cruces griegas, tetrominós en T |
| 33–37 | estrellas islámicas: de ocho puntas con cruces, de seis puntas con hexágonos, de ocho puntas a partir de 4.8.8, de doce puntas a partir de 3.12.12 y de 4.6.12 |
| 38 | celdas de Voronoi, irregulares — las de `randscp`, con un valor propio |
| 39–43 | teselados que nunca se repiten, de las multirrejillas de de Bruijn: los rombos de Penrose, las cometas y dardos de Penrose, los cuadrados y rombos de Ammann–Beenker, y rombos en doce y en siete direcciones |
| 44–45 | por sustitución: el molinete, cuyos triángulos miran en todas las direcciones posibles, y la silla, eles cortadas en eles |

La **pestaña Teselados** de Ayuda → Referencia de fórmulas de usuario tiene una
imagen de cada uno, con el rótulo con que la lista la pestaña Valores, y un
doble clic en una copia el comienzo de su llamada — `randsctile(18, ` — al
portapapeles, lista para la semilla. Las imágenes las dibuja el propio
`randsctile`, ocho unidades de lado alrededor del origen, mediante
`tools/randsctile-thumbnails.cpp`, y se guardan en el ejecutable en lugar de
dibujarse cada vez que se abre la pestaña. Las imágenes guardadas pueden
quedarse desfasadas, así que cada una se guarda con una huella de lo que
dibujaba la función, y un test vuelve a dibujar los mismos puntos y falla,
diciendo qué teselado, hasta que
`cmake --build <build directory> --target randsctile-thumbnails` lo haya vuelto
a dibujar.

Los periódicos son **tablas**, escritas por `tools/randsctile-tables.py`, que
construye cada teselado a partir de su geometría y lo comprueba — veinte mil
puntos al azar, cada uno de los cuales debe caer en exactamente una tesela —
antes de escribirlo. Las estrellas islámicas se hacen a partir de los teselados
arquimedianos: cada polígono de ocho lados o más se convierte en una estrella
con las puntas en los vértices del polígono, y el triángulo fino entre cada
muesca de la estrella y el lado pasa a la tesela del otro lado, así que dos
estrellas que comparten un lado forman un rombo entre ellas y un polígono
pequeño recibe una muesca en cada lado que comparte con una estrella. Las
teselas pueden ser cóncavas — estrellas, cruces, la T — así que un punto se
coloca en su tesela por el número de cruces en lugar de lado a lado. El punto
se divide en periodos enteros en la precisión de trabajo, como `randscq` lo
divide en celdas, así que un teselado periódico es exacto por lejos que esté.

Los que nunca se repiten **no tienen tabla** y se calculan de nuevo en cada
punto. Los rombos de Penrose, Ammann–Beenker y los rombos de doce y de siete
direcciones salen de la **multirrejilla de de Bruijn**: familias de líneas
paralelas en un segundo plano, una tesela dondequiera que se crucen dos,
nombrada por las dos líneas — cuatro números enteros, exactos los pida quien
los pida. El punto se coloca probando los cruces junto a donde cae en ese
plano, primero las familias cuyas líneas pasan más cerca, lo que lo encuentra
en cuatro a ocho intentos de media. Las cometas y dardos salen de los
triángulos de Robinson cortados nivel a nivel de una rueda de diez alrededor
del origen, el molinete de un triángulo rectángulo cortado en cinco, la silla
de una L cortada en cuatro; en cada nivel una línea o dos deciden en qué hijo
está el punto.

**Lo que cuesta**, medido por llamada frente a los 220 ns de `randscq` en esta
máquina: los teselados periódicos de 255 a 380 ns, los que nunca se repiten de
410 a 680. Se calculan en `double` en ambas compilaciones, incluido el lugar
dentro de un periodo — en long double los que nunca se repiten tardaban el
triple — y así las dos compilaciones dibujaron los cuarenta y cinco teselados
iguales píxel a píxel, cuarenta mil puntos de cada uno. Más allá de 2³² celdas
del origen, donde `double` empezaría a perder sus teselas, los teselados que
nunca se repiten se vuelven planos como hace la familia más allá de su rejilla;
los periódicos aguantan hasta donde aguanta el resto de la familia.

La inclinación mide una tesela como `randscp` mide su polígono, respecto al
radio del mayor círculo que cabe en la tesela, así que los contornos son la
tesela encogida, cóncavas incluidas; la roseta gira alrededor del centro de la
tesela.

    z+randsctile(18,13,{0.3,0.3},{1,1})*0.25        pentágonos de El Cairo
    z+randsctile(36,13,{0.3,0.3},{1,1})*0.25        estrellas de doce puntas
    z+randsctile(39,13,{0.3,0.3},{1,1},5)*0.25      Penrose, plegado en cinco

**`fbm(value, seed, [intensity=4], [frequency=8], [octaves=4], [roughness=0.5])`**
— un movimiento browniano fraccionario: ruido de gradiente sumado en octavas,
cada una al doble de frecuencia que la anterior y conservando `roughness` de su
altura — el mismo movimiento que dibujan los modos de coloreado de ese nombre,
ver [Movimiento browniano fraccionario](#movimiento-browniano-fraccionario). Lo
que dibuja es **desgaste** más que un motivo — manchas y abolladuras — porque
ninguna octava es lo bastante grande para verse sola y ninguna lo bastante
pequeña para desaparecer.

Donde `randsc` lee la posición y solo puede leer la posición, **esta lee lo que
escribas delante**. Ese es todo el sentido de tenerla como función:

    fbm(z,7)          se mueve con la órbita
    fbm(x,7)          se queda quieta en el plano
    fbm(z*3+c,7)      lo que sea eso
    z^2+c+fbm(z,7,0.1)*i    el conjunto con un borde rugoso

El valor y la semilla son obligatorios; el resto tiene valores predeterminados.
El movimiento va **de cero a `intensity` y nunca por debajo**, así que sumarlo
a algo no puede llevar eso por debajo de cero. La semilla se lee como la lee
`randsc`, así que el mismo número significa el mismo campo en ambas. `octaves`
se mantiene entre uno y veinticuatro, y solo merece la pena subirlo cuando se
sube `roughness`: con la mitad simple la octava octava aporta una doscientos
cincuentava parte del total.

**No tiene caleidoscopio propio**, a diferencia de la familia `randsc`, porque
ya hay uno con el que componerla: `parchment(z, sectors)` pliega el ángulo de
`z` en `|sectors|` sectores y `parchmenta` refleja las mitades alternas, así
que `fbm(parchmenta(z,6),7)` lo dice, y lo dice donde quien lee la fórmula
puede verlo.

Los modos de coloreado de **Fractal → Coloreado exterior → Movimiento browniano
fraccionario...** aplican el mismo movimiento a un color. Esta lo pone donde
una fórmula puede usarlo.

| | |
| --- | --- |
| `randsc` | suave, curvo |
| `randscq` | duro, regular, cuadrado |
| `randscp` | duro, irregular |
| `randsch` | duro, regular, hexagonal |
| `randsct` | duro, regular, triangular |
| `randsctile` | duro, cualquiera de cuarenta y cinco teselados |

Las seis colocan una celda sobre cada cuadrado unidad del tamaño vigente, así
que `size` significa lo mismo en todas: cambiar una letra cambia la forma de
las celdas y no la escala de la imagen. Un hexágono de circunradio uno cubre
2.6 cuadrados unidad y un triángulo de lado uno cubre 0.43, así que esas dos
rejillas se escalan para coincidir con el resto en lugar de dejarse como
vienen.

Cada pasada tiene un campo propio. La iteración entra en el hash junto con la
posición, así que una fórmula que llama a una de estas una vez por iteración
obtiene un valor nuevo cada vez — incluso con una degradación de uno, que deja
el tamaño como está. La degradación fija lo grandes que son las manchas en una
pasada dada; no es lo que hace distinta la pasada.

Hay un suelo bajo el tamaño. Una celda se encuentra dividiendo la posición por
el tamaño y el resultado tiene que caber en un entero, así que el tamaño no
puede bajar útilmente de más o menos la posición dividida por lo que cabe en un
entero — una milmillonésima de milmillonésima de ella. Una degradación de un
medio llega ahí en unas sesenta pasadas. Más allá no queda estructura de celdas
que resolver y el campo es un único valor plano en todo el plano; sigue
cambiando en cada pasada y sigue difiriendo entre las seis funciones, así que
una fórmula que reste una de otra no se queda en cero iterando hasta el límite
para nada. Para un desvanecimiento que siga siendo una imagen hasta el final,
usa una degradación cercana a uno.

Todas hacen el hash de la posición del punto y de la iteración, nunca de `z`, y
nunca de ningún estado global. Así que la misma imagen vuelve en cada
redibujado, con cualquier número de hilos, y en modo Mandelbrot o Julia por
igual. `randsc` es continuo, así que las dos precisiones coinciden hasta unos
1e-19; los mosaicos son funciones escalonadas y pueden discrepar en un hilo a
lo largo de los bordes de las celdas, lo que es inherente a pedir bordes duros.

Una fórmula que llame a cualquiera de ellas desactiva el seguimiento del borde.
Esa optimización rellena sin calcularla una región que encontró rodeada de un
solo color — cierto para un fractal, falso para un campo de ruido.

El antiguo `rand` no ha cambiado y sigue dependiendo de cuántas veces se ha
llamado, así que la misma posición no se redibuja igual. Usa mejor `randsc`.

### Caleidoscopios

Los dos argumentos que siguen a la degradación pliegan el plano antes de que el
campo se muestree en él. El primero es en cuántas cuñas plegarlo, y `1` — el
predeterminado — lo deja como está; el segundo, `kaled_mirror`, es qué espejo
hace el plegado:

| kaled_mirror | |
| --- | --- |
| `0` (y cualquier otro valor) | la mitad lejana de cada cuña refleja la mitad cercana, así que cada cuña es simétrica respecto a su propia bisectriz |
| `1` | lo mismo al revés, con la mitad cercana reflejando la lejana |

    randsc(13,{0.6,0.6},{1,1},6,0)     seis cuñas, cada una espejo de sí misma
    randsc(13,{0.6,0.6},{1,1},3,1)     tres cuñas, plegadas en el otro sentido

Ambos plegados son continuos donde se juntan las cuñas, así que el ruido sigue
siendo coherente y las dos precisiones siguen coincidiendo; un plegado que se
encontrara consigo mismo de forma desigual se vería como una costura.

Es la única parte de la familia que cuesta trigonometría, y solo un valor de
dos o más llega a ella. Una llamada que no dice nada al respecto paga una
comparación.

### Sacar una imagen del ruido

`randsc` devuelve un valor en `[0, 1)`. Usado solo, la iteración nunca sale de
un bailout de 4, así que nada escapa y la imagen es plana. Multiplícalo, o deja
que perturbe una iteración que sí escapa:

    z+randsc(13,{0.35,0.35},{1,1})*0.25       manchas
    randsc(13,{0.15,0.15},{0.5,0.5})*2.5      movimiento browniano
    randscq(13,{0.15,0.15},{0.5,0.5})*2.5     cuadrados dispersos, cada vez menores
    z+randscp(13,{0.3,0.3},{1,1})*0.25        polígonos irregulares
    z+randsch(13,{0.3,0.3},{1,1})*0.25        un panal
    z+randsct(13,{0.3,0.3},{1,1})*0.25        triángulos
    z^2+c+randsc(13,{0.25,0.25},{1,1})*1.2    el propio conjunto deformado

Una degradación de 0.5 reduce las manchas a la mitad en cada pasada, así que
tras veinte iteraciones son una millonésima de su tamaño y quedan por debajo de
un píxel. Para un fundido lento a lo largo de muchas iteraciones, usa algo
cercano a 1 — `0.97^60` es todavía 0.16.

Una semilla entera es exacta en ambas compilaciones; una fraccionaria se
cuantiza a sus bits iniciales, que coinciden salvo en una semilla de cada
cincuenta millones.

### Todas las pasadas hasta ahora: `selfsim`

Cada pasada de estas funciones es un campo propio — la iteración entra en el
hash — con celdas iguales a la degradación por las de la pasada anterior. Una
fórmula que llama a una en cada pasada lee por tanto cada vez un campo nuevo y
más fino, y un modo de coloreado que lee la última pasada lee solo el último
campo. Con una degradación de 0.75 las celdas de una vista de cinco de ancho
están por debajo de un píxel hacia la decimosexta pasada, y más allá el
interior es nieve: cien pasadas de `randsct(672,,{0.75,0.75},,,{2,2},6)+c` en
modo Julia dejan al 70% de los píxeles interiores vecinos a un dieciseisavo de
la paleta o más unos de otros.

El octavo argumento devuelve en cambio **las pasadas hasta ahora,
promediadas**, que guardan entre sí la relación de las octavas de un movimiento
browniano fraccionario: la pasada *n* pesa d^(nH), donde *d* es la degradación
— la media geométrica de sus dos componentes tomadas sin signo, así que
`0.75+0.75i` es 0.75 — y *H* es el argumento. Los pesos se dividen por su
total, así que la respuesta queda entre los valores promediados: en [0, 1] sin
inclinación.

| selfsim | |
| --- | --- |
| omitido | apagado — la única pasada, hasta el último bit. Lo que recibe una llamada que no nombra nada, y un lugar dejado vacío |
| `0` | la media simple: cada pasada pesa igual, y es lo que da todo *H* cercano a cero — no apagado |
| `1` | el movimiento simple: cada pasada pesa *d* veces la anterior |
| `0.5` | más áspero: las pasadas finas conservan más peso, así que se ve más detalle |
| mayor | más suave: las primeras pasadas se llevan casi todo |
| negativo | las pasadas finas pesan más que las gruesas, y vuelve la nieve |
| complejo, `{1,2}` | la parte real pesa como arriba, la imaginaria gira las pasadas — ver más abajo |

    randsct(672,,{0.75,0.75},,,{2,2},6,1)+c     lo mismo, autosemejante

En esa fórmula la nieve desaparece: los píxeles interiores vecinos quedan de
media a 1.8 celdas de paleta en lugar de 35, y el 0.6% de ellos a más de
dieciséis en lugar del 70%. Lo que queda son los triángulos a todas las escalas
a la vez, los grandes llevando a los pequeños.

* **`selfsim` puede ser complejo.** Entonces d^(nH) también lo es: d^(n·Hr) en
  tamaño y un giro de n·Hi·ln d, así que cada pasada gira eso más que la
  anterior y luego se promedia como promedia un *H* real — una espiral de
  octavas. A 0.75, una parte imaginaria de uno son unos dieciséis grados por
  pasada. Los pesos se dividen por **la suma de sus tamaños**, no por su suma
  compleja: eso mantiene la respuesta dentro del mayor de los valores
  promediados. La suma compleja se acercaría a cero en alguna pasada y haría
  estallar todos los valores en esa pasada para todos los píxeles a la vez, ya
  que los pesos no tienen nada que ver con el píxel. Los giros dan al valor una
  parte imaginaria incluso sin inclinación, así que `imag` y `angle` tienen
  algo que leer. Una degradación de uno tiene logaritmo cero y no gira nada;
  `{1,0}` es `1` hasta el último bit; y `{0,2}` es cada pasada pesada igual, y
  girada.
* **Cero no es apagado.** d^(n·0) vale uno en cada pasada, así que `0` es la
  media simple, y también lo es todo *H* cercano — por ambos lados, y también
  complejo. Al principio el cero era el interruptor, lo que lo convertía en el
  único valor por el que la curva no pasa: `0.000001` era la media simple y `0`
  solo la última pasada, a un cuarto del rango de distancia. Solo la última
  pasada es el extremo lejano de la curva, con *H* tendiendo a menos infinito.
  Para no tener media, omite el argumento, o deja su lugar vacío.
* **Lo que se promedia es el valor inclinado**, así que la inclinación y cada
  bit de `skew_mode` siguen haciendo lo que hacían, y la respuesta es compleja
  cuando la inclinación la hace compleja.
* **Una media junta los valores.** A 0.75 y con *H* de uno la dispersión es
  unos dos quintos de la de una pasada, así que sube la velocidad de coloreado
  para recuperar el contraste.
* **Una degradación de uno** pesa igual cada pasada y da su media simple, que
  se asienta hacia un centro plano a medida que se suman las pasadas. Es lo que
  dice la definición, y se hace como se dice.
* **Las pasadas que una llamada no vio se calculan cuando se le pregunta.** Una
  llamada en una rama que la fórmula solo toma desde la séptima pasada recibe
  igualmente las siete primeras, y la respuesta en una pasada es la misma sea
  cual sea el camino que llegó a ella. Un píxel nuevo se reconoce por su
  posición además de por su pasada: solo la pasada, que es en lo que se basa
  `trap`, permitiría que tal llamada continuara la media que dejó el píxel
  anterior.
* **Nada desborda.** La media se guarda como una sola, y cada pasada toma de
  ella una parte que se deduce de la parte de la pasada anterior, en lugar de
  como una suma y un total — así que un *H* negativo, o una degradación mayor
  que uno, cuyos pesos superarían lo que cabe en el tipo, no tiene dónde
  hacerlo.
* **Cuesta** alrededor de una décima más por llamada con 64 bits de mantisa y
  una séptima con 113. Apagado, nada que pueda medirse.

### Figuras en lugar de ruido

`sierpinskyt`, `sierpinskyc` y `snowflake` están donde está el ruido, y no
tienen nada de aleatorio: la figura está donde parece estar, y el modo
Mandelbrot y el modo Julia dibujan la misma imagen, porque una forma en el
plano está en el mismo sitio en ambos casos.

**Son fractales propios, no campos que multiplicar dentro de uno.** Escrita
sola — `sierpinskyt()` y nada más — cada una dibuja su figura, del mismo modo
que Fractal → Más fórmulas dibuja su **Sierpinski**, su **Sierpinski Carpet** y
su **Koch Snowflake**.

Las tres lo hacen **llevando el punto a su padre**. Cada parte de una de estas
figuras tiene uno: un hueco del triángulo de Sierpinski está dentro de un hueco
mayor un nivel más arriba, un hueco de la alfombra dentro de la celda que se
cortó del mismo modo, un triángulo del copo de nieve sobre un lado del hexágono
del centro o sobre el borde libre de otro triángulo. El paso al padre es un
movimiento del plano — duplicar alejándose de la esquina más cercana en el
triángulo de Sierpinski, ampliar una celda por el número de celdas en la
alfombra. La parte superior no tiene padre dentro de la figura, así que su paso
la lleva fuera del bailout, y una parte *n* niveles más abajo necesita *n*
pasos para llegar allí: **la pasada en que sale un punto es el nivel en que
está**, y el recuento de iteraciones es la imagen.

Lo que vuelve en la pasada en que sale un punto importa tanto como la propia
pasada: **todos los modos de coloreado exterior salvo el recuento de
iteraciones lo leen**. Así que la parte superior de cada figura se lanza fuera
en lugar de simplemente declararla escapada — el triángulo de Sierpinski
duplica su hueco central alejándolo de una esquina y cae fuera del triángulo,
la alfombra escala su celda central respecto a la celda contigua, el copo de
nieve lanza fuera su hexágono central en la dirección hacia la que mira ese
sector. Cada uno lleva consigo de dónde venía, así que `real`, `imag`, `angle`,
el coloreado suave y los demás tienen todos algo con que trabajar. Solo a un
punto que nunca estuvo en la figura se le da un número que representa
«escapado».

Un copo de nieve se lee **desde un hexágono hacia fuera**, y eso es lo que hay
que saber de él. Se descompone exactamente en un hexágono regular en el centro,
seis triángulos del mismo lado que ese hexágono apoyados en sus seis lados,
doce de un tercio de eso sobre sus bordes libres, cuarenta y ocho de un noveno
sobre los suyos, y así para siempre. Tres de los seis son las esquinas del
triángulo del que creció la figura y tres son los primeros salientes puestos en
sus lados, y **nada los distingue**: desde el hexágono hacia fuera, un copo de
nieve tiene simetría de orden seis.

Leído en cambio desde el primer triángulo, cinco octavos de la figura son el
nivel uno — un único triángulo plano que llena la imagen, con el copo de nieve
solo en la franja de alrededor, y nada que decir para el color exterior. Leído
desde el hexágono, los dos primeros niveles son cinco doceavos cada uno y los
demás bajan de forma uniforme. Un triángulo apoyado en el hexágono pasa al
triángulo más grande que el hexágono contiene — esquina sobre esquina, tres
veces el área — así que se amplía por la raíz cuadrada de tres y gira un
doceavo de vuelta; uno apoyado en el borde libre de otro triángulo se amplía
por tres y gira hacia donde mira ese borde.

**Cada paso expande**, y eso es una regla más que una casualidad. El triángulo
de Sierpinski duplica, la alfombra amplía una celda por el número de celdas, el
copo de nieve por tres o por la raíz de tres. Un triángulo apoyado en el
hexágono podría en cambio doblarse sobre el lado en que se apoya, lo que lo
deja sobre uno de los seis de que está hecho el hexágono y es más ordenado;
pero un doblez es una isometría, y una isometría no tiene sensibilidad a dónde
empezó el punto. Escrita dentro de una fórmula mayor, la figura daría entonces
a esa fórmula una imagen rígida de sí misma en cinco doceavos de su área, y la
fórmula dibujaría allí una región plana. Una prueba verifica las tres
expansiones.

El punto viaja de verdad, y eso es lo que las hace utilizables y no solo
correctas: `z` a lo largo del camino es un punto del plano como cualquier otro,
así que los modos de coloreado que lo leen, el coloreado suave y escribir una
figura dentro de una fórmula mayor tienen todos sentido.

Un triángulo de Sierpinski y una alfombra llenan su forma, así que cada punto
de ellos tiene un nivel y sale en él. Un copo de nieve no llena su hexágono:
deja seis esquinas de **fondo**, y ese fondo no forma parte de la figura y no
tiene nivel. Se **gira media vuelta**, lo que lleva fondo a fondo, así que
nunca sale y aparece en el **color interior** — el espacio vacío sigue vacío en
lugar de quitarle a la figura una banda del color exterior.

Girado, y no devuelto donde estaba, porque **un punto fijo es una órbita
detenida y no una órbita acotada**. Quieto, el bucle ejecutaba cada pasada
hasta el límite sin que nada cambiara; los modos de coloreado interior que
comparan una pasada con la anterior no tenían nada que comparar; y
`snowflake()` escrita dentro de una fórmula mayor devolvía a esa fórmula el
punto que ya tenía, mientras que `sierpinskyt()` y `sierpinskyc()` mueven cada
punto que reciben. La media vuelta se prefiere al sexto de vuelta, que la
simetría de la figura también permitiría: ambos llevan fondo a fondo, pero la
media vuelta mantiene el punto dentro de un cuadrado y de un círculo además de
dentro del hexágono. Deja la distancia al centro intacta, así que los modos de
coloreado interior que la leen no se ven afectados, y la imagen que dibuja
`snowflake()` no cambia ni un píxel.

Es el único punto en que estas figuras cuestan algo. A un píxel de fondo se le
hace la pregunta de pertenencia una vez por pasada mientras se lo mire,
mientras que las otras dos figuras y el cuerpo del copo de nieve salen en la
pasada que nombra su nivel. Dibujado directamente con el bucle de iteración, un
copo de nieve entero cuesta unas 2.4 veces lo que costaba cuando al fondo se le
daba una banda propia. El seguimiento del borde queda activado para estas
fórmulas — sus bandas son realmente sólidas, que es la condición para usarlo —
y el fondo es una única región sólida, así que lo que llega a la pantalla es
menos de lo que sugiere esa cifra.

Escritas dentro de una fórmula mayor — `snowflake()^-2+snowflake()` y similares
— las tres **no** se comportan igual, y vale la pena decir la razón con
claridad, porque ningún trabajo sobre la figura la cambiará.

Un triángulo de Sierpinski es un **polvo**: no tiene área, sus puntos nunca
salen, y cada uno de ellos se duplica en cada pasada. Una fórmula construida
sobre él tiene algo caótico debajo en todas partes, a todas las profundidades.
Un copo de nieve es **sólido**: diez doceavos son el hexágono y los seis
triángulos, que salen en la primera y la segunda pasada, así que una fórmula
construida sobre él se queda sin figura casi enseguida y dibuja lo que dibuje
por sí misma. El detalle está a lo largo del borde de Koch, donde los niveles
son profundos, y las regiones planas son los interiores, donde no lo son.

**Lo que sabe una pasada es lo que esa pasada pagó.** El copo de nieve decide
si un punto es fondo bajando por la curva de Koch, y a esa bajada se le permite
**un nivel por pasada**. Antes recorría sus veinticuatro niveles en una sola
evaluación, así que una fórmula escrita alrededor de la figura recibía todo el
borde en su primera pasada y lo dibujaba entero por pocas pasadas que se le
dieran — con `Iteraciones: 2` aún se obtenían todos los salientes de la curva
de Koch, algo que ninguna otra fórmula de XaoS hace y que no es para lo que
sirve un número de iteraciones. Quedarse sin niveles no es lo mismo que estar
fuera: un punto fuera del marco o por encima de la envolvente es fondo y se
sabe ya que lo es, mientras que uno indeciso sigue adelante como parte de la
figura y se le vuelve a preguntar en la pasada siguiente con un nivel más que
gastar. La imagen que dibuja `snowflake()` por sí sola no cambia — un punto de
la figura nunca necesitó la bajada profunda, solo la seguridad de no ser fondo.

Lo que queda es la diferencia honesta entre las dos figuras, y ningún trabajo
sobre ninguna de ellas la cerrará. Un triángulo de Sierpinski es un **polvo**:
sin área, sus puntos nunca salen, cada uno duplicado en cada pasada, así que
una fórmula construida sobre él tiene algo caótico debajo a todas las
profundidades. Un copo de nieve es **sólido**, y diez doceavos de él son el
hexágono y los seis triángulos, que salen en la primera y la segunda pasada.
Medido con un solo zoom — la fracción de píxeles que caen en el borde de una
banda — `sierpinskyt()^-2+sierpinskyt()` llega al 17 por ciento hacia las doce
pasadas y `snowflake()^-2+snowflake()` al 2.5; ambos crecen con el número de
pasadas, que es lo importante, y dónde se detienen es lo que son las figuras.

`radius` significa lo que significa `bailout` y se lee igual: **como el
cuadrado de la distancia**. Lo que dibuja es la forma que dibuja un bailout de
ese número, **inscrita exactamente en ella** — pon el bailout en la forma
correspondiente y el mismo número, y la figura y el bailout quedan
superpuestos:

| | |
| --- | --- |
| `sierpinskyt` | el triángulo del bailout **triángulo −90°**, esquina con esquina |
| `sierpinskyc` | el cuadrado del bailout **cuadrado**, lado con lado |
| `snowflake` | el hexágono del bailout **hexágono 0°**, con sus seis puntas en las seis esquinas del hexágono |

Lo que hay que saber es que un polígono de bailout pone sus **lados** a la raíz
cuadrada del bailout desde el centro — es su apotema, no su radio circunscrito
— así que sus esquinas están más lejos de lo que dice el número: el doble en un
triángulo, y un séptimo más en un hexágono. `sierpinskyt(4)` llega por tanto a
4 con sus esquinas mientras que `sierpinskyc(4)` llega a 2 con sus lados, y
ambos son correctos.

Pon el bailout mayor que el radio y la figura ya no llena la forma: hay un
**margen** entre ambos, y ponlo menor y la figura queda recortada. Ambas cosas
son lo que se pidió y ninguna se corrige. Pero el margen no forma parte de la
figura y no tiene nivel, así que no se le da uno: se **gira según la simetría
propia de la figura** — un tercio de vuelta en un triángulo de Sierpinski, un
cuarto en una alfombra, media vuelta en un copo de nieve — lo que lo mantiene
fuera de la figura y dentro del bailout, así que nunca sale y se dibuja en el
**color interior**, con los modos de coloreado interior trabajando sobre él.
Antes salía en la primera pasada y tomaba la banda más baja del color exterior,
lo que dibujaba todo el margen en un único tono plano y le quitaba una banda a
la figura.

Un giro y no un espejo, aunque un espejo serviría igual para mantenerlo donde
está: lo que una figura devuelve tiene que girar con el punto, para que girar
la imagen gire lo que hay dibujado en ella, y un espejo lo gira en sentido
contrario. Por la misma razón el centro de una alfombra — el bloque de celdas
que encierra el anillo del borde, que no tiene padre y por tanto hay que lanzar
fuera — se lanza en la dirección hacia la que mira su propio **cuarto** del
cuadrado, los cuartos que cortan las diagonales. Lanzado siempre en la misma
dirección, como antes, dibujaba ese centro como una rampa de un lado al otro
sin nada a través; el cuadrado tiene simetría de orden cuatro y así debería
tenerla lo que se dibuja en él. Una prueba verifica que cada figura responde a
un giro de su propio orden con el mismo giro.

**`sierpinskyt([radius=4], [kaleidoscope=1], [kaled_mirror=0])`** — el
triángulo de Sierpinski, en ese triángulo.

**`sierpinskyc([radius=4], [squares=3], [kaleidoscope=1], [kaled_mirror=0])`**
— la alfombra de Sierpinski, en ese cuadrado. El cuadrado se corta en `squares`
por `squares`, el anillo de celdas a lo largo del borde se conserva, todo lo
que ese anillo encierra se desecha, y lo mismo se hace con cada celda
conservada. Así que la imagen es **un cuadrado en el centro y 4×`squares`−4
alrededor**: ocho con tres, que es la alfombra tal como suele dibujarse, doce
con cuatro, dieciséis con cinco. Dos es el único número sin un anillo
propiamente dicho, y ahí lo que se va es la esquina lejana, lo que vuelve a ser
un triángulo de Sierpinski — un cuadrado cortado en cuatro con una esquina
quitada es justo lo que es un triángulo de Sierpinski.

**`snowflake([radius=4], [kaleidoscope=1], [kaled_mirror=0])`** — el copo de
nieve de Koch, con sus puntas en las esquinas de ese hexágono, en bandas por
generación desde el centro hacia fuera: el hexágono en la primera pasada, los
seis triángulos sobre sus lados en la segunda, los doce sobre sus bordes libres
en la tercera, los cuarenta y ocho en la cuarta. Las seis esquinas de fondo que
la figura no cubre no tienen banda alguna — nunca salen, y se dibujan en el
color interior.

Cada argumento tiene un valor predeterminado, así que `snowflake()` es una
llamada, y también lo es `sierpinskyc( ,5)` — una alfombra más calada con el
tamaño predeterminado.

**El caleidoscopio es el de la familia del ruido**, escrito al final en las
tres y con el mismo significado: en cuántas cuñas se corta el plano, y qué
espejo las pliega. `1` — lo que recibe una llamada que no dice nada — no pliega
nada y deja la figura con el número que siempre fue, hasta el último bit.

Lo que pliega es distinto, sin embargo, y tiene que serlo. El ruido pliega la
**posición**, que se queda quieta durante toda la órbita; una figura no tiene
uso para la posición y lee **z**, el punto que va llevando hacia abajo hacia la
parte en que está, así que lo que se pliega es z. La imagen sale igual — la
figura dibujada en una cuña y repetida alrededor del origen — y se mantiene
hasta el fondo, ya que cada paso se da sobre el punto plegado y se vuelve a
plegar en la pasada siguiente. Dentro de la media cuña que el plegado deja
intacta, la figura es la figura.

    sierpinskyt(4,6,0)      el triángulo de Sierpinski seis veces alrededor del origen
    sierpinskyc(4,3,5,1)    la alfombra en cinco cuñas, el espejo cercano
    snowflake(4,5)          cinco cuñas de copo de nieve

Un copo de nieve ya tiene la simetría de orden seis con el espejo a lo largo de
cada bisectriz, así que plegarlo en seis cuñas de ese modo es la identidad
sobre la imagen: `snowflake(4,6,0)` sale **píxel a píxel** igual que
`snowflake(4)` en 90000 píxeles. Pide cinco cuñas, o cuatro, para verlo plegar.

Una cosa que saber sobre la alfombra: lo que recorta lo lanza en una de cuatro
direcciones, los cuartos que cortan las diagonales del cuadrado, y sobre la
propia diagonal la elección es un empate. Un plegado puede llevar todo un rayo
del plano a esa diagonal, y allí los dos lados del empate se lanzan en
direcciones opuestas — una línea finísima a lo largo de la diagonal, que es el
escalón que la alfombra siempre ha tenido y no algo que haga el plegado. Medido
sobre diez mil puntos con seis cuñas: todos los desacuerdos estaban sobre una
diagonal, y los propios puntos plegados coincidían hasta una parte en 10^16.

    sierpinskyt()        con Fractal -> Forma del escape -> triángulo -90
    sierpinskyc()        con forma del bailout cuadrado
    sierpinskyc( ,5)     una alfombra más calada
    snowflake()          con forma del bailout hexágono 0

Ninguna cuesta más que el ruido de al lado, pasada por pasada: medidas frente a
`randsc`, el triángulo de Sierpinski 0.44, la alfombra 0.54 y el copo de nieve
0.50 de su coste con 64 bits de mantisa, y 0.50, 0.87 y 0.65 con 113. Leer un
copo de nieve desde el centro lo hizo más barato además de más agradable a la
vista: seis séptimos de la figura son el hexágono o un triángulo apoyado en él,
y a esos se responde con cuatro multiplicaciones y una comparación sin recorrer
la curva en absoluto. Lo que cuesta un copo de nieve en una imagen entera es el
párrafo de arriba, y depende de cuántas pasadas toma el fondo más que de lo que
cuesta una pasada. No fue gratis, y con 113 bits casi se perdió — una raíz
cuadrada y una división son ahí ambas software, y cualquiera de ellas cuesta
más o menos lo que una figura entera. Así que el recíproco de la raíz cuadrada
del radio se guarda en el punto de llamada en lugar de calcularse en cada
pasada, y la propia raíz se recupera multiplicando el radio por ese recíproco
en lugar de dividiendo; el triángulo de Sierpinski compara sus tres pesos
baricéntricos escalados por un número positivo, lo que no cambia ni cuál es el
mayor ni cuál es negativo; la alfombra cuenta en celdas en lugar de en el
plano, donde la multiplicación por `squares` y la división por él se cancelan;
y el copo de nieve toma el hijo al que pasar del lado bajo cuyo saliente se
encontró, que tenía que encontrar de todos modos. Recorre además la curva de
Koch en `double`, ya que la figura tiene un tamaño fijo y 52 bits de una
fracción de un lado son más de los que pueden usar los 24 niveles que dibuja.
Una prueba verifica las proporciones, así que una reescritura que las pierda lo
dice.

## Paletas

Fractal → Paleta pide un número de algoritmo, y tras el número dice cómo se
llama cada uno. Del 1 al 3 son los de XaoS; del 4 al 7 se rehicieron después de
medir de qué están hechos los tres primeros, y el 8 después de medir unos
cuadros.

| | |
| --- | --- |
| 1 | **Oscuro y color** — un casi negro con un tinte propio y un color, alternados |
| 2 | **Negro, color, blanco** |
| 3 | **Cálidos sobre verde oscuro** — violetas, rojos y naranjas vivos sobre verdes oscuros y verdes azulados, negro o blanco uno de cada tres |
| 4 | **Smog** — oscura y melancólica: grises, ocre, pizarra y violeta polvoriento, con sangre de buey y escarlata entre ellos |
| 5 | **Cálidos sobre la noche** — rosa, rojo, burdeos, marrón, naranja, ocre y amarillo sobre verde azulado, petróleo, azul e índigo, en el esqueleto de la 3 |
| 6 | **Parejas preferidas** — las parejas de colores que los tres primeros ponen juntas con más frecuencia, colocadas enteras entre negro y blanco |
| 7 | **Colores aleatorios** — cada parada un color realmente al azar, y nada más |
| 8 | **Kandinsky** — los degradados de sus cuadros, de lo mate a lo vivo y de vez en cuando casi al gris |

### De qué están hechos los tres primeros

Veinte mil paletas de cada uno, generadas como las genera el programa, cada
color nombrado por el más cercano de tres docenas de colores de referencia en
CIELAB, contando los colores y las parejas de colores:

* **Anclas.** El negro y el blanco son entre un tercio y la mitad de las
  paradas, y los colores se sitúan entre ellos con un periodo de dos o tres
  paradas, así que el patrón se ve en las primeras bandas — que, en una paleta
  de unas tres mil paradas, son todo lo que usa la mayoría de las imágenes.
* **El 1 y el 2 contienen 256 colores y no más.** Toman el byte bajo del
  generador, y el byte bajo de ese generador es por sí solo un generador de
  periodo 256: el rojo decide el verde y el azul, y cada color va siempre
  seguido del mismo. En veinte mil paletas del 2, ningún color fue seguido ni
  una vez por otro distinto. De ahí vienen sus parejas recurrentes — rosa junto
  a amarillo veinticuatro veces más a menudo que por azar, burdeos junto a
  oliva seis, rojo junto a petróleo ocho.
* **El 3 toma de dos ventanas**: colores vivos del violeta pasando por el rojo
  hasta el naranja, y oscuros del verde hasta el verde azulado. Sus parejas más
  comunes son un color cálido vivo sobre un verde oscuro.

### Los cuatro siguientes

Cada uno se midió frente a los tres primeros antes de quedarse: la divergencia
de Jensen-Shannon entre los colores que muestran dos algoritmos, por nombre, y
entre las parejas de vecinos que muestran, sobre cinco mil paletas — 0 para la
misma paleta, 1 para nada en común. Con esa medida el 1 y el 2 distan 0.11, el
1 y el 3 0.16.

* **El 4, Smog**, es oscuro y apagado donde el 1 es vivo sobre negro. Las
  paradas profundas son un color apenas por encima del negro, las claras un
  gris nublado con un tinte y nunca blanco, y una de cada cuatro un gris liso.
  Sus tonos se recorren sobre un anillo de 256 como el 1 y el 2 recorren los
  suyos, así que los vecinos se repiten como los suyos, pero solo sobre los
  tonos del smog — ocre y oliva, y acero, pizarra y un violeta polvoriento — y,
  en un cuarto del anillo, del carmesí al escarlata: sangre de buey cuando es
  profundo y un escarlata apagado cuando es claro, mantenido saturado, porque
  un rojo oscuro que no lo está se lee como marrón. Una parada de cada cinco se
  lee como roja, frente a una de cada treinta en el 1. Dista de 0.43 a 0.56 de
  los tres primeros. Al principio era el mecanismo del 1 sobre otro anillo de
  256 colores, elegido lo más lejos posible del anillo clásico, y medía 0.05
  del 1 — más cerca que el 2: el parecido era el mecanismo, negro la mitad de
  las veces y un color de cualquier sitio la otra mitad, no los colores.
* **El 5, Cálidos sobre la noche**, es el esqueleto del 3 con la temperatura
  invertida: colores cálidos de toda luminosidad, del rosa pasando por rojo,
  burdeos, marrón y naranja hasta ocre y amarillo, sobre colores nocturnos
  mantenidos oscuros, del verde azulado pasando por petróleo y azul hasta
  índigo. Las dos ventanas tienen 125 grados de ancho, cerca de los 135 del 3;
  al principio eran un tercio de eso y se leían como dos colores. El azul
  marino es su color más común, y azul marino con oliva, burdeos o marrón sus
  parejas más comunes.
* **El 6, Parejas preferidas**, toma cada pareja de colores que uno de los tres
  primeros pone uno tras otro al menos tres veces más a menudo que por azar —
  sesenta y seis, con los dos colores tal como salieron de verdad — y coloca
  cada pareja entera, los dos colores tocándose, entre una parada negra y una
  blanca.
* **El 7, Colores aleatorios**, es cada parada un color sacado al azar de la
  parte alta del generador, donde el 1 y el 2 sacan de la baja: setenta y nueve
  mil colores distintos en veinte mil paletas donde el 2 tiene 256, y ninguna
  pareja más de un quinto más común que por azar. Al principio tenía negro y
  blanco una de cada tres paradas, y medía 0.00 del 2: el ritmo de muy oscuro y
  muy claro era todo lo que se veía. Sin él tiene menos recorrido de oscuro a
  claro que los anclados, que es lo que significa aleatorio.

Los cuatro que sustituyen eran un espectro, un duotono, una tríada y una pareja
complementaria, razonados a partir del círculo cromático, y en la pantalla cada
uno era un solo color. Repartían sus tonos a lo largo de toda la paleta, y en
el programa la paleta tiene unas tres mil paradas — de cuatro a nueve en la
prueba que debía detectarlo, por eso no lo detectó — así que una imagen veía
solo el primer tono. Y mantenían cada parada como un color, profundo o pálido,
sin negro y sin blanco, así que nada en ellas tenía un borde. La prueba genera
ahora sus paletas como lo hace el programa, y pide a los siete que la parte que
ve una imagen vaya de oscuro a claro y contenga más de un color; al 8, que
llegó después, se le pide algo menos, y el porqué está más abajo.

Una posición guardada con los números del 4 al 7 antes de este cambio vuelve
con otros colores: el número y la semilla son todo lo que una posición registra
de su paleta, y lo que significa el número ha cambiado. Una guardada con los
números del 1 al 3 no se ve afectada, hasta el último bit, y una guardada con
los números del 4 al 7 nombra un algoritmo que el XaoS original no tiene y
rechazará.

### 8, Kandinsky

Hecha a partir de dos series de imágenes a su manera: cuatro viñetas vivas —
negro, bermellón, amarillo, cerúleo, violeta y rosa sobre papel crema — y una
imagen mate, naranja quemado, ocre, marrón, salvia y petróleo sobre el mismo
papel. Las viñetas dan dieciocho colores, por k-means en CIELAB, y cada uno
tiene una **forma mate**: el color de la imagen mate más cercano en luminosidad
y tono, dejando fuera el croma, ya que el croma es en lo que difieren las dos.
El bermellón se vuelve naranja quemado, el amarillo ocre, el burdeos marrón, el
negro un marrón oscuro, el cerúleo salvia, el azul petróleo. La imagen mate no
tiene gris, así que el gris sigue siendo gris. El aguamarina que parece tener
es un salvia, 128,146,121: se lee como aguamarina solo junto al naranja.

Cada paleta sortea un **ánimo**. Cuatro de cada cinco caen de forma uniforme
entre las formas mates y las vivas; la quinta va más allá de las formas mates
hacia el gris, hasta quince centésimas de su croma, y una paleta de cada veinte
o así es casi gris.

**Las anclas llegan por casualidad.** Una paleta de cada diez puede convertir
alguna parada aquí y allá en negro puro o blanco puro: cada parada después de
la primera con su propia tirada, con una frecuencia que la paleta sortea entre
una parada de cada veinte y una de cada diez, negro o blanco según caiga la
tirada — sin ritmo y sin turnos, a diferencia de las anclas del 1 al 3. Las
otras nueve paletas de cada diez son exactamente lo que eran antes de que
llegaran las anclas: la tirada ya se hacía para cada parada, y se dejaba sin
usar.

Lo que sigue a un color es adonde lo llevan **sus degradados**. Allí donde la
pintura de las viñetas cambia suavemente — a lo largo de ocho píxeles, sin un
paso mayor que un tercio del cambio — se contaron los colores de los dos
extremos, y cada parada sortea la siguiente a partir de esos recuentos: el
papel sigue hacia una aguada de amarillo y luego amarillo, el negro hacia el
bermellón (tres veces de cada diez) y el bermellón hacia el burdeos, el azul
hacia el malva y el malva hacia el rosa, el cerúleo hacia un aguamarina gris.
También se contaron los bordes duros, y se dejaron fuera: la mayoría son papel
junto a papel de otro matiz, y las paletas hechas con ellos salían papel y
gris. Una paleta empieza en el papel o en el negro.

Sobre veinte mil paletas, las de la décima parte superior de ánimos distan 0.06
de las viñetas, con la misma divergencia de arriba, y las de la décima parte
sobre cero 0.07 de la imagen mate; de las otras siete paletas distan de 0.39 a
0.67.

Como un degradado se queda un rato en una parte del círculo cromático, y el 8
no tiene negro y blanco cada pocas paradas a los que volver, las pruebas le
piden menos que a los demás. Donde la paleta se extiende larga, cuarenta
entradas por parada o más, una imagen ve unas diez paradas, y una paleta de
cada cien es entonces un único campo — ocre con salvia, petróleo con bermellón
— sin llegar a un tercio del camino del negro al blanco ni a un segundo color.
Y tres paradas, que es lo que a menudo es una paleta corta de 256 entradas, son
un único degradado.

## Observar la órbita

`trap` y `stripe` guardan un número sobre toda la órbita en lugar de sobre el
punto, que es el tipo de magnitud que los modos de coloreado del motor no
pueden alcanzar: allí un bucle de cálculo está compilado para devolver un
color, y solo el último `z` y el anterior sobreviven para usarse en el
coloreado. Hecho en la fórmula, al motor no le cuesta nada.

Ambas devuelven el valor que recibieron hasta la última pasada que permite el
límite de iteraciones, y en esa lo que han reunido. Así que una fórmula
completa es

    trap(z^2+c, 3)

y nada más: el fractal itera como lo haría, y en la última pasada el valor pasa
a ser la trampa, que luego dibujan los modos de coloreado interior. Un punto
que escapa nunca llega a esa pasada y conserva su color exterior habitual, así
que lo que dibujan estas es el interior.

**`trap(orbit, shape, centre, size)`** — lo más cerca que la órbita llegó nunca
a una forma. `shape` vale 0 por defecto, `centre` el origen, `size` 1.

| shape | |
| --- | --- |
| `0` | el propio centro, un punto |
| `1` | una recta horizontal que pasa por él |
| `2` | una vertical |
| `3` | ambas, una cruz |
| `4` | un anillo de radio `size` |
| `5` | un cuadrado de semilado `size` |
| `6` | un rombo de semidiagonal `size` |

**`stripe(orbit, density)`** — la media de `(sin(density * arg orbit) + 1) / 2`
a lo largo de la órbita. `density` vale 4 por defecto y es cuántas franjas dan
una vuelta; un número entero, o no coinciden donde se cierra la vuelta. Una
media sobre la órbita cambia suavemente con el punto incluso donde el recuento
de iteraciones salta, y eso es lo que dibuja las fibras por las que se conoce
el método.

    trap(z^2+c,0)                    lo cerca que pasó la órbita del origen
    trap(z^2+c,3)                    una cruz, que dibuja rayos
    trap(z^2+c,4,{0,0},{0.5,0})      anillos dentro de la cardioide
    stripe(z^2+c,6)                  seis franjas por vuelta

La magnitud acumulada vive en el punto de llamada, así que dos trampas en una
fórmula guardan cada una la suya y un hilo no puede perturbar a otro. Ambas
desactivan el seguimiento del borde, por la misma razón que el ruido: dos
vecinos que necesitan el mismo número de pasadas pueden haber visto órbitas muy
distintas, y una región rellenada sin calcularla sería errónea.

No son modos de coloreado, y no pueden serlo: los modos viven dentro del bucle
compilado de cada fórmula, que devuelve un color directamente, y hay treinta y
uno, cada uno con su tabla de simetría por posición. Lo que hay aquí funciona
solo para una fórmula de usuario.

## Forma del escape

Cálculo → Radio de escape → Forma del escape. La iteración se detiene cuando
`z` sale de una región, que siempre ha sido un círculo. La forma se ve en las
bandas fuera del conjunto, así que es tanto una herramienta de dibujo como una
numérica.

| forma | escapa cuando |
| --- | --- |
| Círculo (clásico) | \|re\|² + \|im\|² supera el bailout |
| Cuadrado | lo hace cualquiera de las dos componentes sola |
| Rombo | lo hace \|re\| + \|im\| |
| Eje real | la componente real sola |
| Eje imaginario | la componente imaginaria sola |
| Ambos ejes | las dos componentes a la vez |
| Triángulo 0° / 90° / 180° / −90° | el punto sale de un triángulo con esa orientación |
| Hexágono 0° / 90° | el punto sale de un hexágono con esa orientación |
| Octógono | el punto sale de un octógono |

Los polígonos se miden por su **apotema** — la distancia del centro al punto
medio de un lado — así que un polígono y el círculo del mismo bailout se tocan
a lo largo de los lados y no en las esquinas. La orientación es parte de la
forma: un triángulo girado 120° es el mismo triángulo, así que 0°, 90°, 180° y
−90° son cuatro distintos, mientras que un hexágono se repite cada 60° y por
tanto tiene dos, y un octógono cada 45°, lo que no le deja nada que elegir. El
triángulo y el hexágono tienen por eso un submenú propio para la orientación.
Sus normales se calculan una vez por fotograma, así que el bucle interior
cuesta de tres a ocho multiplicaciones-sumas en lugar de un seno y un coseno
por lado.

Frente al círculo, el rombo cambia cerca del 39% de los píxeles y el eje real
el 44%. La silueta del conjunto no se mueve — la decide lo que nunca escapa —
así que mira las bandas.

Se escribe en una posición guardada solo cuando no es el círculo, así que una
posición que no pide una forma sigue pudiendo cargarse en cualquier versión
anterior.

## Movimiento browniano fraccionario

Un diálogo a cada lado — **Fractal → Coloreado interior → Movimiento browniano
fraccionario...** y lo mismo bajo Coloreado exterior — con los tres modos que
ofrece y los cinco números con que se dibujan. Era un submenú, los tres modos
como entradas de elección y los números detrás de un diálogo Ajustes debajo de
ellos; ahora es un único diálogo, marcado cuando uno de los tres modos está en
uso, y Aceptar fija el modo y los números juntos como un solo paso de deshacer.
Una posición sigue guardando el modo como `(incoloring N)` o `(outcoloring N)`
y los números aparte, como `(infbmset ...)` o `(outfbmset ...)`, así que un
archivo de una versión anterior se abre sin cambios.

Lo que hace es que un coloreado limpio parezca **usado**: las bandas conservan
su forma y su orden y dejan de ser perfectas, como si la cosa se hubiera dejado
a la intemperie. No es aleatoriedad píxel a píxel, que se lee como grano y se
probó y se descartó; es un campo de ruido con una escala, así que lo que dibuja
son marcas.

**Por fuera** los tres modos son `fbm + smooth`, `fbm + iter` y `fbm` solo — el
recuento suave desgastado, el recuento por bandas desgastado, y el ruido solo,
que dibuja nubes alrededor del conjunto sin ningún recuento dentro.

**Por dentro** son `fbm`, `fbm + zmag` y `fbm + decomposition`. El primero es
el interesante: el interior es por defecto un único tono plano, y esto le da
una superficie.

### Ruido de gradiente, y ninguna rejilla

El movimiento es **ruido de gradiente**: cada esquina de una retícula recibe
una dirección, y un punto toma de cada esquina la pendiente que esa dirección
da a su distancia. Antes era ruido de valor — un número en cada esquina,
mezclado con una smoothstep — y eso dibujaba una rejilla. Una smoothstep no
tiene pendiente en ninguno de sus extremos, así que el campo se volvía plano a
lo largo de cada línea de la retícula, y como la retícula de cada octava
contenía las líneas de la primera, todas se volvían planas ahí a la vez. A
través de una línea la pendiente era cero; la imagen salía en cuadrados, de
ocho píxeles de lado en `FBM_ERROR.xpf`, como si se hubiera calculado a una
resolución menor que aquella a la que se mostraba. El ruido de gradiente es tan
empinado sobre las líneas como entre ellas — medido, 0.99 con los ajustes
predeterminados y 1.08 con los de esa posición — y la retícula de cada octava
está desplazada respecto a la anterior, en cantidades que ninguna duplicación
vuelve enteras, así que no hay dos que compartan una línea ni una esquina.
`fbm()` en las fórmulas dibuja el mismo movimiento, y `randsc` es una octava de
él.

Lo que se mantuvo: el valor sigue yendo de cero a la intensidad y nunca fuera;
las **marcas tienen el tamaño que tenían**, con la retícula a tres quintos de
la frecuencia para que el ruido de gradiente las dibuje tan grandes como lo
hacía el ruido de valor (medido por cuánto se mueve el campo antes de parecerse
a sí mismo solo la mitad: 0.56 de celda con los ajustes predeterminados, antes
y después); el **contraste**, del que el ruido de gradiente solo tiene siete
décimas — cada octava pasa por una curva que lleva sus valores a la
distribución que tenía el ruido de valor, monótona para que no mueva ningún
contorno, y medido sobre el movimiento con los ajustes predeterminados oscila
0.133 desde su centro frente a 0.132 antes; y el coste, más o menos — los modos
algo más rápidos que antes, `fbm()` cerca de un octavo más lenta. Lo que no:
una posición guardada con un modo fbm conserva sus números, su escala y su
fuerza y **dibuja marcas distintas**.

### Los cinco números

El mismo diálogo que el modo, debajo de él. Se guardan con la posición.

| | |
| --- | --- |
| **intensidad** | cuántas bandas de color mueve el ruido el valor. `4` es el predeterminado; por debajo de uno es un susurro, por encima de diez las bandas dejan de ser bandas |
| **frecuencia** | el tamaño de las marcas, tan grandes como las dibujaría una retícula de ese número de celdas por unidad del plano. **Es la que hay que subir al hacer zoom**: en la primera vista `8` se lee bien, con una amplitud de `0.02` hacen falta unos `200` |
| **octavas** | cuántas se suman, cada una al doble de la frecuencia de la anterior |
| **rugosidad** | lo que cada octava conserva de la altura de la anterior. `0.5` es el movimiento simple; más alto es más arenoso |
| **semilla** | la misma imagen cada vez. Cámbiala para otra del mismo carácter |

**Las octavas y la rugosidad se sostienen mutuamente.** Con una rugosidad de
`0.5` la octava octava aporta una doscientas cincuentava parte del total, así
que las octavas más allá de la cuarta o la quinta no cambian nada visible —
medido, cuatro y dieciocho dieron la misma imagen. Sube la rugosidad y las
octavas finas tienen algo que gastar, y entonces el número empieza a importar.

### El tamaño de las marcas según el recuento

Por fuera, dos campos más al pie del diálogo. **Escalado del ruido** dice qué
decide el tamaño de las marcas: **fijo**, donde empieza, la frecuencia sola,
como arriba; **iteración**, la frecuencia subida con el recuento de
iteraciones; **iteración suavizada**, la frecuencia subida con el recuento
suave. **Factor** dice cuánto: la frecuencia se duplica ese número de veces
cada diez iteraciones, así que a `1` las marcas son de la mitad de tamaño donde
la órbita tardó diez iteraciones más en salir, y de un cuarto donde tardó
veinte. A `0` el recuento no cambia nada. Donde el recuento sube, hacia el
borde del conjunto, las marcas se hacen más finas y siguen el detalle.

Una frecuencia que cambia de un píxel al siguiente no puede tenderse sin más
sobre el plano: el campo se estiraría, tanto más cuanto más lejos del origen.
Así que el movimiento se lee en las dos potencias de dos a cada lado del
recuento, cada una un campo propio a frecuencia fija, y las dos se mezclan como
una textura entre sus mipmaps, devolviendo la mezcla al contraste de un solo
campo. Con **iteración suavizada** las marcas cambian de tamaño sin costuras;
con **iteración** cambian en los bordes de las bandas, cada banda con un grano
propio. Por dentro el recuento es el mismo en todas partes, así que el interior
no tiene ninguno de los dos campos.

Entre dos duplicaciones se leen dos campos, y el coloreado cuesta más o menos
el doble que con **fijo**: 248 ns por píxel frente a 127, con cuatro octavas.
Donde el recuento cae justo en una duplicación, como hace todo recuento entero
con un factor de `10`, se lee un solo campo y el coste es el de **fijo**. Con
**iteración suavizada**, `fbm + iter` y `fbm` los dibuja el bucle suave, que
hace unas pocas iteraciones más allá del bailout para calcular el recuento
suave. Un factor demasiado alto lleva las marcas por debajo de un píxel, donde
se vuelven grano.

Con **fijo** dibuja como antes, bit a bit. Una posición guarda los dos como
`(outfbmscale 'smooth 1.5)`, y solo cuando se han movido, así que una posición
que no los toca sigue abriéndose en versiones anteriores.

### Lo que cuesta, y lo único que no hará

El ruido se calcula **una vez por píxel**, en el coloreado, no una vez por
iteración. Sobre una imagen de 90000 puntos con una amplitud de `0.02` y 400
iteraciones, `fbm + smooth` tardó un quinto más que `smooth` solo con cuatro
octavas y un cuarto más con dieciocho; con 4000 iteraciones, donde la órbita
hace más del trabajo, las mismas cuatro octavas costaron una décima.

Las marcas están ancladas **al plano**, no a la vista: es suciedad sobre la
cosa, no sobre la lente. Así que hacer zoom las amplía, y si vas lo bastante
lejos acabas dentro de una sola marca y la imagen vuelve a estar limpia. Subir
la frecuencia es la respuesta, y por eso la frecuencia es un ajuste y no una
constante. No puede atarse al zoom en su lugar: el motor reutiliza píxeles
calculados a la escala anterior, y un coloreado que se moviera con la vista les
dejaría las marcas equivocadas.

Los modos se numeran **después** de los de color real, así que cada número de
modo que ya contiene una posición guardada sigue significando lo que
significaba.

## Teselado

**Fractal → Coloreado interior → Teselado...** colorea el interior según las
teselas de uno de los cuarenta y cinco teselados de `randsctile`, tendidas
sobre la `z` final, cuarenta por unidad como los cuadrados de `squares`. Su
diálogo tiene cuatro campos.

| campo | |
| --- | --- |
| **Teselado** | qué teselado, con el número con que lo conoce `randsctile` y su nombre: **3: hexágonos**, donde empieza. La pestaña Teselados de Ayuda → Referencia de fórmulas de usuario tiene una imagen de cada uno |
| **Coloreado** | de dónde viene el color, en todas partes: **ángulo de z**, donde empieza; **ángulo en torno al centro de la tesela**, el de la tesela en la que cae el punto, que hace de cada tesela una pequeña rueda de color; **zmag**, el valor del modo de ese nombre; **ninguno**, nada más que las teselas |
| **Teselas** | cómo se separa cada tesela de sus vecinas: **cada una desplazada en una cantidad propia**, a lo largo de la paleta, donde empieza; **sombreadas del centro al borde**, lo mismo, y desplazada un octavo más de los colores de una vuelta del centro al borde, como un bisel; **con el borde marcado**, lo mismo, y la parte exterior de la tesela dibujada en la otra mitad de la paleta; **un color liso cada una**, un color propio, venga de donde venga el color |
| **Intensidad** | cuánto se separan las teselas a lo largo de la paleta: multiplica cada cantidad propia de una tesela, su desplazamiento, el bisel de **sombreadas del centro al borde** y el borde de **con el borde marcado**. En `1`, donde empieza, las teselas son como se describen aquí; en `0` solo queda aquello de lo que viene el color, así que **ángulo de z** es una descomposición simple y **zmag** el zmag simple; en medio, la imagen en `0` se funde con la imagen en `1` en esa medida, de modo que `0.1` trae las teselas una décima parte, y por encima de `1` las cantidades se multiplican y las teselas destacan más. Una paleta de 256 colores, que no tiene colores entre los suyos con los que fundir, también en medio desplaza las teselas esa parte de sus cantidades |

El modo está entre los modos que se repiten y las descomposiciones, y los
ajustes de ambos lo alcanzan. **Frecuencia (modos periódicos)** y **Función
(modos periódicos)** alcanzan las teselas como alcanzan todo argumento
periódico del interior: a frecuencia 2 las teselas tienen la mitad de tamaño.
Los [ajustes de las descomposiciones](#ajustes-de-las-descomposiciones) del
interior alcanzan los dos ángulos, que se cortan y se recorren como los de las
descomposiciones. **Sectores**, **Rotación**, **Repeticiones**, **Sentido** y
**Velocidad** alcanzan a ambos; **Ángulo de** alcanza solo el ángulo de `z`,
porque el de en torno al centro de la tesela no tiene uso para él; **Anillos
por unidad** hace anillos del módulo para el ángulo de `z`, a través de la
frecuencia y la función como en las descomposiciones, y para el ángulo en torno
al centro otros tantos anillos del centro al borde, que siguen la forma de la
tesela. **zmag**, **ninguno** y una tesela lisa no leen ninguno de ellos.

La cantidad propia de una tesela viene del hash que la nombra. Para un ángulo
abarca tanta paleta como una vuelta del ángulo, así que un ajuste de la
descomposición movido deja los colores de las teselas donde estaban; para
**zmag**, **ninguno** y una tesela lisa abarca la paleta entera, así que
**ninguno** con **cada una desplazada en una cantidad propia** es lo mismo que
**un color liso cada una**. **ángulo en torno al centro de la tesela** pone una
vuelta entera de colores en cada tesela, que al tamaño que tienen las teselas
son muchísimos colores en pocos píxeles: una **Velocidad** por debajo de uno, o
una frecuencia por debajo de uno, hace legibles las ruedas.

Con **Coloreado** en **ángulo de z**, **Teselas** en **cada una desplazada en
una cantidad propia** y los ajustes de las descomposiciones donde empiezan, el
modo dibuja lo que dibujaba antes de tener elección, bit a bit. Los teselados
que nunca se repiten quedan planos, sin cantidad propia, más allá de cuatro mil
millones de teselas del origen. Eso vale con **Intensidad** en `1`. En `0` las
teselas aún se ven con **ángulo en torno al centro de la tesela**, cuyo ángulo
se mide desde el centro de cada tesela.

Los teselados periódicos cuestan más o menos lo que los demás modos interiores,
de 120 a 210 ns por píxel; los que nunca se repiten, que se localizan de nuevo
en cada punto, de 250 a 490. **ángulo en torno al centro de la tesela**,
**sombreadas del centro al borde** y **con el borde marcado** miden dónde está
el punto en su tesela, lo que cuesta hasta dos tercios más en los teselados
periódicos y poco más en los demás. Una posición guarda el modo como
`(incoloring 30)`, el teselado como `(intessellation 39)`, solo cuando no es el
de hexágonos, y el coloreado como `(intessellstyle 'tile-angle 'edged)`, solo
cuando se ha movido. La intensidad se guarda como `(intessellintensity 0.5)`,
solo cuando no es uno.

### En el exterior

**Fractal → Coloreado exterior → Teselado...** tiene tres modos, como
Movimiento browniano fraccionario: `tessellation + smooth`, las teselas sobre
el recuento suavizado; `tessellation + iter`, sobre el recuento entero;
`tessellation`, las teselas solas. Las teselas se tienden sobre el lugar del
píxel en el plano, no sobre el punto de escape, que salta de una franja a la
siguiente y solo daría grano, como pasó cuando se probó allí el movimiento.

| campo | |
| --- | --- |
| **Modo** | uno de los tres |
| **Teselado**, **Coloreado**, **Teselas**, **Intensidad** | como en el interior, con ajustes propios. **ángulo de z** es el ángulo del punto de escape y **zmag** el logaritmo de su módulo, como lo suma `iter+log(mag)`; salvo en `tessellation + iter` con escalado **fijo** o **iteración**, ambos se leen unas pocas iteraciones más allá del bailout, donde se calcula el recuento suavizado, y se funden con el mismo punto una iteración después en la medida en que el recuento suavizado ha avanzado hacia la iteración siguiente, de modo que siguen sin cortes a través de los bordes de las iteraciones en lugar de romperse allí |
| **Frecuencia** | teselas por unidad del plano, `8` donde empieza: el tamaño de las teselas, que hay que subir al acercarse, como la frecuencia del movimiento |
| **Escalado de las teselas** y **Factor** | qué da el tamaño a las teselas, como [El tamaño de las marcas según el recuento](#el-tamaño-de-las-marcas-según-el-recuento) se lo da a las marcas del movimiento: **fijo**, la frecuencia sola; **iteración** o **iteración suavizada**, la frecuencia duplicada tantas veces como dice el factor cada diez iteraciones de ese recuento |

Entre dos duplicaciones las teselas se leen en los dos tamaños y los dos
colores se funden, un fundido de un tamaño al siguiente: las teselas no tienen
un valor intermedio hacia el que fundirse, como tiene el movimiento. Una paleta
de 256 colores, que no tiene ninguno entre los suyos, toma el tamaño más
cercano. Entre dos duplicaciones el teselado se lee dos veces, lo que duplica
su coste. Un factor de `1` duplica la frecuencia diez veces cada cien
iteraciones, lo que pronto lleva las teselas por debajo del píxel: unas décimas
son un mejor comienzo. Los [ajustes de las
descomposiciones](#ajustes-de-las-descomposiciones) del exterior alcanzan los
dos ángulos: **Iteraciones por franja** hace franjas del recuento para el
ángulo de z y anillos del centro al borde para el ángulo en torno al centro.

Con **Intensidad** en `0` y **Coloreado** en **ninguno**, `tessellation + iter`
dibuja como `iter` y `tessellation` un solo color. Una posición guarda el modo
como `(outcoloring 21)` a `(outcoloring 23)`, y los ajustes en un solo comando,
`(outtessellation 5 'tile-angle 'shaded 0.6 40 'smooth 0.4)`, solo cuando
alguno se ha movido.

## Más modos de coloreado

No es el cálculo lo que decide el color de un píxel: el color sale de una
fórmula en la `z` final, el parámetro `c` y el recuento de iteraciones, y hay
una fórmula así por cada modo. Veintidós más, en un submenú propio bajo Fractal
→ Coloreado interior y Fractal → Coloreado exterior, por encima de los submenús
de color real.

**Por dentro**, donde la órbita nunca escapó y solo la `z` final y `c` dicen
algo:

| modo | |
| --- | --- |
| `real/mag` | el coseno del ángulo: un sombreado que se cierra sobre sí mismo donde `atan2` tiene una costura |
| `max(\|real\|,\|imag\|)` | la norma cuadrada, así que las bandas son cuadrados donde las de `zmag` son círculos |
| `\|real\|+\|imag\|` | las mismas bandas giradas un octavo de vuelta |
| `min(\|real\|,\|imag\|)` | cerca de cero a lo largo de cualquiera de los ejes, así que rayos |
| `\|z-c\|` | lo lejos del parámetro que se asentó la órbita, lo que distingue los bulbos |
| `\|z*c\|` | los dos módulos uno contra otro |
| `angle(z)-angle(c)` | la dirección de asentamiento medida respecto al parámetro |
| `real*imag` | una silla, cuatro lóbulos alrededor de los ejes |
| `sin(real)*sin(imag)` | una rejilla puesta sobre donde se asentó la órbita |
| `sign(imag)` | dos tonos planos: una descomposición del interior, como la descomposición binaria hace con el exterior |
| `frac(mag)` | curvas de nivel en cada octavo del módulo |
| `log(mag)` | `zmag` con el contraste llevado a los valores pequeños, donde pasa su tiempo una órbita atraída |
| `cos(\|z\|)` | anillos alrededor de donde se asentó la órbita, ocho por unidad del módulo |
| `hexagons` | tres ondas separadas un tercio de vuelta, que sumadas dan una red de hexágonos |
| `spiral` | cinco brazos alrededor de donde se asentó la órbita, enrollados una vez por cada factor e del módulo |
| `cos(\|z-c\|)` | anillos alrededor del parámetro, que distinguen los bulbos como hace `\|z-c\|` |

Once de los modos interiores se repiten — `cos(mag)`, `mag*cos(real^2)`,
`sin(real^2-imag^2)`, `squares`, `sin(real)*sin(imag)`, `frac(mag)`, los cuatro
al pie de la tabla y `tessellation` — y **Frecuencia (modos periódicos)**, el
primer campo del diálogo Velocidad del coloreado interior, dice cada cuánto.
Multiplica el argumento de lo que se repite, así que a 2 los senos giran el
doble de rápido, `squares` dibuja cuadrados de la mitad de tamaño y `frac(mag)`
dieciséis contornos por unidad del módulo donde dibujaba ocho. En `cos(mag)`
estrecha también las bandas que eligen entre el seno y el coseno. En uno, donde
empieza, dibujan como siempre, bit a bit, y los demás modos nunca lo leen. Una
posición lo guarda como `(incolorfrequency 2)` solo cuando no es uno, así que
una posición que no lo toca sigue abriéndose en versiones anteriores.

**Función (modos periódicos)**, debajo de ella, dobla el argumento con una de
las [funciones de velocidad del
coloreado](#funciones-de-velocidad-del-coloreado) antes de que la frecuencia lo
multiplique, de modo que las bandas giran más deprisa donde la función sube
empinada y más despacio donde se aplana: `square` y `cube` las apiñan hacia el
borde y las espacian alrededor de cero, `log` y `atan` al revés, y `sin`, `cos`
y `triangle` repliegan el argumento sobre sí mismo, de modo que las bandas dan
la vuelta y regresan. Alcanza cada función periódica de los once, cuyo
argumento `a` entra como `k·f(a)`, siendo `k` la frecuencia: los senos y
cosenos, el módulo de `frac(mag)`, la parte real y la imaginaria cuyas franjas
forman los cuadrados de `squares` y las baldosas de `tessellation`, y en
`cos(mag)` también las bandas que eligen entre el seno y el coseno. En `spiral`
los cinco brazos siguen siendo cinco, porque un número de vueltas no entero
dejaría una costura donde el ángulo se cierra: la frecuencia y la función
actúan sobre el enrollamiento. Una función periódica dentro del argumento de
otra tendría también el suyo doblado. En `ident`, donde empieza, nada cambia.
Una posición la guarda como `(incolorperiodicfun 6)`, por el lugar de la
función en la lista contando desde cero, solo cuando no es ident.

**Por fuera**, donde lo que hay es el recuento y el punto de escape:

| modo | |
| --- | --- |
| `iter+angle` | la dirección por la que salió, mezclada en el recuento: las bandas adquieren un giro |
| `iter+log(mag)` | cuánto pasó del bailout, un suavizado más barato que el logaritmo de un logaritmo |
| `iter+real*imag` | biomorfos hechos continuos en lugar de por umbral |
| `max(\|real\|,\|imag\|)` | por qué lado salió, sin recuento alguno |
| `iter banded` | el recuento plegado en ocho bandas, que muestra la forma de los conjuntos de nivel en lugar de su número |
| `\|real\|-\|imag\|` | lo desequilibrado que está el punto de escape |

Los números de los modos antiguos no se mueven, así que una posición guardada
sigue cargándose como antes.

La velocidad y el desplazamiento del coloreado llegan ahora a smooth y smooth
log, cosa que nunca hacían: esos dos devolvían un píxel propio en lugar de
pasar por el paso con que terminan los demás modos, así que ambos controles no
hacían nada mientras uno de ellos estaba elegido. Con los valores con que
vienen — velocidad uno, desplazamiento nulo — el paso es la identidad, así que
ninguna imagen guardada antes se mueve.

Smooth y smooth log siguen la forma del escape. Interpolan entre una pasada y
la siguiente preguntando dónde cruzó la órbita el bailout, y antes se lo
preguntaban al módulo fuera cual fuera la forma elegida — así que escapar en un
cuadrado e interpolar en un círculo ponía más de un quinto de los puntos que
escapan, y la mitad con un rombo, fuera de la banda a la que pertenecen. Cada
forma entrega ahora la magnitud que realmente comprueba y el umbral con que la
compara. Para el círculo las dos cosas coinciden, así que nada de lo que se
dibujaba antes se mueve.

También funcionan en modo Newton, donde una órbita converge en lugar de
escapar: la magnitud es entonces cuánto se movió el último paso y el umbral el
límite de convergencia, y la interpolación es la misma pregunta hecha a un
número que baja en lugar de subir. Los dos fractales de Newton, que no tenían
ninguna variante suave, tienen una.

**Convergencia de Newton**, en el menú Cálculo, se lee por fin. El límite
estaba escrito en la iteración como una millonésima, así que el ajuste se
ofrecía, se guardaba y se recargaba mientras nada lo consultaba nunca.

**Smooth y smooth log** se han trasladado a ese submenú en el menú exterior.
Interpolan entre una iteración y la siguiente usando cuánto pasó la órbita del
bailout, lo que necesita una fórmula que escape en el bailout — veinte de las
treinta y una lo hacen, y las otras once o convergen (las Newton) o se detienen
con una prueba propia (las Sierpinski, Koch, el Clock). Para esas el modo no
hace nada en absoluto, y no es algo que la primera lista deba ofrecer. La
fórmula de usuario, que siempre podría haberlo tenido y no lo tenía, ahora lo
tiene.

Lo que estos modos no pueden hacer es nada que necesite la órbita en lugar de
su final: una trampa de órbita, una media de franjas, un exponente de Lyapunov,
la estimación de la distancia. Esas necesitan un valor llevado de una iteración
a la siguiente, lo que en este motor significa una segunda copia compilada de
cada bucle de cálculo. `trap` y `stripe` en una fórmula de usuario lo hacen en
su lugar.

### Ajustes de las descomposiciones

**Fractal → Coloreado interior → Descomposición**, y lo mismo bajo Coloreado
exterior, fija cómo cortan la vuelta las descomposiciones: binary decomposition
y color decomposition en el exterior, `sign(imag)`, `Decomposition-like` y
`fbm + decomposition` en el interior, y los ángulos de los teselados, interior
y exterior. Con cada campo donde empieza — cero, y uno, derecha y uno los tres
últimos — dibujan como siempre, bit a bit.

| campo | |
| --- | --- |
| **Sectores** | en cuántos se corta cada repetición de la vuelta. Cero es lo del modo: los dos semiplanos de las binarias, la vuelta entera de las demás. Las binarias colorean sectores alternos, así que un número impar pasa al par siguiente — dos sectores del mismo color se tocarían. Las de color se convierten en otros tantos escalones planos |
| **Rotación** | cuántos grados se giran los sectores, en sentido antihorario |
| **Ángulo de** | solo en el interior: el ángulo de `z`, de `z-c` o de `z*c`. `z-c` centra los sectores en el parámetro, y `z*c` los gira el ángulo del parámetro |
| **Anillos por unidad** (interior) | anillos del módulo de aquello cuyo ángulo se toma, uno sí y uno no distinguidos: las binarias se vuelven un tablero de ajedrez polar. Pasan por `k·f(a)` como todo argumento periódico del interior, así que los alcanzan la frecuencia y la función de los modos periódicos |
| **Iteraciones por franja** (exterior) | franjas de tantas iteraciones a través de los sectores: con la binary decomposition, la cuadrícula de rayos externos y equipotenciales |
| **Repeticiones** | cuántas veces se repite la vuelta: los sectores, y en las de color la paleta, vuelven a empezar en cada una. Cualquier número, aunque una parte de repetición deja una costura donde el ángulo se cierra |
| **Sentido** | derecha, como fueron siempre; izquierda, al revés; ambidiestro, adelante y atrás en cada repetición — un espejo, sin costura cuando las repeticiones son enteras |
| **Velocidad** | en las de color, cuánta paleta abarca una repetición: por encima de uno, más colores en cada una, no más repeticiones. En las binarias, la separación entre los colores de sectores alternos: en uno como siempre, en cero un solo color, en dos el doble — la binary decomposition lleva los sectores impares del recuento hacia el recuento al revés, `sign(imag)` su segundo tono desde el primero |

Una posición los guarda como `(indecomposition 8 22.5 'z-c 6 2.5 'both 1.5)` y
`(outdecomposition 6 30 2 3 'left 2)` solo cuando han cambiado, así que una
posición que no los toca sigue abriéndose en 2.5. Cambiados, una descomposición
ya no es simétrica respecto al eje real, y XaoS deja de dibujar una mitad como
espejo de la otra.

## Funciones de velocidad del coloreado

**Fractal → Coloreado interior → Velocidad del coloreado**, y lo mismo bajo
Coloreado exterior, es un único diálogo para todo lo que lleva el valor que
entrega un modo de coloreado a un lugar de la paleta. La **función**, elegida
de una lista, dobla el valor antes de que se apliquen la **velocidad** y el
**desplazamiento** del coloreado: la paleta se lee en `f(v) × speed + shift`.
En una columna propia, al lado, está el [cambio de velocidad tras
ITER](#cambio-de-velocidad-tras-iter), que actúa sobre el valor antes de que lo
reciba la función. El interior empieza con la frecuencia y la función de los
modos que se repiten, que actúan sobre el propio modo; véase [Más modos de
coloreado](#más-modos-de-coloreado). Aceptar lo fija todo como un solo paso de
deshacer. Para los modos que cuentan — iter, smooth y los demás — `v` es el
propio recuento de iteraciones, así que una función se puede pensar en
iteraciones. Todas son impares, `f(−v) = −f(v)`, así que una magnitud que
cambia de signo lo conserva.

Eran cuatro entradas del menú: la velocidad, el desplazamiento, un submenú de
las funciones y el diálogo propio del cambio de velocidad. Una posición se
escribe y se lee exactamente como antes, y los comandos que esas entradas
ejecutaban siguen ahí para un script.

`sin`, `cos` y `triangle` dan vueltas: la velocidad del coloreado, aplicada
después de ellas, dice cuántas entradas de la paleta abarca una vuelta, y deja
como estaba cada cuánto vuelve. Eso lo hace la **velocidad inicial** del cambio
de velocidad, aplicada antes: con los dos grupos apagados y el exponente
inicial en uno la función recibe `s0·v`, así que 2 junta las vueltas al doble y
0,5 las separa al doble.

| | |
| --- | --- |
| `ident`, `log`, `atan`, `sin`, `cos`, `square`, `cube`, `square root`, `cube root` | lo que siempre hubo, sin cambiar ni un bit |
| `log log` | `log(1 + log(1 + v))`, para los zooms profundos donde el log solo aún aprieta demasiado las bandas: con un millón de iteraciones log da 13.8 bandas y esta 2.7 |
| `exp` | `e^v − 1`, lo contrario del log: abre los valores pequeños. Para las magnitudes que se quedan pequeñas — las componentes que leen el color real y las descomposiciones; con un recuento de noventa ya pasa de 10³⁸ |
| `triangle` | el seno con las curvas enderezadas: el mismo periodo y la misma altura, la paleta recorrida hacia delante y hacia atrás a paso constante en lugar de demorarse en las vueltas |

Las nuevas van detrás de las antiguas, así que cada número que guarda una
posición sigue nombrando lo que nombraba.

### Cambio de velocidad tras ITER

La parte inferior del mismo diálogo no es una función sino un cambio a la
función que esté elegida: hasta una iteración que indiques, la función sigue
exactamente como era; pasada esa, el recuento que recibe corre más rápido o más
despacio, así que los colores se aceleran o se frenan desde ahí y no antes. Con
`square root` elegida e ITER en 40, los colores son los de la raíz cuadrada
hasta 40 iteraciones y la raíz cuadrada de un recuento que va cogiendo ritmo
después. Contiene un exponente y una velocidad iniciales y luego dos grupos,
cada uno separado por una línea:

| | |
| --- | --- |
| **Exponente inicial** | el exponente desde el principio hasta el primer ITER de un grupo activado: el recuento pasa a ser la potencia simple `v^p0`, la misma estén los grupos activados o no, y con ambos desactivados rige para todo recuento. En uno, donde empieza, el recuento queda intacto hasta el primer ITER |
| **Velocidad inicial** | por lo que se multiplica el ritmo de ese tramo: `s0·v^p0`. En uno, donde empieza, no cambia nada |
| **ITER** | dónde empieza el tramo del grupo, en las unidades que recibe la función. En cero o por debajo empieza al principio — lo que hace de `ident` una potencia con cualquier exponente |
| **Exponente** | el exponente del tramo, desde su ITER hasta el siguiente o en adelante: cómo se curva el tramo. Parte del recuento al que llegó el anterior, `f(ITER)`, y suma `s·ITER·((v/ITER)^p − 1)`, cuyo ritmo en ITER es `s·p`: por encima de uno los colores se aceleran desde ahí, por debajo de uno se frenan |
| **Velocidad** | por lo que se multiplica el ritmo del tramo — lineal donde el exponente es una potencia: con un exponente de uno el tramo es la recta `f(ITER) + s·(v − ITER)`. En uno, donde empieza, no cambia nada |
| **Modo** | **off**, **jump** o **smooth**: si el grupo se usa, y cómo se une su tramo con el anterior |
| **Radio** | solo en smooth, y mostrado solo entonces: cuántas iteraciones ocupa la unión a cada lado de ITER. Cualquier número desde cero, que es el jump |

Con ambos grupos activados y el primer ITER por debajo del segundo, el recuento
es `v^p0` hasta `i1`, sigue desde ahí al ritmo de `p1`, y desde `i2` al ritmo
de `p2` — **cada tramo con su propio exponente y su propia velocidad**, nunca
`p0·p1` ni `p0·p1·p2`. Con `p0` = 2, `p1` = 1 en 50 y `p2` = 2 en 100, todos en
jump y todas las velocidades a uno, el recuento sobre `ident` es `v²` hasta 50;
`v + 2450` de 50 a 100, un exponente de uno que sigue a ritmo uno desde los
2500 a los que llegó el cuadrado; y desde 100 `2550 + 100·((v/100)² − 1)`, de
nuevo al cuadrado desde un ritmo de 2 ahí. Una velocidad de 3 en el primer
grupo hace que el tramo del medio sea `2500 + 3·(v − 50)`. El ITER más bajo va
primero lo tenga el grupo que lo tenga, y **con el mismo ITER el segundo grupo
sustituye al primero**, radio, velocidad y todo; un grupo con un ITER de cero o
menos pone su propio exponente y su propia velocidad en el primer tramo.

**Un grupo desactivado se deja fuera**, con su ITER:

| primer grupo | segundo grupo | el recuento corre a |
| --- | --- | --- |
| activado | activado | `p0` hasta el ITER más bajo, luego el exponente de ese grupo hasta el otro ITER, luego el del otro |
| desactivado | activado | `p0` hasta el ITER del segundo, luego `p2` |
| activado | desactivado | `p0` hasta el ITER del primero, luego `p1` |
| desactivado | desactivado | `p0` para todo recuento, la potencia simple `v^p0` |

Ambos empiezan desactivados, y el exponente inicial y todas las velocidades a
uno, donde al valor no se le hace nada en absoluto, así que cada función da
hasta el último bit lo que daba antes de que existiera tal cosa.

Las dos formas de unir, según cambian el recuento que se entrega a la función
(escritas para `ident`, donde el recuento es lo que ve la paleta, con `f(ITER)`
el recuento al llegar a ITER — `ITER^p0` en el primero — `p` el exponente del
grupo, `s` su velocidad y `r` su radio):

| modo | de ITER − r a ITER + r | más allá |
| --- | --- | --- |
| **jump** | — | `f(ITER) + s·ITER·((v/ITER)^p − 1)`: el ritmo gira en ITER de golpe, lo que se ve como un anillo donde las bandas cambian de ritmo |
| **smooth** | el ritmo va en línea recta desde el que tiene el recuento en `ITER − r` hasta el que tiene en `ITER + r`, a medio camino en ITER, con una parábola alrededor de ITER añadida que devuelve lo que la recta cortó de la esquina del jump | el recuento del jump, hasta la última cifra |

**Smooth** es simétrico respecto a ITER y cambia los colores dentro del radio y
en ningún otro sitio: antes de `ITER − r` y más allá de `ITER + r` la imagen es
la del jump. Donde ambos ritmos son constantes la parábola es nula y la esquina
se corta como la corta una parábola — con el recuento intacto hasta un ITER de
50 y un exponente de cero más allá, un radio de 10 lo deja en 47.5 en ITER, un
cuarto del radio menos, y de nuevo en 50 desde 60; de un ritmo de uno a uno de
tres es dos en ITER. Donde un ritmo se curva, la parábola devuelve lo que la
curva añade: tras `v²` hasta 50 y un exponente de uno, un radio de 20 lleva el
ritmo de 60 en 30 a 1 en 70, 45.5 en ITER, y sigue desde 70 como `v + 2450`,
igual que el jump.

Antes mezclaba los ritmos de los dos tramos tal como cada uno habría seguido:
el ritmo del cuadrado, todavía creciendo, seguía creciendo bien dentro del
redondeo, así que apenas bajaba antes de ITER y estaba a medio bajar solo medio
radio después, y lo que ganaba movía todos los colores más allá del redondeo —
133 bandas con un radio de 20.

La parábola nunca lleva el ritmo por debajo del más lento de los dos extremos.
Donde lo haría — un ritmo lento convertido en uno muy rápido en un radio corto,
como una potencia centésima en diez iteraciones — se mantiene ahí, y el
recuento más allá del redondeo se desplaza en lo que eso retiene. Un radio más
ancho que el primer ITER se detiene antes de cero, y dos que se solaparían se
reducen por igual para que cada uno quede entre sus vecinos — radios de 40 y 40
entre ITER separados 30 pasan a ser 15 cada uno.

Un exponente de 100 es legal e inútil más allá de ITER: unas pocas iteraciones
después, las bandas son más finas que los píxeles y el color es ruido. Con el
segundo grupo en ITER 50 y un exponente de 2 y el primero en 100 y 1.5, ambos
en jump, y un exponente inicial de uno, el recuento sobre `ident` es 30 en 30,
`50·1.5²` = 112.5 en 75, y `200 + 100·(1.5^1.5 − 1)` = 284 en 150: el tramo
desde 100 parte de los 200 a los que llegó el anterior, a un ritmo de 1.5.

Antes había cuatro modos, y las posiciones guardadas antes del radio los
nombran — `outcolorfunafter` de la 1.8.0, y `outcolorfunchange` con dos grupos
y sin radio. Se leen en los dos que hay, con todas las velocidades a uno:

| guardado como | leído como |
| --- | --- |
| off, jump | lo mismo |
| smooth, que redondeaba la esquina en el logaritmo del recuento, un poco a cada lado de ITER por lejos que estuviera | smooth con un radio de tres cuartos de ITER: a menos de cuatro centésimas de ITER de lo que dibujaba hasta un exponente de uno y medio, una décima hasta dos |
| continuous, que seguía desde ITER al ritmo con que llegaba el recuento | smooth a tres cuartos de ITER, con el exponente `0.66 + p/4 + 0.09·p²` — uno en uno, donde continuous tampoco cambiaba nada — a menos de doce centésimas de ITER de él hasta un medio |

Encontrado poniendo las curvas una junto a otra sobre cuatro ITER: una imagen
así sale parecida, no idéntica al píxel. Una posición cuyo cambio de velocidad
no es más que un grupo en jump, con el exponente inicial y todas las
velocidades a uno y el resto como empieza, se guarda como la guardaba la 1.8.0,
`(outcolorfunafter 40 1.5 'jump)`, y se abre en la 1.8.0 y allí dibuja igual.
Cualquier otra cosa se guarda como la guardaba la 2.0, en el orden del diálogo
sin las velocidades, `(outcolorspeedchange 1 40 1.5 'smooth 10 50 2 'off 10)`;
y las velocidades, solo donde alguna no es uno, aparte como
`(outcolorspeedlinear 1 3 1)` — la velocidad inicial y las de los dos grupos —
para que una posición que no las usa siga abriéndose en la 2.0. El cambio de
velocidad solo, como lo tenía su propio diálogo, es el comando
`outcolorspeedgroups`, con sus doce campos en su orden.

### Más allá del final de un int

A la paleta se llega en 256avos de entrada a través de un int, y un valor más
allá de lo que cabe en un int — 8.4 millones de bandas — se convertía de todos
modos, lo cual es indefinido: cada píxel así tomaba una y la misma entrada.
`cube` llega ahí a las 204 iteraciones, así que una imagen coloreada con el
cubo más allá de eso era una mancha plana; `exp` y un cambio de velocidad
grande llegan enseguida. Ahora cae donde dice su valor, dando la vuelta a la
paleta como cualquier valor menor; uno demasiado grande para que sus 256avos
quepan en un double cae en la primera entrada, como siempre hizo un infinito.
Dentro de un int nada se mueve. Fuera, el único cambio a una imagen guardada
antes está en el coloreado suave de Magnet y Magnet2 lejos del conjunto, que
era un único color plano y ahora muestra las bandas que describe ahí el
recuento.

Con ello llegaron dos arreglos menores. Elegir una función de velocidad movía
la marca del menú del coloreado interior en lugar de la del suyo. E
**initstate**, con el que empieza toda posición, no devolvía a sus valores
predeterminados los cinco números del fbm — una posición los guarda solo donde
difieren de ellos, así que una cargada después de otra que los había movido
dibujaba con los de la otra; ahora los restablece, y con ellos desactiva el
cambio de velocidad tras ITER.

### El recuento, entero

El valor que entrega un modo de coloreado perdía cifras de camino a la paleta
en dos sitios. El coloreado suave convertía su recuento en un int de 256avos de
iteración antes de que lo recibieran la velocidad del coloreado, la función o
el cambio de velocidad; y lo que devolvía la función era un float, que más allá
de 65 536 entradas guarda menos de un 256avo de una y más allá de dieciséis
millones ni siquiera una entera. Ninguno se nota con velocidad uno, donde un
256avo de iteración es un 256avo de entrada. Ambos se notan en cuanto algo hace
correr los colores: con un cambio de velocidad que lleva el recuento de cero a
diez mil entradas por iteración, una vista a treinta iteraciones, donde el
ritmo era 2500, salía en franjas de un solo color de ocho píxeles de ancho,
cada una un 256avo de iteración y diez entradas de la siguiente, donde la
paleta debería haber recorrido todas las entradas intermedias.

Ahora el recuento sigue entero, y la función, la velocidad del coloreado y el
cambio de velocidad se calculan en el tipo propio de la compilación, como las
fórmulas — long double en la compilación normal, `__float128` en la quad — y
devuelven lo que calculan. Las dos compilaciones colorean tan parecido como
antes: en las imágenes probadas difieren donde siempre difirieron, en un zoom
profundo cuyas órbitas quad calcula más a fondo, y en ningún otro sitio.

Dos cosas se mantienen como estaban.

* **Un infinito** — el centro de un potencial que se fue al infinito, el fondo
  alrededor de una fórmula que escapa enseguida — cae en la primera entrada,
  como un valor que no es un número. Convertido en 256avos en un float, cada
  valor así salía infinito y caía ahí; qué color debería tomar un infinito no
  lo sabe nadie, y así toma el que siempre tomó. Un valor finito más allá de lo
  que cabe en un float cae donde dice.
* **Un valor en el borde de un 256avo** se queda en él: el valor se redondea a
  un 256avo de un 256avo antes de cortarlo a un 256avo, ya que el coloreado no
  dice nada más fino. La parte real de una órbita de Newton se detiene en una
  raíz, en uno o en menos un medio, lo que pone `iter+real` exactamente en un
  borde así, y las últimas cifras de la convergencia pondrían cada píxel a un
  lado o al otro — un 256avo de entrada, que nadie ve, hasta que un filtro de
  bordes dibuja cada píxel así como borde y las curvas de nivel salen dentadas
  como engranajes. El float redondeaba la mayoría y dejaba muescas; ahora no
  hay ninguna.

Renderizados frente a la compilación anterior a 640×480 — antes del cambio de
la sección siguiente — de los 169 ejemplos 83 dibujaban al píxel como antes y
otros 71 diferían en un solo nivel en como mucho 113 píxeles. De los otros
quince, los de color real y smooth log — `perb4`, `spiral4`, `truecol1`,
`truecol3` — difieren hasta cinco niveles en hasta un séptimo de los píxeles,
al haber vuelto también sus componentes como floats; `edge` ha perdido sus
muescas; y el resto difiere en un puñado de píxeles donde un valor estaba en el
borde de una entrada y un filtro o una paleta dibujada sin mezcla muestra el
escalón. Magnet y Magnet2 en coloreado suave son los que más cambian: lejos del
conjunto su recuento es negativo, lo que el int, leído sin signo, adelantaba
cuatro mil millones de 256avos, y el float dibujaba luego en escalones de una
entrada o dos — manchas planas con bordes duros donde ahora está el degradado
que describe el recuento.

En una imagen que sale igual cuesta del dos al cinco por ciento en la
compilación normal, y del siete al veinte en la quad, donde los logaritmos del
recuento y las funciones se calculan por software. Donde las franjas se vuelven
un degradado el renderizado también tarda más: el relleno por estimación
rellena un píxel solo a partir de vecinos de un único color, y la imagen
escrita es mayor — la vista a treinta iteraciones, a 1200×1200, tardó 1.07
segundos donde sus franjas tardaban 0.30.

### El ritmo del recuento suave

El coloreado suave suma al recuento de iteraciones enteras cuánto de la última
se había recorrido cuando el punto escapó. Lo calculaba donde el bailout cae
entre la magnitud antes de la pasada de escape y después, linealmente en su
logaritmo. Así el recuento sigue sin cortes de una iteración a la siguiente,
pero su ritmo no: de una pasada a la siguiente la magnitud casi se eleva al
cuadrado, así que el ritmo cambia a lo largo de cada iteración y gira de golpe
donde empieza la siguiente. En una vista de 2.6×10⁻⁷ de ancho, donde el
recuento va de 59.94 a 60.02, giraba 4.26 en el borde entre ambas. Con una
velocidad de coloreado de uno eso mueve un color unas centésimas de entrada;
con un cambio de velocidad que corre a diez mil entradas por iteración, son las
bandas cambiando de densidad de golpe a lo largo de cada borde de una
iteración.

Ahora la órbita se lleva más allá del bailout hasta que su módulo al cuadrado
pasa de un millón — tres o cuatro pasadas más para un cuadrado — y la fracción
se toma en el logaritmo del logaritmo, donde elevar al cuadrado es un paso de
longitud uno:

    n + ln(ln T / ln s) / ln(ln m / ln s)

con `s` y `m` el módulo al cuadrado a cada lado de un millón, `T` el bailout y
`n` las pasadas hechas. Las pasadas más allá del bailout se descuentan solas, y
el grado es el de la propia órbita, `ln m / ln s`, así que un cubo se cuenta
como un cubo. El recuento se mantiene a menos de una décima de iteración de lo
que era, y su ritmo sigue: el giro de 4.26 en esa vista sale en 1.0000, y a lo
largo de un rayo desde fuera del conjunto los giros de 1.3 a 1.8 en cada borde
salen en 1.0025.

| | el recuento |
| --- | --- |
| las Mandelbrot, Triceratops, Mandelbar, Lambda, Manowar, Spider, Phoenix y Octo, en el círculo | llevado más allá del bailout |
| una fórmula de usuario que escapa en el círculo | llevado más allá mientras la órbita crece al menos como una potencia de uno y medio, y contado en el bailout donde no |
| las tres Barnsley | en el bailout: multiplican la órbita por `c`, así que su logaritmo crece un paso en lugar de duplicarse, y el recuento en el bailout ya es el correcto para ellas |
| Newton, Magnet, Catseye | en el bailout: convergen, o comprueban su escape a su manera |
| toda forma del escape salvo el círculo | en el bailout |

Las formas se quedan como estaban porque cada una mide la órbita con algo que
gira con ella, mientras que cada pasada duplica su ángulo: llevado más allá con
su propia medida, el grado salía lo bastante inestable como para empeorar los
giros — los del hexágono de 1.28 a entre 1.67 y 3.57 a lo largo de ese rayo — y
los del cuadrado incluso volvían atrás.

De los 169 ejemplos, los 154 que no colorean de forma suave dibujan al píxel
como antes de esto. Los quince que sí cambian hasta trece niveles en hasta
nueve décimas de los píxeles — `eye` hasta 39, donde las entradas vecinas de su
paleta están muy separadas, y `spiral4`, cuyo smooth log sobre una paleta de
color real amplía cualquier cambio en el recuento, en medio de un pequeño
detalle. Frente a la 2.1, con los cambios de la sección anterior, 84 dibujan al
píxel como antes y otros 59 difieren como mucho tres niveles en como mucho 500
píxeles.

Cuesta esas pasadas y dos logaritmos más por cada píxel que escapa: en las
imágenes suaves probadas, del 11 al 28 por ciento en la compilación normal y
del 18 al 20 en la quad, donde un logaritmo en `__float128` tarda 946 ns frente
a 39 en long double; una imagen que no colorea de forma suave no cambia.
Calculados en double, los cuatro logaritmos habrían costado del 8 al 19 por
ciento en la compilación normal y ahorrado el 8 en la quad — un double lleva la
fracción a dieciséis cifras, donde una velocidad de diez mil necesita nueve —
pero el coloreado se calcula entero en el tipo propio de la compilación, como
las fórmulas.

## Zoom

**Zoom por selección** (Vista → Zoom) sustituye el zoom continuo por un
rectángulo: arrastra con el botón izquierdo y suelta para ir allí, como hace
Ultra Fractal. Solo para acercar. El rectángulo se ajusta a la relación de
aspecto de la ventana alargando su lado corto, nunca recortando el largo, y
funciona correctamente con el fractal girado. El modo Julia rápido toma el
arrastre mientras está activo, y lo devuelve al salir.

**Pasos fijos**, en el mismo submenú:

| | |
| --- | --- |
| Acercar 2× | `Ctrl` `+` |
| Alejar 2× | `Ctrl` `-` |
| Acercar 10× | `Ctrl` `Alt` `+` |
| Alejar 10× | `Ctrl` `Alt` `-` |

Escalan respecto al centro de la vista, así que se deshacen exactamente unos a
otros.

## Girar y centrar la vista

**Fractal → Vista** fija el centro, el radio y el ángulo con números. Bajo el
ángulo tiene dos herramientas para esos mismos campos. Rellenan los campos y no
hacen nada más: Aceptar aplica lo que muestran como un solo paso de deshacer,
Cancelar nada, y una posición guardada nunca sabe de ellas.

**Porciones** divide una vuelta completa en tantas iguales, y **−** y **+**
giran el ángulo hacia atrás o hacia delante una: 60 grados con seis, 360 con
una. Con ninguna, los dos quedan desactivados. Los giros se cuentan desde el
ángulo del que partieron, así que seis sextos son una vuelta entera hasta la
última cifra; escribir un ángulo, o cambiar las porciones, vuelve a empezar la
cuenta desde ahí. El número de porciones se conserva mientras XaoS está
abierto, para la próxima vez que se abra el diálogo.

**Centrar en el eje** mueve el centro en línea recta a través de la imagen
(**R**) o en línea recta hacia arriba (**I**) hasta que el origen queda en
medio en ese sentido, con el ángulo del diálogo. Sin giro, es la parte real o
la imaginaria la que pasa a cero y la otra se conserva. Con giro, es lo mismo a
lo largo de los ejes propios de la imagen: R a 90 grados anula la parte
imaginaria, y R y luego I llevan el origen al medio con cualquier ángulo. Un
centro que ya está ahí se deja como está, y una parte que el giro deja como
redondeo se toma por cero, así que el campo dice 0 y no 1e-17.

El giro en sí se calcula en el tipo propio de la compilación, pi incluido. Se
calculaba a partir de un pi double, lo que en la compilación quad significaba
también un seno y un coseno double: un par que es un giro solo hasta dieciséis
cifras, así que el centro de una vista girada se dibujaba a 5e-17 de su
distancia al origen lejos de donde estaba. Por debajo de 1e-16, que es para lo
que está quad, una vista girada 30 grados mostraba otro lugar completamente
distinto. Un cuarto de vuelta es ahora exacto. En la compilación normal una
imagen girada es la misma salvo píxeles en el límite de su precisión: 52 de los
64 ejemplos con un ángulo no cambian, los demás mueven píxeles sueltos donde el
detalle es más fino que la aritmética.

## La misma imagen con cualquier número de hilos

Una imagen salía de vez en cuando con unos pocos píxeles distintos cuando se
dibujaba con varios hilos: `examples/Hubicka/julia3.xpf` renderizado con cuatro
hilos en una máquina ocupada difería en 26 de 3000 renders, hasta en 37 píxeles
a lo largo de sus filamentos, y nunca con un hilo.

Los píxeles venían del relleno por estimación. Rellena un píxel a partir de las
filas y columnas ya calculadas a su alrededor cuando coinciden, y en modo
Mandelbrot usa los mismos vecinos para decidir si merece la pena comprobar la
periodicidad de un píxel. La estimación avanza por pasadas — primero una fila
de cada seis, luego las filas intermedias — y las pasadas corrían en todos los
hilos a la vez sin nada entre ellas: un hilo pasaba a la pasada siguiente
mientras los demás seguían en la anterior, y una fila se marcaba como hecha
antes de calcularse. Así que si un vecino contaba como calculado dependía de lo
lejos que hubieran llegado los otros hilos.

Ahora cada pasada tiene dos fases, primero sus filas y luego sus columnas, y
cada hilo termina una antes de que ninguno empiece la siguiente; lo que calcula
una fase no se usa como vecino dentro de ella. Los vecinos son entonces los
mismos sea cual sea el número de hilos, y los 3000 renders salen todos iguales.
Con un hilo el orden es el de siempre.

Los filtros también salen iguales. **3D** dibuja una franja de columnas por
hilo, y el sombreado de una columna se lleva desde la columna dos a su derecha
en la fila de abajo, así que las dos últimas columnas de una franja leen las
dos primeras de la franja a su derecha; antes cada franja leía ahí copias de su
propia última columna, y la imagen cambiaba a lo largo de los bordes de las
franjas con el número de hilos. Ahora una franja espera, fila a fila, a lo que
la franja de su derecha ha dibujado ahí, y las franjas siguen avanzando lado a
lado. **Campo de estrellas** sacaba sus estrellas de un único generador
aleatorio que todos los hilos hacían avanzar a la vez; ahora cada fila empieza
en el punto de la secuencia al que habría llegado un solo hilo dibujando toda
la imagen, así que las estrellas son las mismas con cualquier número de hilos.
Siguen cambiando de un arranque de XaoS al siguiente, como siempre: el
generador toma su semilla del reloj.

**El antialiasing** nunca fue la causa, aunque tres ejemplos dibujados con él —
`circle`, `helloween` y `warriormask` — salían distintos de un render a otro.
Lo que difería era el mensaje con el nombre de la fórmula: desaparecía tras
cinco segundos del reloj de pared, y un render ralentizado por una máquina
ocupada lo dibujaba en menos fotogramas, o en ninguno. En un render un mensaje
dura ahora cinco segundos del tiempo propio de la animación, sea cual sea la
máquina; en la ventana no cambia nada.

De los 169 ejemplos renderizados con un hilo, 167 dibujan al píxel como los
dibujaba la 1.8.0: `galaxy` toma sus estrellas del reloj, y `heart` conserva
ahora un mensaje que la 1.8.0, cronometrándolo con el reloj de pared, había
quitado antes de dibujar la imagen. Renderizados de ocho en ocho en cuatro
núcleos, una vez con un hilo y tres veces con cuatro, todos dibujan igual salvo
`galaxy`.

## Relación de aspecto

**Vista → Relación de aspecto** cambia el tamaño de la ventana una vez, para
que la imagen tome la forma elegida: 1:1, 4:3, 3:2, 16:10, 16:9, 2:1, 21:9, o
9:16 en vertical. La ventana no queda sujeta a ella: luego se puede arrastrar a
cualquier forma, y la marca del menú está en la forma que tenga la imagen en
ese momento, o en ninguna.

La imagen conserva el número de píxeles que tenía y cambia solo su forma — o,
cuando no cabría en la pantalla, pasa a ser la mayor de esa forma que cabe.
Desde pantalla completa o una ventana maximizada vuelve primero a una ventana
normal, aún con los píxeles que tenía, así que ahí 16:9 es la ventana 16:9 más
ancha que cabe en la pantalla. La forma se cumple al píxel, también en una
pantalla escalada: al 125% un número de píxeles de cada cinco no se puede
tener, y el tamaño se elige entre los que sí.

Para lo que sirve es para el render. Una posición se dibuja de modo que toda la
vista que contiene esté en la imagen, y una imagen de otra forma muestra más
plano por su lado más largo. Así que una vista encuadrada en la ventana normal
4:3 y renderizada a 1920×1080 tiene un tercio más por los lados de lo que hubo
nunca en la pantalla; encuadrada en una ventana 16:9, el render es lo que se
vio.

## Antialiasing en luz lineal

Filtros → Antialiasing en luz lineal, o `-linearaa`. El antialiasing normal
promedia los colores de las muestras tal como están almacenados, lo que
oscurece los bordes porque los valores almacenados no son proporcionales a la
luz. Este convierte a luz lineal, promedia ahí, y vuelve a convertir, así que
un borde conserva el brillo de aquello entre lo que está.

## Idiomas

Vista → Idioma cambia XaoS a otro idioma, tras reiniciar. Seis están completos:
alemán, español, francés, italiano, portugués y ruso — los menús, los diálogos
y los mensajes, la referencia de fórmulas de usuario, y los diálogos Aceptar,
Cancelar, de archivos y de colores que vienen del propio Qt. Estos últimos van
ahora incluidos, mientras que antes se buscaban donde está instalado Qt y en el
idioma del sistema fuera cual fuera el elegido. Los demás idiomas que XaoS
heredó llegan hasta donde los llevaron sus traductores, y más allá están en
inglés.

En la referencia de fórmulas de usuario las descripciones están traducidas y
los nombres y las listas de argumentos no: son la sintaxis en que se escribe
una fórmula.

Esta guía está traducida a los mismos seis idiomas, junto a la guía inglesa
`ms-xaos-guide.md` como `ms-xaos-guide.de.md`, `.es.md`, `.fr.md`, `.it.md`,
`.pt.md` y `.ru.md`, y Ayuda abre la del idioma en que está XaoS. Cada una dice
de qué versión de la guía inglesa se hizo; donde difieren, la referencia es la
inglesa.

## Compilación

    cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build

Los dos ejecutables acaban en `bin/`. Opciones:

| | |
| --- | --- |
| `-DQUAD=OFF` | omite el segundo ejecutable |
| `-DDEEPZOOM=ON` | hace quad el ejecutable principal, y no compila un segundo |
| `-DDIAG=ON` | compila también `XaoS-diag`, que informa en `xaos-diag.txt` del orden de superposición y de cálculo |

`cmake --build build --target deploy` llena `bin/` con el runtime de Qt y del
compilador para que funcione en una máquina sin Qt.

## Pruebas

    ctest --test-dir build

351 pruebas: el analizador de fórmulas en ambas precisiones, la exactitud de
`gamma`, `lambertw`, `erf` y `randsc` frente a referencias exactas donde
existen, los bucles de iteración frente a sumas de verificación registradas,
las funciones de velocidad del coloreado frente a lo que daban antes y a
aquello para lo que sirven, las posiciones guardadas que van y vuelven y cada
render que termina limpiamente, el guardado y restauración de la superposición,
el historial de deshacer, el giro y el centrado del diálogo Vista, la
referencia de la ayuda frente a la propia tabla del analizador, y el `poly`
aleatorio de su pestaña Utilidades frente al polinomio escrito término a
término.

Una de ellas mide velocidad en lugar de respuestas. Medir en nanosegundos
afirmaría en qué máquina se ejecuta la batería, así que en su lugar afirma una
proporción: lo que cuesta una degradación a las funciones de ruido, contado en
multiplicaciones. Debería ser dos, una por componente, porque el tamaño se
multiplica por la degradación una vez por pasada; elevarla en cambio a la
potencia de la pasada da seis o siete y crece con el límite de iteraciones. La
prueba falla con las dos implementaciones a las que esta sustituyó.
