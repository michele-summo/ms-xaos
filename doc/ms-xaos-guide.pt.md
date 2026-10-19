# MS XaoS

> Tradução do [guia em inglês](ms-xaos-guide.md) da versão 2.7. Onde diferirem,
> a referência é o original em inglês.

Um fork do [XaoS](https://github.com/xaos-project/XaoS) 4.3.3. Tudo o que o
original faz, continua a fazer; este guia descreve o que foi acrescentado ou
mudado, e porquê. Versão 2.4.

## Dois executáveis

`XaoS.exe` calcula com 64 bits de mantissa, `XaoS-quad.exe` com 113. Ambos são
compilados a partir das mesmas fontes — a precisão é um tipo escolhido na
compilação, não uma definição — e o quad chama-se **MS XaoS Quad**, para que os
dois se distingam enquanto estão abertos.

O quad faz zoom cerca de quinze ordens de grandeza mais fundo. Ao renderizar a
mesma posição a 1e-18, o executável normal mantém 21% do detalhe e o quad 92%.
Não é o predefinido porque `__float128` é emulado em software: cerca de oito
vezes mais lento nas fórmulas incluídas e quinze nas definidas pelo utilizador.

Uma posição guardada pelo executável quad leva `(precision 113)`. Abri-la com
uma compilação de menor precisão avisa e desenha-a na mesma; abri-la com uma
versão que não conhece o comando é recusado, e é essa a intenção — uma imagem
calculada a 113 bits não pode ser reproduzida a 64, e não a abrir é mais
honesto do que desenhar outra coisa.

## Os números nos diálogos

Um campo que pede um número aceita um escrito de qualquer forma legível, 1e-18
tal como 0.5. Os diálogos de campo único usavam uma caixa giratória, que guarda
duas casas decimais se não lhe disserem outra coisa, por isso 0.01 sobrevivia e
0.001 passava a nada — uma convergência de Newton de um milionésimo não se
podia introduzir de todo.

O que é mostrado é o texto mais curto que se volta a ler como o número exato
guardado, por isso um milionésimo de milionésimo aparece como 1e-18 e não como
9.99999999999999999978e-19, que é o que imprimir todos os algarismos da
compilação faz de um decimal que um float binário não pode representar
exatamente. Uma coordenada que precisa dos vinte e um algarismos continua a
recebê-los: precisar deles é o mesmo que não voltar a ler-se inalterada sem
eles.

A forma exponencial usa-se a partir de um bilião para cima e abaixo de cem mil
milionésimos, e a decimal simples entre os dois — assim vinte é vinte, cem mil
milhões escreve-se por extenso, e 1e-8 não. Deixada a si própria, a impressão
escolheria comparando o expoente com o número de algarismos pedidos, o que não
serve quando os algarismos são escolhidos pela brevidade: basta um algarismo
para voltar a ler-se como vinte, e vinte saía como 2e+01.

## Posições guardadas

A vista é escrita com a precisão completa da compilação — 21 algarismos
significativos normalmente, 36 no quad. Antes era escrita com um número de
algarismos derivado da profundidade do zoom, que ficava aquém e parava nos 17,
por isso uma posição profunda reabria a uma fração de píxel de onde tinha sido
guardada. Com muitas iterações é uma imagem visivelmente diferente.

Tudo o que uma posição pode conter sobrevive agora inalterado a guardar e
carregar: todos os 169 exemplos incluídos e 141 posições geradas, cobrindo cada
fórmula, plano, modo de coloração e filtro, fazem a ida e volta píxel a píxel.

## Anular

Editar → **Anular** (`U`, ou `Ctrl` `Z`) recua um passo e **Refazer** (`Ctrl`
`Y` no Windows) avança de novo. Um passo é qualquer coisa que mude a imagem: um
comando — uma fórmula, uma paleta, um modo de coloração, um diálogo, um zoom
fixo —, um clique do rato com o zoom, a deslocação ou a rotação que põe em
marcha, um zoom por seleção, e as iterações definidas com as setas. Um comando
que não muda nada, como Recalcular, não é um passo, e anular não recua para lá
da imagem com que o XaoS abriu, depois da configuração e da linha de comandos.
Cada passo é guardado como uma posição guardada, por isso anular repõe tudo o
que uma posição contém, tão exatamente como reabri-la; guardam-se os
últimos 256.

Falhava de duas maneiras. O rato guardava o seu passo ao soltar, enquanto o
zoom ainda estava a abrandar, por isso anular trazia uma imagem a meio caminho;
agora o passo guarda-se onde a vista para. E cada comando guardava um passo,
mudasse alguma coisa ou não, enquanto anular recuava dois partindo do princípio
de que tinha guardado o último ele próprio — depois de um Recalcular trazia de
volta a imagem que já estava no ecrã. Agora uma imagem só se guarda onde difere
da última guardada, e anular e refazer movem-se um passo ao longo delas. Um
passo novo depois de anular corta o que se poderia ter refeito, como num
editor.

O que se moveu sozinho desde o último passo — o piloto automático, uma rotação
contínua — guarda-se como passo próprio quando se carrega em Anular, para que
Refazer volte a ele.

Anular encontrou um erro nas posições. Carregar uma com uma paleta incluída
depois de uma paleta personalizada deixava a paleta personalizada marcada como
em uso, por isso a posição voltava a ser guardada como paleta personalizada:
não o que tinha sido carregado, e para anular não o passo a que tinha voltado,
o que custava o refazer. Agora uma paleta incluída limpa a marca, como uma
gerada sempre fez.

## Fórmulas do utilizador

**Referência.** Ajuda → Referência das fórmulas do utilizador lista todas as
funções, variáveis e notações que a linguagem de fórmulas aceita, numa janela
que se pode deixar aberta ao lado da fórmula que se está a escrever. Está
agrupada pelo propósito das funções, com as funções especiais — `erf`, `gamma`,
`lambertw` — numa secção própria, já que nenhuma delas é elementar e não se
recorre a elas por acaso. Um quarto separador lista os números que aparecem
como argumentos e significam algo em particular — contra que forma uma
armadilha mede, com que espelho o caleidoscópio dobra — em vez de os enterrar
na descrição da função que os recebe. Um quinto mostra uma imagem de cada
pavimentação que `randsctile` aceita, já que quarenta e cinco números com uma
linha de palavras cada um dizem pouco das formas.

Um sexto, **Utilitários**, escreve fórmulas em vez de as descrever. Por agora
contém um `poly` aleatório, a partir de quatro campos:

| | |
| --- | --- |
| comprimento mínimo e máximo | m, o número de coeficientes, sorteado entre os dois; o grau é m − 1. São precisos ambos, como números inteiros de 1 a 254 — até onde vai `poly`, já que o parser aceita 255 argumentos por chamada e `z` é um deles |
| módulo mínimo e máximo | o tamanho de cada coeficiente, sorteado uniformemente entre os dois. 0 e 1 quando deixados vazios; no máximo 1000000 |

Cada coeficiente é complexo, `{re,im}`, com o seu argumento sorteado ao acaso
em toda a volta, e cada parte é escrita com quatro algarismos significativos do
seu módulo, para que um intervalo de 0.001 a 0.002 não desapareça no
arredondamento. Não há `c` na fórmula. Uma fórmula do utilizador faz `z`
começar no píxel tanto no modo mandelbrot como no modo julia, por isso o que
desenha é o conjunto de Julia preenchido do próprio polinómio, igual nos dois;
`c` vai onde se quiser, como em `poly(...) + c`.

Cada alteração num campo gera uma fórmula nova, tal como **Gerar outro**, e
**Copiar** põe-na na área de transferência para Fractal → Fórmula do
utilizador. Enquanto um comprimento estiver vazio, ou um campo não for um
número, estiver fora do seu intervalo ou tiver o mínimo acima do máximo, não há
fórmula, apenas uma linha a dizer porquê.

Cada linha da lista de funções diz o que a chamada recebe e por que ordem, com
nomes que dizem para que serve cada lugar — "base, exponent", "z, low, high",
"before, after, n" — em vez de letras, que não diziam qual lugar era qual; os
mesmos nomes são os que a linha por baixo da barra de fórmulas mostra. Uma
convenção atravessa-os todos: `z` é o valor sobre o qual a função trabalha,
seja o que for que lá esteja escrito. A lista e essas contagens são ambas
verificadas por um teste face à própria tabela do parser, por isso uma função
que ganha ou perde um argumento não pode continuar a ser descrita com os que
tinha. Um argumento que pode ser omitido é mostrado entre parênteses retos com
o valor que toma quando é omitido: `z, [power=1], [angle_factor=1]`.

**Os argumentos separam-se com uma vírgula, e com mais nada.** O ponto e
vírgula também servia: o parser transformava-o em vírgula antes de ler o que
quer que fosse, por isso `LOGN(5;Z^2)+C` e `LOGN(5,Z^2)+C` eram uma só fórmula,
e o mesmo valia para as duas partes de um número complexo. Isso acabou. Um
ponto e vírgula é agora recusado onde aparece — *Operador inválido: ;* — e a
vírgula é o separador em todo o lado, como é no resto da sintaxe.

### A sonda da paleta

**Vista → Sonda da paleta (clique num ponto)** transforma o ponteiro numa
sonda: clique em qualquer ponto do fractal e a linha de mensagens diz **de que
lugar da paleta esse ponto foi desenhado**, como um número e não como uma cor —
*Posição na paleta 12.375*, e mais nada.

O editor de paleta dispõe **31 cores**, e a paleta é construída passando de
cada uma à seguinte **oito células de cada vez**. Por isso o número é a célula
a dividir por oito: `12.000` é exatamente a décima segunda cor do editor,
`12.375` está a três oitavos do caminho entre a décima segunda e a décima
terceira. Uma célula é um oitavo de cor, mais fino do que alguém lê no editor.

* Com a **paleta personalizada** em uso, o número é trazido para **`0` a
  `31`**, que é a numeração do próprio editor.
* Com qualquer outra paleta é **deixado onde cai e vai muito além de 31** — uma
  paleta de 65534 células chega a 8191. O número fora do intervalo é o que diz
  que a paleta não é a que o editor desenha.
* Um ponto que tomou a **cor interior** dá *Posição na paleta: nenhuma*: essa
  cor é um único tom liso que nenhuma célula da paleta representa.

Serve menos para saber "que cor é esta" do que **que parte da paleta a imagem
usa realmente**, e a resposta muitas vezes dá que pensar. Com a coloração
`iter` simples a 60 iterações as células vão de 1 a 59 — de `0.125` a `7.375` —
por isso a imagem nunca toca nas cores 8 a 30 do editor e editá-las não muda
nada. Uma coloração que se espalha, `iter + real` por exemplo, percorre em vez
disso o anel inteiro, várias vezes.

Lê o ponto e não o ecrã, por isso é exata: a mesma cor pode estar em duas
células, e uma imagem em cor real mistura duas células vizinhas, por isso
nenhuma cor designa uma célula própria. Tal como o zoom por seleção, a sonda
devolve o ponteiro ao modo Julia rápido enquanto este está ativo.

**A barra.** Fractal → Fórmula do utilizador abre um diálogo com cerca de
metade da largura de um ecrã de 1280×720 — uns seiscentos e cinquenta píxeis,
menos num ecrã mais pequeno — cuja barra passa para uma linha nova quando a
fórmula não cabe numa. Muda de linha **entre nomes, nunca dentro de um**: um
nome ou um número (`randsctile`, `z_p3`, `0.6`) mantém-se inteiro, e uma linha
pode acabar em qualquer outro sítio — depois de um parêntese, de uma vírgula,
de um operador. A quebra de linha do próprio Qt não servia para isto: quebra
onde as regras da prosa o permitem, o que numa fórmula é depois de algumas
vírgulas e pouco mais, e depois corta um nome onde for preciso. Só um nome mais
largo do que a barra inteira é cortado. A barra cresce uma linha de cada vez
com a fórmula, e quando o diálogo já não caberia num ecrã de 720p deixa de
crescer e passa a deslizar. Para tudo o resto continua a ser uma fórmula numa
linha: Enter é OK e não uma linha nova, e numa fórmula colada as quebras de
linha e tabulações passam a espaços, que o parser lê como nada.

**Histórico...** por baixo da barra abre as fórmulas usadas antes — as dez mais
recentes, a mais nova primeiro, cada uma uma vez — que a antiga barra de uma
linha guardava na sua lista pendente. Uma fórmula demasiado comprida para a
janela aparece com o seu início, três pontos e o seu comprimento,
`poly(z, {0.01097,0.00699}, {-0.05... (157 carateres)`, para que duas que
começam da mesma maneira ainda se possam distinguir; o início acaba entre
nomes, como uma linha da barra, e vê-se mais quanto mais larga for a janela.
Faça duplo clique numa, ou selecione-a e prima Usar, e ela vai inteira para a
barra. Uma fórmula é acrescentada quando o diálogo é aceite com ela, e só
então: carregar uma posição que contém uma fórmula do utilizador, reproduzir
uma ou renderizar uma com `-render` executa o mesmo comando `usrform` mas já
não a põe no histórico, onde cada uma delas tomava o lugar de uma que tinha
escrito.

**Seta para cima e seta para baixo** percorrem as mesmas fórmulas sem sair da
barra. A seta para cima, a partir da primeira linha da barra, traz a fórmula
usada antes da que se vê; a seta para baixo, a partir da última linha, volta
para a mais recente e depois para o que a barra tinha antes de começar. As duas
pequenas setas ao lado da barra fazem o mesmo, cada uma a cinzento quando não
há mais nada nesse sentido. Uma fórmula que já está na barra é saltada, e o que
escreve passa a ser o ponto de onde se recomeça.

A última fórmula usada aparece esbatida na barra: inteira enquanto a barra está
vazia, e o resto a seguir ao que escreveu enquanto o que escreve for o seu
início. Tab põe-na inteira na barra. Assim que o que escreve segue outro
caminho, o texto esbatido desaparece, e Tab passa para os botões como antes. A
barra abre vazia, com essa fórmula esbatida, e OK com a barra vazia usa-a; só
da primeiríssima vez, antes de se ter usado qualquer fórmula, abre com a
fórmula em uso, `z^2+c`.

**Vista → Estado** mostra a fórmula da mesma maneira: numa linha se couber na
janela, em duas se não — a segunda por baixo da primeira, mudando de linha
entre nomes — e quando nem duas chegam, a segunda acaba em três pontos e o seu
comprimento, `... (918 carateres)`. As linhas de baixo descem para dar espaço,
e tudo é redisposto quando a janela muda de largura. Antes era uma única linha
que saía da borda da imagem.

**Por baixo da barra há uma linha que diz o que recebe a chamada em que está**,
com o argumento que está a escrever a negrito. Ponha o cursor dentro de
`randsc(` e lê-se *randsc(**seed**, [size=1+i], [degradation=0.5+0.5i],
[kaleidoscope=1], [kaled_mirror=0], [skew=0])*; passe uma vírgula e o negrito
move-se consigo. Não precisa de mais da fórmula do que um parêntese aberto com
um nome à frente, por isso está lá enquanto a fórmula está meio escrita, que é
quando é precisa. Ponha o cursor fora de todas as chamadas e a linha fica
vazia.

Seis das posições incluídas estavam escritas à maneira antiga e foram
reescritas com vírgulas: `circle`, `heart`, `helloween`, `pentafrac`,
`warriormask` e `burnship`, todas em **Ficheiro → Abrir → examples/Malczak**.
**Uma posição sua que tenha um ponto e vírgula na fórmula não carrega** até a
fórmula ser escrita com vírgulas — a mensagem indica o caráter, e a solução é
abrir a fórmula e substituí-los.

**Deixar um lugar vazio.** Um argumento mostrado entre parênteses retos pode
ser omitido a meio de uma chamada tal como no fim, não escrevendo nada entre os
dois separadores: `julian(z, ,3)` dá o primeiro e o terceiro e deixa a função
dizer qual é o segundo. Os espaços não fazem diferença — `f(z, ,5)` e `f(z,,5)`
são a mesma chamada — e um argumento que não está entre parênteses retos tem de
ser escrito, por isso `poly( ,1,2)` e `sin(,)` são recusados em vez de
adivinhados.

O que significa um lugar vazio é com cada função. A maioria toma o valor
predefinido que declara; um coeficiente de `poly` deixado vazio é um termo que
não existe, por isso `poly(z,1, ,1)` é `z^2 + 1`; e um ramo de `ifiter` ou
`ifiterl` deixado vazio repete o anterior, que é a forma de dar a um ramo mais
do que uma passagem no ciclo. `ifiter(f(z), , ,g(z), , , , )` executa `f`
durante três passagens e `g` durante cinco, e não custa mais a avaliar do que
escrevê-las por extenso: para onde leva cada escolha decide-se uma vez, quando
a fórmula é analisada.

Uma fórmula do utilizador ainda por escrever diz `z^2+c`. Antes dizia o burning
ship, que é um fractal próprio e um enigma para encontrar como ponto de
partida.

**`ifiterf(rest, last)`** — avalia `rest` em todas as passagens menos a última
e `last` nessa. A última passagem é a última que o limite de iterações permite:
uma fórmula não tem maneira de saber qual será a passagem que escapa, porque
isso depende do valor que ainda não produziu. Só a escolhida é avaliada.

**`ifiterr(before, after, n)`** — avalia `before` enquanto o número da passagem
estiver abaixo de `n` e `after` a partir de `n`. Só a escolhida é avaliada,
como nas outras duas, embora primeiro tenha sido preciso ensinar o limiar ao
parser: um argumento pode agora ser marcado como lido pelo seletor em vez de
escolhido por ele, e é então avaliado antes de a escolha ser feita.

**As variáveis.** `z` é o valor corrente; `c` é o ponto, que é o píxel no modo
mandelbrot e a constante no modo julia; `x` é a coordenada simples, que é o
píxel em ambos os modos e não muda de passagem para passagem; `n` é o número da
iteração. Na *Inicialização do utilizador* `z` é onde `z` teria começado se não
houvesse inicialização, que é o mesmo valor que `x`, por isso uma inicialização
feita só de `z` diz exatamente o mesmo que não dizer nada.

Essa última parte é nova, e também o é poder dizê-la: `z` não estava registado
em nenhum dos dois parsers por onde passa a inicialização, e `x` no errado, por
isso uma inicialização que nomeasse qualquer um deles era recusada pela thread
que tinha de a calcular -- enquanto o diálogo, cujo parser tem os dois, a
aceitava sem dar sinal de que depois era deitada fora. `n` ficava com o que o
píxel anterior lá tivesse deixado, ou seja, uma inicialização que o lia lia a
ordem pela qual os píxeis eram calculados.

**`p1`, `p2` ... `p9999`** — o valor que `z` tinha numa passagem anterior: `p1`
na passagem antes desta, `p2` na anterior a essa, `p9999` nove mil novecentas e
noventa e nove passagens atrás. `p` é outro nome para `p1`. Antes de haver
tantas passagens valem aquilo de que o histórico parte, que *Fractal → Definir
os valores p na primeira iteração* decide: o ponto a iterar quando está ligado,
zero quando está desligado.

Antes paravam em `p6` e eram guardados deslocando um array de seis posições um
lugar a cada passagem — uma cópia por posição por passagem, paga quer a fórmula
nomeasse alguma quer não. É por isso que paravam em seis: com nove mil seriam
320 KB de cópias por passagem, medidos em vinte e quatro vezes o custo total de
iterar `z^2+c`. O histórico é agora um anel, escrito uma vez por passagem, e só
se leem dele as posições que uma fórmula nomeia realmente, por isso `p9999`
custa o mesmo que `p1` e uma fórmula que não nomeia nenhuma não custa nada. A
janela de referência lista-os entre as variáveis.

**Sufixos numa variável.** `z`, `c`, `x`, `p` e de `p1` a `p9999` aceitam
sufixos que valem pelas chamadas feitas mais vezes sobre eles:

| sufixo | vale por |
| --- | --- |
| `_b` | `bship(…)` |
| `_bi` | `bshipi(…)` |
| `_br` | `bshipr(…)` |
| `_pM` | `parchment(…, M)` |
| `_paM` | `parchmenta(…, M)` |

Leem-se da esquerda para a direita, cada um envolvendo o que os anteriores
fizeram: `z_p3` é `parchment(z,3)`, `c_b_p2` é `parchment(bship(c),2)` e
`p12_p2_p3` é `parchment(parchment(p12,2),3)`. *M* faz parte do nome, por isso
escreve-se em algarismos — um número inteiro, de um ou mais; `z_p0` é recusado
com uma mensagem que o diz. Não em `n`, que é uma contagem e não tem uma
segunda componente sobre a qual trabalhar.

São expandidos nas chamadas antes de a fórmula ser lida, por isso uma variável
com sufixo desenha exatamente o que as chamadas escritas por extenso desenham e
custa o que elas custam; a fórmula é guardada como foi escrita. O que não for
um dos cinco é deixado como está e recusado como o nome desconhecido que é.

**Removidos.** `powi`, `powdc` e `logcn` eram segundos nomes de `pow` e `logn`;
`rad`, `deg` e `sign` estavam listados sem nenhuma implementação por trás, e um
segundo `trunc` era ocultado pelo que funcionava. Uma fórmula que use um alias
removido deve usar o nome que resta — a função é idêntica.

**`erf(z)`** — a função erro no plano complexo. Precisa a cerca de três ulp
dentro de um bailout de dois, onde um fractal realmente itera.

**`poly(z, k1, k2, ..., km)`** — um polinómio em `z`:

    k1*z^(m-1) + k2*z^(m-2) + ... + k(m-1)*z + km

O primeiro coeficiente escrito multiplica a potência mais alta e o último fica
sozinho, por isso a chamada lê-se pela ordem em que se diz o polinómio.
Calculado pela regra de Horner, que são m−1 multiplicações em vez das m(m−1)/2
que elevar cada potência separadamente exigiria, e ainda por cima a mais
precisa das duas.

    poly(z,1,0,0)+c            o Mandelbrot, escrito por extenso
    poly(z,1,0,0,{0.7,0.2})    z^3 + 0.7+0.2i

**`randsc(seed, size, degradation, kaleidoscope, kaled_mirror, skew, skew_mode, selfsim)`**
— ruído coerente sobre o ponto, que dá manchas em vez de neve píxel a píxel.
`size` (predefinição `1+i`) é a largura média de uma mancha ao longo do eixo
real e a sua altura ao longo do imaginário. É **ruído de gradiente**, uma
oitava do movimento que os [modos fbm](#movimento-browniano-fracionário) e
`fbm()` somam. Era ruído de valor, e pela mesma razão que o deles desenhava a
sua grelha: as manchas achatavam-se ao longo de cada uma das suas linhas e
saíam quadradas. As manchas são tão largas como eram — a grelha está a três
quintos do tamanho, e o campo deixa de se parecer consigo próprio metade a onze
dezasseis avos do tamanho, antes e depois — e os valores espalham-se como
antes, 0.2145 a partir do meio, com a curva de contraste do fbm a levar a
distribuição mais estreita do ruído de gradiente de volta à do ruído de valor.
Mas cada valor é diferente, e custa 167 ns por chamada contra 147. Os cinco
mosaicos são a sua grelha e não mudam nem um bit. `degradation` (predefinição
`0.5+0.5i`, que as reduz a metade a cada passagem) encolhe-as à medida que a
iteração avança: o tamanho é multiplicado por ele a cada passagem, componente a
componente, por isso `0.5+0.2i` sobre `1+i` dá `1+i` na primeira passagem,
depois `0.5+0.2i`, depois `0.25+0.04i`. Um zero em qualquer componente de
qualquer dos dois argumentos devolve zero em vez de dividir por zero. Depois
vêm os dois argumentos do caleidoscópio, `skew` (predefinição `0`), `skew_mode`
e `selfsim`, cada um descrito mais abaixo. Só a semente é obrigatória.

### A inclinação, e porque estes campos são difíceis de colorir

O motor colore com as duas componentes da órbita: `real` e `zmag` leem uma,
`imag`, `angle` e `real / imag` a outra. Estes campos devolvem um único número
para uma célula inteira, no eixo real, por isso cada modo desenha **um tom por
célula** — `zmag` e `iter + real` parecem o valor truncado, e os outros três
não têm nada para ler. Medido em novecentos píxeis: de 14 a 126 valores para os
dois primeiros, e exatamente **um** para os restantes.

Variar o valor dentro de uma célula não se consegue com um único número real. O
mesmo número decide se o ponto sai, por isso um valor que se move ao longo de
uma célula faz sair metade da célula e deixa dentro a outra metade, e a célula
sai **cortada em duas por uma linha reta**. Tentou-se duas vezes e cortou os
mosaicos das duas vezes.

`skew` fá-lo antes com a segunda componente, e **roda** o valor em vez de o
escalar. Quanto roda segue **o quanto o ponto está fora do meio da sua
célula**: com `out` a ir de zero no meio a um na borda, o valor é rodado o
dobro do arco-tangente de `skew_re*out + skew_im`.

`out` é medido **na geometria da própria célula**, por isso o que a cor desenha
dentro de uma célula é a forma em que o campo está dividido, não um conjunto de
linhas estendidas por cima — zero no meio de uma célula, um ao longo de toda a
sua borda:

| campo | o que a cor desenha dentro de uma célula |
| --- | --- |
| `randsc` | as suas próprias manchas — a rotação segue o nível, por isso os contornos são os do campo |
| `randscq` | quadrados, em torno do meio do quadrado |
| `randsch` | hexágonos, em torno do meio do hexágono |
| `randsct` | triângulos, em torno do meio do triângulo |
| `randscp` | o seu próprio polígono, encolhido passo a passo — não há duas células com a mesma forma |

A **parte real** da inclinação é esse gradiente ao longo de uma célula. A
**parte imaginária** é uma rotação plana que toda a célula partilha: desloca a
cor de uma célula sem desenhar nada dentro dela, por isso `0.4+0.2i` dá ao
mesmo tempo os contornos da forma e um deslocamento por célula.

### O que a inclinação desenha: `skew_mode`

O sétimo argumento diz **o que varia dentro da célula**. Os números
**somam-se** em vez de se excluírem, por isso uma chamada pode pedir dois e
obter ambos.

| | |
| --- | --- |
| **1** (e `0`) | **a forma** — o próprio contorno da célula, encolhido passo a passo. O que recebe uma chamada que não dá `skew_mode`, e o que se usa quando não se nomeia nem 1 nem 2 |
| **2** | **a roseta** — o ângulo em torno do meio da célula, dobrado em tantas voltas quantas cunhas tem o caleidoscópio, para que cada célula traga a simetria da própria imagem. Raios simples quando nada está dobrado |
| **4** | **ao longo da cunha** — uma rotação que segue onde o ponto está na sua cunha do caleidoscópio, igual em cada cunha |
| **8** | **radial** — o módulo move-se tal como o ângulo |

`3` é a forma e a roseta ao mesmo tempo, o que faz uma espiral. `11` é essas
duas com o módulo também a mover-se.

**`4` segue a geometria do próprio caleidoscópio.** Roda o valor nada nas
bordas de cada cunha e toda a rotação que a parte imaginária dá por si só — o
dobro do seu arco-tangente, 127 graus para `2` — ao longo da linha pelo seu
meio, e entre os dois a rotação da parte imaginária vezes o quanto se
atravessou. Cada cunha recebe o mesmo, por isso um caleidoscópio de seis
continua a ser um caleidoscópio de seis, espelhos incluídos, e a rotação é
contínua através das bordas e dos meios, que é onde a dobra espelha. Não faz
nada quando nada está dobrado, ou quando a parte imaginária da inclinação é
zero. Sem trigonometria, não custa nada a 64 bits de mantissa e um décimo de
chamada a 113.

Ao princípio era **uma rotação diferente em cada cunha** — a parte imaginária
vezes o número da cunha mais um — pensada para que as cópias deixassem de ser
idênticas. É a única coisa que as cópias de um caleidoscópio não podem ser: uma
imagem de simetria seis saía como seis fatias sem relação com um corte reto em
cada junção. O múltiplo também esbarrava na meia volta, por isso quatro fatias
de seis ficavam quase no eixo real, e aí `real / imag` desenhava-as como neve
por muito que `selfsim` fizesse a média do ruído.

O que resta disso é uma propriedade de `real / imag`, que divide pela parte
imaginária de *z*: onde quer que uma rotação a leve através de zero, a
coloração salta de uma ponta da paleta à outra, e os últimos restos das
passagens finas decidem o lado píxel a píxel. Uma rotação que varre a cunha
atravessa aí na maioria das definições. Num `randsct` de simetria seis colorido
por `real / imag`, com `selfsim` a um, uma parte imaginária de `0.25` não
deixou grão nenhum e `2` deixou faixas junto às bordas; uma coloração que não
divide por uma componente não tem nenhum.

**`8` é o que responde a `zmag`.** Uma rotação deixa o módulo onde está, e
`zmag` e o bailout leem o módulo e mais nada — por isso nenhum outro modo os
consegue alcançar. O preço é aquele que a inclinação foi pensada para evitar: a
fuga move-se com o módulo, por isso a figura é cortada onde quer que o bailout
caia dentro de uma célula. Combine-o com um bailout circular e uma coloração
que leia as componentes separadamente, ou aceite o corte de propósito.

Uma inclinação nula deixa cada modo sem efeito, por isso `skew_mode` sozinho
não muda nenhuma imagem.

Primeiro foi escrita uma rampa reta — `skew_re*du + skew_im*dv` ao longo da
célula — e é isso que isto substitui. Coloria, mas desenhava as mesmas faixas
diagonais, verticais ou horizontais em todas as células por igual, fosse qual
fosse a divisão do campo, e as faixas estragavam a imagem.

Em **`0`**, que é o que recebe uma chamada que não o nomeia, a rotação é nula,
a parte imaginária fica em zero, e cada valor é o número que sempre foi **até
ao último bit** — as somas de verificação de antes de o argumento existir
passam inalteradas, e uma posição guardada renderiza-se idêntica. A posição
`talc`, três chamadas de ruído e um `c` sobre 345600 valores, dá a mesma
assinatura que a compilação anterior ao argumento.

**Quanta cor dá.** O índice do motor é `(iter + quantity) * speed + shift`, por
isso uma variação de um na grandeza vale uma faixa a uma velocidade de um. A
variação de `imag` dentro de uma célula:

| inclinação | 0.02 | 0.05 | 0.1 | 0.2 | 0.4 | 0.6 | 1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `randsc`, `randsct` | 0.08 | 0.21 | 0.43 | 0.83 | 1.5 | 1.9 | 2.1 |
| os três mosaicos | 0.06 | 0.14 | 0.28 | 0.53 | 0.96 | 1.2 | 1.4 |

Por isso **de `0.3` a `0.6` é onde aparece uma faixa a uma velocidade de um**,
que é o que faltava a estes campos. Contar valores distintos diria que `0.01`
chega — novecentos — mas mil valores dentro de um centésimo de faixa são todos
uma só cor.

**O que custa.** Escalar movia o módulo do valor, e um bailout redondo olha
exatamente para isso: uma célula cujo nível ficava perto de onde cai o bailout
saía cortada em duas, e os mosaicos perdiam as suas formas. Uma rotação deixa o
módulo onde está, por isso **com um bailout circular a inclinação é grátis com
qualquer intensidade** — em 90000 píxeis, a 0.05, 0.3 e 1, e nos cinco campos:
nem um píxel sai de forma diferente.

Um polígono de bailout é outra coisa, porque não lê o módulo mas as
componentes, e uma rotação leva o valor por um círculo que pode atravessar um
lado. O quanto o pode atravessar é o quanto os cantos da forma sobressaem para
lá dos seus lados, e os números seguem isso — a fração da imagem cuja fuga
muda, nos cinco campos:

| forma do bailout | cantos sobre apótema | inclinação 0.01 | inclinação 0.05 | inclinação 0.3 | inclinação 1 |
| --- | --- | --- | --- | --- | --- |
| círculo | 1.00 | **0%** | **0%** | **0%** | **0%** |
| hexágono | 1.15 | **0%** | 0.0–0.8% | 0.7–2.2% | 0.4–2.4% |
| quadrado | 1.41 | **0%** | 0.0–1.0% | 0.8–2.8% | 1.2–12.8% |
| triângulo | 2.00 | 0.05–1.3% | 0.7–5.0% | 2.6–6.4% | 7.4–21.6% |

O **número** para que o bailout está definido também o move, e pela mesma
razão: o que pode mudar de lado é um valor cujo módulo cai entre o apótema do
polígono e os seus cantos, e onde esse anel fica entre as células é o que o
número decide.

Por isso, com um polígono, os dois desejos puxam em sentidos opostos, e não há
maneira de arranjar isto de outra forma: a fuga lê as duas componentes do valor
e a coloração também, por isso tudo o que dê à cor algo para ler é algo que a
fuga pode ver. **Combine a inclinação com um bailout circular** e não custa
nada com qualquer intensidade. Com um polígono, ou aceita os poucos por cento
ou mantém a inclinação em `0.05` e sobe antes a velocidade de coloração.

O que uma rotação não pode tocar é **`zmag`**, que lê o módulo. Colora com algo
que leia as duas componentes separadamente: `iter + real`, `iter + imag`,
`angle`, `real / imag` no exterior, e `real` ou `real / imag` na coloração
interior.

Daqui seguem-se duas coisas. O valor é complexo enquanto a inclinação não for
zero, por isso `randsc(7,,,,,0.4)*z` roda `z` além de o escalar. E o
caleidoscópio dobra de qualquer maneira: medido em cinco campos, duas, três,
cinco e seis cunhas e ambos os espelhos, uma rotação de uma cunha deixa cada
valor onde estava.

**`randscq(...)`** — o mesmo campo sem a interpolação: um mosaico de células
quadradas planas em vez de manchas. Os mesmos argumentos, o mesmo significado.

**`randscp(...)`** — de novo o mesmo campo, sem as curvas mas não sem a
irregularidade: uma semente é espalhada dentro de cada célula e cada posição
toma o valor da semente mais próxima, o que desenha um diagrama de Voronoi —
polígonos convexos planos de lados retos, não há dois com a mesma forma.

**`randsch(...)`** — outra vez o mesmo em hexágonos: um favo de células planas.
Dos polígonos regulares que pavimentam o plano é o que não tem uma direção
preferida — cada célula tem os mesmos seis vizinhos nos mesmos seis ângulos —
por isso lê-se como um material e não como uma grelha.

**`randsct(...)`** — triângulos equiláteros, de orientação alternada. Um
mosaico triangular tem uma direção preferida, e é isso que se pede ao
escolhê-lo.

**`randsctile(tiling, seed, ...)`** — o mesmo campo sobre qualquer uma de
**quarenta e cinco pavimentações**, que o primeiro argumento escolhe; tudo o
que vem depois é o que a família recebe, um lugar mais à frente, por isso
`selfsim` é o nono. Cada pavimentação está escalada para um ladrilho de área
unitária em média, por isso o primeiro argumento muda a forma das células e não
a sua escala. Um número fora de 1 a 45 desenha zero.

| | |
| --- | --- |
| 1–3 | as pavimentações regulares: quadrados, triângulos, hexágonos |
| 4–11 | as oito arquimedianas, polígonos regulares com todos os vértices iguais: 4.8.8 (octógonos e quadrados), 3.6.3.6, 3.4.6.4, 3.12.12, 4.6.12, 3.3.3.4.4, 3.3.4.3.4 (quadrado achatado), 3.3.3.3.6 (hexagonal achatada) |
| 12–19 | as suas duais, que não são regulares: tetrakis quadrada, rombille (os cubos), deltoidal tri-hexagonal (papagaios), triakis triangular, kisrombille, e os pentágonos prismáticos, do Cairo e em flor |
| 20–32 | tijolos, aparelho flamengo, espinha de peixe, entrançado de cesto, pitagórica (quadrados de dois tamanhos), chevrons, quadrados e losangos, casinhas (um pentágono, em pé e pendurado alternadamente), filas de quadrados e triângulos em dois ritmos, hexágonos entre triângulos, cruzes gregas, tetraminós em T |
| 33–37 | estrelas islâmicas: de oito pontas com cruzes, de seis pontas com hexágonos, de oito pontas a partir de 4.8.8, de doze pontas a partir de 3.12.12 e de 4.6.12 |
| 38 | células de Voronoi, irregulares — as de `randscp`, com um valor próprio |
| 39–43 | pavimentações que nunca se repetem, das multigrelhas de de Bruijn: os losangos de Penrose, os papagaios e dardos de Penrose, os quadrados e losangos de Ammann–Beenker, e losangos em doze e em sete direções |
| 44–45 | por substituição: o cata-vento, cujos triângulos estão virados em todas as direções possíveis, e a cadeira, L cortados em L |

O **separador Pavimentações** de Ajuda → Referência das fórmulas do utilizador
tem uma imagem de cada uma, com a legenda com que o separador Valores a lista,
e um duplo clique numa copia o início da sua chamada — `randsctile(18, ` — para
a área de transferência, pronto para a semente. As imagens são desenhadas pelo
próprio `randsctile`, oito unidades de lado em torno da origem, por
`tools/randsctile-thumbnails.cpp`, e guardadas no executável em vez de serem
desenhadas sempre que se abre o separador. Imagens guardadas podem ficar
desatualizadas, por isso cada uma é guardada com uma impressão digital do que a
função desenhava, e um teste volta a desenhar os mesmos pontos e falha, dizendo
qual pavimentação, até
`cmake --build <build directory> --target randsctile-thumbnails` a ter
desenhado de novo.

As periódicas são **tabelas**, escritas por `tools/randsctile-tables.py`, que
constrói cada pavimentação a partir da sua geometria e a verifica — vinte mil
pontos ao acaso, cada um dos quais tem de cair em exatamente um ladrilho —
antes de a escrever. As estrelas islâmicas são feitas a partir das
pavimentações arquimedianas: cada polígono de oito lados ou mais torna-se uma
estrela com as pontas nos vértices do polígono, e o triângulo fino entre cada
reentrância da estrela e o lado passa para o ladrilho do outro lado desse lado,
por isso duas estrelas que partilham um lado formam um losango entre elas e um
polígono pequeno recebe uma reentrância em cada lado que partilha com uma
estrela. Os ladrilhos podem ser côncavos — estrelas, cruzes, o T — por isso um
ponto é colocado no seu ladrilho pelo número de cruzamentos em vez de lado a
lado. O ponto é dividido em períodos inteiros na precisão de trabalho, como
`randscq` o divide em células, por isso uma pavimentação periódica é exata por
mais longe que esteja.

As que nunca se repetem **não têm tabela** e são calculadas de novo em cada
ponto. Os losangos de Penrose, Ammann–Beenker e os losangos de doze e de sete
direções vêm da **multigrelha de de Bruijn**: famílias de linhas paralelas num
segundo plano, um ladrilho onde quer que duas se cruzem, designado pelas duas
linhas — quatro números inteiros, exatos seja quem for que os peça. O ponto é
colocado experimentando os cruzamentos junto a onde cai nesse plano, primeiro
as famílias cujas linhas passam mais perto, o que o encontra em quatro a oito
tentativas em média. Os papagaios e dardos vêm dos triângulos de Robinson
cortados nível a nível a partir de uma roda de dez em torno da origem, o
cata-vento de um triângulo retângulo cortado em cinco, a cadeira de um L
cortado em quatro; em cada nível uma linha ou duas decidem em que filho está o
ponto.

**O que custa**, medido por chamada face aos 220 ns de `randscq` nesta máquina:
as pavimentações periódicas de 255 a 380 ns, as que nunca se repetem de 410
a 680. São calculadas em `double` nas duas compilações, incluindo o lugar
dentro de um período — em long double as que nunca se repetem demoravam o
triplo — e assim as duas compilações desenharam as quarenta e cinco
pavimentações iguais píxel a píxel, quarenta mil pontos de cada. Para lá de 2³²
células da origem, onde `double` começaria a perder os seus ladrilhos, as
pavimentações que nunca se repetem ficam planas como a família fica para lá da
sua grelha; as periódicas aguentam até onde aguenta o resto da família.

A inclinação mede um ladrilho como `randscp` mede o seu polígono, face ao raio
do maior círculo que o ladrilho contém, por isso os contornos são o ladrilho
encolhido, côncavos incluídos; a roseta roda em torno do meio do ladrilho.

    z+randsctile(18,13,{0.3,0.3},{1,1})*0.25        pentágonos do Cairo
    z+randsctile(36,13,{0.3,0.3},{1,1})*0.25        estrelas de doze pontas
    z+randsctile(39,13,{0.3,0.3},{1,1},5)*0.25      Penrose, dobrado em cinco

**`fbm(value, seed, [intensity=4], [frequency=8], [octaves=4], [roughness=0.5])`**
— um movimento browniano fracionário: ruído de gradiente somado em oitavas,
cada uma com o dobro da frequência da anterior e mantendo `roughness` da sua
altura — o mesmo movimento que os modos de coloração com esse nome desenham,
ver [Movimento browniano fracionário](#movimento-browniano-fracionário). O que
desenha é **desgaste** mais do que um padrão — manchas e amolgadelas — porque
nenhuma oitava é grande o suficiente para se ver sozinha e nenhuma é pequena o
suficiente para desaparecer.

Onde `randsc` lê a posição e só pode ler a posição, **esta lê o que quer que se
escreva à frente dela**. É todo o sentido de a ter como função:

    fbm(z,7)          move-se com a órbita
    fbm(x,7)          fica parada no plano
    fbm(z*3+c,7)      seja isso o que for
    z^2+c+fbm(z,7,0.1)*i    o conjunto com uma borda rugosa

O valor e a semente são obrigatórios; os restantes têm valores predefinidos. O
movimento vai **de zero a `intensity` e nunca abaixo**, por isso somá-lo a algo
não o pode puxar para baixo de zero. A semente é lida como `randsc` lê uma, por
isso o mesmo número significa o mesmo campo em ambos. `octaves` é mantido entre
um e vinte e quatro, e só vale a pena subi-lo quando se sobe `roughness`: com a
metade simples a oitava oitava traz um ducentésimo quinquagésimo do total.

**Não tem caleidoscópio próprio**, ao contrário da família `randsc`, porque já
existe um com que compô-la: `parchment(z, sectors)` dobra o ângulo de `z` em
`|sectors|` setores e `parchmenta` espelha as metades alternadas, por isso
`fbm(parchmenta(z,6),7)` di-lo, e di-lo onde quem lê a fórmula o pode ver.

Os modos de coloração em **Fractal → Coloração exterior → Movimento browniano
fracionário...** aplicam o mesmo movimento a uma cor. Esta põe-no onde uma
fórmula o pode usar.

| | |
| --- | --- |
| `randsc` | suave, curvo |
| `randscq` | duro, regular, quadrado |
| `randscp` | duro, irregular |
| `randsch` | duro, regular, hexagonal |
| `randsct` | duro, regular, triangular |
| `randsctile` | duro, qualquer uma de quarenta e cinco pavimentações |

As seis colocam uma célula sobre cada quadrado unitário do tamanho em vigor,
por isso `size` significa o mesmo em todas: mudar uma letra muda a forma das
células e não a escala da imagem. Um hexágono de raio circunscrito um cobre 2.6
quadrados unitários e um triângulo de lado um cobre 0.43, por isso essas duas
grelhas são escaladas para corresponder às restantes em vez de ficarem como
vêm.

Cada passagem tem um campo próprio. A iteração entra no hash juntamente com a
posição, por isso uma fórmula que chame uma destas uma vez por iteração obtém
um valor novo de cada vez — incluindo com uma degradação de um, que deixa o
tamanho como está. A degradação define o tamanho das manchas numa dada
passagem; não é o que torna a passagem diferente.

Há um limite inferior para o tamanho. Uma célula encontra-se dividindo a
posição pelo tamanho e o resultado tem de caber num inteiro, por isso o tamanho
não pode descer utilmente abaixo de cerca da posição a dividir pelo que cabe
num inteiro — um milésimo de milionésimo de milionésimo dela. Uma degradação de
meio chega aí em cerca de sessenta passagens. Para lá disso não resta estrutura
de células para resolver e o campo é um único valor plano em todo o plano;
continua a mudar em cada passagem e continua a diferir entre as seis funções,
por isso uma fórmula que subtraia uma de outra não fica em zero a iterar até ao
limite por nada. Para um esvanecimento que continue a ser uma imagem até ao
fim, use uma degradação perto de um.

Todas fazem o hash da posição do ponto e da iteração, nunca de `z`, e nunca de
qualquer estado global. Por isso a mesma imagem volta em cada redesenho, com
qualquer número de threads, e no modo Mandelbrot ou Julia por igual. `randsc` é
contínuo, por isso as duas precisões concordam até cerca de 1e-19; os mosaicos
são funções em degrau e podem discordar numa linha fina ao longo das bordas das
células, o que é inerente a pedir bordas duras.

Uma fórmula que chame qualquer uma delas desliga o seguimento da fronteira.
Essa otimização preenche sem a calcular uma região que encontrou rodeada de uma
só cor — verdade para um fractal, falso para um campo de ruído.

O antigo `rand` não mudou e continua a depender de quantas vezes foi chamado,
por isso a mesma posição não se redesenha da mesma maneira. Prefira `randsc`.

### Caleidoscópios

Os dois argumentos depois da degradação dobram o plano antes de o campo ser
amostrado nele. O primeiro é em quantas cunhas o dobrar, e `1` — o predefinido
— deixa-o como está; o segundo, `kaled_mirror`, é qual espelho faz a dobra:

| kaled_mirror | |
| --- | --- |
| `0` (e qualquer outro valor) | a metade afastada de cada cunha espelha a metade próxima, por isso cada cunha é simétrica em relação à sua própria bissetriz |
| `1` | o mesmo ao contrário, com a metade próxima a espelhar a afastada |

    randsc(13,{0.6,0.6},{1,1},6,0)     seis cunhas, cada uma espelho de si mesma
    randsc(13,{0.6,0.6},{1,1},3,1)     três cunhas, dobradas no outro sentido

Ambas as dobras são contínuas onde as cunhas se encontram, por isso o ruído
mantém-se coerente e as duas precisões continuam a concordar; uma dobra que se
encontrasse a si mesma de forma desigual ver-se-ia como uma costura.

É a única parte da família que custa trigonometria, e só um valor de dois ou
mais lá chega. Uma chamada que não diz nada sobre isto paga uma comparação.

### Tirar uma imagem do ruído

`randsc` devolve um valor em `[0, 1)`. Usado sozinho, a iteração nunca sai de
um bailout de 4, por isso nada foge e a imagem é plana. Multiplique-o, ou
deixe-o perturbar uma iteração que foge de facto:

    z+randsc(13,{0.35,0.35},{1,1})*0.25       manchas
    randsc(13,{0.15,0.15},{0.5,0.5})*2.5      movimento browniano
    randscq(13,{0.15,0.15},{0.5,0.5})*2.5     quadrados dispersos, cada vez menores
    z+randscp(13,{0.3,0.3},{1,1})*0.25        polígonos irregulares
    z+randsch(13,{0.3,0.3},{1,1})*0.25        um favo
    z+randsct(13,{0.3,0.3},{1,1})*0.25        triângulos
    z^2+c+randsc(13,{0.25,0.25},{1,1})*1.2    o próprio conjunto deformado

Uma degradação de 0.5 reduz as manchas a metade em cada passagem, por isso ao
fim de vinte iterações são um milionésimo do seu tamanho e ficam abaixo de um
píxel. Para um esvanecimento lento ao longo de muitas iterações, use algo perto
de 1 — `0.97^60` ainda é 0.16.

Uma semente inteira é exata nas duas compilações; uma fracionária é quantizada
aos seus bits iniciais, que concordam exceto em cerca de uma semente em
cinquenta milhões.

### Todas as passagens até agora: `selfsim`

Cada passagem destas funções é um campo próprio — a iteração entra no hash —
com células iguais à degradação vezes as da passagem anterior. Uma fórmula que
chame uma em cada passagem lê por isso de cada vez um campo novo e mais fino, e
um modo de coloração que leia a última passagem lê só o último campo. A uma
degradação de 0.75 as células de uma vista com cinco de largura ficam abaixo de
um píxel por volta da décima sexta passagem, e daí em diante o interior é neve:
cem passagens de `randsct(672,,{0.75,0.75},,,{2,2},6)+c` em modo Julia deixam
70% dos píxeis interiores vizinhos a um dezasseis avos da paleta ou mais uns
dos outros.

O oitavo argumento devolve em vez disso **as passagens até agora, em média**,
relacionadas entre si como as oitavas de um movimento browniano fracionário: a
passagem *n* pesa d^(nH), onde *d* é a degradação — a média geométrica das suas
duas componentes tomadas sem sinal, por isso `0.75+0.75i` é 0.75 — e *H* é o
argumento. Os pesos são divididos pelo seu total, por isso a resposta fica
entre os valores médios: em [0, 1] sem inclinação.

| selfsim | |
| --- | --- |
| omitido | desligado — a única passagem, até ao último bit. O que recebe uma chamada que não nomeia nada, e um lugar deixado vazio |
| `0` | a média simples: cada passagem pesa o mesmo, e é o que dá qualquer *H* perto de zero — não desligado |
| `1` | o movimento simples: cada passagem pesa *d* vezes a anterior |
| `0.5` | mais áspero: as passagens finas guardam mais do seu peso, por isso vê-se mais detalhe |
| maior | mais suave: as primeiras passagens ficam com quase tudo |
| negativo | as passagens finas pesam mais do que as grossas, e volta a neve |
| complexo, `{1,2}` | a parte real pesa como acima, a imaginária roda as passagens — ver abaixo |

    randsct(672,,{0.75,0.75},,,{2,2},6,1)+c     o mesmo, autossemelhante

Nessa fórmula a neve desaparece: os píxeis interiores vizinhos ficam em média a
1.8 células de paleta em vez de 35, e 0.6% deles a mais de dezasseis em vez de
70%. O que resta são os triângulos em todas as escalas ao mesmo tempo, os
grandes a carregar os pequenos.

* **`selfsim` pode ser complexo.** Então d^(nH) também é: d^(n·Hr) em grandeza
  e uma rotação de n·Hi·ln d, por isso cada passagem é rodada isso mais do que
  a anterior e depois feita a média como a faz um *H* real — uma espiral de
  oitavas. A 0.75, uma parte imaginária de um dá cerca de dezasseis graus por
  passagem. Os pesos são divididos pela **soma das suas grandezas**, não pela
  sua soma complexa: isso mantém a resposta dentro do maior dos valores médios.
  A soma complexa aproximar-se-ia de zero nalguma passagem e faria explodir
  todos os valores nessa passagem para todos os píxeis ao mesmo tempo, já que
  os pesos nada têm a ver com o píxel. As rotações dão ao valor uma parte
  imaginária mesmo sem inclinação, por isso `imag` e `angle` têm algo para ler.
  Uma degradação de um tem logaritmo zero e não roda nada; `{1,0}` é `1` até ao
  último bit; e `{0,2}` é cada passagem com o mesmo peso, e rodada.
* **Zero não é desligado.** d^(n·0) vale um em cada passagem, por isso `0` é a
  média simples, e também o é qualquer *H* perto dele — de ambos os lados, e
  também complexo. Ao princípio o zero era o interruptor, o que fazia dele o
  único valor por onde a curva não passa: `0.000001` era a média simples e `0`
  só a última passagem, a um quarto do intervalo de distância. Só a última
  passagem é a ponta afastada da curva, com *H* a tender para menos infinito.
  Para não ter média, omita o argumento, ou deixe o seu lugar vazio.
* **O que se média é o valor inclinado**, por isso a inclinação e cada bit de
  `skew_mode` continuam a fazer o que faziam, e a resposta é complexa quando a
  inclinação a torna complexa.
* **Uma média aproxima os valores.** A 0.75 e com *H* igual a um, a dispersão é
  cerca de dois quintos da de uma passagem, por isso suba a velocidade de
  coloração para recuperar o contraste.
* **Uma degradação de um** pesa cada passagem por igual e dá a sua média
  simples, que assenta num meio plano à medida que as passagens se somam. É o
  que diz a definição, e é feito como se diz.
* **As passagens que uma chamada não viu são calculadas quando lhe perguntam.**
  Uma chamada num ramo que a fórmula só toma a partir da sétima passagem recebe
  na mesma as sete primeiras, e a resposta numa passagem é a mesma qualquer que
  tenha sido o caminho até ela. Um píxel novo reconhece-se pela posição além da
  passagem: só a passagem, que é aquilo por que `trap` se guia, deixaria uma
  chamada destas continuar a média que o píxel anterior deixou.
* **Nada transborda.** A média é guardada como uma só, cada passagem tomando
  dela uma parte que decorre da parte da passagem anterior, em vez de como uma
  soma e um total — por isso um *H* negativo, ou uma degradação acima de um,
  cujos pesos excederiam o que o tipo comporta, não tem onde o fazer.
* **Custa** cerca de um décimo a mais por chamada a 64 bits de mantissa e um
  sétimo a 113. Desligado, nada que se possa medir.

### Figuras em vez de ruído

`sierpinskyt`, `sierpinskyc` e `snowflake` estão onde está o ruído, e nada
nelas é aleatório: a figura está onde parece estar, e o modo Mandelbrot e o
modo Julia desenham a mesma imagem, porque uma forma no plano está no mesmo
sítio em ambos os casos.

**São fractais próprios, não campos para multiplicar dentro de um.** Escrita
sozinha — `sierpinskyt()` e nada mais — cada uma desenha a sua figura, tal como
Fractal → Mais fórmulas desenha o seu **Sierpinski**, **Sierpinski Carpet** e
**Koch Snowflake**.

As três fazem-no **levando o ponto ao seu progenitor**. Cada parte de uma
destas figuras tem um: um buraco do triângulo de Sierpinski está dentro de um
buraco maior um nível acima, um buraco do tapete dentro da célula que foi
cortada da mesma maneira, um triângulo do floco de neve sobre um lado do
hexágono do meio ou sobre a aresta livre de outro triângulo. O passo para o
progenitor é um movimento do plano — duplicar afastando-se do canto mais
próximo no triângulo de Sierpinski, ampliar uma célula pelo número de células
no tapete. A parte de topo não tem progenitor dentro da figura, por isso o seu
passo leva-a para fora do bailout, e uma parte *n* níveis abaixo precisa de *n*
passos para lá chegar: **a passagem em que um ponto sai é o nível em que
está**, e a contagem de iterações é a imagem.

O que volta na passagem em que um ponto sai importa tanto como a própria
passagem: **todos os modos de coloração exterior exceto a contagem de iterações
o leem**. Por isso a parte de topo de cada figura é atirada para fora em vez de
simplesmente declarada saída — o triângulo de Sierpinski duplica o seu buraco
do meio afastando-o de um canto e ele cai fora do triângulo, o tapete escala a
sua célula do meio em relação à célula ao lado, o floco de neve atira o seu
hexágono do meio para fora na direção para onde esse setor está virado. Cada um
leva consigo de onde veio, por isso `real`, `imag`, `angle`, a coloração suave
e os restantes têm todos algo com que trabalhar. Só a um ponto que nunca esteve
na figura é dado um número que representa «saiu».

Um floco de neve lê-se **a partir de um hexágono para fora**, e é isso que há a
saber sobre ele. Decompõe-se exatamente num hexágono regular no meio, seis
triângulos com o mesmo lado desse hexágono apoiados nos seus seis lados, doze
com um terço disso nas suas arestas livres, quarenta e oito com um nono nas
deles, e assim para sempre. Três dos seis são os cantos do triângulo de onde a
figura cresceu e três são as primeiras saliências postas nos seus lados, e
**nada os distingue**: do hexágono para fora, um floco de neve tem simetria de
ordem seis.

Lido antes a partir do primeiro triângulo, cinco oitavos da figura são o nível
um — um único triângulo plano a encher a imagem, com o floco de neve só na
franja à volta, e nada para a cor exterior dizer. Lido a partir do hexágono, os
dois primeiros níveis são cinco doze avos cada e os restantes descem de forma
regular. Um triângulo apoiado no hexágono passa para o maior triângulo que o
hexágono contém — canto sobre canto, três vezes a área — por isso é ampliado
pela raiz quadrada de três e rodado um duodécimo de volta; um apoiado na aresta
livre de outro triângulo é ampliado por três e rodado para a direção para onde
essa aresta está virada.

**Cada passo expande**, e isso é uma regra e não um acaso. O triângulo de
Sierpinski duplica, o tapete amplia uma célula pelo número de células, o floco
de neve por três ou pela raiz de três. Um triângulo apoiado no hexágono podia
em vez disso ser dobrado sobre o lado em que se apoia, o que o põe sobre um dos
seis de que o hexágono é feito e é mais arrumado; mas uma dobra é uma
isometria, e uma isometria não tem sensibilidade a onde o ponto começou.
Escrita dentro de uma fórmula maior, a figura daria então a essa fórmula uma
imagem rígida de si mesma em cinco doze avos da sua área, e a fórmula
desenharia aí uma região plana. Um teste verifica as três expansões.

O ponto viaja de facto, e é isso que as torna utilizáveis e não apenas
corretas: `z` pelo caminho é um ponto do plano como outro qualquer, por isso os
modos de coloração que o leem, a coloração suave e escrever uma figura dentro
de uma fórmula maior têm todos significado.

Um triângulo de Sierpinski e um tapete enchem a sua forma, por isso cada ponto
neles tem um nível e sai nele. Um floco de neve não enche o seu hexágono: deixa
seis cantos de **fundo**, e esse fundo não faz parte da figura e não tem nível.
É **rodado meia volta**, o que leva fundo a fundo, por isso nunca sai e aparece
na **cor interior** — o espaço vazio fica vazio em vez de tirar à figura uma
faixa da cor exterior.

Rodado, e não devolvido onde estava, porque **um ponto fixo é uma órbita parada
e não uma órbita limitada**. Parado, o ciclo corria cada passagem até ao limite
sem nada mudar; os modos de coloração interior que comparam uma passagem com a
anterior não tinham nada para comparar; e `snowflake()` escrita dentro de uma
fórmula maior devolvia a essa fórmula o ponto que ela já tinha, enquanto
`sierpinskyt()` e `sierpinskyc()` movem cada ponto que recebem. A meia volta é
preferida ao sexto de volta, que a simetria da figura também permitiria: ambos
levam fundo a fundo, mas a meia volta mantém o ponto dentro de um quadrado e de
um círculo além de dentro do hexágono. Deixa a distância ao centro inalterada,
por isso os modos de coloração interior que a leem não são afetados, e a imagem
que `snowflake()` desenha não muda um píxel.

É o único sítio em que estas figuras custam alguma coisa. A um píxel de fundo é
feita a pergunta de pertença uma vez por passagem enquanto se olhar para ele,
enquanto as outras duas figuras e o corpo do floco de neve saem na passagem que
o seu nível indica. Desenhado diretamente pelo ciclo de iteração, um floco de
neve inteiro custa cerca de 2.4 vezes o que custava quando ao fundo era dada
uma faixa própria. O seguimento da fronteira fica ligado para estas fórmulas —
as suas faixas são de facto sólidas, que é a condição para o usar — e o fundo é
uma única região sólida, por isso o que chega ao ecrã é menos do que esse
número sugere.

Escritas dentro de uma fórmula maior — `snowflake()^-2+snowflake()` e
semelhantes — as três **não** se comportam da mesma maneira, e vale a pena
dizer a razão com clareza, porque nenhum trabalho sobre a figura a vai mudar.

Um triângulo de Sierpinski é uma **poeira**: não tem área, os seus pontos nunca
saem, e cada um deles é duplicado em cada passagem. Uma fórmula construída
sobre ele tem algo caótico por baixo em toda a parte, em todas as
profundidades. Um floco de neve é **sólido**: dez doze avos dele são o hexágono
e os seis triângulos, que saem na primeira e na segunda passagem, por isso uma
fórmula construída sobre ele fica sem figura quase de imediato e desenha o que
desenhar por si. O detalhe está ao longo da fronteira de Koch, onde os níveis
vão fundo, e as regiões planas são os interiores, onde não vão.

**O que uma passagem sabe é o que essa passagem pagou.** O floco de neve decide
se um ponto é fundo descendo pela curva de Koch, e essa descida tem direito a
**um nível por passagem**. Antes corria os seus vinte e quatro níveis numa só
avaliação, por isso uma fórmula escrita à volta da figura recebia a fronteira
inteira na primeira passagem e desenhava-a por inteiro por poucas passagens que
lhe dessem — com `Iterações: 2` obtinha-se ainda cada saliência da curva de
Koch, coisa que nenhuma outra fórmula do XaoS faz e que não é para isso que
serve um número de iterações. Ficar sem níveis não é o mesmo que estar fora: um
ponto fora da moldura ou acima do invólucro é fundo e sabe-se já que o é,
enquanto um indeciso segue como parte da figura e volta a ser interrogado na
passagem seguinte com mais um nível para gastar. A imagem que `snowflake()`
desenha sozinha não muda — um ponto na figura nunca precisou da descida funda,
só da garantia de não ser fundo.

O que resta é a diferença honesta entre as duas figuras, e nenhum trabalho
sobre nenhuma delas a vai fechar. Um triângulo de Sierpinski é uma **poeira**:
sem área, os seus pontos nunca saem, cada um duplicado em cada passagem, por
isso uma fórmula construída sobre ele tem algo caótico por baixo em todas as
profundidades. Um floco de neve é **sólido**, e dez doze avos dele são o
hexágono e os seis triângulos, que saem na primeira e na segunda passagem.
Medido a um só zoom — a fração de píxeis que ficam na borda de uma faixa —
`sierpinskyt()^-2+sierpinskyt()` chega aos 17 por cento por volta das doze
passagens e `snowflake()^-2+snowflake()` aos 2.5; ambos crescem com o número de
passagens, que é o que interessa, e onde param é o que as figuras são.

`radius` significa o que `bailout` significa e lê-se da mesma maneira: **como o
quadrado da distância**. O que desenha é a forma que um bailout desse número
desenha, **inscrita nela exatamente** — defina o bailout com a forma
correspondente e o mesmo número, e a figura e o bailout ficam sobrepostos:

| | |
| --- | --- |
| `sierpinskyt` | o triângulo do bailout **triângulo −90°**, canto com canto |
| `sierpinskyc` | o quadrado do bailout **quadrado**, lado com lado |
| `snowflake` | o hexágono do bailout **hexágono 0°**, com as suas seis pontas nos seis cantos do hexágono |

O que há a saber é que um polígono de bailout põe os seus **lados** à raiz
quadrada do bailout a partir do centro — é o seu apótema, não o raio
circunscrito — por isso os seus cantos estão mais longe do que o número diz: o
dobro num triângulo, e um sétimo mais num hexágono. `sierpinskyt(4)` chega por
isso a 4 com os cantos enquanto `sierpinskyc(4)` chega a 2 com os lados, e
ambos estão certos.

Defina o bailout maior do que o raio e a figura deixa de encher a forma: há uma
**margem** entre os dois, e defina-o menor e a figura fica cortada. Ambas as
coisas são o que se pediu e nenhuma é corrigida. Mas a margem não faz parte da
figura e não tem nível, por isso não lhe é dado um: é **rodada segundo a
simetria própria da figura** — um terço de volta num triângulo de Sierpinski,
um quarto num tapete, meia volta num floco de neve — o que a mantém fora da
figura e dentro do bailout, por isso nunca sai e é desenhada na **cor
interior**, com os modos de coloração interior a trabalhar sobre ela. Antes
saía na primeira passagem e ficava com a faixa mais baixa da cor exterior, o
que desenhava a margem inteira num único tom plano e tirava uma faixa à figura.

Uma rotação e não um espelho, embora um espelho servisse igualmente para a
manter onde está: o que uma figura devolve tem de rodar com o ponto, para que
rodar a imagem rode o que nela está desenhado, e um espelho roda-o no sentido
contrário. Pela mesma razão o meio de um tapete — o bloco de células que o anel
da borda encerra, que não tem progenitor e por isso tem de ser atirado para
fora — é atirado na direção para onde está virado o seu próprio **quarto** do
quadrado, os quartos cortados pelas diagonais. Atirado sempre na mesma direção,
como era, desenhava esse meio como uma rampa de um lado ao outro sem nada a
atravessá-lo; o quadrado tem simetria de ordem quatro e assim deveria ter o que
nele se desenha. Um teste verifica que cada figura responde a uma rotação da
sua própria ordem com a mesma rotação.

**`sierpinskyt([radius=4], [kaleidoscope=1], [kaled_mirror=0])`** — o triângulo
de Sierpinski, nesse triângulo.

**`sierpinskyc([radius=4], [squares=3], [kaleidoscope=1], [kaled_mirror=0])`**
— o tapete de Sierpinski, nesse quadrado. O quadrado é cortado em `squares` por
`squares`, o anel de células ao longo da borda é mantido, tudo o que esse anel
encerra é deitado fora, e o mesmo é feito a cada célula mantida. Por isso a
imagem é **um quadrado no meio e 4×`squares`−4 à volta**: oito com três, que é
o tapete como costuma ser desenhado, doze com quatro, dezasseis com cinco. Dois
é o único número sem um anel digno desse nome, e aí o que sai é o canto
afastado, o que volta a ser um triângulo de Sierpinski — um quadrado cortado em
quatro com um canto tirado é exatamente o que é um triângulo de Sierpinski.

**`snowflake([radius=4], [kaleidoscope=1], [kaled_mirror=0])`** — o floco de
neve de Koch, com as pontas nos cantos desse hexágono, em faixas por geração do
meio para fora: o hexágono na primeira passagem, os seis triângulos sobre os
seus lados na segunda, os doze sobre as suas arestas livres na terceira, os
quarenta e oito na quarta. Os seis cantos de fundo que a figura não cobre não
têm faixa nenhuma — nunca saem, e são desenhados na cor interior.

Cada argumento tem um valor predefinido, por isso `snowflake()` é uma chamada,
e `sierpinskyc( ,5)` também — um tapete mais rendilhado com o tamanho
predefinido.

**O caleidoscópio é o da família do ruído**, escrito em último lugar em cada
uma das três e com o mesmo significado: em quantas cunhas o plano é cortado, e
qual espelho as dobra. `1` — o que recebe uma chamada que não diz nada — não
dobra nada e deixa a figura com o número que sempre foi, até ao último bit.

O que dobra é diferente, no entanto, e tem de ser. O ruído dobra a **posição**,
que fica parada durante toda a órbita; uma figura não tem uso para a posição e
lê **z**, o ponto que vai levando para baixo em direção à parte em que está,
por isso é z que é dobrado. A imagem sai da mesma maneira — a figura desenhada
numa cunha e repetida à volta da origem — e mantém-se até ao fundo, porque cada
passo é dado sobre o ponto dobrado e dobrado de novo na passagem seguinte.
Dentro da meia cunha que a dobra deixa como está, a figura é a figura.

    sierpinskyt(4,6,0)      o triângulo de Sierpinski seis vezes à volta da origem
    sierpinskyc(4,3,5,1)    o tapete em cinco cunhas, o espelho próximo
    snowflake(4,5)          cinco cunhas de floco de neve

Um floco de neve já tem a simetria de ordem seis com o espelho ao longo de cada
bissetriz, por isso dobrá-lo em seis cunhas dessa maneira é a identidade sobre
a imagem: `snowflake(4,6,0)` sai **píxel a píxel** igual a `snowflake(4)` em
90000 píxeis. Peça cinco cunhas, ou quatro, para o ver dobrar.

Uma coisa a saber sobre o tapete: o que corta, atira-o numa de quatro direções,
os quartos cortados pelas diagonais do quadrado, e sobre a própria diagonal a
escolha é um empate. Uma dobra pode levar todo um raio do plano para essa
diagonal, e aí os dois lados do empate são atirados em direções opostas — uma
linha finíssima ao longo da diagonal, que é o degrau que o tapete sempre teve e
não algo que a dobra faça. Medido em dez mil pontos com seis cunhas: cada
desacordo estava sobre uma diagonal, e os próprios pontos dobrados concordavam
até uma parte em 10^16.

    sierpinskyt()        com Fractal -> Forma da fuga -> triângulo -90
    sierpinskyc()        com forma do bailout quadrado
    sierpinskyc( ,5)     um tapete mais rendilhado
    snowflake()          com forma do bailout hexágono 0

Nenhuma custa mais do que o ruído ao lado, passagem a passagem: medidas face a
`randsc`, o triângulo de Sierpinski 0.44, o tapete 0.54 e o floco de neve 0.50
do seu custo a 64 bits de mantissa, e 0.50, 0.87 e 0.65 a 113. Ler um floco de
neve a partir do meio tornou-o mais barato além de mais bonito: seis sétimos da
figura são o hexágono ou um triângulo apoiado nele, e a esses responde-se com
quatro multiplicações e uma comparação sem percorrer a curva de todo. O que um
floco de neve custa numa imagem inteira é o parágrafo acima, e tem mais a ver
com quantas passagens o fundo leva do que com o que custa uma passagem. Não foi
grátis, e a 113 bits quase se perdeu — uma raiz quadrada e uma divisão são
ambas software aí, e qualquer delas custa mais ou menos o que custa uma figura
inteira. Por isso o recíproco da raiz quadrada do raio é guardado no local de
chamada em vez de ser calculado em cada passagem, e a própria raiz volta
multiplicando o raio por esse recíproco em vez de dividindo; o triângulo de
Sierpinski compara os seus três pesos baricêntricos escalados por um número
positivo, o que não muda nem qual é o maior nem qual é negativo; o tapete conta
em células em vez de no plano, onde a multiplicação por `squares` e a divisão
por ele se anulam; e o floco de neve toma o filho para onde passar a partir da
aresta sob cuja saliência foi encontrado, que tinha de encontrar de qualquer
modo. Percorre além disso a curva de Koch em `double`, já que a figura tem um
tamanho fixo e 52 bits de uma fração de uma aresta são mais do que os 24 níveis
que desenha conseguem usar. Um teste verifica as razões, por isso uma reescrita
que as perca di-lo.

## Paletas

Fractal → Paleta pede um número de algoritmo, e a seguir ao número diz como se
chama cada um. De 1 a 3 são os do XaoS; de 4 a 7 foram refeitos depois de medir
de que são feitos os três primeiros, e o 8 depois de medir pinturas.

| | |
| --- | --- |
| 1 | **Escuro e cor** — um quase preto com uma tinta própria e uma cor, alternados |
| 2 | **Preto, cor, branco** |
| 3 | **Quentes sobre verde-escuro** — violetas, vermelhos e laranjas vivos sobre verdes-escuros e verdes-azulados, preto ou branco de três em três |
| 4 | **Smog** — escura e melancólica: cinzentos, ocre, ardósia e violeta poeirento, com sangue de boi e escarlate pelo meio |
| 5 | **Quentes sobre a noite** — rosa, vermelho, bordô, castanho, laranja, ocre e amarelo sobre verde-azulado, petróleo, azul e índigo, no esqueleto da 3 |
| 6 | **Pares preferidos** — os pares de cores que os três primeiros põem lado a lado mais vezes, postos inteiros entre preto e branco |
| 7 | **Cores aleatórias** — cada paragem uma cor verdadeiramente ao acaso, e mais nada |
| 8 | **Kandinsky** — os degradês das suas pinturas, do mate ao vivo e de vez em quando quase ao cinzento |

### De que são feitos os três primeiros

Vinte mil paletas de cada um, geradas como o programa as gera, cada cor
designada pela mais próxima de três dúzias de cores de referência em CIELAB,
contando as cores e os pares de cores:

* **Âncoras.** O preto e o branco são entre um terço e metade das paragens, e
  as cores ficam entre eles com um período de duas ou três paragens, por isso o
  padrão vê-se nas primeiras faixas — que, numa paleta com cerca de três mil
  paragens, são tudo o que a maioria das imagens usa.
* **O 1 e o 2 contêm 256 cores e não mais.** Usam o byte baixo do gerador, e o
  byte baixo desse gerador é por si só um gerador de período 256: o vermelho
  decide o verde e o azul, e cada cor é sempre seguida pela mesma. Em vinte mil
  paletas do 2, nenhuma cor foi seguida uma única vez por outra diferente. É
  daí que vêm os seus pares recorrentes — rosa ao lado de amarelo vinte e
  quatro vezes mais do que o acaso, bordô ao lado de oliva seis, vermelho ao
  lado de petróleo oito.
* **O 3 tira de duas janelas**: cores vivas do violeta passando pelo vermelho
  até ao laranja, e escuras do verde até ao verde-azulado. Os seus pares mais
  comuns são uma cor quente viva sobre um verde-escuro.

### Os quatro seguintes

Cada um foi medido face aos três primeiros antes de ficar: a divergência de
Jensen-Shannon entre as cores que dois algoritmos mostram, por nome, e entre os
pares de vizinhos que mostram, em cinco mil paletas — 0 para a mesma paleta, 1
para nada em comum. Por essa medida, o 1 e o 2 distam 0.11, o 1 e o 3 0.16.

* **O 4, Smog**, é escuro e gasto onde o 1 é vivo sobre preto. As paragens
  fundas são uma cor mal acima do preto, as claras um cinzento encoberto com
  uma tinta e nunca branco, e uma em cada quatro um cinzento simples. As suas
  tonalidades são percorridas num anel de 256 como o 1 e o 2 percorrem as suas,
  por isso os vizinhos repetem-se como os deles, mas só nas tonalidades do smog
  — ocre e oliva, e aço, ardósia e um violeta poeirento — e, num quarto do
  anel, do carmesim ao escarlate: sangue de boi quando fundo e um escarlate
  esbatido quando claro, mantido saturado, porque um vermelho escuro que não o
  esteja lê-se como castanho. Uma paragem em cada cinco lê-se como vermelha,
  contra uma em trinta no 1. Dista de 0.43 a 0.56 dos três primeiros. Ao
  princípio era o mecanismo do 1 sobre outro anel de 256 cores, escolhido o
  mais longe possível do anel clássico, e media 0.05 do 1 — mais perto do que o
  2: a semelhança era o mecanismo, preto metade das vezes e uma cor de qualquer
  lado na outra metade, não as cores.
* **O 5, Quentes sobre a noite**, é o esqueleto do 3 com a temperatura
  invertida: cores quentes de todas as luminosidades, do rosa passando por
  vermelho, bordô, castanho e laranja até ocre e amarelo, sobre cores noturnas
  mantidas escuras, do verde-azulado passando por petróleo e azul até índigo.
  As duas janelas têm 125 graus de largura, perto dos 135 do 3; ao princípio
  eram um terço disso e liam-se como duas cores. O azul-marinho é a sua cor
  mais comum, e azul-marinho com oliva, bordô ou castanho os seus pares mais
  comuns.
* **O 6, Pares preferidos**, toma cada par de cores que um dos três primeiros
  põe uma depois da outra pelo menos três vezes mais do que o acaso — sessenta
  e seis, com as duas cores tal como saíram de facto — e põe cada par inteiro,
  as duas cores a tocar-se, entre uma paragem preta e uma branca.
* **O 7, Cores aleatórias**, é cada paragem uma cor tirada ao acaso da parte
  alta do gerador, onde o 1 e o 2 tiram da parte baixa: setenta e nove mil
  cores diferentes em vinte mil paletas onde o 2 tem 256, e nenhum par mais de
  um quinto mais comum do que o acaso. Ao princípio tinha preto e branco de
  três em três paragens, e media 0.00 do 2: o ritmo de muito escuro e muito
  claro era tudo o que se via. Sem ele, tem menos amplitude do escuro ao claro
  do que os ancorados, que é o que aleatório quer dizer.

Os quatro que substituem eram um espetro, um duotone, uma tríade e um par
complementar, raciocinados a partir do círculo cromático, e no ecrã cada um era
uma só cor. Espalhavam as suas tonalidades por toda a paleta, e no programa a
paleta tem cerca de três mil paragens — de quatro a nove no teste que devia dar
por isso, e foi por isso que não deu — por isso uma imagem via só a primeira
tonalidade. E mantinham cada paragem uma cor, funda ou pálida, sem preto nem
branco, por isso nada nelas tinha uma borda. O teste gera agora as suas paletas
como o programa faz, e pede às sete que a parte que uma imagem vê vá do escuro
ao claro e contenha mais de uma cor; ao 8, que veio depois, pede um pouco
menos, e o porquê está mais abaixo.

Uma posição guardada com os números de 4 a 7 antes desta mudança volta com
cores diferentes: o número e a semente são tudo o que uma posição regista da
sua paleta, e o que o número significa mudou. Uma guardada com os números de 1
a 3 não é afetada, até ao último bit, e uma guardada com os números de 4 a 7
nomeia um algoritmo que o XaoS original não tem e recusará.

### 8, Kandinsky

Feita a partir de dois conjuntos de imagens à sua maneira: quatro vinhetas
vivas — preto, vermelhão, amarelo, cerúleo, violeta e rosa sobre papel creme —
e uma imagem mate, laranja queimado, ocre, castanho, salva e petróleo sobre o
mesmo papel. As vinhetas dão dezoito cores, por k-means em CIELAB, e cada uma
tem uma **forma mate**: a cor da imagem mate mais próxima em luminosidade e
tonalidade, sem o croma, já que é no croma que as duas diferem. O vermelhão
passa a laranja queimado, o amarelo a ocre, o bordô a castanho, o preto a um
castanho-escuro, o cerúleo a salva, o azul a petróleo. A imagem mate não tem
cinzento, por isso o cinzento fica cinzento. O verde-água que parece ter é um
salva, 128,146,121: lê-se como verde-água só ao lado do laranja.

Cada paleta sorteia um **humor**. Quatro em cada cinco caem uniformemente entre
as formas mates e as vivas; a quinta vai para lá das formas mates em direção ao
cinzento, até quinze centésimos do seu croma, e uma paleta em cada vinte ou
assim é quase cinzenta.

**As âncoras surgem por acaso.** Uma paleta em cada dez pode transformar uma
paragem aqui e ali em preto puro ou branco puro: cada paragem depois da
primeira com o seu próprio lançamento, a uma taxa que a paleta sorteia entre
uma paragem em cada vinte e uma em cada dez, preto ou branco conforme cair o
lançamento — sem ritmo e sem alternância, ao contrário das âncoras de 1 a 3. As
outras nove paletas em cada dez são exatamente o que eram antes de as âncoras
chegarem: o lançamento já era feito para cada paragem, e ficava por usar.

O que se segue a uma cor é para onde a levam **os seus degradês**. Onde quer
que a tinta das vinhetas mude suavemente — ao longo de oito píxeis, sem um
passo maior do que um terço da mudança — foram contadas as cores nas duas
pontas, e cada paragem sorteia a seguinte a partir dessas contagens: o papel
continua numa aguada de amarelo e depois em amarelo, o preto em vermelhão (três
vezes em dez) e o vermelhão em bordô, o azul em malva e o malva em rosa, o
cerúleo num verde-água acinzentado. As arestas duras também foram contadas, e
ficaram de fora: a maioria são papel ao lado de papel de outro tom, e as
paletas feitas com elas saíam papel e cinzento. Uma paleta começa no papel ou
no preto.

Em vinte mil paletas, as do décimo mais alto de humores distam 0.06 das
vinhetas, pela mesma divergência de cima, e as do décimo acima de zero 0.07 da
imagem mate; das outras sete paletas distam de 0.39 a 0.67.

Como um degradê fica algum tempo numa parte do círculo cromático, e o 8 não tem
preto e branco a cada poucas paragens onde se apoiar, os testes pedem-lhe menos
do que aos outros. Onde a paleta é estendida longa, quarenta entradas por
paragem ou mais, uma imagem vê umas dez paragens, e uma paleta em cada cem é
então um único campo — ocre com salva, petróleo com vermelhão — que não chega a
um terço do caminho do preto ao branco nem a uma segunda cor. E três paragens,
que é muitas vezes o que uma paleta curta de 256 entradas é, são um único
degradê.

## Observar a órbita

`trap` e `stripe` guardam um número sobre toda a órbita em vez de sobre o
ponto, que é o tipo de grandeza que os modos de coloração do motor não
conseguem alcançar: aí um ciclo de cálculo é compilado para devolver uma cor, e
só o último `z` e o anterior sobrevivem para serem usados na coloração. Feito
na fórmula, não custa nada ao motor.

Ambas devolvem o valor que receberam até à última passagem que o limite de
iterações permite, e nessa o que recolheram. Por isso uma fórmula inteira é

    trap(z^2+c, 3)

e mais nada: o fractal itera como iteraria, e na última passagem o valor passa
a ser a armadilha, que os modos de coloração interior desenham depois. Um ponto
que foge nunca chega a essa passagem e mantém a sua cor exterior normal, por
isso o que estas desenham é o interior.

**`trap(orbit, shape, centre, size)`** — o mais perto que a órbita alguma vez
chegou de uma forma. `shape` é 0 por predefinição, `centre` a origem, `size` 1.

| shape | |
| --- | --- |
| `0` | o próprio centro, um ponto |
| `1` | uma reta horizontal que passa por ele |
| `2` | uma vertical |
| `3` | ambas, uma cruz |
| `4` | um anel de raio `size` |
| `5` | um quadrado de meio lado `size` |
| `6` | um losango de meia diagonal `size` |

**`stripe(orbit, density)`** — a média de `(sin(density * arg orbit) + 1) / 2`
ao longo da órbita. `density` é 4 por predefinição e é quantas riscas dão uma
volta; um número inteiro, ou não se encontram onde a volta fecha. Uma média
sobre a órbita muda suavemente com o ponto mesmo onde a contagem de iterações
salta, e é isso que desenha as fibras pelas quais o método é conhecido.

    trap(z^2+c,0)                    quão perto da origem passou a órbita
    trap(z^2+c,3)                    uma cruz, que desenha raios
    trap(z^2+c,4,{0,0},{0.5,0})      anéis dentro da cardioide
    stripe(z^2+c,6)                  seis riscas por volta

A grandeza acumulada vive no local de chamada, por isso duas armadilhas numa
fórmula guardam cada uma a sua e uma thread não pode perturbar outra. Ambas
desligam o seguimento da fronteira, pela mesma razão que o ruído: dois vizinhos
que levam o mesmo número de passagens podem ter visto órbitas bem diferentes, e
uma região preenchida sem ser calculada estaria errada.

Não são modos de coloração, e não o podem ser: os modos vivem dentro do ciclo
compilado de cada fórmula, que devolve uma cor diretamente, e há trinta e um
destes, cada um com a sua tabela de simetria por posição. O que está aqui
funciona só para uma fórmula do utilizador.

## Forma da fuga

Cálculo → Raio de fuga → Forma da fuga. A iteração para quando `z` sai de uma
região, que sempre foi um círculo. A forma vê-se nas faixas fora do conjunto,
por isso é tanto uma ferramenta de desenho como uma numérica.

| forma | foge quando |
| --- | --- |
| Círculo (clássico) | \|re\|² + \|im\|² ultrapassa o bailout |
| Quadrado | uma das componentes sozinha o faz |
| Losango | \|re\| + \|im\| o faz |
| Eixo real | a componente real sozinha |
| Eixo imaginário | a componente imaginária sozinha |
| Ambos os eixos | as duas componentes ao mesmo tempo |
| Triângulo 0° / 90° / 180° / −90° | o ponto sai de um triângulo com essa orientação |
| Hexágono 0° / 90° | o ponto sai de um hexágono com essa orientação |
| Octógono | o ponto sai de um octógono |

Os polígonos são medidos pelo seu **apótema** — a distância do centro ao meio
de um lado — por isso um polígono e o círculo do mesmo bailout tocam-se ao
longo dos lados e não nos cantos. A orientação faz parte da forma: um triângulo
rodado 120° é o mesmo triângulo, por isso 0°, 90°, 180° e −90° são quatro
diferentes, enquanto um hexágono se repete a cada 60° e por isso tem dois, e um
octógono a cada 45°, o que não lhe deixa nada para escolher. O triângulo e o
hexágono têm por isso um submenu próprio para a orientação. As suas normais são
calculadas uma vez por fotograma, por isso o ciclo interior custa de três a
oito multiplicações-somas em vez de um seno e um cosseno por lado.

Face ao círculo, o losango muda cerca de 39% dos píxeis e o eixo real 44%. A
silhueta do conjunto não se move — é decidida pelo que nunca foge — por isso
olhe para as faixas.

É escrita numa posição guardada só quando não é o círculo, por isso uma posição
que não pede uma forma continua a poder ser carregada por qualquer versão
anterior.

## Movimento browniano fracionário

Um diálogo de cada lado — **Fractal → Coloração interior → Movimento browniano
fracionário...** e o mesmo em Coloração exterior — com os três modos que
oferece e os cinco números a partir dos quais são desenhados. Era um submenu,
os três modos como entradas de escolha e os números atrás de um diálogo
Definições abaixo deles; agora é um único diálogo, assinalado sempre que um dos
três modos está em uso, e OK define o modo e os números juntos como um único
passo de anular. Uma posição continua a guardar o modo como `(incoloring N)` ou
`(outcoloring N)` e os números à parte, como `(infbmset ...)` ou
`(outfbmset ...)`, por isso um ficheiro de uma versão anterior abre sem
alterações.

O que faz é dar a uma coloração limpa um ar **usado**: as faixas mantêm a sua
forma e a sua ordem e deixam de ser perfeitas, como se a coisa tivesse ficado
ao tempo. Não é aleatoriedade píxel a píxel, que se lê como grão e foi
experimentada e posta de lado; é um campo de ruído com uma escala, por isso o
que desenha são marcas.

**No exterior** os três modos são `fbm + smooth`, `fbm + iter` e `fbm` sozinho
— a contagem suave gasta, a contagem por faixas gasta, e o ruído sozinho, que
desenha nuvens à volta do conjunto sem contagem nenhuma.

**No interior** são `fbm`, `fbm + zmag` e `fbm + decomposition`. O primeiro é o
interessante: o interior é por predefinição um único tom plano, e isto dá-lhe
uma superfície.

### Ruído de gradiente, e nenhuma grelha

O movimento é **ruído de gradiente**: cada canto de uma rede recebe uma
direção, e um ponto toma de cada canto o declive que essa direção dá à sua
distância. Antes era ruído de valor — um número em cada canto, misturado por
uma smoothstep — e isso desenhava uma grelha. Uma smoothstep não tem declive em
nenhuma das pontas, por isso o campo ficava plano ao longo de cada linha da
rede, e como a rede de cada oitava continha as linhas da primeira, todas
ficavam planas aí ao mesmo tempo. Através de uma linha o declive era zero; a
imagem saía em quadrados, com oito píxeis de lado em `FBM_ERROR.xpf`, como se
tivesse sido calculada a uma resolução mais baixa do que aquela a que era
mostrada. O ruído de gradiente é tão íngreme sobre as linhas como entre elas —
medido, 0.99 nas definições predefinidas e 1.08 nas dessa posição — e a rede de
cada oitava está deslocada em relação à anterior, por quantidades que nenhuma
duplicação torna inteiras, por isso duas oitavas nunca partilham uma linha nem
um canto. `fbm()` nas fórmulas desenha o mesmo movimento, e `randsc` é uma
oitava dele.

O que se manteve: o valor continua a ir de zero à intensidade e nunca para
fora; as **marcas têm o tamanho que tinham**, com a rede a três quintos da
frequência para que o ruído de gradiente as desenhe tão grandes como o ruído de
valor desenhava (medido por quanto o campo se move antes de se parecer consigo
mesmo só metade: 0.56 de célula nas definições predefinidas, antes e depois); o
**contraste**, do qual o ruído de gradiente sozinho tem só sete décimos — cada
oitava passa por uma curva que leva os seus valores à distribuição que o ruído
de valor tinha, monótona para que não mova nenhum contorno, e medido sobre o
movimento nas definições predefinidas oscila 0.133 a partir do meio contra
0.132 antes; e o custo, mais ou menos — os modos um pouco mais rápidos do que
antes, `fbm()` cerca de um oitavo mais lenta. O que não se manteve: uma posição
guardada com um modo fbm mantém os seus números, a sua escala e a sua força e
**desenha marcas diferentes**.

### Os cinco números

O mesmo diálogo do modo, abaixo dele. São guardados com a posição.

| | |
| --- | --- |
| **intensidade** | quantas faixas de cor o ruído desloca o valor. `4` é o predefinido; abaixo de um é um sussurro, acima de dez as faixas deixam de ser faixas |
| **frequência** | o tamanho das marcas, tão grandes quanto as desenharia uma rede com esse número de células por unidade do plano. **É esta que se sobe ao fazer zoom**: na primeira vista `8` lê-se bem, com uma amplitude de `0.02` são precisos cerca de `200` |
| **oitavas** | quantas se somam, cada uma com o dobro da frequência da anterior |
| **rugosidade** | o que cada oitava guarda da altura da anterior. `0.5` é o movimento simples; mais alto é mais granuloso |
| **semente** | a mesma imagem sempre. Mude-a para outra com o mesmo carácter |

**As oitavas e a rugosidade sustentam-se uma à outra.** Com uma rugosidade de
`0.5` a oitava oitava traz um ducentésimo quinquagésimo do total, por isso
oitavas para lá da quarta ou da quinta não mudam nada que se veja — medido,
quatro e dezoito deram a mesma imagem. Suba a rugosidade e as oitavas finas têm
algo para gastar, e aí o número começa a contar.

### O tamanho das marcas pela contagem

No exterior, mais dois campos no fundo do diálogo. **Escala do ruído** diz o
que decide o tamanho das marcas: **fixo**, onde começa, só a frequência, como
acima; **iteração**, a frequência subida com a contagem de iterações;
**iteração suavizada**, a frequência subida com a contagem suave. **Fator** diz
quanto: a frequência duplica esse número de vezes a cada dez iterações, por
isso a `1` as marcas têm metade do tamanho onde a órbita levou mais dez
iterações a sair, e um quarto onde levou mais vinte. A `0` a contagem não muda
nada. Onde a contagem sobe, para o bordo do conjunto, as marcas ficam mais
finas e seguem o pormenor.

Uma frequência que muda de um píxel para o seguinte não pode simplesmente ser
estendida sobre o plano: o campo esticar-se-ia, tanto mais quanto mais longe da
origem. Por isso o movimento é lido nas duas potências de dois de cada lado da
contagem, cada uma um campo próprio de frequência fixa, e as duas são
misturadas como uma textura entre os seus mipmaps, trazendo a mistura de volta
ao contraste de um só campo. Com **iteração suavizada** as marcas mudam de
tamanho sem costuras; com **iteração** mudam nos bordos das faixas, cada faixa
com um grão próprio. No interior a contagem é a mesma em todo o lado, por isso
o interior não tem nenhum dos dois campos.

Entre duas duplicações leem-se dois campos, e a coloração custa cerca do dobro
do que com **fixo**: 248 ns por píxel contra 127, com quatro oitavas. Onde a
contagem cai exatamente numa duplicação, como faz qualquer contagem inteira com
um fator de `10`, lê-se um só campo e o custo é o de **fixo**. Com **iteração
suavizada**, `fbm + iter` e `fbm` são desenhados pelo ciclo suave, que faz
algumas iterações para lá do bailout para calcular a contagem suave. Um fator
demasiado alto leva as marcas para baixo de um píxel, onde se tornam grão.

Com **fixo** desenha como antes, bit a bit. Uma posição guarda os dois como
`(outfbmscale 'smooth 1.5)`, e só quando foram mexidos, por isso uma posição
que não lhes toca continua a abrir em versões anteriores.

### O que custa, e a única coisa que não faz

O ruído é calculado **uma vez por píxel**, na coloração, não uma vez por
iteração. Numa imagem de 90000 pontos com uma amplitude de `0.02` e 400
iterações, `fbm + smooth` demorou um quinto mais do que `smooth` simples com
quatro oitavas e um quarto mais com dezoito; a 4000 iterações, onde a órbita
faz mais do trabalho, as mesmas quatro oitavas custaram um décimo.

As marcas estão ancoradas **ao plano**, não à vista: é sujidade na coisa, não
na lente. Por isso o zoom amplia-as, e se for suficientemente longe acaba
dentro de uma só marca e a imagem volta a estar limpa. Subir a frequência é a
resposta, e é por isso que a frequência é uma definição e não uma constante.
Não pode antes ser ligada ao zoom: o motor reutiliza píxeis calculados à escala
anterior, e uma coloração que se movesse com a vista deixá-los-ia com as marcas
erradas.

Os modos são numerados **depois** dos de cor real, por isso cada número de modo
que uma posição guardada já contém continua a significar o que significava.

## Pavimentação

**Fractal → Coloração interior → Pavimentação...** colore o interior pelos
ladrilhos de uma das quarenta e cinco pavimentações de `randsctile`, estendidos
sobre o `z` final, quarenta por unidade como os quadrados de `squares`. O seu
diálogo tem três campos.

| campo | |
| --- | --- |
| **Pavimentação** | qual pavimentação, com o número pelo qual `randsctile` a conhece e o seu nome: **3: hexágonos**, onde começa. O separador Pavimentações de Ajuda → Referência das fórmulas do utilizador tem uma imagem de cada uma |
| **Coloração** | de onde vem a cor, em toda a parte: **ângulo de z**, onde começa; **ângulo em torno do centro do ladrilho**, o do ladrilho em que cai o ponto, que faz de cada ladrilho uma pequena roda de cor; **zmag**, o valor do modo com esse nome; **nenhuma**, nada além dos ladrilhos |
| **Ladrilhos** | como cada ladrilho se destaca dos vizinhos: **cada um deslocado de uma quantidade própria**, ao longo da paleta, onde começa; **sombreados do centro para a borda**, o mesmo, e deslocado mais um oitavo das cores de uma volta do centro para a borda, como um chanfro; **com a borda marcada**, o mesmo, e a parte exterior do ladrilho desenhada na outra metade da paleta; **uma cor lisa cada um**, uma cor sua, venha de onde vier a cor |

O modo está entre os modos que se repetem e as decomposições, e as definições
de ambos chegam-lhe. **Frequência (modos periódicos)** e **Função (modos
periódicos)** chegam aos ladrilhos como chegam a todo o argumento periódico do
interior: na frequência 2 os ladrilhos têm metade do tamanho. As [definições
das decomposições](#definições-das-decomposições) do interior chegam aos dois
ângulos, que são cortados e percorridos como os das decomposições. **Setores**,
**Rotação**, **Repetições**, **Sentido** e **Velocidade** chegam a ambos;
**Ângulo de** chega só ao ângulo de `z`, porque o de em torno do centro do
ladrilho não lhe dá uso; **Anéis por unidade** faz anéis do módulo para o
ângulo de `z`, através da frequência e da função como nas decomposições, e para
o ângulo em torno do centro outros tantos anéis do centro para a borda, que
seguem a forma do ladrilho. **zmag**, **nenhuma** e um ladrilho liso não leem
nenhuma delas.

A quantidade própria de um ladrilho vem do hash que lhe dá o nome. Para um
ângulo abrange tanto da paleta como uma volta do ângulo, de modo que uma
definição da decomposição mudada deixa as cores dos ladrilhos onde estavam;
para **zmag**, **nenhuma** e um ladrilho liso abrange a paleta inteira, de modo
que **nenhuma** com **cada um deslocado de uma quantidade própria** é o mesmo
que **uma cor lisa cada um**. **ângulo em torno do centro do ladrilho** põe uma
volta inteira de cores em cada ladrilho, o que, ao tamanho que os ladrilhos
têm, são muitíssimas cores em poucos píxeis: uma **Velocidade** abaixo de um,
ou uma frequência abaixo de um, torna as rodas legíveis.

Com **Coloração** em **ângulo de z**, **Ladrilhos** em **cada um deslocado de
uma quantidade própria** e as definições das decomposições onde começam, o modo
desenha o que desenhava antes de ter escolha, bit a bit. As pavimentações que
nunca se repetem ficam planas, sem quantidade própria, para lá de quatro mil
milhões de ladrilhos da origem.

As pavimentações periódicas custam mais ou menos o mesmo que os outros modos
interiores, de 120 a 210 ns por píxel; as que nunca se repetem, que são
localizadas de novo em cada ponto, de 250 a 490. **ângulo em torno do centro do
ladrilho**, **sombreados do centro para a borda** e **com a borda marcada**
medem onde está o ponto no seu ladrilho, o que custa até dois terços mais nas
pavimentações periódicas e pouco mais nas outras. Uma posição guarda o modo
como `(incoloring 30)`, a pavimentação como `(intessellation 39)`, só quando
não é a dos hexágonos, e a coloração como
`(intessellstyle 'tile-angle 'edged)`, só quando foi mudada.

## Mais modos de coloração

Não é o cálculo que decide a cor de um píxel: a cor vem de uma fórmula no `z`
final, no parâmetro `c` e na contagem de iterações, e há uma fórmula dessas por
cada modo. Mais vinte e dois, num submenu próprio em Fractal → Coloração
interior e Fractal → Coloração exterior, acima dos submenus de cor real.

**No interior**, onde a órbita nunca fugiu e só o `z` final e `c` dizem alguma
coisa:

| modo | |
| --- | --- |
| `real/mag` | o cosseno do ângulo: um sombreado que se fecha sobre si mesmo onde `atan2` tem uma costura |
| `max(\|real\|,\|imag\|)` | a norma quadrada, por isso as faixas são quadrados onde as de `zmag` são círculos |
| `\|real\|+\|imag\|` | as mesmas faixas rodadas um oitavo de volta |
| `min(\|real\|,\|imag\|)` | perto de zero ao longo de qualquer dos eixos, por isso raios |
| `\|z-c\|` | quão longe do parâmetro a órbita assentou, o que distingue os bolbos |
| `\|z*c\|` | os dois módulos um contra o outro |
| `angle(z)-angle(c)` | a direção em que assentou, medida em relação ao parâmetro |
| `real*imag` | uma sela, quatro lóbulos à volta dos eixos |
| `sin(real)*sin(imag)` | uma grelha posta sobre onde a órbita assentou |
| `sign(imag)` | dois tons planos: uma decomposição do interior, como a decomposição binária faz ao exterior |
| `frac(mag)` | curvas de nível a cada oitavo do módulo |
| `log(mag)` | `zmag` com o contraste levado para os valores pequenos, onde uma órbita atraída passa o seu tempo |
| `cos(\|z\|)` | anéis em torno de onde a órbita assentou, oito por unidade do módulo |
| `hexagons` | três ondas a um terço de volta umas das outras, que somadas dão uma rede de hexágonos |
| `spiral` | cinco braços em torno de onde a órbita assentou, enrolados uma vez por cada fator e do módulo |
| `cos(\|z-c\|)` | anéis em torno do parâmetro, que distinguem os bolbos como faz `\|z-c\|` |

Onze dos modos interiores repetem-se — `cos(mag)`, `mag*cos(real^2)`,
`sin(real^2-imag^2)`, `squares`, `sin(real)*sin(imag)`, `frac(mag)`, os quatro
no fundo da tabela e `tessellation` — e **Frequência (modos periódicos)**, o
primeiro campo do diálogo Velocidade da coloração do interior, diz com que
frequência. Multiplica o argumento do que se repete, por isso a 2 os senos
rodam ao dobro da velocidade, `squares` desenha quadrados com metade do tamanho
e `frac(mag)` dezasseis contornos por unidade do módulo onde desenhava oito. Em
`cos(mag)` estreita também as bandas que escolhem entre o seno e o cosseno. A
um, onde começa, desenham como sempre desenharam, bit a bit, e os outros modos
nunca a leem. Uma posição guarda-a como `(incolorfrequency 2)` só quando não é
um, por isso uma posição que não lhe mexe continua a abrir nas versões
anteriores.

**Função (modos periódicos)**, por baixo dela, dobra o argumento com uma das
[funções de velocidade da coloração](#funções-de-velocidade-da-coloração) antes
de a frequência o multiplicar, de modo que as bandas dão a volta mais depressa
onde a função sobe a pique e mais devagar onde se aplana: `square` e `cube`
apertam-nas em direção à borda e afastam-nas em torno de zero, `log` e `atan`
ao contrário, e `sin`, `cos` e `triangle` dobram o argumento sobre si mesmo, de
modo que as bandas viram e voltam. Chega a cada função periódica dos onze, cujo
argumento `a` entra como `k·f(a)`, sendo `k` a frequência: os senos e cossenos,
o módulo de `frac(mag)`, a parte real e a imaginária cujas faixas fazem os
quadrados de `squares` e os mosaicos de `tessellation`, e em `cos(mag)` também
as bandas que escolhem entre o seno e o cosseno. Em `spiral` os cinco braços
continuam cinco, porque um número de voltas não inteiro deixaria uma costura
onde o ângulo fecha: a frequência e a função agem sobre o enrolamento. Uma
função periódica dentro do argumento de outra teria também o seu dobrado. Em
`ident`, onde começa, nada muda. Uma posição guarda-a como
`(incolorperiodicfun 6)`, pelo lugar da função na lista contando a partir de
zero, só quando não é ident.

**No exterior**, onde o que há é a contagem e o ponto de fuga:

| modo | |
| --- | --- |
| `iter+angle` | a direção por onde saiu, misturada na contagem: as faixas ganham uma torção |
| `iter+log(mag)` | quanto passou do bailout, uma suavização mais barata do que o logaritmo de um logaritmo |
| `iter+real*imag` | biomorfos tornados contínuos em vez de por limiar |
| `max(\|real\|,\|imag\|)` | por que lado saiu, sem contagem nenhuma |
| `iter banded` | a contagem dobrada em oito faixas, mostrando a forma dos conjuntos de nível em vez do seu número |
| `\|real\|-\|imag\|` | quão desequilibrado é o ponto de fuga |

Os números dos modos antigos não mudam, por isso uma posição guardada continua
a carregar como antes.

A velocidade e o deslocamento da coloração chegam agora a smooth e smooth log,
o que nunca acontecia: esses dois devolviam um píxel próprio em vez de passarem
pelo passo com que os outros modos terminam, por isso ambos os controlos não
faziam nada enquanto um deles estava escolhido. Com os valores com que vêm —
velocidade um, deslocamento nulo — o passo é a identidade, por isso nenhuma
imagem guardada antes se mexe.

Smooth e smooth log seguem a forma da fuga. Interpolam entre uma passagem e a
seguinte perguntando onde a órbita atravessou o bailout, e antes perguntavam-no
ao módulo fosse qual fosse a forma escolhida — por isso fugir num quadrado e
interpolar num círculo punha mais de um quinto dos pontos que fogem, e metade
com um losango, fora da faixa a que pertencem. Cada forma entrega agora a
grandeza que de facto testa e o limiar com que a compara. Para o círculo as
duas coisas são a mesma, por isso nada do que se desenhava antes se mexe.

Funcionam também em modo Newton, onde uma órbita converge em vez de fugir: a
grandeza é então quanto o último passo se moveu e o limiar o limite de
convergência, e a interpolação é a mesma pergunta feita a um número que desce
em vez de subir. Os dois fractais de Newton, que não tinham variante suave
nenhuma, passam a ter uma.

**Convergência de Newton**, no menu Cálculo, é finalmente lida. O limite estava
escrito na iteração como um milionésimo, por isso a definição era oferecida,
guardada e recarregada sem que nada a consultasse.

**Smooth e smooth log** passaram para esse submenu no menu exterior. Interpolam
entre uma iteração e a seguinte usando quanto a órbita passou do bailout, o que
precisa de uma fórmula que fuja no bailout — vinte das trinta e uma fazem-no, e
as outras onze ou convergem (as Newton) ou param num teste próprio (as
Sierpinski, Koch, o Clock). Para essas o modo não faz nada, e não é algo que a
primeira lista deva oferecer. A fórmula do utilizador, que sempre o poderia ter
tido e não tinha, tem-no agora.

O que estes modos não conseguem fazer é tudo o que precise da órbita em vez do
seu fim: uma armadilha de órbita, uma média de riscas, um expoente de Lyapunov,
a estimativa da distância. Essas precisam de um valor levado de uma iteração
para a seguinte, o que neste motor significa uma segunda cópia compilada de
cada ciclo de cálculo. `trap` e `stripe` numa fórmula do utilizador fazem-no em
vez disso.

### Definições das decomposições

**Fractal → Coloração interior → Decomposição**, e o mesmo em Coloração
exterior, define como as decomposições cortam a volta: binary decomposition e
color decomposition no exterior, `sign(imag)`, `Decomposition-like` e
`fbm + decomposition` no interior, e os ângulos de `tessellation`. Com cada
campo onde começa — zero, e um, direita e um nos três últimos — desenham como
sempre desenharam, bit a bit.

| campo | |
| --- | --- |
| **Setores** | em quantos se corta cada repetição da volta. Zero é o do modo: os dois semiplanos das binárias, a volta inteira das outras. As binárias colorem setores alternados, por isso um número ímpar passa ao par seguinte — dois setores da mesma cor tocar-se-iam. As de cor tornam-se outros tantos degraus planos |
| **Rotação** | quantos graus os setores rodam, no sentido anti-horário |
| **Ângulo de** | só no interior: o ângulo de `z`, de `z-c` ou de `z*c`. `z-c` centra os setores no parâmetro, e `z*c` roda-os pelo ângulo do parâmetro |
| **Anéis por unidade** (interior) | anéis do módulo daquilo cujo ângulo se toma, um sim um não distinguidos: as binárias tornam-se um tabuleiro de xadrez polar. Passam por `k·f(a)` como todo o argumento periódico do interior, por isso chegam-lhes a frequência e a função dos modos periódicos |
| **Iterações por faixa** (exterior) | faixas de tantas iterações através dos setores: com a binary decomposition, a grelha de raios externos e equipotenciais |
| **Repetições** | quantas vezes a volta se repete: os setores, e nas de cor a paleta, recomeçam em cada uma. Qualquer número, mas uma parte de repetição deixa uma costura onde o ângulo fecha |
| **Sentido** | direita, como sempre foram; esquerda, ao contrário; ambidestro, para a frente e para trás em cada repetição — um espelho, sem costura quando as repetições são inteiras |
| **Velocidade** | nas de cor, quanta paleta abrange uma repetição: acima de um, mais cores em cada uma, não mais repetições. Nas binárias, o afastamento entre as cores de setores alternados: em um como sempre, em zero uma só cor, em dois o dobro — a binary decomposition leva os setores ímpares da contagem para a contagem ao contrário, `sign(imag)` o seu segundo tom a partir do primeiro |

Uma posição guarda-as como `(indecomposition 8 22.5 'z-c 6 2.5 'both 1.5)` e
`(outdecomposition 6 30 2 3 'left 2)` só quando mudaram, por isso uma posição
que não lhes mexe continua a abrir na 2.5. Mudadas, uma decomposição deixa de
ser simétrica em relação ao eixo real, e o XaoS deixa de desenhar uma metade
como espelho da outra.

## Funções de velocidade da coloração

**Fractal → Coloração interior → Velocidade da coloração**, e o mesmo em
Coloração exterior, é um único diálogo para tudo o que leva o valor que um modo
de coloração entrega a um lugar na paleta. A **função**, escolhida de uma
lista, dobra o valor antes de se aplicarem a **velocidade** e o
**deslocamento** da coloração: a paleta é lida em `f(v) × speed + shift`. Numa
coluna própria, ao lado, está a [mudança de velocidade após
ITER](#mudança-de-velocidade-após-iter), que age sobre o valor antes de a
função o receber. O do interior começa pela frequência e pela função dos modos
que se repetem, que agem sobre o próprio modo; ver [Mais modos de
coloração](#mais-modos-de-coloração). OK define tudo como um único passo de
anular. Para os modos que contam — iter, smooth e os restantes — `v` é a
própria contagem de iterações, por isso uma função pode ser pensada em
iterações. Todas são ímpares, `f(−v) = −f(v)`, por isso uma grandeza que muda
de sinal mantém-no.

Eram quatro entradas do menu: a velocidade, o deslocamento, um submenu das
funções e o diálogo próprio da mudança de velocidade. Uma posição é escrita e
lida exatamente como antes, e os comandos que essas entradas executavam
continuam lá para um script.

`sin`, `cos` e `triangle` dão voltas: a velocidade da coloração, aplicada
depois delas, diz quantas entradas da paleta uma volta abrange, e deixa como
estava a frequência com que volta. Isso faz a **velocidade inicial** da mudança
de velocidade, aplicada antes: com os dois grupos desligados e o expoente
inicial a um a função recebe `s0·v`, por isso 2 junta as voltas para o dobro e
0,5 afasta-as para o dobro.

| | |
| --- | --- |
| `ident`, `log`, `atan`, `sin`, `cos`, `square`, `cube`, `square root`, `cube root` | o que sempre houve, inalterado até ao último bit |
| `log log` | `log(1 + log(1 + v))`, para os zooms profundos onde o log sozinho ainda aperta demasiado as faixas: a um milhão de iterações o log dá 13.8 faixas e esta 2.7 |
| `exp` | `e^v − 1`, o contrário do log: abre os valores pequenos. Para as grandezas que ficam pequenas — as componentes que a cor real e as decomposições leem; a uma contagem de noventa já passa de 10³⁸ |
| `triangle` | o seno com as curvas endireitadas: o mesmo período e a mesma altura, a paleta percorrida para a frente e para trás a passo constante em vez de se demorar nas viragens |

As novas vêm depois das antigas, por isso cada número que uma posição guardada
contém continua a designar o que designava.

### Mudança de velocidade após ITER

A parte de baixo do mesmo diálogo não é uma função mas uma alteração à função
que estiver escolhida: até uma iteração que indique, a função continua
exatamente como era; a partir daí, a contagem que recebe corre mais depressa ou
mais devagar, por isso as cores aceleram ou abrandam dali em diante e não
antes. Com `square root` escolhida e ITER a 40, as cores são as da raiz
quadrada até 40 iterações e a raiz quadrada de uma contagem que vai ganhando
ritmo depois. Contém um expoente e uma velocidade iniciais e depois dois
grupos, cada um separado por uma linha:

| | |
| --- | --- |
| **Expoente inicial** | o expoente desde o início até ao primeiro ITER de um grupo ligado: a contagem passa a ser a potência simples `v^p0`, a mesma estejam os grupos ligados ou não, e com ambos desligados vale para toda a contagem. Em um, onde começa, a contagem fica intacta até ao primeiro ITER |
| **Velocidade inicial** | aquilo por que se multiplica o ritmo desse troço: `s0·v^p0`. Em um, onde começa, não muda nada |
| **ITER** | onde começa o troço do grupo, nas unidades que a função recebe. Em zero ou abaixo começa no início — o que faz de `ident` uma potência com qualquer expoente |
| **Expoente** | o expoente do troço, do seu ITER até ao seguinte ou daí em diante: como o troço se curva. Parte da contagem a que o anterior chegou, `f(ITER)`, e soma `s·ITER·((v/ITER)^p − 1)`, cujo ritmo em ITER é `s·p`: acima de um as cores aceleram dali, abaixo de um abrandam |
| **Velocidade** | aquilo por que se multiplica o ritmo do troço — linear onde o expoente é uma potência: com um expoente de um o troço é a reta `f(ITER) + s·(v − ITER)`. Em um, onde começa, não muda nada |
| **Modo** | **off**, **jump** ou **smooth**: se o grupo é usado, e como o seu troço se liga ao anterior |
| **Raio** | só em smooth, e mostrado só então: quantas iterações a ligação ocupa de cada lado de ITER. Qualquer número a partir de zero, que é o jump |

Com ambos os grupos ligados e o primeiro ITER abaixo do segundo, a contagem é
`v^p0` até `i1`, continua dali ao ritmo de `p1`, e a partir de `i2` ao ritmo de
`p2` — **cada troço com o seu próprio expoente e a sua própria velocidade**,
nunca `p0·p1` nem `p0·p1·p2`. Com `p0` = 2, `p1` = 1 em 50 e `p2` = 2 em 100,
todos em jump e todas as velocidades a um, a contagem em `ident` é `v²` até 50;
`v + 2450` de 50 a 100, um expoente de um que continua a ritmo um a partir dos
2500 a que o quadrado chegou; e a partir de 100 `2550 + 100·((v/100)² − 1)`, de
novo ao quadrado a partir de um ritmo de 2 aí. Uma velocidade de 3 no primeiro
grupo faz o troço do meio `2500 + 3·(v − 50)`. O ITER mais baixo vai primeiro
seja qual for o grupo que o tem, e **no mesmo ITER o segundo grupo substitui o
primeiro**, raio, velocidade e tudo; um grupo com ITER em zero ou abaixo põe o
seu próprio expoente e a sua própria velocidade no primeiro troço.

**Um grupo desligado fica de fora**, com o seu ITER:

| primeiro grupo | segundo grupo | a contagem corre a |
| --- | --- | --- |
| ligado | ligado | `p0` até ao ITER mais baixo, depois o expoente desse grupo até ao outro ITER, depois o do outro |
| desligado | ligado | `p0` até ao ITER do segundo, depois `p2` |
| ligado | desligado | `p0` até ao ITER do primeiro, depois `p1` |
| desligado | desligado | `p0` para toda a contagem, a potência simples `v^p0` |

Ambos começam desligados, e o expoente inicial e todas as velocidades a um,
onde nada é feito ao valor, por isso cada função dá até ao último bit o que
dava antes de existir tal coisa.

As duas maneiras de ligar, conforme mudam a contagem entregue à função
(escritas para `ident`, onde a contagem é o que a paleta vê, com `f(ITER)` a
contagem ao chegar a ITER — `ITER^p0` no primeiro — `p` o expoente do grupo,
`s` a sua velocidade e `r` o seu raio):

| modo | de ITER − r a ITER + r | daí em diante |
| --- | --- | --- |
| **jump** | — | `f(ITER) + s·ITER·((v/ITER)^p − 1)`: o ritmo vira em ITER de uma vez, o que se vê como um anel onde as faixas mudam de ritmo |
| **smooth** | o ritmo vai em linha reta do que a contagem tem em `ITER − r` ao que tem em `ITER + r`, a meio caminho em ITER, com uma parábola à volta de ITER somada que devolve o que a reta cortou ao canto do jump | a contagem do jump, até ao último dígito |

**Smooth** é simétrico em relação a ITER e muda as cores dentro do raio e em
mais lado nenhum: antes de `ITER − r` e depois de `ITER + r` a imagem é a do
jump. Onde ambos os ritmos são constantes a parábola é nula e o canto é cortado
como uma parábola o corta — com a contagem intacta até um ITER de 50 e um
expoente de zero depois, um raio de 10 dá 47.5 em ITER, um quarto do raio a
menos, e de novo 50 a partir de 60; de um ritmo de um para um de três dá dois
em ITER. Onde um ritmo se curva, a parábola devolve o que a curva acrescenta:
depois de `v²` até 50 e de um expoente de um, um raio de 20 leva o ritmo de 60
em 30 a 1 em 70, 45.5 em ITER, e continua a partir de 70 como `v + 2450`, como
faz o jump.

Antes misturava os ritmos dos dois troços tal como cada um teria continuado: o
ritmo do quadrado, ainda a crescer, continuava a crescer bem dentro do
arredondamento, por isso mal descia antes de ITER e estava a meio da descida só
meio raio depois, e o que ganhava movia todas as cores para lá do
arredondamento — 133 faixas com um raio de 20.

A parábola nunca leva o ritmo abaixo do mais lento dos dois extremos. Onde o
faria — um ritmo lento tornado muito rápido num raio curto, como uma centésima
potência em dez iterações — é mantida aí, e a contagem para lá do
arredondamento é deslocada pelo que isso retém. Um raio mais largo do que o
primeiro ITER para antes de zero, e dois que se sobreporiam são reduzidos por
igual para que cada um fique entre os seus vizinhos — raios de 40 e 40 entre
ITER a 30 de distância passam a 15 cada.

Um expoente de 100 é legal e inútil para lá de ITER: poucas iterações depois,
as faixas são mais finas do que os píxeis e a cor é ruído. Com o segundo grupo
em ITER 50 e expoente 2 e o primeiro em 100 e 1.5, ambos em jump, e um expoente
inicial de um, a contagem em `ident` é 30 em 30, `50·1.5²` = 112.5 em 75, e
`200 + 100·(1.5^1.5 − 1)` = 284 em 150: o troço a partir de 100 parte dos 200 a
que o anterior chegou, a um ritmo de 1.5.

Antes havia quatro modos, e as posições guardadas antes do raio nomeiam-nos — o
`outcolorfunafter` da 1.8.0, e `outcolorfunchange` com dois grupos e sem raio.
São lidos nos dois que existem, com todas as velocidades a um:

| guardado como | lido como |
| --- | --- |
| off, jump | o mesmo |
| smooth, que arredondava o canto no logaritmo da contagem, um pouco de cada lado de ITER por mais longe que estivesse | smooth com um raio de três quartos de ITER: a menos de quatro centésimos de ITER do que desenhava até um expoente de um e meio, um décimo até dois |
| continuous, que continuava a partir de ITER ao ritmo com que a contagem chegava | smooth a três quartos de ITER, com o expoente `0.66 + p/4 + 0.09·p²` — um em um, onde continuous também não mudava nada — a menos de doze centésimos de ITER dele até um meio |

Encontrado pondo as curvas lado a lado em quatro ITER: uma imagem destas sai
parecida, não igual ao píxel. Uma posição cuja mudança de velocidade não passa
de um grupo em jump, com o expoente inicial e todas as velocidades a um e o
resto como começa, é guardada como a 1.8.0 a guardava,
`(outcolorfunafter 40 1.5 'jump)`, e abre na 1.8.0 e aí desenha igual. Tudo o
resto é guardado como a 2.0 o guardava, na ordem do diálogo sem as velocidades,
`(outcolorspeedchange 1 40 1.5 'smooth 10 50 2 'off 10)`; e as velocidades, só
onde alguma não é um, à parte como `(outcolorspeedlinear 1 3 1)` — a velocidade
inicial e as dos dois grupos — para que uma posição que não as use continue a
abrir na 2.0. A mudança de velocidade sozinha, como o seu próprio diálogo a
tinha, é o comando `outcolorspeedgroups`, com os seus doze campos pela sua
ordem.

### Para lá do fim de um int

A paleta é alcançada em 256 avos de entrada através de um int, e um valor para
lá do que um int comporta — 8.4 milhões de faixas — era convertido na mesma, o
que é indefinido: cada píxel desses ficava com uma e a mesma entrada. `cube`
chega lá às 204 iterações, por isso uma imagem colorida com o cubo para lá
disso era uma mancha plana; `exp` e uma grande mudança de velocidade chegam lá
de imediato. Agora cai onde o seu valor diz, dando a volta à paleta como
qualquer valor menor; um grande demais para que os seus 256 avos caibam num
double cai na primeira entrada, como um infinito sempre fez. Dentro de um int
nada se mexe. Fora dele, a única mudança a uma imagem guardada antes está na
coloração suave de Magnet e Magnet2 longe do conjunto, que era uma única cor
plana e agora mostra as faixas que a contagem aí descreve.

Vieram com isso duas reparações mais pequenas. Escolher uma função de
velocidade mudava o visto no menu da coloração interior em vez de no seu. E
**initstate**, com que começa toda a posição, não repunha os cinco números do
fbm nos seus valores predefinidos — uma posição só os guarda onde diferem
deles, por isso uma carregada depois de outra que os tinha mexido desenhava com
os da outra; agora repõe-nos, e com eles desliga a mudança de velocidade após
ITER.

### A contagem, inteira

O valor que um modo de coloração entrega perdia dígitos a caminho da paleta em
dois sítios. A coloração suave fazia da sua contagem um int de 256 avos de
iteração antes de a velocidade da coloração, a função ou a mudança de
velocidade a receberem; e o que a função devolvia era um float, que para lá de
65 536 entradas guarda menos de um 256 avo de uma e para lá de dezasseis
milhões nem sequer uma inteira. Nenhum se vê a velocidade um, onde um 256 avo
de iteração é um 256 avo de entrada. Ambos se veem assim que algo faz as cores
correr: com uma mudança de velocidade que leva a contagem de zero a dez mil
entradas por iteração, uma vista a trinta iterações, onde o ritmo era 2500,
saía em riscas de uma só cor com oito píxeis de largura, cada uma a um 256 avo
de iteração e a dez entradas da seguinte, onde a paleta devia ter percorrido
todas as entradas intermédias.

Agora a contagem continua inteira, e a função, a velocidade da coloração e a
mudança de velocidade são calculadas no tipo próprio da compilação, como as
fórmulas — long double na compilação normal, `__float128` na quad — e devolvem
o que calculam. As duas compilações colorem tão parecido como antes: nas
imagens experimentadas diferem onde sempre diferiram, num zoom profundo cujas
órbitas a quad calcula mais longe, e em mais lado nenhum.

Duas coisas ficam como estavam.

* **Um infinito** — o centro de um potencial que foi para o infinito, o fundo à
  volta de uma fórmula que foge de imediato — cai na primeira entrada, como um
  valor que não é número. Feito em 256 avos num float, cada valor desses saía
  infinito e caía aí; que cor deve tomar um infinito ninguém sabe, e assim toma
  a que sempre tomou. Um valor finito para lá do que um float comporta cai onde
  diz.
* **Um valor na borda de um 256 avo** fica nela: o valor é arredondado a um 256
  avo de um 256 avo antes de ser cortado a um 256 avo, já que a coloração não
  diz nada mais fino. A parte real de uma órbita de Newton para numa raiz, em
  um ou em menos meio, o que põe `iter+real` exatamente numa borda dessas, e os
  últimos dígitos da convergência poriam cada píxel de um lado ou do outro — um
  256 avo de entrada, que ninguém vê, até que um filtro de bordas desenha cada
  píxel desses como borda e as curvas de nível saem dentadas como engrenagens.
  O float arredondava a maior parte e deixava entalhes; agora não há nenhum.

Renderizados face à compilação anterior a 640×480 — antes da mudança da secção
seguinte — dos 169 exemplos 83 desenhavam ao píxel como antes e mais 71
diferiam num só nível em no máximo 113 píxeis. Dos outros quinze, os de cor
real e smooth log — `perb4`, `spiral4`, `truecol1`, `truecol3` — diferem até
cinco níveis em até um sétimo dos píxeis, por as suas componentes também terem
voltado como floats; `edge` perdeu os entalhes; e o resto difere numa mão-cheia
de píxeis onde um valor estava na borda de uma entrada e um filtro ou uma
paleta desenhada sem mistura mostra o degrau. Magnet e Magnet2 em coloração
suave são os que mais mudam: longe do conjunto a sua contagem é negativa, o que
o int, lido sem sinal, avançava quatro mil milhões de 256 avos, e o float
desenhava depois em degraus de uma entrada ou duas — manchas planas com arestas
duras onde agora está o degradê que a contagem descreve.

Numa imagem que sai igual custa de dois a cinco por cento na compilação normal,
e de sete a vinte na quad, onde os logaritmos da contagem e as funções são
calculados em software. Onde as riscas se tornam um degradê a renderização
também demora mais: o preenchimento por estimativa só preenche um píxel a
partir de vizinhos de uma única cor, e a imagem escrita é maior — a vista a
trinta iterações, a 1200×1200, demorou 1.07 segundos onde as suas riscas
demoravam 0.30.

### O ritmo da contagem suave

A coloração suave soma à contagem de iterações inteiras quanto da última o
ponto já tinha percorrido quando fugiu. Calculava-o onde o bailout cai entre a
grandeza antes da passagem de fuga e depois, linearmente no seu logaritmo.
Assim a contagem continua sem quebras de uma iteração para a seguinte, mas o
seu ritmo não: de uma passagem para a seguinte a grandeza é quase elevada ao
quadrado, por isso o ritmo muda ao longo de cada iteração e vira de uma vez
onde começa a seguinte. Numa vista com 2.6×10⁻⁷ de largura, onde a contagem vai
de 59.94 a 60.02, virava 4.26 na borda entre as duas. A uma velocidade de
coloração de um isso move uma cor alguns centésimos de entrada; com uma mudança
de velocidade a correr a dez mil entradas por iteração, são as faixas a mudar
de densidade de uma vez ao longo de cada borda de uma iteração.

Agora a órbita é levada para lá do bailout até o seu módulo ao quadrado passar
de um milhão — mais três ou quatro passagens para um quadrado — e a fração é
tomada no logaritmo do logaritmo, onde uma elevação ao quadrado é um passo de
comprimento um:

    n + ln(ln T / ln s) / ln(ln m / ln s)

com `s` e `m` o módulo ao quadrado de cada lado de um milhão, `T` o bailout e
`n` as passagens feitas. As passagens para lá do bailout descontam-se a si
mesmas, e o grau é o da própria órbita, `ln m / ln s`, por isso um cubo é
contado como um cubo. A contagem fica a menos de um décimo de iteração do que
era, e o seu ritmo continua: a viragem de 4.26 nessa vista sai a 1.0000, e ao
longo de um raio vindo de fora do conjunto as viragens de 1.3 a 1.8 em cada
borda saem a 1.0025.

| | a contagem |
| --- | --- |
| as Mandelbrot, Triceratops, Mandelbar, Lambda, Manowar, Spider, Phoenix e Octo, no círculo | levada para lá do bailout |
| uma fórmula do utilizador que foge no círculo | levada para lá enquanto a órbita cresce pelo menos como uma potência de um e meio, e contada no bailout onde não |
| as três Barnsley | no bailout: multiplicam a órbita por `c`, de modo que o seu logaritmo cresce um passo em vez de duplicar, e a contagem no bailout já é a certa para elas |
| Newton, Magnet, Catseye | no bailout: convergem, ou testam a sua fuga à sua maneira |
| todas as formas da fuga exceto o círculo | no bailout |

As formas ficam como estavam porque cada uma mede a órbita com algo que roda
com ela, enquanto cada passagem lhe duplica o ângulo: levado adiante com a sua
própria medida, o grau saía instável o bastante para piorar as viragens — as do
hexágono de 1.28 para entre 1.67 e 3.57 ao longo desse raio — e as do quadrado
chegavam a voltar para trás.

Dos 169 exemplos, os 154 que não colorem de forma suave desenham ao píxel como
antes disto. Os quinze que o fazem mudam até treze níveis em até nove décimos
dos píxeis — `eye` até 39, onde as entradas vizinhas da sua paleta estão muito
afastadas, e `spiral4`, cujo smooth log numa paleta de cor real amplia qualquer
mudança na contagem, no meio de um pequeno pormenor. Face à 2.1, com as
mudanças da secção anterior, 84 desenham ao píxel como antes e mais 59 diferem
no máximo três níveis em no máximo 500 píxeis.

Custa essas passagens e mais dois logaritmos por cada píxel que foge: nas
imagens suaves experimentadas, de 11 a 28 por cento na compilação normal e de
18 a 20 na quad, onde um logaritmo em `__float128` demora 946 ns contra 39 em
long double; uma imagem que não colore de forma suave fica inalterada.
Calculados em double, os quatro logaritmos teriam custado de 8 a 19 por cento
na compilação normal e poupado 8 na quad — um double leva a fração a dezasseis
dígitos, onde uma velocidade de dez mil precisa de nove — mas a coloração é
calculada toda no tipo próprio da compilação, como as fórmulas.

## Zoom

**Zoom por seleção** (Vista → Zoom) substitui o zoom contínuo por um retângulo:
arraste com o botão esquerdo e solte para ir para lá, como faz o Ultra Fractal.
Só para aproximar. O retângulo é forçado à proporção da janela alongando o seu
lado curto, nunca cortando o longo, e funciona corretamente com o fractal
rodado. O modo Julia rápido fica com o arrastamento enquanto está ativo, e
devolve-o à saída.

**Passos fixos**, no mesmo submenu:

| | |
| --- | --- |
| Aproximar 2× | `Ctrl` `+` |
| Afastar 2× | `Ctrl` `-` |
| Aproximar 10× | `Ctrl` `Alt` `+` |
| Afastar 10× | `Ctrl` `Alt` `-` |

Escalam em relação ao centro da vista, por isso desfazem-se exatamente uns aos
outros.

## Rodar e centrar a vista

**Fractal → Vista** define o centro, o raio e o ângulo por números. Por baixo
do ângulo tem duas ferramentas para esses mesmos campos. Preenchem os campos e
não fazem mais nada: OK aplica o que mostram como um único passo de anular,
Cancelar nada disso, e uma posição guardada nunca sabe delas.

**Fatias** corta uma volta completa em tantas iguais, e **−** e **+** rodam o
ângulo para trás ou para a frente uma: 60 graus para seis, 360 para uma. Com
nenhuma, os dois ficam desativados. As voltas contam-se a partir do ângulo de
onde partiram, por isso seis sextos são uma volta inteira até ao último dígito;
escrever um ângulo, ou mudar as fatias, recomeça a contagem a partir daí. O
número de fatias é guardado enquanto o XaoS está aberto, para a próxima vez que
o diálogo for aberto.

**Centrar no eixo** move o centro a direito através da imagem (**R**) ou a
direito para cima (**I**) até a origem ficar no meio nesse sentido, ao ângulo
do diálogo. Sem rotação, é a parte real ou a imaginária que passa a zero e a
outra que se mantém. Com rotação, é o mesmo ao longo dos eixos próprios da
imagem: R a 90 graus anula a parte imaginária, e R e depois I levam a origem ao
meio a qualquer ângulo. Um centro que já lá está é deixado como está, e uma
parte que a rotação deixa como arredondamento é tomada por zero, por isso o
campo diz 0 e não 1e-17.

A própria rotação é calculada no tipo próprio da compilação, pi incluído. Era
calculada a partir de um pi double, o que na compilação quad significava também
um seno e um cosseno double: um par que é uma rotação só até dezasseis dígitos,
por isso o centro de uma vista rodada era desenhado a 5e-17 da sua distância à
origem longe de onde estava. Abaixo de 1e-16, que é para onde a quad serve, uma
vista rodada 30 graus mostrava um sítio completamente diferente. Um quarto de
volta é agora exato. Na compilação normal uma imagem rodada é a mesma exceto
píxeis no limite da sua precisão: 52 dos 64 exemplos com um ângulo ficam
inalterados, os restantes movem píxeis dispersos onde o pormenor é mais fino do
que a aritmética.

## A mesma imagem com qualquer número de threads

Uma imagem saía de vez em quando com alguns píxeis diferentes quando era
desenhada com várias threads: `examples/Hubicka/julia3.xpf` renderizado com
quatro threads numa máquina ocupada diferia em 26 de 3000 renderizações, até 37
píxeis ao longo dos seus filamentos, e nunca com uma thread.

Os píxeis vinham do preenchimento por estimativa. Preenche um píxel a partir
das linhas e colunas já calculadas à sua volta quando concordam, e em modo
Mandelbrot usa os mesmos vizinhos para decidir se vale a pena verificar a
periodicidade de um píxel. A estimativa avança por passagens — primeiro uma
linha em cada seis, depois as linhas intermédias — e as passagens corriam em
todas as threads ao mesmo tempo sem nada entre elas: uma thread passava à
passagem seguinte enquanto as outras ainda estavam na anterior, e uma linha era
marcada como feita antes de ser calculada. Por isso o facto de um vizinho
contar como calculado dependia de quão longe as outras threads tinham chegado.

Agora cada passagem tem duas fases, primeiro as suas linhas e depois as suas
colunas, e cada thread acaba uma antes de alguma começar a seguinte; o que uma
fase calcula não é usado como vizinho dentro dela. Os vizinhos são então os
mesmos seja qual for o número de threads, e as 3000 renderizações saem todas
iguais. Com uma thread a ordem é a de sempre.

Os filtros também saem iguais. **3D** desenha uma faixa de colunas por thread,
e o sombreado de uma coluna é levado a partir da coluna duas à direita na linha
de baixo, por isso as duas últimas colunas de uma faixa leem as duas primeiras
da faixa à sua direita; antes cada faixa lia aí cópias da sua própria última
coluna, e a imagem mudava ao longo das bordas das faixas com o número de
threads. Agora uma faixa espera, linha a linha, pelo que a faixa à sua direita
desenhou aí, e as faixas continuam a correr lado a lado. **Campo de estrelas**
tirava as suas estrelas de um único gerador aleatório que todas as threads
faziam avançar ao mesmo tempo; agora cada linha começa no ponto da sequência a
que teria chegado uma só thread a desenhar a imagem inteira, por isso as
estrelas são as mesmas com qualquer número de threads. Continuam a mudar de um
arranque do XaoS para o seguinte, como sempre: o gerador tira a semente do
relógio.

**O antialiasing** nunca foi a causa, embora três exemplos desenhados com ele —
`circle`, `helloween` e `warriormask` — saíssem diferentes de uma renderização
para a seguinte. O que diferia era a mensagem com o nome da fórmula:
desaparecia ao fim de cinco segundos do relógio de parede, e uma renderização
abrandada por uma máquina ocupada desenhava-a em menos fotogramas, ou em
nenhum. Numa renderização uma mensagem dura agora cinco segundos do tempo
próprio da animação, seja qual for a máquina; na janela nada muda.

Dos 169 exemplos renderizados com uma thread, 167 desenham ao píxel como a
1.8.0 os desenhava: `galaxy` tira as suas estrelas do relógio, e `heart` mantém
agora uma mensagem que a 1.8.0, cronometrando-a pelo relógio de parede, tinha
retirado antes de a imagem ser desenhada. Renderizados oito de cada vez em
quatro núcleos, uma vez com uma thread e três vezes com quatro, todos desenham
igual exceto `galaxy`.

## Proporção

**Vista → Proporção** redimensiona a janela uma vez, para que a imagem tome a
forma escolhida: 1:1, 4:3, 3:2, 16:10, 16:9, 2:1, 21:9, ou 9:16 na vertical. A
janela não fica presa a ela: pode depois ser arrastada para qualquer forma, e o
visto no menu está na forma que a imagem tiver no momento, ou em nenhuma.

A imagem mantém o número de píxeis que tinha e muda só a forma — ou, quando não
caberia no ecrã, passa a ser a maior dessa forma que cabe. A partir de ecrã
inteiro ou de uma janela maximizada volta primeiro a uma janela normal, ainda
com os píxeis que tinha, por isso 16:9 aí é a janela 16:9 mais larga que o ecrã
comporta. A forma é cumprida ao píxel, também num ecrã com escala: a 125% uma
contagem de píxeis em cada cinco não é possível, e o tamanho é escolhido entre
as que são.

Aquilo para que serve é a renderização. Uma posição é desenhada de forma que
toda a vista que contém fique na imagem, e uma imagem de outra forma mostra
mais plano do seu lado mais comprido. Por isso uma vista enquadrada na janela
normal 4:3 e renderizada a 1920×1080 tem um terço a mais dos lados do que
alguma vez esteve no ecrã; enquadrada numa janela 16:9, a renderização é o que
se viu.

## Antialiasing em luz linear

Filtros → Antialiasing em luz linear, ou `-linearaa`. O antialiasing normal faz
a média das cores das amostras tal como estão guardadas, o que escurece as
bordas porque os valores guardados não são proporcionais à luz. Este converte
para luz linear, faz a média aí, e converte de volta, por isso uma borda mantém
o brilho daquilo entre o que está.

## Idiomas

Vista → Idioma muda o XaoS para outro idioma, depois de reiniciar. Seis estão
completos: alemão, espanhol, francês, italiano, português e russo — os menus,
os diálogos e as mensagens, a referência das fórmulas do utilizador, e os
diálogos OK, Cancelar, de ficheiros e de cores que vêm do próprio Qt. Estes
últimos vêm agora incluídos, quando antes eram procurados onde o Qt está
instalado e no idioma do sistema fosse qual fosse o escolhido. Os outros
idiomas que o XaoS herdou vão até onde os seus tradutores os levaram, e para lá
disso estão em inglês.

Na referência das fórmulas do utilizador as descrições estão traduzidas e os
nomes e as listas de argumentos não: são a sintaxe em que se escreve uma
fórmula.

Este guia está traduzido nos mesmos seis idiomas, ao lado do guia inglês
`ms-xaos-guide.md` como `ms-xaos-guide.de.md`, `.es.md`, `.fr.md`, `.it.md`,
`.pt.md` e `.ru.md`, e Ajuda abre o do idioma em que o XaoS está. Cada um diz
de que versão do guia inglês foi feito; onde diferem, a referência é o inglês.

## Compilação

    cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build

Os dois executáveis vão parar a `bin/`. Opções:

| | |
| --- | --- |
| `-DQUAD=OFF` | salta o segundo executável |
| `-DDEEPZOOM=ON` | torna quad o executável principal, e não compila um segundo |
| `-DDIAG=ON` | compila também `XaoS-diag`, que regista em `xaos-diag.txt` a ordem de sobreposição e de cálculo |

`cmake --build build --target deploy` enche `bin/` com o runtime do Qt e do
compilador para que corra numa máquina sem Qt.

## Testes

    ctest --test-dir build

351 testes: o analisador de fórmulas nas duas precisões, a exatidão de `gamma`,
`lambertw`, `erf` e `randsc` face a referências exatas onde existem, os ciclos
de iteração face a somas de verificação registadas, as funções de velocidade da
coloração face ao que davam antes e àquilo para que servem, as posições
guardadas que vão e voltam e cada renderização que termina de forma limpa, a
gravação e reposição da sobreposição, o histórico de anular, a rotação e a
centragem do diálogo Vista, a referência da ajuda face à própria tabela do
analisador, e o `poly` aleatório do seu separador Utilitários face ao polinómio
escrito termo a termo.

Um deles mede velocidade em vez de respostas. Cronometrar em nanossegundos
afirmaria em que máquina a bateria está a correr, por isso afirma antes uma
razão: quanto custa uma degradação às funções de ruído, contado em
multiplicações. Deveria ser dois, um por componente, porque o tamanho é
multiplicado pela degradação uma vez por passagem; elevá-la antes à potência da
passagem dá seis ou sete e cresce com o limite de iterações. O teste falha nas
duas implementações que esta substituiu.
