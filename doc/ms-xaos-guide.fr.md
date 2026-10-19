# MS XaoS

> Traduction du [guide en anglais](ms-xaos-guide.md) de la version 2.8. En cas
> de différence, l'original anglais fait foi.

Un fork de [XaoS](https://github.com/xaos-project/XaoS) 4.3.3. Tout ce que fait
l'original, il le fait encore ; ce guide décrit ce qui a été ajouté ou modifié,
et pourquoi. Version 2.4.

## Deux exécutables

`XaoS.exe` calcule avec 64 bits de mantisse, `XaoS-quad.exe` avec 113. Les deux
sont compilés à partir des mêmes sources — la précision est un type choisi à la
compilation, pas un réglage — et celui en quad s'appelle **MS XaoS Quad**, pour
que les deux se distinguent une fois lancés.

Quad zoome environ quinze ordres de grandeur plus profond. En rendant la même
position à 1e-18, l'exécutable ordinaire garde 21 % de son détail, le quad 92
%. Ce n'est pas le choix par défaut parce que `__float128` est émulé de façon
logicielle : environ huit fois plus lent sur les formules intégrées et quinze
sur les formules utilisateur.

Une position enregistrée par l'exécutable quad contient `(precision 113)`.
L'ouvrir avec une version moins précise le signale et la dessine quand même ;
l'ouvrir avec une version qui ne connaît pas la commande est refusé, et c'est
voulu — une image calculée sur 113 bits ne peut pas être reproduite sur 64, et
ne pas l'ouvrir est plus honnête que de dessiner autre chose.

## Les nombres dans les dialogues

Un champ qui demande un nombre en accepte un écrit de n'importe quelle façon
lisible, 1e-18 aussi bien que 0.5. Les dialogues à un seul champ utilisaient
une boîte de sélection numérique, qui garde deux décimales sauf indication
contraire, si bien que 0.01 survivait et que 0.001 devenait nul — une
convergence de Newton d'un millionième ne pouvait pas du tout être saisie.

Ce qui est affiché est le texte le plus court qui se relit comme le nombre
exact conservé, donc un millionième de millionième apparaît comme 1e-18 et non
comme 9.99999999999999999978e-19, ce que fait l'impression de tous les chiffres
de la version d'un décimal qu'un flottant binaire ne peut pas représenter
exactement. Une coordonnée qui a besoin des vingt et un chiffres les reçoit
toujours : en avoir besoin revient à ne pas se relire à l'identique sans eux.

La forme exponentielle est utilisée à partir de mille milliards et en dessous
de cent milliardièmes, et la forme décimale simple entre les deux — ainsi vingt
est vingt, cent milliards s'écrit en entier, et 1e-8 non. Livrée à elle-même,
l'impression choisirait en comparant l'exposant au nombre de chiffres demandés,
ce qui ne sert à rien quand les chiffres sont choisis pour leur brièveté : un
seul chiffre suffit pour se relire comme vingt, et vingt sortait comme 2e+01.

## Positions enregistrées

La vue est écrite avec toute la précision de la version — 21 chiffres
significatifs d'ordinaire, 36 en quad. Elle était écrite avec un nombre de
chiffres déduit de la profondeur du zoom, qui sous-estimait et plafonnait à 17,
si bien qu'une position profonde se rouvrait à une fraction de pixel de
l'endroit où elle avait été enregistrée. Avec beaucoup d'itérations, c'est une
image visiblement différente.

Tout ce qu'une position peut contenir survit désormais sans changement à un
enregistrement et un chargement : les 169 exemples fournis et 141 positions
générées couvrant chaque formule, plan, mode de coloration et filtre font
l'aller-retour au pixel près.

## Annuler

Edition → **Annuler** (`U`, ou `Ctrl` `Z`) revient d'un pas en arrière et
**Rétablir** (`Ctrl` `Y` sous Windows) avance de nouveau. Un pas est tout ce
qui change l'image : une commande — une formule, une palette, un mode de
coloration, un dialogue, un zoom fixe —, un appui de la souris avec le zoom, le
déplacement ou la rotation qu'il lance, un zoom par sélection, et les
itérations réglées avec les flèches. Une commande qui ne change rien, comme
Recalculer, n'est pas un pas, et annuler ne remonte pas au-delà de l'image avec
laquelle XaoS s'est ouvert, après la configuration et la ligne de commande.
Chaque pas est conservé comme une position enregistrée, donc annuler rétablit
tout ce qu'une position contient, aussi exactement que la rouvrir ; les 256
derniers sont conservés.

Cela allait mal de deux façons. La souris conservait son pas au relâchement,
pendant que le zoom ralentissait encore, donc annuler ramenait une image à
mi-chemin ; le pas est désormais conservé là où la vue s'arrête. Et chaque
commande conservait un pas, qu'elle change quelque chose ou non, tandis
qu'annuler reculait de deux en supposant avoir conservé le dernier lui-même —
après un Recalculer, il ramenait l'image déjà à l'écran. Désormais une image
n'est conservée que là où elle diffère de la dernière conservée, et annuler et
rétablir avancent d'un pas parmi elles. Un nouveau pas après une annulation
coupe ce qui aurait pu être rétabli, comme dans un éditeur.

Ce qui a bougé de lui-même depuis le dernier pas — le pilote automatique, une
rotation continue — est conservé comme un pas à part quand on appuie sur
Annuler, pour que Rétablir y revienne.

Annuler a trouvé une erreur dans les positions. En charger une avec une palette
intégrée après une palette personnelle laissait la palette personnelle marquée
comme utilisée, si bien que la position était réenregistrée comme palette
personnelle : pas ce qui avait été chargé, et pour annuler pas le pas auquel il
était revenu, ce qui coûtait le rétablissement. Désormais une palette intégrée
efface la marque, comme l'a toujours fait une palette générée.

## Formules utilisateur

**Référence.** Aide → Référence des formules utilisateur liste chaque fonction,
variable et notation qu'accepte le langage des formules, dans une fenêtre qu'on
peut laisser ouverte à côté de la formule en cours d'écriture. Elle est
regroupée selon l'usage des fonctions, avec les fonctions spéciales — `erf`,
`gamma`, `lambertw` — dans une section à part, puisqu'aucune n'est élémentaire
et qu'on n'y recourt pas par hasard. Un quatrième onglet liste les nombres qui
apparaissent comme arguments et ont un sens particulier — par rapport à quelle
forme mesure un piège, avec quel miroir le kaléidoscope replie — au lieu de les
enterrer dans la description de la fonction qui les prend. Un cinquième montre
une image de chaque pavage que prend `randsctile`, puisque quarante-cinq
nombres avec une ligne de mots chacun en disent peu sur les formes.

Un sixième, **Outils**, écrit des formules au lieu de les décrire. Pour
l'instant il contient un `poly` aléatoire, à partir de quatre champs :

| | |
| --- | --- |
| longueur minimale et maximale | m, le nombre de coefficients, tiré entre les deux ; le degré est m − 1. Les deux sont nécessaires, en nombres entiers de 1 à 254 — aussi loin que va `poly`, l'analyseur prenant 255 arguments par appel et `z` étant l'un d'eux |
| module minimal et maximal | la taille de chaque coefficient, tirée uniformément entre les deux. 0 et 1 si laissés vides ; au plus 1000000 |

Chaque coefficient est complexe, `{re,im}`, avec son argument tiré au hasard
sur tout le tour, et chaque partie est écrite avec quatre chiffres
significatifs de son module, pour qu'un intervalle de 0.001 à 0.002 ne
disparaisse pas à l'arrondi. Il n'y a pas de `c` dans la formule. Une formule
utilisateur fait partir `z` du pixel en mode mandelbrot comme en mode julia,
donc ce qu'elle trace est l'ensemble de Julia rempli du polynôme lui-même,
identique dans les deux ; `c` se place où on le veut, comme dans
`poly(...) + c`.

Chaque modification d'un champ génère une nouvelle formule, comme **En générer
un autre**, et **Copier** la met dans le presse-papiers pour Fractale → Formule
utilisateur. Tant qu'une longueur est vide, ou qu'un champ n'est pas un nombre,
est hors de son intervalle ou a son minimum au-dessus de son maximum, il n'y a
pas de formule, seulement une ligne qui dit pourquoi.

Chaque ligne de la liste des fonctions dit ce que prend l'appel et dans quel
ordre, avec des noms qui disent à quoi sert chaque place — "base, exponent",
"z, low, high", "before, after, n" — plutôt qu'avec des lettres, qui ne
disaient pas quelle place était laquelle ; ce sont les mêmes noms que montre la
ligne sous la barre de formule. Une convention les traverse tous : `z` est la
valeur sur laquelle travaille une fonction, quoi qu'on y ait réellement écrit.
La liste et ces nombres d'arguments sont vérifiés par un test par rapport à la
table propre de l'analyseur, si bien qu'une fonction qui gagne ou perd un
argument ne peut pas continuer à être décrite avec les anciens. Un argument
qu'on peut omettre est montré entre crochets avec la valeur qu'il prend alors :
`z, [power=1], [angle_factor=1]`.

**Les arguments sont séparés par une virgule, et par rien d'autre.** Un
point-virgule servait aussi : l'analyseur le changeait en virgule avant de lire
quoi que ce soit, donc `LOGN(5;Z^2)+C` et `LOGN(5,Z^2)+C` étaient une seule
formule, et il en allait de même des deux parties d'un nombre complexe. C'est
fini. Un point-virgule est désormais refusé là où il se trouve — *Opérateur
invalide : ;* — et la virgule est le séparateur partout, comme dans tout le
reste de la syntaxe.

### La sonde de palette

**Vue → Sonde de palette (cliquer sur un point)** transforme le pointeur en
sonde : cliquez n'importe où sur la fractale et la ligne de messages dit **de
quel endroit de la palette ce point a été tracé**, comme un nombre plutôt que
comme une couleur — *Place dans la palette 12.375*, et rien d'autre.

L'éditeur de palette pose **31 couleurs**, et la palette se construit en allant
de chacune à la suivante **huit cellules à la fois**. Le nombre est donc la
cellule divisée par huit : `12.000` est exactement la douzième couleur de
l'éditeur, `12.375` est aux trois huitièmes du chemin de la douzième à la
treizième. Une cellule est un huitième de couleur, plus fin que ce que
quiconque lit dans l'éditeur.

* Avec la **palette personnelle** en usage, le nombre est ramené entre **`0` et
  `31`**, la numérotation de l'éditeur lui-même.
* Avec toute autre palette, il est **laissé où il tombe et va bien au-delà de
  31** — une palette de 65534 cellules atteint 8191. C'est le nombre hors
  limites qui dit que la palette n'est pas celle que trace l'éditeur.
* Un point qui a pris la **couleur intérieure** donne *Place dans la palette :
  aucune* : cette couleur est un ton uni qu'aucune cellule de la palette ne
  représente.

Elle sert moins à savoir « quelle est cette couleur » qu'à voir **quelle part
de la palette l'image utilise réellement**, et la réponse est souvent
édifiante. Avec la coloration `iter` simple à 60 itérations, les cellules vont
de 1 à 59 — de `0.125` à `7.375` —, donc l'image ne touche jamais les couleurs
8 à 30 de l'éditeur et les modifier ne change rien. Une coloration qui s'étale,
`iter + real` par exemple, parcourt au contraire tout l'anneau, plusieurs fois.

Elle lit le point et non l'écran, elle est donc exacte : la même couleur peut
se trouver dans deux cellules, et une image en vraies couleurs mélange deux
cellules voisines, donc aucune couleur ne désigne une cellule à elle. Comme le
zoom par sélection, la sonde rend le pointeur au mode Julia rapide tant que
celui-ci est actif.

**La barre.** Fractale → Formule utilisateur ouvre un dialogue d'environ la
moitié de la largeur d'un écran de 1280×720 — quelque six cent cinquante
pixels, moins sur un écran plus petit — dont la barre passe à la ligne quand la
formule ne tient pas sur une seule. Elle passe à la ligne **entre les noms,
jamais à l'intérieur d'un nom** : un nom ou un nombre (`randsctile`, `z_p3`,
`0.6`) reste entier, et une ligne peut finir partout ailleurs — après une
parenthèse, une virgule, un opérateur. Le retour à la ligne propre à Qt ne
convenait pas : il coupe là où le permettent les règles de la prose, ce qui
dans une formule revient à après quelques virgules et guère plus, puis coupe un
nom là où il le faut. Seul un nom plus large que toute la barre est jamais
coupé. La barre grandit d'une ligne à la fois avec la formule, et dès que le
dialogue ne tiendrait plus sur un écran 720p, elle cesse de grandir et défile.
Pour tout le reste, cela reste une formule sur une ligne : Entrée vaut OK et
non une nouvelle ligne, et dans une formule collée les sauts de ligne et
tabulations deviennent des espaces, que l'analyseur lit comme rien.

**Historique...** sous la barre ouvre les formules utilisées auparavant — les
dix plus récentes, la plus récente d'abord, chacune une fois —, que l'ancienne
barre d'une ligne gardait dans sa liste déroulante. Une formule trop longue
pour la fenêtre est listée par son début, trois points et sa longueur,
`poly(z, {0.01097,0.00699}, {-0.05... (157 caractères)`, pour que deux qui
commencent pareil se distinguent encore ; le début s'arrête entre des noms,
comme une ligne de la barre, et on en voit d'autant plus que la fenêtre est
large. Double-cliquez sur une formule, ou sélectionnez-la et appuyez sur
Utiliser, et elle va entière dans la barre. Une formule est ajoutée quand le
dialogue est validé avec elle, et seulement alors : charger une position qui
contient une formule utilisateur, en rejouer une ou en rendre une avec
`-render` exécute la même commande `usrform` mais ne la pousse plus dans
l'historique, où chacune prenait la place d'une que vous aviez tapée.

**Flèche haut et flèche bas** parcourent les mêmes formules sans quitter la
barre. Flèche haut, depuis la première ligne de la barre, amène la formule
utilisée avant celle affichée ; flèche bas, depuis la dernière ligne, revient
vers la plus récente puis à ce que la barre contenait avant que vous ne
commenciez. Les deux petites flèches à côté de la barre font de même, chacune
grisée quand il n'y a plus rien de ce côté. Une formule déjà dans la barre est
sautée, et ce que vous tapez devient le point d'où l'on repart.

La dernière formule utilisée apparaît en pâle dans la barre : entière tant que
la barre est vide, et sa suite après ce que vous avez tapé tant que ce que vous
tapez en est le début. Tab la met entière dans la barre. Dès que ce que vous
tapez prend un autre chemin, le texte pâle disparaît, et Tab passe aux boutons
comme avant. La barre s'ouvre vide, avec cette formule en pâle, et OK sur une
barre vide l'utilise ; seule la toute première fois, avant qu'aucune formule
n'ait été utilisée, elle s'ouvre sur la formule en cours, `z^2+c`.

**Vue → Affichage état** montre la formule de la même façon : sur une ligne si
elle tient dans la fenêtre, sur deux sinon — la seconde sous la première,
passant à la ligne entre les noms — et quand deux ne suffisent pas non plus, la
seconde se termine par trois points et sa longueur, `... (918 caractères)`. Les
lignes en dessous descendent pour faire de la place, et tout est recomposé
quand la fenêtre change de largeur. C'était une seule ligne qui débordait du
bord de l'image.

**Sous la barre, une ligne dit ce que prend l'appel dans lequel vous êtes**,
avec l'argument que vous écrivez en gras. Placez le curseur dans `randsc(` et
on lit *randsc(**seed**, [size=1+i], [degradation=0.5+0.5i], [kaleidoscope=1],
[kaled_mirror=0], [skew=0])* ; dépassez une virgule et le gras vous suit. Elle
n'a besoin de la formule que d'une parenthèse ouverte précédée d'un nom, elle
est donc là quand la formule est à moitié écrite, c'est-à-dire quand on en a
besoin. Placez le curseur hors de tout appel et la ligne est vide.

Six des positions fournies étaient écrites à l'ancienne et ont été réécrites
avec des virgules : `circle`, `heart`, `helloween`, `pentafrac`, `warriormask`
et `burnship`, toutes sous **Fichier → Ouvrir → examples/Malczak**. **Une
position à vous qui a un point-virgule dans sa formule ne se chargera pas**
tant que la formule ne sera pas écrite avec des virgules — le message nomme le
caractère, et la solution est d'ouvrir la formule et de les remplacer.

**Laisser une place vide.** Un argument montré entre crochets peut être omis au
milieu d'un appel aussi bien qu'à la fin, en n'écrivant rien entre les deux
séparateurs : `julian(z, ,3)` donne le premier et le troisième et laisse la
fonction dire ce qu'est le deuxième. Les espaces ne changent rien — `f(z, ,5)`
et `f(z,,5)` sont le même appel — et un argument qui n'est pas entre crochets
doit quand même être écrit, donc `poly( ,1,2)` et `sin(,)` sont refusés au lieu
d'être devinés.

Ce que signifie une place vide regarde la fonction. La plupart prennent la
valeur par défaut qu'elles déclarent ; un coefficient de `poly` laissé vide est
un terme absent, donc `poly(z,1, ,1)` vaut `z^2 + 1` ; et une branche de
`ifiter` ou `ifiterl` laissée vide répète la précédente, c'est ainsi qu'on
donne à une branche plus d'une passe dans le cycle.
`ifiter(f(z), , ,g(z), , , , )` exécute `f` pendant trois passes et `g` pendant
cinq, et ne coûte pas plus à évaluer que de les écrire en entier : où mène
chaque choix est fixé une fois, quand la formule est analysée.

Une formule utilisateur pas encore écrite dit `z^2+c`. Elle disait le burning
ship, qui est une fractale à part entière et une énigme comme point de départ.

**`ifiterf(rest, last)`** — évalue `rest` à chaque passe sauf la dernière et
`last` à celle-ci. La dernière passe est la dernière que permet la limite
d'itérations : une formule n'a aucun moyen de savoir quelle passe sera celle
qui s'échappe, cela dépendant de la valeur qu'elle n'a pas encore produite.
Seule celle choisie est évaluée.

**`ifiterr(before, after, n)`** — évalue `before` tant que le numéro de passe
est inférieur à `n` et `after` à partir de `n`. Seule celle choisie est
évaluée, comme pour les deux autres, même s'il a d'abord fallu apprendre le
seuil à l'analyseur : un argument peut maintenant être marqué comme lu par le
sélecteur plutôt que choisi par lui, et il est alors évalué avant que le choix
soit fait.

**Les variables.** `z` est la valeur courante ; `c` est le point, qui est le
pixel en mode mandelbrot et la constante en mode julia ; `x` est la coordonnée
simple, qui est le pixel dans les deux modes et ne change pas d'une passe à
l'autre ; `n` est le numéro d'itération. Dans l'*Initialisation utilisateur*,
`z` est là où `z` aurait commencé sans initialisation, c'est-à-dire la même
valeur que `x`, donc une initialisation faite de `z` seul dit exactement ce que
ne rien dire dit.

Cette dernière partie est nouvelle, et pouvoir l'écrire l'est aussi : `z`
n'était enregistré dans aucun des deux analyseurs par lesquels passe
l'initialisation, et `x` dans le mauvais, donc une initialisation qui nommait
l'un ou l'autre était refusée par le thread qui devait la calculer -- tandis
que le dialogue, dont l'analyseur a les deux, l'acceptait sans montrer qu'elle
était ensuite jetée. `n` restait à ce que le pixel précédent y avait laissé,
autrement dit une initialisation qui le lisait lisait l'ordre dans lequel les
pixels étaient calculés.

**`p1`, `p2` ... `p9999`** — la valeur qu'avait `z` à une passe antérieure :
`p1` à la passe d'avant celle-ci, `p2` à celle d'avant, `p9999` neuf mille neuf
cent quatre-vingt-dix-neuf passes en arrière. `p` est un autre nom pour `p1`.
Avant qu'il y ait eu autant de passes, ils valent ce par quoi l'historique
commence, que décide *Fractale → Fixer les valeurs p à la première itération* :
le point itéré quand c'est activé, zéro sinon.

Ils s'arrêtaient à `p6` et étaient conservés en décalant d'un cran un tableau
de six places à chaque passe — une copie par place par passe, payée que la
formule en nomme ou non. C'est pourquoi ils s'arrêtaient à six : à neuf mille,
cela aurait fait 320 Ko de copie par passe, mesurés à vingt-quatre fois le coût
total d'itérer `z^2+c`. L'historique est maintenant un anneau, écrit une fois
par passe, et seules les places qu'une formule nomme réellement y sont lues,
donc `p9999` coûte ce que coûte `p1` et une formule qui n'en nomme aucun ne
coûte rien. La fenêtre de référence les liste parmi les variables.

**Suffixes sur une variable.** `z`, `c`, `x`, `p` et de `p1` à `p9999` prennent
des suffixes qui valent pour les appels faits le plus souvent sur eux :

| suffixe | vaut pour |
| --- | --- |
| `_b` | `bship(…)` |
| `_bi` | `bshipi(…)` |
| `_br` | `bshipr(…)` |
| `_pM` | `parchment(…, M)` |
| `_paM` | `parchmenta(…, M)` |

Ils se lisent de gauche à droite, chacun enveloppant ce qu'ont produit les
précédents : `z_p3` est `parchment(z,3)`, `c_b_p2` est `parchment(bship(c),2)`
et `p12_p2_p3` est `parchment(parchment(p12,2),3)`. *M* fait partie du nom, il
s'écrit donc en chiffres — un nombre entier, un ou plus ; `z_p0` est refusé
avec un message qui le dit. Pas sur `n`, qui est un compte et n'a pas de
seconde composante sur laquelle l'un d'eux pourrait travailler.

Ils sont développés en appels avant que la formule soit lue, donc une variable
suffixée trace exactement ce que tracent les appels écrits en entier et coûte
ce qu'ils coûtent ; la formule est enregistrée telle qu'elle a été écrite. Tout
ce qui n'est pas l'un des cinq est laissé tel quel et refusé comme le nom
inconnu qu'il est.

**Supprimés.** `powi`, `powdc` et `logcn` étaient des seconds noms de `pow` et
`logn` ; `rad`, `deg` et `sign` étaient listés sans aucune implémentation
derrière, et un second `trunc` était masqué par celui qui fonctionnait. Une
formule qui utilise un alias supprimé doit utiliser le nom restant — la
fonction est identique.

**`erf(z)`** — la fonction d'erreur dans le plan complexe. Précise à environ
trois ulp à l'intérieur d'un bailout de deux, là où une fractale itère
réellement.

**`poly(z, k1, k2, ..., km)`** — un polynôme en `z` :

    k1*z^(m-1) + k2*z^(m-2) + ... + k(m-1)*z + km

Le premier coefficient écrit multiplie la plus haute puissance et le dernier
est seul, donc l'appel se lit dans l'ordre où l'on énonce le polynôme. Calculé
par la méthode de Horner, soit m−1 multiplications au lieu des m(m−1)/2 qu'il
faudrait pour élever chaque puissance séparément, et de surcroît la plus
précise des deux.

    poly(z,1,0,0)+c            le Mandelbrot, écrit en entier
    poly(z,1,0,0,{0.7,0.2})    z^3 + 0.7+0.2i

**`randsc(seed, size, degradation, kaleidoscope, kaled_mirror, skew, skew_mode, selfsim)`**
— un bruit cohérent sur le point, qui donne des taches plutôt que de la neige
pixel par pixel. `size` (par défaut `1+i`) est la largeur moyenne d'une tache
le long de l'axe réel et sa hauteur le long de l'axe imaginaire. C'est un
**bruit de gradient**, une octave du mouvement que somment les [modes
fbm](#mouvement-brownien-fractionnaire) et `fbm()`. C'était un bruit de valeur,
et pour la même raison que le leur il traçait son réseau : les taches
devenaient plates le long de chacune de ses lignes et sortaient carrées. Les
taches sont aussi larges qu'avant — le réseau est aux trois cinquièmes de la
taille, et le champ cesse de se ressembler à moitié aux onze seizièmes de la
taille, avant comme après — et les valeurs s'étalent comme avant, 0.2145 depuis
le milieu, la courbe de contraste du fbm ramenant l'étalement plus étroit du
bruit de gradient sur celui du bruit de valeur. Mais chaque valeur est
différente, et un appel coûte 167 ns contre 147. Les cinq mosaïques sont leur
réseau et sont inchangées au bit près. `degradation` (par défaut `0.5+0.5i`,
qui les divise par deux à chaque passe) les rétrécit au fil de l'itération : la
taille est multipliée par lui à chaque passe, composante par composante, donc
`0.5+0.2i` sur `1+i` donne `1+i` à la première passe, puis `0.5+0.2i`, puis
`0.25+0.04i`. Un zéro dans une composante de l'un ou l'autre argument renvoie
zéro au lieu de diviser par zéro. Viennent ensuite les deux arguments du
kaléidoscope, `skew` (par défaut `0`), `skew_mode` et `selfsim`, chacun décrit
plus bas. Seule la valeur initiale est obligatoire.

### L'inclinaison, et pourquoi ces champs sont difficiles à colorer

Le moteur colore avec les deux composantes de l'orbite : `real` et `zmag` en
lisent une, `imag`, `angle` et `real / imag` l'autre. Ces champs renvoient un
seul nombre pour toute une cellule, sur l'axe réel, si bien que chaque mode
trace **un ton par cellule** — `zmag` et `iter + real` ressemblent à la valeur
tronquée, et les trois autres n'ont rien du tout à lire. Mesuré sur neuf cents
pixels : de 14 à 126 valeurs pour les deux premiers, et exactement **une** pour
les autres.

Faire varier la valeur à l'intérieur d'une cellule ne peut pas se faire avec un
seul nombre réel. Le même nombre décide si le point s'échappe, donc une valeur
qui bouge à travers une cellule fait sortir la moitié de la cellule et laisse
l'autre moitié dedans, et la cellule sort **coupée en deux par une ligne
droite**. Cela a été essayé deux fois et a coupé les mosaïques les deux fois.

`skew` le fait plutôt avec la seconde composante, et il **tourne** la valeur au
lieu de la mettre à l'échelle. De combien il tourne suit **l'écart du point au
milieu de sa cellule** : avec `out` allant de zéro au milieu à un au bord, la
valeur est tournée de deux fois l'arc tangente de `skew_re*out + skew_im`.

`out` est mesuré **dans la géométrie propre de la cellule**, donc ce que la
couleur trace dans une cellule est la forme en laquelle le champ est découpé,
et non un ensemble de lignes posées en travers — zéro au milieu d'une cellule,
un tout le long de son bord :

| champ | ce que la couleur trace dans une cellule |
| --- | --- |
| `randsc` | ses propres taches — la rotation suit le niveau, donc les contours sont ceux du champ |
| `randscq` | des carrés, autour du milieu du carré |
| `randsch` | des hexagones, autour du milieu de l'hexagone |
| `randsct` | des triangles, autour du milieu du triangle |
| `randscp` | son propre polygone, rétréci pas à pas — il n'y a pas deux cellules de même forme |

La **partie réelle** de l'inclinaison est ce gradient à travers une cellule. La
**partie imaginaire** est une rotation uniforme que partage toute la cellule :
elle décale la couleur d'une cellule sans rien y tracer, donc `0.4+0.2i` donne
à la fois les contours de la forme et un décalage par cellule.

### Ce que trace l'inclinaison : `skew_mode`

Le septième argument dit **ce qui varie à l'intérieur de la cellule**. Les
nombres **s'additionnent** au lieu de s'exclure, donc un appel peut en demander
deux et obtenir les deux.

| | |
| --- | --- |
| **1** (et `0`) | **la forme** — le contour propre de la cellule, rétréci pas à pas. Ce que reçoit un appel qui ne donne pas `skew_mode`, et ce qui est utilisé quand ni 1 ni 2 n'est nommé |
| **2** | **la rosette** — l'angle autour du milieu de la cellule, replié en autant de tours que le kaléidoscope a de coins, pour que chaque cellule porte la symétrie de l'image elle-même. De simples rayons quand rien n'est replié |
| **4** | **à travers le coin** — une rotation qui suit la place du point dans son coin du kaléidoscope, la même dans chaque coin |
| **8** | **radial** — le module bouge aussi bien que l'angle |

`3` est la forme et la rosette à la fois, ce qui fait une spirale. `11` est ces
deux-là avec le module qui bouge aussi.

**`4` suit la géométrie propre du kaléidoscope.** Il tourne la valeur de rien
sur les bords de chaque coin et de toute la rotation que donne à elle seule la
partie imaginaire — deux fois son arc tangente, 127 degrés pour `2` — le long
de la ligne par son milieu, et entre les deux de la rotation de la partie
imaginaire fois la part traversée. Chaque coin reçoit la même chose, donc un
kaléidoscope à six reste un kaléidoscope à six, miroirs compris, et la rotation
est continue à travers les bords et les milieux, là où le pli reflète. Il ne
fait rien quand rien n'est replié, ou quand la partie imaginaire de
l'inclinaison est nulle. Pas de trigonométrie ; il ne coûte rien à 64 bits de
mantisse et un dixième d'appel à 113.

C'était d'abord **une rotation différente dans chaque coin** — la partie
imaginaire fois le numéro du coin plus un — destinée à ce que les copies
cessent d'être identiques. C'est la seule chose que les copies d'un
kaléidoscope ne peuvent pas être : une image d'ordre six sortait en six
tranches sans rapport, avec une coupure droite à chaque jointure. Le multiple
butait aussi sur le demi-tour, si bien que quatre tranches sur six se
trouvaient presque sur l'axe réel, et là `real / imag` les traçait comme de la
neige, quelle que soit la moyenne que `selfsim` faisait du bruit.

Ce qu'il en reste est une propriété de `real / imag`, qui divise par la partie
imaginaire de *z* : partout où une rotation la fait passer par zéro, la
coloration saute d'un bout de la palette à l'autre, et les derniers restes des
passes fines décident du côté pixel par pixel. Une rotation qui balaie le coin
y passe pour la plupart des réglages. Sur un `randsct` d'ordre six coloré par
`real / imag`, avec `selfsim` à un, une partie imaginaire de `0.25` ne laissait
aucun grain et `2` en laissait des bandes près des bords ; une coloration qui
ne divise pas par une composante n'en a pas.

**`8` est celui qui répond à `zmag`.** Une rotation laisse le module où il est,
et `zmag` et le bailout lisent le module et rien d'autre — aucun autre mode ne
peut donc les atteindre. Le prix est celui que l'inclinaison devait éviter :
l'échappement bouge avec le module, donc la figure est coupée partout où le
bailout tombe dans une cellule. Associez-le à un bailout circulaire et à une
coloration qui lit les composantes séparément, ou acceptez la coupure à
dessein.

Une inclinaison nulle laisse chaque mode sans effet, donc `skew_mode` seul ne
change aucune image.

On avait d'abord écrit une rampe droite — `skew_re*du + skew_im*dv` à travers
la cellule — et c'est ce que ceci remplace. Elle colorait, mais traçait les
mêmes bandes diagonales, verticales ou horizontales dans chaque cellule de la
même façon, quel que soit le découpage du champ, et les bandes gâchaient
l'image.

À **`0`**, ce que reçoit un appel qui ne le nomme pas, la rotation est nulle,
la partie imaginaire reste à zéro, et chaque valeur est le nombre qu'elle a
toujours été, **au bit près** — les sommes de contrôle de référence d'avant
l'existence de l'argument passent inchangées, et une position enregistrée se
rend à l'identique. La position `talc`, trois appels de bruit et un `c` sur
345600 valeurs, donne la même signature que la version d'avant l'argument.

**Combien de couleur elle donne.** L'indice du moteur est
`(iter + quantity) * speed + shift`, donc un écart de un dans la grandeur vaut
une bande à une vitesse de un. L'écart de `imag` à l'intérieur d'une cellule :

| inclinaison | 0.02 | 0.05 | 0.1 | 0.2 | 0.4 | 0.6 | 1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `randsc`, `randsct` | 0.08 | 0.21 | 0.43 | 0.83 | 1.5 | 1.9 | 2.1 |
| les trois mosaïques | 0.06 | 0.14 | 0.28 | 0.53 | 0.96 | 1.2 | 1.4 |

Donc **de `0.3` à `0.6` apparaît une bande à une vitesse de un**, ce qui
manquait à ces champs. Compter les valeurs distinctes dirait que `0.01` suffit
— neuf cents valeurs — mais mille valeurs dans un centième de bande sont toutes
une seule couleur.

**Ce qu'elle coûte.** La mise à l'échelle déplaçait le module de la valeur, et
un bailout rond regarde exactement cela : une cellule dont le niveau était
proche de l'endroit où tombe le bailout sortait coupée en deux, et les
mosaïques perdaient leurs formes. Une rotation laisse le module où il est, donc
**avec un bailout circulaire l'inclinaison est gratuite à toute intensité** —
sur 90000 pixels, à 0.05, 0.3 et 1, et sur les cinq champs : pas un pixel ne
s'échappe différemment.

Un polygone de bailout est une autre affaire, car il ne lit pas le module mais
les composantes, et une rotation promène la valeur sur un cercle qui peut
traverser un côté. De combien il peut le traverser, c'est de combien les coins
de la forme dépassent ses côtés, et les chiffres le suivent — la part de
l'image dont l'échappement change, sur les cinq champs :

| forme du bailout | coins sur apothème | inclinaison 0.01 | inclinaison 0.05 | inclinaison 0.3 | inclinaison 1 |
| --- | --- | --- | --- | --- | --- |
| cercle | 1.00 | **0%** | **0%** | **0%** | **0%** |
| hexagone | 1.15 | **0%** | 0.0–0.8% | 0.7–2.2% | 0.4–2.4% |
| carré | 1.41 | **0%** | 0.0–1.0% | 0.8–2.8% | 1.2–12.8% |
| triangle | 2.00 | 0.05–1.3% | 0.7–5.0% | 2.6–6.4% | 7.4–21.6% |

Le **nombre** auquel est réglé le bailout le déplace aussi, et pour la même
raison : ce qui peut changer de côté est une valeur dont le module tombe entre
l'apothème du polygone et ses coins, et la place de cet anneau parmi les
cellules, c'est le nombre qui la décide.

Donc avec un polygone les deux souhaits tirent en sens contraire, et il n'y a
pas moyen d'arranger cela autrement : l'échappement lit les deux composantes de
la valeur, et la coloration aussi, donc tout ce qui donne à la couleur quelque
chose à lire est quelque chose que l'échappement peut voir. **Associez
l'inclinaison à un bailout circulaire** et elle ne coûte rien à aucune
intensité. Avec un polygone, soit vous acceptez ces quelques pour cent, soit
vous gardez l'inclinaison à `0.05` et augmentez plutôt la vitesse de
coloration.

Ce qu'une rotation ne peut pas toucher, c'est **`zmag`**, qui lit le module.
Colorez avec quelque chose qui lit les deux composantes séparément :
`iter + real`, `iter + imag`, `angle`, `real / imag` à l'extérieur, et `real`
ou `real / imag` dans la coloration intérieure.

Deux choses en découlent. La valeur est complexe tant que l'inclinaison n'est
pas nulle, donc `randsc(7,,,,,0.4)*z` tourne `z` en plus de le mettre à
l'échelle. Et le kaléidoscope replie dans les deux cas : mesuré sur cinq
champs, deux, trois, cinq et six coins et les deux miroirs, une rotation d'un
coin laisse chaque valeur où elle était.

**`randscq(...)`** — le même champ sans l'interpolation : une mosaïque de
cellules carrées plates au lieu de taches. Mêmes arguments, même sens.

**`randscp(...)`** — encore le même champ, les courbes ôtées mais pas
l'irrégularité : un germe est semé dans chaque cellule et chaque position prend
la valeur du germe le plus proche, ce qui trace un diagramme de Voronoï — des
polygones convexes plats à bords droits, pas deux de même forme.

**`randsch(...)`** — la même chose sur des hexagones : un nid d'abeilles de
cellules plates. Des polygones réguliers qui pavent le plan, c'est celui qui
n'a pas de direction privilégiée — chaque cellule a les mêmes six voisines aux
mêmes six angles — donc il se lit comme une matière et non comme une grille.

**`randsct(...)`** — des triangles équilatéraux, d'orientation alternée. Une
mosaïque triangulaire a bien une direction privilégiée, et c'est ce qu'on
demande en la choisissant.

**`randsctile(tiling, seed, ...)`** — le même champ sur l'un de **quarante-cinq
pavages**, que choisit le premier argument ; tout ce qui suit est ce que prend
la famille, décalé d'une place, donc `selfsim` est le neuvième. Chaque pavage
est mis à l'échelle d'une tuile d'aire unité en moyenne, donc le premier
argument change la forme des cellules et non leur taille. Un nombre hors de 1 à
45 trace zéro.

| | |
| --- | --- |
| 1–3 | les pavages réguliers : carrés, triangles, hexagones |
| 4–11 | les huit archimédiens, polygones réguliers aux sommets tous semblables : 4.8.8 (octogones et carrés), 3.6.3.6, 3.4.6.4, 3.12.12, 4.6.12, 3.3.3.4.4, 3.3.4.3.4 (carré adouci), 3.3.3.3.6 (hexagonal adouci) |
| 12–19 | leurs duaux, qui ne sont pas réguliers : tétrakis carré, rhombille (les cubes), deltoïdal trihexagonal (cerfs-volants), triakis triangulaire, kisrhombille, et les pentagones prismatiques, du Caire et en fleur |
| 20–32 | briques, appareil flamand, arête de poisson, vannerie, pythagoricien (carrés de deux tailles), chevrons, carrés et losanges, maisons (un pentagone, debout et suspendu tour à tour), rangées de carrés et de triangles en deux rythmes, hexagones parmi des triangles, croix grecques, tétrominos en T |
| 33–37 | étoiles islamiques : à huit branches avec des croix, à six branches avec des hexagones, à huit branches tirées de 4.8.8, à douze branches tirées de 3.12.12 et de 4.6.12 |
| 38 | cellules de Voronoï, irrégulières — celles de `randscp`, avec une valeur propre |
| 39–43 | pavages qui ne se répètent jamais, tirés des multigrilles de de Bruijn : les losanges de Penrose, les cerfs-volants et fléchettes de Penrose, les carrés et losanges d'Ammann–Beenker, et des losanges dans douze et dans sept directions |
| 44–45 | par substitution : le moulinet, dont les triangles sont tournés dans toutes les directions possibles, et la chaise, des L découpés en L |

L'**onglet Pavages** d'Aide → Référence des formules utilisateur a une image de
chacun, légendée comme la liste l'onglet Valeurs, et un double-clic sur l'une
copie le début de son appel — `randsctile(18, ` — dans le presse-papiers, prêt
pour la valeur initiale. Les images sont tracées par `randsctile` lui-même,
huit unités de côté autour de l'origine, par `tools/randsctile-thumbnails.cpp`,
et conservées dans l'exécutable au lieu d'être tracées à chaque ouverture de
l'onglet. Des images conservées peuvent devenir périmées, donc chacune est
stockée avec une empreinte de ce que traçait la fonction, et un test retrace
les mêmes points et échoue, en nommant le pavage, tant que
`cmake --build <build directory> --target randsctile-thumbnails` ne l'a pas
retracé.

Les périodiques sont des **tables**, écrites par `tools/randsctile-tables.py`,
qui construit chaque pavage à partir de sa géométrie et le vérifie — vingt
mille points au hasard, dont chacun doit tomber dans exactement une tuile —
avant de l'écrire. Les étoiles islamiques sont tirées des pavages archimédiens
: chaque polygone de huit côtés ou plus devient une étoile dont les pointes
sont sur les sommets du polygone, et le mince triangle entre chaque encoche de
l'étoile et le côté va à la tuile de l'autre côté, donc deux étoiles partageant
un côté forment un losange entre elles et un petit polygone prend une encoche
sur chaque côté qu'il partage avec une étoile. Les tuiles peuvent être concaves
— étoiles, croix, le T — donc un point est placé dans sa tuile par le nombre de
croisements plutôt que côté par côté. Le point est découpé en périodes entières
dans la précision de travail, comme `randscq` le découpe en cellules, donc un
pavage périodique est exact aussi loin qu'on aille.

Ceux qui ne se répètent jamais **n'ont pas de table** et sont recalculés en
chaque point. Les losanges de Penrose, Ammann–Beenker et les losanges à douze
et à sept directions viennent de la **multigrille de de Bruijn** : des familles
de droites parallèles dans un second plan, une tuile partout où deux se
croisent, désignée par les deux droites — quatre nombres entiers, exacts quel
que soit celui qui les demande. Le point est placé en essayant les croisements
voisins de l'endroit où il tombe dans ce plan, d'abord les familles dont les
droites passent le plus près, ce qui le trouve en quatre à huit essais en
moyenne. Les cerfs-volants et fléchettes viennent des triangles de Robinson
découpés niveau par niveau dans une roue de dix autour de l'origine, le
moulinet d'un triangle rectangle découpé en cinq, la chaise d'un L découpé en
quatre ; à chaque niveau une droite ou deux décident dans quel enfant se trouve
le point.

**Ce qu'il coûte**, mesuré par appel face aux 220 ns de `randscq` sur cette
machine : les pavages périodiques de 255 à 380 ns, ceux qui ne se répètent
jamais de 410 à 680. Ils sont calculés en `double` dans les deux versions,
place dans une période comprise — en long double ceux qui ne se répètent jamais
prenaient trois fois plus de temps — et ainsi les deux versions ont tracé les
quarante-cinq pavages à l'identique au pixel près, quarante mille points
chacun. Au-delà de 2³² cellules de l'origine, où `double` commencerait à perdre
leurs tuiles, les pavages qui ne se répètent jamais deviennent plats comme la
famille au-delà de sa grille ; les périodiques tiennent aussi loin que le reste
de la famille.

L'inclinaison mesure une tuile comme `randscp` mesure son polygone, par rapport
au rayon du plus grand cercle que contient la tuile, donc les contours sont la
tuile rétrécie, concaves comprises ; la rosette tourne autour du milieu de la
tuile.

    z+randsctile(18,13,{0.3,0.3},{1,1})*0.25        pentagones du Caire
    z+randsctile(36,13,{0.3,0.3},{1,1})*0.25        étoiles à douze branches
    z+randsctile(39,13,{0.3,0.3},{1,1},5)*0.25      Penrose, replié en cinq

**`fbm(value, seed, [intensity=4], [frequency=8], [octaves=4], [roughness=0.5])`**
— un mouvement brownien fractionnaire : du bruit de gradient additionné en
octaves, chacune à deux fois la fréquence de la précédente et gardant
`roughness` de sa hauteur — le même mouvement que tracent les modes de
coloration de ce nom, voir [Mouvement brownien
fractionnaire](#mouvement-brownien-fractionnaire). Ce qu'il trace est de
l'**usure** plutôt qu'un motif — taches et bosses — parce qu'aucune octave
n'est assez grande pour se voir seule et aucune assez petite pour disparaître.

Là où `randsc` lit la position et ne peut lire que la position, **celle-ci lit
ce que vous écrivez devant elle**. C'est tout l'intérêt de l'avoir sous forme
de fonction :

    fbm(z,7)          bouge avec l'orbite
    fbm(x,7)          reste immobile sur le plan
    fbm(z*3+c,7)      quoi que ce soit
    z^2+c+fbm(z,7,0.1)*i    l'ensemble avec un bord rugueux

La valeur et la valeur initiale sont obligatoires ; les autres ont des valeurs
par défaut. Le mouvement va **de zéro à `intensity` et jamais en dessous**,
donc l'ajouter à quelque chose ne peut pas faire passer cela sous zéro. La
valeur initiale se lit comme `randsc` en lit une, donc le même nombre désigne
le même champ dans les deux. `octaves` est maintenu entre un et vingt-quatre,
et il ne vaut la peine de l'augmenter que si l'on augmente `roughness` : avec
la moitié simple, la huitième octave porte un deux cent cinquantième du tout.

**Il n'a pas de kaléidoscope propre**, contrairement à la famille `randsc`,
parce qu'il y en a déjà un avec lequel le composer : `parchment(z, sectors)`
replie l'angle de `z` en `|sectors|` secteurs et `parchmenta` reflète une
moitié sur deux, donc `fbm(parchmenta(z,6),7)` le dit, et le dit là où le
lecteur de la formule peut le voir.

Les modes de coloration sous **Fractale → Coloriage extérieur → Mouvement
brownien fractionnaire...** appliquent le même mouvement à une couleur.
Celui-ci le place là où une formule peut s'en servir.

| | |
| --- | --- |
| `randsc` | doux, courbe |
| `randscq` | dur, régulier, carré |
| `randscp` | dur, irrégulier |
| `randsch` | dur, régulier, hexagonal |
| `randsct` | dur, régulier, triangulaire |
| `randsctile` | dur, n'importe lequel de quarante-cinq pavages |

Les six posent une cellule sur chaque carré unité de la taille en vigueur, donc
`size` veut dire la même chose partout : changer une lettre change la forme des
cellules et non l'échelle de l'image. Un hexagone de rayon circonscrit un
couvre 2.6 carrés unité et un triangle de côté un en couvre 0.43, donc ces deux
grilles sont mises à l'échelle pour correspondre aux autres au lieu d'être
laissées telles quelles.

Chaque passe a son propre champ. L'itération entre dans le hachage avec la
position, donc une formule qui appelle l'une d'elles une fois par itération
obtient une nouvelle valeur à chaque fois — même avec une dégradation de un,
qui laisse la taille telle quelle. La dégradation fixe la taille des taches à
une passe donnée ; ce n'est pas ce qui rend la passe différente.

Il y a un plancher sous la taille. Une cellule se trouve en divisant la
position par la taille et le résultat doit tenir dans un entier, donc la taille
ne peut utilement descendre sous environ la position divisée par ce que
contient un entier — un milliardième de milliardième d'elle. Une dégradation
d'un demi y arrive en une soixantaine de passes. Au-delà, il ne reste aucune
structure de cellules à résoudre et le champ est une seule valeur plate sur
tout le plan ; il change encore à chaque passe et diffère encore entre les six
fonctions, donc une formule qui en soustrait une d'une autre ne s'arrête pas
sur zéro en itérant jusqu'à la limite pour rien. Pour un fondu qui reste une
image jusqu'au bout, utilisez une dégradation proche de un.

Toutes hachent la position du point et l'itération, jamais `z`, et jamais un
état global. Donc la même image revient à chaque nouveau tracé, avec n'importe
quel nombre de threads, et en mode Mandelbrot comme en mode Julia. `randsc` est
continu, donc les deux précisions concordent à environ 1e-19 ; les mosaïques
sont des fonctions en escalier et peuvent différer d'un cheveu le long des
bords des cellules, ce qui est inhérent au fait de demander des bords nets.

Une formule qui appelle l'une d'elles désactive la détection de frontière.
Cette optimisation remplit sans la calculer une région qu'elle a trouvée
entourée d'une seule couleur — vrai pour une fractale, faux pour un champ de
bruit.

L'ancien `rand` est inchangé et dépend toujours du nombre de fois qu'il a été
appelé, donc la même position ne se retrace pas de la même façon. Préférez
`randsc`.

### Kaléidoscopes

Les deux arguments qui suivent la dégradation replient le plan avant que le
champ y soit échantillonné. Le premier est le nombre de coins en lesquels le
replier, `1` — la valeur par défaut — le laissant tel quel ; le second,
`kaled_mirror`, est le miroir qui fait le pli :

| kaled_mirror | |
| --- | --- |
| `0` (et toute autre valeur) | la moitié éloignée de chaque coin reflète la moitié proche, donc chaque coin est symétrique par rapport à sa propre bissectrice |
| `1` | la même chose dans l'autre sens, la moitié proche reflétant l'éloignée |

    randsc(13,{0.6,0.6},{1,1},6,0)     six coins, chacun miroir de lui-même
    randsc(13,{0.6,0.6},{1,1},3,1)     trois coins, repliés dans l'autre sens

Les deux plis sont continus là où les coins se rejoignent, donc le bruit reste
cohérent et les deux précisions continuent de concorder ; un pli qui se
rejoindrait lui-même de façon inégale se verrait comme une couture.

C'est la seule partie de la famille qui coûte de la trigonométrie, et seule une
valeur de deux ou plus l'atteint. Un appel qui n'en dit rien paie une
comparaison.

### Tirer une image du bruit

`randsc` renvoie une valeur dans `[0, 1)`. Utilisé seul, l'itération ne quitte
jamais un bailout de 4, donc rien ne s'échappe et l'image est plate.
Multipliez-le, ou laissez-le perturber une itération qui s'échappe vraiment :

    z+randsc(13,{0.35,0.35},{1,1})*0.25       des taches
    randsc(13,{0.15,0.15},{0.5,0.5})*2.5      mouvement brownien
    randscq(13,{0.15,0.15},{0.5,0.5})*2.5     carrés épars, de plus en plus petits
    z+randscp(13,{0.3,0.3},{1,1})*0.25        polygones irréguliers
    z+randsch(13,{0.3,0.3},{1,1})*0.25        un nid d'abeilles
    z+randsct(13,{0.3,0.3},{1,1})*0.25        des triangles
    z^2+c+randsc(13,{0.25,0.25},{1,1})*1.2    l'ensemble lui-même déformé

Une dégradation de 0.5 divise les taches par deux à chaque passe, donc après
vingt itérations elles font un millionième de leur taille et sont plus petites
qu'un pixel. Pour un fondu lent sur un long calcul, prenez quelque chose proche
de 1 — `0.97^60` vaut encore 0.16.

Une valeur initiale entière est exacte dans les deux versions ; une
fractionnaire est quantifiée à ses bits de tête, qui concordent sauf pour
environ une valeur sur cinquante millions.

### Toutes les passes jusqu'ici : `selfsim`

Chaque passe de ces fonctions est un champ à part — l'itération entre dans le
hachage — dont les cellules valent la dégradation fois celles de la passe
précédente. Une formule qui en appelle une à chaque passe lit donc à chaque
fois un champ nouveau et plus fin, et un mode de coloration qui lit la dernière
passe ne lit que le dernier champ. Avec une dégradation de 0.75, les cellules
d'une vue large de cinq passent sous le pixel vers la seizième passe, et
au-delà l'intérieur est de la neige : cent passes de
`randsct(672,,{0.75,0.75},,,{2,2},6)+c` en mode Julia écartent 70 % des pixels
intérieurs voisins d'un seizième de la palette ou plus.

Le huitième argument renvoie à la place **les passes jusqu'ici, moyennées**,
qui sont entre elles comme les octaves d'un mouvement brownien fractionnaire :
la passe *n* pèse d^(nH), où *d* est la dégradation — la moyenne géométrique de
ses deux composantes prises sans leur signe, donc `0.75+0.75i` vaut 0.75 — et
*H* est l'argument. Les poids sont divisés par leur total, donc la réponse
reste parmi les valeurs moyennées : dans [0, 1] sans inclinaison.

| selfsim | |
| --- | --- |
| omis | désactivé — la seule passe, au bit près. Ce que reçoit un appel qui ne nomme rien, et une place laissée vide |
| `0` | la moyenne simple : chaque passe pèse autant, et c'est ce que donne tout *H* proche de zéro — pas désactivé |
| `1` | le mouvement simple : chaque passe pèse *d* fois la précédente |
| `0.5` | plus rugueux : les passes fines gardent plus de leur poids, donc plus de détail apparaît |
| plus grand | plus lisse : les premières passes prennent presque tout |
| négatif | les passes fines l'emportent sur les grossières, et la neige revient |
| complexe, `{1,2}` | la partie réelle pondère comme ci-dessus, la partie imaginaire tourne les passes — voir plus bas |

    randsct(672,,{0.75,0.75},,,{2,2},6,1)+c     le même, autosimilaire

Sur cette formule la neige disparaît : les pixels intérieurs voisins sont en
moyenne à 1.8 cellule de palette au lieu de 35, et 0.6 % d'entre eux à plus de
seize au lieu de 70 %. Ce qui reste, ce sont les triangles à toutes les
échelles à la fois, les grands portant les petits.

* **`selfsim` peut être complexe.** Alors d^(nH) l'est aussi : d^(n·Hr) en
  grandeur et une rotation de n·Hi·ln d, donc chaque passe est tournée d'autant
  de plus que la précédente puis moyennée comme le fait un *H* réel — une
  spirale d'octaves. À 0.75, une partie imaginaire de un fait environ seize
  degrés par passe. Les poids sont divisés par **la somme de leurs grandeurs**,
  non par leur somme complexe : cela garde la réponse dans les limites de la
  plus grande des valeurs moyennées. La somme complexe approcherait de zéro à
  une certaine passe et ferait exploser toutes les valeurs à cette passe pour
  tous les pixels à la fois, les poids n'ayant rien à voir avec le pixel. Les
  rotations donnent à la valeur une partie imaginaire même sans inclinaison,
  donc `imag` et `angle` ont quelque chose à lire. Une dégradation de un a un
  logarithme nul et ne tourne rien ; `{1,0}` est `1` au bit près ; et `{0,2}`
  est chaque passe pondérée également, et tournée.
* **Zéro n'est pas désactivé.** d^(n·0) vaut un à chaque passe, donc `0` est la
  moyenne simple, et tout *H* proche aussi — des deux côtés, et complexe
  également. Zéro était d'abord l'interrupteur, ce qui en faisait la seule
  valeur par laquelle la courbe ne passe pas : `0.000001` était la moyenne
  simple et `0` la dernière passe seule, à un quart de l'intervalle d'écart. La
  dernière passe seule est l'extrémité de la courbe, *H* tendant vers moins
  l'infini. Pour ne pas avoir de moyenne, omettez l'argument, ou laissez sa
  place vide.
* **C'est la valeur inclinée qui est moyennée**, donc l'inclinaison et chaque
  bit de `skew_mode` continuent de faire ce qu'ils faisaient, et la réponse est
  complexe quand l'inclinaison la rend telle.
* **Une moyenne rapproche les valeurs.** À 0.75 et *H* égal à un, l'écart est
  d'environ deux cinquièmes de celui d'une seule passe, donc augmentez la
  vitesse de coloration pour retrouver le contraste.
* **Une dégradation de un** pondère chaque passe également et donne leur
  moyenne simple, qui se pose vers un milieu plat à mesure que les passes
  s'additionnent. C'est ce que dit la définition, et c'est fait comme il est
  dit.
* **Les passes qu'un appel n'a pas vues sont calculées quand on l'interroge.**
  Un appel dans une branche que la formule ne prend qu'à partir de la septième
  passe reçoit quand même les sept premières, et la réponse à une passe est la
  même quel que soit le chemin qui y a mené. Un nouveau pixel se reconnaît à sa
  position autant qu'à sa passe : la passe seule, sur laquelle se fonde `trap`,
  permettrait à un tel appel de poursuivre la moyenne laissée par le pixel
  précédent.
* **Rien ne déborde.** La moyenne est tenue comme une seule valeur, chaque
  passe en prenant une part qui découle de la part de la passe précédente,
  plutôt que comme une somme et un total — donc un *H* négatif, ou une
  dégradation supérieure à un, dont les poids dépasseraient ce que contient le
  type, n'a pas où le faire.
* **Cela coûte** environ un dixième de plus par appel à 64 bits de mantisse et
  un septième à 113. Désactivé, rien de mesurable.

### Des figures au lieu du bruit

`sierpinskyt`, `sierpinskyc` et `snowflake` se tiennent là où se tient le
bruit, et rien en elles n'est aléatoire : la figure est là où elle semble être,
et le mode Mandelbrot et le mode Julia tracent la même image, une forme du plan
étant au même endroit dans les deux cas.

**Ce sont des fractales à part entière, pas des champs à multiplier dans une
autre.** Écrite seule — `sierpinskyt()` et rien d'autre — chacune trace sa
figure, comme Fractale → Autres formules trace son **Sierpinski**, son
**Sierpinski Carpet** et son **Koch Snowflake**.

Toutes trois le font en **portant le point vers son parent**. Chaque partie de
l'une de ces figures en a un : un trou du triangle de Sierpinski est dans un
trou plus grand un niveau au-dessus, un trou du tapis dans la cellule découpée
de la même façon, un triangle du flocon sur un côté de l'hexagone du milieu ou
sur le bord libre d'un autre triangle. Le pas vers le parent est un mouvement
du plan — un doublement en s'éloignant du coin le plus proche pour le triangle,
un agrandissement d'une cellule par le nombre de cellules pour le tapis. La
partie supérieure n'a pas de parent dans la figure, donc son pas la porte hors
du bailout, et une partie *n* niveaux plus bas met *n* pas à y arriver : **la
passe où un point sort est le niveau où il se trouve**, et le nombre
d'itérations est l'image.

Ce qui revient à la passe où un point sort compte autant que la passe elle-même
: **tous les modes de coloration extérieure sauf le nombre d'itérations le
lisent**. Donc la partie supérieure de chaque figure est jetée dehors au lieu
d'être simplement déclarée sortie — le triangle double son trou central en
l'éloignant d'un coin et il atterrit hors du triangle, le tapis met à l'échelle
sa cellule centrale par rapport à la cellule voisine, le flocon jette son
hexagone central dans la direction vers laquelle regarde ce secteur. Chacun
emporte avec lui d'où il venait, donc `real`, `imag`, `angle`, la coloration
lissée et le reste ont tous de quoi travailler. Seul un point qui n'a jamais
été dans la figure reçoit un nombre qui signifie « sorti ».

Un flocon se lit **à partir d'un hexagone vers l'extérieur**, et c'est ce qu'il
faut savoir de lui. Il se décompose exactement en un hexagone régulier au
milieu, six triangles de même côté que cet hexagone posés sur ses six côtés,
douze d'un tiers de cela sur leurs bords libres, quarante-huit d'un neuvième
sur les leurs, et ainsi de suite sans fin. Trois des six sont les coins du
triangle dont la figure est née et trois sont les premières bosses posées sur
ses côtés, et **rien ne les distingue** : de l'hexagone vers l'extérieur, un
flocon a une symétrie d'ordre six.

Lu plutôt à partir du premier triangle, cinq huitièmes de la figure sont le
niveau un — un seul triangle plat remplissant l'image, avec le flocon seulement
dans la frange autour, et rien à dire pour la couleur extérieure. Lu à partir
de l'hexagone, les deux premiers niveaux font cinq douzièmes chacun et les
autres décroissent régulièrement. Un triangle posé sur l'hexagone passe sur le
plus grand triangle que contient l'hexagone — coin sur coin, trois fois l'aire
— donc il est agrandi de la racine carrée de trois et tourné d'un douzième de
tour ; un triangle posé sur le bord libre d'un autre est agrandi de trois et
tourné vers la direction où regarde ce bord.

**Chaque pas dilate**, et c'est une règle plutôt qu'un hasard. Le triangle
double, le tapis agrandit une cellule par le nombre de cellules, le flocon de
trois ou de la racine de trois. Un triangle posé sur l'hexagone pourrait plutôt
être rabattu à plat sur le côté où il se pose, ce qui le pose sur l'un des six
dont l'hexagone est fait et serait plus propre ; mais un rabattement est une
isométrie, et une isométrie est insensible à l'endroit d'où le point est parti.
Écrite dans une formule plus grande, la figure donnerait alors à cette formule
une image rigide d'elle-même sur cinq douzièmes de son aire, et la formule y
tracerait une région plate. Un test vérifie les trois dilatations.

Le point voyage vraiment, et c'est ce qui les rend utilisables et pas seulement
correctes : `z` en chemin est un point du plan comme un autre, donc les modes
de coloration qui le lisent, la coloration lissée et l'écriture d'une figure
dans une formule plus grande ont tous un sens.

Un triangle de Sierpinski et un tapis remplissent leur forme, donc chaque point
de l'un d'eux a un niveau et sort à ce niveau. Un flocon ne remplit pas son
hexagone : il laisse six coins de **fond**, et ce fond ne fait pas partie de la
figure et n'a pas de niveau. Il est **tourné d'un demi-tour**, ce qui envoie du
fond sur du fond, donc il ne sort jamais et apparaît dans la **couleur
intérieure** — l'espace vide reste vide au lieu de prendre à la figure une
bande de la couleur extérieure.

Tourné, et non rendu là où il était, car **un point fixe est une orbite arrêtée
et non une orbite bornée**. Immobile, la boucle exécutait chaque passe jusqu'à
la limite sans que rien change ; les modes de coloration intérieure qui
comparent une passe à la précédente n'avaient rien à comparer ; et
`snowflake()` écrite dans une formule plus grande rendait à cette formule le
point qu'elle tenait déjà, alors que `sierpinskyt()` et `sierpinskyc()`
déplacent chaque point qu'on leur donne. Le demi-tour est préféré au sixième de
tour, que la symétrie de la figure permettrait aussi : tous deux envoient du
fond sur du fond, mais le demi-tour garde le point dans un carré et un cercle
autant que dans l'hexagone. Il laisse intacte la distance au centre, donc les
modes de coloration intérieure qui la lisent ne sont pas touchés, et l'image
que trace `snowflake()` est inchangée au pixel près.

C'est le seul endroit où ces figures coûtent quelque chose. On pose à un pixel
de fond la question d'appartenance une fois par passe aussi longtemps qu'on le
regarde, alors que les deux autres figures et le corps du flocon sortent à la
passe que nomme leur niveau. Tracé directement par la boucle d'itération, un
flocon entier coûte environ 2.4 fois ce qu'il coûtait quand le fond recevait
une bande à lui. La détection de frontière reste activée pour ces formules —
leurs bandes sont vraiment pleines, ce qui est la condition pour l'employer —
et le fond est une seule région pleine, donc ce qui arrive à l'écran est moins
que ce que suggère ce chiffre.

Écrites dans une formule plus grande — `snowflake()^-2+snowflake()` et autres —
les trois ne se comportent **pas** de la même façon, et la raison mérite d'être
dite clairement, car aucun travail sur la figure ne la changera.

Un triangle de Sierpinski est une **poussière** : il n'a pas d'aire, ses points
ne sortent jamais, et chacun d'eux est doublé à chaque passe. Une formule bâtie
dessus a quelque chose de chaotique en dessous partout, à toutes les
profondeurs. Un flocon est **plein** : dix douzièmes de lui sont l'hexagone et
les six triangles, qui sortent à la première et à la deuxième passe, donc une
formule bâtie dessus manque de figure presque tout de suite et trace ce qu'elle
trace d'elle-même. Le détail est le long du bord de Koch, où les niveaux vont
loin, et les régions plates sont les intérieurs, où ils ne vont pas.

**Ce que sait une passe, c'est ce que cette passe a payé.** Le flocon décide si
un point est du fond en descendant la courbe de Koch, et cette descente a droit
à **un niveau par passe**. Elle parcourait auparavant ses vingt-quatre niveaux
en une seule évaluation, donc une formule écrite autour de la figure recevait
tout le bord à sa première passe et le traçait en entier si peu de passes qu'on
lui donne — à `Itérations : 2` on obtenait encore chaque bosse de la courbe de
Koch, ce qu'aucune autre formule de XaoS ne fait et à quoi ne sert pas un
nombre d'itérations. Manquer de niveaux n'est pas la même chose qu'être dehors
: un point hors du cadre ou au-dessus de l'enveloppe est du fond et on le sait
déjà, tandis qu'un point indécis est porté plus loin comme partie de la figure
et interrogé de nouveau à la passe suivante avec un niveau de plus à dépenser.
L'image que trace `snowflake()` elle-même est inchangée — un point de la figure
n'a jamais eu besoin de la descente profonde, seulement de l'assurance de ne
pas être du fond.

Ce qui reste, c'est la différence honnête entre les deux figures, et aucun
travail sur l'une ou l'autre ne la comblera. Un triangle de Sierpinski est une
**poussière** : pas d'aire, ses points ne sortent jamais, chacun doublé à
chaque passe, donc une formule bâtie dessus a quelque chose de chaotique en
dessous à toutes les profondeurs. Un flocon est **plein**, et dix douzièmes de
lui sont l'hexagone et les six triangles, qui sortent à la première et à la
deuxième passe. Mesuré à un seul zoom — la part des pixels situés sur un bord
de bande — `sierpinskyt()^-2+sierpinskyt()` atteint 17 pour cent vers douze
passes et `snowflake()^-2+snowflake()` 2.5 ; tous deux croissent avec le nombre
de passes, ce qui est le but, et là où ils s'arrêtent, c'est ce que sont les
figures.

`radius` veut dire ce que veut dire `bailout` et se lit de la même façon :
**comme le carré de la distance**. Ce qu'il trace est la forme que trace un
bailout de ce nombre, **exactement inscrite dedans** — réglez le bailout sur la
forme correspondante et le même nombre, et la figure et le bailout se
superposent :

| | |
| --- | --- |
| `sierpinskyt` | le triangle du bailout **triangle −90°**, coin pour coin |
| `sierpinskyc` | le carré du bailout **carré**, côté pour côté |
| `snowflake` | l'hexagone du bailout **hexagone 0°**, ses six pointes sur les six coins de l'hexagone |

Ce qu'il faut savoir, c'est qu'un polygone de bailout place ses **côtés** à la
racine carrée du bailout depuis le centre — c'est son apothème, pas son rayon
circonscrit — donc ses coins sont plus loin que ne le dit le nombre : deux fois
plus loin pour un triangle, et d'un septième pour un hexagone. `sierpinskyt(4)`
atteint donc 4 par ses coins tandis que `sierpinskyc(4)` atteint 2 par ses
côtés, et les deux ont raison.

Réglez le bailout plus grand que le rayon et la figure ne remplit plus la forme
: il y a une **marge** entre les deux, et réglez-le plus petit et la figure est
coupée. Les deux sont ce qui a été demandé et aucun n'est corrigé. Mais la
marge ne fait pas partie de la figure et n'a pas de niveau, donc on ne lui en
donne pas : elle est **tournée selon la symétrie propre de la figure** — un
tiers de tour pour un triangle, un quart pour un tapis, un demi-tour pour un
flocon — ce qui la garde hors de la figure et dans le bailout, donc elle ne
sort jamais et est tracée dans la **couleur intérieure**, les modes de
coloration intérieure travaillant sur elle. Elle sortait auparavant à la
première passe et prenait la bande la plus basse de la couleur extérieure, ce
qui traçait toute la marge d'un seul ton plat et ôtait une bande à la figure.

Une rotation et non un miroir, bien qu'un miroir servirait aussi bien à la
garder où elle est : ce qu'une figure rend doit tourner avec le point, pour que
tourner l'image tourne ce qui y est tracé, et un miroir le tourne dans l'autre
sens. Pour la même raison, le milieu d'un tapis — le bloc de cellules que
l'anneau du bord enferme, qui n'a pas de parent et doit donc être jeté dehors —
est jeté dans la direction vers laquelle regarde son propre **quart** du carré,
les quarts découpés par les diagonales. Jeté toujours dans la même direction,
comme avant, il traçait ce milieu comme une rampe d'un côté à l'autre sans rien
en travers ; le carré a une symétrie d'ordre quatre, et ce qu'on y trace
devrait l'avoir aussi. Un test vérifie que chaque figure répond à une rotation
de son propre ordre par la même rotation.

**`sierpinskyt([radius=4], [kaleidoscope=1], [kaled_mirror=0])`** — le triangle
de Sierpinski, dans ce triangle.

**`sierpinskyc([radius=4], [squares=3], [kaleidoscope=1], [kaled_mirror=0])`**
— le tapis de Sierpinski, dans ce carré. Le carré est découpé en `squares` sur
`squares`, l'anneau de cellules le long du bord est gardé, tout ce que cet
anneau enferme est jeté, et on fait de même pour chaque cellule gardée. L'image
est donc **un carré au milieu et 4×`squares`−4 autour** : huit pour trois, ce
qui est le tapis tel qu'on le dessine d'habitude, douze pour quatre, seize pour
cinq. Deux est le seul nombre sans véritable anneau, et là c'est le coin
éloigné qui s'en va, ce qui redonne un triangle de Sierpinski — un carré
découpé en quatre dont on retire un coin, c'est exactement ce qu'est un
triangle de Sierpinski.

**`snowflake([radius=4], [kaleidoscope=1], [kaled_mirror=0])`** — le flocon de
Koch, ses pointes sur les coins de cet hexagone, en bandes par génération du
milieu vers l'extérieur : l'hexagone à la première passe, les six triangles sur
ses côtés à la deuxième, les douze sur leurs bords libres à la troisième, les
quarante-huit à la quatrième. Les six coins de fond que la figure ne couvre pas
n'ont aucune bande — ils ne sortent jamais, et sont tracés dans la couleur
intérieure.

Chaque argument a une valeur par défaut, donc `snowflake()` est un appel, et
`sierpinskyc( ,5)` aussi — un tapis plus ajouré à la taille par défaut.

**Le kaléidoscope est celui de la famille du bruit**, écrit en dernier pour
chacune des trois et avec le même sens : en combien de coins le plan est
découpé, et quel miroir les replie. `1` — ce que reçoit un appel qui ne dit
rien — ne replie rien et laisse la figure le nombre qu'elle a toujours été, au
bit près.

Ce qu'il replie est pourtant différent, et doit l'être. Le bruit replie la
**position**, qui reste immobile toute l'orbite ; une figure n'a que faire de
la position et lit **z**, le point qu'elle porte vers le bas jusqu'à la partie
où il se trouve, donc c'est z qui est replié. L'image sort de la même manière —
la figure tracée dans un coin et répétée autour de l'origine — et elle tient
jusqu'en bas, chaque pas étant fait sur le point replié et replié de nouveau à
la passe suivante. Dans le demi-coin que le pli laisse tel quel, la figure est
la figure.

    sierpinskyt(4,6,0)      le triangle six fois autour de l'origine
    sierpinskyc(4,3,5,1)    le tapis en cinq coins, le miroir proche
    snowflake(4,5)          cinq coins de flocon

Un flocon a déjà la symétrie d'ordre six avec le miroir le long de chaque
bissectrice, donc le replier ainsi en six coins est l'identité sur l'image :
`snowflake(4,6,0)` sort **pixel pour pixel** identique à `snowflake(4)` sur
90000 pixels. Demandez cinq coins, ou quatre, pour le voir se replier.

Une chose à savoir sur le tapis : ce qu'il découpe, il le jette dans l'une de
quatre directions, les quarts découpés par les diagonales du carré, et sur une
diagonale même le choix est une égalité. Un pli peut envoyer toute une
demi-droite du plan sur cette diagonale, et là les deux côtés de l'égalité sont
jetés dans des directions opposées — un trait fin le long de la diagonale, qui
est la marche que le tapis a toujours eue et non quelque chose que fait le pli.
Mesuré sur dix mille points à six coins : chaque désaccord se trouvait sur une
diagonale, et les points repliés eux-mêmes concordaient à une partie sur 10^16
près.

    sierpinskyt()        avec Fractale -> Forme de l'échappement -> triangle -90
    sierpinskyc()        avec forme du bailout carré
    sierpinskyc( ,5)     un tapis plus ajouré
    snowflake()          avec forme du bailout hexagone 0

Aucune ne coûte plus que le bruit à côté, passe pour passe : mesurés par
rapport à `randsc`, le triangle 0.44, le tapis 0.54 et le flocon 0.50 de son
coût à 64 bits de mantisse, et 0.50, 0.87 et 0.65 à 113. Lire un flocon depuis
son milieu l'a rendu moins cher en plus de plus beau : six septièmes de la
figure sont l'hexagone ou un triangle posé dessus, et on y répond par quatre
multiplications et une comparaison sans parcourir la courbe du tout. Ce que
coûte un flocon sur une image entière, c'est le paragraphe ci-dessus, et cela
tient au nombre de passes que prend le fond plutôt qu'au coût d'une passe. Ce
n'était pas gratuit, et à 113 bits ce fut presque perdu — une racine carrée et
une division y sont toutes deux logicielles, et chacune coûte à peu près une
figure entière. Donc l'inverse de la racine carrée du rayon est gardé au point
d'appel au lieu d'être recalculé à chaque passe, et la racine elle-même revient
en multipliant le rayon par cet inverse plutôt qu'en divisant ; le triangle
compare ses trois poids barycentriques mis à l'échelle par un nombre positif,
ce qui ne change ni lequel est le plus grand ni lequel est négatif ; le tapis
compte en cellules plutôt que dans le plan, où la multiplication par `squares`
et la division par lui s'annulent ; et le flocon prend l'enfant vers lequel
passer à partir du bord sous la bosse duquel il a été trouvé, qu'il devait de
toute façon trouver. Il parcourt en outre la courbe de Koch en `double`,
puisque la figure a une taille fixe et que 52 bits d'une fraction d'un bord
sont plus que ce que peuvent utiliser les 24 niveaux qu'il trace. Un test
vérifie les rapports, donc une réécriture qui les perd le signale.

## Palettes

Fractale → Couleurs demande un numéro d'algorithme, et dit après le numéro
comment s'appelle chacun. 1 à 3 sont ceux de XaoS ; 4 à 7 ont été refaits après
avoir mesuré de quoi sont faits les trois premiers, et 8 après avoir mesuré des
tableaux.

| | |
| --- | --- |
| 1 | **Sombre et couleur** — un presque-noir avec une teinte à lui et une couleur, à tour de rôle |
| 2 | **Noir, couleur, blanc** |
| 3 | **Chaud sur vert sombre** — violets, rouges et orangés vifs sur des verts sombres et des bleu-vert, noir ou blanc un sur trois |
| 4 | **Smog** — sombre et mélancolique : gris, ocre, ardoise et violet poussiéreux, avec du sang-de-bœuf et de l'écarlate parmi eux |
| 5 | **Chaud sur nuit** — rose, rouge, bordeaux, brun, orange, ocre et jaune sur bleu-vert, pétrole, bleu et indigo, dans le squelette de la 3 |
| 6 | **Paires favorites** — les paires de couleurs que les trois premiers placent le plus souvent côte à côte, posées entières entre noir et blanc |
| 7 | **Couleurs aléatoires** — chaque étape une couleur vraiment au hasard, et rien d'autre |
| 8 | **Kandinsky** — les dégradés de ses tableaux, du mat au vif et de temps en temps presque jusqu'au gris |

### De quoi sont faits les trois premiers

Vingt mille palettes de chacun, faites comme le programme les fait, chaque
couleur nommée d'après la plus proche de trois douzaines de couleurs de
référence en CIELAB, en comptant les couleurs et les paires de couleurs :

* **Ancres.** Le noir et le blanc forment d'un tiers à la moitié des étapes, et
  les couleurs se placent entre eux avec une période de deux ou trois étapes,
  donc le motif se voit dès les premières bandes — qui, dans une palette longue
  d'environ trois mille étapes, sont tout ce qu'utilisent la plupart des
  images.
* **1 et 2 contiennent 256 couleurs et pas plus.** Ils prennent l'octet bas du
  générateur, et l'octet bas de ce générateur est à lui seul un générateur de
  période 256 : le rouge décide du vert et du bleu, et chaque couleur est
  toujours suivie de la même. Sur vingt mille palettes de 2, aucune couleur n'a
  été suivie une seule fois d'une autre. De là viennent leurs paires
  récurrentes — rose à côté de jaune vingt-quatre fois plus souvent que le
  hasard, bordeaux à côté d'olive six fois, rouge à côté de pétrole huit fois.
* **3 puise dans deux fenêtres** : des couleurs vives du violet en passant par
  le rouge jusqu'à l'orange, et des sombres du vert jusqu'au bleu-vert. Ses
  paires les plus courantes sont une couleur chaude et vive sur un vert sombre.

### Les quatre suivants

Chacun a été mesuré par rapport aux trois premiers avant d'être gardé : la
divergence de Jensen-Shannon entre les couleurs que montrent deux algorithmes,
par nom, et entre les paires de voisines qu'ils montrent, sur cinq mille
palettes — 0 pour la même palette, 1 pour rien en commun. 1 et 2 sont à 0.11
l'un de l'autre selon cette mesure, 1 et 3 à 0.16.

* **4, Smog**, est sombre et éteint là où 1 est vif sur noir. Les étapes
  profondes sont une couleur à peine au-dessus du noir, les claires un gris
  couvert avec une teinte et jamais blanc, et une étape sur quatre un gris
  simple. Ses teintes sont parcourues sur un anneau de 256 comme 1 et 2
  parcourent les leurs, donc les voisines reviennent comme chez eux, mais
  seulement parmi les teintes du smog — ocre et olive, et acier, ardoise et un
  violet poussiéreux — et, sur un quart de l'anneau, du cramoisi à l'écarlate :
  sang-de-bœuf quand c'est profond et un écarlate atténué quand c'est clair,
  gardé saturé, car un rouge sombre qui ne l'est pas se lit comme du brun. Une
  étape sur cinq se lit comme rouge, contre une sur trente pour 1. Il est à
  0.43 à 0.56 des trois premiers. Ce fut d'abord le mécanisme de 1 sur un autre
  anneau de 256 couleurs, choisi aussi loin de l'anneau classique qu'un anneau
  peut l'être, et il mesurait 0.05 de 1 — plus près que 2 : la ressemblance
  tenait au mécanisme, noir la moitié du temps et une couleur de n'importe où
  l'autre moitié, non aux couleurs.
* **5, Chaud sur nuit**, est le squelette de 3 avec la température inversée :
  des couleurs chaudes de toute luminosité, du rose en passant par le rouge, le
  bordeaux, le brun et l'orange jusqu'à l'ocre et au jaune, sur des couleurs de
  nuit gardées sombres, du bleu-vert en passant par le pétrole et le bleu
  jusqu'à l'indigo. Les deux fenêtres font 125 degrés de large, près des 135 de
  3 ; elles en faisaient d'abord le tiers et se lisaient comme deux couleurs.
  Le bleu marine est sa couleur la plus courante, et bleu marine avec olive,
  bordeaux ou brun ses paires les plus courantes.
* **6, Paires favorites**, prend chaque paire de couleurs que l'un des trois
  premiers met l'une après l'autre au moins trois fois plus souvent que le
  hasard — soixante-six, avec les deux couleurs telles qu'elles sont réellement
  sorties — et pose chaque paire entière, les deux couleurs se touchant, entre
  une étape noire et une blanche.
* **7, Couleurs aléatoires**, fait de chaque étape une couleur tirée au hasard
  du haut du générateur, là où 1 et 2 tirent du bas : soixante-dix-neuf mille
  couleurs différentes dans vingt mille palettes là où 2 en a 256, et aucune
  paire plus d'un cinquième plus fréquente que le hasard. Il avait d'abord du
  noir et du blanc une étape sur trois, et mesurait 0.00 de 2 : le rythme du
  très sombre et du très clair était tout ce qu'on voyait. Sans lui, il a moins
  d'étendue du sombre au clair que ceux qui sont ancrés, ce qui est le sens
  même d'aléatoire.

Les quatre qu'ils remplacent étaient un spectre, une bichromie, une triade et
une paire complémentaire, raisonnés à partir du cercle chromatique, et à
l'écran chacun était une seule couleur. Ils répartissaient leurs teintes sur
toute la palette, et la palette fait dans le programme environ trois mille
étapes — quatre à neuf dans le test censé s'en apercevoir, ce qui explique
qu'il ne l'ait pas fait — si bien qu'une image ne voyait que la première
teinte. Et ils gardaient chaque étape comme une couleur, profonde ou pâle, sans
noir ni blanc, donc rien en eux n'avait de bord. Le test fait maintenant ses
palettes comme le programme, et demande aux sept que la partie que voit une
image aille du sombre au clair et contienne plus d'une couleur ; on demande un
peu moins à 8, venu ensuite, et la raison est plus bas.

Une position enregistrée avec les numéros 4 à 7 avant ce changement revient
avec d'autres couleurs : le numéro et la valeur initiale sont tout ce qu'une
position enregistre de sa palette, et ce que signifie le numéro a changé. Une
position enregistrée avec les numéros 1 à 3 n'est pas touchée, au bit près, et
une enregistrée avec les numéros 4 à 7 nomme un algorithme que le XaoS
d'origine n'a pas et refusera.

### 8, Kandinsky

Fait à partir de deux séries d'images dans sa manière : quatre vignettes vives
— noir, vermillon, jaune, bleu céruléen, violet et rose sur papier crème — et
une image mate, orange brûlé, ocre, brun, sauge et pétrole sur le même papier.
Les vignettes donnent dix-huit couleurs, par k-moyennes en CIELAB, et chacune a
une **forme mate** : la couleur de l'image mate la plus proche en luminosité et
en teinte, sans la saturation, puisque c'est par la saturation que les deux
diffèrent. Le vermillon devient orange brûlé, le jaune ocre, le bordeaux brun,
le noir un brun foncé, le céruléen sauge, le bleu pétrole. L'image mate n'a pas
de gris, donc le gris reste gris. Le bleu-vert d'eau qu'elle semble avoir est
un sauge, 128,146,121 : il ne se lit comme tel qu'à côté de l'orange.

Chaque palette tire une **humeur**. Quatre sur cinq tombent également entre les
formes mates et les vives ; la cinquième dépasse les formes mates vers le gris,
jusqu'à quinze centièmes de leur saturation, et une palette sur vingt environ
est presque grise.

**Les ancres viennent par hasard.** Une palette sur dix peut changer çà et là
une étape en noir pur ou en blanc pur : chaque étape après la première par son
propre tirage, à un taux que la palette tire entre une étape sur vingt et une
sur dix, noir ou blanc selon le tirage — sans rythme et sans alternance,
contrairement aux ancres de 1 à 3. Les neuf autres palettes sur dix sont
exactement ce qu'elles étaient avant l'arrivée des ancres : le tirage était
déjà fait pour chaque étape, et laissé inutilisé.

Ce qui suit une couleur, c'est là où la mènent **ses dégradés**. Partout où la
peinture des vignettes change doucement — sur huit pixels, sans pas plus grand
qu'un tiers du changement — on a compté les couleurs aux deux extrémités, et
chaque étape tire la suivante de ces comptes : le papier se poursuit en un
lavis de jaune puis en jaune, le noir en vermillon (trois fois sur dix) et le
vermillon en bordeaux, le bleu en mauve et le mauve en rose, le céruléen en un
bleu-vert d'eau grisé. Les bords nets ont été comptés aussi, et écartés : la
plupart sont du papier à côté de papier d'une autre nuance, et les palettes
qu'on en faisait sortaient papier et gris. Une palette commence sur le papier
ou sur le noir.

Sur vingt mille palettes, celles du dixième supérieur des humeurs sont à 0.06
des vignettes, selon la même divergence que plus haut, et celles du dixième
au-dessus de zéro à 0.07 de l'image mate ; elles sont à 0.39 à 0.67 des sept
autres palettes.

Comme un dégradé reste un moment dans une partie du cercle chromatique, et que
8 n'a pas de noir et de blanc toutes les quelques étapes sur lesquels se
replier, les tests lui demandent moins qu'aux autres. Là où la palette est
étalée longue, quarante entrées par étape et plus, une image voit une dizaine
d'étapes, et une palette sur cent est alors un seul champ — ocre avec du sauge,
pétrole avec du vermillon — sans atteindre un tiers du chemin du noir au blanc
ni une seconde couleur. Et trois étapes, ce qu'est souvent une palette courte
de 256 entrées, forment un seul dégradé.

## Observer l'orbite

`trap` et `stripe` gardent un nombre sur toute une orbite plutôt que sur le
point, ce qui est le genre de grandeur que les modes de coloration du moteur ne
peuvent atteindre : une boucle de calcul y est compilée pour rendre une
couleur, et seuls le dernier `z` et celui d'avant survivent pour servir à la
coloration. Fait dans la formule, cela ne coûte rien au moteur.

Toutes deux rendent la valeur qu'elles ont reçue jusqu'à la dernière passe que
permet la limite d'itérations, et à celle-là ce qu'elles ont recueilli. Donc
une formule entière est

    trap(z^2+c, 3)

et rien d'autre : la fractale itère comme elle le ferait, et à la dernière
passe la valeur devient le piège, que tracent ensuite les modes de coloration
intérieure. Un point qui s'échappe n'atteint jamais cette passe et garde sa
couleur extérieure ordinaire, donc ce qu'elles tracent, c'est l'intérieur.

**`trap(orbit, shape, centre, size)`** — le plus près que l'orbite soit jamais
venue d'une forme. `shape` vaut 0 par défaut, `centre` l'origine, `size` 1.

| shape | |
| --- | --- |
| `0` | le centre lui-même, un point |
| `1` | une droite horizontale qui le traverse |
| `2` | une verticale |
| `3` | les deux, une croix |
| `4` | un anneau de rayon `size` |
| `5` | un carré de demi-côté `size` |
| `6` | un losange de demi-diagonale `size` |

**`stripe(orbit, density)`** — la moyenne de
`(sin(density * arg orbit) + 1) / 2` le long de l'orbite. `density` vaut 4 par
défaut et c'est le nombre de rayures qui font un tour ; un nombre entier, sinon
elles ne se rejoignent pas là où le tour se ferme. Une moyenne sur l'orbite
change doucement avec le point même là où le nombre d'itérations saute, et
c'est ce qui trace les fibres qui ont fait connaître la méthode.

    trap(z^2+c,0)                    à quelle distance de l'origine l'orbite est passée
    trap(z^2+c,3)                    une croix, qui trace des rayons
    trap(z^2+c,4,{0,0},{0.5,0})      des anneaux dans la cardioïde
    stripe(z^2+c,6)                  six rayures par tour

La grandeur cumulée vit au point d'appel, donc deux pièges dans une formule
gardent chacun la leur et un thread ne peut en déranger un autre. Les deux
désactivent la détection de frontière, pour la même raison que le bruit : deux
voisins qui prennent le même nombre de passes peuvent tout de même avoir vu des
orbites bien différentes, et une région remplie sans être calculée serait
fausse.

Ce ne sont pas des modes de coloration, et ils ne peuvent pas l'être : les
modes vivent dans la boucle compilée de chaque formule, qui rend directement
une couleur, et il y en a trente et une, chacune avec sa table de symétrie par
position. Ce qui est ici ne fonctionne que pour une formule utilisateur.

## Forme de l'échappement

Calcul → Valeur d'échappement → Forme de l'échappement. L'itération s'arrête
quand `z` quitte une région, qui a toujours été un cercle. La forme se voit
dans les bandes hors de l'ensemble, c'est donc autant un outil de dessin qu'un
outil numérique.

| forme | s'échappe quand |
| --- | --- |
| Cercle (classique) | \|re\|² + \|im\|² dépasse le bailout |
| Carré | l'une des deux composantes seule le fait |
| Losange | \|re\| + \|im\| le fait |
| Axe réel | la composante réelle seule |
| Axe imaginaire | la composante imaginaire seule |
| Les deux axes | les deux composantes à la fois |
| Triangle 0° / 90° / 180° / −90° | le point quitte un triangle de cette orientation |
| Hexagone 0° / 90° | le point quitte un hexagone de cette orientation |
| Octogone | le point quitte un octogone |

Les polygones sont mesurés par leur **apothème** — la distance du centre au
milieu d'un côté — donc un polygone et le cercle du même bailout se touchent le
long des côtés et non aux coins. L'orientation fait partie de la forme : un
triangle tourné de 120° est le même triangle, donc 0°, 90°, 180° et −90° en
sont quatre différents, alors qu'un hexagone se répète tous les 60° et en a
donc deux, et un octogone tous les 45°, ce qui ne lui laisse rien à choisir. Le
triangle et l'hexagone ont donc un sous-menu à eux pour l'orientation. Leurs
normales sont calculées une fois par image, donc la boucle intérieure coûte de
trois à huit multiplications-additions au lieu d'un sinus et d'un cosinus par
côté.

Par rapport au cercle, le losange change environ 39 % des pixels et l'axe réel
44 %. La silhouette de l'ensemble ne bouge pas — c'est ce qui ne s'échappe
jamais qui la décide — donc regardez les bandes.

Elle n'est écrite dans une position enregistrée que lorsque ce n'est pas le
cercle, donc une position qui ne demande pas de forme reste chargeable par
toute version antérieure.

## Mouvement brownien fractionnaire

Une boîte de dialogue à lui de chaque côté — **Fractale → Coloriage intérieur →
Mouvement brownien fractionnaire...** et la même sous Coloriage extérieur —
contenant les trois modes qu'elle offre et les cinq nombres d'après lesquels
ils sont tracés. C'était un sous-menu, les trois modes en entrées de choix et
les nombres derrière une boîte de dialogue Paramètres en dessous ; c'est
maintenant une seule boîte de dialogue, cochée dès que l'un des trois modes est
utilisé, et OK règle le mode et les nombres ensemble en une seule étape
d'annulation. Une position enregistre toujours le mode comme `(incoloring N)`
ou `(outcoloring N)` et les nombres à part, comme `(infbmset ...)` ou
`(outfbmset ...)`, de sorte qu'un fichier d'une version antérieure s'ouvre sans
changement.

Ce qu'il fait, c'est donner à une coloration propre un air **usé** : les bandes
gardent leur forme et leur ordre et cessent d'être parfaites, comme si la chose
avait été laissée aux intempéries. Ce n'est pas un hasard pixel par pixel, qui
se lit comme du grain et a été essayé puis écarté ; c'est un champ de bruit
avec une échelle, donc ce qu'il trace, ce sont des marques.

**À l'extérieur** les trois modes sont `fbm + smooth`, `fbm + iter` et `fbm`
seul — le compte lissé usé, le compte en bandes usé, et le bruit seul, qui
trace des nuages autour de l'ensemble sans aucun compte dedans.

**À l'intérieur** ce sont `fbm`, `fbm + zmag` et `fbm + decomposition`. Le
premier est l'intéressant : l'intérieur est par défaut un seul ton plat, et
ceci lui donne une surface.

### Bruit de gradient, et aucune grille

Le mouvement est un **bruit de gradient** : chaque coin d'un réseau reçoit une
direction, et un point prend de chaque coin la pente que cette direction donne
à sa distance. C'était auparavant un bruit de valeur — un nombre à chaque coin,
fondu par une smoothstep — et cela traçait une grille. Une smoothstep n'a de
pente à aucune de ses extrémités, donc le champ devenait plat le long de chaque
ligne du réseau, et comme le réseau de chaque octave contenait les lignes de la
première, toutes y devenaient plates ensemble. En travers d'une ligne, la pente
était nulle ; l'image sortait en carrés, de huit pixels de côté dans
`FBM_ERROR.xpf`, comme si elle avait été calculée à une résolution plus basse
que celle de l'affichage. Le bruit de gradient est aussi raide sur les lignes
qu'entre elles — mesuré, 0.99 aux réglages par défaut et 1.08 à ceux de cette
position — et le réseau de chaque octave est décalé par rapport au précédent,
de quantités qu'aucun doublement ne rend jamais entières, donc deux octaves ne
partagent jamais une ligne ni un coin. `fbm()` dans les formules trace le même
mouvement, et `randsc` en est une octave.

Ce qui a été gardé : la valeur va toujours de zéro à l'intensité et jamais
au-delà ; les **marques ont la taille qu'elles avaient**, le réseau tournant
aux trois cinquièmes de la fréquence pour que le bruit de gradient les trace
aussi grandes que le faisait le bruit de valeur (mesuré par la distance dont le
champ se déplace avant de ne plus se ressembler qu'à moitié : 0.56 cellule aux
réglages par défaut, avant et après) ; le **contraste**, dont le bruit de
gradient seul n'a que sept dixièmes — chaque octave passe par une courbe qui
ramène ses valeurs sur la distribution du bruit de valeur, monotone pour ne
déplacer aucun contour, et mesuré sur le mouvement aux réglages par défaut il
oscille de 0.133 autour de son milieu contre 0.132 avant ; et le coût, à peu
près — les modes un peu plus rapides qu'avant, `fbm()` environ un huitième plus
lente. Ce qui ne l'a pas été : une position enregistrée avec un mode fbm garde
ses nombres, son échelle et sa force et **trace d'autres marques**.

### Les cinq nombres

La même boîte de dialogue que le mode, en dessous. Ils sont enregistrés avec la
position.

| | |
| --- | --- |
| **intensité** | de combien de bandes de couleur le bruit déplace la valeur. `4` est la valeur par défaut ; en dessous de un c'est un murmure, au-dessus de dix les bandes cessent d'être des bandes |
| **fréquence** | la taille des marques, aussi grandes que les tracerait un réseau d'autant de cellules par unité du plan. **C'est celle-ci qu'on augmente en zoomant** : à la première vue `8` se lit bien, à une étendue de `0.02` il faut environ `200` |
| **octaves** | combien on en additionne, chacune au double de la fréquence de la précédente |
| **rugosité** | ce que chaque octave garde de la hauteur de la précédente. `0.5` est le mouvement simple ; plus haut, c'est plus granuleux |
| **valeur initiale** | la même image à chaque fois. Changez-la pour une autre du même caractère |

**Les octaves et la rugosité se soutiennent l'une l'autre.** Avec une rugosité
de `0.5`, la huitième octave porte un deux cent cinquantième du tout, donc les
octaves au-delà de la quatrième ou de la cinquième ne changent rien de visible
— mesuré, quatre et dix-huit ont donné la même image. Augmentez la rugosité et
les octaves fines ont de quoi dépenser, et alors le nombre commence à compter.

### La taille des marques selon le compte

À l'extérieur, deux champs de plus au bas de la boîte de dialogue. **Mise à
l'échelle du bruit** dit ce qui décide de la taille des marques : **fixe**, où
elle commence, la fréquence seule, comme ci-dessus ; **itération**, la
fréquence relevée avec le nombre d'itérations ; **itération lissée**, la
fréquence relevée avec le compte lissé. **Facteur** dit de combien : la
fréquence double autant de fois toutes les dix itérations, si bien qu'à `1` les
marques sont deux fois plus petites là où l'orbite a mis dix itérations de plus
à sortir, et quatre fois là où elle en a mis vingt. À `0` le compte ne change
rien. Là où le compte monte, vers le bord de l'ensemble, les marques s'affinent
et suivent le détail.

Une fréquence qui change d'un pixel au suivant ne peut pas simplement être
posée sur le plan : le champ s'étirerait, d'autant plus qu'on s'éloigne de
l'origine. Le mouvement est donc lu aux deux puissances de deux de part et
d'autre du compte, chacune un champ à part à fréquence fixe, et les deux sont
mêlés comme une texture entre ses mipmaps, le mélange ramené au contraste d'un
seul champ. Avec **itération lissée** les marques changent de taille sans
couture ; avec **itération** elles en changent aux bords des bandes, chaque
bande avec un grain à elle. À l'intérieur le compte est partout le même, si
bien que l'intérieur n'a aucun des deux champs.

Entre deux doublements on lit deux champs, et la coloration coûte à peu près le
double de **fixe** : 248 ns par pixel contre 127, à quatre octaves. Là où le
compte tombe exactement sur un doublement, comme tout compte entier à un
facteur de `10`, on lit un seul champ et le coût est celui de **fixe**. Avec
**itération lissée**, `fbm + iter` et `fbm` sont tracés par la boucle lissée,
qui fait quelques itérations au-delà du bailout pour calculer le compte lissé.
Un facteur trop élevé fait passer les marques sous le pixel, où elles
deviennent du grain.

Avec **fixe** elle trace comme avant, au bit près. Une position enregistre les
deux sous la forme `(outfbmscale 'smooth 1.5)`, et seulement quand ils ont
bougé, si bien qu'une position qui n'y touche pas s'ouvre encore dans les
versions précédentes.

### Ce que cela coûte, et la seule chose que cela ne fera pas

Le bruit est calculé **une fois par pixel**, dans la coloration, et non une
fois par itération. Sur une image de 90000 points à une étendue de `0.02` et
400 itérations, `fbm + smooth` a pris un cinquième de plus que `smooth` seul à
quatre octaves et un quart de plus à dix-huit ; à 4000 itérations, où l'orbite
fait davantage du travail, les mêmes quatre octaves coûtaient un dixième.

Les marques sont ancrées **au plan**, non à la vue : c'est de la saleté sur la
chose, pas sur l'objectif. Donc zoomer les agrandit, et si l'on va assez loin
on finit dans une seule marque et l'image redevient propre. Augmenter la
fréquence est la réponse, et c'est pourquoi la fréquence est un réglage et non
une constante. On ne peut pas la lier plutôt au zoom : le moteur réutilise des
pixels calculés à l'échelle précédente, et une coloration qui bougerait avec la
vue leur laisserait les mauvaises marques.

Les modes sont numérotés **après** ceux en couleurs vraies, donc chaque numéro
de mode que contient déjà une position enregistrée veut toujours dire ce qu'il
voulait dire.

## Pavage

**Fractale → Coloriage intérieur → Pavage...** colore l'intérieur selon les
tuiles de l'un des quarante-cinq pavages de `randsctile`, posées sur le `z`
final, quarante par unité comme les carrés de `squares`. Son dialogue a quatre
champs.

| champ | |
| --- | --- |
| **Pavage** | quel pavage, avec le nombre sous lequel `randsctile` le connaît et son nom : **3: hexagones**, où il commence. L'onglet Pavages d'Aide → Référence des formules utilisateur a une image de chacun |
| **Coloriage** | d'où vient la couleur, partout : **angle de z**, où il commence ; **angle autour du centre de la tuile**, celui de la tuile où tombe le point, qui fait de chaque tuile une petite roue de couleur ; **zmag**, la valeur du mode de ce nom ; **aucun**, rien que les tuiles |
| **Tuiles** | comment chaque tuile se détache de ses voisines : **chacune décalée d'une quantité propre**, le long de la palette, où il commence ; **ombrées du centre au bord**, de même, et décalée d'un huitième des couleurs d'un tour de plus du centre au bord, comme un biseau ; **avec le bord marqué**, de même, et la partie extérieure de la tuile tracée dans l'autre moitié de la palette ; **une couleur unie chacune**, une couleur à elle, d'où que vienne la couleur |
| **Intensité** | à quelle distance les tuiles sont placées le long de la palette : elle multiplie chaque quantité propre d'une tuile, son décalage, le biseau de **ombrées du centre au bord** et le bord de **avec le bord marqué**. À `1`, là où elle commence, les tuiles sont comme décrites ici ; à `0` il ne reste que ce dont vient la couleur, de sorte que **angle de z** est une décomposition simple et **zmag** le zmag simple ; entre les deux les tuiles se voient à peine, et au-dessus de `1` davantage |

Le mode se tient entre les modes qui se répètent et les décompositions, et les
réglages des deux l'atteignent. **Fréquence (modes périodiques)** et **Fonction
(modes périodiques)** atteignent les tuiles comme elles atteignent tout
argument périodique de l'intérieur : à la fréquence 2 les tuiles sont deux fois
plus petites. Les [réglages des décompositions](#réglages-des-décompositions)
de l'intérieur atteignent les deux angles, qui sont découpés et parcourus comme
ceux des décompositions. **Secteurs**, **Rotation**, **Répétitions**, **Sens**
et **Vitesse** atteignent les deux ; **Angle de** n'atteint que l'angle de `z`,
celui autour du centre de la tuile n'en ayant pas l'usage ; **Anneaux par
unité** fait des anneaux du module pour l'angle de `z`, à travers la fréquence
et la fonction comme pour les décompositions, et pour l'angle autour du centre
autant d'anneaux du centre au bord, qui suivent la forme de la tuile. **zmag**,
**aucun** et une tuile unie n'en lisent aucun.

La quantité propre d'une tuile vient du hachage qui la nomme. Pour un angle
elle couvre autant de la palette qu'un tour de l'angle, si bien qu'un réglage
de la décomposition déplacé laisse les couleurs des tuiles où elles étaient ;
pour **zmag**, **aucun** et une tuile unie elle couvre la palette entière, si
bien qu'**aucun** avec **chacune décalée d'une quantité propre** revient à
**une couleur unie chacune**. **angle autour du centre de la tuile** met un
tour entier de couleurs dans chaque tuile, ce qui, à la taille qu'ont les
tuiles, fait beaucoup de couleurs pour peu de pixels : une **Vitesse**
inférieure à un, ou une fréquence inférieure à un, rend les roues lisibles.

Avec **Coloriage** sur **angle de z**, **Tuiles** sur **chacune décalée d'une
quantité propre** et les réglages des décompositions là où ils commencent, le
mode dessine ce qu'il dessinait avant d'avoir le choix, au bit près. Les
pavages qui ne se répètent jamais deviennent plats, sans quantité propre,
au-delà de quatre milliards de tuiles de l'origine. Cela vaut avec
**Intensité** à `1`. À `0` les tuiles se voient encore avec **angle autour du
centre de la tuile**, dont l'angle est mesuré depuis le centre de chaque tuile.

Les pavages périodiques coûtent à peu près ce que coûtent les autres modes
intérieurs, de 120 à 210 ns par pixel ; ceux qui ne se répètent jamais,
localisés à nouveau en chaque point, de 250 à 490. **angle autour du centre de
la tuile**, **ombrées du centre au bord** et **avec le bord marqué** mesurent
où le point se tient dans sa tuile, ce qui coûte jusqu'à deux tiers de plus sur
les pavages périodiques et peu de plus sur les autres. Une position enregistre
le mode sous la forme `(incoloring 30)`, le pavage sous la forme
`(intessellation 39)`, seulement quand ce n'est pas celui des hexagones, et le
coloriage sous la forme `(intessellstyle 'tile-angle 'edged)`, seulement quand
il a été déplacé. L'intensité est enregistrée sous la forme
`(intessellintensity 0.5)`, seulement quand elle n'est pas un.

### À l'extérieur

**Fractale → Coloriage extérieur → Pavage...** a trois modes, comme Mouvement
brownien fractionnaire : `tessellation + smooth`, les tuiles portées par le
compte lissé ; `tessellation + iter`, par le compte entier ; `tessellation`,
les tuiles seules. Les tuiles sont posées sur la place du pixel dans le plan,
non sur le point d'échappement, qui saute d'une bande à la suivante et ne
donnerait que du grain, comme lorsque le mouvement y fut essayé.

| champ | |
| --- | --- |
| **Mode** | l'un des trois |
| **Pavage**, **Coloriage**, **Tuiles**, **Intensité** | comme à l'intérieur, avec des réglages à eux. **angle de z** est l'angle du point d'échappement et **zmag** le logarithme de son module, tel que `iter+log(mag)` l'ajoute ; avec **itération lissée** l'angle est celui de quelques itérations après le bailout, là où le compte lissé est calculé |
| **Fréquence** | tuiles par unité du plan, `8` au départ : la taille des tuiles, qu'il faut augmenter en zoomant, comme la fréquence du mouvement |
| **Mise à l'échelle des tuiles** et **Facteur** | ce qui donne leur taille aux tuiles, comme [La taille des marques selon le compte](#la-taille-des-marques-selon-le-compte) la donne aux marques du mouvement : **fixe**, la fréquence seule ; **itération** ou **itération lissée**, la fréquence doublée autant de fois que le dit le facteur toutes les dix itérations de ce compte |

Entre deux doublements les tuiles sont lues aux deux tailles et les deux
couleurs sont fondues, un fondu d'une taille à la suivante : les tuiles n'ont
pas de valeur médiane vers laquelle se fondre, comme en a le mouvement. Une
palette de 256 couleurs, qui n'en a aucune entre les siennes, prend la taille
la plus proche. Entre deux doublements le pavage est lu deux fois, ce qui en
double le coût. Un facteur de `1` double la fréquence dix fois toutes les cent
itérations, ce qui amène vite les tuiles sous le pixel : quelques dixièmes sont
un meilleur départ. Les [réglages des
décompositions](#réglages-des-décompositions) de l'extérieur atteignent les
deux angles : **Itérations par bande** fait des bandes du compte pour l'angle
de z et des anneaux du centre au bord pour l'angle autour du centre.

Avec **Intensité** à `0` et **Coloriage** sur **aucun**, `tessellation + iter`
dessine comme `iter` et `tessellation` une seule couleur. Une position
enregistre le mode sous la forme `(outcoloring 21)` à `(outcoloring 23)`, et
les réglages en une seule commande,
`(outtessellation 5 'tile-angle 'shaded 0.6 40 'smooth 0.4)`, seulement quand
l'un d'eux a été déplacé.

## Plus de modes de coloration

Ce n'est pas le calcul qui décide de la couleur d'un pixel : la couleur vient
d'une formule portant sur le `z` final, le paramètre `c` et le nombre
d'itérations, et il y a une telle formule par mode. Vingt-deux de plus, dans un
sous-menu à eux sous Fractale → Coloriage intérieur et Fractale → Coloriage
extérieur, au-dessus des sous-menus en vraies couleurs.

**À l'intérieur**, où l'orbite ne s'est jamais échappée et où seuls le `z`
final et `c` disent quelque chose :

| mode | |
| --- | --- |
| `real/mag` | le cosinus de l'angle : un ombrage qui se referme sur lui-même là où `atan2` a une couture |
| `max(\|real\|,\|imag\|)` | la norme carrée, donc les bandes sont des carrés là où celles de `zmag` sont des cercles |
| `\|real\|+\|imag\|` | les mêmes bandes tournées d'un huitième de tour |
| `min(\|real\|,\|imag\|)` | proche de zéro le long de l'un ou l'autre axe, donc des rayons |
| `\|z-c\|` | à quelle distance du paramètre l'orbite s'est posée, ce qui distingue les bulbes |
| `\|z*c\|` | les deux modules l'un contre l'autre |
| `angle(z)-angle(c)` | la direction où elle s'est posée, mesurée par rapport au paramètre |
| `real*imag` | une selle, quatre lobes autour des axes |
| `sin(real)*sin(imag)` | une grille posée sur l'endroit où l'orbite s'est posée |
| `sign(imag)` | deux tons plats : une décomposition de l'intérieur, comme la décomposition binaire le fait pour l'extérieur |
| `frac(mag)` | des courbes de niveau à chaque huitième du module |
| `log(mag)` | `zmag` avec le contraste reporté sur les petites valeurs, où une orbite attirée passe son temps |
| `cos(\|z\|)` | des anneaux autour de l'endroit où l'orbite s'est posée, huit par unité du module |
| `hexagons` | trois ondes à un tiers de tour l'une de l'autre, qui s'additionnent en un réseau d'hexagones |
| `spiral` | cinq bras autour de l'endroit où l'orbite s'est posée, enroulés une fois par facteur e du module |
| `cos(\|z-c\|)` | des anneaux autour du paramètre, qui distinguent les bulbes comme le fait `\|z-c\|` |

Onze des modes intérieurs se répètent — `cos(mag)`, `mag*cos(real^2)`,
`sin(real^2-imag^2)`, `squares`, `sin(real)*sin(imag)`, `frac(mag)`, les quatre
au bas du tableau et `tessellation` — et **Fréquence (modes périodiques)**, le
premier champ de la boîte de dialogue Vitesse de coloration du côté intérieur,
dit à quelle fréquence. Elle multiplie l'argument de ce qui se répète, donc à 2
les sinus tournent deux fois plus vite, `squares` trace des carrés deux fois
plus petits et `frac(mag)` seize contours par unité du module là où il en
traçait huit. Dans `cos(mag)` elle resserre aussi les bandes qui choisissent
entre le sinus et le cosinus. À un, où elle commence, ils tracent comme ils
l'ont toujours fait, au bit près, et les autres modes ne la lisent jamais. Une
position l'enregistre comme `(incolorfrequency 2)` seulement quand elle n'est
pas un, donc une position qui n'y touche pas s'ouvre toujours dans les versions
antérieures.

**Fonction (modes périodiques)**, en dessous, courbe l'argument par l'une des
[fonctions de vitesse de coloration](#fonctions-de-vitesse-de-coloration) avant
que la fréquence le multiplie, de sorte que les bandes reviennent plus vite là
où la fonction monte raide et plus lentement là où elle s'aplatit : `square` et
`cube` les serrent vers le bord et les écartent autour de zéro, `log` et `atan`
à l'inverse, et `sin`, `cos` et `triangle` replient l'argument sur lui-même, de
sorte que les bandes font demi-tour et reviennent. Elle atteint chaque fonction
périodique des onze, dont l'argument `a` entre comme `k·f(a)`, `k` étant la
fréquence : les sinus et cosinus, le module de `frac(mag)`, les parties réelle
et imaginaire dont les bandes font les carrés de `squares` et les tuiles de
`tessellation`, et dans `cos(mag)` aussi les bandes qui choisissent entre le
sinus et le cosinus. Dans `spiral` les cinq bras restent cinq, car un nombre de
tours non entier laisserait une couture là où l'angle se referme : la fréquence
et la fonction agissent sur l'enroulement. Une fonction périodique dans
l'argument d'une autre aurait aussi le sien courbé. À `ident`, où elle
commence, rien ne change. Une position l'enregistre comme
`(incolorperiodicfun 6)`, par la place de la fonction dans la liste en comptant
à partir de zéro, seulement quand elle n'est pas ident.

**À l'extérieur**, où ce qu'il y a, c'est le compte et le point de sortie :

| mode | |
| --- | --- |
| `iter+angle` | la direction par laquelle elle est sortie, mêlée au compte : les bandes prennent une torsion |
| `iter+log(mag)` | de combien elle a dépassé le bailout, un lissage moins cher que le logarithme d'un logarithme |
| `iter+real*imag` | des biomorphes rendus continus plutôt que par seuil |
| `max(\|real\|,\|imag\|)` | par quel côté elle est sortie, sans aucun compte |
| `iter banded` | le compte replié en huit bandes, montrant la forme des ensembles de niveau plutôt que leur numéro |
| `\|real\|-\|imag\|` | à quel point le point de sortie est déséquilibré |

Les numéros des anciens modes ne bougent pas, donc une position enregistrée
continue de se charger comme avant.

La vitesse et le décalage de coloration atteignent maintenant smooth et smooth
log, ce qu'ils ne faisaient jamais : ces deux-là rendaient un pixel à eux au
lieu de passer par l'étape par laquelle finissent les autres modes, donc les
deux réglages ne faisaient rien du tout tant que l'un d'eux était choisi. Aux
valeurs avec lesquelles ils arrivent — vitesse un, décalage nul — l'étape est
l'identité, donc aucune image enregistrée auparavant ne bouge.

Smooth et smooth log suivent la forme de l'échappement. Ils interpolent entre
une passe et la suivante en demandant où l'orbite a franchi le bailout, et ils
le demandaient auparavant au module quelle que soit la forme choisie — si bien
que s'échapper sur un carré et interpoler sur un cercle mettait plus d'un
cinquième des points qui s'échappent, et la moitié avec un losange, hors de la
bande à laquelle ils appartiennent. Chaque forme transmet maintenant la
grandeur qu'elle teste réellement et le seuil auquel elle la compare. Pour le
cercle les deux sont la même chose, donc rien de ce qui était tracé avant ne
bouge.

Ils fonctionnent aussi en mode Newton, où une orbite converge au lieu de
s'échapper : la grandeur est alors de combien le dernier pas a bougé et le
seuil la limite de convergence, et l'interpolation est la même question posée à
un nombre qui descend au lieu de monter. Les deux fractales de Newton, qui
n'avaient aucune variante lissée, en ont une.

**Convergence de Newton**, dans le menu Calcul, est enfin lue. La limite était
écrite dans l'itération comme un millionième, donc le réglage était proposé,
enregistré et rechargé alors que rien ne le consultait jamais.

**Smooth et smooth log** sont passés dans ce sous-menu du menu extérieur. Ils
interpolent entre une itération et la suivante d'après le dépassement du
bailout par l'orbite, ce qui demande une formule qui s'échappe sur le bailout —
vingt des trente et une le font, et les onze autres soit convergent (les
Newton), soit s'arrêtent sur un test à elles (les Sierpinski, Koch, le Clock).
Pour celles-là le mode ne fait rien du tout, et ce n'est pas quelque chose que
la première liste devrait proposer. La formule utilisateur, qui aurait toujours
pu l'avoir et ne l'avait pas, l'a maintenant.

Ce que ces modes ne peuvent pas faire, c'est tout ce qui demande l'orbite
plutôt que sa fin : un piège d'orbite, une moyenne de rayures, un exposant de
Lyapunov, l'estimation de distance. Il faut pour cela une valeur portée d'une
itération à la suivante, ce qui dans ce moteur veut dire une seconde copie
compilée de chaque boucle de calcul. `trap` et `stripe` dans une formule
utilisateur le font à la place.

### Réglages des décompositions

**Fractale → Coloriage intérieur → Décomposition**, et de même sous Coloriage
extérieur, règle la façon dont les décompositions découpent le tour : binary
decomposition et color decomposition à l'extérieur, `sign(imag)`,
`Decomposition-like` et `fbm + decomposition` à l'intérieur, et les angles des
pavages, intérieur et extérieur. Avec chaque champ là où il commence — zéro, et
un, droite et un pour les trois derniers —, elles dessinent comme elles l'ont
toujours fait, au bit près.

| champ | |
| --- | --- |
| **Secteurs** | en combien chaque répétition du tour est découpée. Zéro est celui du mode : les deux demi-plans des binaires, le tour entier des autres. Les binaires colorent des secteurs alternés, donc un nombre impair passe au pair suivant — deux secteurs d'une même couleur se toucheraient. Celles en couleur deviennent autant de paliers plats |
| **Rotation** | de combien de degrés les secteurs sont tournés, dans le sens antihoraire |
| **Angle de** | à l'intérieur seulement : l'angle de `z`, de `z-c` ou de `z*c`. `z-c` centre les secteurs sur le paramètre, et `z*c` les tourne de l'angle du paramètre |
| **Anneaux par unité** (intérieur) | des anneaux du module de ce dont on prend l'angle, un sur deux distingué : les binaires deviennent un damier polaire. Ils passent par `k·f(a)` comme tout argument périodique de l'intérieur, donc la fréquence et la fonction des modes périodiques les atteignent |
| **Itérations par bande** (extérieur) | des bandes d'autant d'itérations à travers les secteurs : avec la binary decomposition, la grille de rayons externes et d'équipotentielles |
| **Répétitions** | combien de fois le tour se répète : les secteurs, et pour celles en couleurs la palette, recommencent dans chacune. N'importe quel nombre, mais une fraction laisse une couture là où l'angle se referme |
| **Sens** | droite, comme elles allaient toujours ; gauche, dans l'autre sens ; ambidextre, aller et retour dans chaque répétition — un miroir, sans couture quand les répétitions sont entières |
| **Vitesse** | pour celles en couleurs, quelle part de la palette couvre une répétition : au-dessus de un, plus de couleurs dans chacune, pas plus de répétitions. Pour les binaires, l'écart entre les couleurs de secteurs alternés : à un comme toujours, à zéro une seule couleur, à deux le double — la binary decomposition porte les secteurs impairs du compte vers le compte à rebours, `sign(imag)` sa deuxième teinte depuis la première |

Une position les enregistre comme
`(indecomposition 8 22.5 'z-c 6 2.5 'both 1.5)` et
`(outdecomposition 6 30 2 3 'left 2)` seulement quand ils ont changé, donc une
position qui n'y touche pas s'ouvre toujours dans 2.5. Changés, une
décomposition n'est plus symétrique par rapport à l'axe réel, et XaoS cesse de
dessiner une moitié comme le miroir de l'autre.

## Fonctions de vitesse de coloration

**Fractale → Coloriage intérieur → Vitesse de coloration**, et la même chose
sous Coloriage extérieur, est une seule boîte de dialogue pour tout ce qui mène
la valeur que transmet un mode de coloration à une place dans la palette. La
**fonction**, choisie dans une liste, courbe la valeur avant que la **vitesse**
et le **décalage** de coloration soient appliqués : la palette est lue à
`f(v) × speed + shift`. Dans une colonne à lui, à côté, se trouve le
[changement de vitesse après ITER](#changement-de-vitesse-après-iter), qui agit
sur la valeur avant que la fonction la reçoive. Celle du côté intérieur
commence par la fréquence et la fonction des modes qui se répètent, qui
agissent sur le mode lui-même ; voir [Plus de modes de
coloration](#plus-de-modes-de-coloration). OK règle le tout en une seule étape
d'annulation. Pour les modes qui comptent — iter, smooth et les autres — `v`
est le nombre d'itérations lui-même, donc une fonction peut se penser en
itérations. Toutes sont impaires, `f(−v) = −f(v)`, donc une grandeur qui change
de signe le garde.

C'étaient quatre entrées du menu : la vitesse, le décalage, un sous-menu des
fonctions et la boîte de dialogue propre au changement de vitesse. Une position
est écrite et lue exactement comme avant, et les commandes que ces entrées
lançaient sont toujours là pour un script.

`sin`, `cos` et `triangle` tournent : la vitesse de coloration, appliquée après
elles, dit combien d'entrées de la palette couvre un tour, et laisse tel quel
la fréquence à laquelle il revient. C'est le rôle de la **vitesse initiale** du
changement de vitesse, appliquée avant : avec les deux groupes éteints et
l'exposant initial à un, la fonction reçoit `s0·v`, donc 2 rapproche les tours
deux fois plus et 0,5 les écarte deux fois plus.

| | |
| --- | --- |
| `ident`, `log`, `atan`, `sin`, `cos`, `square`, `cube`, `square root`, `cube root` | ce qu'il y a toujours eu, inchangé au bit près |
| `log log` | `log(1 + log(1 + v))`, pour les zooms profonds où log seul serre encore trop les bandes : à un million d'itérations log donne 13.8 bandes et celle-ci 2.7 |
| `exp` | `e^v − 1`, l'inverse de log : elle ouvre les petites valeurs. Pour les grandeurs qui restent petites — les composantes que lisent les vraies couleurs et les décompositions ; à un compte de quatre-vingt-dix elle dépasse déjà 10³⁸ |
| `triangle` | le sinus aux courbes redressées : même période et même hauteur, la palette parcourue en avant et en arrière à pas régulier au lieu de s'attarder dans les virages |

Les nouvelles viennent après les anciennes, donc chaque numéro que contient une
position enregistrée nomme toujours ce qu'il nommait.

### Changement de vitesse après ITER

La partie inférieure de la même boîte de dialogue n'est pas une fonction mais
une modification de la fonction choisie : jusqu'à une itération que vous
indiquez, la fonction continue exactement comme elle était ; au-delà, le compte
qu'elle reçoit court plus vite ou plus lentement, donc les couleurs accélèrent
ou ralentissent à partir de là et pas avant. Avec `square root` choisie et ITER
à 40, les couleurs sont celles de la racine carrée jusqu'à 40 itérations, et
ensuite la racine carrée d'un compte qui prend de la vitesse. Elle contient un
exposant et une vitesse initiaux puis deux groupes, chacun séparé par un trait
:

| | |
| --- | --- |
| **Exposant initial** | l'exposant du début jusqu'au premier ITER d'un groupe activé : le compte devient la puissance simple `v^p0`, la même que les groupes soient activés ou non, et avec les deux désactivés il vaut pour tout compte. À un, où il commence, le compte est intact jusqu'au premier ITER |
| **Vitesse initiale** | ce par quoi est multiplié le rythme de ce tronçon : `s0·v^p0`. À un, où elle commence, elle ne change rien |
| **ITER** | là où commence le tronçon du groupe, dans les unités que reçoit la fonction. À zéro ou en dessous il commence au début — ce qui fait de `ident` une puissance d'exposant quelconque |
| **Exposant** | l'exposant du tronçon, de son ITER jusqu'au suivant ou au-delà : comment le tronçon se courbe. Il part du compte auquel était arrivé le précédent, `f(ITER)`, et ajoute `s·ITER·((v/ITER)^p − 1)`, dont le rythme à ITER est `s·p` : au-dessus de un les couleurs accélèrent à partir de là, en dessous de un elles ralentissent |
| **Vitesse** | ce par quoi est multiplié le rythme du tronçon — linéaire là où l'exposant est une puissance : à un exposant de un, le tronçon est la droite `f(ITER) + s·(v − ITER)`. À un, où elle commence, elle ne change rien |
| **Mode** | **off**, **jump** ou **smooth** : si le groupe est utilisé, et comment son tronçon se raccorde au précédent |
| **Rayon** | en smooth seulement, et affiché seulement alors : combien d'itérations prend le raccord de chaque côté d'ITER. Tout nombre à partir de zéro, qui est le jump |

Avec les deux groupes activés et le premier ITER sous le second, le compte vaut
`v^p0` jusqu'à `i1`, continue de là au rythme de `p1`, et à partir de `i2` au
rythme de `p2` — **chaque tronçon avec son propre exposant et sa propre
vitesse**, jamais `p0·p1` ni `p0·p1·p2`. Avec `p0` = 2, `p1` = 1 à 50 et `p2` =
2 à 100, tous en jump et toutes les vitesses à un, le compte sur `ident` vaut
`v²` jusqu'à 50 ; `v + 2450` de 50 à 100, un exposant de un qui continue au
rythme un depuis les 2500 auxquels était arrivé le carré ; et à partir de 100
`2550 + 100·((v/100)² − 1)`, de nouveau au carré à partir d'un rythme de 2 à
cet endroit. Une vitesse de 3 sur le premier groupe fait du tronçon du milieu
`2500 + 3·(v − 50)`. L'ITER le plus bas passe en premier quel que soit le
groupe qui le porte, et **au même ITER le second groupe remplace le premier**,
rayon, vitesse et tout ; un groupe dont l'ITER est à zéro ou en dessous met son
propre exposant et sa propre vitesse sur le premier tronçon.

**Un groupe désactivé est laissé de côté**, avec son ITER :

| premier groupe | second groupe | le compte court à |
| --- | --- | --- |
| activé | activé | `p0` jusqu'à l'ITER le plus bas, puis l'exposant de ce groupe jusqu'à l'autre ITER, puis celui de l'autre |
| désactivé | activé | `p0` jusqu'à l'ITER du second, puis `p2` |
| activé | désactivé | `p0` jusqu'à l'ITER du premier, puis `p1` |
| désactivé | désactivé | `p0` pour tout compte, la puissance simple `v^p0` |

Les deux commencent désactivés, et l'exposant initial et toutes les vitesses à
un, où l'on ne fait rien du tout à la valeur, donc chaque fonction donne au bit
près ce qu'elle donnait avant qu'une telle chose existe.

Les deux façons de raccorder, selon la manière dont elles changent le compte
transmis à la fonction (écrites pour `ident`, où le compte est ce que voit la
palette, avec `f(ITER)` le compte en arrivant à ITER — `ITER^p0` au premier —
`p` l'exposant du groupe, `s` sa vitesse et `r` son rayon) :

| mode | de ITER − r à ITER + r | au-delà |
| --- | --- | --- |
| **jump** | — | `f(ITER) + s·ITER·((v/ITER)^p − 1)` : le rythme tourne à ITER d'un coup, ce qui se voit comme un anneau là où les bandes changent de rythme |
| **smooth** | le rythme va en ligne droite de celui qu'a le compte à `ITER − r` à celui qu'il a à `ITER + r`, à mi-chemin à ITER, avec une parabole autour d'ITER ajoutée qui rend ce que la droite a coupé du coin du jump | le compte du jump, au chiffre près |

**Smooth** est symétrique autour d'ITER et change les couleurs dans le rayon et
nulle part ailleurs : avant `ITER − r` et après `ITER + r` l'image est celle du
jump. Là où les deux rythmes sont réguliers, la parabole est nulle et le coin
est coupé comme le coupe une parabole — avec le compte intact jusqu'à un ITER
de 50 et un exposant de zéro au-delà, un rayon de 10 le met à 47.5 à ITER, un
quart du rayon de moins, et de nouveau à 50 à partir de 60 ; d'un rythme de un
à un rythme de trois, c'est deux à ITER. Là où un rythme se courbe, la parabole
rend ce que la courbure ajoute : après `v²` jusqu'à 50 et un exposant de un, un
rayon de 20 mène le rythme de 60 à 30 à 1 à 70, 45.5 à ITER, et continue à
partir de 70 comme `v + 2450`, comme le fait jump.

Il mélangeait auparavant les rythmes des deux tronçons comme chacun aurait
continué : le rythme du carré, toujours croissant, continuait de croître loin
dans l'arrondi, donc il baissait à peine avant ITER et n'était à mi-descente
qu'un demi-rayon après, et ce qu'il gagnait déplaçait chaque couleur au-delà de
l'arrondi — de 133 bandes avec un rayon de 20.

La parabole ne mène jamais le rythme sous le plus lent des deux bouts. Là où
elle le ferait — un rythme lent rendu très rapide sur un rayon court, comme une
puissance centième sur dix itérations — il est tenu là, et le compte au-delà de
l'arrondi est décalé de ce que cela retient. Un rayon plus large que le premier
ITER s'arrête avant zéro, et deux qui se chevaucheraient sont réduits de même
pour que chacun reste entre ses voisins — des rayons de 40 et 40 entre des ITER
distants de 30 deviennent 15 chacun.

Un exposant de 100 est permis et inutile au-delà d'ITER : quelques itérations
plus loin, les bandes sont plus fines que les pixels et la couleur est du
bruit. Avec le second groupe à ITER 50 et un exposant de 2 et le premier à 100
et 1.5, tous deux en jump, et un exposant initial de un, le compte sur `ident`
vaut 30 à 30, `50·1.5²` = 112.5 à 75, et `200 + 100·(1.5^1.5 − 1)` = 284 à 150
: le tronçon à partir de 100 part des 200 auxquels était arrivé le précédent, à
un rythme de 1.5.

Il y avait auparavant quatre modes, et les positions enregistrées avant le
rayon les nomment — `outcolorfunafter` de la 1.8.0, et `outcolorfunchange` avec
deux groupes et sans rayon. Ils sont lus dans les deux qui existent, avec
toutes les vitesses à un :

| enregistré comme | lu comme |
| --- | --- |
| off, jump | la même chose |
| smooth, qui arrondissait le coin dans le logarithme du compte, un peu de chaque côté d'ITER si loin qu'il soit | smooth avec un rayon de trois quarts d'ITER : à moins de quatre centièmes d'ITER de ce qu'il traçait jusqu'à un exposant de un et demi, un dixième jusqu'à deux |
| continuous, qui continuait à partir d'ITER au rythme avec lequel arrivait le compte | smooth à trois quarts d'ITER, avec l'exposant `0.66 + p/4 + 0.09·p²` — un à un, où continuous ne changeait rien non plus — à moins de douze centièmes d'ITER de lui jusqu'à un demi |

Trouvé en mettant les courbes côte à côte sur quatre ITER : une telle image
sort semblable, pas identique au pixel près. Une position dont le changement de
vitesse n'est rien de plus qu'un groupe en jump, l'exposant initial et toutes
les vitesses à un et le reste comme au départ, est enregistrée comme
l'enregistrait la 1.8.0, `(outcolorfunafter 40 1.5 'jump)`, et s'ouvre dans la
1.8.0 et y trace la même chose. Tout le reste est enregistré comme
l'enregistrait la 2.0, dans l'ordre de la boîte de dialogue sans les vitesses,
`(outcolorspeedchange 1 40 1.5 'smooth 10 50 2 'off 10)` ; et les vitesses,
seulement quand l'une d'elles n'est pas un, à part comme
`(outcolorspeedlinear 1 3 1)` — la vitesse initiale et celles des deux groupes
— pour qu'une position qui ne s'en sert pas s'ouvre toujours dans la 2.0. Le
changement de vitesse seul, tel que l'avait sa propre boîte de dialogue, est la
commande `outcolorspeedgroups`, avec ses douze champs dans leur ordre.

### Au-delà de la fin d'un int

La palette est atteinte en 256es d'entrée à travers un int, et une valeur
au-delà de ce que contient un int — 8.4 millions de bandes — était convertie
malgré tout, ce qui est indéfini : chaque pixel de ce genre prenait une seule
et même entrée. `cube` y arrive à 204 itérations, donc une image colorée au
cube au-delà était une tache plate ; `exp` et un fort changement de vitesse y
arrivent tout de suite. Elle tombe maintenant là où dit sa valeur, enroulée
autour de la palette comme toute valeur plus petite ; une valeur trop grande
pour que ses 256es tiennent dans un double tombe sur la première entrée, comme
l'a toujours fait un infini. À l'intérieur d'un int rien ne bouge. Au-delà, le
seul changement à une image enregistrée auparavant est dans la coloration
lissée de Magnet et Magnet2 loin de l'ensemble, qui était une seule couleur
plate et montre maintenant les bandes que le compte y décrit.

Deux réparations plus petites sont venues avec. Choisir une fonction de vitesse
déplaçait la coche dans le menu du coloriage intérieur au lieu du sien. Et
**initstate**, par lequel commence toute position, ne remettait pas à leurs
valeurs par défaut les cinq nombres du fbm — une position ne les enregistre que
là où ils en diffèrent, donc une position chargée après une autre qui les avait
changés traçait avec ceux de l'autre ; il les réinitialise maintenant, et
désactive avec eux le changement de vitesse après ITER.

### Le compte, entier

La valeur que transmet un mode de coloration perdait des chiffres en route vers
la palette à deux endroits. La coloration lissée faisait de son compte un int
de 256es d'itération avant que la vitesse de coloration, la fonction ou le
changement de vitesse ne le reçoivent ; et ce que rendait la fonction était un
float, qui au-delà de 65 536 entrées garde moins d'un 256e d'entrée et au-delà
de seize millions pas même une entière. Ni l'un ni l'autre ne se voit à une
vitesse de un, où un 256e d'itération est un 256e d'entrée. Les deux se voient
dès que quelque chose fait courir les couleurs : avec un changement de vitesse
portant le compte de zéro à dix mille entrées par itération, une vue à trente
itérations, où le rythme était de 2500, sortait en rayures d'une seule couleur
de huit pixels de large, chacune à un 256e d'itération et dix entrées de la
suivante, là où la palette aurait dû parcourir chaque entrée intermédiaire.

Le compte continue maintenant entier, et la fonction, la vitesse de coloration
et le changement de vitesse sont calculés dans le type propre de la version,
comme les formules — long double dans la version ordinaire, `__float128` en
quad — et rendent ce qu'ils calculent. Les deux versions colorent aussi
semblablement qu'avant : sur les images essayées, elles diffèrent là où elles
ont toujours différé, dans un zoom profond dont quad calcule les orbites plus
loin, et nulle part ailleurs.

Deux choses sont gardées comme elles étaient.

* **Un infini** — le milieu d'un potentiel parti à l'infini, le fond autour
  d'une formule qui s'échappe tout de suite — tombe sur la première entrée,
  comme une valeur qui n'est pas un nombre. Mise en 256es dans un float, chaque
  valeur de ce genre sortait autrefois infinie et tombait là ; quelle couleur
  doit prendre un infini, nul ne le sait, et ainsi il prend celle qu'il a
  toujours prise. Une valeur finie au-delà de ce que contient un float tombe là
  où elle dit.
* **Une valeur sur le bord d'un 256e** y reste : la valeur est arrondie à un
  256e de 256e avant d'être coupée à un 256e, puisque la coloration ne dit rien
  de plus fin. La partie réelle d'une orbite de Newton s'arrête sur une racine,
  à un ou à moins un demi, ce qui met `iter+real` exactement sur un tel bord,
  et les derniers chiffres de la convergence mettraient chaque pixel d'un côté
  ou de l'autre — un 256e d'entrée, que personne ne voit, jusqu'à ce qu'un
  filtre de contours trace chaque pixel de ce genre comme un bord et que les
  courbes de niveau sortent dentées comme des engrenages. Le float en
  arrondissait la plupart et laissait des encoches ; il n'y en a plus aucune.

Rendus par rapport à la version précédente en 640×480 — avant le changement de
la section suivante — sur les 169 exemples, 83 se traçaient au pixel près comme
avant et 71 de plus différaient d'un seul niveau sur au plus 113 pixels. Des
quinze autres, ceux en vraies couleurs et en smooth log — `perb4`, `spiral4`,
`truecol1`, `truecol3` — diffèrent jusqu'à cinq niveaux sur jusqu'à un septième
des pixels, leurs composantes étant elles aussi revenues en floats ; `edge` a
perdu ses encoches ; et le reste diffère sur une poignée de pixels où une
valeur était sur le bord d'une entrée et où un filtre ou une palette tracée
sans fondu montre la marche. Magnet et Magnet2 en coloration lissée changent le
plus : loin de l'ensemble, leur compte est négatif, ce que l'int, lu sans
signe, avançait de quatre milliards de 256es, et que le float traçait ensuite
par marches d'une entrée ou deux — des taches plates aux bords nets là où il y
a maintenant le dégradé que décrit le compte.

Sur une image qui sort identique, cela coûte de deux à cinq pour cent dans la
version ordinaire, et de sept à vingt en quad, où les logarithmes du compte et
les fonctions sont calculés de façon logicielle. Là où les rayures deviennent
un dégradé, le rendu prend aussi plus de temps : l'estimation par blocs ne
remplit un pixel qu'à partir de voisins d'une seule couleur, et l'image écrite
est plus grande — la vue à trente itérations, en 1200×1200, a pris 1.07 seconde
là où ses rayures prenaient 0.30.

### Le rythme du compte lissé

La coloration lissée ajoute au compte des itérations entières la part de la
dernière que le point avait parcourue en s'échappant. Elle la calculait là où
le bailout tombe entre la grandeur avant la passe de sortie et après,
linéairement dans son logarithme. Le compte continue ainsi sans rupture d'une
itération à la suivante, mais pas son rythme : d'une passe à la suivante la
grandeur est presque élevée au carré, donc le rythme change le long de chaque
itération et tourne d'un coup là où commence la suivante. Dans une vue large de
2.6×10⁻⁷, où le compte va de 59.94 à 60.02, il tournait de 4.26 au bord entre
les deux. À une vitesse de coloration de un, cela déplace une couleur de
quelques centièmes d'entrée ; avec un changement de vitesse à dix mille entrées
par itération, ce sont les bandes qui changent de densité d'un coup le long de
chaque bord d'une itération.

L'orbite est maintenant poursuivie au-delà du bailout jusqu'à ce que son module
au carré dépasse un million — trois ou quatre passes de plus pour un carré — et
la fraction est prise dans le logarithme du logarithme, où une élévation au
carré est un pas de longueur un :

    n + ln(ln T / ln s) / ln(ln m / ln s)

avec `s` et `m` le module au carré de part et d'autre d'un million, `T` le
bailout et `n` les passes faites. Les passes au-delà du bailout se retirent
d'elles-mêmes, et le degré est celui de l'orbite elle-même, `ln m / ln s`, donc
un cube est compté comme un cube. Le compte reste à moins d'un dixième
d'itération de ce qu'il était, et son rythme continue : le virage de 4.26 dans
cette vue sort à 1.0000, et le long d'un rayon venant de l'extérieur de
l'ensemble les virages de 1.3 à 1.8 à chaque bord sortent à 1.0025.

| | le compte |
| --- | --- |
| les Mandelbrot, Triceratops, Mandelbar, Lambda, Manowar, Spider, Phoenix et Octo, sur le cercle | poursuivi au-delà du bailout |
| une formule utilisateur qui s'échappe sur le cercle | poursuivi tant que l'orbite croît au moins comme une puissance un et demi, et compté au bailout sinon |
| les trois Barnsley | au bailout : elles multiplient l'orbite par `c`, de sorte que son logarithme croît d'un pas au lieu de doubler, et le compte au bailout est déjà le bon pour elles |
| Newton, Magnet, Catseye | au bailout : elles convergent, ou testent leur sortie à leur façon |
| toute forme de l'échappement sauf le cercle | au bailout |

Les formes restent comme elles étaient parce que chacune mesure l'orbite par
quelque chose qui tourne avec elle, alors que chaque passe double son angle :
poursuivi avec sa propre mesure, le degré sortait assez instable pour aggraver
les virages — ceux de l'hexagone de 1.28 à entre 1.67 et 3.57 le long de ce
rayon — et ceux du carré allaient même jusqu'à revenir en arrière.

Sur les 169 exemples, les 154 qui ne colorent pas en lissé se tracent au pixel
près comme avant ceci. Les quinze qui le font changent jusqu'à treize niveaux
sur jusqu'à neuf dixièmes des pixels — `eye` jusqu'à 39, là où les entrées
voisines de sa palette sont très éloignées, et `spiral4`, dont le smooth log
sur une palette en vraies couleurs grossit tout changement du compte, au milieu
d'un petit détail. Par rapport à la 2.1, avec les changements de la section
précédente, 84 se tracent au pixel près comme avant et 59 de plus diffèrent
d'au plus trois niveaux sur au plus 500 pixels.

Cela coûte ces passes et deux logarithmes de plus pour chaque pixel qui
s'échappe : sur les images lissées essayées, de 11 à 28 pour cent dans la
version ordinaire et de 18 à 20 en quad, où un logarithme en `__float128` prend
946 ns contre 39 en long double ; une image qui ne colore pas en lissé est
inchangée. Calculés en double, les quatre logarithmes auraient coûté de 8 à 19
pour cent dans la version ordinaire et épargné 8 en quad — un double porte la
fraction à seize chiffres, là où une vitesse de dix mille en demande neuf —
mais la coloration est calculée de bout en bout dans le type propre de la
version, comme les formules.

## Zoomer

**Le zoom par sélection** (Vue → Zoomer) remplace le zoom continu par un
rectangle : tirez avec le bouton gauche et relâchez pour y aller, comme dans
Ultra Fractal. En avant seulement. Le rectangle est ramené aux proportions de
la fenêtre en allongeant son petit côté, jamais en rognant le grand, et il
fonctionne correctement avec la fractale tournée. Le mode Julia rapide prend le
tracé pendant qu'il tourne, et le rend en sortant.

**Pas fixes**, dans le même sous-menu :

| | |
| --- | --- |
| Zoomer 2× | `Ctrl` `+` |
| Dézoomer 2× | `Ctrl` `-` |
| Zoomer 10× | `Ctrl` `Alt` `+` |
| Dézoomer 10× | `Ctrl` `Alt` `-` |

Ils mettent à l'échelle autour du centre de la vue, donc ils s'annulent
exactement les uns les autres.

## Tourner et centrer la vue

**Fractale → Vue** règle le centre, le rayon et l'angle par des nombres. Sous
l'angle, elle a deux outils pour ces mêmes champs. Ils remplissent les champs
et ne font rien d'autre : OK applique ce qu'ils affichent en une seule étape
d'annulation, Annuler rien de tout cela, et une position enregistrée n'en sait
jamais rien.

**Parts** découpe un tour complet en autant de parts égales, et **−** et **+**
tournent l'angle d'une part en arrière ou en avant : 60 degrés pour six, 360
pour une. Avec aucune, les deux sont désactivés. Les tours sont comptés depuis
l'angle d'où ils sont partis, donc six sixièmes font un tour entier jusqu'au
dernier chiffre ; taper un angle, ou changer le nombre de parts, recommence le
compte à partir de là. Le nombre de parts est gardé tant que XaoS tourne, pour
la prochaine ouverture de la boîte de dialogue.

**Centrer sur l'axe** déplace le centre tout droit à travers l'image (**R**) ou
tout droit vers le haut (**I**) jusqu'à ce que l'origine soit au milieu dans ce
sens, à l'angle de la boîte de dialogue. Sans rotation, c'est la partie réelle
ou imaginaire qui devient nulle et l'autre qui est gardée. Avec rotation, c'est
la même chose le long des axes propres de l'image : R à 90 degrés efface la
partie imaginaire, et R puis I amènent l'origine au milieu à n'importe quel
angle. Un centre qui y est déjà est laissé tel quel, et une partie que la
rotation laisse comme arrondi est prise pour zéro, donc le champ affiche 0 et
non 1e-17.

La rotation elle-même est calculée dans le type propre de la version, pi
compris. Elle était calculée à partir d'un pi double, ce qui dans la version
quad voulait dire aussi un sinus et un cosinus double : une paire qui n'est une
rotation qu'à seize chiffres près, donc le centre d'une vue tournée était tracé
à 5e-17 de sa distance à l'origine à côté de là où il était. En dessous de
1e-16, là où quad sert, une vue tournée de 30 degrés montrait un tout autre
endroit. Un quart de tour est maintenant exact. Dans la version ordinaire, une
image tournée est la même sauf des pixels à la limite de sa précision : 52 des
64 exemples ayant un angle sont inchangés, les autres déplacent des pixels
épars là où le détail est plus fin que l'arithmétique.

## La même image avec n'importe quel nombre de threads

Une image sortait de temps en temps différente de quelques pixels quand elle
était tracée avec plusieurs threads : `examples/Hubicka/julia3.xpf` rendu avec
quatre threads sur une machine chargée différait dans 26 rendus sur 3000, de
jusqu'à 37 pixels le long de ses filaments, et jamais avec un seul thread.

Les pixels venaient de l'estimation par blocs. Elle remplit un pixel à partir
des lignes et des colonnes déjà calculées autour quand elles concordent, et en
mode Mandelbrot elle se sert des mêmes voisins pour décider s'il vaut la peine
de tester la périodicité d'un pixel. L'estimation procède par passes — d'abord
une ligne sur six, puis les lignes intermédiaires — et les passes tournaient
dans tous les threads à la fois sans rien entre elles : un thread passait à la
passe suivante pendant que les autres en étaient encore à la précédente, et une
ligne était marquée faite avant d'être calculée. Donc qu'un voisin compte comme
calculé dépendait de l'avancement des autres threads.

Chaque passe a maintenant deux phases, d'abord ses lignes puis ses colonnes, et
chaque thread finit l'une avant qu'aucun ne commence la suivante ; ce qu'une
phase calcule ne sert pas de voisin en son sein. Les voisins sont alors les
mêmes quel que soit le nombre de threads, et les 3000 rendus sortent tous
identiques. Avec un seul thread, l'ordre est celui de toujours.

Les filtres sortent identiques aussi. **3D** trace une bande de colonnes par
thread, et l'ombrage d'une colonne est reporté depuis la colonne deux crans à
sa droite dans la ligne du dessous, donc les deux dernières colonnes d'une
bande lisent les deux premières de la bande à sa droite ; chaque bande y lisait
auparavant des copies de sa propre dernière colonne, et l'image changeait le
long des bords des bandes avec le nombre de threads. Maintenant une bande
attend, ligne par ligne, ce que la bande à sa droite y a tracé, et les bandes
tournent toujours côte à côte. **Champ d'étoiles** tirait ses étoiles d'un seul
générateur aléatoire que tous les threads faisaient avancer à la fois ;
maintenant chaque ligne part de l'endroit de la séquence qu'aurait atteint un
seul thread traçant toute l'image, donc les étoiles sont les mêmes avec
n'importe quel nombre de threads. Elles changent toujours d'un lancement de
XaoS au suivant, comme elles l'ont toujours fait : le générateur prend sa
valeur initiale à l'horloge.

**L'anticrénelage** n'a jamais été en cause, bien que trois exemples tracés
avec lui — `circle`, `helloween` et `warriormask` — sortent différents d'un
rendu à l'autre. Ce qui différait, c'était le message portant le nom de la
formule : il disparaissait après cinq secondes de l'horloge murale, et un rendu
ralenti par une machine chargée le traçait sur moins d'images, ou sur aucune.
Dans un rendu, un message dure maintenant cinq secondes du temps propre de
l'animation, quelle que soit la machine ; dans la fenêtre rien ne change.

Sur les 169 exemples rendus avec un seul thread, 167 se tracent au pixel près
comme la 1.8.0 les traçait : `galaxy` tient ses étoiles de l'horloge, et
`heart` garde maintenant un message que la 1.8.0, le chronométrant à l'horloge
murale, avait retiré avant que l'image soit tracée. Rendus huit à la fois sur
quatre cœurs, une fois avec un thread et trois fois avec quatre, tous se
tracent de même sauf `galaxy`.

## Proportions

**Vue → Proportions** redimensionne la fenêtre une fois, pour que l'image
prenne la forme choisie : 1:1, 4:3, 3:2, 16:10, 16:9, 2:1, 21:9, ou 9:16
debout. La fenêtre n'y est pas tenue : on peut ensuite la tirer à n'importe
quelle forme, et la coche du menu est sur la forme qu'a l'image à ce moment, ou
sur aucune.

L'image garde le nombre de pixels qu'elle avait et ne change que de forme — ou,
quand elle ne tiendrait pas à l'écran, devient la plus grande de cette forme
qui tient. Depuis le plein écran ou une fenêtre agrandie, elle revient d'abord
à une fenêtre ordinaire, toujours avec les pixels qu'elle avait, donc 16:9 y
est la fenêtre 16:9 la plus large que contient l'écran. La forme est respectée
au pixel près, sur un écran mis à l'échelle aussi : à 125 %, un nombre de
pixels sur cinq ne peut pas être obtenu, et la taille est choisie parmi ceux
qui le peuvent.

Ce à quoi cela sert, c'est le rendu. Une position est tracée de sorte que toute
la vue qu'elle contient soit dans l'image, et une image d'une autre forme
montre davantage du plan sur son côté le plus long. Donc une vue cadrée dans la
fenêtre ordinaire 4:3 et rendue en 1920×1080 a un tiers de plus sur les côtés
que ce qui a jamais été à l'écran ; cadrée dans une fenêtre 16:9, le rendu est
ce qu'on a vu.

## Anticrénelage en lumière linéaire

Filtres → Anticrénelage en lumière linéaire, ou `-linearaa`. L'anticrénelage
ordinaire fait la moyenne des couleurs des échantillons telles qu'elles sont
stockées, ce qui assombrit les bords parce que les valeurs stockées ne sont pas
proportionnelles à la lumière. Celui-ci convertit en lumière linéaire, y fait
la moyenne, et reconvertit, donc un bord garde la luminosité de ce entre quoi
il se trouve.

## Langues

Vue → Choisir la langue fait passer XaoS dans une autre langue, après un
redémarrage. Six sont complètes : allemand, espagnol, français, italien,
portugais et russe — les menus, les boîtes de dialogue et les messages, la
référence des formules utilisateur, et les boîtes OK, Annuler, de fichiers et
de couleurs qui viennent de Qt lui-même. Ces dernières sont maintenant
intégrées, alors qu'on les cherchait auparavant là où Qt est installé, et dans
la langue du système quelle que soit celle choisie. Les autres langues dont
XaoS a hérité vont aussi loin que leurs traducteurs les ont menées, et sont en
anglais au-delà.

Dans la référence des formules utilisateur, les descriptions sont traduites et
les noms et les listes d'arguments ne le sont pas : ils sont la syntaxe dans
laquelle s'écrit une formule.

Ce guide est traduit dans les mêmes six langues, à côté du guide anglais
`ms-xaos-guide.md`, sous les noms `ms-xaos-guide.de.md`, `.es.md`, `.fr.md`,
`.it.md`, `.pt.md` et `.ru.md`, et Aide ouvre celui de la langue dans laquelle
est XaoS. Chacun indique la version du guide anglais à partir de laquelle il a
été fait ; là où ils diffèrent, c'est l'anglais qui fait foi.

## Compilation

    cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build

Les deux exécutables arrivent dans `bin/`. Options :

| | |
| --- | --- |
| `-DQUAD=OFF` | ne pas faire le second exécutable |
| `-DDEEPZOOM=ON` | rendre quad l'exécutable principal, et n'en compiler aucun second |
| `-DDIAG=ON` | compiler aussi `XaoS-diag`, qui consigne dans `xaos-diag.txt` l'ordre de superposition et de calcul |

`cmake --build build --target deploy` remplit `bin/` avec l'environnement
d'exécution de Qt et du compilateur pour qu'il tourne sur une machine sans Qt.

## Tests

    ctest --test-dir build

351 tests : l'analyseur de formules dans les deux précisions, l'exactitude de
`gamma`, `lambertw`, `erf` et `randsc` par rapport à des références exactes là
où il en existe, les boucles d'itération par rapport à des sommes de contrôle
enregistrées, les fonctions de vitesse de coloration par rapport à ce qu'elles
donnaient avant et à ce à quoi elles servent, les positions enregistrées en
aller-retour et chaque rendu qui se termine proprement, l'enregistrement et la
restauration de la superposition, l'historique d'annulation, la rotation et le
centrage de la boîte de dialogue Vue, la référence de l'aide par rapport à la
table propre de l'analyseur, et le `poly` aléatoire de son onglet Utilitaires
par rapport au polynôme écrit terme à terme.

L'un d'eux mesure une vitesse plutôt que des réponses. Chronométrer en
nanosecondes affirmerait sur quelle machine tourne la batterie, donc il affirme
plutôt un rapport : ce qu'une dégradation coûte aux fonctions de bruit, compté
en multiplications. Cela devrait être deux, une par composante, parce que la
taille est multipliée par la dégradation une fois par passe ; l'élever plutôt à
la puissance de la passe donne six ou sept et croît avec la limite
d'itérations. Le test échoue sur les deux implémentations que celle-ci a
remplacées.
