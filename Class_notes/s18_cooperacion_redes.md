# Sesión 18: Cooperación en redes. Notas formales

> **Cómo leer estas notas.** Siguen el orden del programa: primero se observa
> en el cuaderno (`Code/s18_cooperacion_redes.jl`), después se formaliza aquí.
> Lo que no tiene marca es **ruta básica**. Las secciones y demostraciones
> marcadas **[Avanzado]** son opcionales.
>
> Repaso recomendado: [L.B4](repaso/algebra_lineal.md#lb4) (matriz de adyacencia
> y grado) y [C.B4](repaso/calculo.md#cb4) (función logística).
>
> Lecturas: Nowak & May (1992). Complementaria: Slanina (2014), sec. 8.1. Ruta
> avanzada: Santos & Pacheco (2005); Ohtsuki et al. (2006); de consulta:
> Schecter & Gintis (2016).

---

## 1. El dilema del prisionero y su extinción en población bien mezclada

En el **dilema del prisionero**, cada agente elige entre **cooperar** ($C$,
$s_i = 1$) o **traicionar** ($D$, $s_i = 0$). Dos cooperadores reciben la
recompensa $R$; dos traidores reciben el castigo $P$; y cuando un traidor se
enfrenta a un cooperador, el traidor obtiene la tentación $T$ mientras el
cooperador recibe el pago del incauto $S$, con

$$
T > R > P \ge S .
$$

En una población bien mezclada con fracción $x$ de cooperadores ([sesión 17](s17_replicador_goodwin.md)),
los pagos esperados son $\pi_C(x) = Rx + S(1-x)$ y $\pi_D(x) = Tx + P(1-x)$.
Como $T > R$ y $P \ge S$ (con al menos una desigualdad estricta), se cumple
$\pi_C(x) - \pi_D(x) < 0$ para todo $x \in (0, 1]$, por lo que en la dinámica
del replicador $\dot x = x(1-x)(\pi_C - \pi_D) < 0$ para todo $x \in (0, 1)$ y
el único punto fijo estable es la extinción total de la cooperación ($x^\ast =
0$).

---

## 2. El juego espacial de Nowak y May: reciprocidad en retículas

Nowak y May (1992) estudiaron la versión débil de un solo parámetro ($R = 1$,
$T = b > 1$, $P = S = 0$):

$$
A_{\text{NM}} = \begin{pmatrix} 1 & 0 \cr b & 0 \end{pmatrix}, \qquad b > 1 .
$$

En una red $G = (V, E)$ con matriz de adyacencia $A_{ij}$
([L.B4](repaso/algebra_lineal.md#lb4)), en cada ronda:

1. Cada nodo $i$ juega con todos sus vecinos $j \in \mathcal N_i$ y acumula
   el pago total (en la convención original de Nowak y May también juega
   consigo mismo, lo que desplaza los umbrales de $b$ pero no cambia la
   formación de clústeres; el cuaderno usa sólo a los vecinos):
   $$
   \Pi_i = \sum_{j \in \mathcal N_i} \bigl[s_i s_j + b(1 - s_i)s_j\bigr]
   = \bigl[s_i + b(1 - s_i)\bigr] k_i^C ,
   $$
   donde $k_i^C = \sum_{j \in \mathcal N_i} s_j$ es el número de vecinos
   cooperadores de $i$.
2. **Actualización por imitación:** Cada nodo compara su pago con el de sus
   vecinos:
   - *Imitación determinista del mejor* (Nowak & May, 1992; la retícula del
     cuaderno): $i$ adopta la estrategia del nodo con mayor pago en
     $\lbrace i\rbrace \cup \mathcal N_i$.
   - *Imitación proporcional* (Santos & Pacheco, 2005; la red BA del
     cuaderno): $i$ elige un vecino $j$ al azar y, si $\Pi_j > \Pi_i$, copia su
     estrategia con probabilidad $(\Pi_j - \Pi_i)/(b \max(k_i, k_j))$.
   - *Imitación de Fermi (logit por pares, Szabó & Tőke, 1998):* $i$ elige un
     vecino $j$ al azar y copia su estrategia con probabilidad logística
     ([C.B4](repaso/calculo.md#cb4)):
     $$
     \mathbb P(s_i \leftarrow s_j) = \sigma\bigl(\beta(\Pi_j - \Pi_i)\bigr)
     = \frac{1}{1 + e^{-\beta(\Pi_j - \Pi_i)}} .
     $$

**Estabilidad de clústeres de cooperadores en una cuadrícula.** ¿Por qué no se
extinguen los cooperadores en la retícula cuando $b > 1$? En una población bien
mezclada, un cooperador y un traidor enfrentan la misma fracción $x$ de
cooperadores, por lo que $\Pi_D = b x > x = \Pi_C$. En cambio, en una retícula
2D (por ejemplo, con $k = 4$ vecinos en cruz), consideremos un bloque compacto
y grande de cooperadores rodeado de traidores:

- Un cooperador en la frontera plana del bloque tiene $3$ vecinos cooperadores y
  $1$ vecino traidor, obteniendo un pago $\Pi_C = 3$; un cooperador del
  interior obtiene $4$.
- Un traidor adyacente a esa frontera plana sólo toca a $1$ cooperador (y a $3$
  traidores), por lo que su pago es $\Pi_D = b$.

Bajo imitación del mejor vecino, el traidor de la frontera se convierte en
cooperador cuando $3 > b$, así que el bloque **se expande** por sus caras planas
si $b < 3$; el cooperador de la frontera ve a un vecino interior con pago $4$,
así que el bloque **resiste** mientras $b < 4$. En vecindad de Moore ($8$
vecinos) y un bloque de $3 \times 3$, el cooperador del centro obtiene $8$, el
del borde plano $5$, y el traidor pegado a ese borde toca a $3$ cooperadores
($3b$): el bloque crece por sus caras si $5 > 3b$, es decir, $b < 5/3$, y
resiste mientras $3b < 8$. Como los pagos de los cooperadores son enteros y los
de los traidores son $b$ por un entero, la fracción estacionaria de
cooperadores sólo puede cambiar por escalones en cocientes $b = m/n$ de números
de vecinos, donde distintas geometrías de clústeres ganan o pierden
estabilidad (**reciprocidad espacial**).

---

## 3. Cooperación en redes libres de escala

En una red libre de escala de Barabási-Albert ([sesión 10](s10_redes_modelos.md)),
la distribución de grados $P(k) \sim k^{-3}$ produce concentradores (*hubs*) con
decenas de conexiones. Santos y Pacheco (2005) descubrieron que esta
heterogeneidad de grados promueve la cooperación de forma aún más potente que
una retícula regular, mediante un bucle de retroalimentación asimétrico:

- **Si un *hub* de grado $k_H \gg 1$ es cooperador** y una fracción $q$ de sus
  vecinos coopera, el *hub* acumula un pago $\Pi_H^C = q k_H$, mucho mayor que
  el pago máximo $\le b k_{\text{hoja}}$ de cualquiera de sus vecinos de bajo
  grado. Los vecinos imitan al *hub*, aumentando $q$ y reforzando aún más el
  pago del *hub*.
- **Si el *hub* es traidor**, obtiene inicialmente un pago alto $\Pi_H^D = b q
  k_H$. Precisamente por su éxito, sus vecinos copian su estrategia traidora;
  pero a medida que los vecinos del *hub* se vuelven traidores, $q \to 0$ y el
  pago del *hub* traidor se desploma a $\Pi_H^D \to 0$, haciéndolo vulnerable a
  ser reconvertido por cualquier vecino cooperador conectado a otro grupo.

Así, un *hub* traidor destruye la fuente de su propio pago, mientras que un
*hub* cooperador estabiliza a su vecindario.

---

## 4. Selección débil y la regla $b/c > k$ de reciprocidad de red [Avanzado]

Consideremos el **juego de donación** (caso particular del dilema del
prisionero): un cooperador paga un costo $c > 0$ por cada vecino para otorgarle
un beneficio $b > c$, mientras que un traidor no paga nada ni otorga beneficios.
La matriz de pagos por interacción es

$$
\begin{array}{c|cc}
 & C & D \cr
\hline
C & b - c & -c \cr
D & b & 0
\end{array}
$$

En una red regular de grado $k$, un agente con estrategia $s_i \in \lbrace 0,
1\rbrace$ y fracción local $q_i = k_i^C / k$ de vecinos cooperadores recibe un
pago medio por vecino $\pi_i = b q_i - c s_i$. Su aptitud reproductiva
(*fitness*) bajo **selección débil** ($w \ll 1$) es $f_i = 1 + w\pi_i$.

Con la regla de actualización **muerte-nacimiento** (*Death-Birth*): en cada
pasos muere un individuo $i$ al azar y sus $k$ vecinos compiten por ocupar la
vacante con probabilidad proporcional a su aptitud $f_j$.

**Resultado** (Ohtsuki, Hauert, Lieberman & Nowak, 2006). En una red grande
de grado $k$ bajo actualización muerte-nacimiento y selección débil
($w \to 0$, $N \gg k$), la probabilidad de fijación $\rho_C$ de un único
mutante cooperador supera a la de una mutación neutra ($\rho_C > 1/N$) si y
sólo si

$$
\frac{b}{c} > k .
$$

Ohtsuki et al. lo obtuvieron con la aproximación de pares y lo confirmaron por
simulación. Para grafos regulares finitos el umbral exacto bajo selección
débil es

$$
\left(\frac{b}{c}\right)^{\ast} = \frac{N - 2}{N/k - 2} ,
$$

que tiende a $k$ cuando $N \to \infty$ (Taylor, Day & Wild, 2007, para grafos
homogéneos; Chen, 2013, con una demostración rigurosa para grafos regulares).
Con $N = 36$, como en el cuaderno, el umbral es $34/16 \approx 2.1$ en el anillo
($k = 2$) y $34/7 \approx 4.9$ con $k = 4$.

*Idea del mecanismo (correlaciones locales).* Para que un cooperador tenga mayor
aptitud esperada que un traidor se requiere que los cooperadores tengan, en
promedio sobre la red, más vecinos cooperadores que los traidores:

$$
\mathbb E[\pi_C] - \mathbb E[\pi_D]
= b\bigl(q_{C\mid C} - q_{C\mid D}\bigr) - c > 0 ,
$$

donde $q_{C\mid C}$ es la probabilidad condicional de que un vecino de un
cooperador sea cooperador, y $q_{C\mid D}$ la de que un vecino de un traidor sea
cooperador. En el límite de selección débil ($w \to 0$), las correlaciones
espaciales están gobernadas por el modelo de votante neutro ([sesión 15](s15_votante_kirman.md)):
en la aproximación de pares, $q_{C\mid C} - q_{C\mid D} = 1/(k - 1)$, y al
tener en cuenta la competencia entre los vecinos de la vacante bajo la regla
muerte-nacimiento se obtiene la condición $b/c > k$ (la demostración completa
está en el suplemento de Ohtsuki et al., 2006).

Cuanto más escasas son las conexiones locales ($k$ pequeño), más fuertes son los
lazos de agrupamiento local y más fácil resulta que emerja la cooperación.
Cuando $k \to \infty$ (población bien mezclada), el umbral $b/c > k$ se vuelve
inalcanzable y la cooperación desaparece.

---

## 5. Resumen formal

| Modelo / Concepto | Resultado clave |
|---|---|
| Dilema del prisionero bien mezclado | $\pi_C(x) < \pi_D(x)$ para todo $x$; la cooperación se extingue ($x^\ast = 0$) |
| Juego espacial de Nowak y May | Pago $\Pi_i = [s_i + b(1-s_i)]k_i^C$; en retícula 2D sobreviven clústeres de cooperadores con escalones en $b = m/n$ |
| Redes libres de escala (Barabási-Albert) | Los *hubs* cooperadores estabilizan su vecindario; los *hubs* traidores destruyen su propio pago |
| Regla de reciprocidad de red | Bajo selección débil (muerte-nacimiento), la cooperación es favorecida si $b/c > k$ ($N$ grande); en un grafo regular finito, si $b/c > (N - 2)/(N/k - 2)$ |

---

## 6. Ejercicios formales

**Ruta básica**

- **F1.** En una cuadrícula 2D con vecindad de von Neumann ($k = 4$ vecinos:
  arriba, abajo, izquierda, derecha) bajo el dilema del prisionero de Nowak y
  May ($R = 1, T = b > 1, P = S = 0$, jugando sólo con los $4$ vecinos),
  calcula el pago de un cooperador aislado rodeado de traidores y el de los
  traidores vecinos, y explica por qué el cooperador aislado desaparece en una
  ronda para todo $b > 1$.
- **F2.** En la misma cuadrícula con $k = 4$ vecinos, considera un cuadrado de
  $2 \times 2$ cooperadores rodeado de traidores. Calcula el pago de cada
  cooperador del cuadrado y el pago de un traidor adyacente, y encuentra para
  qué valores de $b > 1$ sobrevive el cuadrado bajo imitación del mejor vecino.
- **F3.** En un grafo regular donde cada agente tiene $k = 4$ vecinos y juega el
  juego de donación con costo $c = 1$, ¿cuál es el beneficio mínimo $b$ necesario
  para que la reciprocidad de red favorezca la cooperación según la regla de
  Ohtsuki et al. ($N$ grande)? ¿Y si la red es un anillo ($k = 2$)?

**Ruta avanzada**

- **F4.** En una cuadrícula con vecindad de Moore ($8$ vecinos), considera un
  bloque de $3 \times 3$ cooperadores rodeado de traidores. Calcula los pagos
  del cooperador central, de un cooperador en la cara plana, de un cooperador en
  la esquina y de un traidor pegado a la cara plana, y deduce los umbrales de
  $b$ donde el bloque se expande ($b < 5/3$) y donde resiste ($b < 8/3$).
- **F5.** Compara el umbral de tamaño finito $(N-2)/(N/k - 2)$ con el límite
  $k$ para $N = 36$ y $N = 360$ cuando $k = 4$. ¿A qué tamaño $N$ difiere el
  umbral finito en menos de $5\%$ de $k = 4$?

---

## Referencias

- Axelrod, R. (1984). *The Evolution of Cooperation*. Basic Books.
- Chen, Y.-T. (2013). Sharp benefit-to-cost rules for the evolution of
  cooperation on regular graphs. *The Annals of Applied Probability*, 23(2),
  637–664. https://doi.org/10.1214/12-AAP849
- Nowak, M. A., & May, R. M. (1992). Evolutionary games and spatial chaos.
  *Nature*, 359(6398), 826–829. https://doi.org/10.1038/359826a0
- Ohtsuki, H., Hauert, C., Lieberman, E., & Nowak, M. A. (2006). A simple rule
  for the evolution of cooperation on graphs and social networks. *Nature*,
  441(7092), 502–505. https://doi.org/10.1038/nature04605
- Santos, F. C., & Pacheco, J. M. (2005). Scale-free networks provide a
  unifying framework for the emergence of cooperation. *Physical Review
  Letters*, 95(9), 098104. https://doi.org/10.1103/PhysRevLett.95.098104
- Schecter, S., & Gintis, H. (2016). *Game Theory in Action: An Introduction to
  Classical and Evolutionary Game Theory*. Princeton University Press.
- Slanina, F. (2014). *Essentials of Econophysics Modelling*. Oxford
  University Press.
- Szabó, G., & Tőke, C. (1998). Evolutionary prisoner's dilemma game on a
  square lattice. *Physical Review E*, 58(1), 69–73.
  https://doi.org/10.1103/PhysRevE.58.69
- Taylor, P. D., Day, T., & Wild, G. (2007). Evolution of cooperation in a
  finite homogeneous graph. *Nature*, 447(7143), 469–472.
  https://doi.org/10.1038/nature05784
