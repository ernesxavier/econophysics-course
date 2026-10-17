# Sesiones 10 y 11: Redes: estructura y modelos. Notas formales

> **Cómo leer estas notas.** Siguen el orden del programa: primero se observa
> en el cuaderno (`Code/s10_redes_modelos.jl`), después se formaliza aquí. Lo que no
> tiene marca es **ruta básica**. Las secciones y demostraciones marcadas
> **[Avanzado]** son opcionales.
>
> Repaso recomendado: [L.B4](repaso/algebra_lineal.md#lb4) (matriz de
> adyacencia), [P.B1](repaso/probabilidad.md#pb1) (distribuciones),
> [F.B4](repaso/fisica_estadistica.md#fb4) (transiciones de fase),
> [P.I5](repaso/probabilidad.md#pi5) (leyes de potencias),
> [L.I4](repaso/algebra_lineal.md#li4) (potencias de matrices y
> centralidad); para la ruta avanzada, [L.A1](repaso/algebra_lineal.md#la1)
> (Laplaciano), [L.A5](repaso/algebra_lineal.md#la5) (Perron-Frobenius) y
> [F.A5](repaso/fisica_estadistica.md#fa5) (percolación).
>
> Lecturas: Newman (2018), caps. 6–7 y 10–11; Jackson (2008), cap. 2. Ruta
> avanzada: Newman (2018), caps. 12, 13 y 15; Jackson (2008), caps. 3–5.

---

## 1. Representación

### 1.1 Matriz de adyacencia y grados

Una **red** (o grafo) tiene $n$ nodos y un conjunto de enlaces entre pares
de nodos. En estas notas las redes son **simples**: no dirigidas, sin pesos,
sin enlaces de un nodo consigo mismo y con a lo más un enlace por par. Se
representan con la matriz de adyacencia $A$, simétrica, con $A_{ij} = 1$ si
$i$ y $j$ están conectados y $A_{ii} = 0$.

- **Grado**: $k_i = \sum_j A_{ij}$, o $k = A\mathbf 1$.
- **Número de enlaces**: $m = \tfrac12 \sum_i k_i$, porque cada enlace se
  cuenta desde sus dos extremos. En consecuencia, $\langle k\rangle = 2m/n$.
- **Densidad**: $m / \binom{n}{2}$, la fracción de pares conectados.
- **Distribución de grados**: $p_k$, la fracción de nodos con grado $k$. Su
  CCDF es $P(K \ge k) = \sum_{j \ge k} p_j$.

### 1.2 Caminos, distancias y componentes

Un **camino** de longitud $\ell$ es una sucesión de $\ell$ enlaces
consecutivos. La **distancia** $d_{ij}$ es la longitud del camino más corto
entre $i$ y $j$ (infinita si no hay camino). Una **componente conexa** es
un conjunto máximo de nodos con caminos entre todos sus pares.

**Proposición (caminos).** $[A^\ell]_{ij}$ es el número de caminos de
longitud $\ell$ de $i$ a $j$ (los caminos pueden repetir nodos).

*Demostración.* Por inducción sobre $\ell$. La fórmula del producto de
matrices da

$$
\left[A^{\ell+1}\right]_{ij} = \sum_r \left[A^\ell\right]_{ir} A_{rj} ,
$$

que suma, sobre cada nodo $r$ vecino de $j$, los caminos de longitud $\ell$
que llegan de $i$ a $r$. $\square$

**Corolario (triángulos).** Un triángulo es un camino cerrado de longitud 3.
Cada triángulo da 6 caminos cerrados de longitud 3: uno por cada vértice de
salida y cada sentido. Por eso

$$
\text{triángulos} = \frac{1}{6}\thinspace\mathrm{tr}(A^3),
\qquad
t_i = \frac{1}{2}[A^3]_{ii}
$$

donde $t_i$ es el número de triángulos que pasan por el nodo $i$.

Ejemplo en Julia, con la red de cuatro nodos del repaso L.B4:

```julia
using LinearAlgebra
A = [0 1 1 0; 1 0 1 0; 1 1 0 1; 0 0 1 0]
k = A * ones(Int, 4)      # grados: [2, 2, 3, 1]
tr(A^3) ÷ 6               # 1 triángulo
```

---

## 2. Medidas

### 2.1 Centralidades

| Centralidad | Definición | Idea |
|---|---|---|
| Grado | $k_i$ | Contactos directos |
| Cercanía | $c_i = (n-1) / \sum_{j \ne i} d_{ij}$ | Inverso de la distancia promedio a los demás |
| Intermediación | $b_i = \sum_{s < t;\ s, t \ne i} \sigma_{st}(i)/\sigma_{st}$ | Fracción de caminos más cortos que pasan por $i$ |
| Eigenvector | $Ax = \lambda_1 x$ | Conectado con nodos importantes |

En la intermediación, $\sigma_{st}$ es el número de caminos más cortos entre
$s$ y $t$, y $\sigma_{st}(i)$ el número de ellos que pasan por $i$. Se suele
normalizar dividiendo entre el número de pares $\binom{n-1}{2}$, de modo que
$0 \le b_i \le 1$; así lo hace `betweenness_centrality` de `Graphs.jl`.

**Centralidad de eigenvector.** La idea "un nodo es importante si sus
vecinos lo son" se escribe $x_i = \lambda^{-1}\sum_j A_{ij}x_j$, es decir,
$Ax = \lambda x$. Hace falta que $x$ tenga todas sus componentes positivas
para interpretarlo como importancia.

**[Avanzado] Proposición.** Si la red es conexa, el eigenvector del mayor
eigenvalor $\lambda_1$ de $A$ tiene todas sus componentes del mismo signo
(y se puede tomar positivo), y es el único con esa propiedad.

*Demostración.* Una red conexa tiene $A$ no negativa e irreducible; el
teorema de Perron-Frobenius ([L.A5](repaso/algebra_lineal.md#la5)) da un
eigenvalor real simple $\lambda_1 = r(A)$ con eigenvector positivo. Como
$A$ es simétrica, los eigenvectores de eigenvalores distintos son
ortogonales: ninguno puede ser ortogonal a un vector positivo y tener
todas sus componentes positivas. $\square$

Si la red no es conexa, la centralidad de eigenvector es cero fuera de la
componente con mayor $\lambda_1$: conviene calcularla por componente.

### 2.2 Agrupamiento

El **coeficiente de agrupamiento local** del nodo $i$, con $k_i \ge 2$, es
la fracción de pares de vecinos de $i$ que están conectados entre sí:

$$
C_i = \frac{t_i}{\binom{k_i}{2}} = \frac{[A^3]_{ii}}{k_i(k_i - 1)} .
$$

Para $k_i < 2$ se define $C_i = 0$. El **agrupamiento global**
(transitividad) cuenta triángulos sobre tripletas conectadas, es decir,
caminos de longitud 2:

$$
C = \frac{3 \times \text{triángulos}}{\text{tripletas conectadas}}
= \frac{\mathrm{tr}(A^3)}{\sum_i k_i(k_i - 1)} .
$$

Es la probabilidad de que dos vecinos de un nodo estén conectados. Es la
que calcula `global_clustering_coefficient`. No es igual al promedio de los
$C_i$: el global da más peso a los nodos de grado alto.

### 2.3 Las familias de Florencia

En la red de matrimonios de Padgett y Ansell (1993), sin los Pucci (que no
tienen enlaces), hay $n = 15$ familias, $m = 20$ enlaces y 3 triángulos.

| Familia | Grado | Cercanía | Intermediación | Eigenvector |
|---|---|---|---|---|
| Médici | 6 | 0.56 | 0.52 | 0.43 |
| Guadagni | 4 | 0.47 | 0.25 | 0.29 |
| Albizzi | 3 | 0.48 | 0.21 | 0.24 |
| Strozzi | 4 | 0.44 | 0.10 | 0.36 |

Las cuatro centralidades ponen a los Médici en primer lugar, pero sólo la
intermediación marca una diferencia de 2 a 1 con el segundo lugar. Los
Strozzi tienen tanto grado como los Guadagni y una centralidad de
eigenvector alta (están conectados con familias conectadas), pero casi no
son intermediarios: sus vecinos se conocen entre sí.

---

## 3. Erdős-Rényi

### 3.1 Definición y grados

**Modelo** $G(n, p)$: cada uno de los $\binom{n}{2}$ pares se conecta con
probabilidad $p$, independientemente. El grado de un nodo es una suma de
$n - 1$ variables de Bernoulli independientes:

$$
K \sim \mathrm{Binomial}(n-1, p),
\qquad
\langle k\rangle = (n-1)p \equiv c .
$$

**Proposición.** Si $n \to \infty$ con $c$ fijo, $p_k \to e^{-c} c^k / k!$
(Poisson).

*Demostración.* Con $p = c/(n-1)$,

$$
\binom{n-1}{k}p^k(1-p)^{n-1-k}
= \frac{(n-1)(n-2)\cdots(n-k)}{(n-1)^k}\thinspace\frac{c^k}{k!}
\left(1 - \frac{c}{n-1}\right)^{n-1-k} .
$$

El primer cociente tiende a 1 porque tiene $k$ factores, cada uno con
límite 1, y el último factor tiende a $e^{-c}$. $\square$

La Poisson tiene todos sus momentos finitos y una cola que cae más rápido
que cualquier potencia: el grado máximo en una red de $n$ nodos crece como
$\ln n / \ln\ln n$, muy despacio. En el cuaderno, con $n = 10^4$ y $c = 4$,
el grado máximo es 14.

### 3.2 Agrupamiento y distancias

**Agrupamiento.** Dos vecinos de un nodo están conectados con probabilidad
$p$, independientemente de todo lo demás, así que

$$
C = p = \frac{c}{n-1} \to 0 .
$$

**Distancias (heurística).** Desde un nodo, hay del orden de $c$ nodos a
distancia 1, $c^2$ a distancia 2 y, en general, $c^\ell$ a distancia $\ell$,
mientras no se agoten los nodos. Se alcanza a toda la red cuando
$c^\ell \approx n$:

$$
\ell \approx \frac{\ln n}{\ln c} .
$$

Con $n = 10^4$ y $c = 4$ da 6.6; la simulación da 6.8. La distancia crece
como el **logaritmo** del tamaño: es el efecto de "mundo pequeño".

### 3.3 El componente gigante

**Proposición.** Sea $S$ la fracción de nodos en el componente más grande.
Cuando $n \to \infty$, $S$ es la mayor solución en $[0, 1]$ de

$$
S = 1 - e^{-cS} .
$$

Si $c \le 1$, la única solución es $S = 0$; si $c > 1$, hay una solución
positiva: aparece un **componente gigante**.

*Argumento.* Sea $u = 1 - S$ la probabilidad de que un nodo **no** esté en
el componente gigante. Eso ocurre si ninguno de sus vecinos está en él. En
una red grande y dispersa, los vecinos de un nodo casi no están conectados
entre sí por otro lado (no hay ciclos cortos), así que cada vecino cuenta de
forma independiente. Condicionando en el grado $k$:

$$
u = \sum_k p_k u^k = \sum_k e^{-c}\frac{(cu)^k}{k!} = e^{-c(1-u)} .
$$

Con $u = 1 - S$ queda $S = 1 - e^{-cS}$. $\square$

**[Avanzado] El umbral.** Sea $g(S) = 1 - e^{-cS}$. Siempre $g(0) = 0$, y $g$
es creciente y cóncava con $g(1) < 1$. Una segunda intersección con la recta
$S$ existe si y sólo si la pendiente en el origen es mayor que 1:
$g'(0) = c > 1$. Cerca del umbral, desarrollando
$e^{-cS} \approx 1 - cS + c^2S^2/2$,

$$
S \approx \frac{2(c - 1)}{c^2} \approx 2(c - 1),
$$

así que el parámetro de orden se anula linealmente: es una transición de
fase continua ([F.B4](repaso/fisica_estadistica.md#fb4)), la transición de
**percolación** ([F.A5](repaso/fisica_estadistica.md#fa5)).

---

## 4. Watts-Strogatz

**Modelo.** Se parte de un anillo de $n$ nodos, cada uno conectado con sus
$K$ vecinos más cercanos ($K$ par, $K/2$ de cada lado). Cada enlace se
recablea con probabilidad $p$: uno de sus extremos se cambia por un nodo
elegido al azar, sin repetir enlaces.

**Agrupamiento del anillo.** Para $p = 0$,

$$
C(0) = \frac{3(K - 2)}{4(K - 1)},
$$

que tiende a $3/4$ para $K$ grande. Con $K = 10$ da $2/3$, como en el
cuaderno.

**[Avanzado]** Para $p > 0$, un triángulo del anillo sobrevive si sus tres
enlaces sobreviven, con probabilidad $(1-p)^3$. Despreciando los triángulos
nuevos que crean los atajos (Barrat y Weigt, 2000):

$$
C(p) \approx C(0)\thinspace(1 - p)^3 .
$$

Con $p = 0.1$ da $0.486$, y la simulación da $0.476$.

**Distancias.** En el anillo, $\ell \approx n/(2K)$: crece linealmente con
$n$. Con atajos, la distancia pasa a crecer como $\ln n$ en cuanto hay del
orden de un atajo por cada "barrio", es decir, cuando $pnK/2 \gg 1$. Por eso
basta un $p$ pequeño. Hay un rango amplio de $p$ (en el cuaderno, de 0.01 a
0.1) con agrupamiento alto, como el anillo, y distancias logarítmicas, como
Erdős-Rényi: el régimen de **mundo pequeño**.

---

## 5. Barabási-Albert

### 5.1 Definición

**Modelo.** Se parte de una red pequeña. En cada paso llega un nodo nuevo
con $m$ enlaces, que se conectan a nodos existentes con probabilidad
proporcional a su grado:

$$
\Pi(i) = \frac{k_i}{\sum_j k_j} .
$$

Después de $t$ pasos hay $n \approx t$ nodos y $\sum_j k_j \approx 2mt$, así
que $\langle k\rangle = 2m$.

### 5.2 Distribución de grados (campo medio)

**Proposición (aproximación continua).** El nodo que llegó en el paso $t_i$
tiene, en el paso $t$, grado

$$
k_i(t) = m\left(\frac{t}{t_i}\right)^{1/2},
$$

y la distribución de grados tiene cola $P(K \ge k) \approx (m/k)^2$.

*Demostración.* En cada paso, el nodo $i$ recibe en promedio
$m \Pi(i) = m k_i/(2mt) = k_i/(2t)$ enlaces nuevos. Tratando $k_i$ como
continua,

$$
\frac{dk_i}{dt} = \frac{k_i}{2t},
\qquad k_i(t_i) = m
\quad\Longrightarrow\quad
k_i(t) = m\sqrt{t/t_i} .
$$

Un nodo tiene grado $k_i(t) \ge k$ si llegó antes de $t_i = t m^2/k^2$. Como
los tiempos de llegada están repartidos uniformemente entre 0 y $t$, la
fracción de nodos con grado al menos $k$ es $t_i/t = (m/k)^2$. $\square$

El resultado exacto (Dorogovtsev, Mendes y Samukhin, 2000; Krapivsky, Redner
y Leyvraz, 2000) es

$$
p_k = \frac{2m(m+1)}{k(k+1)(k+2)},
\qquad
P(K \ge k) = \frac{m(m+1)}{k(k+1)}, \qquad k \ge m .
$$

La CCDF se obtiene sumando la serie telescópica
$\frac{1}{j(j+1)(j+2)} = \frac12\left[\frac{1}{j(j+1)} - \frac{1}{(j+1)(j+2)}\right]$.
Las dos expresiones dan una cola de Pareto con exponente $\alpha = 2$ para
la CCDF ($p_k \sim k^{-3}$; repaso [P.I5](repaso/probabilidad.md#pi5)).

**Consecuencias.**

- El **primer momento** es finito: $\langle k\rangle = 2m$.
- El **segundo momento** diverge logarítmicamente con el tamaño: el nodo
  más antiguo tiene grado $k_{\max} \approx m\sqrt{n}$ y
  $\langle k^2\rangle = \sum_k k^2 p_k \approx 2m(m+1)\ln(k_{\max}/m) \approx m(m+1)\ln n$.
  Esto es lo que hace a la red robusta ante fallas al azar (sección 7).
- **Sólo la combinación** de crecimiento y conexión preferencial produce la
  cola. Con conexión uniforme ($\Pi(i) = 1/t$) los grados quedan con cola
  exponencial; sin crecimiento, la red termina siendo casi completa.
- El agrupamiento tiende a cero al crecer $n$, aunque más despacio que en
  Erdős-Rényi.

**Relación con la sesión 6.** La regla "el grado esperado crece en
proporción al grado" es crecimiento proporcional (Gibrat). Lo que la
convierte en una ley de potencias, en lugar de una lognormal, es la llegada
continua de nodos nuevos con grado $m$: el mismo papel que hace la barrera
inferior en Gibrat con barrera. Es el proceso de Yule y Simon.

---

## 6. Una red de correlaciones

A partir de la matriz de correlación $C$ de los rendimientos (sesión 7) se
define la **red de umbral**

$$
A_{ij} = \mathbb 1\left[\rho_{ij} > u\right], \qquad i \ne j .
$$

A diferencia del árbol de expansión mínima, que siempre tiene $n - 1$
enlaces y ningún ciclo, la red de umbral conserva los triángulos. Con las
95 acciones del S&P 100 (2015–2026) y $u = 0.5$:

| Red | Enlaces | $\langle k\rangle$ | $k_{\max}$ | Agrupamiento | Distancia promedio | Componente gigante |
|---|---|---|---|---|---|---|
| Mercado | 455 | 9.6 | 37 | 0.59 | 2.84 | 86% |
| Erdős-Rényi | 455 | 9.6 | 19 | 0.09 | 2.24 | 100% |
| Barabási-Albert | 450 | 9.5 | 30 | 0.16 | 2.23 | 100% |

**Por qué el agrupamiento es alto.** La correlación tiende a ser
"transitiva": si los vectores de rendimientos estandarizados $\hat z_i$,
$\hat z_j$ y $\hat z_l$ tienen norma 1, entonces
$\Vert \hat z_j - \hat z_l\Vert \le \Vert \hat z_j - \hat z_i\Vert + \Vert \hat z_i - \hat z_l\Vert$
([L.I5](repaso/algebra_lineal.md#li5)), y con $d = \sqrt{2(1-\rho)}$ eso
acota $\rho_{jl}$ por abajo cuando $\rho_{ij}$ y $\rho_{il}$ son altas. Por
ejemplo, si $\rho_{ij} = \rho_{il} = 0.75$, entonces $d_{jl} \le 2\sqrt{0.5}$
y $\rho_{jl} \ge 0$. La geometría no garantiza triángulos, pero los
favorece; los sectores hacen el resto.

Las redes de umbral fueron estudiadas por Boginski, Butenko y Pardalos
(2005) y por Onnela, Kaski y Kertész (2004), que compararon su agrupamiento
con el del árbol de expansión mínima.

---

## 7. Robustez y percolación

### 7.1 El problema

Se quita una fracción $f$ de los nodos, con todos sus enlaces, y se mide la
fracción $S(f)$ de los $n$ nodos originales que sigue en el componente
gigante. En **fallas aleatorias** los nodos se eligen al azar. En un
**ataque** se quitan en orden decreciente de grado.

En el cuaderno, con $n = 10^4$ y $\langle k\rangle = 4$:

| Red | Fallas: $S$ baja de 1% en | Ataque: se rompe en |
|---|---|---|
| Erdős-Rényi | $f \approx 0.77$ | $f \approx 0.37$ |
| Barabási-Albert | $f \approx 0.83$ | $f \approx 0.15$ |

Albert, Jeong y Barabási (2000) observaron este contraste en internet y en
la web: robustas ante fallas, frágiles ante ataques.

### 7.2 [Avanzado] El criterio de Molloy y Reed

Considera una red aleatoria con distribución de grados $p_k$ y por lo demás
al azar (el modelo de configuración de la sección 9).

**Proposición.** Hay componente gigante si y sólo si

$$
\kappa \equiv \frac{\langle k^2\rangle}{\langle k\rangle} > 2 .
$$

*Argumento.* Exploramos la red desde un nodo, siguiendo enlaces. Al llegar a
un vecino por un enlace, la probabilidad de llegar a un nodo de grado $k$ es
proporcional a $k p_k$, porque un nodo con $k$ enlaces tiene $k$ maneras de
ser alcanzado. Desde ahí puede seguir por $k - 1$ enlaces nuevos (el
**grado excedente**). El número promedio de enlaces nuevos es

$$
\sum_k (k - 1)\frac{k p_k}{\langle k\rangle}
= \frac{\langle k^2\rangle - \langle k\rangle}{\langle k\rangle}
= \kappa - 1 .
$$

La exploración es un proceso de ramificación: crece sin límite (componente
gigante) si cada rama genera en promedio más de una rama nueva, es decir, si
$\kappa - 1 > 1$. $\square$

Para Erdős-Rényi, $\langle k^2\rangle = c^2 + c$, así que $\kappa = c + 1$ y
el criterio da $c > 1$, como en la sección 3.3.

### 7.3 [Avanzado] Umbral de fallas aleatorias

**Proposición** (Cohen et al., 2000). Si se quita al azar una fracción $f$
de los nodos de una red con parámetro $\kappa$, el componente gigante
desaparece en

$$
f_c = 1 - \frac{1}{\kappa - 1} .
$$

*Demostración.* Sea $q = 1 - f$. Cada vecino de un nodo que sobrevive
sobrevive con probabilidad $q$, así que el nuevo grado $K'$ de ese nodo es
binomial con parámetros $K$ y $q$:
$\langle k'\rangle = q\langle k\rangle$ y
$\langle k'^2\rangle = q^2\langle k^2\rangle + q(1-q)\langle k\rangle$.
Entonces

$$
\kappa' = \frac{\langle k'^2\rangle}{\langle k'\rangle} = q\kappa + 1 - q .
$$

La condición $\kappa' = 2$ da $q(\kappa - 1) = 1$. $\square$

**Consecuencias.**

- Erdős-Rényi: $f_c = 1 - 1/c$. Con $c = 4$, $f_c = 0.75$; la simulación da
  0.77.
- Barabási-Albert: $\kappa \approx \frac{m+1}{2}\ln n$ crece sin límite, así
  que $f_c \to 1$: en una red infinita **no hay umbral**. Con $n = 10^3$,
  $10^4$ y $10^5$ (y $m = 2$), $\kappa$ vale 11.1, 14.9 y 19.0, y $f_c$ vale
  0.90, 0.93 y 0.94. Cerca de $f_c$ el componente gigante es diminuto, por
  eso la simulación lo "pierde" antes.
- Un **ataque** quita los nodos de mayor grado: corta la cola de $p_k$, y
  $\langle k^2\rangle$ cae muy rápido. Para Barabási-Albert basta quitar
  una fracción pequeña de los nodos (Callaway et al., 2000; Cohen et al.,
  2001).

---

## 8. [Avanzado] El Laplaciano

El Laplaciano es $L = D - A$, con $D = \mathrm{diag}(k_1, \dots, k_n)$
([L.A1](repaso/algebra_lineal.md#la1)).

**Proposición.**

1. $x^\top L x = \sum_{\text{enlaces } (i,j)} (x_i - x_j)^2 \ge 0$, así que
   los eigenvalores son reales y no negativos:
   $0 = \lambda_1 \le \lambda_2 \le \dots \le \lambda_n$.
2. La multiplicidad del eigenvalor 0 es el número de componentes conexas.

*Demostración.* (1) $x^\top L x = \sum_i k_i x_i^2 - \sum_{i,j} A_{ij}x_ix_j$.
Cada enlace $(i, j)$ aporta $x_i^2 + x_j^2$ al primer término y $2x_ix_j$ al
segundo, es decir, $(x_i - x_j)^2$. (2) $Lx = 0$ implica $x^\top L x = 0$,
así que $x_i = x_j$ en cada enlace: $x$ es constante en cada componente.
Los vectores indicadores de las componentes son una base de ese espacio. $\square$

**Conectividad algebraica.** Por (1), $\lambda_2 = \min x^\top L x$ sobre los
$x$ de norma 1 ortogonales a $\mathbf 1$. Si la red se puede partir en dos
grupos con pocos enlaces entre ellos, el vector que vale $+1$ en un grupo y
$-1$ en el otro (centrado y normalizado) da un valor pequeño. Por eso el
eigenvector de $\lambda_2$ (**vector de Fiedler**) sugiere una partición:
los nodos con componente positiva de un lado, los negativos del otro.

En el cuaderno, la red de umbral del mercado con $u = 0.5$ tiene 11
componentes (10 son acciones aisladas o pares) y su Laplaciano tiene 11
ceros. En el componente gigante, el vector de Fiedler separa a las
acciones defensivas (consumo básico, servicios públicos, buena parte de
salud) de las cíclicas, igual que el segundo eigenvector de la matriz de
correlación en la sesión 7.

---

## 9. [Avanzado] Modelo de configuración y asortatividad

### 9.1 El modelo de configuración

**Modelo.** Se fija una sucesión de grados $k_1, \dots, k_n$. Cada nodo
recibe $k_i$ "medios enlaces" y se emparejan todos al azar. Es la red más
aleatoria posible con esos grados; Erdős-Rényi es el caso de grados Poisson.

La probabilidad de un enlace entre $i$ y $j$ es $p_{ij} \approx k_ik_j/(2m)$,
y el agrupamiento es

$$
C = \frac{1}{n}\frac{\left(\langle k^2\rangle - \langle k\rangle\right)^2}{\langle k\rangle^3} .
$$

Tiende a cero si $\langle k^2\rangle$ es finito, pero puede ser apreciable
con grados muy desiguales y $n$ pequeño. Con los grados del mercado
($u = 0.5$) la fórmula da 0.35 y la red recableada del cuaderno tiene 0.31,
contra 0.59 de la red real. El resto del agrupamiento es estructura de
sectores, no consecuencia de los grados.

En la práctica, la red se muestrea recableando la real: se eligen dos
enlaces $(a, b)$ y $(c, d)$ y se cambian por $(a, d)$ y $(c, b)$, si eso no
crea enlaces repetidos ni lazos. Cada cambio conserva todos los grados.

### 9.2 Asortatividad

La **asortatividad** $r$ es el coeficiente de correlación de Pearson entre
los grados de los dos extremos de un enlace, calculado sobre todos los
enlaces en los dos sentidos (Newman, 2002). Va de $-1$ a 1.

- $r > 0$: los nodos de grado alto se conectan entre sí. Típico de redes
  sociales. El mercado tiene $r = 0.12$: un núcleo financiero denso.
- $r < 0$: los concentradores se conectan con nodos de grado bajo. Típico
  de redes tecnológicas y biológicas, y de redes simples con grados muy
  desiguales: si sólo puede haber un enlace por par, los concentradores no
  encuentran suficientes concentradores. Por eso la red recableada tiene
  $r = -0.23$.

---

## 10. Resumen formal

| Concepto | Resultado |
|---|---|
| Caminos y triángulos | $[A^\ell]_{ij}$ = caminos de longitud $\ell$; triángulos $= \mathrm{tr}(A^3)/6$ |
| Agrupamiento | $C_i = [A^3]_{ii}/(k_i(k_i-1))$; global $= \mathrm{tr}(A^3)/\sum_i k_i(k_i-1)$ |
| Erdős-Rényi | Grados Poisson($c$); $C = c/(n-1)$; $\ell \approx \ln n/\ln c$; gigante si $c > 1$, con $S = 1 - e^{-cS}$ |
| Watts-Strogatz | $C(p) \approx \frac{3(K-2)}{4(K-1)}(1-p)^3$; distancias logarítmicas con pocos atajos |
| Barabási-Albert | $P(K \ge k) = m(m+1)/(k(k+1))$; $\langle k^2\rangle \sim \ln n$ |
| Molloy-Reed | Gigante si $\kappa = \langle k^2\rangle/\langle k\rangle > 2$ |
| Fallas aleatorias | $f_c = 1 - 1/(\kappa - 1)$; libre de escala: $f_c \to 1$ |
| Laplaciano | Ceros = componentes; vector de Fiedler = partición |

---

## 11. Ejercicios formales

**Ruta básica**

- **F1.** Para la red de cuatro nodos de la sección 1.2, calcula a mano
  $A^2$ y $A^3$. Verifica que $[A^2]_{ii} = k_i$ y obtén $C_i$ para cada
  nodo.
- **F2.** Una red de Erdős-Rényi tiene $n = 10^6$ nodos y $c = 10$. ¿Cuál es
  su agrupamiento? ¿Cuál es la distancia típica entre dos nodos?
- **F3.** En una red de Barabási-Albert con $m = 3$, ¿qué fracción de los
  nodos tiene grado 30 o más? ¿Y grado 300 o más? Compara con una Poisson
  de media 6.
- **F4.** ¿Qué centralidad usarías para identificar a un banco cuya quiebra
  cortaría el crédito entre dos regiones? ¿Y para identificar al banco con
  más contrapartes?

**Ruta avanzada**

- **F5.** Demuestra que en una estrella de $n$ nodos (un centro conectado
  con $n - 1$ hojas) el centro tiene intermediación normalizada 1 y las
  hojas 0. ¿Cuánto vale el agrupamiento global?
- **F6.** Deriva $C(0) = 3(K-2)/(4(K-1))$ para el anillo de Watts-Strogatz.
  Pista: cuenta cuántos pares de vecinos de un nodo están a distancia menor
  o igual a $K/2$ sobre el anillo.
- **F7.** Generaliza la sección 5.2 a una conexión preferencial con
  "atractivo inicial": $\Pi(i) \propto k_i + a$. ¿Qué exponente tiene la
  cola?
- **F8.** Usa el criterio de Molloy y Reed para una red con
  $p_k \propto k^{-\gamma}$ y $k \ge 1$. ¿Para qué valores de $\gamma$ hay
  componente gigante sin importar cuántos nodos se quiten al azar (en el
  límite $n \to \infty$)?
- **F9.** Demuestra que $\lambda_2 > 0$ si y sólo si la red es conexa, y
  calcula el espectro del Laplaciano de un anillo de $n$ nodos.

---

## Referencias

- Albert, R., Jeong, H., & Barabási, A.-L. (2000). Error and attack
  tolerance of complex networks. *Nature*, 406(6794), 378–382.
  https://doi.org/10.1038/35019019
- Barabási, A.-L., & Albert, R. (1999). Emergence of scaling in random
  networks. *Science*, 286(5439), 509–512.
  https://doi.org/10.1126/science.286.5439.509
- Barrat, A., & Weigt, M. (2000). On the properties of small-world network
  models. *The European Physical Journal B*, 13(3), 547–560.
  https://doi.org/10.1007/s100510050067
- Boginski, V., Butenko, S., & Pardalos, P. M. (2005). Statistical analysis
  of financial networks. *Computational Statistics & Data Analysis*, 48(2),
  431–443. https://doi.org/10.1016/j.csda.2004.02.004
- Callaway, D. S., Newman, M. E. J., Strogatz, S. H., & Watts, D. J. (2000).
  Network robustness and fragility: Percolation on random graphs. *Physical
  Review Letters*, 85(25), 5468–5471.
  https://doi.org/10.1103/PhysRevLett.85.5468
- Cohen, R., Erez, K., ben-Avraham, D., & Havlin, S. (2000). Resilience of
  the Internet to random breakdowns. *Physical Review Letters*, 85(21),
  4626–4628. https://doi.org/10.1103/PhysRevLett.85.4626
- Cohen, R., Erez, K., ben-Avraham, D., & Havlin, S. (2001). Breakdown of
  the Internet under intentional attack. *Physical Review Letters*, 86(16),
  3682–3685. https://doi.org/10.1103/PhysRevLett.86.3682
- Dorogovtsev, S. N., Mendes, J. F. F., & Samukhin, A. N. (2000). Structure
  of growing networks with preferential linking. *Physical Review Letters*,
  85(21), 4633–4636. https://doi.org/10.1103/PhysRevLett.85.4633
- Jackson, M. O. (2008). *Social and Economic Networks*. Princeton
  University Press.
- Krapivsky, P. L., Redner, S., & Leyvraz, F. (2000). Connectivity of
  growing random networks. *Physical Review Letters*, 85(21), 4629–4632.
  https://doi.org/10.1103/PhysRevLett.85.4629
- Molloy, M., & Reed, B. (1995). A critical point for random graphs with a
  given degree sequence. *Random Structures & Algorithms*, 6(2–3), 161–180.
  https://doi.org/10.1002/rsa.3240060204
- Newman, M. E. J. (2002). Assortative mixing in networks. *Physical Review
  Letters*, 89(20), 208701. https://doi.org/10.1103/PhysRevLett.89.208701
- Newman, M. (2018). *Networks* (2.ª ed.). Oxford University Press.
  https://doi.org/10.1093/oso/9780198805090.001.0001
- Onnela, J.-P., Kaski, K., & Kertész, J. (2004). Clustering and
  information in correlation based financial networks. *The European
  Physical Journal B*, 38(2), 353–362.
  https://doi.org/10.1140/epjb/e2004-00128-7
- Padgett, J. F., & Ansell, C. K. (1993). Robust action and the rise of the
  Medici, 1400–1434. *American Journal of Sociology*, 98(6), 1259–1319.
  https://doi.org/10.1086/230190
- Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
  'small-world' networks. *Nature*, 393(6684), 440–442.
  https://doi.org/10.1038/30918
