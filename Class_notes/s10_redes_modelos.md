# Sesión 10: Redes: medidas y modelos. Notas formales

> **Cómo usar estas notas.** La clase combina el cuaderno
> (`Code/s10_redes_modelos.jl`) y el pizarrón. El **plan de clase** de abajo
> dice qué se proyecta, qué secciones de estas notas se escriben en el
> pizarrón, cuándo se trabaja en parejas y qué queda para leer. En el
> pizarrón van las definiciones, las proposiciones y los resultados en
> recuadro; las demostraciones, sólo donde el plan lo indica. El texto que las
> rodea se explica de palabra. Como todo está aquí, en clase no hace falta
> copiar. El cuerpo de las notas es la **ruta básica**. Las demostraciones y
> los temas de la ruta avanzada, opcionales, están al final, en **Para
> profundizar**, y no se ven en clase; donde estaban queda una línea que los
> enlaza.
>
> Repaso recomendado: [L.B4](repaso/algebra_lineal.md#lb4) (matriz de
> adyacencia), [P.B1](repaso/probabilidad.md#pb1) (distribuciones),
> [C.B2](repaso/calculo.md#cb2) (escalas log-log),
> [F.B4](repaso/fisica_estadistica.md#fb4) (transiciones de fase),
> [P.I5](repaso/probabilidad.md#pi5) (leyes de potencias),
> [L.I4](repaso/algebra_lineal.md#li4) (potencias de matrices y centralidad),
> [S.I3](repaso/estadistica.md#si3) (comparar datos con una distribución);
> para la ruta avanzada, [L.A5](repaso/algebra_lineal.md#la5)
> (Perron-Frobenius) y [F.A5](repaso/fisica_estadistica.md#fa5) (percolación).
>
> Lecturas: Newman (2018), cap. 7. Complementaria: Newman (2018), caps. 6 y 11;
> Requião da Cunha (2021), sec. 8.1.

**Plan de clase** (120 minutos)

| Minutos | Dónde | Qué | Secciones |
|---|---|---|---|
| 0–10 | Proyector | El fenómeno: ¿por qué ganaron los Médici? | Cuaderno, sección 0 |
| 10–20 | Pizarra | Repaso: la matriz de adyacencia; grados y $m = \tfrac12\sum_i k_i$; $[A^\ell]_{ij}$ cuenta caminos, con la demostración, y los triángulos | 1.1–1.2 |
| 20–30 | Pizarra | Las cuatro centralidades (la tabla) y el agrupamiento local y global | 2.1–2.2 |
| 30–42 | Pizarra | Erdős-Rényi: grados binomiales y el límite de Poisson (enunciado; es el ejercicio 6 del cuaderno); $C = p$; $\ell \approx \ln n/\ln c$; el componente gigante, con el argumento | 3.1–3.3 |
| 42–45 | Pizarra | Watts-Strogatz: el agrupamiento del anillo (enunciado; es el ejercicio 7) y el mundo pequeño, de palabra | 4 |
| 45–55 | Pizarra | Barabási-Albert: la conexión preferencial; el grado por campo medio, con la demostración, y la cola con $\alpha = 2$ | 5.1–5.2 |
| 55–100 | Parejas | Laboratorio: medir la red de Florencia y comparar los tres modelos; ejercicios de código 1 y 2 | Cuaderno, secciones 2 y 3 y ejercicios 1–2 |
| 100–120 | Pizarra | Cierre: puesta en común de las preguntas 3 a 5 y la tabla del resumen | Cuaderno, preguntas 3–5; notas, 6 |
| Después | Lectura | Las familias de Florencia (2.3), el resultado exacto de Barabási-Albert y sus consecuencias (5.2) y Para profundizar | |

**Pregunta de la sesión.** ¿Por qué ganaron los Médici? ¿Qué distingue la
posición de un nodo en una red, y qué redes producen nodos así?

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

**[Avanzado]** Proposición: en [Para profundizar, P.1](#pp1).

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
\boxed{\ \ell \approx \frac{\ln n}{\ln c}\ }
$$

Con $n = 10^4$ y $c = 4$ da 6.6; la simulación da 6.8. La distancia crece
como el **logaritmo** del tamaño: es el efecto de "mundo pequeño".

### 3.3 El componente gigante

**Proposición.** Sea $S$ la fracción de nodos en el componente más grande.
Cuando $n \to \infty$, $S$ es la mayor solución en $[0, 1]$ de

$$
\boxed{\ S = 1 - e^{-cS}\ }
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

**[Avanzado]** El umbral: en [Para profundizar, P.2](#pp2).

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

**[Avanzado]** El agrupamiento para p > 0: en [Para profundizar, P.3](#pp3).

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
\boxed{\ k_i(t) = m\sqrt{t/t_i}\ }
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
  Esto es lo que hace a la red robusta ante fallas al azar
  ([sesión 11](s11_redes_mercado.md)).
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

## 6. Resumen formal

| Concepto | Resultado |
|---|---|
| Caminos y triángulos | $[A^\ell]_{ij}$ = caminos de longitud $\ell$; triángulos $= \mathrm{tr}(A^3)/6$ |
| Agrupamiento | $C_i = [A^3]_{ii}/(k_i(k_i-1))$; global $= \mathrm{tr}(A^3)/\sum_i k_i(k_i-1)$ |
| Erdős-Rényi | Grados Poisson($c$); $C = c/(n-1)$; $\ell \approx \ln n/\ln c$; gigante si $c > 1$, con $S = 1 - e^{-cS}$ |
| Watts-Strogatz | $C(p) \approx \frac{3(K-2)}{4(K-1)}(1-p)^3$; distancias logarítmicas con pocos atajos |
| Barabási-Albert | $P(K \ge k) = m(m+1)/(k(k+1))$; $\langle k^2\rangle \sim \ln n$ |

---

## 7. Ejercicios formales

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

---

## Para profundizar

Material opcional de la ruta avanzada: demostraciones y temas que no se ven
en clase. Cada apartado indica de qué sección de las notas viene.

<a id="pp1"></a>

### P.1 Proposición (sección 2.1)

**Proposición.** Si la red es conexa, el eigenvector del mayor
eigenvalor $\lambda_1$ de $A$ tiene todas sus componentes del mismo signo
(y se puede tomar positivo), y es el único con esa propiedad.

*Demostración.* Una red conexa tiene $A$ no negativa e irreducible; el
teorema de Perron-Frobenius ([L.A5](repaso/algebra_lineal.md#la5)) da un
eigenvalor real simple $\lambda_1 = r(A)$ con eigenvector positivo. Como
$A$ es simétrica, los eigenvectores de eigenvalores distintos son
ortogonales: ninguno puede ser ortogonal a un vector positivo y tener
todas sus componentes positivas. $\square$

<a id="pp2"></a>

### P.2 El umbral (sección 3.3)

**El umbral.** Sea $g(S) = 1 - e^{-cS}$. Siempre $g(0) = 0$, y $g$
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

<a id="pp3"></a>

### P.3 El agrupamiento para p > 0 (sección 4)

Para $p > 0$, un triángulo del anillo sobrevive si sus tres
enlaces sobreviven, con probabilidad $(1-p)^3$. Despreciando los triángulos
nuevos que crean los atajos (Barrat y Weigt, 2000):

$$
C(p) \approx C(0)\thinspace(1 - p)^3 .
$$

Con $p = 0.1$ da $0.486$, y la simulación da $0.476$.

---

## Referencias

- Barabási, A.-L., & Albert, R. (1999). Emergence of scaling in random
  networks. *Science*, 286(5439), 509–512.
  https://doi.org/10.1126/science.286.5439.509
- Barrat, A., & Weigt, M. (2000). On the properties of small-world network
  models. *The European Physical Journal B*, 13(3), 547–560.
  https://doi.org/10.1007/s100510050067
- Dorogovtsev, S. N., Mendes, J. F. F., & Samukhin, A. N. (2000). Structure
  of growing networks with preferential linking. *Physical Review Letters*,
  85(21), 4633–4636. https://doi.org/10.1103/PhysRevLett.85.4633
- Jackson, M. O. (2008). *Social and Economic Networks*. Princeton
  University Press.
- Krapivsky, P. L., Redner, S., & Leyvraz, F. (2000). Connectivity of
  growing random networks. *Physical Review Letters*, 85(21), 4629–4632.
  https://doi.org/10.1103/PhysRevLett.85.4629
- Newman, M. (2018). *Networks* (2.ª ed.). Oxford University Press.
  https://doi.org/10.1093/oso/9780198805090.001.0001
- Padgett, J. F., & Ansell, C. K. (1993). Robust action and the rise of the
  Medici, 1400–1434. *American Journal of Sociology*, 98(6), 1259–1319.
  https://doi.org/10.1086/230190
- Requião da Cunha, C. (2021). *Introduction to Econophysics: Contemporary
  Approaches with Python Simulations*. CRC Press.
  https://doi.org/10.1201/9781003127956
- Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
  'small-world' networks. *Nature*, 393(6684), 440–442.
  https://doi.org/10.1038/30918
