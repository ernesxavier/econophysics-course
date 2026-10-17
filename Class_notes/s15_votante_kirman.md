# Sesión 15: Votante y Kirman: imitación y herding. Notas formales

> **Cómo leer estas notas.** Siguen el orden del programa: primero se observa
> en el cuaderno (`Code/s15_votante_kirman.jl`), después se formaliza aquí. Lo que no
> tiene marca es **ruta básica**. Las secciones y demostraciones marcadas
> **[Avanzado]** son opcionales.
>
> Repaso recomendado: [E.B2](repaso/procesos_estocasticos.md#eb2) (caminata
> aleatoria), [E.B3](repaso/procesos_estocasticos.md#eb3) (promedios),
> [L.B4](repaso/algebra_lineal.md#lb4) (matriz de adyacencia y grado),
> [C.B4](repaso/calculo.md#cb4) (función logística),
> [F.B3](repaso/fisica_estadistica.md#fb3) (temperatura y ruido),
> [F.B4](repaso/fisica_estadistica.md#fb4) (transiciones de fase),
> [E.I2](repaso/procesos_estocasticos.md#ei2) (martingalas),
> [E.I3](repaso/procesos_estocasticos.md#ei3) (cadenas de Markov y estados
> absorbentes), [F.I3](repaso/fisica_estadistica.md#fi3) (Ising) y
> [F.I4](repaso/fisica_estadistica.md#fi4) (campo medio); para la ruta
> avanzada, [E.A2](repaso/procesos_estocasticos.md#ea2) (lema de Itô),
> [E.A3](repaso/procesos_estocasticos.md#ea3) (ecuación maestra),
> [E.A4](repaso/procesos_estocasticos.md#ea4) (Fokker-Planck) y
> [E.A6](repaso/procesos_estocasticos.md#ea6) (ruina del jugador y tiempos de
> salida).
>
> Lecturas: Sen & Chakrabarti (2014), cap. 3. Complementaria: Kirman (2011),
> caps. 2–3. Ruta avanzada: Slanina (2014), cap. 3; Castellano, Fortunato &
> Loreto (2009); Kirman (1993); Sood & Redner (2005).

---

## 1. El modelo de votante

### 1.1 Definición y cadena de Markov en la red completa

$N$ agentes ocupan los nodos de un grafo conexo $G = (V, E)$. Cada agente
$i \in V$ tiene una opinión binaria $s_i \in \lbrace 0, 1\rbrace$. En cada paso
elemental ($\Delta t = 1/N$ de barrido):

1. Se elige un agente $i$ al azar con probabilidad uniforme $1/N$.
2. Se elige un vecino $j$ de $i$ al azar con probabilidad $1/k_i$, donde $k_i$
   es el grado de $i$ ([L.B4](repaso/algebra_lineal.md#lb4)).
3. El agente $i$ copia la opinión de $j$: $s_i \leftarrow s_j$.

Los dos estados de consenso, $s = \mathbf 0$ y $s = \mathbf 1$, son
**estados absorbentes** ([E.I3](repaso/procesos_estocasticos.md#ei3)): una vez
que todos los agentes coinciden, ningún paso posterior puede cambiar ninguna
opinión.

**Red completa.** Cuando todos son vecinos de todos ($k_i = N - 1$), el estado
macroscópico queda determinado por el número $k = \sum_{i=1}^N s_i \in
\lbrace 0, 1, \dots, N\rbrace$ de agentes con opinión $1$ (o la fracción
$x = k/N$). Para que $k$ aumente en una unidad debe elegirse un agente en $0$
(con probabilidad $(N-k)/N$) y que el vecino elegido esté en $1$, lo que
ocurre con probabilidad $k/(N-1)$. Para que disminuya en una unidad debe
ocurrir lo simétrico. Por tanto, las
probabilidades de transición por paso elemental son

$$
p_k = \mathbb P(k \to k+1) = \frac{k(N-k)}{N(N-1)}, \qquad
q_k = \mathbb P(k \to k-1) = \frac{k(N-k)}{N(N-1)},
$$

y $\mathbb P(k \to k) = 1 - 2p_k$.

### 1.2 Martingala y probabilidad de victoria

Como $p_k = q_k$ para todo $k$, el cambio esperado en un paso es nulo:

$$
\mathbb E\left[k_{n+1} - k_n \mid k_n = k\right] = (+1)p_k + (-1)q_k = 0 .
$$

La fracción $x_n = k_n/N$ es una **martingala** acotada
([E.I2](repaso/procesos_estocasticos.md#ei2)) con estados absorbentes en $0$ y
$1$.

**Proposición** (probabilidad de consenso en $1$). Partiendo de $k_0$ agentes
con opinión $1$, la probabilidad $v_{k_0}$ de que toda la población termine en
el consenso $1$ es

$$
v_{k_0} = \frac{k_0}{N} = x_0 .
$$

*Demostración.* Condicionando en el primer paso desde $k \in
\lbrace 1, \dots, N-1\rbrace$,

$$
v_k = p_k v_{k+1} + q_k v_{k-1} + (1 - p_k - q_k) v_k .
$$

Como $p_k = q_k > 0$, esto se simplifica a $v_k = \frac12(v_{k+1} + v_{k-1})$,
es decir, $v_{k+1} - v_k = v_k - v_{k-1}$: la sucesión $v_k$ es lineal en $k$.
Con las condiciones de frontera absorbentes $v_0 = 0$ y $v_N = 1$, se obtiene
$v_k = k/N$ (la ruina del jugador simétrica,
[E.A6](repaso/procesos_estocasticos.md#ea6)). Equivalentemente, por el teorema
de paro opcional para la martingala $x_n$ en el tiempo de absorción $T$,
$\mathbb E[x_T] = x_0$; como $x_T \in \lbrace 0, 1\rbrace$,
$\mathbb E[x_T] = 1\cdot v_{k_0} + 0\cdot(1 - v_{k_0}) = x_0$. $\square$

### 1.3 Tiempo esperado hasta el consenso en la red completa [Avanzado]

Sea $\tau_k$ el número esperado de **barridos** (donde $1$ barrido $= N$ pasos
elementales) hasta alcanzar cualquiera de los dos consensos ($k = 0$ o
$k = N$), partiendo de $k$ agentes en $1$.

**Proposición.** Para $N$ grande y $x = k/N \in (0, 1)$, el tiempo esperado
hasta el consenso (en barridos) es

$$
\tau(x) = -(N - 1)\bigl[x\ln x + (1 - x)\ln(1 - x)\bigr] .
$$

En particular, partiendo de una población dividida en partes iguales ($x_0 =
1/2$), $\tau(1/2) = (N - 1)\ln 2 \approx 0.693\thinspace(N - 1)$ barridos.

*Demostración.* En un paso elemental ($\Delta t = 1/N$ de barrido), condicionando
en el resultado del paso:

$$
\tau_k = \frac{1}{N} + p_k \tau_{k+1} + q_k \tau_{k-1} + (1 - 2p_k)\tau_k ,
$$

con $\tau_0 = \tau_N = 0$. Reordenando y sustituyendo $p_k =
\frac{k(N-k)}{N(N-1)}$:

$$
\tau_{k+1} - 2\tau_k + \tau_{k-1} = -\frac{1}{N p_k}
= -\frac{N - 1}{k(N - k)} .
$$

Escribiendo $\tau_k = \tau(x)$ con $x = k/N$ y aproximando la segunda diferencia
por $\tau_{k+1} - 2\tau_k + \tau_{k-1} \approx \frac{1}{N^2}\tau''(x)$ para $N$
grande, obtenemos la ecuación diferencial ordinaria

$$
\frac{1}{N^2}\tau''(x) = -\frac{N - 1}{N^2 x(1 - x)}
\quad\Longrightarrow\quad
\tau''(x) = -(N - 1)\left(\frac{1}{x} + \frac{1}{1 - x}\right) .
$$

Integrando dos veces se obtiene $\tau(x) = -(N - 1)[x\ln x + (1-x)\ln(1-x)] +
C_1 x + C_2$. Las condiciones de frontera $\tau(0) = \tau(1) = 0$ implican
$C_1 = C_2 = 0$. $\square$

**Fórmula exacta.** La ecuación de diferencias se resuelve sin aproximar,
sumando dos veces y usando $\frac{1}{j(N-j)} = \frac{1}{N}\left(\frac{1}{j} + \frac{1}{N-j}\right)$:

$$
\tau_k = \frac{N-1}{N}\left[(N-k)\sum_{j=1}^{k-1}\frac{1}{N-j} + k\sum_{j=k}^{N-1}\frac{1}{j}\right] .
$$

Con $\sum_{j \le n} 1/j \approx \ln n + \gamma_E$ se recupera la fórmula
continua. La diferencia es pequeña: con $N = 80$ y $k = 24$, la exacta da
47.8 barridos y la continua 48.3. El cuaderno usa la exacta.

Nótese que $S(x) = -x\ln x - (1-x)\ln(1-x)$ es la entropía de Shannon de una
variable binaria: el tiempo esperado hasta el consenso es proporcional a la
incertidumbre inicial sobre qué opinión ganará.

### 1.4 Efecto de la topología de la red y leyes de escala

En un grafo arbitrario con matriz de adyacencia $A_{ij}$ y grados $k_i = \sum_j
A_{ij}$ ([L.B4](repaso/algebra_lineal.md#lb4)), condicionando en el estado
actual $s$, la probabilidad de que el nodo $i$ sea elegido y cambie de estado
en un paso elemental hace que

$$
\mathbb E\left[s_i(t + \Delta t) - s_i(t) \mid s(t)\right]
= \frac{1}{N}\left(\frac{1}{k_i}\sum_{j=1}^N A_{ij} s_j(t) - s_i(t)\right) .
$$

Si la red no es regular ($k_i$ varía entre nodos), la suma simple $\sum_i s_i$
**ya no es una martingala**, porque los nodos de alto grado aparecen como
vecinos con mayor frecuencia que como nodos actualizados.

**Proposición** (martingala ponderada por el grado; Suchecki, Eguíluz & San
Miguel, 2005; Sood & Redner, 2005). En cualquier grafo conexo no dirigido, la
opinión media ponderada por el grado

$$
\bar s_w = \frac{\sum_{i=1}^N k_i s_i}{\sum_{i=1}^N k_i}
$$

es una martingala. Por tanto, la probabilidad de alcanzar el consenso $1$ desde
una configuración inicial $s(0)$ es $\bar s_w(0)$.

*Demostración.* Multiplicando por $k_i$ y sumando sobre $i$:

$$
\mathbb E\left[\sum_{i=1}^N k_i s_i(t + \Delta t) - \sum_{i=1}^N k_i s_i(t) \ \middle|\ s(t)\right]
= \frac{1}{N}\sum_{i=1}^N \sum_{j=1}^N A_{ij} s_j(t) - \frac{1}{N}\sum_{i=1}^N k_i s_i(t) = 0 ,
$$

porque $A_{ij} = A_{ji}$ y $\sum_i A_{ij} = k_j$. $\square$

Un concentrador (*hub*) con grado $k_i$ grande tiene peso proporcional a $k_i$
en el destino final de la población. Además, la dimensión y la heterogeneidad
de grados cambian el escalamiento del tiempo medio al consenso $\langle\tau\rangle$
(en barridos) con el tamaño $N$:

| Topología de la red | Tiempo al consenso $\langle\tau\rangle$ (barridos) | Mecanismo físico |
|---|---|---|
| Anillo / red regular 1D | $\mathcal O(N^2)$ | Fronteras de dominio que difunden como caminatas aleatorias hasta aniquilarse |
| Cuadrícula 2D | $\mathcal O(N \ln N)$ | Dimensión crítica $d_c = 2$: los dominios engrosan (*coarsening*) sólo logarítmicamente |
| Redes regulares $d \ge 3$, Erdős-Rényi, completa | $\mathcal O(N)$ | Los dominios no engrosan; el sistema sigue mezclado y llega al consenso por una fluctuación global, como en la red completa |
| Libre de escala ($\gamma = 3$, Barabási-Albert) | $\mathcal O(N / \ln N)$ | Los *hubs* sincronizan rápidamente a sus vecinos (Sood & Redner, 2005) |

---

## 2. El modelo de hormigas de Kirman

### 2.1 Tasas de transición: cambio autónomo e imitación

Kirman (1993) propuso un modelo estocástico inspirado en experimentos con
hormigas que eligen entre dos fuentes de alimento idénticas, el mismo
fenómeno que Becker (1991) describió con dos restaurantes. Después se aplicó
a inversionistas que alternan entre reglas fundamentalistas y chartistas;
Alfarano, Lux y Wagner (2005) estiman una versión asimétrica con datos
financieros. A
diferencia del modelo de votante, cada agente tiene una propensión positiva
$a > 0$ a cambiar de opinión por cuenta propia, evitando que los estados
$k = 0$ y $k = N$ sean absorbentes ([F.B3](repaso/fisica_estadistica.md#fb3)).

Sea $k \in \lbrace 0, 1, \dots, N\rbrace$ el número de agentes en el estado
$1$. En cada paso elemental se elige un agente al azar:

- Con peso $N - 1$ imita a uno de los demás $N - 1$ agentes elegido al azar.
- Con peso $2a$ decide de forma autónoma y elige $0$ o $1$ con igual
  probabilidad (peso $a$ para cada opción).

Si el agente elegido estaba en $0$ (probabilidad $(N-k)/N$), pasa a $1$ con
probabilidad $(a + k)/(N - 1 + 2a)$. Si estaba en $1$ (probabilidad $k/N$),
pasa a $0$ con probabilidad $(a + N - k)/(N - 1 + 2a)$. Las probabilidades de
transición en un paso son

$$
W(k \to k+1) = \frac{N-k}{N}\thinspace\frac{a + k}{N - 1 + 2a}, \qquad
W(k \to k-1) = \frac{k}{N}\thinspace\frac{a + N - k}{N - 1 + 2a} .
$$

Cuando $a \to 0$ se recupera exactamente el modelo de votante en la red
completa.

### 2.2 Ecuación maestra y distribución Beta-binomial exacta [Avanzado]

Como el proceso es un proceso de nacimiento y muerte irreducible en el espacio
finito $\lbrace 0, 1, \dots, N\rbrace$, su distribución estacionaria $P_k$ es
única y satisface la condición de **balance detallado**
([E.A3](repaso/procesos_estocasticos.md#ea3)):

$$
P_k W(k \to k+1) = P_{k+1} W(k+1 \to k), \qquad k = 0, 1, \dots, N-1 .
$$

**Proposición.** La distribución estacionaria exacta del modelo de Kirman para
cualquier $N \ge 2$ y $a > 0$ es la distribución **Beta-binomial** simétrica:

$$
P_k = \binom{N}{k}\frac{\mathrm{B}(k + a,\ N - k + a)}{\mathrm{B}(a, a)},
\qquad k = 0, 1, \dots, N ,
$$

donde $\mathrm{B}(u, v) = \frac{\Gamma(u)\Gamma(v)}{\Gamma(u+v)}$ es la función
Beta.

*Demostración.* Del balance detallado, el cociente entre probabilidades
consecutivas es

$$
\frac{P_{k+1}}{P_k} = \frac{W(k \to k+1)}{W(k+1 \to k)}
= \frac{N - k}{k + 1}\thinspace\frac{a + k}{a + N - 1 - k} .
$$

Iterando desde $k = 0$ y usando $\Gamma(z+1) = z\Gamma(z)$,

$$
P_k = P_0 \prod_{j=0}^{k-1} \frac{N - j}{j + 1}
\prod_{j=0}^{k-1} \frac{a + j}{a + N - 1 - j}
= P_0 \binom{N}{k} \frac{\Gamma(k+a)/\Gamma(a)}{\Gamma(N+a)/\Gamma(N-k+a)} .
$$

Normalizando con la identidad de Pólya $\sum_{k=0}^N \binom{N}{k}
\mathrm{B}(k+a, N-k+a) = \mathrm{B}(a, a)$, se obtiene la fórmula enunciada.
$\square$

**Caso crítico $a = 1$.** Cuando $a = 1$, el factor $(a + k)/(a + N - 1 - k)$
en el cociente $P_{k+1}/P_k$ se cancela con $(N - k)/(k + 1)$ al escribirlo en
términos de factoriales, o bien observando que $\binom{N}{k}\mathrm{B}(k+1,
N-k+1) = \frac{N!}{k!(N-k)!}\frac{k!(N-k)!}{(N+1)!} = \frac{1}{N+1}$. Así, en
$a = 1$ la distribución estacionaria es **exactamente uniforme** para todo $N$:

$$
P_k = \frac{1}{N + 1}, \qquad k = 0, 1, \dots, N \quad (a = 1) .
$$

### 2.3 Límite continuo: Fokker-Planck y distribución $\mathrm{Beta}(a,a)$ [Avanzado]

Sea $x = k/N \in [0, 1]$ la fracción de agentes en el estado $1$, y midamos el
tiempo en unidades de $N(N - 1 + 2a)/2 \approx N^2/2$ pasos elementales. En un
paso elemental, el cambio $\Delta x = \pm 1/N$ tiene momentos condicionales

$$
\mathbb E[\Delta x \mid x] = \frac{W(k \to k+1) - W(k \to k-1)}{N}
= \frac{2a(1/2 - x)}{N(N - 1 + 2a)} ,
$$

$$
\mathbb E\left[(\Delta x)^2 \mid x\right] = \frac{W(k \to k+1) + W(k \to k-1)}{N^2}
= \frac{a + 2N x(1 - x)}{N^2(N - 1 + 2a)} .
$$

En el límite $N \to \infty$ con $a > 0$ fijo (desarrollando hasta orden
$1/N^2$, expansión de Kramers-Moyal,
[E.A3](repaso/procesos_estocasticos.md#ea3)), la densidad $p(x, t)$ obedece la
**ecuación de Fokker-Planck** de la difusión $dx = A(x)\thinspace dt +
B(x)\thinspace dW_t$ ([E.A4](repaso/procesos_estocasticos.md#ea4); difusión de
Wright-Fisher / Jacobi):

$$
\frac{\partial p(x,t)}{\partial t}
= -\frac{\partial}{\partial x}\bigl[A(x)p(x,t)\bigr] +
\frac12\frac{\partial^2}{\partial x^2}\bigl[B(x)^2 p(x,t)\bigr],
$$

con deriva e intensidad de ruido

$$
A(x) = a\left(\frac12 - x\right), \qquad B(x)^2 = x(1 - x) .
$$

**Proposición.** La densidad estacionaria de la ecuación de Fokker-Planck en
$(0, 1)$ es la distribución $\mathrm{Beta}(a, a)$:

$$
p_{\text{est}}(x) = \frac{1}{\mathrm{B}(a, a)}\thinspace x^{a - 1}(1 - x)^{a - 1},
\qquad x \in (0, 1) .
$$

*Demostración.* Con flujo nulo en las fronteras ($J(x) = A(x)p(x) -
\frac12\frac{d}{dx}[B(x)^2 p(x)] = 0$),

$$
\frac{d}{dx}\bigl[x(1-x)p(x)\bigr] = 2a\left(\frac12 - x\right)p(x)
= a(1 - 2x)p(x) .
$$

Expandiendo la derivada izquierda: $(1 - 2x)p(x) + x(1-x)p'(x) = a(1-2x)p(x)$,
de donde

$$
\frac{p'(x)}{p(x)} = \frac{(a - 1)(1 - 2x)}{x(1 - x)}
= \frac{a - 1}{x} - \frac{a - 1}{1 - x} .
$$

Integrando: $\ln p(x) = (a-1)\ln x + (a-1)\ln(1-x) + \text{cte}$, es decir,
$p(x) \propto x^{a-1}(1-x)^{a-1}$. $\square$

### 2.4 Transición entre unimodalidad y bimodalidad: Kirman frente a Ising

La forma de $p_{\text{est}}(x) \propto [x(1-x)]^{a-1}$ cambia cualitativamente
en $a = 1$ ([F.B4](repaso/fisica_estadistica.md#fb4)):

- **$a > 1$ (régimen unimodal):** El exponente $a - 1 > 0$ hace que
  $p_{\text{est}}(0) = p_{\text{est}}(1) = 0$ y que la densidad alcance un
  máximo único en $x = 1/2$. La decisión autónoma domina sobre la imitación.
- **$a = 1$ (umbral uniforme):** $p_{\text{est}}(x) = 1$ en $[0, 1]$.
- **$a < 1$ (régimen bimodal de herding):** El exponente $a - 1 < 0$ hace que
  $p_{\text{est}}(x) \to \infty$ cuando $x \to 0^+$ o $x \to 1^-$ (integrable
  para todo $a > 0$). El sistema pasa largos periodos cerca de $x \approx 0$ o
  $x \approx 1$ y realiza transiciones abruptas entre ambos consensos.

**Comparación con el modelo de Ising ([sesión 13](s13_logit_ising.md) y [sesión 14](s14_campo_medio.md), [F.I3](repaso/fisica_estadistica.md#fi3)–[F.I4](repaso/fisica_estadistica.md#fi4)):**

| Propiedad | Ising en campo medio (sesión 14) | Modelo de Kirman |
|---|---|---|
| Regla individual | No lineal: probabilidad logit $\sigma(2\beta J m)$ ([C.B4](repaso/calculo.md#cb4)) | Lineal en la fracción $x = k/N$ |
| Deriva determinista $A(x)$ | Tiene tres ceros cuando $\beta J > 1$ (doble pozo) | Lineal: $A(x) = a(1/2 - x)$, siempre apunta a $x = 1/2$ |
| Intensidad de ruido $B(x)^2$ | Constante de orden $1/N$ lejos de los bordes | Multiplicativa: $B(x)^2 = x(1-x)$, se anula en $x = 0, 1$ |
| Origen de la bimodalidad | Bifurcación determinista de la deriva (ruptura de simetría) | **Bimodalidad inducida por el ruido**: el ruido se apaga en los bordes |
| Tiempo de cruce entre modas | Exponencial en $N$: $\sim e^{N\Delta f}$ (barrera de energía libre) | Polinomial en $N$: las olas de opinión continúan aun con $N$ grande |

En otras palabras, en el modelo de Kirman la tendencia determinista $A(x)$
*siempre* empuja hacia el centro $x = 1/2$; la bimodalidad cuando $a < 1$ no
proviene de un potencial de doble pozo determinista, sino de que la varianza de
la imitación $B(x)^2 = x(1-x)$ se anula en $x = 0$ y $x = 1$. Cuando casi todos
los agentes coinciden, casi no hay pares discrepantes y el sistema queda
atrapado temporalmente cerca de la frontera.

---

## 3. Resumen formal

| Modelo / Concepto | Resultado clave |
|---|---|
| Modelo de votante (red completa) | $p_k = q_k = \frac{k(N-k)}{N(N-1)}$; $x_n = k_n/N$ es martingala; $\mathbb P(\text{gana } 1) = x_0$ |
| Tiempo al consenso (red completa) | $\tau(x_0) = -(N-1)[x_0\ln x_0 + (1-x_0)\ln(1-x_0)]$ barridos; $\tau(1/2) = (N-1)\ln 2$ |
| Votante en red arbitraria | $\bar s_w = (\sum_i k_i s_i)/(\sum_i k_i)$ es martingala; $\langle\tau\rangle \sim N^2$ en 1D, $N\ln N$ en 2D, $N/\ln N$ en BA |
| Modelo de Kirman (discreto) | $W(k\to k+1) = \frac{N-k}{N}\frac{a+k}{N-1+2a}$; estacionaria exacta $\mathrm{BetaBinomial}(N, a, a)$ |
| Modelo de Kirman (continuo) | Fokker-Planck con $A(x) = a(1/2 - x)$ y $B(x)^2 = x(1-x)$; estacionaria $\mathrm{Beta}(a,a) \propto [x(1-x)]^{a-1}$ |
| Umbral de herding en Kirman | Bimodal si $a < 1$, exactamente uniforme si $a = 1$, unimodal si $a > 1$ |

---

## 4. Ejercicios formales

**Ruta básica**

- **F1.** En una red completa de $N = 50$ agentes, $15$ empiezan con la opinión
  $1$ y $35$ con la opinión $0$. ¿Cuál es la probabilidad exacta de que toda la
  red termine en la opinión $1$? ¿Cuántos barridos en promedio tardará en
  alcanzarse un consenso (cualquiera de los dos)?
- **F2.** En una red estrella con un nodo central conectado a $N - 1$ hojas (y
  sin aristas entre las hojas), el nodo central empieza con opinión $1$ y todas
  las hojas empiezan con opinión $0$. Usa la martingala ponderada por el grado
  $\bar s_w$ para calcular la probabilidad exacta de que toda la red termine en
  la opinión $1$.
- **F3.** En el modelo de Kirman con $N = 2$ agentes, escribe explícitamente
  las probabilidades de transición entre los tres estados $k \in
  \lbrace 0, 1, 2\rbrace$ y verifica que $P_0 = P_1 = P_2 = 1/3$ cuando
  $a = 1$.

**Ruta avanzada**

- **F4.** Demuestra que para $N = 3$ agentes en la red completa, resolviendo el
  sistema lineal exacto $\tau_1 = \frac13 + p_1 \tau_2 + (1 - 2p_1)\tau_1$ (por
  simetría $\tau_1 = \tau_2$, con $\tau_0 = \tau_3 = 0$), el tiempo esperado
  exacto hasta el consenso desde $k_0 = 1$ es $\tau_1 = 1$ barrido ($3$ pasos
  elementales), y compáralo con la aproximación continua
  $-(N-1)[x\ln x + (1-x)\ln(1-x)]$ evaluada en $x = 1/3$.
- **F5.** Usando que la distribución estacionaria continua del modelo de Kirman
  es $\mathrm{Beta}(a, a)$ con densidad $p(x) = \frac{1}{\mathrm{B}(a,a)}
  x^{a-1}(1-x)^{a-1}$, demuestra que $\mathbb E[x] = 1/2$ y que la varianza es
  $\operatorname{Var}(x) = \frac{1}{4(2a + 1)}$. Usa la identidad
  $\mathrm{B}(u+1, v) = \frac{u}{u+v}\mathrm{B}(u, v)$ y verifica los límites
  $a \to 0$, $a = 1$ y $a \to \infty$.
- **F6.** Aplica la fórmula de Itô ([E.A2](repaso/procesos_estocasticos.md#ea2))
  con el cambio de variable $y = \arcsin(2x - 1) \in [-\pi/2, \pi/2]$ a la EDE
  del modelo de Kirman $dx = a(1/2 - x)\thinspace dt + \sqrt{x(1-x)}\thinspace dW_t$ y encuentra
  la difusión con ruido aditivo constante que satisface $y_t$. ¿Cómo se
  interpreta el signo de la deriva efectiva de $y_t$ cerca de las fronteras
  $y = \pm\pi/2$ según sea $a < 1$ o $a > 1$?

---

## Referencias

- Alfarano, S., Lux, T., & Wagner, F. (2005). Estimation of agent-based models:
  The case of an asymmetric herding model. *Computational Economics*, 26(1),
  19–49. https://doi.org/10.1007/s10614-005-6415-1
- Becker, G. S. (1991). A note on restaurant pricing and other examples of
  social influences on price. *Journal of Political Economy*, 99(5),
  1109–1116. https://doi.org/10.1086/261791
- Castellano, C., Fortunato, S., & Loreto, V. (2009). Statistical physics of
  social dynamics. *Reviews of Modern Physics*, 81(2), 591–646.
  https://doi.org/10.1103/RevModPhys.81.591
- Clifford, P., & Sudbury, A. (1973). A model for spatial conflict.
  *Biometrika*, 60(3), 581–588. https://doi.org/10.1093/biomet/60.3.581
- Föllmer, H. (1974). Random economies with many interacting agents. *Journal
  of Mathematical Economics*, 1(1), 51–62.
  https://doi.org/10.1016/0304-4068(74)90035-4
- Holley, R. A., & Liggett, T. M. (1975). Ergodic theorems for weakly
  interacting infinite systems and the voter model. *The Annals of
  Probability*, 3(4), 643–663. https://doi.org/10.1214/aop/1176996306
- Kirman, A. (1993). Ants, rationality, and recruitment. *Quarterly Journal of
  Economics*, 108(1), 137–156. https://doi.org/10.2307/2118498
- Kirman, A. (2011). *Complex Economics: Individual and Collective
  Rationality*. Routledge.
- Sen, P., & Chakrabarti, B. K. (2014). *Sociophysics: An Introduction*.
  Oxford University Press.
- Slanina, F. (2014). *Essentials of Econophysics Modelling*. Oxford
  University Press.
- Sood, V., & Redner, S. (2005). Voter model on heterogeneous graphs.
  *Physical Review Letters*, 94(17), 178701.
  https://doi.org/10.1103/PhysRevLett.94.178701
- Suchecki, K., Eguíluz, V. M., & San Miguel, M. (2005). Voter model dynamics
  in complex networks: Role of dimensionality, disorder, and degree
  distribution. *Physical Review E*, 72(3), 036132.
  https://doi.org/10.1103/PhysRevE.72.036132
