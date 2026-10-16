# Sesión 6: Máxima entropía y modelos cinéticos de intercambio. Notas formales

> **Cómo leer estas notas.** Siguen el orden del programa: primero se observa
> en el cuaderno (`Code/s06_intercambio.jl`), después se formaliza aquí. Lo
> que no tiene marca es **ruta básica**. Las secciones y demostraciones
> marcadas **[Avanzado]** son opcionales.
>
> Repaso recomendado: [P.B3](repaso/probabilidad.md#pb3) (histogramas),
> [C.B6](repaso/calculo.md#cb6) (la integral como área),
> [F.B1–B3](repaso/fisica_estadistica.md#fb1) (de lo micro a lo macro,
> entropía, temperatura), [F.I1](repaso/fisica_estadistica.md#fi1)
> (Boltzmann), [F.I6](repaso/fisica_estadistica.md#fi6) (modelos de
> intercambio), [S.I3](repaso/estadistica.md#si3) (comparar datos con una
> distribución); para la ruta avanzada, [C.A1](repaso/calculo.md#ca1)
> (multiplicadores de Lagrange), [C.A4](repaso/calculo.md#ca4) (función
> gamma) y [P.A5](repaso/probabilidad.md#pa5) (desigualdades).
>
> Lecturas: Yakovenko & Rosser (2009); Requião da Cunha, sec. 8.3; Slanina,
> sec. 7.2. Ruta avanzada: Chakrabarti et al. (2013), caps. 2–5.

---

## 1. Medir la desigualdad

### 1.1 Curva de Lorenz

Sea $X \ge 0$ un ingreso con distribución acumulada $F$ y media
$\mu = \mathbb{E}[X] < \infty$.

**Definición.** La curva de Lorenz es

$$
L(p) = \frac{1}{\mu}\int_0^{F^{-1}(p)} x\thinspace dF(x), \qquad 0 \le p \le 1,
$$

la fracción del ingreso total que recibe la fracción $p$ de menores
ingresos. Para una muestra ponderada (como la ENIGH, con factores de
expansión $w_k$) se ordenan los ingresos de menor a mayor y se grafican los
pares $\left(\sum_{k \le i} w_k / \sum_k w_k,\ \sum_{k \le i} w_k x_k / \sum_k w_k x_k\right)$.

**Propiedades.** $L(0) = 0$, $L(1) = 1$, $L$ es creciente y **convexa** (su
pendiente en $p$ es $F^{-1}(p)/\mu$, que crece con $p$), y $L(p) \le p$. La
igualdad $L(p) = p$ sólo ocurre si todos tienen el mismo ingreso.

### 1.2 Índice de Gini

**Definición.** $G = 1 - 2\int_0^1 L(p)\thinspace dp$: el doble del área
entre la diagonal y la curva de Lorenz ([C.B6](repaso/calculo.md#cb6)).

**Proposición.** Si $X$ y $X'$ son independientes con la misma distribución,

$$
G = \frac{\mathbb{E}\lvert X - X'\rvert}{2\mu} .
$$

Es decir, el Gini es la diferencia promedio entre dos personas al azar,
dividida entre el doble del ingreso promedio. De aquí se sigue que $G$ **no
cambia si se multiplica a todos los ingresos por la misma constante**: no
depende de las unidades ni del nivel general de ingreso.

*Demostración* **[Avanzado]**. Usa
$\lvert a - b\rvert = a + b - 2\min(a, b)$, así que
$\mathbb{E}\lvert X - X'\rvert = 2\mu - 2\thinspace\mathbb{E}[\min(X, X')]$. Por otro
lado, integrando por partes,
$\int_0^1 L(p)\thinspace dp = \frac{1}{\mu}\mathbb{E}\left[X(1 - F(X))\right]$,
porque cada $x$ cuenta para todas las $p$ entre $F(x)$ y 1. Para una
distribución continua, $1 - F(X)$ es la probabilidad de que $X'$ sea mayor
que $X$, así que $\mathbb{E}[X(1 - F(X))]$ es la mitad de
$\mathbb{E}[\min(X, X')]$ (el mínimo es $X$ o $X'$, con la misma probabilidad).
Sustituyendo en la definición de $G$ queda el resultado. $\square$

**Ejemplo: la exponencial.** Si $X$ es exponencial con media $T$, la
diferencia $X - X'$ tiene distribución de Laplace con
$\mathbb{E}\lvert X - X'\rvert = T$, así que $G = T/(2T) = 1/2$ **para todo**
$T$.

**En la ENIGH 2024.** Con el ingreso corriente mensual y el factor de
expansión, $G = 0.401$ por hogar y $0.429$ por persona (ingreso per cápita
del hogar, ponderado por número de integrantes). El 50% de hogares con menor
ingreso recibe el 23% del ingreso total.

---

## 2. Intercambio con dinero conservado

### 2.1 El modelo

**Modelo (Dragulescu y Yakovenko, 2000).** $N$ agentes tienen dinero
$m_1, \dots, m_N \ge 0$, con total $M = \sum_i m_i$ fijo. En cada paso se
eligen dos agentes $i \ne j$ al azar y se reparten su dinero:

$$
m_i' = \varepsilon\thinspace(m_i + m_j), \qquad m_j' = (1 - \varepsilon)(m_i + m_j),
$$

con $\varepsilon$ uniforme en $[0, 1]$ e independiente de todo lo demás. El
dinero total se conserva en cada transacción y nadie puede tener dinero
negativo.

**Proposición.** Para $N$ grande, la distribución estacionaria del dinero
de un agente es la exponencial

$$
p(m) = \frac{1}{T}\thinspace e^{-m/T}, \qquad T = \frac{M}{N} .
$$

Damos tres argumentos de distinta naturaleza.

### 2.2 Primer argumento: conteo (entropía)

Divide el dinero en unidades: el agente $i$ tiene $k_i$ unidades y
$\sum_i k_i = K$. Una **distribución** es la lista $n_0, n_1, n_2, \dots$, donde
$n_k$ es el número de agentes con $k$ unidades. El número de asignaciones de
dinero a agentes compatibles con esa distribución es

$$
W = \frac{N!}{n_0!\thinspace n_1!\thinspace n_2! \cdots} .
$$

Si las transacciones exploran al azar todas las asignaciones con el mismo
total, la distribución que veremos es la que tiene **más asignaciones**: la
que maximiza $\ln W$, la entropía ([F.B2](repaso/fisica_estadistica.md#fb2)),
sujeta a $\sum_k n_k = N$ y $\sum_k k\thinspace n_k = K$. Con la aproximación
de Stirling y multiplicadores de Lagrange (sección 2.3) el resultado es
$n_k \propto e^{-\beta k}$: la **distribución de Boltzmann**
([F.I1](repaso/fisica_estadistica.md#fi1)), con el dinero en el papel de la
energía y $T = 1/\beta$ en el de la temperatura.

### 2.3 [Avanzado] Segundo argumento: máxima entropía continua

**Problema.** Maximizar $S[p] = -\int_0^\infty p(m)\ln p(m)\thinspace dm$ sujeto a
$\int_0^\infty p\thinspace dm = 1$ y $\int_0^\infty m\thinspace p\thinspace dm = T$.

**Solución con multiplicadores de Lagrange** ([C.A1](repaso/calculo.md#ca1)).
Al variar $p$ en cada punto $m$,

$$
\frac{\delta}{\delta p(m)}\left[S - \alpha\int p - \beta\int m\thinspace p\right] = -\ln p(m) - 1 - \alpha - \beta m = 0 ,
$$

de donde $p(m) = e^{-1-\alpha} e^{-\beta m}$. La normalización da
$e^{-1-\alpha} = \beta$ y la media da $\beta = 1/T$. La entropía máxima es
$S = 1 + \ln T$.

**Demostración de que es el máximo.** Sea $q(m) = e^{-m/T}/T$ y $p$
cualquier densidad en $[0, \infty)$ con media $T$. Por la desigualdad de
Gibbs, $\int p \ln(p/q) \ge 0$ (es consecuencia de la desigualdad de Jensen,
[P.A5](repaso/probabilidad.md#pa5)), así que

$$
S[p] = -\int p\ln p \le -\int p \ln q = \int p\left(\ln T + \frac{m}{T}\right) dm = \ln T + 1 = S[q] ,
$$

con igualdad sólo si $p = q$. $\square$

En el cuaderno, la entropía de la distribución simulada (con $T = 1$) sube
desde valores negativos hasta acercarse a 1 en unas pocas transacciones por
agente.

### 2.4 [Avanzado] Tercer argumento: la exponencial es estacionaria

**Proposición.** Si $m_i$ y $m_j$ son independientes y exponenciales con
media $T$, después de una transacción $m_i'$ y $m_j'$ también son
independientes y exponenciales con media $T$.

*Demostración.* $s = m_i + m_j$ tiene densidad gamma $s\thinspace e^{-s/T}/T^2$.
Como $m_i' = \varepsilon s$ con $\varepsilon$ uniforme e independiente de $s$,
su densidad es

$$
f(x) = \int_x^\infty \frac{1}{s}\cdot\frac{s\thinspace e^{-s/T}}{T^2}\thinspace ds = \frac{1}{T}\thinspace e^{-x/T} .
$$

La independencia sale de calcular la densidad conjunta de
$(\varepsilon s, (1-\varepsilon)s)$ con el cambio de variables
$(s, \varepsilon) \mapsto (x, y) = (\varepsilon s, (1-\varepsilon)s)$, cuyo
jacobiano vale $s$. Dividiendo la densidad de $(s, \varepsilon)$ entre él, la
densidad conjunta resulta $e^{-(x+y)/T}/T^2$, un producto. $\square$

**Universalidad.** Otras reglas que conservan el dinero, como transferir un
monto fijo si el que paga lo tiene, llevan a la misma exponencial (en el
cuaderno, más despacio). Lo esencial es la conservación y el azar, no el
detalle de la regla. Lo que rompe la exponencial es permitir deudas, cobrar
impuestos o, como veremos ahora, ahorrar.

---

## 3. Ahorro uniforme

### 3.1 El modelo

**Modelo (Chakraborti y Chakrabarti, 2000).** Cada agente guarda una fracción
$\lambda \in [0, 1)$ de su dinero y sólo pone en juego el resto:

$$
m_i' = \lambda m_i + \varepsilon(1-\lambda)(m_i + m_j), \qquad
m_j' = \lambda m_j + (1-\varepsilon)(1-\lambda)(m_i + m_j).
$$

Con $\lambda = 0$ es el modelo de la sección 2.

### 3.2 La varianza estacionaria

Tomemos dinero promedio 1 y supongamos, como aproximación de campo medio,
que $m_i$ y $m_j$ son independientes con la distribución estacionaria. Sea
$\mu_2 = \mathbb{E}[m^2]$.

**Proposición.** En el estado estacionario,

$$
\mathrm{Var}(m) = \frac{1-\lambda}{1+2\lambda} .
$$

*Demostración.* $\mathbb{E}[m_i'] = \lambda + (1-\lambda)\thinspace\mathbb{E}[\varepsilon]\cdot 2 = 1$,
así que la media se conserva. Para el segundo momento usamos
$\mathbb{E}[\varepsilon] = 1/2$, $\mathbb{E}[\varepsilon^2] = 1/3$,
$\mathbb{E}[m_i(m_i + m_j)] = \mu_2 + 1$ y
$\mathbb{E}[(m_i + m_j)^2] = 2\mu_2 + 2$:

$$
\mathbb{E}[m_i'^2] = \lambda^2\mu_2 + \lambda(1-\lambda)(\mu_2 + 1) + \tfrac{2}{3}(1-\lambda)^2(\mu_2 + 1) .
$$

Imponiendo $\mathbb{E}[m_i'^2] = \mu_2$ y agrupando, el coeficiente de $\mu_2$
es $(1-\lambda)(1+2\lambda)/3$ y el término independiente es
$(1-\lambda)(2+\lambda)/3$. Entonces $\mu_2 = (2+\lambda)/(1+2\lambda)$ y
$\mathrm{Var}(m) = \mu_2 - 1 = (1-\lambda)/(1+2\lambda)$. $\square$

Con $\lambda = 0$ la varianza es 1, la de la exponencial con media 1. Al
aumentar $\lambda$ la varianza baja hasta 0: con ahorro total no hay
intercambio.

### 3.3 La aproximación gamma

**Aproximación (Patriarca, Chakraborti y Kaski, 2004).** La distribución
estacionaria es muy cercana a una gamma con media 1,

$$
p(m) = \frac{n^n}{\Gamma(n)}\thinspace m^{n-1} e^{-nm}, \qquad n = 1 + \frac{3\lambda}{1-\lambda} .
$$

La forma $n$ es la que iguala la varianza de la gamma, $1/n$, con la de la
sección 3.2: $n = (1+2\lambda)/(1-\lambda) = 1 + 3\lambda/(1-\lambda)$. No es
exacta (los momentos de orden 3 y mayores difieren un poco), pero en el
cuaderno la forma ajustada por máxima verosimilitud coincide con la fórmula
para todo $\lambda$ entre 0 y 0.9.

**Interpretación.** Una gamma con forma entera $n$ es la distribución de la
energía de un gas ideal en $2n$ dimensiones: el ahorro actúa como si cada
agente repartiera su dinero entre más "grados de libertad". Para $n > 1$ la
densidad vale cero en $m = 0$: casi nadie queda en la ruina.

**Gini de la gamma** **[Avanzado]**. Con la fórmula de la sección 1.2 y la
función gamma ([C.A4](repaso/calculo.md#ca4)),

$$
G(n) = \frac{\Gamma(n + 1/2)}{n\thinspace\Gamma(n)\sqrt{\pi}} ,
$$

que vale $1/2$ para $n = 1$ y decrece hacia 0 cuando $n \to \infty$.

---

## 4. Ahorro heterogéneo y la cola de Pareto

**Modelo (Chatterjee, Chakrabarti y Manna, 2004).** Cada agente tiene su
propio ahorro $\lambda_i$, fijo, tomado de una distribución $\rho(\lambda)$
en $[0, 1)$:

$$
m_i' = \lambda_i m_i + \varepsilon\left[(1-\lambda_i)m_i + (1-\lambda_j)m_j\right],
\qquad
m_j' = \lambda_j m_j + (1-\varepsilon)\left[(1-\lambda_i)m_i + (1-\lambda_j)m_j\right].
$$

**Proposición (campo medio).** El dinero promedio del agente $i$ en el
estado estacionario es

$$
\langle m_i\rangle = \frac{C}{1 - \lambda_i}, \qquad C = \langle (1-\lambda) m\rangle ,
$$

donde $C$ es el promedio sobre todos los agentes de lo que cada uno pone en
juego.

*Demostración.* Cuando $i$ participa en una transacción, en promedio
conserva $\lambda_i m_i$ y recibe la mitad de lo que se pone en juego,
$\frac{1}{2}\left[(1-\lambda_i)m_i + (1-\lambda_j)m_j\right]$. Promediando
sobre su pareja $j$ e imponiendo que el promedio no cambie,
$\langle m_i\rangle = \lambda_i\langle m_i\rangle + \frac{1}{2}(1-\lambda_i)\langle m_i\rangle + \frac{1}{2}C$,
de donde $(1-\lambda_i)\langle m_i\rangle = C$. $\square$

**Consecuencia: la cola.** Los agentes que ahorran casi todo son los ricos.
Si cerca de $\lambda = 1$ la densidad de ahorros se comporta como
$\rho(\lambda) \approx k(1-\lambda)^{a-1}$, entonces

$$
\mathbb{P}\left(\langle m\rangle > x\right) = \mathbb{P}\left(1-\lambda < \frac{C}{x}\right) \approx \frac{k}{a}\left(\frac{C}{x}\right)^{a},
$$

una cola de **Pareto con exponente $\alpha = a$** (Mohanty, 2006). Con
$\lambda_i$ uniforme ($a = 1$) se obtiene la ley de Zipf, el resultado
original de Chatterjee, Chakrabarti y Manna.

En el cuaderno, $\lambda_i = 1 - U_i^{1/a}$ con $U_i$ uniforme, que tiene
exactamente $\rho(\lambda) = a(1-\lambda)^{a-1}$.

**Cuidados.**

1. El argumento es para el dinero **promedio** de cada agente; el dinero en
   un instante fluctúa alrededor de él.
2. Los agentes con $\lambda_i$ cercano a 1 tardan un tiempo del orden de
   $1/(1-\lambda_i)$ en llegar a su promedio. Con $a = 1$, además,
   $\langle 1/(1-\lambda)\rangle$ diverge: con $N$ finito, unos cuantos agentes
   acumulan una parte enorme del total. Por eso las simulaciones con $N$ y
   tiempo finitos dan exponentes estimados algo mayores que $a$.

---

## 5. Comparación con la ENIGH 2024

Datos: 91,414 hogares que representan 38.8 millones; ingreso corriente
mensual dividido entre su media ponderada.

| Estadística | ENIGH (hogares) | Exponencial | Gamma ajustada | Lognormal ajustada |
|---|---|---|---|---|
| Parámetro | | | $n = 1.94$ ($\lambda = 0.24$) | |
| Gini | 0.40 | 0.50 | 0.38 | |
| Distancia KS a los datos | | 0.163 | 0.059 | 0.012 |

La gamma se ajustó por máxima verosimilitud ponderada; la distancia de
Kolmogorov-Smirnov ([S.I3](repaso/estadistica.md#si3)) es la máxima
diferencia entre la distribución acumulada ponderada de los datos y la del
modelo. La cola, estimada con el estimador de Hill ponderado sobre el 10%,
5% y 1% de hogares con mayor ingreso, tiene exponente entre 2.7 y 3.0.

**Lectura.**

1. La exponencial falla en la parte baja: predice demasiados hogares con
   ingreso cercano a cero.
2. Con ahorro uniforme $\lambda \approx 0.24$ el cuerpo y el Gini mejoran
   mucho, pero la cola de una gamma cae exponencialmente.
3. Con ahorro heterogéneo $a \approx 2.5$ la cola se parece a la de la
   encuesta, pero el cuerpo empeora.
4. La lognormal ajusta mejor que todas. Es lo que predice el crecimiento
   multiplicativo de la sesión 4, otro mecanismo distinto.

**Tres advertencias.**

- Los modelos describen **riqueza** que se intercambia y se conserva; la
  ENIGH mide **ingreso**, un flujo. La analogía es útil, no exacta.
- Las encuestas subestiman los ingresos más altos: los hogares más ricos
  casi no aparecen o declaran de menos. La cola real es más pesada que la de
  la ENIGH; para medirla se combinan encuestas con registros fiscales.
- Que un modelo ajuste no prueba su mecanismo
  ([S.A2](repaso/estadistica.md#sa2)).

---

## 6. Resumen formal

| Modelo | Regla | Resultado |
|---|---|---|
| Dragulescu-Yakovenko | $m_i' = \varepsilon(m_i + m_j)$ | Exponencial; máxima entropía con media fija; $G = 1/2$ |
| Chakraborti-Chakrabarti | Ahorro uniforme $\lambda$ | $\mathrm{Var}(m) = (1-\lambda)/(1+2\lambda)$; gamma con $n = 1 + 3\lambda/(1-\lambda)$ |
| Chatterjee-Chakrabarti-Manna | Ahorro $\lambda_i$ con densidad $\approx k(1-\lambda)^{a-1}$ cerca de 1 | $\langle m_i\rangle \propto 1/(1-\lambda_i)$; cola de Pareto con $\alpha = a$ |
| Gini | $1 - 2\int_0^1 L$ | $\mathbb{E}\lvert X - X'\rvert / (2\mu)$; invariante de escala |

---

## 7. Ejercicios formales

**Ruta básica**

- **F1.** Calcula la curva de Lorenz y el Gini de una población donde el 90%
  gana 1 y el 10% gana 10.
- **F2.** Para la exponencial con media $T$, calcula la mediana y la
  fracción de agentes con menos de la mitad del promedio. Compara con la
  ENIGH (mediana/media = 0.76).
- **F3.** ¿Qué valor de $\lambda$ da Gini de 0.3 con la aproximación gamma?
  (Resuélvelo numéricamente.)
- **F4.** En el modelo con ahorro uniforme, ¿cuánto vale la desviación
  estándar del dinero si $\lambda = 0.5$? Compárala con la de la exponencial.

**Ruta avanzada**

- **F5.** Calcula la curva de Lorenz de una Pareto con exponente $\alpha > 1$
  y demuestra que su Gini es $1/(2\alpha - 1)$.
- **F6.** Demuestra que $0 \le G < 1$ para cualquier ingreso con media
  finita. Para una población de $N$ personas donde una sola recibe todo el
  ingreso, calcula $G$ con la fórmula de la sección 1.2 y comprueba que vale
  $(N-1)/N$.
- **F7.** Calcula $\mathbb{E}[m^3]$ estacionario en el modelo con ahorro
  uniforme y compáralo con el de la gamma de Patriarca. ¿Para qué $\lambda$ es
  mayor la diferencia?
- **F8.** En el modelo heterogéneo con $\lambda_i$ uniforme en $[0, \lambda_{\max}]$
  y $\lambda_{\max} < 1$, ¿hay cola de Pareto? Usa la proposición de la
  sección 4.

---

## Referencias

- Chakrabarti, B. K., Chakraborti, A., Chakravarty, S. R., & Chatterjee, A.
  (2013). *Econophysics of Income and Wealth Distributions*. Cambridge
  University Press. https://doi.org/10.1017/CBO9781139004169
- Chakraborti, A., & Chakrabarti, B. K. (2000). Statistical mechanics of
  money: How saving propensity affects its distribution. *The European
  Physical Journal B*, 17(1), 167–170. https://doi.org/10.1007/s100510070173
- Chatterjee, A., Chakrabarti, B. K., & Manna, S. S. (2004). Pareto law in a
  kinetic model of market with random saving propensity. *Physica A*,
  335(1–2), 155–163. https://doi.org/10.1016/j.physa.2003.11.014
- Dragulescu, A., & Yakovenko, V. M. (2000). Statistical mechanics of money.
  *The European Physical Journal B*, 17(4), 723–729.
  https://doi.org/10.1007/s100510070114
- INEGI (2025). *Encuesta Nacional de Ingresos y Gastos de los Hogares
  (ENIGH) 2024, nueva serie: datos abiertos*. Instituto Nacional de
  Estadística y Geografía.
- Mohanty, P. K. (2006). Generic features of the wealth distribution in
  ideal-gas-like markets. *Physical Review E*, 74(1), 011117.
  https://doi.org/10.1103/PhysRevE.74.011117
- Patriarca, M., Chakraborti, A., & Kaski, K. (2004). Statistical model with
  a standard Γ distribution. *Physical Review E*, 70(1), 016104.
  https://doi.org/10.1103/PhysRevE.70.016104
- Requião da Cunha, C. (2021). *Introduction to Econophysics: Contemporary
  Approaches with Python Simulations*. CRC Press.
  https://doi.org/10.1201/9781003127956
- Slanina, F. (2014). *Essentials of Econophysics Modelling*. Oxford
  University Press.
- Yakovenko, V. M., & Rosser, J. B., Jr. (2009). Colloquium: Statistical
  mechanics of money, wealth, and income. *Reviews of Modern Physics*,
  81(4), 1703–1725. https://doi.org/10.1103/RevModPhys.81.1703
