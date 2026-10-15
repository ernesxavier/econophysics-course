# Repaso: probabilidad

> Cada concepto tiene una clave (por ejemplo **[P.B3](#pb3)**) y dice en qué semanas
> se usa. El índice [README.md](README.md) indica qué repasar antes de cada
> sesión. Los fragmentos de Julia corren en el entorno del curso
> (`julia --project=Code`).
>
> - **Básico**: lo necesario para la ruta básica.
> - **Intermedio**: lo necesario para leer las notas formales.
> - **Avanzado**: lo necesario para la ruta avanzada.

---

## Nivel básico

<a id="pb1"></a>

### P.B1 Variable aleatoria y su distribución
*Semanas 1, 2, 6, 8.*

Una **variable aleatoria** $X$ asigna un número a cada resultado de un
experimento: el retorno de mañana, la riqueza de un agente, el grado de un
nodo.

- **Discreta** (toma valores contables): se describe con la
  **función de masa** $p(k) = \mathbb{P}(X = k)$, con $\sum_k p(k) = 1$.
- **Continua**: se describe con la **densidad** $f(x) \ge 0$, con
  $\int f(x)\thinspace dx = 1$ y $\mathbb{P}(a < X \le b) = \int_a^b f(x)\thinspace dx$.
  Es decir, la probabilidad es el **área bajo la densidad** (repaso de la
  integral como área: [C.B6](calculo.md#cb6)). Ojo: $f(x)$ no es una probabilidad; puede ser
  mayor que 1.

Dos funciones que valen para ambos casos:

- **Distribución acumulada (CDF)**: $F(x) = \mathbb{P}(X \le x)$. Crece de 0
  a 1.
- **Distribución acumulada complementaria (CCDF)**:
  $\bar F(x) = \mathbb{P}(X > x) = 1 - F(x)$. Es la herramienta principal
  para estudiar colas (semana 2).

Para una variable continua, $f(x) = F'(x)$.

<a id="pb2"></a>

### P.B2 Media, varianza y desviación estándar
*Semanas 1 en adelante.*

$$
\mu = \mathbb{E}[X] = \sum_k k\thinspace p(k) \ \ \text{o}\ \  \int x f(x)\thinspace dx,
\qquad
\sigma^2 = \mathrm{Var}[X] = \mathbb{E}[(X-\mu)^2] = \mathbb{E}[X^2] - \mu^2.
$$

La **desviación estándar** $\sigma$ tiene las mismas unidades que $X$; en
finanzas se llama **volatilidad**.

Propiedades que usaremos todo el semestre (para constantes $a, b$):

$$
\mathbb{E}[aX + b] = a\thinspace\mathbb{E}[X] + b,
\qquad
\mathrm{Var}[aX + b] = a^2\thinspace\mathrm{Var}[X],
\qquad
\mathbb{E}[X + Y] = \mathbb{E}[X] + \mathbb{E}[Y] \text{ siempre}.
$$

La media de una suma es la suma de las medias **aunque las variables
dependan entre sí**. Con la varianza no pasa lo mismo: sólo si $X$ e $Y$ son
**independientes** se cumple

$$
\mathrm{Var}[X + Y] = \mathrm{Var}[X] + \mathrm{Var}[Y].
$$

Con $n$ variables independientes de varianza $\sigma^2$ cada una, la varianza
de su suma es $n\sigma^2$ y su desviación estándar $\sigma\sqrt{n}$. El caso
con variables correlacionadas está en [P.I1](#pi1).

<a id="pb3"></a>

### P.B3 Histogramas y densidad empírica
*Semanas 1, 2, 6.*

Un **histograma** divide el eje en intervalos (*bins*) y cuenta cuántos
datos caen en cada uno. Para compararlo con una densidad teórica hay que
**normalizarlo**: dividir cada conteo entre $N \cdot \Delta x$ (número de
datos por ancho del bin), de modo que el área total sea 1.

Cuidados:

- El número de bins cambia el aspecto. Con pocos datos, pocos bins.
- Las colas casi no tienen datos y en escala lineal **se ven vacías**. En
  escala logarítmica del eje vertical aparecen ([sesión 1](../s01_hechos_estilizados.md), hecho 1).
- Para leyes de potencias no conviene el histograma con bins iguales: mejor
  la CCDF empírica ([P.I5](#pi5)).

```julia
using Plots, Distributions, Random
x = rand(Xoshiro(1), Normal(0, 1), 10_000)
histogram(x; normalize = :pdf, bins = 60, yscale = :log10)
plot!(-4:0.01:4, z -> pdf(Normal(0, 1), z); lw = 2)
```

<a id="pb4"></a>

### P.B4 Distribución normal y estandarización
*Semanas 1, 2, 3, 10.*

$$
f(x) = \frac{1}{\sigma\sqrt{2\pi}}\exp\negthinspace\left(-\frac{(x-\mu)^2}{2\sigma^2}\right),
\qquad X \sim \mathcal{N}(\mu, \sigma^2).
$$

- **Estandarizar**: $Z = (X - \mu)/\sigma \sim \mathcal{N}(0, 1)$.
  Su CDF se denota $\Phi$.
- **Regla 68–95–99.7**: $\mathbb{P}(|Z| < 1) \approx 0.68$,
  $\mathbb{P}(|Z| < 2) \approx 0.95$, $\mathbb{P}(|Z| < 3) \approx 0.997$.
- Las colas de la normal decaen como $e^{-x^2/2}$: **más rápido que
  cualquier potencia**. Por eso un evento de $5\sigma$ es casi imposible
  bajo la normal (probabilidad $\approx 6 \times 10^{-7}$).
- La suma de normales independientes es normal, con medias y varianzas
  sumadas.

<a id="pb5"></a>

### P.B5 Muestras, estimadores y ley de los grandes números
*Semanas 1 en adelante.*

Los datos $x_1, \dots, x_N$ son una **muestra**; $\mu$ y $\sigma$ son
propiedades de la distribución que no conocemos. Las estimamos con

$$
\hat\mu = \frac{1}{N}\sum_{i=1}^N x_i,
\qquad
\hat\sigma^2 = \frac{1}{N-1}\sum_{i=1}^N (x_i - \hat\mu)^2.
$$

**Ley de los grandes números.** Si los $x_i$ son independientes con media
$\mu$, entonces $\hat\mu \to \mu$ cuando $N \to \infty$.

**Error estándar.** $\mathrm{Var}[\hat\mu] = \sigma^2/N$, así que el error
típico de la media es $\sigma/\sqrt{N}$. Para reducir el error a la mitad
hacen falta **cuatro veces** más datos. Lo mismo vale para promediar
simulaciones ([E.B3](procesos_estocasticos.md#eb3)).

**Cuidado:** si la varianza es infinita (colas con $\alpha < 2$, ver [P.I5](#pi5)),
la media muestral converge muy lentamente o no converge.

---

## Nivel intermedio

<a id="pi1"></a>

### P.I1 Independencia, covarianza y correlación
*Semanas 1, 3, 5, 6.*

$X$ e $Y$ son **independientes** si
$\mathbb{P}(X \in A, Y \in B) = \mathbb{P}(X \in A)\thinspace\mathbb{P}(Y \in B)$
para todo $A, B$. Equivalentemente, la densidad conjunta se factoriza:
$f(x, y) = f_X(x) f_Y(y)$.

$$
\mathrm{Cov}(X, Y) = \mathbb{E}[(X - \mu_X)(Y - \mu_Y)] = \mathbb{E}[XY] - \mu_X\mu_Y,
\qquad
\rho_{XY} = \frac{\mathrm{Cov}(X, Y)}{\sigma_X\thinspace\sigma_Y} \in [-1, 1].
$$

- Independientes $\Rightarrow$ $\mathrm{Cov} = 0$. **El recíproco es
  falso**: si $X \sim \mathcal{N}(0,1)$ e $Y = X^2$, $\mathrm{Cov}(X, Y) = 0$
  pero $Y$ depende totalmente de $X$. La correlación sólo mide dependencia
  **lineal**.
- **Varianza de una suma**:
  $\mathrm{Var}\big[\sum_i X_i\big] = \sum_i \mathrm{Var}[X_i] + 2\sum_{i<j}\mathrm{Cov}(X_i, X_j)$.
  Si son incorreladas, las varianzas se suman. De aquí sale la regla
  $\sigma\sqrt{\Delta t}$ ([sesión 1](../s01_hechos_estilizados.md)) y la idea de diversificación.
- **Matriz de correlación** de $N$ variables: $C_{ij} = \rho_{X_i X_j}$. Es
  simétrica, con unos en la diagonal (semana 5).

<a id="pi2"></a>

### P.I2 Esperanza condicional
*Semanas 1, 3, 14.*

$\mathbb{E}[X \mid Y]$ es la **mejor predicción** de $X$ cuando se conoce
$Y$ (la que minimiza el error cuadrático medio). Es una variable aleatoria:
depende del valor que tome $Y$.

Propiedades clave:

- **Ley de esperanzas iteradas**: $\mathbb{E}\big[\mathbb{E}[X \mid Y]\big] = \mathbb{E}[X]$.
- **Sacar lo conocido**: $\mathbb{E}[g(Y)\thinspace X \mid Y] = g(Y)\thinspace\mathbb{E}[X \mid Y]$.
- Si $X$ e $Y$ son independientes, $\mathbb{E}[X \mid Y] = \mathbb{E}[X]$:
  conocer $Y$ no ayuda.

En finanzas se condiciona sobre la información disponible hasta hoy,
$\mathcal F_t$. Ejemplo: la **varianza condicional**
$\mathrm{Var}[r_{t+1} \mid \mathcal F_t]$ es la volatilidad que esperamos
mañana dado lo que sabemos hoy; los modelos GARCH la modelan ([E.I4](procesos_estocasticos.md#ei4)).

<a id="pi3"></a>

### P.I3 Teorema del límite central
*Semanas 1, 2, 4.*

**Enunciado.** Si $X_1, \dots, X_n$ son i.i.d. con media $\mu$ y varianza
$\sigma^2 < \infty$, y $S_n = X_1 + \cdots + X_n$, entonces

$$
\frac{S_n - n\mu}{\sigma\sqrt{n}} \xrightarrow{d} \mathcal{N}(0, 1).
$$

En palabras: la suma de muchas contribuciones pequeñas e independientes es
aproximadamente normal, **sin importar** la distribución de cada una.

Las tres hipótesis importan, y la econofísica se interesa justo en cuando
fallan:

| Hipótesis | Si falla… | Dónde lo vemos |
|---|---|---|
| Independencia | Convergencia más lenta o a otra cosa | Sesión 1, normalidad agregacional |
| Varianza finita | El límite es una distribución estable, no normal | Semana 2 ([P.A2](#pa2)) |
| Suma (no producto) | Un producto de factores tiende a la **lognormal** | Semana 4 ([P.I6](#pi6)) |

```julia
using Distributions, Statistics, Random
rng = Xoshiro(1)
S = [sum(rand(rng, Exponential(1.0), 50)) for _ in 1:100_000]  # sumas de 50 exponenciales
z = (S .- 50) ./ sqrt(50)
mean(z), std(z), kurtosis(z)   # ≈ 0, 1 y algo cercano a 6/50 = 0.12 (la curtosis estimada es ruidosa)
```

<a id="pi4"></a>

### P.I4 Catálogo de distribuciones del curso
*Semanas 1–7.*

| Distribución | Densidad o masa | Media | Varianza | Julia | Dónde |
|---|---|---|---|---|---|
| Binomial $(n, p)$ | $\binom{n}{k}p^k(1-p)^{n-k}$ | $np$ | $np(1-p)$ | `Binomial(n, p)` | Caminata aleatoria (sem. 3) |
| Uniforme $(a, b)$ | $\frac{1}{b-a}$ en $[a,b]$ | $\frac{a+b}{2}$ | $\frac{(b-a)^2}{12}$ | `Uniform(a, b)` | Ahorro aleatorio (sem. 6) |
| Normal $(\mu, \sigma^2)$ | ver [P.B4](#pb4) | $\mu$ | $\sigma^2$ | `Normal(μ, σ)` | Referencia general |
| Exponencial $(\lambda)$ | $\lambda e^{-\lambda x}$, $x \ge 0$ | $1/\lambda$ | $1/\lambda^2$ | `Exponential(1/λ)` | Dinero conservado (sem. 6) |
| Gamma $(k, \theta)$ | $\frac{x^{k-1}e^{-x/\theta}}{\Gamma(k)\theta^k}$ | $k\theta$ | $k\theta^2$ | `Gamma(k, θ)` | Modelo con ahorro (sem. 6) |
| Lognormal $(\mu, \sigma)$ | $\frac{1}{x\sigma\sqrt{2\pi}}e^{-(\ln x-\mu)^2/2\sigma^2}$ | $e^{\mu+\sigma^2/2}$ | $(e^{\sigma^2}-1)e^{2\mu+\sigma^2}$ | `LogNormal(μ, σ)` | Gibrat (sem. 4) |
| Pareto $(\alpha, x_m)$ | $\frac{\alpha x_m^\alpha}{x^{\alpha+1}}$, $x \ge x_m$ | $\frac{\alpha x_m}{\alpha-1}$ si $\alpha>1$ | finita sólo si $\alpha>2$ | `Pareto(α, xm)` | Colas, riqueza (sem. 2, 6, 7) |
| $t$ de Student $(\nu)$ | $\propto (1 + x^2/\nu)^{-(\nu+1)/2}$ | $0$ si $\nu>1$ | $\frac{\nu}{\nu-2}$ si $\nu>2$ | `TDist(ν)` | Modelo simple de retornos con colas pesadas |

**Ojo con Julia:** `Exponential(θ)` usa la **escala** $\theta = 1/\lambda$,
no la tasa. `Normal(μ, σ)` recibe la desviación estándar, no la varianza.

Notas:

- La **exponencial** es la distribución de máxima entropía con media fija
  en $[0, \infty)$; por eso aparece cuando se conserva el dinero (sem. 6).
  Su CCDF es $e^{-\lambda x}$: una **recta en escala semilog**.
- La **Pareto** tiene CCDF $(x_m/x)^\alpha$: una **recta en escala log-log**.
- La **$t$ de Student** tiene colas $\sim x^{-\nu}$ (en la CCDF) y se
  vuelve normal cuando $\nu \to \infty$.

<a id="pi5"></a>

### P.I5 Leyes de potencias y colas pesadas
*Semanas 2, 4, 6, 7, 8.*

Una variable tiene **cola de ley de potencias** con exponente $\alpha$ si

$$
\bar F(x) = \mathbb{P}(X > x) \sim C\thinspace x^{-\alpha} \quad (x \to \infty)
\qquad\Longleftrightarrow\qquad
f(x) \sim \alpha C\thinspace x^{-(\alpha+1)}.
$$

- **Exponente de la CCDF vs. de la densidad.** Si la CCDF decae con
  $\alpha$, la densidad decae con $\alpha + 1$. Muchos artículos reportan el
  de la densidad; revisa siempre cuál usan.
- **Recta en log-log**: $\ln \bar F(x) \approx \ln C - \alpha \ln x$. La
  pendiente es $-\alpha$.
- **Momentos**: $\mathbb{E}[X^n]$ es finito sólo si $n < \alpha$. Con
  $\alpha \approx 3$ (retornos) la varianza existe pero la curtosis no; con
  $\alpha \approx 1$ (ciudades, Zipf) ni la varianza existe.
- **Invarianza de escala**: $\bar F(\lambda x)/\bar F(x) \to \lambda^{-\alpha}$,
  independiente de $x$. No hay un "tamaño típico".
- **Cola pesada** se usa en sentido amplio para cualquier cola más lenta que
  exponencial (incluye la lognormal, que no es ley de potencias pero en un
  rango finito puede parecerlo).

```julia
using Distributions, Random, Plots
x = sort(rand(Xoshiro(1), Pareto(2.5, 1.0), 10_000); rev = true)
ccdf_emp = (1:length(x)) ./ length(x)       # P(X ≥ x_(i))
plot(x, ccdf_emp; xscale = :log10, yscale = :log10, seriestype = :scatter, ms = 1)
```

<a id="pi6"></a>

### P.I6 Transformación de variables
*Semanas 4, 6.*

Si $Y = g(X)$ con $g$ creciente y derivable, la densidad de $Y$ es

$$
f_Y(y) = f_X\big(g^{-1}(y)\big)\thinspace\left|\frac{d g^{-1}(y)}{dy}\right|.
$$

**Ejemplo: la lognormal.** Si $\ln X \sim \mathcal{N}(\mu, \sigma^2)$,
entonces $X = e^{\ln X}$ tiene densidad
$f_X(x) = \frac{1}{x}\thinspace\varphi\big(\frac{\ln x - \mu}{\sigma}\big)/\sigma$.
Como $\ln(a_1 a_2 \cdots a_n) = \sum \ln a_i$, el TLC aplicado a los
logaritmos dice que **un producto de muchos factores aleatorios es
aproximadamente lognormal**. Esa es la ley de Gibrat (semana 4).

**Ejemplo: de Pareto a exponencial.** Si $X \sim \text{Pareto}(\alpha, x_m)$,
entonces $\ln(X/x_m) \sim \text{Exponencial}(\alpha)$. Esta transformación
está detrás del estimador de Hill ([P.A3](#pa3)).

---

## Nivel avanzado

<a id="pa1"></a>

### P.A1 Funciones generadoras y cumulantes
*Semanas 1, 2.*

- **Función generadora de momentos**: $M(s) = \mathbb{E}[e^{sX}]$. Sus
  derivadas en 0 dan los momentos: $M^{(n)}(0) = \mathbb{E}[X^n]$. No existe
  si las colas son pesadas.
- **Función característica**: $\varphi(s) = \mathbb{E}[e^{isX}]$. Existe
  siempre. Para la normal, $\varphi(s) = e^{i\mu s - \sigma^2 s^2/2}$.
- **Cumulantes** $c_n$: coeficientes de $\ln M(s) = \sum_n c_n s^n/n!$.
  Los primeros son $c_1 = \mu$, $c_2 = \sigma^2$, $c_3 = \gamma\sigma^3$ y
  $c_4 = \kappa\sigma^4$, con $\gamma$ la asimetría y $\kappa$ la curtosis en exceso. La normal
  tiene $c_n = 0$ para $n \ge 3$.
- **Aditividad**: si $X$ e $Y$ son independientes,
  $M_{X+Y} = M_X M_Y$, así que $c_n(X+Y) = c_n(X) + c_n(Y)$. Consecuencia:
  la curtosis de una suma de $n$ i.i.d. es $\kappa/n$ ([sesión 1](../s01_hechos_estilizados.md)).

<a id="pa2"></a>

### P.A2 Distribuciones estables y TLC generalizado
*Semana 2.*

$X$ es **estable** si la suma de dos copias independientes tiene la misma
distribución que $X$, salvo escala y traslación. La normal es estable; la
de Cauchy ($f(x) = \frac{1}{\pi(1+x^2)}$) también.

Las estables se indexan por $\alpha \in (0, 2]$:

- $\alpha = 2$: normal.
- $\alpha < 2$: colas $\bar F(x) \sim x^{-\alpha}$, varianza infinita. La
  suma de $n$ copias escala como $n^{1/\alpha}$ (más rápido que $\sqrt{n}$).

**TLC generalizado.** Si los $X_i$ son i.i.d. con colas
$\sim x^{-\alpha}$ y $\alpha < 2$, la suma, normalizada con $n^{1/\alpha}$,
converge a una estable de índice $\alpha$. Si $\alpha > 2$, converge a la
normal. Mandelbrot (1963) propuso estables para los precios del algodón;
hoy se sabe que los retornos tienen $\alpha \approx 3$, así que agregan a la
normal (lentamente).

<a id="pa3"></a>

### P.A3 Máxima verosimilitud
*Semana 2.*

Dada una muestra i.i.d. y un modelo con densidad $f(x; \theta)$, el
**estimador de máxima verosimilitud** maximiza

$$
\ell(\theta) = \sum_{i=1}^N \ln f(x_i; \theta).
$$

**Ejemplo: exponente de Pareto.** Con $f(x) = \alpha x_m^\alpha x^{-(\alpha+1)}$
para $x \ge x_m$:

$$
\ell(\alpha) = N\ln\alpha + N\alpha\ln x_m - (\alpha+1)\sum_i \ln x_i,
\qquad
\ell'(\alpha) = 0 \quad\Longrightarrow\quad
\hat\alpha = \frac{N}{\sum_{i=1}^N \ln(x_i/x_m)}.
$$

El error estándar es $\hat\alpha/\sqrt{N}$.

El **estimador de Hill** es este mismo cálculo aplicado sólo a los $k$ datos
más grandes, tomando $x_m$ como el dato número $k+1$. El problema práctico
es elegir $k$ (o $x_{\min}$): el método de Clauset et al. (2009) lo elige
minimizando la distancia de Kolmogorov-Smirnov ([P.A4](#pa4)).

<a id="pa4"></a>

### P.A4 Prueba de Kolmogorov-Smirnov
*Semana 2.*

Compara la CDF empírica $F_N(x) = \frac{1}{N}\sum_{i=1}^N \mathbf{1}(x_i \le x)$ (la fracción de datos menores o iguales a $x$) con una CDF
teórica $F$ mediante

$$
D = \sup_x |F_N(x) - F(x)|.
$$

Bajo la hipótesis de que los datos vienen de $F$, $\sqrt{N}D$ tiene una
distribución conocida que no depende de $F$. Si los parámetros de $F$ se
estimaron con los mismos datos, esa distribución ya no aplica y el
$p$-valor se obtiene por simulación (*bootstrap*), como en Clauset et al.

<a id="pa5"></a>

### P.A5 Dos desigualdades útiles
*Semanas 1, 4, 6.*

- **Jensen.** Si $g$ es convexa, $g(\mathbb{E}[X]) \le \mathbb{E}[g(X)]$.
  Con $g(x) = x^2$ da $\mathrm{Var}[X] \ge 0$; con $g = -\ln$, que también es
  convexa, da $\mathbb{E}[\ln X] \le \ln\mathbb{E}[X]$: la tasa de crecimiento
  típica es menor que la del promedio ([E.I6](procesos_estocasticos.md#ei6)).
- **Chebyshev.** $\mathbb{P}(|X - \mu| \ge k\sigma) \le 1/k^2$ para
  cualquier distribución con varianza finita. Es una cota muy holgada: para
  $k = 5$ da 4%, cuando la normal da $6 \times 10^{-7}$. Los datos
  financieros están en medio.

---

## Autodiagnóstico

Intenta responder antes de abrir la respuesta.

<details>
<summary><b>1 (B).</b> Si <i>X</i> tiene media 2 y varianza 9, ¿cuánto valen la media y la desviación estándar de <i>Y</i> = 3 − 2<i>X</i>?</summary>

$\mathbb{E}[Y] = 3 - 2\cdot 2 = -1$. $\mathrm{Var}[Y] = 4 \cdot 9 = 36$, así
que $\sigma_Y = 6$.
</details>

<details>
<summary><b>2 (B).</b> Una densidad vale 3 en algún punto. ¿Es un error?</summary>

No. La densidad no es una probabilidad; sólo su integral debe valer 1. Por
ejemplo, la uniforme en $[0, 1/3]$ vale 3 en todo su soporte.
</details>

<details>
<summary><b>3 (B).</b> Con 2 500 datos, el error estándar de la media es 0.02. ¿Cuántos datos necesitas para bajarlo a 0.01?</summary>

El error escala como $1/\sqrt{N}$: para dividirlo entre 2 hay que
multiplicar $N$ por 4, es decir, 10 000 datos.
</details>

<details>
<summary><b>4 (I).</b> Dos activos tienen la misma volatilidad σ y correlación ρ. ¿Cuál es la volatilidad de un portafolio 50/50?</summary>

$\mathrm{Var} = \tfrac14\sigma^2 + \tfrac14\sigma^2 + 2\cdot\tfrac14\rho\sigma^2 = \tfrac{1+\rho}{2}\sigma^2$,
así que la volatilidad es $\sigma\sqrt{(1+\rho)/2}$. Con $\rho = 1$ no hay
diversificación; con $\rho = 0$ se reduce por $1/\sqrt2$.
</details>

<details>
<summary><b>5 (I).</b> En una gráfica log-log, la CCDF de los tamaños de ciudad es una recta de pendiente −1. ¿Cuál es el exponente de la densidad? ¿Tiene media finita?</summary>

$\alpha = 1$ en la CCDF, así que la densidad decae como $x^{-2}$. La media
requiere $\alpha > 1$: en el límite $\alpha = 1$ diverge (logarítmicamente).
</details>

<details>
<summary><b>6 (I).</b> ¿Por qué el TLC lleva a la lognormal cuando el crecimiento es multiplicativo?</summary>

Si $x_n = x_0 \prod a_i$, entonces $\ln x_n = \ln x_0 + \sum \ln a_i$. El
TLC aplica a la suma de logaritmos, así que $\ln x_n$ es aproximadamente
normal y $x_n$ es aproximadamente lognormal.
</details>

<details>
<summary><b>7 (A).</b> Muestra que la suma de dos normales independientes es normal usando funciones características.</summary>

$\varphi_{X+Y}(s) = \varphi_X(s)\varphi_Y(s) = e^{i\mu_1 s - \sigma_1^2 s^2/2}e^{i\mu_2 s - \sigma_2^2 s^2/2} = e^{i(\mu_1+\mu_2)s - (\sigma_1^2+\sigma_2^2)s^2/2}$,
que es la función característica de $\mathcal{N}(\mu_1+\mu_2, \sigma_1^2+\sigma_2^2)$.
</details>

<details>
<summary><b>8 (A).</b> Si <i>X</i> ~ Pareto(α, 1), ¿qué distribución tiene ln <i>X</i>? Úsalo para justificar el estimador de máxima verosimilitud de α.</summary>

$\mathbb{P}(\ln X > y) = \mathbb{P}(X > e^y) = e^{-\alpha y}$: exponencial con
tasa $\alpha$ y media $1/\alpha$. El estimador natural de $1/\alpha$ es la
media de los $\ln x_i$, de donde $\hat\alpha = N/\sum \ln x_i$, que coincide
con el de máxima verosimilitud.
</details>

---

## Para profundizar

- Requião da Cunha (2021), cap. 2: probabilidad con orientación a
  econofísica.
- Mantegna & Stanley (1999), caps. 3–4: TLC y distribuciones estables.
- Clauset, Shalizi & Newman (2009): leyes de potencias en datos empíricos.
- Gardiner (2021), cap. 2: tratamiento más formal.
