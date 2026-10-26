# Sesión 13: Logit, Boltzmann e Ising. Notas formales

> **Cómo usar estas notas.** La clase combina el cuaderno
> (`Code/s13_logit_ising.jl`) y el pizarrón. El **plan de clase** de abajo
> dice qué se proyecta, qué secciones de estas notas se escriben en el
> pizarrón, cuándo se trabaja en parejas y qué queda para leer. En el
> pizarrón van las definiciones, las proposiciones y los resultados en
> recuadro; las demostraciones, sólo donde el plan lo indica. El texto que las
> rodea se explica de palabra. Como todo está aquí, en clase no hace falta
> copiar. El cuerpo de las notas es la **ruta básica**. Esta sesión no tiene
> desarrollos avanzados aparte: la ruta avanzada son los ejercicios F4 a F6,
> y el campo medio formal está en las notas de la
> [sesión 14](s14_campo_medio.md).
>
> Repaso recomendado: [C.B4](repaso/calculo.md#cb4) (función logística),
> [F.B3](repaso/fisica_estadistica.md#fb3) (temperatura),
> [F.B4](repaso/fisica_estadistica.md#fb4) (transiciones de fase),
> [F.I1](repaso/fisica_estadistica.md#fi1) (Boltzmann),
> [F.I2](repaso/fisica_estadistica.md#fi2) (logit = Boltzmann),
> [F.I3](repaso/fisica_estadistica.md#fi3) (Ising),
> [F.I5](repaso/fisica_estadistica.md#fi5) (Metropolis); para la ruta
> avanzada, [E.A3](repaso/procesos_estocasticos.md#ea3) (ecuación maestra).
>
> Lecturas: Sen & Chakrabarti (2014), cap. 2. Ruta avanzada: Brock & Durlauf
> (2001).

**Plan de clase** (120 minutos)

| Minutos | Dónde | Qué | Secciones |
|---|---|---|---|
| 0–10 | Proyector | El fenómeno: mismos incentivos, resultados distintos | Cuaderno, sección 0 |
| 10–22 | Pizarra | Repaso: la función logística; utilidad aleatoria y la distribución de Gumbel; la fórmula logit de McFadden (enunciado; es el ejercicio 6 del cuaderno) y el caso binario | 1.1 |
| 22–30 | Pizarra | Logit = Boltzmann: el diccionario entre física estadística y elección discreta | 1.2 |
| 30–40 | Pizarra | El modelo de Ising: la energía, la magnetización y la ruptura espontánea de simetría; la solución exacta de la red cuadrada (enunciado) | 2.1–2.2 |
| 40–52 | Pizarra | Metropolis y baño térmico; el baño térmico es la regla logit; balance detallado, con la demostración | 2.3 |
| 52–100 | Parejas | Laboratorio: logit y Boltzmann, Ising y el algoritmo de Metropolis; la transición con el deslizador; ejercicios de código 1 y 2 | Cuaderno, secciones 2 a 4 y ejercicios 1–2 |
| 100–120 | Pizarra | Cierre: puesta en común de las preguntas 3 a 5; la lectura económica de la dinámica logit | Cuaderno, preguntas 3–5; notas, 2.3 |
| Después | Lectura | La fórmula *log-sum* (1.2), el resumen (3) y los ejercicios F4 a F6 | |

**Pregunta de la sesión.** Dos poblaciones con los mismos incentivos, ¿por
qué terminan eligiendo distinto?

---

## 1. El modelo logit

### 1.1 Utilidad aleatoria

Un agente elige una opción de un conjunto finito. La opción $i$ le da
utilidad $U_i = u_i + \varepsilon_i$, donde $u_i$ es la parte que observa el
analista y $\varepsilon_i$ un choque que el analista no observa. El agente
elige la opción con mayor $U_i$.

**Distribución de Gumbel.** Una variable tiene distribución de Gumbel con
escala $1/\beta$ (y posición 0) si su CDF es

$$
F(x) = \exp\left(-e^{-\beta x}\right), \qquad x \in \mathbb R .
$$

Es la distribución límite del máximo de muchas variables con colas de tipo
exponencial (teoría de valores extremos).

**Proposición** (McFadden). Si los $\varepsilon_i$ son independientes con
distribución de Gumbel de escala $1/\beta$,

$$
\boxed{\ \mathbb P(\text{elige } i) = \frac{e^{\beta u_i}}{\sum_j e^{\beta u_j}}\ }
$$

*Demostración.* La opción $i$ gana si $\varepsilon_j < u_i - u_j + \varepsilon_i$
para todo $j \ne i$. Condicionando en $\varepsilon_i = x$, con densidad
$f(x) = \beta e^{-\beta x}F(x)$,

$$
\mathbb P(i) = \int_{-\infty}^{\infty} \beta e^{-\beta x} e^{-e^{-\beta x}}
\prod_{j \ne i} \exp\left(-e^{-\beta(x + u_i - u_j)}\right) dx
= \int_{-\infty}^{\infty} \beta e^{-\beta x}
\exp\left(-e^{-\beta x}\thinspace\frac{S}{e^{\beta u_i}}\right) dx ,
$$

con $S = \sum_j e^{\beta u_j}$ (el término $j = i$ es el primer factor). Con
el cambio $y = e^{-\beta x}$, $dy = -\beta y\thinspace dx$:

$$
\mathbb P(i) = \int_0^\infty \exp\left(-y\thinspace\frac{S}{e^{\beta u_i}}\right) dy
= \frac{e^{\beta u_i}}{S} . \qquad \square
$$

**Caso binario.** Con dos opciones,
$\mathbb P(A) = \sigma\left(\beta(u_A - u_B)\right)$, con
$\sigma(x) = 1/(1 + e^{-x})$ ([C.B4](repaso/calculo.md#cb4)).

### 1.2 Logit = Boltzmann

La distribución de Boltzmann a temperatura $T = 1/\beta$
([F.I1](repaso/fisica_estadistica.md#fi1)) asigna al estado $i$, de energía
$E_i$, probabilidad $e^{-\beta E_i}/Z$. Es la fórmula logit con
$E_i = -u_i$. El diccionario:

| Física estadística | Elección discreta |
|---|---|
| Estado $i$ | Opción $i$ |
| Energía $E_i$ | Menos la utilidad, $-u_i$ |
| Temperatura $T = 1/\beta$ | Escala de los choques no observados |
| Función de partición $Z$ | Denominador $\sum_j e^{\beta u_j}$ |
| $\ln Z/\beta$ | Excedente esperado del consumidor (fórmula *log-sum*), salvo una constante |

La última línea dice que $\mathbb E\left[\max_i U_i\right] = \beta^{-1}\ln\sum_j e^{\beta u_j}$
más la media de una Gumbel de posición 0 ($\gamma_E/\beta$, con
$\gamma_E \approx 0.577$ la constante de Euler-Mascheroni). Es la energía
libre con el signo cambiado, que veremos en la sesión 14.

---

## 2. El modelo de Ising

### 2.1 Definición

$N$ agentes en los nodos de una red eligen $s_i \in \lbrace -1, +1\rbrace$. La
energía de una configuración $s$ es

$$
H(s) = -J\sum_{\langle i,j\rangle} s_is_j - h\sum_i s_i ,
$$

donde $\langle i,j\rangle$ recorre los pares de vecinos, $J > 0$ premia
coincidir y $h$ es un incentivo común ([F.I3](repaso/fisica_estadistica.md#fi3)). En equilibrio, la configuración
aparece con probabilidad de Boltzmann
$\mathbb P(s) = e^{-\beta H(s)}/Z$. La **magnetización**
$m = \frac1N\sum_i s_i$ es la adopción neta.

**Simetría.** Con $h = 0$, $H(s) = H(-s)$: las dos opciones son
equivalentes y $\mathbb E[m] = 0$ para cualquier $N$ finito. La **ruptura
espontánea de simetría** sólo se ve en el límite $N \to \infty$: debajo de
$T_c$, el sistema se queda en uno de los dos consensos y nunca pasa al otro
([F.B4](repaso/fisica_estadistica.md#fb4)).
En una simulación finita, el sistema salta entre ambos de vez en cuando;
por eso se mide $|m|$.

### 2.2 La solución exacta en la red cuadrada

En la red cuadrada infinita con $h = 0$ (Onsager, 1944; Yang, 1952):

$$
T_c = \frac{2J}{\ln\left(1 + \sqrt 2\right)} \approx 2.269\thinspace J,
\qquad
m(T) = \left[1 - \sinh^{-4}(2J/T)\right]^{1/8} \quad (T < T_c),
$$

y $m = 0$ para $T \ge T_c$. Cerca de $T_c$, $m \sim (T_c - T)^{1/8}$; en la
sesión 14 compararemos esta transición exacta con la aproximación de campo
medio.

### 2.3 Dinámicas: Metropolis y baño térmico

Ambas son cadenas de Markov que cambian un agente a la vez
([F.I5](repaso/fisica_estadistica.md#fi5)). Sea
$h_i^{\text{loc}} = J\sum_{j \sim i} s_j + h$ el campo que siente $i$ y
$\Delta E = 2s_ih_i^{\text{loc}}$ el cambio de energía si $i$ cambia de opción.

- **Metropolis**: cambiar con probabilidad $\min\lbrace 1, e^{-\beta\Delta E}\rbrace$.
- **Baño térmico** (Glauber, 1963): elegir $s_i = +1$ con probabilidad
  $\sigma\left(2\beta h_i^{\text{loc}}\right)$, sin importar el valor anterior.
  Es exactamente la regla **logit** ([F.I2](repaso/fisica_estadistica.md#fi2)): la ventaja en utilidad de $+1$ sobre
  $-1$ es $2h_i^{\text{loc}}$.

**Proposición.** Las dos dinámicas cumplen **balance detallado** respecto a
la distribución de Boltzmann: si $s'$ difiere de $s$ sólo en el agente $i$,

$$
\boxed{\ \mathbb P(s)\thinspace W(s \to s') = \mathbb P(s')\thinspace W(s' \to s)\ }
$$

Por lo tanto la distribución de Boltzmann es estacionaria
([E.A3](repaso/procesos_estocasticos.md#ea3)).

*Demostración.* Sea $\Delta E = H(s') - H(s)$, así que
$\mathbb P(s')/\mathbb P(s) = e^{-\beta\Delta E}$. Metropolis: si
$\Delta E > 0$, $W(s \to s') = e^{-\beta\Delta E}$ y $W(s' \to s) = 1$; si
$\Delta E \le 0$, al revés. En ambos casos
$W(s \to s')/W(s' \to s) = e^{-\beta\Delta E}$. Baño térmico: si $s'_i = +1$,
$W(s \to s') = \sigma(2\beta h)$ y $W(s' \to s) = \sigma(-2\beta h)$, con
$h = h_i^{\text{loc}}$, y $\sigma(x)/\sigma(-x) = e^x$ da
$e^{2\beta h} = e^{-\beta\Delta E}$. (Los factores $1/N$ por elegir al agente
$i$ se cancelan.) $\square$

**Interpretación económica.** Una población de agentes que revisan su
elección de vez en cuando con una regla logit, mirando lo que hacen sus
vecinos, termina distribuida según el modelo de Ising. Blume (1993) lo
estudió como dinámica de juegos de coordinación.

---

## 3. Resumen formal

| Concepto | Resultado |
|---|---|
| Logit | Choques de Gumbel con escala $1/\beta$ dan $\mathbb P(i) = e^{\beta u_i}/\sum_j e^{\beta u_j}$ |
| Logit = Boltzmann | $E_i = -u_i$, $T = 1/\beta$, $Z$ = denominador del logit |
| Ising en red cuadrada | $T_c = 2J/\ln(1+\sqrt2)$; $m \sim (T_c - T)^{1/8}$ |
| Dinámicas | Metropolis y baño térmico (logit) cumplen balance detallado |

---

## 4. Ejercicios formales

**Ruta básica**

- **F1.** Un consumidor elige entre tres marcas con utilidades 2, 1 y 0 y
  $\beta = 1$. Calcula las probabilidades logit. ¿Qué pasa si se agrega una
  cuarta marca idéntica a la segunda?
- **F2.** En una red cuadrada con $J = 1$ y $h = 0$, un agente con $s_i = +1$
  tiene tres vecinos en $+1$ y uno en $-1$. Calcula $\Delta E$ si cambia a
  $-1$, la probabilidad de aceptar el cambio en Metropolis y la probabilidad
  de quedar en $-1$ con el baño térmico para $\beta = 0.5$.
- **F3.** En una simulación de un modelo de Ising finito con $h = 0$, ¿por qué
  se mide el promedio de $|m|$ en lugar del promedio de $m$ para detectar la
  transición de fase?

**Ruta avanzada**

- **F4.** Demuestra que si $\varepsilon_A$ y $\varepsilon_B$ son Gumbel
  independientes de escala $1/\beta$, entonces $\varepsilon_A - \varepsilon_B$
  tiene CDF $\sigma(\beta x)$.
- **F5.** Demuestra que $\ln Z/\beta$ es, salvo una constante, la utilidad
  esperada del agente: $\mathbb E\left[\max_i U_i\right]$.
- **F6.** Verifica el balance detallado del baño térmico para un modelo de
  Ising en una red de 3 agentes en triángulo, escribiendo las 8
  configuraciones y sus probabilidades.

---

## Referencias

- Blume, L. E. (1993). The statistical mechanics of strategic interaction.
  *Games and Economic Behavior*, 5(3), 387–424.
  https://doi.org/10.1006/game.1993.1023
- Brock, W. A., & Durlauf, S. N. (2001). Discrete choice with social
  interactions. *Review of Economic Studies*, 68(2), 235–260.
  https://doi.org/10.1111/1467-937X.00168
- Glaeser, E. L., Sacerdote, B., & Scheinkman, J. A. (1996). Crime and
  social interactions. *Quarterly Journal of Economics*, 111(2), 507–548.
  https://doi.org/10.2307/2946686
- Glauber, R. J. (1963). Time-dependent statistics of the Ising model.
  *Journal of Mathematical Physics*, 4(2), 294–307.
  https://doi.org/10.1063/1.1703954
- McFadden, D. (1974). Conditional logit analysis of qualitative choice
  behavior. En P. Zarembka (Ed.), *Frontiers in Econometrics*
  (pp. 105–142). Academic Press.
- Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., &
  Teller, E. (1953). Equation of state calculations by fast computing
  machines. *Journal of Chemical Physics*, 21(6), 1087–1092.
  https://doi.org/10.1063/1.1699114
- Onsager, L. (1944). Crystal statistics. I. A two-dimensional model with
  an order-disorder transition. *Physical Review*, 65(3–4), 117–149.
  https://doi.org/10.1103/PhysRev.65.117
- Sen, P., & Chakrabarti, B. K. (2014). *Sociophysics: An Introduction*.
  Oxford University Press.
- Yang, C. N. (1952). The spontaneous magnetization of a two-dimensional
  Ising model. *Physical Review*, 85(5), 808–816.
  https://doi.org/10.1103/PhysRev.85.808
