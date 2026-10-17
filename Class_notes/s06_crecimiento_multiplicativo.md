# Sesión 6: Crecimiento multiplicativo. Notas formales

> **Cómo leer estas notas.** Siguen el orden del programa: primero se observa
> en el cuaderno (`Code/s06_crecimiento_multiplicativo.jl`), después se
> formaliza aquí. Lo que no tiene marca es **ruta básica**. Las secciones y
> demostraciones marcadas **[Avanzado]** son opcionales.
>
> Repaso recomendado: [C.B1](repaso/calculo.md#cb1) (logaritmos),
> [P.I3](repaso/probabilidad.md#pi3) (TLC), [P.I6](repaso/probabilidad.md#pi6)
> (lognormal), [E.I6](repaso/procesos_estocasticos.md#ei6) (Gibrat y Kesten);
> para la ruta avanzada, [E.A4](repaso/procesos_estocasticos.md#ea4)
> (Fokker-Planck) y [E.A5](repaso/procesos_estocasticos.md#ea5) (exponente de
> Kesten).
>
> Lecturas: Slanina, secs. 7.1 y 7.3; Gabaix (2009), sec. 3. Ruta avanzada:
> Gabaix (1999); Sornette & Cont (1997).

---

## 1. La ley de Gibrat

### 1.1 Definición

**Definición (ley de Gibrat, 1931).** El tamaño $x_t$ de una unidad (empresa,
ciudad) sigue la ley de Gibrat si

$$
x_{t+1} = a_t\thinspace x_t,
$$

con factores de crecimiento $a_t > 0$ i.i.d. e **independientes del tamaño**
$x_t$. Equivalentemente, la tasa de crecimiento $\ln(x_{t+1}/x_t) = \ln a_t$
tiene la misma distribución para unidades grandes y pequeñas: misma media y
misma dispersión.

### 1.2 Consecuencia: la lognormal

**Proposición.** Si $\ln a_t$ tiene media $m$ y varianza $s^2 < \infty$,
entonces

$$
\ln x_t = \ln x_0 + \sum_{k=0}^{t-1} \ln a_k
$$

es una caminata aleatoria con $\mathbb{E}[\ln x_t] = \ln x_0 + mt$ y
$\mathrm{Var}[\ln x_t] = s^2 t$. Para $t$ grande, $x_t$ es aproximadamente
lognormal.

*Demostración.* El logaritmo convierte el producto en suma
([C.B1](repaso/calculo.md#cb1)); los sumandos son i.i.d., así que la media y
la varianza se suman, y el TLC da la normalidad de $\ln x_t$. $\square$

**Dos consecuencias.**

1. **No hay distribución estacionaria.** La varianza de $\ln x_t$ crece sin
   límite: la distribución de tamaños se ensancha para siempre. En cambio,
   la distribución de tamaños de ciudades cambia muy poco entre censos.
2. **La CCDF se curva en log-log.** Una lognormal no es una ley de potencias,
   aunque con dispersión grande lo parezca en un rango limitado (sesión 3,
   sección 3).

**Promedio contra trayectoria típica.** $\mathbb{E}[x_t] = x_0\thinspace\mathbb{E}[a]^t$,
pero la mediana es $x_0 e^{mt}$ con $m = \mathbb{E}[\ln a] \le \ln\mathbb{E}[a]$
por la desigualdad de Jensen ([P.A5](repaso/probabilidad.md#pa5)). El promedio
puede crecer mientras la unidad típica se encoge
([E.I6](repaso/procesos_estocasticos.md#ei6)).

### 1.3 Prueba con los municipios de México

Con la población de los 2 456 municipios presentes en los censos de 2010 y
2020 (INEGI), sea $g_i = \ln(P_{i,2020}/P_{i,2010})$ la tasa de crecimiento.
La ley de Gibrat implica dos predicciones:

| Predicción | Cómo se mide | Resultado |
|---|---|---|
| La dispersión de $g$ no depende del tamaño | Desviación estándar de $g$ por decil de $P_{2010}$ | Entre 10.5% y 17.1%, sin tendencia: **se cumple** aproximadamente |
| La media de $g$ no depende del tamaño | Pendiente de $g$ contra $\ln P_{2010}$ | $\hat b = 0.026$, $R^2 = 0.082$: **no se cumple** del todo |

Los municipios más pequeños crecieron menos, o se despoblaron: es la
urbanización. La desviación de Gibrat es sistemática pero explica poco de la
variación ($R^2$ de 8%). Gibrat es una buena primera aproximación, no una ley
exacta.

---

## 2. Gibrat con una barrera inferior

### 2.1 El modelo

$$
x_{t+1} = \max\left(a_t\thinspace x_t,\ x_{\min}\right), \qquad \ln a_t \sim \mathcal{N}(m, s^2),\quad m < 0 .
$$

La deriva negativa empuja hacia la barrera; el ruido aleja de ella. La
barrera impide que los tamaños se dispersen hacia cero.

### 2.2 El exponente

**Proposición.** En la aproximación continua, la distribución estacionaria es
una Pareto con $x_{\min}$ y exponente

$$
\alpha = -\frac{2m}{s^2} .
$$

*Demostración* **[Avanzado]**. Sea $y = \ln(x/x_{\min}) \ge 0$. En tiempo
continuo, $y$ es una difusión con deriva $m$ y varianza $s^2$ por unidad de
tiempo, reflejada en $y = 0$. Su densidad cumple la ecuación de Fokker-Planck
([E.A4](repaso/procesos_estocasticos.md#ea4))

$$
\frac{\partial p}{\partial t} = -m\frac{\partial p}{\partial y} + \frac{s^2}{2}\frac{\partial^2 p}{\partial y^2}.
$$

En el estado estacionario con flujo nulo en la barrera,
$m\thinspace p - \frac{s^2}{2}p' = 0$, de donde

$$
p(y) = \frac{2\lvert m\rvert}{s^2}\thinspace e^{-2\lvert m\rvert y/s^2}, \qquad y \ge 0 .
$$

Es decir, $\mathbb{P}(Y > y) = e^{-\alpha y}$ con $\alpha = 2\lvert m\rvert/s^2$.
Como $Y = \ln(X/x_{\min})$,
$\mathbb{P}(X > x) = \mathbb{P}\big(Y > \ln(x/x_{\min})\big) = (x/x_{\min})^{-\alpha}$. $\square$

**Idea clave.** Una distribución **exponencial en el logaritmo** es una
**ley de potencias** en la variable original
([P.I6](repaso/probabilidad.md#pi6)). La difusión con deriva hacia una
barrera produce exponenciales; el crecimiento multiplicativo las traslada a
leyes de potencias.

**Interpretación.** Más ruido ($s$ mayor) o menos encogimiento ($\lvert m\rvert$
menor) dan $\alpha$ menor: colas más pesadas.

---

## 3. Procesos de Kesten

**Definición.** Un **proceso de Kesten** es

$$
x_{t+1} = a_t\thinspace x_t + b_t,
$$

con $(a_t, b_t)$ i.i.d., $a_t > 0$ y $b_t > 0$. El término aditivo $b_t$ cumple
el papel de la barrera: impide que $x$ llegue a cero.

**Teorema (Kesten, 1973; Goldie, 1991; enunciado).** Supón que
$\mathbb{E}[\ln a] < 0$, que $\mathbb{P}(a > 1) > 0$ y que la distribución de
$\ln a$ no está concentrada en una red de puntos. Si existe $\alpha > 0$ con

$$
\mathbb{E}\left[a^{\alpha}\right] = 1,
$$

y $\mathbb{E}[a^\alpha \ln^{+} a]$ y $\mathbb{E}[b^\alpha]$ son finitos, entonces
$x_t$ converge a una distribución estacionaria con
$\mathbb{P}(X > x) \sim C x^{-\alpha}$.

**Argumento heurístico.** Para $x$ grande, $b$ es despreciable y $X \approx aX$
en distribución. Si $\mathbb{P}(X > x) \approx Cx^{-\alpha}$, entonces
$\mathbb{P}(aX > x) = \mathbb{E}\big[\mathbb{P}(X > x/a)\big] \approx C x^{-\alpha}\mathbb{E}[a^\alpha]$,
y la cola se conserva sólo si $\mathbb{E}[a^\alpha] = 1$.

**Proposición.** Si $\ln a \sim \mathcal{N}(m, s^2)$ con $m < 0$, el exponente
de Kesten es $\alpha = -2m/s^2$, el mismo que el de la barrera.

*Demostración.* $\mathbb{E}[a^\alpha] = \mathbb{E}[e^{\alpha\ln a}] = e^{\alpha m + \alpha^2 s^2/2}$.
Igualar a 1 da $\alpha m + \alpha^2 s^2/2 = 0$, con solución positiva
$\alpha = -2m/s^2$. $\square$

**Muestras finitas.** La ley de potencias de Kesten es **asintótica**:
aparece lejos de la escala típica $b/(1 - \mathbb{E}[a])$. En el laboratorio
(10 000 unidades, 2000 periodos, exponente estimado con el método de Clauset):

| $\alpha$ teórico | Barrera: $\hat\alpha$ | Kesten: $\hat\alpha$ |
|---|---|---|
| 0.50 | 0.51 | 0.49 |
| 1.00 | 0.99 | 1.00 |
| 1.50 | 1.51 | 1.37 |
| 2.00 | 2.01 | 1.72 |
| 3.00 | 3.08 | 2.46 |

Con la barrera, la Pareto empieza exactamente en $x_{\min}$ y la estimación
es excelente. Con Kesten, la estimación se queda corta cuando $\alpha$ crece:
un exponente medido con datos finitos puede diferir del asintótico.

---

## 4. La ley de Zipf como consecuencia: el argumento de Gabaix

**Planteamiento** (Gabaix, 1999). Midamos cada ciudad como fracción de la
población total, $x_i = S_i/\sum_j S_j$, y supongamos crecimiento de Gibrat con
una barrera $x_{\min}$. Como las fracciones suman 1, el tamaño **promedio**
$\bar x$ es constante.

**Proposición.** Si la distribución estacionaria es una Pareto con mínimo
$x_{\min}$ y media $\bar x$, entonces

$$
\alpha = \frac{1}{1 - x_{\min}/\bar x}.
$$

*Demostración.* La media de una Pareto con $\alpha > 1$ es
$\bar x = \frac{\alpha}{\alpha - 1}x_{\min}$ (sesión 3, sección 1.3).
Despejando $\alpha$ se obtiene la fórmula. $\square$

**Consecuencia.** Cuando la barrera es muy pequeña comparada con el tamaño
promedio ($x_{\min}/\bar x \to 0$), $\alpha \to 1$: la ley de Zipf. Las
localidades de México, con $\hat\alpha = 0.99 \pm 0.02$ (sesión 3), son
consistentes con este argumento.

**Cuidados.**

- Con $\alpha$ cercano a 1, la media de una muestra finita la dominan las
  pocas unidades más grandes y converge muy lentamente; en simulaciones con
  pocas ciudades, el exponente estimado se aleja de la fórmula.
- Que un mecanismo **pueda** producir Zipf no prueba que sea el que opera:
  la sección 1.3 muestra que el crecimiento de los municipios no es
  exactamente proporcional, y otros modelos también generan Zipf.

---

## 5. [Avanzado] GARCH como proceso de Kesten

En GARCH(1,1), con $r_t = \sigma_t\varepsilon_t$,

$$
\sigma_{t+1}^2 = \left(a\varepsilon_t^2 + b\right)\sigma_t^2 + \omega ,
$$

que es un proceso de Kesten con factor $a\varepsilon_t^2 + b$ y término
aditivo $\omega$. Si $\kappa > 0$ resuelve
$\mathbb{E}\big[(a\varepsilon^2 + b)^{\kappa}\big] = 1$, entonces
$\mathbb{P}(\sigma^2 > x) \sim C x^{-\kappa}$, así que $\sigma$ tiene exponente
$2\kappa$. Como $\varepsilon$ es normal (de colas ligeras), el producto
$r = \sigma\varepsilon$ hereda la cola de $\sigma$:

$$
\mathbb{P}(\lvert r\rvert > x) \sim C' x^{-2\kappa}.
$$

Con los parámetros estimados en la sesión 5:

| Activo | $a$ | $b$ | $2\kappa$ (GARCH) | $\alpha$ de los datos |
|---|---|---|---|---|
| IPC (BMV) | 0.075 | 0.9125 | 6.0 | 2.8 |
| USD/MXN | 0.120 | 0.8575 | 4.9 | 3.1 |
| S&P 500 | 0.125 | 0.8575 | 4.2 | 3.1 |
| Bitcoin | 0.110 | 0.8475 | 7.5 | 4.0 |

Las colas de GARCH con innovaciones normales son más ligeras que las reales.
Los exponentes estimados en las simulaciones de la sesión 5 (entre 3.3 y 4.4)
quedaron por debajo de $2\kappa$ por el mismo efecto de muestra finita de la
sección 3.

---

## 6. Resumen formal

| Modelo | Ecuación | Resultado |
|---|---|---|
| Gibrat | $x_{t+1} = a_t x_t$ | $\ln x_t$ caminata aleatoria; lognormal; sin estado estacionario |
| Gibrat con barrera | $x_{t+1} = \max(a_t x_t, x_{\min})$ | Pareto con $\alpha = -2m/s^2$ |
| Kesten | $x_{t+1} = a_t x_t + b_t$ | $\mathbb{P}(X > x) \sim Cx^{-\alpha}$ con $\mathbb{E}[a^\alpha] = 1$ |
| Gabaix | Fracciones del total con barrera | $\alpha = 1/(1 - x_{\min}/\bar x) \to 1$ |
| GARCH | $\sigma_{t+1}^2 = (a\varepsilon_t^2 + b)\sigma_t^2 + \omega$ | Cola de $r$ con exponente $2\kappa$ |

---

## 7. Ejercicios formales

**Ruta básica**

- **F1.** Una empresa crece 30% o se encoge 25% cada año, con igual
  probabilidad. ¿Crece su tamaño promedio? ¿Y su tamaño típico? Calcula
  $\mathbb{E}[a]$ y $\mathbb{E}[\ln a]$.
- **F2.** En el modelo con barrera, ¿qué valores de $m$ y $s$ dan la ley de
  Zipf? Da dos combinaciones distintas.
- **F3.** Con el argumento de Gabaix, ¿qué exponente predices si la barrera
  es la quinta parte del tamaño promedio?

**Ruta avanzada**

- **F4.** Para el modelo con barrera, calcula la media estacionaria de $x$
  en función de $m$, $s$ y $x_{\min}$. ¿Cuándo es infinita?
- **F5.** Para Kesten con $a$ que vale $\lambda$ o $\mu$ con probabilidad $1/2$
  ($0 < \lambda < 1 < \mu$), encuentra la ecuación que determina $\alpha$ y
  resuélvela para $\lambda = 0.5$, $\mu = 1.5$.
- **F6.** Demuestra que si $X$ tiene cola de potencias con exponente $\alpha$
  y $\varepsilon$ es independiente con todos sus momentos finitos, entonces
  $\mathbb{P}(\lvert X\varepsilon\rvert > x) \sim \mathbb{E}\lvert\varepsilon\rvert^{\alpha}\thinspace\mathbb{P}(\lvert X\rvert > x)$
  (lema de Breiman), al menos de manera heurística.
- **F7.** Calcula numéricamente $2\kappa$ para un GARCH(1,1) con innovaciones
  $t$ de Student con $\nu$ grados de libertad. ¿Puede el modelo reproducir
  $\alpha \approx 3$?

---

## Referencias

- Gabaix, X. (1999). Zipf's law for cities: An explanation. *The Quarterly
  Journal of Economics*, 114(3), 739–767. https://doi.org/10.1162/003355399556133
- Gabaix, X. (2009). Power laws in economics and finance. *Annual Review of
  Economics*, 1, 255–294. https://doi.org/10.1146/annurev.economics.050708.142940
- Gibrat, R. (1931). *Les inégalités économiques*. Sirey.
- Goldie, C. M. (1991). Implicit renewal theory and tails of solutions of
  random equations. *The Annals of Applied Probability*, 1(1), 126–166.
  https://doi.org/10.1214/aoap/1177005985
- INEGI (2011, 2021). *Censos de Población y Vivienda 2010 y 2020: Principales
  resultados por localidad (ITER)*. Instituto Nacional de Estadística y
  Geografía.
- Kesten, H. (1973). Random difference equations and renewal theory for
  products of random matrices. *Acta Mathematica*, 131, 207–248.
  https://doi.org/10.1007/BF02392040
- Slanina, F. (2014). *Essentials of Econophysics Modelling*. Oxford
  University Press.
- Sornette, D., & Cont, R. (1997). Convergent multiplicative processes
  repelled from zero: Power laws and truncated power laws. *Journal de
  Physique I*, 7(3), 431–444. https://doi.org/10.1051/jp1:1997169
