# Sesión 5: Langevin, GBM y GARCH. Notas formales

> **Cómo leer estas notas.** Siguen el orden del programa: primero se observa
> en el cuaderno (`Code/s05_gbm_garch.jl`), después se formaliza aquí. El
> cuerpo de las notas es la **ruta básica**. Las demostraciones y los temas de
> la ruta avanzada, opcionales, están al final, en **Para profundizar**; donde
> estaban queda una línea que los enlaza.
>
> Repaso recomendado: [C.B3](repaso/calculo.md#cb3) (la derivada como tasa de
> cambio), [C.I4](repaso/calculo.md#ci4) (ecuaciones diferenciales),
> [E.I1](repaso/procesos_estocasticos.md#ei1) (autocorrelación),
> [E.I4](repaso/procesos_estocasticos.md#ei4) (GARCH),
> [E.I5](repaso/procesos_estocasticos.md#ei5) (Langevin y Euler-Maruyama); para
> la ruta avanzada, [E.A1](repaso/procesos_estocasticos.md#ea1) a
> [E.A4](repaso/procesos_estocasticos.md#ea4) y
> [S.A1](repaso/estadistica.md#sa1) (máxima verosimilitud).
>
> Lecturas: Requião da Cunha (2021), sec. 2.3.5. Complementaria: Slanina
> (2014), sec. 2.2. Ruta avanzada: Mantegna & Stanley (1999), cap. 10;
> Requião da Cunha (2021), secs. 3.3 y 4.3; Gardiner (consulta).

---

## 1. Ecuaciones diferenciales estocásticas

### 1.1 La ecuación de Langevin

Una **ecuación diferencial estocástica** (EDE) es

$$
dX_t = A(X_t)\thinspace dt + B(X_t)\thinspace dW_t,
$$

con $W_t$ un movimiento browniano. Se interpreta como la ecuación integral
$X_t = X_0 + \int_0^t A(X_s)\thinspace ds + \int_0^t B(X_s)\thinspace dW_s$, donde la
segunda integral es la integral de Itô. Langevin (1908) la introdujo para una
partícula browniana: la tendencia es la fricción y el ruido, los choques
moleculares.

### 1.2 El método de Euler-Maruyama

Con paso $\Delta t$ y $Z_k$ i.i.d. $\mathcal{N}(0, 1)$,

$$
X_{k+1} = X_k + A(X_k)\thinspace\Delta t + B(X_k)\sqrt{\Delta t}\thinspace Z_k .
$$

El término de ruido es $B\thinspace\Delta W$, con
$\Delta W = W_{t+\Delta t} - W_t \sim \mathcal{N}(0, \Delta t)$: escala con
$\sqrt{\Delta t}$.

**[Avanzado]** Convergencia: en [Para profundizar, P.1](#pp1).

## 2. Ornstein-Uhlenbeck

**Definición.** $dX_t = -\theta X_t\thinspace dt + \sigma\thinspace dW_t$, con
$\theta, \sigma > 0$ (Uhlenbeck & Ornstein, 1930).

**Proposición.** Con $X_0 = x_0$,

$$
\mathbb{E}[X_t] = x_0 e^{-\theta t}, \qquad
\mathrm{Var}[X_t] = \frac{\sigma^2}{2\theta}\left(1 - e^{-2\theta t}\right),
$$

y $X_t$ es normal. La distribución estacionaria es
$\mathcal{N}\big(0, \sigma^2/(2\theta)\big)$ y la autocorrelación
estacionaria es $e^{-\theta\tau}$.

**[Avanzado]** Demostración: en [Para profundizar, P.2](#pp2).

**Relación con AR(1).** Muestreado cada $\Delta t$, el proceso cumple
exactamente $X_{t+\Delta t} = \phi X_t + \eta_t$ con $\phi = e^{-\theta\Delta t}$
y $\eta_t$ normal independiente: es un AR(1)
([E.I1](repaso/procesos_estocasticos.md#ei1)). El tiempo $1/\theta$ es la
escala característica en la que el proceso olvida su estado.

## 3. El movimiento browniano geométrico

**Definición.** $dP_t = \mu P_t\thinspace dt + \sigma P_t\thinspace dW_t$.

**Proposición.**
$P_t = P_0\exp\big((\mu - \sigma^2/2)t + \sigma W_t\big)$. En consecuencia,
$\ln P_t$ es normal con media $\ln P_0 + (\mu - \sigma^2/2)t$ y varianza
$\sigma^2 t$, el precio promedio es $\mathbb{E}[P_t] = P_0 e^{\mu t}$ y la
mediana es $P_0 e^{(\mu - \sigma^2/2)t}$.

*Demostración.* Por el lema de Itô con $f(P) = \ln P$
([E.A2](repaso/procesos_estocasticos.md#ea2)),
$d\ln P_t = (\mu - \sigma^2/2)\thinspace dt + \sigma\thinspace dW_t$. Integrar da la
fórmula. La media usa $\mathbb{E}[e^{\sigma W_t}] = e^{\sigma^2 t/2}$; la
mediana, que la mediana de $\sigma W_t$ es 0. $\square$

**Consecuencia para los hechos estilizados.** Los retornos logarítmicos
$\ln(P_{t+\Delta t}/P_t)$ son i.i.d.
$\mathcal{N}\big((\mu - \sigma^2/2)\Delta t,\ \sigma^2\Delta t\big)$: no tienen
colas pesadas ni agrupamiento de volatilidad. El modelo sólo reproduce la
ausencia de autocorrelación lineal.

---

## 4. GARCH(1,1)

### 4.1 Definición y propiedades

**Definición** (Bollerslev, 1986). Con $\varepsilon_t$ i.i.d.
$\mathcal{N}(0, 1)$, $\omega > 0$ y $a, b \ge 0$,

$$
r_t = \sigma_t\varepsilon_t, \qquad \sigma_t^2 = \omega + a\thinspace r_{t-1}^2 + b\thinspace\sigma_{t-1}^2 .
$$

**Proposición.** Sea $\mathcal F_{t-1}$ la información hasta $t - 1$.

1. $\mathbb{E}[r_t \mid \mathcal F_{t-1}] = 0$: los retornos son una
   diferencia de martingala y no tienen autocorrelación lineal.
2. $\mathrm{Var}[r_t \mid \mathcal F_{t-1}] = \sigma_t^2$: la varianza
   condicional depende del pasado.
3. Si $a + b < 1$, el proceso es estacionario con varianza incondicional
   $\bar\sigma^2 = \omega/(1 - a - b)$.

*Demostración.* (1) y (2): $\sigma_t$ es conocido en $\mathcal F_{t-1}$ y
$\varepsilon_t$ es independiente del pasado. (3) Tomando esperanzas y usando
$\mathbb{E}[r_{t-1}^2] = \mathbb{E}[\sigma_{t-1}^2]$, la varianza estacionaria
$v$ cumple $v = \omega + (a + b)v$. $\square$

**Persistencia.** Como
$\mathbb{E}[\sigma_{t+h}^2 \mid \mathcal F_t] - \bar\sigma^2 = (a+b)^{h-1}(\sigma_{t+1}^2 - \bar\sigma^2)$,
un choque de volatilidad decae geométricamente con razón $a + b$, y su
**vida media** es $\ln(1/2)/\ln(a+b)$.

### 4.2 Curtosis

**Proposición.** Si $(a+b)^2 + 2a^2 < 1$,

$$
\frac{\mathbb{E}[r_t^4]}{\mathbb{E}[r_t^2]^2} = \frac{3\big(1 - (a+b)^2\big)}{1 - (a+b)^2 - 2a^2} > 3 .
$$

**[Avanzado]** Demostración: en [Para profundizar, P.3](#pp3).

GARCH produce colas más pesadas que la normal aunque $\varepsilon_t$ sea
normal. Si $(a+b)^2 + 2a^2 \ge 1$, el cuarto momento es infinito.

**[Avanzado]** Autocorrelación de $r^2$: en [Para profundizar, P.4](#pp4).

### 4.3 Estimación

Si se fija la varianza de largo plazo en la varianza muestral
$\hat\sigma^2$, es decir, $\omega = \hat\sigma^2(1 - a - b)$ (en inglés,
*variance targeting*), la log-verosimilitud normal es, salvo una constante,

$$
\ell(a, b) = -\frac{1}{2}\sum_{t=1}^{T}\left(\ln\sigma_t^2 + \frac{r_t^2}{\sigma_t^2}\right),
$$

con $\sigma_t^2$ calculada con la recursión. El cuaderno la maximiza en una
cuadrícula de valores de $(a, b)$ ([S.A1](repaso/estadistica.md#sa1)).

### 4.4 Resultados del laboratorio

GARCH(1,1) ajustado a cada activo y comparado con 20 series simuladas del
modelo ajustado, del mismo largo que los datos:

| Activo | $a$ | $b$ | Vida media (días) | ACF de $\lvert r\rvert$, $\tau = 1$: datos / GARCH | ACF de $\lvert r\rvert$, $\tau = 100$: datos / GARCH | $\alpha$ de la cola: datos / GARCH |
|---|---|---|---|---|---|---|
| IPC (BMV) | 0.075 | 0.9125 | 55 | 0.21 / 0.20 | 0.079 / 0.040 | 2.82 / 4.36 |
| USD/MXN | 0.120 | 0.8575 | 30 | 0.28 / 0.27 | 0.046 / 0.012 | 3.10 / 3.84 |
| S&P 500 | 0.125 | 0.8575 | 39 | 0.29 / 0.31 | 0.084 / 0.027 | 3.12 / 3.29 |
| Bitcoin | 0.110 | 0.8475 | 16 | 0.21 / 0.20 | 0.024 / −0.007 | 4.04 / 4.31 |

(El exponente $\alpha$ se estima como en la sesión 3; el de GARCH es la mediana
de las 20 simulaciones.)

**Conclusión.** GARCH reproduce la ausencia de autocorrelación lineal y el
agrupamiento de volatilidad **de corto plazo**. No reproduce la **memoria
larga**: a 100 días, la autocorrelación de $\lvert r\rvert$ del modelo es la
mitad o menos de la observada, porque decae como $(a+b)^\tau$. Sus colas son
pesadas, pero en general menos que las reales.

**Una pista para la sesión 6.** La recursión
$\sigma_{t+1}^2 = \omega + (a\varepsilon_t^2 + b)\thinspace\sigma_t^2$ es un
**proceso de Kesten**: un factor multiplicativo aleatorio más un término
aditivo. Por eso la cola de $\sigma^2$, y con ella la de $r$, es una ley de
potencias, con un exponente que depende de $a$ y $b$.

---

## 5. [Avanzado] Fokker-Planck y Black-Scholes

El desarrollo de esta sección está en [Para profundizar, P.5](#pp5).


---

## 6. Resumen formal

| Proceso | Ecuación | Resultado clave |
|---|---|---|
| Ornstein-Uhlenbeck | $dX = -\theta X\thinspace dt + \sigma\thinspace dW$ | Estacionaria $\mathcal{N}(0, \sigma^2/2\theta)$; AR(1) al muestrear |
| Browniano geométrico | $dP = \mu P\thinspace dt + \sigma P\thinspace dW$ | Mediana crece a $\mu - \sigma^2/2$; retornos i.i.d. normales |
| GARCH(1,1) | $\sigma_t^2 = \omega + a r_{t-1}^2 + b\sigma_{t-1}^2$ | $\bar\sigma^2 = \omega/(1-a-b)$; curtosis > 3; memoria corta |

---

## 7. Ejercicios formales

**Ruta básica**

- **F1.** En Euler-Maruyama con $\Delta t = 0.01$ y $B(x) = 0.5$, ¿qué
  desviación estándar tiene el término de ruido de un paso? ¿Y el de la
  tendencia si $A(x) = -2x$ y $x = 1$?
- **F2.** En un browniano geométrico con $\mu = 0.08$ y $\sigma = 0.3$
  anuales, ¿a qué tasa crece la mediana del precio? ¿En cuántos años se
  duplica el precio típico?
- **F3.** Un GARCH(1,1) tiene $a = 0.1$, $b = 0.85$ y volatilidad diaria de
  largo plazo de 1%. Calcula $\omega$, la vida media de un choque y la
  curtosis.

**Ruta avanzada**

- **F4.** Demuestra que en el estado estacionario de Ornstein-Uhlenbeck,
  $\mathrm{Cov}(X_t, X_{t+\tau}) = \frac{\sigma^2}{2\theta}e^{-\theta\tau}$.
- **F5.** Completa la demostración de la proposición de la sección 5.2.
- **F6.** Para GARCH(1,1), encuentra la región de $(a, b)$ en la que la
  varianza es finita pero el cuarto momento no. ¿Dónde caen los parámetros
  estimados en el laboratorio?

---

## Para profundizar

Material opcional de la ruta avanzada: demostraciones y temas que no se ven
en clase. Cada apartado indica de qué sección de las notas viene.

<a id="pp1"></a>

### P.1 Convergencia (sección 1.2)

**Convergencia.** Bajo condiciones de regularidad sobre $A$ y $B$,
el error en las trayectorias es de orden $\Delta t^{1/2}$ (convergencia
**fuerte**), y el error en promedios $\mathbb{E}[g(X_T)]$, de orden
$\Delta t$ (convergencia **débil**). Para estudiar distribuciones, como en la
sección 5 del cuaderno, importa la segunda.

<a id="pp2"></a>

### P.2 Demostración (sección 2)

*Demostración.* Por el lema de Itô aplicado a
$e^{\theta t}X_t$, $d(e^{\theta t}X_t) = \sigma e^{\theta t}\thinspace dW_t$, así que

$$
X_t = x_0 e^{-\theta t} + \sigma\int_0^t e^{-\theta(t-s)}\thinspace dW_s .
$$

La integral de Itô de una función determinista es normal con media cero y,
por la isometría de Itô, varianza
$\sigma^2\int_0^t e^{-2\theta(t-s)}\thinspace ds = \frac{\sigma^2}{2\theta}(1 - e^{-2\theta t})$. $\square$

<a id="pp3"></a>

### P.3 Demostración (sección 4.2)

*Demostración.* Escribimos
$\sigma_{t+1}^2 = \omega + (a\varepsilon_t^2 + b)\sigma_t^2$. Como
$\varepsilon_t$ es independiente de $\sigma_t$ y
$\mathbb{E}[(a\varepsilon^2 + b)^2] = 3a^2 + 2ab + b^2 = (a+b)^2 + 2a^2$, en el
estado estacionario $m_4 = \mathbb{E}[\sigma^4]$ cumple

$$
m_4 = \omega^2 + 2\omega(a+b)\bar\sigma^2 + \big((a+b)^2 + 2a^2\big)m_4 .
$$

Con $\bar\sigma^2 = \omega/(1-a-b)$, el término independiente es
$\omega^2(1 + a + b)/(1 - a - b)$. Despejando, dividiendo entre $\bar\sigma^4$
y usando $\mathbb{E}[r^4] = 3m_4$, se obtiene la fórmula. $\square$

<a id="pp4"></a>

### P.4 Autocorrelación de $r^2$ (sección 4.2)

**Autocorrelación de $r^2$.** Cuando el cuarto momento existe,
$\rho_{r^2}(\tau) = \rho_1 (a+b)^{\tau-1}$ con
$\rho_1 = a(1 - b^2 - ab)/(1 - b^2 - 2ab)$. Decae exponencialmente.

<a id="pp5"></a>

### P.5 Fokker-Planck y Black-Scholes (sección 5)


#### 5.1 Fokker-Planck

Para $dX = A(X)\thinspace dt + B(X)\thinspace dW$, la densidad $p(x, t)$ cumple

$$
\frac{\partial p}{\partial t} = -\frac{\partial}{\partial x}\big[A(x)\thinspace p\big] + \frac{1}{2}\frac{\partial^2}{\partial x^2}\big[B(x)^2\thinspace p\big]
$$

([E.A4](repaso/procesos_estocasticos.md#ea4)). Se obtiene como la ecuación de difusión de la
[sesión 4](s04_caminata_markov.md): con el desarrollo de Kramers-Moyal de la
ecuación maestra, cortado en el segundo orden.

**Ejemplo: Ornstein-Uhlenbeck.** Con $A = -\theta x$ y $B = \sigma$, una
solución estacionaria con flujo nulo cumple
$\theta x\thinspace p + \frac{\sigma^2}{2}p' = 0$, de donde
$p(x) \propto e^{-\theta x^2/\sigma^2}$: una normal con varianza
$\sigma^2/(2\theta)$, como en la sección 2.

#### 5.2 Black-Scholes como ecuación del calor

Si la acción sigue un browniano geométrico y la tasa libre de riesgo es
$r_f$, el precio $V(S, t)$ de una opción europea con vencimiento $T$ cumple
(Black & Scholes, 1973)

$$
\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2\frac{\partial^2 V}{\partial S^2} + r_f S\frac{\partial V}{\partial S} - r_f V = 0 .
$$

**Proposición.** Con $k = 2r_f/\sigma^2$, el cambio de variables

$$
S = K e^{x}, \qquad t = T - \frac{2\tau}{\sigma^2}, \qquad
V = K\exp\left(-\frac{k-1}{2}x - \frac{(k+1)^2}{4}\tau\right)u(x, \tau)
$$

transforma la ecuación de Black-Scholes en la ecuación del calor
$\partial_\tau u = \partial_x^2 u$.

*Idea de la demostración.* El cambio $S = Ke^x$ elimina los coeficientes
variables ($S\partial_S = \partial_x$); el cambio de tiempo invierte el sentido
(la condición es final, en $T$) y absorbe el factor $\sigma^2/2$; el factor
exponencial elimina los términos de primer orden y el término $-r_f V$. Los
detalles son un cálculo directo con la regla de la cadena. $\square$

La valuación de opciones se reduce así a un problema de difusión: el mismo
de la caminata aleatoria de Bachelier ([sesión 4](s04_caminata_markov.md)).

---

## Referencias

- Black, F., & Scholes, M. (1973). The pricing of options and corporate
  liabilities. *Journal of Political Economy*, 81(3), 637–654.
  https://doi.org/10.1086/260062
- Bollerslev, T. (1986). Generalized autoregressive conditional
  heteroskedasticity. *Journal of Econometrics*, 31(3), 307–327.
  https://doi.org/10.1016/0304-4076(86)90063-1
- Engle, R. F. (1982). Autoregressive conditional heteroscedasticity with
  estimates of the variance of United Kingdom inflation. *Econometrica*,
  50(4), 987–1007. https://doi.org/10.2307/1912773
- Gardiner, C. (2021). *Elements of Stochastic Methods*. AIP Publishing.
- Mantegna, R. N., & Stanley, H. E. (1999). *An Introduction to Econophysics:
  Correlations and Complexity in Finance*. Cambridge University Press.
  https://doi.org/10.1017/CBO9780511755767
- Requião da Cunha, C. (2021). *Introduction to Econophysics: Contemporary
  Approaches with Python Simulations*. CRC Press.
  https://doi.org/10.1201/9781003127956
- Slanina, F. (2014). *Essentials of Econophysics Modelling*. Oxford
  University Press.
- Uhlenbeck, G. E., & Ornstein, L. S. (1930). On the theory of the Brownian
  motion. *Physical Review*, 36(5), 823–841.
  https://doi.org/10.1103/PhysRev.36.823
