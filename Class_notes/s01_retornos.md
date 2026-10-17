# Sesiones 1 y 2: Econofísica y hechos estilizados. Notas formales

> **Cómo leer estas notas.** Siguen el orden del programa: primero se observa
> en el cuaderno (`Code/s01_retornos.jl`), después se formaliza
> aquí. Lo que no tiene marca es **ruta básica**. Las secciones y
> demostraciones marcadas **[Avanzado]** son opcionales.
>
> Repaso recomendado: [C.B1](repaso/calculo.md#cb1) (logaritmos),
> [P.B1–B5](repaso/probabilidad.md#pb1) (variables aleatorias, media,
> varianza, normal), [E.B1](repaso/procesos_estocasticos.md#eb1) (procesos
> estocásticos), [S.B1–B3](repaso/estadistica.md#sb1) y
> [S.B5](repaso/estadistica.md#sb5) (estimación, pruebas y permutaciones),
> [P.I1–I3](repaso/probabilidad.md#pi1) (correlación, esperanza condicional,
> TLC), [E.I1–I2](repaso/procesos_estocasticos.md#ei1) (autocorrelación y
> martingalas), [C.I1](repaso/calculo.md#ci1) (Taylor),
> [S.I4](repaso/estadistica.md#si4) (datos dependientes); para la ruta
> avanzada, [P.A1](repaso/probabilidad.md#pa1) (cumulantes),
> [F.A2](repaso/fisica_estadistica.md#fa2) (universalidad) y
> [S.A4](repaso/estadistica.md#sa4) (memoria larga).
>
> Lecturas: Cont (2001), secciones 1–3; Mantegna & Stanley, caps. 1–2, 5 y 7;
> Slanina 1.1–1.2; Requião da Cunha 2.4.

---

## 1. Precios y retornos

### 1.1 Definiciones

Sea $P_t > 0$ el precio de cierre de un activo en el día $t = 0, 1, \dots, N$.

**Definición (retorno simple).**

$$
R_t = \frac{P_t - P_{t-1}}{P_{t-1}}.
$$

**Definición (retorno logarítmico).**

$$
r_t = \ln\frac{P_t}{P_{t-1}} = \ln P_t - \ln P_{t-1} = \ln(1 + R_t).
$$

Ambos son adimensionales: no dependen de la moneda ni del nivel del precio,
lo que permite comparar activos distintos.

**Relación entre ambos.** Por el desarrollo de Taylor de $\ln(1+x)$ alrededor
de $0$,

$$
r_t = R_t - \frac{R_t^2}{2} + \frac{R_t^3}{3} - \cdots,
$$

de modo que $r_t \approx R_t$ cuando $|R_t| \ll 1$. Con retornos diarios
típicos ($|R_t| \sim 0.01$, es decir, 1%) la diferencia es del orden de $10^{-5}$. Además
$r_t \le R_t$ siempre, porque $\ln(1+x) \le x$.

### 1.2 Aditividad temporal

**Proposición.** El retorno logarítmico en un horizonte de $\Delta t$ días es
la suma de los retornos diarios:

$$
r_t(\Delta t) \equiv \ln\frac{P_{t+\Delta t}}{P_t}
= \sum_{s=1}^{\Delta t} r_{t+s}.
$$

*Demostración.* Suma telescópica:
$\sum_{s=1}^{\Delta t} (\ln P_{t+s} - \ln P_{t+s-1}) = \ln P_{t+\Delta t} - \ln P_t$. $\square$

El retorno simple, en cambio, se **compone**:
$1 + R_t(\Delta t) = \prod_{s=1}^{\Delta t}(1 + R_{t+s})$.

Esta es la razón por la que trabajamos con $r_t$: convierte la pregunta "¿cómo
es la distribución de los retornos a $\Delta t$ días?" en una pregunta sobre
**sumas de variables aleatorias**, para la cual tenemos el teorema del límite
central (sección 6).

### 1.3 Escalas de tiempo y volatilidad anualizada

Si los retornos diarios fueran independientes con varianza $\sigma^2$, por la
proposición anterior

$$
\mathrm{Var}[r_t(\Delta t)] = \Delta t\thinspace\sigma^2
\quad\Longrightarrow\quad
\sigma(\Delta t) = \sigma\sqrt{\Delta t}.
$$

Con unos 252 días hábiles al año, la **volatilidad anualizada** se reporta
como $\sigma_{\text{anual}} = \sigma_{\text{diaria}}\sqrt{252}$. La regla de
la raíz cuadrada es exacta sólo bajo independencia (más precisamente, bajo
incorrelación; ver ejercicio F2), y es una de las primeras cosas que la
econofísica pone a prueba.

---

## 2. Los retornos como variable aleatoria

### 2.1 Momentos

Tratamos $r_t$ como realización de una variable aleatoria $X$ con media
$\mu = \mathbb{E}[X]$ y varianza $\sigma^2 = \mathbb{E}[(X-\mu)^2]$.

**Definición (momentos estandarizados).**

$$
\text{asimetría: } \gamma = \mathbb{E}\negthinspace\left[\left(\frac{X-\mu}{\sigma}\right)^3\right],
\qquad
\text{curtosis en exceso: } \kappa = \mathbb{E}\negthinspace\left[\left(\frac{X-\mu}{\sigma}\right)^4\right] - 3.
$$

- $\gamma$ mide la asimetría: $\gamma < 0$ significa que las caídas grandes
  son más frecuentes que las subidas grandes.
- $\kappa$ mide el **peso de las colas** relativo a la normal. Se resta 3
  para que la normal tenga $\kappa = 0$ (ver 2.3). Una distribución con
  $\kappa > 0$ se llama **leptocúrtica**.

Ambas son invariantes ante cambios de escala y de origen ($X \mapsto aX + b$,
$a > 0$), así que no dependen de la volatilidad del activo.

### 2.2 Estimadores muestrales

Con una muestra $r_1, \dots, r_N$:

$$
\hat\mu = \frac{1}{N}\sum_{t=1}^N r_t,
\qquad
\hat\sigma^2 = \frac{1}{N-1}\sum_{t=1}^N (r_t - \hat\mu)^2,
\qquad
\hat\kappa = \frac{1}{N}\sum_{t=1}^N \left(\frac{r_t - \hat\mu}{\hat\sigma}\right)^4 - 3.
$$

El factor $N-1$ hace que $\hat\sigma^2$ sea insesgado. Con $N \sim 6000$ días
la diferencia entre dividir por $N$ o por $N-1$ es despreciable. (Las
funciones de `StatsBase.jl` usan $N$ dentro de $\hat\kappa$.)

**Advertencia.** $\hat\kappa$ es muy sensible a los valores extremos: un solo
día con $z = 10$ aporta $10^4/N$ a la suma. Si la distribución tiene colas
de ley de potencias con exponente $\alpha \le 4$, el cuarto momento **no
existe** y $\hat\kappa$ no converge al crecer $N$. Por eso en la sesión 3
pasaremos a medir las colas directamente.

### 2.3 Estandarización y la normal de referencia

**Definición.** El retorno estandarizado es $z_t = (r_t - \hat\mu)/\hat\sigma$.
Tiene media $0$ y varianza $1$; $z = k$ significa "un movimiento de $k$
desviaciones estándar".

**Proposición.** Si $Z \sim \mathcal{N}(0,1)$, entonces
$\mathbb{E}[Z^4] = 3$, es decir, $\kappa = 0$.

*Demostración* **[Avanzado]**. Sea $\varphi(z) = e^{-z^2/2}/\sqrt{2\pi}$, que
cumple $\varphi'(z) = -z\thinspace\varphi(z)$. Integrando por partes con $u = z^3$,
$dv = z\varphi(z)\thinspace dz$, $v = -\varphi(z)$:

$$
\mathbb{E}[Z^4] = \int_{-\infty}^{\infty} z^3 \cdot z\varphi(z)\thinspace dz
= \big[-z^3\varphi(z)\big]_{-\infty}^{\infty} + 3\int_{-\infty}^{\infty} z^2\varphi(z)\thinspace dz
= 0 + 3\thinspace\mathbb{E}[Z^2] = 3. \quad\square
$$

En general, $\mathbb{E}[Z^{2n}] = (2n-1)!! = 1\cdot 3 \cdot 5 \cdots (2n-1)$.

---

## 3. Series de tiempo: estacionariedad y autocorrelación

### 3.1 Procesos estocásticos

Una serie de retornos $\lbrace r_t\rbrace_{t \in \mathbb{Z}}$ es un **proceso
estocástico**: una familia de variables aleatorias indexada por el tiempo.
Observamos una sola trayectoria, así que para estimar cualquier cosa
necesitamos suponer que las propiedades estadísticas no cambian en el tiempo.

**Definición (estacionariedad débil).** $\lbrace r_t\rbrace$ es débilmente estacionario
si

1. $\mathbb{E}[r_t] = \mu$ no depende de $t$;
2. $\mathrm{Var}[r_t] = \sigma^2 < \infty$ no depende de $t$;
3. $\mathrm{Cov}(r_t, r_{t+\tau}) = C(\tau)$ depende sólo del rezago $\tau$.

Los **precios** claramente no son estacionarios (crecen y su varianza
aumenta con el horizonte); los **retornos** son estacionarios en buena
aproximación. Esta es otra razón para trabajar con $r_t$.

### 3.2 Función de autocorrelación

**Definición.** Para un proceso débilmente estacionario, la función de
autocorrelación (ACF) es

$$
\rho(\tau) = \frac{C(\tau)}{C(0)} = \mathrm{corr}(r_t, r_{t+\tau}),
\qquad \tau = 0, 1, 2, \dots
$$

Cumple $\rho(0) = 1$, $|\rho(\tau)| \le 1$ y $\rho(-\tau) = \rho(\tau)$.

**Estimador muestral.**

$$
\hat\rho(\tau) =
\frac{\sum_{t=1}^{N-\tau} (r_t - \hat\mu)(r_{t+\tau} - \hat\mu)}
     {\sum_{t=1}^{N} (r_t - \hat\mu)^2}.
$$

**Resultado (Bartlett).** Si los $r_t$ son independientes e idénticamente
distribuidos (i.i.d.) con varianza finita, entonces para cada $\tau \ge 1$
fijo,

$$
\sqrt{N}\thinspace\hat\rho(\tau) \xrightarrow{d} \mathcal{N}(0, 1)
\qquad (N \to \infty).
$$

Por eso la franja $\pm 1.96/\sqrt{N}$ del cuaderno es un intervalo de 95%
bajo la hipótesis nula de independencia: una barra fuera de la franja es
evidencia (al 5%) de correlación en ese rezago.

**Dos cuidados.**

- Con 30 rezagos, esperamos que en promedio $0.05 \times 30 = 1.5$ barras
  salgan de la franja **por azar**. Una barra aislada fuera de la franja no
  es, por sí sola, evidencia fuerte.
- La banda supone i.i.d. Con agrupamiento de volatilidad (sección 5) los
  $r_t$ son incorrelados pero **no** independientes, y la varianza verdadera
  de $\hat\rho(\tau)$ es mayor que $1/N$. La franja es, entonces, demasiado
  estrecha: exagera la significancia.

### 3.3 Incorrelación no es independencia

**Definición.** $X$ e $Y$ son **independientes** si
$\mathbb{E}[f(X)\thinspace g(Y)] = \mathbb{E}[f(X)]\thinspace\mathbb{E}[g(Y)]$ para todo par de
funciones (acotadas) $f, g$. Son **incorreladas** si eso se cumple sólo para
$f(x) = x$, $g(y) = y$.

La independencia implica incorrelación, pero no al revés. El ejemplo clásico:
$X \sim \mathcal{N}(0,1)$ y $Y = X^2$ son incorreladas
($\mathrm{Cov}(X, X^2) = \mathbb{E}[X^3] = 0$), aunque $Y$ está completamente
determinada por $X$.

En finanzas, la pregunta relevante es si $\rho_{f(r)}(\tau) = 0$ para
funciones no lineales como $f(r) = |r|$ o $f(r) = r^2$. Los datos dicen que
no (sección 5).

---

## 4. Hipótesis de mercado eficiente

### 4.1 Formulación como martingala

Sea $\mathcal F_t$ la información disponible al cierre del día $t$ (todos
los precios hasta $t$, como mínimo).

**Definición.** El precio (descontado) sigue una **martingala** respecto a
$\lbrace\mathcal F_t\rbrace$ si

$$
\mathbb{E}[P_{t+1} \mid \mathcal F_t] = P_t.
$$

En términos de retornos, la versión que usaremos es que el exceso de retorno
es una **diferencia de martingala**:

$$
\mathbb{E}[r_{t+1} - \mu \mid \mathcal F_t] = 0.
$$

**Hipótesis de mercado eficiente (forma débil).** Los precios incorporan toda
la información contenida en los precios pasados. Formalmente: con la
información pasada no se puede predecir el retorno futuro mejor que con su
media incondicional.

El argumento económico es de arbitraje: si existiera una regla del tipo "si
hoy subió, mañana sube", los agentes comprarían hoy, el precio subiría hoy y
la regularidad desaparecería.

### 4.2 Consecuencia: ausencia de autocorrelación lineal

**Proposición.** Si $\lbrace r_t - \mu\rbrace$ es una diferencia de martingala con
varianza finita, entonces $\rho(\tau) = 0$ para todo $\tau \ge 1$.

*Demostración.* Sea $x_t = r_t - \mu$. Como $x_t$ es conocido en
$\mathcal F_{t+\tau-1}$ y por la ley de esperanzas iteradas,

$$
\mathbb{E}[x_t\thinspace x_{t+\tau}]
= \mathbb{E}\big[\thinspace x_t\thinspace\mathbb{E}[x_{t+\tau} \mid \mathcal F_{t+\tau-1}]\thinspace\big]
= \mathbb{E}[x_t \cdot 0] = 0. \quad\square
$$

**Lo que la proposición NO dice.** La martingala restringe sólo la
**media condicional** del retorno. No dice nada sobre la **varianza
condicional** $\mathbb{E}[x_{t+1}^2 \mid \mathcal F_t]$, que puede depender
del pasado. Por eso la eficiencia es perfectamente compatible con el
agrupamiento de volatilidad: la dirección es impredecible, el tamaño no.

**Observación sobre el rezago 1.** En el cuaderno, $\hat\rho(1) \approx 0.08$
en el IPC y $\approx -0.10$ en el S&P 500. Una correlación $\rho$ explica una
fracción $\rho^2$ de la varianza del día siguiente: aquí, alrededor del 1%.
Para que la regularidad sea explotable, la ganancia esperada debe superar los
costos de transacción. Por eso se habla de eficiencia *en el sentido de que
no hay ganancias sistemáticas*, no de que $\rho(1)$ sea exactamente cero.

---

## 5. Los cuatro hechos estilizados

Un **hecho estilizado** es una propiedad estadística que se repite, con
variaciones menores, en muchos activos, mercados y periodos (Cont, 2001). Es
cualitativa: dice "las colas son más pesadas que las de la normal", no "la
curtosis vale 7.3".

### Hecho 1: colas pesadas

**Enunciado.** La distribución incondicional de los retornos tiene colas más
pesadas que la normal: $\kappa > 0$, y los eventos extremos son órdenes de
magnitud más frecuentes de lo que predice una normal con la misma varianza.

**Cómo se cuantifica.** Si $z \sim \mathcal{N}(0,1)$, la probabilidad de un
movimiento mayor a $k$ desviaciones estándar es

$$
p_k = \mathbb{P}(|Z| > k) = 2\big(1 - \Phi(k)\big),
$$

donde $\Phi$ es la función de distribución normal estándar. En $N$ días, el
número esperado de eventos es $N p_k$ y el **periodo de retorno** (tiempo
promedio entre eventos) es $1/p_k$ días.

| $k$ | $p_k$ | Periodo de retorno (normal) |
|---|---|---|
| 3 | $2.7 \times 10^{-3}$ | ≈ 1.5 años |
| 4 | $6.3 \times 10^{-5}$ | ≈ 63 años |
| 5 | $5.7 \times 10^{-7}$ | ≈ 7 000 años |
| 6 | $2.0 \times 10^{-9}$ | ≈ 2 millones de años |

(Periodos calculados con 252 días hábiles por año.) En los datos del curso,
eventos con $|z| > 5$ ocurren varias veces por década.

**Qué tan pesadas.** Una distribución tiene **cola de ley de potencias** con
exponente $\alpha$ si

$$
\mathbb{P}(|r| > x) \sim C\thinspace x^{-\alpha} \qquad (x \to \infty).
$$

Para los retornos financieros se observa $\alpha$ entre 2 y 5, con valores
típicos cercanos a 3 (la "ley cúbica inversa"). Consecuencia: los momentos
$\mathbb{E}[|r|^n]$ con $n \ge \alpha$ divergen. Si $\alpha \approx 3$, la
varianza existe pero la curtosis no (ver la advertencia de 2.2). La
estimación de $\alpha$ es el tema de la sesión 3.

### Hecho 2: ausencia de autocorrelación lineal

**Enunciado.** $\rho_r(\tau) \approx 0$ para todo $\tau \ge 1$, salvo quizás
una correlación pequeña en $\tau = 1$, en escalas de tiempo diarias. (En
escalas intradía de pocos minutos sí aparecen correlaciones, por efectos de
microestructura.)

**Interpretación.** Es la huella empírica de la hipótesis de mercado
eficiente en su forma débil (proposición 4.2).

### Hecho 3: agrupamiento de volatilidad

**Enunciado.** La autocorrelación de los retornos absolutos es positiva y
decae lentamente:

$$
\rho_{|r|}(\tau) > 0 \quad \text{para } \tau \text{ de días a meses},
\qquad
\rho_{|r|}(\tau) \sim A\thinspace\tau^{-\beta}, \quad \beta \in [0.2, 0.4].
$$

Lo mismo ocurre con $r_t^2$. Un decaimiento como ley de potencias, en lugar
de exponencial $e^{-\tau/\tau_0}$, indica **memoria larga**: no hay una escala
de tiempo característica a partir de la cual la volatilidad "olvide".

**Un modelo mínimo.** Considera

$$
r_t = \sigma_t\thinspace\varepsilon_t,
$$

donde $\varepsilon_t$ son i.i.d. con media 0 y varianza 1, y $\lbrace\sigma_t\rbrace$
es un proceso positivo (la volatilidad) independiente de los $\varepsilon$.

**Proposición.** En este modelo:

1. $\mathrm{Cov}(r_t, r_{t+\tau}) = 0$ para $\tau \ge 1$ (no hay
   autocorrelación lineal).
2. $\mathrm{Cov}(r_t^2, r_{t+\tau}^2) = \mathrm{Cov}(\sigma_t^2, \sigma_{t+\tau}^2)$:
   si la volatilidad es persistente, los cuadrados están correlacionados.

*Demostración.* Por independencia entre $\lbrace\sigma\rbrace$ y $\lbrace\varepsilon\rbrace$ y
entre los $\varepsilon$:

1. $\mathbb{E}[r_t r_{t+\tau}] = \mathbb{E}[\sigma_t\sigma_{t+\tau}]\thinspace\mathbb{E}[\varepsilon_t]\thinspace\mathbb{E}[\varepsilon_{t+\tau}] = 0$,
   y $\mathbb{E}[r_t] = 0$.
2. $\mathbb{E}[r_t^2 r_{t+\tau}^2] = \mathbb{E}[\sigma_t^2\sigma_{t+\tau}^2]\thinspace\mathbb{E}[\varepsilon_t^2]\thinspace\mathbb{E}[\varepsilon_{t+\tau}^2] = \mathbb{E}[\sigma_t^2\sigma_{t+\tau}^2]$
   y $\mathbb{E}[r_t^2] = \mathbb{E}[\sigma_t^2]$. $\square$

Este es el ejemplo de 3.3 en versión financiera: incorrelación sin
independencia. Los modelos GARCH (sesión 5) son una forma concreta de
especificar $\sigma_t$, con la diferencia de que $\sigma_t$ depende de los
retornos pasados.

**Proposición [Avanzado].** Si además $\varepsilon_t \sim \mathcal{N}(0,1)$,

$$
\mathbb{E}[r^4]/\mathbb{E}[r^2]^2 = 3\thinspace\frac{\mathbb{E}[\sigma^4]}{\mathbb{E}[\sigma^2]^2} \ \ge\  3,
$$

con igualdad sólo si $\sigma$ es constante.

*Demostración.* $\mathbb{E}[r^4] = \mathbb{E}[\sigma^4]\thinspace\mathbb{E}[\varepsilon^4] = 3\thinspace\mathbb{E}[\sigma^4]$
y $\mathbb{E}[r^2] = \mathbb{E}[\sigma^2]$. La desigualdad es
$\mathbb{E}[\sigma^4] \ge \mathbb{E}[\sigma^2]^2$, es decir,
$\mathrm{Var}[\sigma^2] \ge 0$. $\square$

**Moraleja:** una mezcla de normales con varianzas distintas tiene colas más
pesadas que una normal. Los hechos 1 y 3 no son independientes entre sí:
**las fluctuaciones de la volatilidad generan colas pesadas**.

### Hecho 4: normalidad agregacional

**Enunciado.** Al aumentar el horizonte $\Delta t$, la distribución de
$r_t(\Delta t)$ se parece cada vez más a la normal. Su curtosis disminuye,
pero más lentamente de lo que predice la independencia.

**El punto de referencia: teorema del límite central.** Si
$r_1, r_2, \dots$ son i.i.d. con media $\mu$ y varianza $\sigma^2 < \infty$,
entonces

$$
\frac{r_t(\Delta t) - \Delta t\thinspace\mu}{\sigma\sqrt{\Delta t}}
\xrightarrow{d} \mathcal{N}(0, 1)
\qquad (\Delta t \to \infty).
$$

**Proposición.** Si los retornos diarios son i.i.d. con curtosis en exceso
$\kappa_1$ (finita), la curtosis en exceso de la suma de $n$ de ellos es

$$
\kappa_n = \frac{\kappa_1}{n}.
$$

*Demostración* **[Avanzado]**. Usamos **cumulantes**: $c_m$ es el coeficiente
de $s^m/m!$ en $\ln\mathbb{E}[e^{sX}]$. Se cumple $c_2 = \sigma^2$ y
$c_4 = \kappa\thinspace\sigma^4$. Como el logaritmo de la función generadora de una
suma de independientes es la suma de los logaritmos, los cumulantes son
aditivos: $c_m(S_n) = n\thinspace c_m(X)$. Entonces

$$
\kappa_n = \frac{c_4(S_n)}{c_2(S_n)^2} = \frac{n\thinspace c_4}{n^2 c_2^2} = \frac{\kappa_1}{n}. \quad\square
$$

**Contraste con los datos.** En el cuaderno, la curtosis empírica a
$\Delta t$ días queda **por encima** de la curva $\kappa_1/\Delta t$. Las
dos hipótesis de la proposición pueden fallar:

- **Independencia.** Con agrupamiento de volatilidad, los días turbulentos
  vienen juntos y la suma "hereda" la turbulencia, así que la curtosis baja
  más despacio.
- **Curtosis finita.** Si $\alpha \le 4$, $\kappa_1$ no existe y la
  proposición no aplica. Si $\alpha < 2$ ni siquiera aplicaría el teorema del
  límite central clásico (sesión 3, **[Avanzado]**: distribuciones estables).

---

## 6. Modelos nulos

Un patrón empírico sólo es interesante si no aparece en una versión de los
datos donde, por construcción, no debería existir.

**Modelo nulo 1: caminata gaussiana.** $\tilde r_t \overset{\text{i.i.d.}}{\sim} \mathcal{N}(\hat\mu, \hat\sigma^2)$.
Conserva la media y la varianza; destruye las colas pesadas y toda estructura
temporal.

**Modelo nulo 2: barajado (permutación).** $\tilde r_t = r_{\pi(t)}$, donde
$\pi$ es una permutación aleatoria uniforme de $\lbrace1, \dots, N\rbrace$. Conserva
**exactamente** la distribución empírica (mismo histograma, misma
$\hat\kappa$) y destruye el orden temporal: bajo la permutación,
$\mathbb{E}[\hat\rho(\tau)] = -(N-\tau)/\big(N(N-1)\big) \approx -1/N \approx 0$ para toda
función de los retornos, incluida $|r|$.

| | Colas pesadas | ACF de $r$ ≈ 0 | ACF de $\lvert r\rvert$ > 0 |
|---|---|---|---|
| Datos reales | Sí | Sí | Sí |
| Barajado | Sí | Sí | **No** |
| Caminata gaussiana | **No** | Sí | **No** |

**Conclusión.** Las colas pesadas son una propiedad de la **distribución
marginal**; el agrupamiento de volatilidad es una propiedad del **orden
temporal**. Un modelo de mercado tiene que reproducir ambas, y el barajado
nos dice que son ingredientes separables.

La idea de comparar contra datos barajados reaparece en la sesión 7: lo que
sobresale del espectro de la matriz de correlación barajada es señal.

---

## 7. Universalidad

Los cuatro activos del curso (IPC, USD/MXN, S&P 500 y BTC) tienen
volatilidades anuales entre 12% y 66% y mecanismos de mercado muy distintos.
Aun así, los cuatro cumplen los cuatro hechos estilizados.

En física estadística, **universalidad** significa que sistemas con detalles
microscópicos distintos comparten el mismo comportamiento macroscópico
(por ejemplo, los mismos exponentes críticos). La hipótesis de trabajo de la
econofísica es análoga: los hechos estilizados no dependen de los detalles
institucionales de cada mercado, sino de la forma en que interactúan muchos
agentes heterogéneos. Si esto es cierto, modelos simples de agentes deberían
reproducirlos. Esa es la pregunta que cerramos en la sesión 20.

---

## 8. Resumen formal

| Hecho | Enunciado | Modelo nulo que lo destruye |
|---|---|---|
| Colas pesadas | $\kappa > 0$; $\mathbb{P}(\lvert r\rvert > x) \sim x^{-\alpha}$, $\alpha \approx 3$ | Caminata gaussiana |
| Sin autocorrelación lineal | $\rho_r(\tau) \approx 0$, $\tau \ge 1$ | (Ninguno: también lo cumplen) |
| Agrupamiento de volatilidad | $\rho_{\lvert r\rvert}(\tau) \sim \tau^{-\beta}$, $\beta \in [0.2, 0.4]$ | Barajado y caminata gaussiana |
| Normalidad agregacional | $\kappa_{\Delta t} \to 0$, pero más lento que $\kappa_1/\Delta t$ | — |

Relaciones entre ellos:

- Mercado eficiente (martingala) $\Rightarrow$ hecho 2, pero no restringe la
  varianza condicional, así que permite el hecho 3.
- Volatilidad fluctuante $\Rightarrow$ colas pesadas (hecho 3 $\Rightarrow$
  contribuye al hecho 1).
- Dependencia en la volatilidad $\Rightarrow$ convergencia lenta a la normal
  (hecho 3 $\Rightarrow$ hecho 4 más lento).

---

## 9. Ejercicios formales

**Ruta básica**

- **F1.** Muestra que si una acción sube 50% un día y baja 50% al siguiente,
  el retorno simple total no es 0, pero la suma de los retornos simples sí.
  Calcula los retornos logarítmicos y verifica la aditividad.
- **F2.** Usando $\mathrm{Var}[\sum X_i] = \sum_i \mathrm{Var}[X_i] + 2\sum_{i<j}\mathrm{Cov}(X_i, X_j)$,
  muestra que la regla $\sigma(\Delta t) = \sigma\sqrt{\Delta t}$ requiere
  sólo incorrelación, no independencia. ¿Qué pasa si $\rho(1) > 0$?
- **F3.** Con la tabla de periodos de retorno: si en 25 años de datos hubo 8
  días con $|z| > 5$, ¿cuántas veces más frecuentes son que lo que predice la
  normal?

**Ruta avanzada**

- **F4.** Demuestra que $\mathbb{E}[Z^6] = 15$ para $Z \sim \mathcal{N}(0,1)$.
- **F5.** Sea $r = \sigma\varepsilon$ con $\varepsilon \sim \mathcal{N}(0,1)$ y
  $\sigma$ igual a $\sigma_1$ o $\sigma_2$ con probabilidad $1/2$ cada uno.
  Calcula $\kappa$ en función de $\sigma_1/\sigma_2$. ¿Cuánto vale con
  $\sigma_2 = 2\sigma_1$?
- **F6.** Sea $\mathbb{P}(X > x) = x^{-\alpha}$ para $x \ge 1$ (Pareto).
  Muestra que $\mathbb{E}[X^n] < \infty$ si y sólo si $n < \alpha$.
- **F7.** Generaliza la proposición del hecho 4 a retornos con varianza
  finita pero autocorrelación lineal $\rho(\tau)$: calcula
  $\mathrm{Var}[r_t(\Delta t)]$ y muestra cómo cambia la regla de la raíz
  cuadrada.

---

## Referencias

- Cont, R. (2001). Empirical properties of asset returns: stylized facts and
  statistical issues. *Quantitative Finance*, 1(2), 223–236.
  https://doi.org/10.1080/713665670
- Mantegna, R. N., & Stanley, H. E. (1999). *An Introduction to Econophysics:
  Correlations and Complexity in Finance*. Cambridge University Press.
  https://doi.org/10.1017/CBO9780511755767
- Slanina, F. (2014). *Essentials of Econophysics Modelling*. Oxford
  University Press.
- Requião da Cunha, C. (2021). *Introduction to Econophysics: Contemporary
  Approaches with Python Simulations*. CRC Press.
  https://doi.org/10.1201/9781003127956
