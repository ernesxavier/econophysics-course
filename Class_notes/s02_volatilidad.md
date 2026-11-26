# Sesión 2: Volatilidad y modelos nulos. Notas formales

> **Cómo usar estas notas.** La clase combina el cuaderno
> (`Code/s02_volatilidad.jl`) y el pizarrón. El **plan de clase** de abajo
> dice qué se proyecta, qué secciones de estas notas se escriben en el
> pizarrón, cuándo se trabaja en parejas y qué queda para leer. En el
> pizarrón van las definiciones, las proposiciones y los resultados en
> recuadro; las demostraciones, sólo donde el plan lo indica. El texto que las
> rodea se explica de palabra. Como todo está aquí, en clase no hace falta
> copiar. El cuerpo de las notas es la **ruta básica**. Las demostraciones y
> los temas de la ruta avanzada, opcionales, están al final, en **Para
> profundizar**, y no se ven en clase; donde estaban queda una línea que los
> enlaza. La función de autocorrelación y la hipótesis de mercado eficiente
> están en las [notas de la sesión 1](s01_retornos.md).
>
> Las secciones marcadas con ★ son el **núcleo** que se ve en clase; las
> demás se leen por cuenta propia y se evalúan igual.
>
> Repaso recomendado: [P.B3–B4](repaso/probabilidad.md#pb3) (histogramas y la
> normal), [S.B5](repaso/estadistica.md#sb5) (permutaciones),
> [P.I1](repaso/probabilidad.md#pi1) (independencia y correlación),
> [P.I3](repaso/probabilidad.md#pi3) (teorema del límite central),
> [E.I1](repaso/procesos_estocasticos.md#ei1) (autocorrelación),
> [S.I4](repaso/estadistica.md#si4) (datos dependientes); para la ruta
> avanzada, [P.A1](repaso/probabilidad.md#pa1) (cumulantes),
> [F.A2](repaso/fisica_estadistica.md#fa2) (universalidad) y
> [S.A4](repaso/estadistica.md#sa4) (memoria larga).
>
> Lecturas: Cont (2001), secs. 1–3. Complementaria: Mantegna & Stanley (1999),
> caps. 5 y 7; Requião da Cunha (2021), sec. 2.4.

**Plan de clase** (120 minutos: 105 de núcleo y 15 de margen). El núcleo:
incorrelación no es independencia, los hechos 3 y 4 y los modelos nulos.

| Minutos | Dónde | Qué | Secciones |
|---|---|---|---|
| 0–8 | Proyector | El fenómeno: los días malos llegan juntos | Cuaderno, sección 0 |
| 8–15 | Pizarra | Incorrelación no es independencia: el ejemplo de $X$ y $X^2$ | 1 |
| 15–30 | Pizarra | Hecho 3: el enunciado; el modelo $r_t = \sigma_t\varepsilon_t$ y su proposición, con la demostración; la volatilidad fluctuante genera colas pesadas | 2 |
| 30–40 | Pizarra | Hecho 4: el teorema del límite central y $\kappa_n = \kappa_1/n$ (enunciado); por qué los datos quedan por encima | 3 |
| 40–45 | Pizarra | Los dos modelos nulos y su tabla | 4 |
| 45–90 | Parejas | Laboratorio: autocorrelación de $\lvert r\rvert$, agregación y modelos nulos; ejercicio de código 1 | Cuaderno, secciones 2 a 5 y ejercicio 1 |
| 90–105 | Pizarra | Cierre: puesta en común de las preguntas 3 a 5 y las relaciones entre los cuatro hechos | Cuaderno, preguntas 3–5; notas, 6 |
| 105–120 | Margen | Si alcanza, en este orden: la universalidad (cuaderno, sección 6; notas, 5) y la pregunta 6; el ejercicio de código 2 | |
| Después | Lectura | Lo que no entró del margen; universalidad (5) y Para profundizar | |

Lo del **margen** se ve sólo si sobra tiempo, en el orden indicado; lo
que no se vea pasa a la lectura y se evalúa igual.

**Pregunta de la sesión.** Si la dirección del mercado no tiene memoria, ¿por
qué los días malos llegan juntos?

---

## 1. Incorrelación no es independencia ★

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
no (sección 2).

---

## 2. Hecho 3: agrupamiento de volatilidad ★

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

$$
\boxed{\ \mathrm{Cov}(r_t, r_{t+\tau}) = 0,
\qquad
\mathrm{Cov}(r_t^2, r_{t+\tau}^2) = \mathrm{Cov}(\sigma_t^2, \sigma_{t+\tau}^2)\ }
$$

Este es el ejemplo de la sección 1 en versión financiera: incorrelación sin
independencia. Los modelos GARCH (sesión 5) son una forma concreta de
especificar $\sigma_t$, con la diferencia de que $\sigma_t$ depende de los
retornos pasados.

**[Avanzado]** Proposición: en [Para profundizar, P.1](#pp1).

**Moraleja:** una mezcla de normales con varianzas distintas tiene colas más
pesadas que una normal. Los hechos 1 ([sesión 1](s01_retornos.md)) y 3 no son
independientes entre sí: **las fluctuaciones de la volatilidad generan colas
pesadas**.

## 3. Hecho 4: normalidad agregacional ★

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
\boxed{\ \kappa_n = \frac{\kappa_1}{n}\ }
$$

**[Avanzado]** Demostración: en [Para profundizar, P.2](#pp2).

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

## 4. Modelos nulos ★

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

## 5. Universalidad

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

## 6. Resumen formal

| Hecho | Enunciado | Modelo nulo que lo destruye |
|---|---|---|
| Agrupamiento de volatilidad | $\rho_{\lvert r\rvert}(\tau) \sim \tau^{-\beta}$, $\beta \in [0.2, 0.4]$ | Barajado y caminata gaussiana |
| Normalidad agregacional | $\kappa_{\Delta t} \to 0$, pero más lento que $\kappa_1/\Delta t$ | — |

Relaciones con los hechos de la sesión 1:

- Mercado eficiente (martingala) $\Rightarrow$ hecho 2, pero no restringe la
  varianza condicional, así que permite el hecho 3.
- Volatilidad fluctuante $\Rightarrow$ colas pesadas (hecho 3 $\Rightarrow$
  contribuye al hecho 1).
- Dependencia en la volatilidad $\Rightarrow$ convergencia lenta a la normal
  (hecho 3 $\Rightarrow$ hecho 4 más lento).

---

## 7. Ejercicios formales

**Ruta básica**

- **F1.** Explica por qué barajar los retornos conserva exactamente la
  curtosis muestral $\hat\kappa$, pero destruye la autocorrelación de
  $\lvert r\rvert$.
- **F2.** Si los retornos diarios fueran i.i.d. con curtosis en exceso
  $\kappa_1 = 9$, ¿qué curtosis tendrían los retornos de 20 días? Compáralo con
  la curtosis a 20 días que mediste en el cuaderno y explica la diferencia.
- **F3.** En el modelo $r_t = \sigma_t\varepsilon_t$, supón que $\sigma_t$
  alterna entre periodos calmados y turbulentos que duran semanas. Sin hacer
  cuentas, ¿qué signo tiene $\mathrm{Cov}(\lvert r_t\rvert, \lvert r_{t+1}\rvert)$?
  ¿Y $\mathrm{Cov}(r_t, r_{t+1})$?

**Ruta avanzada**

- **F4.** Sea $r = \sigma\varepsilon$ con $\varepsilon \sim \mathcal{N}(0,1)$ y
  $\sigma$ igual a $\sigma_1$ o $\sigma_2$ con probabilidad $1/2$ cada uno.
  Calcula $\kappa$ en función de $\sigma_1/\sigma_2$. ¿Cuánto vale con
  $\sigma_2 = 2\sigma_1$?
- **F5.** Generaliza la proposición del hecho 4 a retornos con varianza
  finita pero autocorrelación lineal $\rho(\tau)$: calcula
  $\mathrm{Var}[r_t(\Delta t)]$ y muestra cómo cambia la regla de la raíz
  cuadrada.
- **F6.** Demuestra el valor esperado de $\hat\rho(\tau)$ bajo la permutación
  que se usa en la sección 4. Usa que, si $x_i = r_i - \hat\mu$, entonces
  $\sum_i x_i = 0$.

---

## Para profundizar

Material opcional de la ruta avanzada: demostraciones y temas que no se ven
en clase. Cada apartado indica de qué sección de las notas viene.

<a id="pp1"></a>

### P.1 Proposición (sección 2)

**Proposición.** Si además $\varepsilon_t \sim \mathcal{N}(0,1)$,

$$
\mathbb{E}[r^4]/\mathbb{E}[r^2]^2 = 3\thinspace\frac{\mathbb{E}[\sigma^4]}{\mathbb{E}[\sigma^2]^2} \ \ge\  3,
$$

con igualdad sólo si $\sigma$ es constante.

*Demostración.* $\mathbb{E}[r^4] = \mathbb{E}[\sigma^4]\thinspace\mathbb{E}[\varepsilon^4] = 3\thinspace\mathbb{E}[\sigma^4]$
y $\mathbb{E}[r^2] = \mathbb{E}[\sigma^2]$. La desigualdad es
$\mathbb{E}[\sigma^4] \ge \mathbb{E}[\sigma^2]^2$, es decir,
$\mathrm{Var}[\sigma^2] \ge 0$. $\square$

<a id="pp2"></a>

### P.2 Demostración (sección 3)

*Demostración.* Usamos **cumulantes**: $c_m$ es el coeficiente
de $s^m/m!$ en $\ln\mathbb{E}[e^{sX}]$. Se cumple $c_2 = \sigma^2$ y
$c_4 = \kappa\thinspace\sigma^4$. Como el logaritmo de la función generadora de una
suma de independientes es la suma de los logaritmos, los cumulantes son
aditivos: $c_m(S_n) = n\thinspace c_m(X)$. Entonces

$$
\kappa_n = \frac{c_4(S_n)}{c_2(S_n)^2} = \frac{n\thinspace c_4}{n^2 c_2^2} = \frac{\kappa_1}{n}. \quad\square
$$

---

## Referencias

- Cont, R. (2001). Empirical properties of asset returns: stylized facts and
  statistical issues. *Quantitative Finance*, 1(2), 223–236.
  https://doi.org/10.1080/713665670
- Mantegna, R. N., & Stanley, H. E. (1999). *An Introduction to Econophysics:
  Correlations and Complexity in Finance*. Cambridge University Press.
  https://doi.org/10.1017/CBO9780511755767
- Requião da Cunha, C. (2021). *Introduction to Econophysics: Contemporary
  Approaches with Python Simulations*. CRC Press.
  https://doi.org/10.1201/9781003127956
