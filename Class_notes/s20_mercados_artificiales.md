# Sesión 20: Mercados artificiales y expectativas heterogéneas. Notas formales

> **Cómo leer estas notas.** Siguen el orden del programa: primero se observa
> en el cuaderno (`Code/s20_mercados_artificiales.jl`), después se formaliza
> aquí. El cuerpo de las notas es la **ruta básica**. Las demostraciones y los
> temas de la ruta avanzada, opcionales, están al final, en **Para
> profundizar**; donde estaban queda una línea que los enlaza.
>
> Repaso recomendado: [E.B3](repaso/procesos_estocasticos.md#eb3)
> (estacionariedad y promedios temporales),
> [S.B5](repaso/estadistica.md#sb5) (pruebas de permutación),
> [E.I1](repaso/procesos_estocasticos.md#ei1) (autocorrelación),
> [E.I4](repaso/procesos_estocasticos.md#ei4) (ARCH y GARCH),
> [C.I5](repaso/calculo.md#ci5) (puntos fijos y estabilidad en 1D),
> [S.I2](repaso/estadistica.md#si2) (prueba de Ljung-Box),
> [S.I4](repaso/estadistica.md#si4) (datos dependientes); para la ruta
> avanzada, [L.A3](repaso/algebra_lineal.md#la3) (jacobiano y estabilidad),
> [C.A5](repaso/calculo.md#ca5) (mapas discretos, bifurcaciones y caos),
> [S.A3](repaso/estadistica.md#sa3) (calibración por momentos simulados) y
> [S.A4](repaso/estadistica.md#sa4) (memoria larga).
>
> Lecturas: LeBaron (2006). Complementaria: Hommes (2013), cap. 1; Lux &
> Marchesi (1999). Ruta avanzada: Hommes (2013), caps. 2 y 6; Slanina (2014),
> sec. 3.2; Giardina & Bouchaud (2003); Brock & Hommes (1998). De consulta:
> Hommes (2006); Franke & Westerhoff (2012).

---

## 1. De los hechos estilizados a los mercados de agentes heterogéneos

### 1.1 El problema que quedó abierto en las sesiones 1 y 2

En las [sesiones 1](s01_retornos.md) y [2](s02_volatilidad.md) vimos que los retornos logarítmicos
diarios $r_t = \ln(P_t/P_{t-1})$ de cualquier activo financiero líquido cumplen
tres regularidades universales (Cont, 2001):

1. **Ausencia de autocorrelación lineal:** $\hat\rho_r(\tau) \approx 0$ para
   $\tau \ge 1$. Predecir el signo del retorno de mañana con una combinación
   lineal de retornos pasados es casi imposible.
2. **Colas pesadas (leptocurtosis):** la curtosis en exceso
   $\kappa = \mathbb E[z_t^4] - 3$ es muy positiva y la distribución acumulada
   complementaria (CCDF) $\mathbb P(|z_t| > x)$ cae mucho más despacio que la
   normal estándar.
3. **Agrupamiento de volatilidad:** la autocorrelación de la magnitud de los
   retornos, $\hat\rho_{|r|}(\tau) = \mathrm{corr}(|r_t|, |r_{t+\tau}|)$, es
   positiva y decae lentamente durante semanas o meses.

En la [sesión 5](s05_gbm_garch.md) vimos que los procesos GARCH
describen estadísticamente esa varianza cambiante, pero no explican **qué
mecanismo microeconómico** la produce ni por qué la volatilidad estalla incluso
en días sin noticias macroeconómicas relevantes (Cutler, Poterba y Summers,
1989).

### 1.2 Las piezas microscópicas reunidas

Los **mercados financieros artificiales** y los **modelos de expectativas
heterogéneas** (Brock y Hommes, 1998; Lux y Marchesi, 1999; Giardina y Bouchaud,
2003; LeBaron, 2006) cierran el ciclo del curso uniendo las herramientas de los
módulos anteriores:

- **Reglas simples de predicción** en lugar de expectativas racionales únicas:
  conviven **fundamentalistas** (que apuestan a que el precio regresará a su
  valor fundamental $p^\ast$) y **chartistas** o **seguidores de tendencia** (que
  extrapolan las subidas o caídas recientes).
- **Elección discreta logit** ([sesión 13](s13_logit_ising.md)): los agentes cambian de
  regla con probabilidad proporcional a $e^{\beta \pi_h}$, donde $\pi_h$ es la
  ganancia reciente de la regla $h$ y $\beta$ es la intensidad de elección.
- **Contagio e imitación** ([sesión 15](s15_votante_kirman.md)): el optimismo o
  pesimismo de los chartistas se propaga socialmente.
- **Retroalimentación ecológica** ([sesiones 17](s17_replicador_goodwin.md) y
  [19](s19_juego_minoria.md)): el éxito de una estrategia depende de cuántos
  otros agentes la están usando.

---

## 2. El modelo de Brock y Hommes (1998): selección evolutiva de expectativas

### 2.1 Equilibrio temporal de mercado con creencias heterogéneas

Consideremos un mercado con un activo libre de riesgo que paga un factor bruto
constante $R = 1 + r > 1$ por periodo, y un activo riesgoso con precio $p_t$ y
dividendo $y_t$ con media constante $\bar y$.

Si los inversionistas tienen utilidad media-varianza con aversión absoluta al
riesgo $a > 0$ y creen que la varianza condicional del exceso de retorno es
$\sigma^2$, la demanda de acciones de un agente que usa la regla predictiva $h$
en el periodo $t$ es:

$$
z_{h,t} = \frac{\mathbb E_{h,t}[p_{t+1} + y_{t+1}] - R\thinspace p_t}{a\thinspace \sigma^2} .
$$

Con oferta neta cero de acciones ($\sum_{h=1}^H n_{h,t} z_{h,t} = 0$, donde
$n_{h,t}$ es la fracción de agentes que usa la regla $h$ en el periodo $t$), el
precio que vacía el mercado satisface:

$$
R\thinspace p_t = \sum_{h=1}^H n_{h,t}\thinspace \mathbb E_{h,t}[p_{t+1} + y_{t+1}] .
$$

**Valor fundamental y desviaciones.** Si todos los agentes fueran racionales e
idénticos, el único precio estacionario sin burbujas sería el **precio
fundamental** $p^\ast$ tal que $R\thinspace p^\ast = p^\ast + \bar y$, es decir:

$$
p^\ast = \frac{\bar y}{R - 1} = \frac{\bar y}{r} .
$$

Definiendo la **desviación respecto al fundamental** $x_t = p_t - p^\ast$, y
suponiendo que las reglas $h = 1, \dots, H$ difieren únicamente en su pronóstico
$f_{h,t} = \mathbb E_{h,t}[x_{t+1}]$ de la desviación futura, la ecuación de
precios se simplifica a:

$$
\boxed{R\thinspace x_t = \sum_{h=1}^H n_{h,t}\thinspace f_{h,t}}
$$

(o $R\thinspace x_t = \sum_{h=1}^H n_{h,t}\thinspace f_{h,t} + \varepsilon_t$ en presencia de
demanda ruidosa $\varepsilon_t$).

### 2.2 Fundamentalistas contra seguidores de tendencia

El caso canónico de Brock y Hommes (1998) enfrenta dos reglas de predicción:

1. **Fundamentalistas ($h = 1$):** creen que el precio regresará al valor
   fundamental $p^\ast$ ($x = 0$), por lo que pronostican
   $$
   f_{1,t} = 0 .
   $$
   Obtener y procesar la información fundamental tiene un costo fijo $C \ge 0$
   por periodo.
2. **Chartistas / seguidores de tendencia ($h = 2$):** extrapolan la última
   desviación observada con un coeficiente de tendencia $g > 0$:
   $$
   f_{2,t} = g\thinspace x_{t-1} .
   $$
   Mirar el gráfico de precios no tiene costo de información ($C_2 = 0$).

Con estas dos reglas, el precio en el periodo $t$ queda determinado por la
fracción de chartistas $n_{2,t}$:

$$
R\thinspace x_t = n_{2,t}\thinspace g\thinspace x_{t-1} .
$$

### 2.3 Ganancias realizadas y actualización evolutiva logit

¿Cómo deciden los agentes qué regla usar en el periodo $t$? Mirando qué regla
ganó más dinero en el periodo anterior.

El exceso de retorno por acción entre $t-1$ y $t$ es:

$$
R_{e,t} = p_t + y_t - R\thinspace p_{t-1} = x_t - R\thinspace x_{t-1} + \delta_t ,
$$

donde $\delta_t = y_t - \bar y$ es la sorpresa del dividendo. La demanda que
colocó la regla $h$ en $t-1$ (cuando pronosticaba $f_{h,t-1}$ para $x_t$) fue
$z_{h,t-1} = (f_{h,t-1} - R\thinspace x_{t-1}) / (a\sigma^2)$. Reescalando unidades para
hacer $a\sigma^2 = 1$, la **ganancia neta realizada** de la regla $h$ al cierre
del periodo $t-1$ es:

$$
\pi_{h,t-1} = (x_{t-1} - R\thinspace x_{t-2})\thinspace (f_{h,t-2} - R\thinspace x_{t-2}) - C_h .
$$

Para nuestras dos reglas ($f_{1,t-2} = 0$ con costo $C$, y
$f_{2,t-2} = g\thinspace x_{t-3}$ con costo $0$):

$$
\begin{aligned}
\pi_{1,t-1} &= -R\thinspace x_{t-2}\thinspace (x_{t-1} - R\thinspace x_{t-2}) - C , \cr
\pi_{2,t-1} &= (g\thinspace x_{t-3} - R\thinspace x_{t-2})\thinspace (x_{t-1} - R\thinspace x_{t-2}) .
\end{aligned}
$$

Las fracciones de agentes se actualizan mediante el modelo **logit multinomial**
([sesión 13](s13_logit_ising.md)) con **intensidad de elección** $\beta \ge 0$:

$$
n_{h,t} = \frac{e^{\beta \pi_{h,t-1}}}{\sum_{k=1}^H e^{\beta \pi_{k,t-1}}}
\quad\Longrightarrow\quad
n_{2,t} = \frac{1}{1 + \exp\big(\beta\thinspace (\pi_{1,t-1} - \pi_{2,t-1})\big)} .
$$

- Si $\beta = 0$, los agentes no reaccionan a las ganancias y se dividen a
  mitades iguales ($n_{1,t} = n_{2,t} = 1/2$).
- Si $\beta \to \infty$, todos los agentes saltan instantáneamente a la regla
  que tuvo mayor ganancia en el paso previo.

### 2.4 El ciclo endógeno de calma, burbuja y corrección

Cuando $R < g < R^2$ (por ejemplo, con $R = 1.10$, $g = 1.20 < R^2 = 1.21$ y
costo fundamentalista $C = 1.0$, como en el cuaderno), la
interacción entre $x_t$ y $n_{2,t}$ genera una narrativa económica transparente:

1. **Calma engañosa cerca del fundamental:** Cuando $x_t \approx 0$, ambas
   reglas ganan casi lo mismo antes de costos ($\pi_1 \approx -C$,
   $\pi_2 \approx 0$). Como investigar fundamentales cuesta $C > 0$, los
   agentes abandonan la regla fundamentalista y se pasan a la regla chartista
   gratuita ($n_{2,t}$ sube hacia $1/(1 + e^{-\beta C})$).
2. **Despegue de la burbuja:** En cuanto la fracción de chartistas supera el
   umbral $n_{2,t} > R/g$, el factor de amplificación $n_{2,t}\thinspace g / R$ se vuelve
   mayor que $1$. Cualquier pequeña perturbación se amplifica: la tendencia se
   autorrealiza, los chartistas ganan dinero gracias a sus propias compras y
   atraen aún más seguidores ($n_{2,t} \to 1$).
3. **Freno y reversión:** Cuando $n_{2,t} \approx 1$, el precio crece al ritmo
   $x_{t-1} \approx (g/R)\thinspace x_{t-2}$. Como $g/R < R$ (porque $g < R^2$), la
   burbuja crece más despacio que el rendimiento libre de riesgo $R$: el exceso
   de retorno realizado $x_{t-1} - R\thinspace x_{t-2} \approx (g/R - R)\thinspace x_{t-2}$ cambia
   de signo y la apuesta fundamentalista obtiene una ganancia cuadrática
   $-R\thinspace x_{t-2}(x_{t-1} - R\thinspace x_{t-2}) - C = R(R - g/R)\thinspace x_{t-2}^2 - C > 0$. Los
   agentes regresan de golpe a la regla fundamentalista ($n_{2,t} \downarrow$) y
   el precio colapsa hacia $x = 0$, reiniciando el ciclo de manera caótica.

---

## 3. Fundamentalistas, chartistas y contagio estocástico

### 3.1 De Lux-Marchesi a un modelo estimable: Franke y Westerhoff (2012)

En un mercado con un número finito $N$ de agentes y órdenes estocásticas, el
cambio entre estrategias no ocurre sólo de manera determinista: incluye
**contagio social (herding)** como en el modelo de Kirman
([sesión 15](s15_votante_kirman.md)) y fluctuaciones de demanda. Lux y Marchesi
(1999) construyeron un modelo de agentes con optimistas, pesimistas y
fundamentalistas que se contagian entre sí, y mostraron que el ir y venir entre
estrategias genera colas pesadas y agrupamiento de volatilidad; Giardina y
Bouchaud (2003) llegaron a conclusiones parecidas con un mercado de agentes que
eligen entre reglas según su desempeño.

Franke y Westerhoff (2012) destilaron esos ingredientes en un modelo de pocas
ecuaciones que puede estimarse con datos. El cuaderno simula su versión
**DCA-HPM** (enfoque de elección discreta con *herding*, predisposición y
desalineamiento), con los parámetros que estimaron para el S&P 500. El
log-precio $p_t$ responde al exceso de demanda de una fracción $n^f_t$ de
fundamentalistas y una fracción $n^c_t = 1 - n^f_t$ de chartistas:

$$
p_t = p_{t-1} + \mu\bigl(n^f_{t-1} d^f_{t-1} + n^c_{t-1} d^c_{t-1}\bigr) ,
$$

con demandas medias

$$
d^f_t = \phi\thinspace(p^\ast - p_t) + \sigma_f\thinspace\varepsilon^f_t ,
\qquad
d^c_t = \chi\thinspace(p_t - p_{t-1}) + \sigma_c\thinspace\varepsilon^c_t ,
\qquad
\varepsilon^f_t, \varepsilon^c_t \sim \mathcal{N}(0, 1) .
$$

- **Fundamentalistas:** apuestan a la reversión hacia el log-precio
  fundamental $p^\ast$ con velocidad $\phi > 0$.
- **Chartistas:** siguen la tendencia reciente con intensidad $\chi > 0$, con
  una demanda más ruidosa ($\sigma_c > \sigma_f$).
- **Elección de estrategia:** la fracción de fundamentalistas es logit
  ([sesión 13](s13_logit_ising.md)) en el **atractivo relativo** $a_t$ del
  fundamentalismo frente al chartismo:

$$
n^f_t = \frac{1}{1 + e^{-\beta a_{t-1}}} ,
\qquad
a_t = \alpha_0 + \alpha_n\thinspace(n^f_t - n^c_t) + \alpha_p\thinspace(p_t - p^\ast)^2 .
$$

Los tres términos de $a_t$ son las tres fuerzas del modelo:

1. **Predisposición** ($\alpha_0$): preferencia de fondo por una de las dos
   estrategias.
2. **Contagio de manada** ($\alpha_n > 0$): la estrategia que ya es mayoritaria
   resulta más atractiva, lo que mantiene a la población en una misma fase
   durante un tiempo.
3. **Ancla fundamental (desalineamiento)** ($\alpha_p > 0$): cuando el precio se
   aleja demasiado del fundamental, crece el temor a que la burbuja reviente y
   los agentes regresan al fundamentalismo.

Los valores estimados son $\mu = 0.01$, $\beta = 1$, $\phi = 0.12$,
$\chi = 1.5$, $\sigma_f = 0.758$, $\sigma_c = 2.087$, $\alpha_0 = -0.327$,
$\alpha_n = 1.79$ y $\alpha_p = 18.43$.

### 3.2 Por qué este mecanismo genera los hechos estilizados

Escribimos $z_t = n^c_t$ para la fracción de chartistas. La varianza condicional
del retorno $r_t = p_t - p_{t-1}$ dado el estado del mercado en $t-1$ es

$$
\mathrm{Var}(r_t \mid \mathcal{F}_{t-1}) = \mu^2\bigl[(1 - z_{t-1})^2\thinspace\sigma_f^2 + z_{t-1}^2\thinspace\sigma_c^2\bigr] .
$$

Como $\sigma_c > \sigma_f$ y la fracción de chartistas $z_t$ cambia lentamente
en el tiempo:

1. **Agrupamiento de volatilidad ($\hat\rho_{\lvert r\rvert}(\tau) > 0$):** Cuando
   el mercado entra en una fase especulativa ($z_t$ alto), la varianza
   condicional permanece alta durante muchos pasos consecutivos; cuando el
   desalineamiento $(p_t - p^\ast)^2$ devuelve a los agentes al
   fundamentalismo ($z_t$ bajo), el mercado entra en una etapa de baja
   volatilidad. La persistencia de $\lvert r_t\rvert$ es el reflejo directo de
   la persistencia de la composición poblacional $z_t$, que el contagio
   $\alpha_n$ refuerza.
2. **Colas pesadas ($\kappa > 0$):** Una mezcla de normales con varianza
   cambiante (baja la mayor parte del tiempo y alta durante las fases
   chartistas) produce una distribución incondicional leptocúrtica
   ([sesión 5](s05_gbm_garch.md)). En el cuaderno la curtosis
   simulada es positiva pero bastante menor que la del S&P 500.
3. **Ausencia de autocorrelación lineal ($\hat\rho_r(\tau) \approx 0$):** El
   término de tendencia de los chartistas entra multiplicado por
   $\mu\chi = 0.015$, y la reversión fundamentalista es lenta
   ($\mu\phi = 0.0012$ por día), de modo que la parte predecible del retorno es
   pequeña frente al ruido y los retornos tienen autocorrelación lineal cercana
   a cero.

---

## 4. Cierre del ciclo: qué hechos estilizados explica cada modelo

Al llegar al final del curso podemos evaluar los modelos estudiados desde la
sesión 1 hasta la 20 frente a los tres hechos estilizados fundamentales de los
retornos financieros:

| Modelo | Sesión | $\hat\rho_r(\tau) \approx 0$ (sin memoria lineal) | Colas pesadas ($\kappa \gg 0$) | Agrupamiento de volatilidad ($\hat\rho_{\lvert r\rvert}(\tau) > 0$) | ¿El mecanismo es microeconómico? |
|---|---|---|---|---|---|
| **Caminata aleatoria gaussiana (GBM)** | 1, 3 | Sí | No ($\kappa = 0$) | No ($\hat\rho_{\lvert r\rvert}(\tau) = 0$) | No (choques externos i.i.d.) |
| **Ruido i.i.d. con colas pesadas ($t$-Student / Lévy)** | 2 | Sí | Sí | No (días independientes) | No (distribución impuesta) |
| **GARCH(1,1)** | 3 | Sí | Sí | Sí (decaimiento exponencial) | No (ecuación estadística para $\sigma_t^2$) |
| **Juego de la minoría básico ($N$ fijo)** | 13 | Depende de $\alpha$: sin predictibilidad ($H = 0$) sólo si $\alpha < \alpha_c$ | No (fluctuaciones casi gaussianas) | No (aparece en variantes gran canónicas con número variable de agentes activos; Challet, Marsili & Zhang, 2005) | Sí (competencia inductiva por estar en la minoría) |
| **Fundamentalistas y chartistas (Lux-Marchesi, Franke-Westerhoff)** | 14 | Sí | Sí, aunque menores que en los datos | Sí (persistencia larga por cambio de régimen en $z_t$) | **Sí** (selección evolutiva + contagio entre reglas) |

---

## 5. [Avanzado] Bifurcaciones y caos en el modelo de Brock-Hommes

El desarrollo de esta sección está en [Para profundizar, P.1](#pp1).


---

## Tabla de conceptos clave

| Concepto | Fórmula / condición | Interpretación económica |
|---|---|---|
| **Desviación respecto al fundamental** | $x_t = p_t - p^\ast$, con $p^\ast = \bar y / (R - 1)$ | Burbuja positiva ($x_t > 0$) o subvaluación ($x_t < 0$) |
| **Vaciado de mercado (Brock-Hommes)** | $R\thinspace x_t = \sum_{h=1}^H n_{h,t}\thinspace f_{h,t}$ | El precio de hoy promedia los pronósticos ponderados por su cuota de mercado |
| **Selección evolutiva logit** | $n_{h,t} = e^{\beta \pi_{h,t-1}} / \sum_k e^{\beta \pi_{k,t-1}}$ | $\beta$ mide qué tan rápido migran los agentes hacia la regla más rentable |
| **Umbral de inestabilidad ($\beta_c$)** | $g\thinspace n_2^\ast(\beta_c) = R \implies \beta_c = \frac{1}{C}\ln\left(\frac{R}{g-R}\right)$ | Si ser fundamentalista cuesta $C > 0$, la calma hace migrar a los agentes al chartismo y desestabiliza $x^\ast = 0$ |
| **Agrupamiento de volatilidad endógeno** | $\mathrm{Var}(r_t \mid \mathcal{F}_{t-1}) = \mu^2\bigl[(1-z_{t-1})^2\sigma_f^2 + z_{t-1}^2\sigma_c^2\bigr]$ | Las olas de chartismo ($z_t \uparrow$) crean ráfagas de alta volatilidad y colas pesadas |

---

## Ejercicios

### Ruta básica

**B1. Fracción logit de chartistas (`fraccion_logit_dos`).** Si en un mercado
conviven fundamentalistas con ganancia neta $\pi_1$ y chartistas con ganancia
neta $\pi_2$, escribe la expresión de la fracción de chartistas $n_2$ bajo
elección logit con intensidad $\beta$. ¿Cuánto vale $n_2$ cuando $\beta = 0$?
¿Y cuando $\pi_1 = -C < 0$ y $\pi_2 = 0$?

**B2. Paso de precio en Brock-Hommes (`paso_precio_bh`).** Con factor libre de
riesgo $R = 1.01$, fundamentalistas ($f_{1,t} = 0$) y seguidores de tendencia
($f_{2,t} = g\thinspace x_{t-1}$ con $g = 1.2$), si ayer la desviación fue $x_{t-1} = 2.0$
y hoy el 85% de los agentes son chartistas ($n_{2,t} = 0.85$), ¿cuál es la
nueva desviación $x_t$? ¿La burbuja creció o se contrajo?

**B3. Diagnóstico de hechos estilizados.** Explica por qué barajar al azar el
orden temporal de los retornos simulados de un mercado artificial conserva la
curtosis en exceso $\kappa$, pero destruye por completo la autocorrelación de la
volatilidad $\hat\rho_{|r|}(\tau)$.

### Ruta avanzada

**A1. Intensidad crítica de bifurcación (`beta_critico_bh`).** Demuestra que
cuando $R < g < 2R$ y $C > 0$, la condición $g\thinspace n_2^\ast(\beta_c) = R$ implica
$\beta_c = \frac{1}{C}\ln\left(\frac{R}{g - R}\right)$. Evalúa $\beta_c$ para
$R = 1.01$, $g = 1.20$ y $C = 1.0$. ¿Qué ocurre con $\beta_c$ cuando el costo
de información $C$ tiende a cero?

**A2. Contrarios y duplicación de periodo.** Supón que el segundo grupo de
agentes no sigue la tendencia sino que es **contrario**: $f_{2,t} = g\thinspace x_{t-1}$
con $-2R < g < -R$ (por ejemplo, $R = 1.10$ y $g = -1.35$), y que los
fundamentalistas pagan $C = 1$. Calcula el eigenvalor $\lambda_1 = g\thinspace n_2^\ast / R$
en el origen, muestra que $\lambda_1 > -1$ cuando $\beta = 0$ y encuentra el
valor de $\beta$ en que cruza $\lambda_1 = -1$. Explica por qué en ese cruce
nace un ciclo de **periodo 2** donde la desviación $x_t$ alterna de signo en
cada paso. ¿Qué pasaría con $C = 0$?

---

## Para profundizar

Material opcional de la ruta avanzada: demostraciones y temas que no se ven
en clase. Cada apartado indica de qué sección de las notas viene.

<a id="pp1"></a>

### P.1 Bifurcaciones y caos en el modelo de Brock-Hommes (sección 5)


#### 5.1 Puntos fijos y estabilidad lineal del equilibrio fundamental

Estudiemos el sistema determinista de Brock y Hommes (1998) de la sección 2 con
fundamentalistas ($f_{1,t} = 0$, costo $C \ge 0$) y seguidores de tendencia
($f_{2,t} = g\thinspace x_{t-1}$, costo $0$). Sustituyendo la fracción logit $n_{2,t}$ en
$R\thinspace x_t = n_{2,t}\thinspace g\thinspace x_{t-1}$, obtenemos el sistema dinámico en tiempo discreto
([C.A5](repaso/calculo.md#ca5)):

$$
\begin{aligned}
x_t &= \frac{g}{R}\thinspace n_{2,t}\thinspace x_{t-1} , \cr
n_{2,t} &= \frac{1}{1 + \exp\Big(\beta\thinspace \big[-C - g\thinspace x_{t-3}\thinspace (x_{t-1} - R\thinspace x_{t-2})\big]\Big)} .
\end{aligned}
$$

**Proposición 5.1 (Estabilidad local del estado fundamental).** El estado
estacionario fundamental $x^\ast = 0$ existe para todos los valores de los
parámetros. En él, la fracción estacionaria de chartistas es

$$
n_2^\ast = \frac{1}{1 + e^{-\beta C}} \in \left[\frac{1}{2},\thinspace  1\right) .
$$

Además:
1. Si $0 \le g < R$, el punto fijo fundamental $x^\ast = 0$ es globalmente estable
   para todo $\beta \ge 0$.
2. Si $R < g < 2R$ y $C > 0$, existe una intensidad de elección crítica

   $$
   \boxed{\beta_c = \frac{1}{C}\thinspace \ln\left(\frac{R}{g - R}\right)}
   $$

   tal que $x^\ast = 0$ es localmente asintóticamente estable para $\beta < \beta_c$
   y se vuelve **inestable** para $\beta > \beta_c$ (donde $n_2^\ast > R/g$).

*Demostración.* Cerca de $x_t = x_{t-1} = x_{t-2} = x_{t-3} = 0$, el término
cuadrático $g\thinspace x_{t-3}\thinspace (x_{t-1} - R\thinspace x_{t-2})$ es de segundo orden
$\mathcal{O}(x^2)$. Por tanto, la linealización de $n_{2,t}$ alrededor del
origen evalúa $n_{2,t}$ en $n_2^\ast = 1/(1 + e^{-\beta C}) + \mathcal{O}(x^2)$, y
la ecuación linealizada para $x_t$ es simplemente:

$$
x_t = \left(\frac{g\thinspace n_2^\ast}{R}\right) x_{t-1} + \mathcal{O}(x^3) .
$$

Los eigenvalores del jacobiano en el origen son $\lambda_1 = g\thinspace n_2^\ast / R$ y
$\lambda_2 = \lambda_3 = 0$. El punto fijo $x^\ast = 0$ es localmente estable si y
sólo si $|\lambda_1| < 1$, es decir, $g\thinspace n_2^\ast < R$.
- Como $n_2^\ast < 1$, si $0 \le g < R$ siempre se cumple $g\thinspace n_2^\ast < R$.
- Si $R < g < 2R$ y $C > 0$, para $\beta = 0$ tenemos $n_2^\ast = 1/2 < R/g$, por
  lo que el origen es estable. Al crecer $\beta$, $n_2^\ast(\beta)$ crece
  monótonamente desde $1/2$ hacia $1$ y cruza el valor crítico $R/g$ exactamente
  cuando $1/(1 + e^{-\beta_c C}) = R/g$, de donde $e^{-\beta_c C} = (g - R)/R$,
  es decir, $\beta_c = \frac{1}{C}\ln\big(\frac{R}{g - R}\big)$. $\square$

#### 5.2 La bifurcación pitchfork y la ruta al caos

¿Qué ocurre cuando $\beta$ supera el valor crítico $\beta_c$ en el régimen
canónico $R < g < R^2$ con $C > 0$ (por ejemplo, $R = 1.10$, $g = 1.20$,
$C = 1.0$, donde $\beta_c = \ln(11) \approx 2.40$)?

1. **Primera bifurcación ($\beta = \beta_c$, pitchfork):** El origen $x^\ast = 0$
   se vuelve inestable y nacen dos puntos fijos simétricos no fundamentales
   $\pm x^\ast \ne 0$ en los cuales la fracción de chartistas se equilibra en
   $n_2^\ast = R/g$:
   $$
   n_2^\ast = \frac{R}{g} \quad\Longleftrightarrow\quad
   (x^\ast)^2 = \frac{C - \frac{1}{\beta}\ln\left(\frac{R}{g - R}\right)}{g\thinspace (R - 1)}
   = \frac{C\thinspace \big(1 - \beta_c / \beta\big)}{g\thinspace (R - 1)} > 0 \quad (\text{para } \beta > \beta_c) .
   $$
2. **Después: ciclos y caos.** Al seguir aumentando $\beta$, los equilibrios
   $\pm x^\ast$ también pierden estabilidad, aparecen ciclos cuasiperiódicos
   y, para $\beta$ mayor, **caos determinista**: las burbujas crecen hasta que
   volver al fundamentalismo se vuelve rentable y la órbita regresa cerca de
   $x = 0$, una y otra vez, sin repetirse. Es lo que muestra el diagrama de
   bifurcación del cuaderno. El cuaderno también estima el **exponente de
   Lyapunov** $\lambda$: si $\lambda > 0$, dos trayectorias que empiezan casi
   iguales se separan exponencialmente. La clasificación de estas
   bifurcaciones y la ruta al caos están en Brock y Hommes (1998) y en Hommes
   (2013); Brock y Hommes (1997) la explican para $\beta$ grande con órbitas
   homoclínicas.

---

## Referencias

- Brock, W. A., & Hommes, C. H. (1997). A rational route to randomness.
  *Econometrica*, 65(5), 1059–1095. https://doi.org/10.2307/2171879
- Brock, W. A., & Hommes, C. H. (1998). Heterogeneous beliefs and routes to
  chaos in a simple asset pricing model. *Journal of Economic Dynamics and
  Control*, 22(8–9), 1235–1274. https://doi.org/10.1016/S0165-1889(98)00011-6
- Challet, D., Marsili, M., & Zhang, Y.-C. (2005). *Minority Games: Interacting
  Agents in Financial Markets*. Oxford University Press.
- Cont, R. (2001). Empirical properties of asset returns: stylized facts and
  statistical issues. *Quantitative Finance*, 1(2), 223–236.
  https://doi.org/10.1080/713665670
- Cutler, D. M., Poterba, J. M., & Summers, L. H. (1989). What moves stock
  prices? *The Journal of Portfolio Management*, 15(3), 4–12.
  https://doi.org/10.3905/jpm.1989.409212
- Franke, R., & Westerhoff, F. (2012). Structural stochastic volatility in asset
  pricing dynamics: Estimation and model contest. *Journal of Economic Dynamics
  and Control*, 36(8), 1193–1211. https://doi.org/10.1016/j.jedc.2011.10.004
- Giardina, I., & Bouchaud, J.-P. (2003). Bubbles, crashes and intermittency in
  agent based market models. *The European Physical Journal B*, 31(3), 421–437.
  https://doi.org/10.1140/epjb/e2003-00050-6
- Hommes, C. H. (2006). Heterogeneous agent models in economics and finance. En
  L. Tesfatsion & K. L. Judd (Eds.), *Handbook of Computational Economics*
  (Vol. 2, pp. 1109–1186). Elsevier.
  https://doi.org/10.1016/S1574-0021(05)02023-X
- Hommes, C. H. (2013). *Behavioral Rationality and Heterogeneous Expectations
  in Complex Economic Systems*. Cambridge University Press.
  https://doi.org/10.1017/CBO9781139094276
- LeBaron, B. (2006). Agent-based financial markets: Matching stylized facts
  with style. En D. Colander (Ed.), *Post Walrasian Macroeconomics*
  (pp. 221–236). Cambridge University Press.
  https://doi.org/10.1017/CBO9780511617751.013
- Lux, T., & Marchesi, M. (1999). Scaling and criticality in a stochastic
  multi-agent model of a financial market. *Nature*, 397(6719), 498–500.
  https://doi.org/10.1038/17290
- Slanina, F. (2014). *Essentials of Econophysics Modelling*. Oxford
  University Press.
