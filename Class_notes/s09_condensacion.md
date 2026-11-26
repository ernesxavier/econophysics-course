# Sesión 9: Condensación de riqueza y redistribución. Notas formales

> **Cómo usar estas notas.** La clase combina el cuaderno
> (`Code/s09_condensacion.jl`) y el pizarrón. El **plan de clase** de abajo
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
> Las secciones marcadas con ★ son el **núcleo** que se ve en clase; las
> demás se leen por cuenta propia y se evalúan igual.
>
> Repaso recomendado: [E.B3](repaso/procesos_estocasticos.md#eb3)
> (promediar simulaciones), [F.B4](repaso/fisica_estadistica.md#fb4)
> (transiciones de fase), [E.I6](repaso/procesos_estocasticos.md#ei6)
> (Gibrat y Kesten), [F.I6](repaso/fisica_estadistica.md#fi6) (modelos de
> intercambio); para la ruta avanzada,
> [E.A4](repaso/procesos_estocasticos.md#ea4) (Fokker-Planck) y
> [C.A4](repaso/calculo.md#ca4) (función gamma).
>
> Lecturas: Chakrabarti et al. (2013), sec. 4.5. Ruta avanzada: Bouchaud &
> Mézard (2000); Ribeiro (2020), de consulta.

**Plan de clase** (120 minutos: 105 de núcleo y 15 de margen). El núcleo:
el modelo de Bouchaud y Mézard, su distribución estacionaria, la
condensación y el efecto de un impuesto redistribuido.

| Minutos | Dónde | Qué | Secciones |
|---|---|---|---|
| 0–8 | Proyector | El problema: inversión sin intercambio | Cuaderno, sección 0 |
| 8–19 | Pizarra | El modelo de Bouchaud y Mézard: inversión más intercambio; la riqueza relativa; sin intercambio, Gibrat: un agente acaba con casi todo | 1.1–1.2 |
| 19–28 | Pizarra | La distribución estacionaria: la gamma inversa con $\mu = 1 + 2J/\sigma^2$ (enunciado) y su lectura; la tabla del Gini | 2.1, 2.3 |
| 28–37 | Pizarra | La razón de participación $Y_2$ como parámetro de orden; condensación y exponente de la cola: por qué en campo medio no hay condensación | 3.1–3.2 |
| 37–45 | Pizarra | En una red sí hay transición; el impuesto redistribuido equivale a más intercambio, con la demostración | 3.3, 4.1 |
| 45–90 | Parejas | Laboratorio: el modelo, la transición de condensación y el efecto de un impuesto; ejercicio de código 1 | Cuaderno, secciones 1 a 4 y ejercicio 1 |
| 90–105 | Pizarra | Cierre: puesta en común de las preguntas 3 a 5 y la tabla del impuesto en la red | Cuaderno, preguntas 3–5; notas, 4.1 |
| 105–120 | Margen | Si alcanza: el ejercicio de código 2 (el exponente de Pareto) y su comparación con la cola de la simulación | |
| Después | Lectura | Lo que no entró del margen; las convenciones de Itô y Stratonovich (1.1), los resultados en la red (3.3), el resumen (5) y Para profundizar | |

Lo del **margen** se ve sólo si sobra tiempo, en el orden indicado; lo
que no se vea pasa a la lectura y se evalúa igual.

**Pregunta de la sesión.** Con las mismas oportunidades para todos, ¿por qué
unos pocos acumulan casi toda la riqueza, y qué lo evita?

---

## 1. El modelo de Bouchaud y Mézard ★

### 1.1 Definición ★

**Modelo.** $N$ agentes con riquezas $W_1, \dots, W_N > 0$ evolucionan según

$$
dW_i = \sigma W_i\thinspace dB_i + J\left(\bar W - W_i\right)dt,
\qquad \bar W = \frac{1}{N}\sum_{j} W_j ,
$$

con movimientos brownianos $B_i$ independientes, en la convención de Itô.

- El primer término es la **inversión**: un rendimiento aleatorio
  proporcional a la riqueza, con volatilidad $\sigma$. Sin él, el modelo es
  un intercambio puro.
- El segundo es el **intercambio** en campo medio: cada agente recibe
  $J\bar W\thinspace dt$ de los demás y gasta $JW_i\thinspace dt$.

El intercambio conserva la riqueza total, y el término de inversión tiene
media cero, así que $\mathbb{E}[\bar W]$ es constante. Para $N$ grande,
$\bar W$ es prácticamente constante (ley de los grandes números) y basta
estudiar la **riqueza relativa** $w_i = W_i / \bar W$, con promedio 1:

$$
dw = \sigma w\thinspace dB + J(1 - w)\thinspace dt .
$$

**Sobre las convenciones.** Bouchaud y Mézard (2000) escriben
$dW_i/dt = \eta_i W_i + \dots$ con un ruido $\eta_i$ de varianza
$2\sigma_{BM}^2$, en la convención de Stratonovich. Su $\sigma_{BM}^2$
equivale a $\sigma^2/2$ en estas notas. Todas las fórmulas de aquí están en
la convención de Itô; para traducir resultados del artículo, cambia
$\sigma_{BM}^2$ por $\sigma^2/2$.

**Escala.** Si se cambia la unidad de tiempo, $J$ y $\sigma^2$ se
multiplican por el mismo factor. Las propiedades de largo plazo sólo
dependen del cociente $J/\sigma^2$; en el cuaderno, $\sigma = 1$.

### 1.2 Sin intercambio: Gibrat ★

Con $J = 0$, cada $W_i$ es un movimiento browniano geométrico:
$\ln W_i(t) = \ln W_i(0) - \sigma^2 t/2 + \sigma B_i(t)$
([E.I6](repaso/procesos_estocasticos.md#ei6)). Los logaritmos se dispersan
con varianza $\sigma^2 t$, sin límite.

**Consecuencia (heurística).** Con $N$ fijo, la diferencia entre el mayor y
el segundo mayor de los $\ln W_i$ crece como $\sigma\sqrt{t}$, así que el
cociente entre la riqueza del segundo y la del primero tiende a cero: **un
solo agente acaba con casi toda la riqueza**, aunque todos tengan las mismas
oportunidades. En el cuaderno, con $N = 1000$ y $\sigma = 1$, el más rico
tiene casi el 100% hacia $t = 100$.

---

## 2. La distribución estacionaria

### 2.1 Solución de Fokker-Planck ★

**Proposición.** Para $J > 0$, la densidad estacionaria de la riqueza
relativa es la **gamma inversa**

$$
p(w) = \frac{(\mu - 1)^{\mu}}{\Gamma(\mu)}\thinspace\frac{e^{-(\mu-1)/w}}{w^{1+\mu}},
\qquad w > 0,
\qquad \boxed{\ \mu = 1 + \frac{2J}{\sigma^2}\ }
$$

Tiene media 1 y cola de Pareto, $\mathbb{P}(w > x) \sim C x^{-\mu}$, con
exponente $\mu$.

**[Avanzado]** Demostración: en [Para profundizar, P.1](#pp1).

**Interpretación.**

- Más intercambio ($J$ mayor) o inversiones menos riesgosas ($\sigma$
  menor) dan $\mu$ mayor: menos desigualdad.
- Cerca de $w = 0$ el factor $e^{-(\mu-1)/w}$ anula la densidad: nadie
  queda en la ruina, porque el intercambio siempre trae algo.

### 2.2 [Avanzado] Comprobación con Kesten

El desarrollo de esta subsección está en [Para profundizar, P.2](#pp2).


### 2.3 El índice de Gini ★

**Proposición.** El Gini de la gamma inversa con exponente $\mu > 1$ es igual
al de una gamma con forma $\mu - 1$ (sesión 8):

$$
G(\mu) = \frac{\Gamma(\mu - 1/2)}{(\mu - 1)\thinspace\Gamma(\mu - 1)\sqrt{\pi}} .
$$

**[Avanzado]** Demostración: en [Para profundizar, P.3](#pp3).

| $\mu$ | 1.2 | 1.5 | 2 | 3 | 5 |
|---|---|---|---|---|---|
| $G$ | 0.80 | 0.64 | 0.50 | 0.375 | 0.27 |

---

## 3. Condensación ★

### 3.1 El parámetro de orden ★

**Definición.** La **razón de participación inversa** es

$$
Y_2 = \sum_{i=1}^{N} f_i^2, \qquad f_i = \frac{W_i}{\sum_j W_j} .
$$

Cumple $1/N \le Y_2 \le 1$: vale $1/N$ si todos tienen lo mismo y 1 si uno
tiene todo. $1/Y_2$ es el número efectivo de agentes que tienen la riqueza.

**Definición.** Hay **condensación** si $\mathbb{E}[Y_2]$ tiende a un valor
positivo cuando $N \to \infty$: unos pocos agentes tienen una fracción
finita de la riqueza total, por grande que sea la población.
$\mathbb{E}[Y_2]$ es el **parámetro de orden**
([F.B4](repaso/fisica_estadistica.md#fb4)).

### 3.2 Condensación y exponente de la cola ★

Supongamos que las riquezas son independientes con cola de Pareto de
exponente $\mu$. El mayor de $N$ valores es del orden de $N^{1/\mu}$
(sesión 3).

- **Si $\mu < 1$**, la media es infinita y la suma de los $N$ valores también
  es del orden de $N^{1/\mu}$: el mayor tiene una fracción finita del total.
  Hay condensación; de hecho $\mathbb{E}[Y_2] \to 1 - \mu$ (resultado citado por
  Bouchaud y Mézard, 2000).
- **Si $1 < \mu < 2$**, la suma crece como $N$ y la suma de cuadrados está
  dominada por el mayor, del orden de $N^{2/\mu}$. Entonces
  $Y_2 \sim N^{2/\mu - 2}$, que tiende a cero, pero **muy despacio** cuando
  $\mu$ está cerca de 1.
- **Si $\mu > 2$**, la varianza es finita y $Y_2 \approx \mathbb{E}[w^2]/N \sim 1/N$.

**En campo medio**, $\mu = 1 + 2J/\sigma^2 > 1$ para todo $J > 0$: **no hay
condensación**. Pero con $J$ pequeño, $Y_2$ baja muy lentamente con $N$,
como se ve en el cuaderno (con $J = 0.05$, de 0.10 a 0.06 al pasar de 1,000
a 20,000 agentes).

### 3.3 Condensación en una red ★

Si cada agente sólo intercambia con sus $c_i$ socios en una red,

$$
dW_i = \sigma W_i\thinspace dB_i + \frac{J}{c_i}\sum_{j \sim i}\left(W_j - W_i\right)dt ,
$$

la riqueza de un agente afortunado sólo se derrama a sus pocos socios. El
argumento de campo medio ya no vale, y aparece una **transición de fase**
en un valor finito $J_c$: para $J < J_c$ la riqueza se condensa.

**Resultados.** Para una red aleatoria con 4 socios en promedio:

- Bouchaud y Mézard (2000), midiendo el exponente de la cola, encontraron
  la transición en $J/\sigma_{BM}^2 \approx 0.3$, es decir,
  $J_c \approx 0.15\thinspace\sigma^2$ en nuestra convención.
- En las simulaciones del curso (`Code/simular_condensacion.jl`, con $N$ de
  1,000 a 20,000), $Y_2$ no cambia con $N$ para $J \le 0.1$ y baja con $N$
  para $J \ge 0.4$ (con $\sigma = 1$). Entre ambos valores el comportamiento
  no es monótono: con estos tamaños no se puede ubicar $J_c$ con más
  precisión.

Es una lección general sobre transiciones de fase: en sistemas finitos se
ven suavizadas, y ubicarlas con precisión exige tamaños cada vez mayores.

---

## 4. Impuestos y redistribución

### 4.1 Impuesto a la riqueza redistribuido ★

Un impuesto con tasa $\varphi$ sobre la riqueza, repartido por igual entre
todos, agrega a $dW_i$ el término $\varphi(\bar W - W_i)\thinspace dt$.

**Proposición.** En campo medio, el impuesto redistribuido equivale a más
intercambio:

$$
\boxed{\ \mu = 1 + \frac{2(J + \varphi)}{\sigma^2}\ }
$$

*Demostración.* El término del impuesto tiene la misma forma que el del
intercambio, así que basta sustituir $J$ por $J + \varphi$ en la sección 2.
$\square$

**En la red** la diferencia es importante: el impuesto llega a **todos**,
mientras que el intercambio sólo llega a los socios. Simulaciones de la red
con $J = 0.1$ (fase condensada), de 1,000 a 20,000 agentes
(`Code/simular_condensacion.jl`, datos en `Data/condensacion_impuesto.csv`;
mediana de 4 simulaciones):

| $\varphi$ | $Y_2$, $N = 1{,}000$ | $Y_2$, $N = 5{,}000$ | $Y_2$, $N = 20{,}000$ |
|---|---|---|---|
| 0 | 0.104 | 0.105 | 0.098 |
| 0.025 | 0.064 | 0.039 | 0.039 |
| 0.05 | 0.047 | 0.025 | 0.022 |
| 0.1 | 0.029 | 0.013 | 0.0094 |
| 0.2 | 0.012 | 0.0051 | 0.0026 |

Sin impuesto, $Y_2$ no cambia con $N$. Con $\varphi = 0.05$, $Y_2$ es la
cuarta parte que sin impuesto con 20,000 agentes, y con $\varphi \ge 0.1$
baja claramente al crecer $N$: la economía sale de la fase condensada.

### 4.2 [Avanzado] Cobrar no es lo mismo que redistribuir

El desarrollo de esta subsección está en [Para profundizar, P.4](#pp4).


---

## 5. Resumen formal

| Situación | Resultado |
|---|---|
| $J = 0$ (Gibrat) | $\mathrm{Var}(\ln W) = \sigma^2 t$; un agente acaba con casi todo |
| Campo medio, $J > 0$ | Gamma inversa con $\mu = 1 + 2J/\sigma^2$; Gini de una gamma con forma $\mu - 1$ |
| Parámetro de orden | $Y_2 = \sum_i f_i^2$; condensación si $\mathbb{E}[Y_2] \not\to 0$ |
| Pareto con $\mu < 1$ | Condensación, $\mathbb{E}[Y_2] \to 1 - \mu$ |
| Red con 4 socios | Condensación para $J < J_c$, con $J_c \approx 0.15\thinspace\sigma^2$ según Bouchaud y Mézard |
| Impuesto redistribuido | En campo medio, $J \to J + \varphi$; en la red, elimina la condensación |

---

## 6. Ejercicios formales

**Ruta básica**

- **F1.** Con $\sigma^2 = 0.04$ por año (rendimientos con 20% de
  volatilidad anual), ¿qué $J$ da $\mu = 1.5$? Interpreta: ¿qué fracción de
  su riqueza tendría que gastar cada agente al año?
- **F2.** Calcula $Y_2$ para una población donde 10 agentes tienen cada uno
  el 5% de la riqueza y los demás se reparten el resto por igual. ¿Qué pasa
  cuando $N \to \infty$?
- **F3.** Con la tabla de la sección 2.3, ¿qué tasa de impuesto
  redistribuido hace falta para bajar el Gini de 0.8 a 0.5 si $\sigma = 1$
  y $J = 0.1$?

**Ruta avanzada**

- **F4.** Verifica que la densidad de la sección 2.1 integra 1 y tiene
  media 1, y calcula su varianza. ¿Para qué $\mu$ es infinita?
- **F5.** Para $1 < \mu < 2$, justifica con más cuidado que
  $Y_2 \sim N^{2/\mu - 2}$. ¿Cuántos agentes harían falta para que $Y_2$ baje
  a la mitad de su valor con $N = 1000$, si $\mu = 1.1$?
- **F6.** Escribe la ecuación de Fokker-Planck del modelo con impuesto
  redistribuido y con un ingreso fijo $I$ igual para todos
  ($dW_i = \dots + I\thinspace dt$). ¿Cambia el exponente de la cola?
- **F7.** Demuestra la afirmación de la sección 4.2: un impuesto
  proporcional no redistribuido no cambia la distribución de la riqueza
  relativa.

---

## Para profundizar

Material opcional de la ruta avanzada: demostraciones y temas que no se ven
en clase. Cada apartado indica de qué sección de las notas viene.

<a id="pp1"></a>

### P.1 Demostración (sección 2.1)

*Demostración.* La difusión tiene deriva $A(w) = J(1-w)$ y
coeficiente de ruido $B(w) = \sigma w$. Su solución estacionaria con flujo
nulo es ([E.A4](repaso/procesos_estocasticos.md#ea4))

$$
p(w) \propto \frac{1}{B(w)^2}\exp\left(\int^w \frac{2A(y)}{B(y)^2}\thinspace dy\right)
= \frac{1}{\sigma^2 w^2}\exp\left(\frac{2J}{\sigma^2}\int^w \left(\frac{1}{y^2} - \frac{1}{y}\right) dy\right).
$$

La integral vale $-1/w - \ln w$, así que
$p(w) \propto w^{-2}\thinspace w^{-2J/\sigma^2}\thinspace e^{-(2J/\sigma^2)/w} = w^{-1-\mu} e^{-(\mu-1)/w}$.
Con el cambio $u = (\mu-1)/w$, la normalización y la media se reducen a
integrales de la función gamma ([C.A4](repaso/calculo.md#ca4)), que dan la
constante y media 1. Para $w$ grande, $e^{-(\mu-1)/w} \to 1$ y queda
$p(w) \sim w^{-1-\mu}$; integrando, la cola de la CCDF tiene exponente
$\mu$. $\square$

<a id="pp2"></a>

### P.2 Comprobación con Kesten (sección 2.2)


La simulación del cuaderno avanza en pasos $dt$ con la regla

$$
w_{t+dt} = a_t\thinspace w_t + J\thinspace dt,
\qquad a_t = (1 - J\thinspace dt)\thinspace e^{\sigma\sqrt{dt}\thinspace\xi_t - \sigma^2 dt/2},
$$

con $\xi_t$ normal estándar (y una normalización final que no cambia la
forma de la distribución). Es un proceso de Kesten (sesión 6), cuyo
exponente resuelve $\mathbb{E}[a^\alpha] = 1$. Como
$\mathbb{E}\left[e^{\alpha\sigma\sqrt{dt}\xi}\right] = e^{\alpha^2\sigma^2 dt/2}$,

$$
\mathbb{E}[a^\alpha] = (1 - J\thinspace dt)^\alpha\thinspace e^{\alpha(\alpha - 1)\sigma^2 dt/2} = 1
\quad\Longrightarrow\quad
-\alpha J + \frac{\alpha(\alpha-1)\sigma^2}{2} = 0
$$

a primer orden en $dt$, de donde $\alpha = 1 + 2J/\sigma^2 = \mu$. Los dos
caminos dan el mismo exponente.

<a id="pp3"></a>

### P.3 Demostración (sección 2.3)

*Demostración.* El Gini no depende de la escala, así que
basta tomar $X = 1/Y$ con $Y$ gamma de forma $\mu$ y escala 1. Sean $Y, Y'$
independientes. Entonces

$$
\mathbb{E}\lvert X - X'\rvert = \mathbb{E}\left[\frac{\lvert Y - Y'\rvert}{Y Y'}\right].
$$

Multiplicar la densidad conjunta $y^{\mu-1}y'^{\mu-1}e^{-y-y'}/\Gamma(\mu)^2$
por $1/(yy')$ da la densidad de dos gammas independientes de forma $\mu - 1$,
multiplicada por $(\Gamma(\mu-1)/\Gamma(\mu))^2$. Así,
$\mathbb{E}\lvert X - X'\rvert = \left(\Gamma(\mu-1)/\Gamma(\mu)\right)^2\mathbb{E}\lvert Z - Z'\rvert$,
con $Z, Z'$ gammas de forma $\mu - 1$. Como $\mathbb{E}[X] = \Gamma(\mu-1)/\Gamma(\mu)$ y
$\mathbb{E}[Z] = \mu - 1 = \Gamma(\mu)/\Gamma(\mu-1)$,

$$
G_X = \frac{\mathbb{E}\lvert X - X'\rvert}{2\thinspace\mathbb{E}[X]}
= \frac{\Gamma(\mu-1)}{\Gamma(\mu)}\cdot\frac{\mathbb{E}\lvert Z - Z'\rvert}{2}
= \frac{\mathbb{E}\lvert Z - Z'\rvert}{2\thinspace\mathbb{E}[Z]} = G_Z .
$$

El resultado sale de la fórmula de la sesión 8 para el Gini de la gamma. $\square$

<a id="pp4"></a>

### P.4 Cobrar no es lo mismo que redistribuir (sección 4.2)


- Un impuesto a la riqueza que **no** se redistribuye, por sí solo, quita a
  todos la misma fracción: la riqueza relativa no cambia, y la distribución
  tampoco.
- Bouchaud y Mézard (2000) analizan un caso más general, con impuesto al
  ingreso $\varphi_I$ y al capital $\varphi_C$, de los que se redistribuyen las
  fracciones $f_I$ y $f_C$. Encuentran que el impuesto al ingreso aumenta
  $\mu$ (reduce la desigualdad), pero un impuesto al capital combinado con
  uno al ingreso puede **reducir** $\mu$ si se redistribuye poco: sólo si
  $f_C > f_I\varphi_I/(1 + \varphi_I)$ el impuesto al capital reduce la
  desigualdad.

---

## Referencias

- Bouchaud, J.-P., & Mézard, M. (2000). Wealth condensation in a simple
  model of economy. *Physica A*, 282(3–4), 536–545.
  https://doi.org/10.1016/S0378-4371(00)00205-3
- Chakrabarti, B. K., Chakraborti, A., Chakravarty, S. R., & Chatterjee, A.
  (2013). *Econophysics of Income and Wealth Distributions*. Cambridge
  University Press. https://doi.org/10.1017/CBO9781139004169
- Ribeiro, M. B. (2020). *Income Distribution Dynamics of Economic
  Systems: An Econophysical Approach*. Cambridge University Press.
  https://doi.org/10.1017/9781316136119
