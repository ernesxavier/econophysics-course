# Sesión 17: Replicador y ciclo de Goodwin. Notas formales

> **Cómo leer estas notas.** Siguen el orden del programa: primero se observa
> en el cuaderno (`Code/s17_replicador_goodwin.jl`), después se formaliza aquí.
> Lo que no tiene marca es **ruta básica**. Las secciones y demostraciones
> marcadas **[Avanzado]** son opcionales.
>
> Repaso recomendado: [C.B3](repaso/calculo.md#cb3) (la derivada como tasa de
> cambio), [C.I1](repaso/calculo.md#ci1) (desarrollo de Taylor),
> [C.I3](repaso/calculo.md#ci3) (funciones de varias variables y derivadas
> parciales), [C.I4](repaso/calculo.md#ci4) (ecuaciones diferenciales en una
> dimensión, método de Euler y retrato de fase) y
> [L.I1](repaso/algebra_lineal.md#li1) (eigenvalores y eigenvectores); para la
> ruta avanzada, [C.A2](repaso/calculo.md#ca2) (sistemas de EDO en dos
> dimensiones) y [L.A3](repaso/algebra_lineal.md#la3) (jacobiano y estabilidad
> lineal).
>
> Lecturas: Requião da Cunha (2021), cap. 7. Ruta avanzada: Goodwin (1967); de
> consulta: Strogatz (2019), caps. 5–6; Schecter & Gintis (2016).

---

## 1. Juegos evolutivos y la ecuación del replicador

### 1.1 Juegos simétricos de $2 \times 2$ y pagos esperados

Consideremos una población grande de agentes que interactúan en pares y eligen
entre dos estrategias puras, $E_1$ y $E_2$. Los pagos de una interacción se
resumen en la matriz de pagos del jugador fila frente al jugador columna:

$$
\begin{array}{c|cc}
 & E_1 & E_2 \cr
\hline
E_1 & a & b \cr
E_2 & c & d
\end{array}
\qquad\Longleftrightarrow\qquad
A = \begin{pmatrix} a & b \cr c & d \end{pmatrix} .
$$

Si la población está **bien mezclada** (cada agente se empareja al azar con
cualquier otro) y una fracción $x \in [0, 1]$ usa la estrategia $E_1$ mientras
que la fracción $1 - x$ usa $E_2$, los pagos esperados de cada estrategia son
funciones lineales de $x$:

$$
\pi_1(x) = a x + b(1 - x), \qquad
\pi_2(x) = c x + d(1 - x) .
$$

El pago medio de toda la población es

$$
\bar\pi(x) = x\pi_1(x) + (1 - x)\pi_2(x) .
$$

### 1.2 Derivación de la ecuación del replicador por imitación

En lugar de suponer que los agentes calculan equilibrios de Nash de antemano,
la teoría de juegos evolutivos modela el **aprendizaje social por imitación**
(Taylor & Jonker, 1978; Schlag, 1998).

Supongamos que, con tasa constante por unidad de tiempo, cada agente elige al
azar un modelo $j$ de la población y compara su pago actual $\pi_i$ con el pago
$\pi_j$. Si $\pi_j > \pi_i$, el agente $i$ adopta la estrategia de $j$ con
probabilidad proporcional a la diferencia de pagos:

$$
\mathbb P(i \text{ imita a } j \mid \pi_j > \pi_i) = \kappa(\pi_j - \pi_i) ,
$$

donde $\kappa > 0$ fija la escala temporal. Un agente que usa $E_2$ (fracción
$1 - x$) encuentra a un agente que usa $E_1$ con probabilidad $x$, y a la
inversa. Como entre esos dos grupos siempre uno tiene mayor pago que el otro,
el flujo neto de agentes entre $E_2$ y $E_1$ da la tasa de cambio de $x$:

$$
\dot x = x(1 - x)\thinspace\kappa\bigl[\pi_1(x) - \pi_2(x)\bigr] .
$$

Reescalando el tiempo ($\kappa = 1$) y observando que $\pi_1 - \bar\pi =
\pi_1 - [x\pi_1 + (1-x)\pi_2] = (1 - x)(\pi_1 - \pi_2)$, obtenemos las dos
formas equivalentes de la **ecuación del replicador**
([C.I4](repaso/calculo.md#ci4)):

$$
\dot x = x\bigl[\pi_1(x) - \bar\pi(x)\bigr]
= x(1 - x)\bigl[\pi_1(x) - \pi_2(x)\bigr] .
$$

La primera expresión dice que la tasa de crecimiento relativa $\dot x / x$ de
una estrategia es igual al exceso de su pago sobre el promedio de la población.
La segunda muestra que la ventaja relativa entre estrategias es lineal en $x$:

$$
\Delta\pi(x) = \pi_1(x) - \pi_2(x) = (b - d) + \bigl[(a - c) - (b - d)\bigr]x .
$$

### 1.3 Los tres regímenes de los juegos de $2 \times 2$

Los estados homogéneos $x^\ast = 0$ (todos usan $E_2$) y $x^\ast = 1$ (todos
usan $E_1$) son **siempre puntos fijos** de la ecuación del replicador, porque
sin diversidad no hay estrategia alternativa que copiar. Un punto fijo interior
$x^\ast \in (0, 1)$ existe si y sólo si $\pi_1(x^\ast) = \pi_2(x^\ast)$, lo que
requiere que $a - c$ y $b - d$ tengan signos opuestos:

$$
x^\ast = \frac{d - b}{(a - c) + (d - b)} .
$$

Sea $f(x) = x(1 - x)\Delta\pi(x)$. Derivando en los tres posibles puntos fijos
([C.I4](repaso/calculo.md#ci4)):

$$
f'(0) = b - d, \qquad
f'(1) = c - a, \qquad
f'(x^\ast) = x^\ast(1 - x^\ast)\bigl[(a - c) - (b - d)\bigr] .
$$

Dependiendo de los signos de $a - c$ (ventaja de $E_1$ cuando todos usan $E_1$)
y $b - d$ (ventaja de $E_1$ cuando todos usan $E_2$), surgen tres clases
clásicas de interacción estratégica:

| Juego | Condición sobre pagos | Signo de $\Delta\pi(x)$ | Puntos fijos estables | Interpretación |
|---|---|---|---|---|
| **Dilema del prisionero** ($E_1 = C$, $E_2 = D$) | $c > a$ y $d > b$ ($T > R > P > S$) | $\Delta\pi(x) < 0$ para todo $x$ | $x^\ast = 0$ (todos traicionan) | Traicionar domina estrictamente; la cooperación se extingue ([sesión 18](s18_cooperacion_redes.md)) |
| **Coordinación** ($E_1 = A$, $E_2 = B$) | $a > c$ y $d > b$ | Negativo si $x < x^\ast$, positivo si $x > x^\ast$ | $x^\ast = 0$ y $x^\ast = 1$ (biestabilidad) | Convenciones o estándares tecnológicos; $x^\ast$ interior es un umbral (separatriz) |
| **Halcón-Paloma** ($E_1 = H$, $E_2 = P$) | $a < c$ y $b > d$ | Positivo si $x < x^\ast$, negativo si $x > x^\ast$ | $x^\ast \in (0, 1)$ único estable | Coexistencia polimórfica; cada estrategia invade cuando es rara |

**El juego Halcón-Paloma** (Maynard Smith & Price, 1973). Dos agentes compiten
por un recurso de valor $V > 0$. Un **Halcón** ($H$) siempre pelea; una
**Paloma** ($P$) amenaza pero se retira si el oponente pelea. Si dos halcones
se enfrentan, cada uno gana el recurso con probabilidad $1/2$ pero sufre un
costo de lesión $C/2$, con $C > V$; si dos palomas se encuentran, comparten el
recurso sin pelear ($V/2$ cada una); si un halcón encuentra a una paloma, el
halcón se lleva $V$ y la paloma $0$:

$$
A_{\text{HP}} = \begin{pmatrix} (V - C)/2 & V \cr 0 & V/2 \end{pmatrix} .
$$

Con una fracción $x$ de halcones, los pagos esperados son $\pi_H(x) =
x\frac{V-C}{2} + (1-x)V$ y $\pi_P(x) = (1-x)\frac{V}{2}$, cuya diferencia es

$$
\Delta\pi(x) = \pi_H(x) - \pi_P(x) = \frac{V - C x}{2} .
$$

La ecuación del replicador toma la forma explícita

$$
\dot x = \frac12\thinspace x(1 - x)(V - C x) .
$$

Como $C > V > 0$, existe un único equilibrio interior en

$$
x^\ast = \frac{V}{C} \in (0, 1) ,
$$

con $f'(x^\ast) = -\frac{C}{2}x^\ast(1 - x^\ast) = -\frac{V(C - V)}{2C} < 0$.
Por tanto, $x^\ast = V/C$ es **globalmente asintóticamente estable** en $(0,
1)$ (y constituye una estrategia evolutivamente estable mixta): cuanto más caro
es el conflicto ($C$ grande frente a $V$), menor es la fracción de agresivos que
puede sostener la población.

---

## 2. De Lotka-Volterra al ciclo distributivo de Goodwin

### 2.1 Derivación económica del modelo de Goodwin

Goodwin (1967) formuló en ecuaciones diferenciales la explicación de Marx sobre
el ciclo económico impulsado por el conflicto distributivo entre salarios y
ganancias, obteniendo exactamente las ecuaciones de depredador-presa de Lotka
(1925) y Volterra (1926).

Consideremos una economía cerrada con producción $Y$, capital $K$, empleo $L$ y
fuerza laboral disponible $N$ (todas funciones del tiempo $t$):

1. **Tecnología y crecimiento demográfico.** La relación capital-producto
   $\sigma = K/Y > 0$ es constante. La productividad del trabajo $a(t) = Y/L$ y
   la oferta laboral $N(t)$ crecen exponencialmente a tasas constantes
   $\alpha > 0$ y $\beta > 0$:
   $$
   \frac{\dot a}{a} = \alpha, \qquad \frac{\dot N}{N} = \beta .
   $$
2. **Las dos variables de estado.** Definimos la **tasa de empleo** $v$ y la
   **participación salarial** $u$ en el ingreso nacional (con salario real $w$):
   $$
   v = \frac{L}{N} \in (0, 1), \qquad
   u = \frac{w L}{Y} = \frac{w}{a} \in (0, 1) .
   $$
3. **Inversión y dinámica del empleo.** Los trabajadores consumen su salario y
   los capitalistas reinvierten todas las ganancias $\Pi = Y - wL = (1 - u)Y$.
   Si el capital se deprecia a tasa $\delta_K \ge 0$, entonces
   $\dot K = (1 - u)Y - \delta_K K$. Como $Y = K/\sigma$, la tasa de
   crecimiento del producto es $\frac{\dot Y}{Y} = \frac{\dot K}{K} =
   \frac{1 - u}{\sigma} - \delta_K$. Usando $v = \frac{Y}{a N}$ y la derivada
   logarítmica ([C.B3](repaso/calculo.md#cb3)),
   $$
   \frac{\dot v}{v} = \frac{\dot Y}{Y} - \frac{\dot a}{a} - \frac{\dot N}{N}
   = \underbrace{\left(\frac{1}{\sigma} - \delta_K - \alpha - \beta\right)}_{\alpha_G}
   - \underbrace{\frac{1}{\sigma}}_{\beta_G} u .
   $$
4. **Curva de Phillips real y dinámica salarial.** Cuando el empleo $v$ se
   acerca al pleno empleo, el poder de negociación de los trabajadores sube y
   el salario real crece a tasa lineal $\frac{\dot w}{w} = -\gamma_0 + \delta_G
   v$ (con $\gamma_0, \delta_G > 0$). Como $u = w/a$,
   $$
   \frac{\dot u}{u} = \frac{\dot w}{w} - \frac{\dot a}{a}
   = -\underbrace{(\gamma_0 + \alpha)}_{\gamma_G} + \delta_G v .
   $$

### 2.2 El sistema de Lotka-Volterra y su punto fijo

Reuniendo ambas ecuaciones (con parámetros positivos $\alpha_G, \beta_G,
\gamma_G, \delta_G > 0$), obtenemos el sistema de **Goodwin – Lotka-Volterra**
([C.A2](repaso/calculo.md#ca2)):

$$
\begin{aligned}
\dot v &= v\left(\alpha_G - \beta_G u\right), \cr
\dot u &= u\left(-\gamma_G + \delta_G v\right) .
\end{aligned}
$$

El diccionario entre ecología y macroeconomía es directo:

| Ecología (Lotka-Volterra) | Macroeconomía (Goodwin, 1967) | Mecanismo dinámico |
|---|---|---|
| Población de **presas** $v$ | **Tasa de empleo** $v = L/N$ | Si $u$ es baja, las ganancias $(1-u)Y$ son altas, se acumula capital y crece el empleo |
| Población de **depredadores** $u$ | **Participación salarial** $u = wL/Y$ | Si $v$ es alta, escasea la mano de obra, los salarios suben más rápido que la productividad y crece $u$ |
| Depredación $-\beta_G v u$ | Compresión de ganancias (*profit squeeze*) | Salarios altos reducen la reinversión y frenan la creación de empleo |
| Hambre del depredador $-\gamma_G u$ | Rezago salarial con alto desempleo | Con bajo empleo, la productividad crece más rápido que el salario real |

Las **ceroclinas** ($\dot v = 0$ y $\dot u = 0$) se cruzan en el origen
$(0, 0)$ y en un único **punto fijo interior**:

$$
v^\ast = \frac{\gamma_G}{\delta_G}, \qquad
u^\ast = \frac{\alpha_G}{\beta_G} .
$$

### 2.3 Linealización con el jacobiano y periodo del ciclo [Avanzado]

Sea $F(v, u) = \bigl(v(\alpha_G - \beta_G u),\ u(-\gamma_G + \delta_G v)\bigr)$.
La matriz jacobiana ([C.I3](repaso/calculo.md#ci3),
[L.A3](repaso/algebra_lineal.md#la3)) en cualquier punto $(v, u)$ es

$$
J(v, u) = \begin{pmatrix}
\alpha_G - \beta_G u & -\beta_G v \cr
\delta_G u & -\gamma_G + \delta_G v
\end{pmatrix} .
$$

- **En el origen $(0, 0)$:** $J(0, 0) = \operatorname{diag}(\alpha_G, -\gamma_G)$
  tiene eigenvalores $\lambda_1 = \alpha_G > 0$ y $\lambda_2 = -\gamma_G < 0$
  ([L.I1](repaso/algebra_lineal.md#li1), $\det J = -\alpha_G\gamma_G < 0$): el
  origen es un **punto silla** inestable.
- **En el equilibrio interior $(v^\ast, u^\ast)$:** Los términos diagonales se
  anulan porque $\alpha_G - \beta_G u^\ast = 0$ y $-\gamma_G + \delta_G v^\ast =
  0$, de modo que
  $$
  J(v^\ast, u^\ast) = \begin{pmatrix}
  0 & -\beta_G v^\ast \cr
  \delta_G u^\ast & 0
  \end{pmatrix}
  = \begin{pmatrix}
  0 & -\frac{\beta_G \gamma_G}{\delta_G} \cr
  \frac{\alpha_G \delta_G}{\beta_G} & 0
  \end{pmatrix} .
  $$

La traza y el determinante en $(v^\ast, u^\ast)$ son

$$
\operatorname{tr} J = 0, \qquad
\det J = \alpha_G \gamma_G > 0 .
$$

El polinomio característico $\lambda^2 + \alpha_G \gamma_G = 0$ tiene dos
**eigenvalores imaginarios puros**:

$$
\lambda_{1,2} = \pm i\omega, \qquad \omega = \sqrt{\alpha_G \gamma_G} .
$$

Por tanto, para oscilaciones pequeñas alrededor de $(v^\ast, u^\ast)$, las
desviaciones $(\delta v, \delta u) = (v - v^\ast, u - u^\ast)$ describen elipses
recorridas en **sentido antihorario** en el plano con $v$ en el eje horizontal
y $u$ en el vertical (cuando $v > v^\ast$, $\dot u > 0$; cuando $u > u^\ast$,
$\dot v < 0$). Si se dibuja $u$ en el eje horizontal, como en muchos textos
sobre Goodwin, el mismo ciclo se ve en sentido horario. El periodo es

$$
T_{\text{ciclo}} = \frac{2\pi}{\omega} = \frac{2\pi}{\sqrt{\alpha_G \gamma_G}} .
$$

La tasa de empleo $v(t)$ adelanta en un cuarto de ciclo ($\pi/2$ radianes) a la
participación salarial $u(t)$.

### 2.4 Cantidad conservada y órbitas cerradas exactas [Avanzado]

Como los eigenvalores lineales tienen parte real nula ($\operatorname{tr} J =
0$), el análisis lineal por sí solo no descarta que los términos no lineales
formen una espiral lenta. Para demostrar que las trayectorias son órbitas
cerradas exactas para cualquier amplitud en el primer cuadrante, separamos
variables dividiendo una ecuación entre la otra:

$$
\frac{du}{dv} = \frac{u(-\gamma_G + \delta_G v)}{v(\alpha_G - \beta_G u)}
\quad\Longrightarrow\quad
\left(\beta_G - \frac{\alpha_G}{u}\right) du +
\left(\delta_G - \frac{\gamma_G}{v}\right) dv = 0 .
$$

**Proposición.** La función separable

$$
V(v, u) = \left(\delta_G v - \gamma_G \ln v\right) +
\left(\beta_G u - \alpha_G \ln u\right)
$$

es una **integral primera** (constante de movimiento) del sistema de Goodwin en
$(v, u) \in (0, \infty)^2$, es decir, $\frac{d}{dt}V(v(t), u(t)) = 0$, y todas
las trayectorias distintas del punto fijo $(v^\ast, u^\ast)$ son curvas cerradas
periódicas.

*Demostración.* Derivando respecto al tiempo mediante la regla de la cadena:

$$
\dot V = \left(\delta_G - \frac{\gamma_G}{v}\right)\dot v +
\left(\beta_G - \frac{\alpha_G}{u}\right)\dot u
= (\delta_G v - \gamma_G)(\alpha_G - \beta_G u) +
(\beta_G u - \alpha_G)(-\gamma_G + \delta_G v) = 0 .
$$

Además, la función $g_1(v) = \delta_G v - \gamma_G \ln v$ es estrictamente
convexa ($g_1''(v) = \gamma_G / v^2 > 0$), alcanza su único mínimo global en
$v^\ast = \gamma_G / \delta_G$ y diverge a $+\infty$ cuando $v \to 0^+$ o
$v \to \infty$; lo mismo ocurre con $g_2(u) = \beta_G u - \alpha_G \ln u$ en
$u^\ast = \alpha_G / \beta_G$. Así, $V(v, u) = g_1(v) + g_2(u)$ es una función
estrictamente convexa con mínimo único en $(v^\ast, u^\ast)$, por lo que sus
curvas de nivel $V(v, u) = C > V(v^\ast, u^\ast)$ son curvas cerradas y
acotadas que rodean al equilibrio. $\square$

> **Nota numérica.** Si se integra el sistema de Goodwin con el método de Euler
> explícito ([C.I4](repaso/calculo.md#ci4)), $z_{n+1} = z_n + \Delta t F(z_n)$
> con $z = (v, u)$, el desarrollo de Taylor de segundo orden
> ([C.I1](repaso/calculo.md#ci1)) da
> $V(z_{n+1}) = V(z_n) + \Delta t\thinspace \nabla V \cdot F + \tfrac12 \Delta t^2 F^\top H F + \mathcal O(\Delta t^3)$,
> donde $\nabla V \cdot F = \dot V = 0$ y la hessiana
> $H = \operatorname{diag}(\gamma_G / v^2, \alpha_G / u^2)$ es definida positiva.
> Por tanto $V$ crece en cada paso: el error de truncamiento inyecta "energía"
> y la trayectoria se abre en una espiral artificial. En cambio, el **método de
> Euler simpléctico** (actualizar primero $v_{n+1}$ y usar ese nuevo valor para
> actualizar $u_{n+1}$) o el método de Runge-Kutta de cuarto orden mantienen
> $V$ prácticamente constante y las órbitas cerradas durante cientos de ciclos.

---

## 3. Resumen formal

| Modelo / Concepto | Resultado clave |
|---|---|
| Ecuación del replicador ($2 \times 2$) | $\dot x = x(\pi_1 - \bar\pi) = x(1-x)(\pi_1(x) - \pi_2(x))$ |
| Juego Halcón-Paloma ($C > V > 0$) | $\dot x = \frac12 x(1-x)(V - Cx)$; único punto fijo estable $x^\ast = V/C$ (coexistencia mixta) |
| Ciclo de Goodwin – Lotka-Volterra | $\dot v = v(\alpha_G - \beta_G u)$, $\dot u = u(-\gamma_G + \delta_G v)$; punto fijo $(v^\ast, u^\ast) = (\gamma_G/\delta_G, \alpha_G/\beta_G)$ |
| Estabilidad y periodo de Goodwin | $\operatorname{tr} J(v^\ast, u^\ast) = 0$, $\det J = \alpha_G\gamma_G > 0$; eigenvalores $\pm i\sqrt{\alpha_G\gamma_G}$, periodo $T = 2\pi/\sqrt{\alpha_G\gamma_G}$ |
| Integral primera de Goodwin | $V(v, u) = (\delta_G v - \gamma_G \ln v) + (\beta_G u - \alpha_G \ln u)$ cumple $\dot V = 0$ (órbitas cerradas) |

---

## 4. Ejercicios formales

**Ruta básica**

- **F1.** En el juego Halcón-Paloma con recurso $V = 4$ y costo de pelea
  $C = 10$, escribe la matriz de pagos $2 \times 2$, calcula la diferencia de
  pagos $\pi_H(x) - \pi_P(x)$ y encuentra la fracción de equilibrio $x^\ast$ de
  halcones. ¿Qué ocurre con $x(t)$ si la población empieza con $x_0 = 0.8$?
- **F2.** Considera un juego de coordinación con matriz de pagos $a = 3$,
  $b = 0$, $c = 0$, $d = 2$. Encuentra los tres puntos fijos de la ecuación del
  replicador $\dot x = x(1-x)\bigl[\pi_1(x) - \pi_2(x)\bigr]$ y determina hacia
  dónde converge la población si empieza en $x_0 = 0.35$ y si empieza en
  $x_0 = 0.45$.
- **F3.** En el modelo de Goodwin con parámetros $\alpha_G = 0.24$,
  $\beta_G = 0.40$, $\gamma_G = 0.36$ y $\delta_G = 0.40$: (a) calcula la tasa
  de empleo $v^\ast$ y la participación salarial $u^\ast$ de equilibrio; (b) si
  en un instante dado $v = 0.90$ y $u = 0.50$, determina los signos de $\dot v$
  y $\dot u$ e interpreta económicamente en qué fase del ciclo se encuentra la
  economía.

**Ruta avanzada**

- **F4.** Demuestra que si se transforma la fracción $x \in (0, 1)$ de la
  ecuación del replicador mediante el log-cociente de momios (*log-odds*)
  $y = \ln\bigl(\frac{x}{1-x}\bigr)$, entonces $\dot y = \pi_1(x) - \pi_2(x)$.
  Usa este resultado para resolver explícitamente $x(t)$ cuando la diferencia de
  pagos es constante, $\pi_C(x) - \pi_D(x) = -c < 0$ (juego de donación en
  población bien mezclada), y comprueba que $x(t)$ decae hacia $0$.
- **F5.** Calcula la matriz jacobiana del sistema de Goodwin evaluada en el
  equilibrio interior $(v^\ast, u^\ast)$ para los parámetros del ejercicio F3
  ($\alpha_G = 0.24$, $\gamma_G = 0.36$). Obtén sus eigenvalores y calcula el
  periodo aproximado $T_{\text{ciclo}}$ (en años) de las oscilaciones pequeñas.
- **F6.** Demuestra que para oscilaciones pequeñas $(\delta v, \delta u) =
  (v - v^\ast, u - u^\ast)$ alrededor del punto fijo de Goodwin, el desarrollo
  de Taylor de segundo orden ([C.I1](repaso/calculo.md#ci1)) de la cantidad
  conservada $V(v, u) - V(v^\ast, u^\ast)$ produce la elipse
  $$
  V(v, u) - V(v^\ast, u^\ast) \approx
  \frac{\delta_G^2}{2\gamma_G}(\delta v)^2 +
  \frac{\beta_G^2}{2\alpha_G}(\delta u)^2 .
  $$

---

## Referencias

- Goodwin, R. M. (1967). A growth cycle. En C. H. Feinstein (Ed.), *Socialism,
  Capitalism and Economic Growth: Essays Presented to Maurice Dobb* (pp. 54–58).
  Cambridge University Press.
- Hofbauer, J., & Sigmund, K. (1998). *Evolutionary Games and Population
  Dynamics*. Cambridge University Press.
  https://doi.org/10.1017/CBO9781139173179
- Lotka, A. J. (1925). *Elements of Physical Biology*. Williams & Wilkins.
- Maynard Smith, J., & Price, G. R. (1973). The logic of animal conflict.
  *Nature*, 246(5427), 15–18. https://doi.org/10.1038/246015a0
- Requião da Cunha, C. (2021). *Introduction to Econophysics: Contemporary
  Approaches with Python Simulations*. CRC Press.
  https://doi.org/10.1201/9781003127956
- Schecter, S., & Gintis, H. (2016). *Game Theory in Action: An Introduction to
  Classical and Evolutionary Game Theory*. Princeton University Press.
- Schlag, K. H. (1998). Why imitate, and if so, how? A boundedly rational
  approach to multi-armed bandits. *Journal of Economic Theory*, 78(1),
  130–156. https://doi.org/10.1006/jeth.1997.2347
- Strogatz, S. H. (2019). *Nonlinear Dynamics and Chaos* (2.ª ed.). CRC Press.
- Taylor, P. D., & Jonker, L. B. (1978). Evolutionary stable strategies and
  game dynamics. *Mathematical Biosciences*, 40(1–2), 145–156.
  https://doi.org/10.1016/0025-5564(78)90077-9
- Volterra, V. (1926). Fluctuations in the abundance of a species considered
  mathematically. *Nature*, 118(2972), 558–560.
  https://doi.org/10.1038/118558a0
