### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 34de0b78-ca60-11f1-8542-bf8eb1a504d5
begin
    # Usa el entorno del curso (Code/Project.toml) en lugar del gestor de
    # paquetes propio de Pluto: un solo entorno para todo el repositorio.
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, Statistics, Random, Distributions, Graphs, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

# ╔═╡ 34de0a74-ca60-11f1-9cf0-6136d3253a85
md"""
# Sesión 12: Juegos evolutivos y cooperación

**Objetivos de la sesión**
1. Usar una **ecuación diferencial ordinaria (EDO)** como regla de
   actualización paso a paso (método de Euler) e interpretar su retrato de fase.
2. Simular la **dinámica del replicador** a partir de la imitación de pagos y
   encontrar el equilibrio mixto del juego **Halcón-Paloma**.
3. Conectar el sistema depredador-presa de **Lotka-Volterra** con el **ciclo
   distributivo de Goodwin** entre tasa de empleo y participación salarial.
4. Simular el **dilema del prisionero en redes** y comparar la supervivencia de
   la cooperación en una **retícula 2D** frente a una red libre de escala de
   **Barabási-Albert**.
5. **[Avanzado]** Analizar la estabilidad de puntos fijos con el **jacobiano**
   y verificar la condición ``b/c > k`` de **reciprocidad de red**.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s12_juegos_evolutivos.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): C.B3 (la derivada como tasa
de cambio), C.B4 (función logística), L.B4 (matriz de adyacencia), C.I1
(desarrollo de Taylor), C.I3 (funciones de varias variables), C.I4 (EDO en una
dimensión y retrato de fase), L.I1 (eigenvalores). Ruta avanzada: C.A2
(sistemas de EDO en dos dimensiones), L.A3 (jacobiano y estabilidad lineal).
"""

# ╔═╡ 34de0aa6-ca60-11f1-b2bd-213ba88088be
TableOfContents(title = "Contenido")

# ╔═╡ 34de0ab0-ca60-11f1-a9d5-d16579bd03b2
md"""
## 0. El fenómeno: cooperación sin autoridad y ciclos sin choques

En la sesión 11 los agentes copiaban la opinión de un vecino sin preguntar qué
tan bien le iba. Hoy damos el siguiente paso: cada agente elige una
**estrategia**, recibe un **pago** según lo que hacen los demás y tiende a
imitar a quienes obtienen mejores resultados.

Esta idea conecta dos fenómenos clásicos de las ciencias sociales y la biología:

1. **Cooperación entre egoístas** (Axelrod, 1984; Nowak y May, 1992). En el
   **dilema del prisionero**, traicionar siempre deja un pago individual mayor
   que cooperar, sin importar qué haga el otro. Si todos interactúan al azar con
   todos, los traidores explotan a los cooperadores hasta extinguirlos. Pero si
   los agentes interactúan únicamente con sus vecinos en el espacio o en una
   red, los cooperadores forman **agrupaciones compactas** (*clústeres*) que se
   apoyan entre sí y resisten la invasión.
2. **Ciclos económicos sin choques externos** (Goodwin, 1967). Cuando el empleo
   es alto, los salarios suben y comprimen las ganancias; con menos ganancias,
   cae la inversión y sube el desempleo; el desempleo modera los salarios,
   recupera las ganancias y relanza el empleo. La economía oscila perpetuamente
   con la misma matemática que las poblaciones de presas y depredadores de
   Lotka y Volterra.
"""

# ╔═╡ 34de0aba-ca60-11f1-80cf-7f5030690333
md"""
Las dos poblaciones empezaron con el mismo 50% de cooperadores (azul) y 50% de
traidores (naranja), y con la misma tentación a traicionar (``b = 1.35``). Al
emparejarse al azar, la cooperación desaparece; en la retícula 2D, sobrevive en
islas compactas de cooperadores.
"""

# ╔═╡ 34de0ac4-ca60-11f1-b4d1-3b65cc417c97
md"""
## 1. Repaso: una ecuación diferencial como regla de actualización

Una **ecuación diferencial ordinaria (EDO)** ``\dot{x} = dx/dt = f(x)`` dice a
qué velocidad cambia ``x`` cuando su valor actual es ``x`` (repaso C.B3 y
C.I4). En la computadora la leemos directamente como una **regla de
actualización** en pasos pequeños ``\Delta t`` (**método de Euler**):

```math
x_{t + \Delta t} = x_t + f(x_t)\,\Delta t .
```

Sin necesidad de resolver integrales, el signo de ``f(x)`` nos dice todo el
comportamiento cualitativo (**retrato de fase**):
- Donde ``f(x) > 0``, ``x`` **crece** (se mueve a la derecha).
- Donde ``f(x) < 0``, ``x`` **decrece** (se mueve a la izquierda).
- Donde ``f(x^\ast) = 0``, ``x^\ast`` es un **punto fijo** (equilibrio): es
  **estable** si la curva cruza el cero bajando (``f'(x^\ast) < 0``) e
  **inestable** si cruza subiendo (``f'(x^\ast) > 0``).
"""

# ╔═╡ 34de0acc-ca60-11f1-aebb-a3da9de9ac44
# Un paso del método de Euler: nuevo valor = valor actual + tasa de cambio × paso de tiempo
paso_euler(f, x, Δt) = x + f(x) * Δt

# ╔═╡ 34de0acc-ca60-11f1-a6f2-25a850deb231
# Si f(x) = x(1 - x)(0.6 - x), el punto fijo x* = 0.6 atrae desde ambos lados
paso_euler(x -> x * (1 - x) * (0.6 - x), 0.2, 0.1),
paso_euler(x -> x * (1 - x) * (0.6 - x), 0.8, 0.1)

# ╔═╡ 34de0ad8-ca60-11f1-bb4c-3bca4d55456b
md"""
## 2. Dinámica del replicador y el juego Halcón-Paloma

Supongamos una población grande donde una fracción ``x \in [0, 1]`` usa la
estrategia 1 y ``1 - x`` usa la estrategia 2. Al emparejarse al azar, sus pagos
esperados son ``\pi_1(x)`` y ``\pi_2(x)``, y el pago promedio de la población
es ``\bar{\pi}(x) = x\pi_1(x) + (1-x)\pi_2(x)``.

**De la imitación a la EDO.** Si cada agente se compara al azar con otro y,
cuando el otro gana más, lo imita con probabilidad proporcional a la diferencia
de pagos, la fracción ``x`` obedece la **ecuación del replicador** (Taylor y
Jonker, 1978; repaso C.I4):

```math
\dot{x} = x\bigl[\pi_1(x) - \bar{\pi}(x)\bigr] = x(1 - x)\bigl[\pi_1(x) - \pi_2(x)\bigr] .
```

El factor ``x(1 - x)`` refleja que para que haya conversiones deben encontrarse
las dos estrategias; el factor ``\pi_1(x) - \pi_2(x)`` dice hacia cuál fluye la
población.

### El juego Halcón-Paloma (Maynard Smith y Price, 1973)

Dos agentes compiten por un recurso de valor ``V > 0`` (un mercado, un
territorio). Cada uno puede actuar como **Halcón** (agresivo: pelea siempre) o
como **Paloma** (pacífico: comparte con otra paloma y se retira ante un
halcón). Si dos halcones pelean, se reparten el recurso pero sufren un costo de
conflicto ``C > V``:

| Jugador fila \ columna | Halcón (fracción ``x``) | Paloma (fracción ``1 - x``) |
|---|---|---|
| **Halcón** | ``(V - C) / 2`` | ``V`` |
| **Paloma** | ``0`` | ``V / 2`` |

Los pagos esperados y su diferencia son:

```math
\pi_H(x) = x\,\frac{V - C}{2} + (1 - x)V, \qquad
\pi_P(x) = (1 - x)\frac{V}{2}
\quad\Longrightarrow\quad
\pi_H(x) - \pi_P(x) = \frac{V - C x}{2} .
```

Sustituyendo en la ecuación del replicador:

```math
\dot{x} = \frac{1}{2}\,x(1 - x)(V - C x) .
```

- Si casi todos son palomas (``x`` pequeño), ser halcón paga ``V > V/2`` y los
  halcones invaden (``\dot{x} > 0``).
- Si casi todos son halcones (``x`` grande), las peleas cuestan ``(V-C)/2 < 0``
  y conviene ser paloma (``\dot{x} < 0``).
- Las dos estrategias coexisten en el **equilibrio interior estable**
  ``x^\ast = V / C``.

### Laboratorio: Halcón-Paloma por EDO y por simulación de agentes

Fijamos el valor del recurso en ``V = 2``. Mueve el costo de pelear ``C > V`` y
la fracción inicial ``x_0`` de halcones:

Costo del conflicto ``C`` (con ``V = 2``): $(@bind C_hp Slider([2.5, 3.0, 4.0, 5.0, 8.0], default = 4.0, show_value = true))
  Fracción inicial de halcones ``x_0``: $(@bind x0_hp Slider(0.1:0.1:0.9, default = 0.1, show_value = true))
"""

# ╔═╡ 34de0aec-ca60-11f1-b206-c7b11fa25adf
md"""
Observa dos resultados clave:
1. **La ecuación diferencial es el promedio del modelo de agentes:** la
   simulación microscópica de 250 agentes que se comparan en pares (curva
   naranja) sigue de cerca la solución de Euler de la ecuación del replicador
   (curva azul) y fluctúa alrededor de ``x^\ast = V / C``.
2. **Efecto del costo del conflicto:** cuando pelear se vuelve más destructivo
   (``C`` sube de 2.5 a 8), el equilibrio ``x^\ast = V/C`` baja de 80% a 25% de
   halcones. Ningún halcón ni paloma cambia de personalidad: cambia la
   **proporción** en la población porque los halcones se dañan entre sí.
"""

# ╔═╡ 34de0af6-ca60-11f1-9944-2744d9670242
md"""
## 3. Lotka-Volterra y el ciclo distributivo de Goodwin

En un juego de ``2 \times 2`` sólo hay una variable independiente ``x(t)``, y
en una dimensión una trayectoria continua no puede oscilar (siempre se detiene
en un punto fijo). En cuanto dos variables se retroalimentan con signos
opuestos, pueden aparecer **ciclos endógenos**.

Goodwin (1967) mostró que la interacción entre la **tasa de empleo**
``v = L / N`` y la **participación de los salarios en el ingreso nacional**
``u = wL / Y`` tiene exactamente la estructura del sistema depredador-presa de
**Lotka-Volterra** (repaso C.A2):

```math
\begin{aligned}
\dot{v} &= v(\alpha_G - \beta_G u) \qquad &&\text{(la inversión de ganancias crea empleo: } v \text{ hace el papel de la presa)}, \\
\dot{u} &= u(-\gamma_G + \delta_G v) \qquad &&\text{(el alto empleo empuja los salarios: } u \text{ hace el papel del depredador)}.
\end{aligned}
```

**El mecanismo económico en cuatro pasos:**
1. **Auge de inversión** (``u < u^\ast``, ``v`` subiendo): los salarios absorben
   poco del producto, las ganancias ``(1 - u)Y`` son altas, las empresas
   invierten y el empleo ``v`` crece.
2. **Compresión de ganancias** (*profit squeeze*, ``v > v^\ast``, ``u``
   subiendo): con casi pleno empleo, los trabajadores negocian aumentos
   salariales por encima del crecimiento de la productividad y sube ``u``.
3. **Recesión** (``u > u^\ast``, ``v`` bajando): la alta participación salarial
   reduce las ganancias y la inversión; la producción crece menos que la fuerza
   laboral y cae el empleo ``v``.
4. **Recuperación de márgenes** (``v < v^\ast``, ``u`` bajando): el desempleo
   frena los salarios reales mientras la productividad sigue creciendo; ``u``
   baja y el ciclo vuelve a empezar.

El único punto fijo interior está donde se anulan ambos paréntesis:

```math
v^\ast = \frac{\gamma_G}{\delta_G}, \qquad u^\ast = \frac{\alpha_G}{\beta_G} .
```

### Laboratorio: el ciclo de Goodwin en el tiempo y en el plano de fase

Con ``\alpha_G = 0.24``, ``\gamma_G = 0.36`` y ``\delta_G = 0.40``, el empleo de
equilibrio es ``v^\ast = 0.90`` (90%). Mueve la tasa de empleo inicial ``v_0``
y la sensibilidad de la inversión ``\beta_G``:

Empleo inicial ``v_0``: $(@bind v0_gw Slider([0.82, 0.85, 0.88, 0.94, 0.96], default = 0.85, show_value = true))
  Sensibilidad a los salarios ``\beta_G``: $(@bind β_gw Slider([0.35, 0.40, 0.45], default = 0.40, show_value = true))
"""

# ╔═╡ 34de0b0a-ca60-11f1-bca5-33dffc73fc7b
md"""
En la gráfica izquierda, el pico de la tasa de empleo ``v(t)`` (azul) ocurre
siempre **un cuarto de ciclo antes** que el pico de la participación salarial
``u(t)`` (naranja). En la gráfica derecha, la economía gira en sentido
**antihorario** alrededor del punto fijo ``(v^\ast, u^\ast)`` sin caer en él
ni alejarse: con ``v > v^\ast`` la participación salarial sube, y con
``u > u^\ast`` el empleo baja. Cuanto más lejos del equilibrio empieza ``v_0``,
más amplia es la órbita. (Muchos textos dibujan ``u`` en el eje horizontal y
``v`` en el vertical; con esos ejes el mismo ciclo se ve en sentido horario.)
"""

# ╔═╡ 34de0b0a-ca60-11f1-9afb-53e07cdd8b41
md"""
## 4. Dilema del prisionero en redes: retícula frente a Barabási-Albert

Volvamos a la paradoja de la sección 0. En el **dilema del prisionero débil**
de Nowak y May (1992), cada encuentro entre dos agentes da los siguientes
pagos:
- Si ambos **cooperan** (``C``), cada uno recibe ``R = 1``.
- Si un **traidor** (``D``) se enfrenta a un cooperador, el traidor recibe la
  **tentación** ``T = b > 1`` y el cooperador recibe ``S = 0``.
- Si ambos **traicionan**, cada uno recibe ``P = 0``.

En una población bien mezclada, un traidor gana ``b x`` y un cooperador gana
``x``; como ``b > 1``, ``\dot{x} = (1 - b)x^2(1 - x) < 0`` y los cooperadores
siempre se extinguen.

**¿Qué cambia en una red?** En cada ronda:
1. Cada nodo ``i`` juega con todos sus vecinos en la red y suma sus pagos
   ``\Pi_i``.
2. Cada nodo imita a vecinos con mayor pago. En la retícula usamos la regla
   de Nowak y May: adoptar la estrategia del mayor pago en el vecindario
   (incluido el propio). En la red de Barabási-Albert usamos la regla de Santos
   y Pacheco (2005): elegir un vecino ``j`` al azar y, si ``\Pi_j > \Pi_i``,
   copiar su estrategia con probabilidad
   ``(\Pi_j - \Pi_i) / (b \max(k_i, k_j))``.

En ambos casos cada agente juega sólo con sus vecinos, no consigo mismo.

### Laboratorio: clústeres en la retícula 2D y concentradores en Barabási-Albert

Mueve la tentación a traicionar ``b > 1`` y observa primero cómo evoluciona la
retícula 2D de ``36 \times 36`` agentes:

Tentación de traicionar ``b``: $(@bind b_pd Slider([1.05, 1.15, 1.25, 1.35, 1.55, 1.75, 1.95], default = 1.35, show_value = true))
"""

# ╔═╡ 34de0b1e-ca60-11f1-a4f6-992b2ffc5b2c
md"""
Al inicio, los traidores invaden a los cooperadores aislados y la cooperación
cae en las primeras generaciones. Pero en la retícula 2D, los cooperadores que
quedan juntos forman **bloques compactos**: un cooperador en el interior o en la
frontera plana de un bloque recibe el apoyo de muchos vecinos cooperadores,
mientras que un traidor en la orilla sólo puede explotar a los pocos
cooperadores que toca.

¿Cómo depende la cooperación de la estructura de la red en todo el rango de
``b``? Comparamos tres topologías con ``N = 400`` agentes y el mismo grado
medio (``\langle k \rangle = 8``). La dinámica en la red de Barabási-Albert es
lenta, así que la dejamos correr hasta 5 000 generaciones (o hasta que toda la
red coopere o traicione):
1. **Mezcla al azar** (en cada ronda se sortean 8 rivales al azar);
2. **Retícula 2D** (8 vecinos de Moore fijos);
3. **Red de Barabási-Albert** (libre de escala, ``m = 4 \Rightarrow \langle k \rangle \approx 8``, Santos y Pacheco, 2005).
"""

# ╔═╡ 34de0b32-ca60-11f1-b8d1-0d9743fd2efa
md"""
La comparación revela dos mecanismos distintos de cooperación en redes:
- **En la retícula 2D (azul)**, la cooperación sobrevive por **reciprocidad
  espacial** para ``b`` moderado, cae por escalones y desaparece hacia
  ``b \approx 1.6``. Los escalones sólo pueden estar en cocientes ``b = m/n`` de
  enteros entre 1 y 8 (números de vecinos cooperadores), como ``4/3`` o
  ``3/2``, porque ahí cambia qué comparación de pagos gana. A partir de
  ``b = 5/3`` un traidor que toca a 3 cooperadores (pago ``3b``) ya gana más que
  el cooperador del borde plano de un bloque de ``3 \times 3`` (pago 5), y los
  bloques dejan de crecer por sus caras.
- **En la red de Barabási-Albert (naranja)**, la cooperación resiste valores de
  ``b`` mucho mayores (Santos y Pacheco, 2005). La razón está en los
  **concentradores (*hubs*)**: un *hub* cooperador rodeado de cooperadores suma
  decenas de pagos ``R = 1`` y se vuelve el nodo más exitoso de la red,
  convirtiendo a sus vecinos; en cambio, si un *hub* es traidor, sus vecinos
  copian su traición, dejan de cooperar con él y el pago del *hub* traidor se
  desploma a cero.
"""

# ╔═╡ 34de0b32-ca60-11f1-829f-cbff897663c8
md"""
## 5. [Avanzado] Estabilidad con el jacobiano y la regla ``b/c > k``

### 5.1 El jacobiano del ciclo de Goodwin y su cantidad conservada

Linealizando el sistema de Goodwin alrededor del punto fijo interior
``(v^\ast, u^\ast) = (\gamma_G/\delta_G,\, \alpha_G/\beta_G)``, la matriz
**jacobiana** (repaso L.A3) es:

```math
J(v^\ast, u^\ast) =
\begin{pmatrix}
\dfrac{\partial \dot{v}}{\partial v} & \dfrac{\partial \dot{v}}{\partial u} \\[8pt]
\dfrac{\partial \dot{u}}{\partial v} & \dfrac{\partial \dot{u}}{\partial u}
\end{pmatrix}_{(v^\ast, u^\ast)}
=
\begin{pmatrix}
0 & -\beta_G v^\ast \\
\delta_G u^\ast & 0
\end{pmatrix} .
```

Su traza es ``\mathrm{tr}(J) = 0`` y su determinante es
``\det(J) = \beta_G \delta_G v^\ast u^\ast = \alpha_G \gamma_G > 0``. Por tanto,
los eigenvalores son **imaginarios puros** (polinomio característico
``\lambda^2 + \alpha_G \gamma_G = 0``):

```math
\lambda_{1,2} = \pm i\omega, \qquad \omega = \sqrt{\alpha_G \gamma_G},
\qquad T_{\text{ciclo}} = \frac{2\pi}{\sqrt{\alpha_G \gamma_G}} .
```

Con los parámetros de la sección 3 (``\alpha_G = 0.24``, ``\gamma_G = 0.36``),
el periodo teórico para oscilaciones pequeñas es
``T_{\text{ciclo}} = 2\pi / \sqrt{0.0864} \approx 21.4`` años.

Con parte real nula, la linealización **no basta** para decidir la
estabilidad: los términos no lineales podrían convertir el centro lineal en una
espiral lenta. Lo que garantiza las órbitas cerradas es que el sistema conserva
exactamente la función tipo energía (repaso C.A2):

```math
V(v, u) = (\delta_G v - \gamma_G \ln v) + (\beta_G u - \alpha_G \ln u),
\qquad \frac{dV}{dt} = 0 ,
```

cuyas curvas de nivel cerradas son las órbitas del plano de fase.

### 5.2 La condición ``b/c > k`` de reciprocidad de red

En el **juego de donación**, cooperar cuesta ``c > 0`` por vecino y otorga un
beneficio ``b > c`` al vecino. Ohtsuki, Hauert, Lieberman y Nowak (2006)
demostraron que bajo selección débil (aptitud ``f_i = 1 + w\pi_i`` con
``w \ll 1``) y actualización **muerte-nacimiento** (un nodo al azar es
reemplazado por uno de sus ``k`` vecinos proporcionalmente a su aptitud), la
probabilidad ``\rho_C`` de que un solo mutante cooperador tome toda la red
supera al azar neutro ``1/N`` cuando la razón beneficio-costo supera el grado
``k`` de la red:

```math
\frac{b}{c} > k .
```

La regla vale para redes grandes (``N \gg k``). En una red regular finita el
umbral exacto es algo mayor (Taylor, Day y Wild, 2007; Chen, 2013):

```math
\left(\frac{b}{c}\right)^{\ast} = \frac{N - 2}{N/k - 2}
\;\xrightarrow{\;N \to \infty\;}\; k .
```

La gráfica muestra ``N \cdot \rho_C`` simulado (30 000 ensayos por punto, con
``N = 36`` y ``w = 0.03``) en un anillo (``k = 2``, umbral finito
``34/16 \approx 2.1``) y en una red regular aleatoria de grado ``k = 4``
(umbral finito ``34/7 \approx 4.9``). Cada curva cruza la deriva neutra
``N \rho_C = 1`` cerca de su umbral finito, un poco a la derecha de
``b/c = k``.
"""

# ╔═╡ 34de0b46-ca60-11f1-9959-85a5afc38314
md"""
## Resumen

| Modelo | Ecuación / Regla microscópica | Comportamiento colectivo |
|---|---|---|
| **Dinámica del replicador (Halcón-Paloma)** | ``\dot{x} = x(1-x)(\pi_H - \pi_P) = \tfrac{1}{2}x(1-x)(V - Cx)`` | Equilibrio interior estable ``x^\ast = V/C``: coexistencia de agresivos y pacíficos |
| **Goodwin – Lotka-Volterra** | ``\dot{v} = v(\alpha_G - \beta_G u)``, ``\dot{u} = u(-\gamma_G + \delta_G v)`` | Órbitas cerradas alrededor de ``(\gamma_G/\delta_G,\, \alpha_G/\beta_G)``; el empleo adelanta a los salarios |
| **Dilema del prisionero en redes** | Cooperar (``R=1, S=0``) vs. traicionar (``T=b>1, P=0``) con vecinos | En mezcla global ``x \to 0``; en retícula 2D sobreviven clústeres; en Barabási-Albert los *hubs* sostienen la cooperación; bajo selección débil emerge si ``b/c > k`` |
"""

# ╔═╡ 34de0b50-ca60-11f1-a089-3fbb5729b76a
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `paso_replicador(x, π1, π2, Δt)`, que devuelva el nuevo valor de
la fracción `x` después de un paso de tiempo `Δt` con el método de Euler
aplicado a la ecuación del replicador ``\dot{x} = x(1 - x)(\pi_1 - \pi_2)``.
"""

# ╔═╡ 34de0b50-ca60-11f1-8145-d54bf95b14cc
function paso_replicador(x, π1, π2, Δt)
    missing  # reemplaza esta línea
end

# ╔═╡ 34de0b5a-ca60-11f1-8997-e1325165499b
md"""
**2.** Escribe `equilibrio_goodwin(α, β, γ, δ)`, que devuelva la tupla
`(v_est, u_est)` con la tasa de empleo ``v^\ast`` y la participación salarial
``u^\ast`` en el punto fijo interior del modelo de Goodwin
``\dot{v} = v(\alpha - \beta u)``, ``\dot{u} = u(-\gamma + \delta v)``.
"""

# ╔═╡ 34de0b5a-ca60-11f1-80b1-d1e5f47953a2
function equilibrio_goodwin(α, β, γ, δ)
    missing  # reemplaza esta línea
end

# ╔═╡ 34de0b64-ca60-11f1-a56d-f59833580b58
md"""
### Ruta básica: interpretación

3. En el juego Halcón-Paloma con ``V = 2``, explica por qué la fracción de
   halcones de equilibrio disminuye cuando aumenta el costo de pelear ``C``, y
   por qué nunca llega a cero mientras ``C`` sea finito.
4. En el laboratorio de Halcón-Paloma, si la población empieza con ``x_0 = 0``
   (100% palomas), ¿qué predice exactamente la ecuación diferencial ``\dot{x}``?
   ¿Qué ocurriría en cambio si aparece un solo halcón por mutación?
5. En el ciclo de Goodwin, explica con intuición económica por qué el máximo de
   la tasa de empleo ``v(t)`` ocurre **antes** que el máximo de la participación
   salarial ``u(t)``, y no al mismo tiempo.
6. ¿Por qué un cooperador aislado rodeado de traidores en la retícula 2D
   desaparece en una sola generación, mientras que un bloque cuadrado de
   ``3 \times 3`` cooperadores logra sobrevivir cuando ``b = 1.35``?
7. Explica por qué un concentrador (*hub*) traidor en una red de
   Barabási-Albert termina perdiendo su alto pago inicial cuando sus vecinos lo
   imitan, mientras que un *hub* cooperador refuerza su propio pago.

### Ruta avanzada

8. Verifica que para la matriz jacobiana ``J(v^\ast, u^\ast)`` del modelo de
   Goodwin, ``\mathrm{tr}(J) = 0`` y ``\det(J) = \alpha_G \gamma_G``. Calcula el
   periodo de las oscilaciones pequeñas cuando ``\alpha_G = 0.24`` y
   ``\gamma_G = 0.36``.
9. Deriva respecto al tiempo la función
   ``V(v, u) = \delta_G v - \gamma_G \ln v + \beta_G u - \alpha_G \ln u`` usando
   la regla de la cadena y demuestra que ``\dot{V} = 0`` a lo largo de cualquier
   trayectoria de Goodwin.
10. En una cuadrícula con vecindad de Moore (8 vecinos), calcula para qué valor
    crítico de ``b`` un traidor en la orilla plana de un bloque de ``3 \times 3``
    cooperadores (que toca a 3 cooperadores) obtiene un pago mayor que un
    cooperador del borde del bloque (que toca a 5 cooperadores). Compáralo con
    la caída de la curva azul en la sección 4.
11. Según la condición ``b/c > k`` de Ohtsuki et al. (2006), ¿en cuál red es más
    fácil que prospere un mutante cooperador bajo selección débil: en un anillo
    (``k = 2``), en una cuadrícula (``k = 4``) o en una red aleatoria densa con
    ``k = 20``? Explica intuitivamente por qué tener demasiados vecinos dificulta
    la reciprocidad local.

## Para la próxima sesión

- Lectura: Requião da Cunha (2021), capítulo 7; Slanina (2014), sección 8.1;
  Nowak y May (1992). Para la ruta avanzada, Goodwin (1967); Ohtsuki et al.
  (2006); Santos y Pacheco (2005). De consulta: Strogatz (2019) y Schecter y
  Gintis (2016).
- En la semana 13 estudiaremos qué pasa cuando ganar depende de estar en la
  **minoría**: el problema del bar El Farol y el **juego de la minoría** como
  modelo de eficiencia y volatilidad en un mercado financiero.
"""

# ╔═╡ 34de0b70-ca60-11f1-9fcd-ef8d2e035d94
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 34de0b78-ca60-11f1-afda-e9cd8a803efe
begin
    const AZUL       = "#2a78d6"   # serie principal, cooperadores / halcones / empleo
    const NARANJA    = "#eb6834"   # segunda serie o teoría, traidores / salarios
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const BLANCO     = "#ffffff"
end;

# ╔═╡ 34de0b82-ca60-11f1-a60f-79fdc16f43b4
"""
Simula el juego Halcón-Paloma con recurso `V` y costo `C` desde `x0`:
1) Integración de Euler de la EDO del replicador;
2) Simulación microscópica de `N` agentes que se emparejan e imitan según la diferencia de pagos.
"""
function simular_halcon_paloma(V, C, x0; N = 250, T = 18.0, Δt = 0.05, rng = Xoshiro(42))
    pasos = round(Int, T / Δt)
    t_edo = collect(0:pasos) .* Δt
    f(x) = 0.5 * x * (1 - x) * (V - C * x)

    x_edo = zeros(pasos + 1)
    x_edo[1] = x0
    x0_alt = x0 < 0.5 ? 0.92 : 0.08
    x_alt = zeros(pasos + 1)
    x_alt[1] = x0_alt
    for k in 1:pasos
        x_edo[k+1] = x_edo[k] + f(x_edo[k]) * Δt
        x_alt[k+1] = x_alt[k] + f(x_alt[k]) * Δt
    end

    # Simulación basada en agentes: imitación proporcional a la diferencia de pagos
    k_h = clamp(round(Int, x0 * N), 1, N - 1)
    escala = max(V, (C - V) / 2)
    pasos_por_unidad = round(Int, N * escala)
    pasos_abm = round(Int, T * pasos_por_unidad)
    cada = max(1, pasos_por_unidad ÷ 4)
    t_abm = Float64[0.0]
    x_abm = Float64[k_h / N]
    for p in 1:pasos_abm
        x = k_h / N
        π_H = x * (V - C) / 2 + (1 - x) * V
        π_P = (1 - x) * V / 2
        if rand(rng) < 2 * x * (1 - x)
            Δπ = π_H - π_P
            if rand(rng) < abs(Δπ) / (2 * escala)
                k_h += Δπ > 0 ? 1 : -1
                k_h = clamp(k_h, 1, N - 1)
            end
        end
        if p % cada == 0
            push!(t_abm, p / pasos_por_unidad)
            push!(x_abm, k_h / N)
        end
    end
    (t_edo = t_edo, x_edo = x_edo, x_alt = x_alt, x0_alt = x0_alt, t_abm = t_abm, x_abm = x_abm)
end

# ╔═╡ 34de0ae2-ca60-11f1-b8fb-496d030b090e
sim_hp = simular_halcon_paloma(2.0, C_hp, x0_hp; N = 250, T = 18.0, rng = Xoshiro(42));

# ╔═╡ 34de0ae2-ca60-11f1-ba00-55d677ebabd8
let
    V = 2.0
    x_est = V / C_hp
    xs = range(0.0, 1.0, length = 200)
    f_hp(x) = 0.5 * x * (1 - x) * (V - C_hp * x)

    p1 = plot(xs, f_hp.(xs); color = AZUL, lw = 2.4, label = "ẋ = ½ x(1−x)(V − Cx)",
              xlabel = "fracción de halcones x", ylabel = "tasa de cambio ẋ",
              legend = :outerbottom, title = "Retrato de fase (equilibrio x* = V/C = $(round(x_est, digits=2)))")
    hline!(p1, [0.0]; color = GRIS, ls = :dot, lw = 1.2, label = "")
    scatter!(p1, [0.0, 1.0], [0.0, 0.0]; color = BLANCO, msc = GRIS, ms = 5, msw = 1.5,
             label = "Inestables (0 y 1)")
    scatter!(p1, [x_est], [0.0]; color = NARANJA, ms = 6, msw = 0,
             label = "Estable (x* = $(round(x_est, digits=2)))")

    p2 = plot(sim_hp.t_edo, sim_hp.x_edo; color = AZUL, lw = 2.5,
              label = "Replicador (Euler, desde x₀ = $x0_hp)",
              xlabel = "tiempo", ylabel = "fracción de halcones x(t)",
              ylims = (-0.02, 1.02), legend = :outerbottom,
              title = "EDO frente a N = 250 agentes que imitan pagos")
    plot!(p2, sim_hp.t_edo, sim_hp.x_alt; color = AZUL, lw = 1.4, ls = :dash,
          label = "Replicador (desde x₀ = $(sim_hp.x0_alt))")
    plot!(p2, sim_hp.t_abm, sim_hp.x_abm; color = NARANJA, lw = 1.5, alpha = 0.85,
          label = "Simulación de agentes (N = 250)")
    hline!(p2, [x_est]; color = GRIS, ls = :dot, lw = 1.5, label = "Equilibrio V / C")

    plot(p1, p2; layout = (1, 2), size = (900, 390), bottom_margin = 4Plots.mm)
end

# ╔═╡ 34de0b8c-ca60-11f1-9669-f945d2d90a72
"""
Integra el sistema de Goodwin – Lotka-Volterra con RK4: a diferencia de Euler, la
cantidad conservada V(v, u) casi no se desvía y las órbitas se ven cerradas.
"""
function integrar_goodwin(α, β, γ, δ, v0, u0; T = 45.0, Δt = 0.04)
    pasos = round(Int, T / Δt)
    ts = collect(0:pasos) .* Δt
    vs = zeros(pasos + 1)
    us = zeros(pasos + 1)
    vs[1], us[1] = v0, u0
    f(v, u) = (v * (α - β * u), u * (-γ + δ * v))
    for k in 1:pasos
        v, u = vs[k], us[k]
        k1v, k1u = f(v, u)
        k2v, k2u = f(v + 0.5Δt * k1v, u + 0.5Δt * k1u)
        k3v, k3u = f(v + 0.5Δt * k2v, u + 0.5Δt * k2u)
        k4v, k4u = f(v + Δt * k3v, u + Δt * k3u)
        vs[k+1] = v + (Δt / 6) * (k1v + 2k2v + 2k3v + k4v)
        us[k+1] = u + (Δt / 6) * (k1u + 2k2u + 2k3u + k4u)
    end
    (t = ts, v = vs, u = us)
end

# ╔═╡ 34de0b8c-ca60-11f1-bdd9-cb49ddb09e55
"Simula la órbita principal y una órbita de referencia del ciclo de Goodwin."
function simular_goodwin(α, β, γ, δ; v0 = 0.85, u0 = α / β, T = 45.0)
    orb = integrar_goodwin(α, β, γ, δ, v0, u0; T = T)
    ref = integrar_goodwin(α, β, γ, δ, 0.88, α / β; T = T)
    (t = orb.t, v = orb.v, u = orb.u, v_ref = ref.v, u_ref = ref.u,
     v_est = γ / δ, u_est = α / β)
end

# ╔═╡ 34de0b00-ca60-11f1-a9fb-498261e3d40e
sim_gw = simular_goodwin(0.24, β_gw, 0.36, 0.40; v0 = v0_gw, u0 = 0.24 / β_gw, T = 45.0);

# ╔═╡ 34de0b00-ca60-11f1-9586-271aca99373c
let
    g = sim_gw
    p1 = plot(g.t, 100 .* g.v; color = AZUL, lw = 2.3, label = "Tasa de empleo v(t)",
              xlabel = "tiempo (años)", ylabel = "porcentaje (%)",
              ylims = (42, 102), legend = :outerbottom, legend_columns = 2,
              title = "Series de tiempo: el empleo adelanta a los salarios")
    plot!(p1, g.t, 100 .* g.u; color = NARANJA, lw = 2.3, ls = :dash,
          label = "Participación salarial u(t)")
    hline!(p1, [100g.v_est]; color = AZUL, alpha = 0.35, ls = :dot, label = "")
    hline!(p1, [100g.u_est]; color = NARANJA, alpha = 0.35, ls = :dot, label = "")

    p2 = plot(100 .* g.v_ref, 100 .* g.u_ref; color = GRIS_CLARO, lw = 1.4, ls = :dash,
              label = "Referencia (v₀ = 88%)",
              xlabel = "tasa de empleo v (%)", ylabel = "participación salarial u (%)",
              legend = :outerbottom, legend_columns = 2,
              title = "Plano de fase (v, u): órbita cerrada")
    vline!(p2, [100g.v_est]; color = GRIS, ls = :dot, lw = 1, label = "Ceroclinas (v*, u*)")
    hline!(p2, [100g.u_est]; color = GRIS, ls = :dot, lw = 1, label = "")
    plot!(p2, 100 .* g.v, 100 .* g.u; color = AZUL, lw = 2.5,
          label = "Órbita actual (v₀ = $(round(Int, 100v0_gw))%)")
    scatter!(p2, [100g.v_est], [100g.u_est]; color = NARANJA, ms = 6, msw = 0,
             label = "Centro ($(round(Int, 100g.v_est))%, $(round(Int, 100g.u_est))%)")
    scatter!(p2, [100g.v[1]], [100g.u[1]]; color = AZUL, ms = 5, msw = 0, label = "")

    plot(p1, p2; layout = (1, 2), size = (900, 400), bottom_margin = 4Plots.mm)
end

# ╔═╡ 34de0b96-ca60-11f1-8eab-e3b8dd56c358
"""
Simula el dilema del prisionero espacial de Nowak y May en una cuadrícula L × L
con vecindad de Moore (8 vecinos, frontera periódica) y tentación `b` (R = 1, P = S = 0).
Compara con la misma población cuando los 8 rivales y el modelo imitado se eligen al azar.
"""
function simular_dilema_reticula(L, b; pasos = 50, rng = Xoshiro(12))
    tablero = Int.(rand(rng, L, L) .< 0.5)
    pagos = zeros(Float64, L, L)
    nuevo = similar(tablero)
    historia = [mean(tablero)]

    for _ in 1:pasos
        for r in 1:L, c in 1:L
            coop_vec = 0
            for dr in -1:1, dc in -1:1
                (dr == 0 && dc == 0) && continue
                coop_vec += tablero[mod1(r + dr, L), mod1(c + dc, L)]
            end
            pagos[r, c] = (tablero[r, c] == 1 ? 1.0 : b) * coop_vec
        end
        for r in 1:L, c in 1:L
            mejor_pago = pagos[r, c]
            mejor_est = tablero[r, c]
            for dr in -1:1, dc in -1:1
                (dr == 0 && dc == 0) && continue
                nr, nc = mod1(r + dr, L), mod1(c + dc, L)
                if pagos[nr, nc] > mejor_pago
                    mejor_pago = pagos[nr, nc]
                    mejor_est = tablero[nr, nc]
                end
            end
            nuevo[r, c] = mejor_est
        end
        tablero .= nuevo
        push!(historia, mean(tablero))
    end

    # Comparación en mezcla al azar con el mismo tamaño N = L^2 y 8 interacciones por ronda
    N = L * L
    rng_m = Xoshiro(99)
    est_m = Int.(rand(rng_m, N) .< 0.5)
    pagos_m = zeros(Float64, N)
    historia_m = [mean(est_m)]
    for _ in 1:pasos
        for i in 1:N
            coop_v = 0
            for _ in 1:8
                coop_v += est_m[rand(rng_m, 1:N)]
            end
            pagos_m[i] = (est_m[i] == 1 ? 1.0 : b) * coop_v
        end
        nuevo_m = similar(est_m)
        for i in 1:N
            mejor_p = pagos_m[i]
            mejor_s = est_m[i]
            for _ in 1:8
                j = rand(rng_m, 1:N)
                if pagos_m[j] > mejor_p
                    mejor_p = pagos_m[j]
                    mejor_s = est_m[j]
                end
            end
            nuevo_m[i] = mejor_s
        end
        est_m .= nuevo_m
        push!(historia_m, mean(est_m))
    end
    (tablero = tablero, historia = historia,
     tablero_mezcla = reshape(est_m, L, L), historia_mezcla = historia_m)
end

# ╔═╡ 34de0b14-ca60-11f1-892a-83546bf6513d
sim_pd_ret = simular_dilema_reticula(36, b_pd; pasos = 50, rng = Xoshiro(12));

# ╔═╡ 34de0b14-ca60-11f1-936f-f53921f91168
let
    r = sim_pd_ret
    paleta = cgrad([NARANJA, AZUL])
    p1 = heatmap(r.tablero; color = paleta, clims = (0, 1), aspect_ratio = :equal,
                 colorbar = false, axis = false, ticks = false, grid = false,
                 title = "Retícula 2D (b = $b_pd): $(round(Int, 100r.historia[end]))% coopera")

    t = 0:length(r.historia)-1
    p2 = plot(t, 100 .* r.historia; color = AZUL, lw = 2.4,
              label = "Retícula 2D (clústeres locales)",
              xlabel = "generación", ylabel = "cooperadores (%)",
              ylims = (-2, 102), legend = :outerbottom, legend_columns = 2,
              title = "Evolución temporal de la cooperación (b = $b_pd)")
    plot!(p2, t, 100 .* r.historia_mezcla; color = NARANJA, lw = 2.0, ls = :dash,
          label = "Interacción al azar (mezcla global)")
    plot(p1, p2; layout = Plots.grid(1, 2; widths = [0.44, 0.56]), size = (900, 390), bottom_margin = 4Plots.mm)
end

# ╔═╡ 34de0ba2-ca60-11f1-9c22-439a45260574
"Ejemplo fijo de la sección 0: comparación entre interacción al azar y retícula 2D con b = 1.35."
ejemplo_dilema_espacial = let
    s = simular_dilema_reticula(36, 1.35; pasos = 50, rng = Xoshiro(12))
    (mezcla = (tablero = s.tablero_mezcla, coop = s.historia_mezcla[end]),
     reticula = (tablero = s.tablero, coop = s.historia[end]))
end;

# ╔═╡ 34de0ab0-ca60-11f1-af6f-0d45955b6dcf
let
    mezcla, reticula = ejemplo_dilema_espacial
    paleta = cgrad([NARANJA, AZUL])
    p1 = heatmap(mezcla.tablero; color = paleta, clims = (0, 1), aspect_ratio = :equal,
                 colorbar = false, axis = false, ticks = false, grid = false,
                 title = "Interacción al azar: $(round(Int, 100mezcla.coop))% coopera")
    p2 = heatmap(reticula.tablero; color = paleta, clims = (0, 1), aspect_ratio = :equal,
                 colorbar = false, axis = false, ticks = false, grid = false,
                 title = "Retícula 2D (Nowak y May): $(round(Int, 100reticula.coop))% coopera")
    plot(p1, p2; layout = (1, 2), size = (820, 390))
end

# ╔═╡ 34de0ba2-ca60-11f1-960a-513ef7bcd020
"""
Simula el dilema del prisionero en un grafo arbitrario `g` con tentación `b`
bajo la regla de imitación proporcional de Santos y Pacheco (2005). Se detiene
antes de `pasos` generaciones si toda la red coopera o toda traiciona.
"""
function dilema_en_grafo(g, b, rng; pasos = 5000)
    N = nv(g)
    s = Int.(rand(rng, N) .< 0.5)
    adj = [neighbors(g, i) for i in 1:N]
    grados = [length(adj[i]) for i in 1:N]
    pagos = zeros(Float64, N)
    nuevo = similar(s)
    for _ in 1:pasos
        0 < sum(s) < N || break
        for i in 1:N
            c_vec = 0
            for j in adj[i]
                c_vec += s[j]
            end
            pagos[i] = (s[i] == 1 ? 1.0 : b) * c_vec
        end
        for i in 1:N
            vec = adj[i]
            if isempty(vec)
                nuevo[i] = s[i]
                continue
            end
            j = vec[rand(rng, 1:length(vec))]
            if pagos[j] > pagos[i]
                prob = (pagos[j] - pagos[i]) / (b * max(grados[i], grados[j]))
                nuevo[i] = rand(rng) < prob ? s[j] : s[i]
            else
                nuevo[i] = s[i]
            end
        end
        s .= nuevo
    end
    mean(s)
end

# ╔═╡ 34de0baa-ca60-11f1-98bd-b1db3a2a9ee9
"Compara la fracción final de cooperadores en función de b entre mezcla al azar, retícula 2D y Barabási-Albert."
function comparar_topologias_pd(N = 400; bs = 1.05:0.10:1.95, M = 14)
    L = round(Int, sqrt(N))
    m_ret, se_ret = Float64[], Float64[]
    m_ba, se_ba = Float64[], Float64[]
    m_mez = Float64[]

    for (ib, b) in enumerate(bs)
        vals_ret = [simular_dilema_reticula(L, b; pasos = 45, rng = Xoshiro(200 + 30ib + m)).historia[end] for m in 1:M]
        push!(m_ret, mean(vals_ret))
        push!(se_ret, std(vals_ret) / sqrt(M))

        vals_ba = Float64[]
        for m in 1:M
            g_ba = barabasi_albert(N, 4; seed = 500 + m)
            push!(vals_ba, dilema_en_grafo(g_ba, b, Xoshiro(900 + 30ib + m); pasos = 5000))
        end
        push!(m_ba, mean(vals_ba))
        push!(se_ba, std(vals_ba) / sqrt(M))

        push!(m_mez, simular_dilema_reticula(L, b; pasos = 35, rng = Xoshiro(50 + ib)).historia_mezcla[end])
    end
    (bs = collect(bs), reticula = m_ret, se_ret = se_ret, ba = m_ba, se_ba = se_ba, mezcla = m_mez)
end

# ╔═╡ 34de0b28-ca60-11f1-b8f8-29a3d079547a
barrido_pd = comparar_topologias_pd(400; bs = 1.05:0.10:1.95, M = 14);

# ╔═╡ 34de0b28-ca60-11f1-947b-e7ea80d0681b
let
    b = barrido_pd
    p = plot(b.bs, 100 .* b.reticula; yerror = 100 .* b.se_ret, color = AZUL, lw = 2.4,
             marker = :circle, ms = 4.5, msc = AZUL, label = "Retícula 2D (clústeres espaciales)",
             xlabel = "tentación de traicionar b", ylabel = "cooperadores finales (%)",
             ylims = (-3, 103), legend = :outerbottom, legend_columns = 3, size = (780, 380),
             title = "Cooperación frente a la tentación b: retícula 2D vs. Barabási-Albert")
    plot!(p, b.bs, 100 .* b.ba; yerror = 100 .* b.se_ba, color = NARANJA, lw = 2.4,
          marker = :square, ms = 4.5, msc = NARANJA, label = "Barabási-Albert (concentradores)")
    plot!(p, b.bs, 100 .* b.mezcla; color = GRIS, lw = 1.8, ls = :dash,
          marker = :utriangle, ms = 4, msw = 0, label = "Mezcla al azar")
    p
end

# ╔═╡ 34de0baa-ca60-11f1-83f8-47e4af3c38a9
"""
Simula el proceso muerte-nacimiento (Death-Birth, Ohtsuki et al., 2006) en un
grafo regular `g` de grado `k` con razón beneficio/costo `r = b/c` y selección débil `w`,
partiendo de un único mutante cooperador. Devuelve N · ρ_C.
"""
function fijacion_db(g, r, w; ensayos = 4_000, rng = Xoshiro(1))
    N = nv(g)
    adj = [neighbors(g, i) for i in 1:N]
    k = length(adj[1])
    s = zeros(Int, N)
    exitos = 0
    for e in 1:ensayos
        fill!(s, 0)
        s[mod1(e, N)] = 1
        n_coop = 1
        pasos = 0
        while 0 < n_coop < N && pasos < 40 * N * N
            pasos += 1
            i = rand(rng, 1:N)
            vec_i = adj[i]
            c_en_vec = 0
            for j in vec_i
                c_en_vec += s[j]
            end
            (c_en_vec == 0 && s[i] == 0) && continue
            (c_en_vec == k && s[i] == 1) && continue

            suma_f = 0.0
            suma_fc = 0.0
            for j in vec_i
                cj = 0
                for m in adj[j]
                    cj += s[m]
                end
                fj = 1.0 + w * (r * (cj / k) - s[j])
                suma_f += fj
                s[j] == 1 && (suma_fc += fj)
            end
            nuevo = rand(rng) * suma_f < suma_fc ? 1 : 0
            if nuevo != s[i]
                n_coop += nuevo - s[i]
                s[i] = nuevo
            end
        end
        exitos += (n_coop == N)
    end
    N * (exitos / ensayos)
end

# ╔═╡ 34de0bb4-ca60-11f1-94ec-574c04a04c0f
"N · ρ_C simulado (con su error estándar) bajo actualización muerte-nacimiento en un anillo (k = 2) y una red regular aleatoria (k = 4)."
curva_reciprocidad_red = let
    N = 36
    w = 0.03
    ensayos = 30_000
    bc = collect(1.0:0.5:7.0)
    anillo = watts_strogatz(N, 2, 0.0; seed = 1)
    reg4 = random_regular_graph(N, 4; seed = 4)
    sim_k2 = [fijacion_db(anillo, r, w; ensayos, rng = Xoshiro(100 + round(Int, 10r))) for r in bc]
    sim_k4 = [fijacion_db(reg4, r, w; ensayos, rng = Xoshiro(200 + round(Int, 10r))) for r in bc]
    # Error estándar de N · ρ̂ con ρ̂ = éxitos / ensayos (binomial)
    se(x) = N * sqrt((x / N) * (1 - x / N) / ensayos)
    umbral(k) = (N - 2) / (N / k - 2)
    (bc = bc, sim_k2 = sim_k2, sim_k4 = sim_k4, se_k2 = se.(sim_k2), se_k4 = se.(sim_k4),
     umbral_k2 = umbral(2), umbral_k4 = umbral(4))
end;

# ╔═╡ 34de0b3c-ca60-11f1-aa8a-7d09c5440f22
let
    d = curva_reciprocidad_red
    p = plot(d.bc, d.sim_k2; yerror = d.se_k2, color = AZUL, lw = 2.0, msc = AZUL,
             marker = :circle, ms = 4.5, label = "Anillo (k = 2)",
             xlabel = "razón beneficio / costo (b / c)",
             ylabel = "probabilidad relativa de fijación N · ρ_C",
             legend = :outerbottom, legend_columns = 2, size = (780, 380),
             title = "Regla b/c > k de Ohtsuki et al.: fijación de un mutante cooperador (N = 36)")
    plot!(p, d.bc, d.sim_k4; yerror = d.se_k4, color = NARANJA, lw = 2.0, msc = NARANJA,
          marker = :square, ms = 4.5, label = "Red regular aleatoria (k = 4)")
    hline!(p, [1.0]; color = GRIS, ls = :dot, lw = 1.5, label = "Deriva neutra (N · ρ_C = 1)")
    vline!(p, [2.0, 4.0]; color = GRIS, alpha = 0.5, ls = :dash, lw = 1.2,
           label = "Regla para N grande (b/c = k)")
    vline!(p, [d.umbral_k2]; color = AZUL, alpha = 0.5, ls = :dot, lw = 1.5,
           label = "Umbral finito (N − 2)/(N/k − 2), k = 2")
    vline!(p, [d.umbral_k4]; color = NARANJA, alpha = 0.5, ls = :dot, lw = 1.5,
           label = "Umbral finito, k = 4")
    p
end

# ╔═╡ 34de0bbe-ca60-11f1-bb01-91106b7a7d94
"Muestra si la respuesta de un ejercicio es correcta."
function verificar(f, esperado; rtol = 1e-8, pista = md"")
    caja(tipo, titulo, texto) = Markdown.MD(Markdown.Admonition(tipo, titulo, [texto]))
    resultado = try
        f()
    catch e
        return caja("danger", "Tu código produce un error", md"$(sprint(showerror, e))")
    end
    if resultado === missing || (resultado isa AbstractArray && any(ismissing, resultado))
        caja("warning", "Pendiente", md"Reemplaza `missing` por tu código.")
    elseif resultado isa AbstractArray && esperado isa AbstractArray && size(resultado) == size(esperado) &&
           all(isapprox.(resultado, esperado; rtol, atol = 1e-12))
        caja("correct", "¡Correcto!", md"Bien hecho.")
    elseif resultado isa typeof(esperado) && all(isapprox.(resultado, esperado; rtol, atol = 1e-12))
        caja("correct", "¡Correcto!", md"Bien hecho.")
    else
        caja("danger", "Todavía no", pista)
    end
end

# ╔═╡ 34de0b5a-ca60-11f1-b029-3fe3db9d28c6
verificar(() -> [paso_replicador(0.5, 2.0, 1.0, 0.1), paso_replicador(0.2, 1.0, 3.0, 0.05), paso_replicador(1.0, 4.0, 1.0, 0.1)],
          [0.5 + 0.25 * 1.0 * 0.1, 0.2 + 0.16 * (-2.0) * 0.05, 1.0];
          pista = md"Usa `x + x * (1 - x) * (π1 - π2) * Δt`.")

# ╔═╡ 34de0b64-ca60-11f1-9f3e-373562f15ad8
verificar(() -> (r = equilibrio_goodwin(0.24, 0.40, 0.36, 0.40); r === missing ? missing : collect(r)),
          [0.9, 0.6];
          pista = md"De ``\alpha - \beta u^\ast = 0`` y ``-\gamma + \delta v^\ast = 0`` se obtiene `(γ / δ, α / β)`.")

# ╔═╡ Cell order:
# ╟─34de0a74-ca60-11f1-9cf0-6136d3253a85
# ╟─34de0aa6-ca60-11f1-b2bd-213ba88088be
# ╟─34de0ab0-ca60-11f1-a9d5-d16579bd03b2
# ╟─34de0ab0-ca60-11f1-af6f-0d45955b6dcf
# ╟─34de0aba-ca60-11f1-80cf-7f5030690333
# ╟─34de0ac4-ca60-11f1-b4d1-3b65cc417c97
# ╠═34de0acc-ca60-11f1-aebb-a3da9de9ac44
# ╠═34de0acc-ca60-11f1-a6f2-25a850deb231
# ╟─34de0ad8-ca60-11f1-bb4c-3bca4d55456b
# ╠═34de0ae2-ca60-11f1-b8fb-496d030b090e
# ╟─34de0ae2-ca60-11f1-ba00-55d677ebabd8
# ╟─34de0aec-ca60-11f1-b206-c7b11fa25adf
# ╟─34de0af6-ca60-11f1-9944-2744d9670242
# ╠═34de0b00-ca60-11f1-a9fb-498261e3d40e
# ╟─34de0b00-ca60-11f1-9586-271aca99373c
# ╟─34de0b0a-ca60-11f1-bca5-33dffc73fc7b
# ╟─34de0b0a-ca60-11f1-9afb-53e07cdd8b41
# ╠═34de0b14-ca60-11f1-892a-83546bf6513d
# ╟─34de0b14-ca60-11f1-936f-f53921f91168
# ╟─34de0b1e-ca60-11f1-a4f6-992b2ffc5b2c
# ╠═34de0b28-ca60-11f1-b8f8-29a3d079547a
# ╟─34de0b28-ca60-11f1-947b-e7ea80d0681b
# ╟─34de0b32-ca60-11f1-b8d1-0d9743fd2efa
# ╟─34de0b32-ca60-11f1-829f-cbff897663c8
# ╟─34de0b3c-ca60-11f1-aa8a-7d09c5440f22
# ╟─34de0b46-ca60-11f1-9959-85a5afc38314
# ╟─34de0b50-ca60-11f1-a089-3fbb5729b76a
# ╠═34de0b50-ca60-11f1-8145-d54bf95b14cc
# ╟─34de0b5a-ca60-11f1-b029-3fe3db9d28c6
# ╟─34de0b5a-ca60-11f1-8997-e1325165499b
# ╠═34de0b5a-ca60-11f1-80b1-d1e5f47953a2
# ╟─34de0b64-ca60-11f1-9f3e-373562f15ad8
# ╟─34de0b64-ca60-11f1-a56d-f59833580b58
# ╟─34de0b70-ca60-11f1-9fcd-ef8d2e035d94
# ╠═34de0b78-ca60-11f1-8542-bf8eb1a504d5
# ╠═34de0b78-ca60-11f1-afda-e9cd8a803efe
# ╠═34de0b82-ca60-11f1-a60f-79fdc16f43b4
# ╠═34de0b8c-ca60-11f1-9669-f945d2d90a72
# ╠═34de0b8c-ca60-11f1-bdd9-cb49ddb09e55
# ╠═34de0b96-ca60-11f1-8eab-e3b8dd56c358
# ╠═34de0ba2-ca60-11f1-9c22-439a45260574
# ╠═34de0ba2-ca60-11f1-960a-513ef7bcd020
# ╠═34de0baa-ca60-11f1-98bd-b1db3a2a9ee9
# ╠═34de0baa-ca60-11f1-83f8-47e4af3c38a9
# ╠═34de0bb4-ca60-11f1-94ec-574c04a04c0f
# ╠═34de0bbe-ca60-11f1-bb01-91106b7a7d94
