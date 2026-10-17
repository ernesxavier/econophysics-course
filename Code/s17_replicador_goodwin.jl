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

# ╔═╡ d6098612-ca7e-11f1-8eb0-8107b7571c11
begin
    # Usa el entorno del curso (Code/Project.toml) en lugar del gestor de
    # paquetes propio de Pluto: un solo entorno para todo el repositorio.
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, Statistics, Random, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

# ╔═╡ d6098540-ca7e-11f1-9151-17c3274092ef
md"""
# Sesión 17: Replicador y ciclo de Goodwin

**Objetivos de la sesión**
1. Usar una **ecuación diferencial ordinaria (EDO)** como regla de
   actualización paso a paso (método de Euler) e interpretar su retrato de fase.
2. Simular la **dinámica del replicador** a partir de la imitación de pagos y
   encontrar el equilibrio mixto del juego **Halcón-Paloma**.
3. Conectar el sistema depredador-presa de **Lotka-Volterra** con el **ciclo
   distributivo de Goodwin** entre tasa de empleo y participación salarial.
4. **[Avanzado]** Analizar la estabilidad de puntos fijos con el **jacobiano**
   y verificar la **cantidad conservada** de Goodwin.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s17_replicador_goodwin.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): C.B3 (la derivada como tasa
de cambio), C.I1 (desarrollo de Taylor), C.I3 (funciones de varias variables),
C.I4 (EDO en una dimensión y retrato de fase), L.I1 (eigenvalores). Ruta
avanzada: C.A2 (sistemas de EDO en dos dimensiones), L.A3 (jacobiano y
estabilidad lineal).
"""

# ╔═╡ d609857c-ca7e-11f1-82c4-69d50d94bd5d
TableOfContents(title = "Contenido")

# ╔═╡ d6098586-ca7e-11f1-9f1b-d13a44fda248
md"""
## 0. El fenómeno: ¿por qué no todos son agresivos? ¿Y por qué la economía oscila sin choques?

En la sesión 15 los agentes copiaban la opinión de un vecino sin preguntar qué
tan bien le iba. Hoy damos el siguiente paso: cada agente elige una
**estrategia**, recibe un **pago** según lo que hacen los demás y tiende a
imitar a quienes obtienen mejores resultados.

Esta idea explica dos fenómenos clásicos:

1. **Coexistencia de agresivos y pacíficos** (Maynard Smith y Price, 1973). Si
   ser agresivo (*halcón*) permite arrebatar recursos a los pacíficos
   (*palomas*), ¿por qué no todos se vuelven agresivos? Porque cuando casi todos
   son halcones, pelean entre sí y el costo del conflicto supera al premio.
2. **Ciclos económicos sin choques externos** (Goodwin, 1967). Cuando el empleo
   es alto, los salarios suben y comprimen las ganancias; con menos ganancias,
   cae la inversión y sube el desempleo; el desempleo modera los salarios,
   recupera las ganancias y relanza el empleo. La economía oscila perpetuamente
   con la misma matemática que las poblaciones de presas y depredadores de
   Lotka y Volterra.
"""

# ╔═╡ d6098590-ca7e-11f1-9674-bdf05e7b2f82
md"""
En la izquierda, sin importar si la población empieza con 10% o con 92% de
agresivos, converge al mismo equilibrio mixto. En la derecha, empleo y salarios
se persiguen en un ciclo perpetuo sin que ocurra ningún choque exógeno.
"""

# ╔═╡ d609859a-ca7e-11f1-ab80-e7bc0b9dd447
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

# ╔═╡ d609859a-ca7e-11f1-8750-1f67c3c055c7
# Un paso del método de Euler: nuevo valor = valor actual + tasa de cambio × paso de tiempo
paso_euler(f, x, Δt) = x + f(x) * Δt

# ╔═╡ d60985a6-ca7e-11f1-bddc-63d975b81747
# Si f(x) = x(1 - x)(0.6 - x), el punto fijo x* = 0.6 atrae desde ambos lados
paso_euler(x -> x * (1 - x) * (0.6 - x), 0.2, 0.1),
paso_euler(x -> x * (1 - x) * (0.6 - x), 0.8, 0.1)

# ╔═╡ d60985ae-ca7e-11f1-96b8-39e55239dfa1
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

# ╔═╡ d60985c2-ca7e-11f1-a2cf-c7e6844615e2
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

# ╔═╡ d60985cc-ca7e-11f1-8c05-19ac41debd09
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

# ╔═╡ d60985e0-ca7e-11f1-84eb-d71476a373ab
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

# ╔═╡ d60985e0-ca7e-11f1-8d26-3bc153968f4c
md"""
## 4. [Avanzado] Estabilidad con el jacobiano y la cantidad conservada de Goodwin

Linealizando el sistema de Goodwin alrededor del punto fijo interior
``(v^\ast, u^\ast) = (\gamma_G/\delta_G,\, \alpha_G/\beta_G)`` mediante el
desarrollo de Taylor en dos variables (repaso C.I1 y C.I3), la matriz
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
los eigenvalores (repaso L.I1) son **imaginarios puros** (polinomio
característico ``\lambda^2 + \alpha_G \gamma_G = 0``):

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
"""

# ╔═╡ d60985ea-ca7e-11f1-bc1f-c901504332f0
md"""
## Resumen

| Modelo | Ecuación / Regla microscópica | Comportamiento colectivo |
|---|---|---|
| **Dinámica del replicador (Halcón-Paloma)** | ``\dot{x} = x(1-x)(\pi_H - \pi_P) = \tfrac{1}{2}x(1-x)(V - Cx)`` | Equilibrio interior estable ``x^\ast = V/C``: coexistencia de agresivos y pacíficos |
| **Goodwin – Lotka-Volterra** | ``\dot{v} = v(\alpha_G - \beta_G u)``, ``\dot{u} = u(-\gamma_G + \delta_G v)`` | Órbitas cerradas alrededor de ``(\gamma_G/\delta_G,\, \alpha_G/\beta_G)``; el empleo adelanta a los salarios |
"""

# ╔═╡ d60985f4-ca7e-11f1-a16e-238a90e7c64b
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `paso_replicador(x, π1, π2, Δt)`, que devuelva el nuevo valor de
la fracción `x` después de un paso de tiempo `Δt` con el método de Euler
aplicado a la ecuación del replicador ``\dot{x} = x(1 - x)(\pi_1 - \pi_2)``.
"""

# ╔═╡ d60985f4-ca7e-11f1-a89d-d7a32949ce4b
function paso_replicador(x, π1, π2, Δt)
    missing  # reemplaza esta línea
end

# ╔═╡ d60985fe-ca7e-11f1-a94e-df54157eac9a
md"""
**2.** Escribe `equilibrio_goodwin(α, β, γ, δ)`, que devuelva la tupla
`(v_est, u_est)` con la tasa de empleo ``v^\ast`` y la participación salarial
``u^\ast`` en el punto fijo interior del modelo de Goodwin
``\dot{v} = v(\alpha - \beta u)``, ``\dot{u} = u(-\gamma + \delta v)``.
"""

# ╔═╡ d60985fe-ca7e-11f1-a3ce-a1f2d6ae027c
function equilibrio_goodwin(α, β, γ, δ)
    missing  # reemplaza esta línea
end

# ╔═╡ d609860a-ca7e-11f1-b29b-45854126a07d
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

### Ruta avanzada

6. Verifica que para la matriz jacobiana ``J(v^\ast, u^\ast)`` del modelo de
   Goodwin, ``\mathrm{tr}(J) = 0`` y ``\det(J) = \alpha_G \gamma_G``. Calcula el
   periodo de las oscilaciones pequeñas cuando ``\alpha_G = 0.24`` y
   ``\gamma_G = 0.36``.
7. Deriva respecto al tiempo la función
   ``V(v, u) = \delta_G v - \gamma_G \ln v + \beta_G u - \alpha_G \ln u`` usando
   la regla de la cadena y demuestra que ``\dot{V} = 0`` a lo largo de cualquier
   trayectoria de Goodwin.

## Para la próxima sesión

- Lectura: Requião da Cunha (2021), capítulo 7. Para la ruta avanzada, Goodwin
  (1967); de consulta, Strogatz (2019).
- En la **sesión 18** (`Code/s18_cooperacion_redes.jl`) veremos qué pasa cuando
  traicionar domina estrictamente a cooperar (**dilema del prisionero**), pero
  los agentes interactúan en una **red**: reciprocidad espacial en retículas y
  el papel de los concentradores en redes de Barabási-Albert.
"""

# ╔═╡ d6098612-ca7e-11f1-bdde-a97d8c9b068f
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ d609861c-ca7e-11f1-9c55-83195afc1dc7
begin
    const AZUL       = "#2a78d6"   # serie principal, halcones / empleo
    const NARANJA    = "#eb6834"   # segunda serie o teoría, salarios
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const BLANCO     = "#ffffff"
end;

# ╔═╡ d6098626-ca7e-11f1-a304-bf7f062c07cb
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

# ╔═╡ d60985b8-ca7e-11f1-9af4-1f1fb3b7b55d
sim_hp = simular_halcon_paloma(2.0, C_hp, x0_hp; N = 250, T = 18.0, rng = Xoshiro(42));

# ╔═╡ d60985b8-ca7e-11f1-9446-df311c3ecc53
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

# ╔═╡ d6098626-ca7e-11f1-9ed0-395070c035ec
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

# ╔═╡ d6098630-ca7e-11f1-9ff3-dd2755fbacdb
"Simula la órbita principal y una órbita de referencia del ciclo de Goodwin."
function simular_goodwin(α, β, γ, δ; v0 = 0.85, u0 = α / β, T = 45.0)
    orb = integrar_goodwin(α, β, γ, δ, v0, u0; T = T)
    ref = integrar_goodwin(α, β, γ, δ, 0.88, α / β; T = T)
    (t = orb.t, v = orb.v, u = orb.u, v_ref = ref.v, u_ref = ref.u,
     v_est = γ / δ, u_est = α / β)
end

# ╔═╡ d6098590-ca7e-11f1-be87-d5955be78822
let
    hp = simular_halcon_paloma(2.0, 4.0, 0.1; N = 250, T = 18.0, rng = Xoshiro(42))
    p1 = plot(hp.t_edo, 100 .* hp.x_edo; color = AZUL, lw = 2.4,
              label = "Desde 10% de halcones", xlabel = "tiempo",
              ylabel = "halcones en la población (%)", ylims = (0, 100), legend = :topright,
              title = "Halcón-Paloma: coexistencia en V/C = 50%")
    plot!(p1, hp.t_edo, 100 .* hp.x_alt; color = NARANJA, lw = 2.4, ls = :dash,
          label = "Desde 92% de halcones")
    hline!(p1, [50]; color = GRIS, ls = :dot, label = "")

    gw = simular_goodwin(0.24, 0.40, 0.36, 0.40; v0 = 0.84, u0 = 0.60, T = 45.0)
    p2 = plot(gw.t, 100 .* gw.v; color = AZUL, lw = 2.2, label = "Empleo v(t)",
              xlabel = "tiempo (años)", ylabel = "porcentaje (%)", ylims = (42, 102),
              legend = :bottomright, title = "Ciclo de Goodwin: oscilación sin choques")
    plot!(p2, gw.t, 100 .* gw.u; color = NARANJA, lw = 2.2, ls = :dash,
          label = "Salarios / PIB u(t)")
    plot(p1, p2; layout = (1, 2), size = (880, 350), bottom_margin = 4Plots.mm)
end

# ╔═╡ d60985cc-ca7e-11f1-9991-7131fc208b43
sim_gw = simular_goodwin(0.24, β_gw, 0.36, 0.40; v0 = v0_gw, u0 = 0.24 / β_gw, T = 45.0);

# ╔═╡ d60985d8-ca7e-11f1-bfa7-bdfbf8ffbcd9
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

# ╔═╡ d6098630-ca7e-11f1-9d99-e3cf169dab27
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

# ╔═╡ d60985f4-ca7e-11f1-9210-559d60c0761e
verificar(() -> [paso_replicador(0.5, 2.0, 1.0, 0.1), paso_replicador(0.2, 1.0, 3.0, 0.05), paso_replicador(1.0, 4.0, 1.0, 0.1)],
          [0.5 + 0.25 * 1.0 * 0.1, 0.2 + 0.16 * (-2.0) * 0.05, 1.0];
          pista = md"Usa `x + x * (1 - x) * (π1 - π2) * Δt`.")

# ╔═╡ d609860a-ca7e-11f1-9abe-3dfb2536a760
verificar(() -> (r = equilibrio_goodwin(0.24, 0.40, 0.36, 0.40); r === missing ? missing : collect(r)),
          [0.9, 0.6];
          pista = md"De ``\alpha - \beta u^\ast = 0`` y ``-\gamma + \delta v^\ast = 0`` se obtiene `(γ / δ, α / β)`.")

# ╔═╡ Cell order:
# ╟─d6098540-ca7e-11f1-9151-17c3274092ef
# ╟─d609857c-ca7e-11f1-82c4-69d50d94bd5d
# ╟─d6098586-ca7e-11f1-9f1b-d13a44fda248
# ╟─d6098590-ca7e-11f1-be87-d5955be78822
# ╟─d6098590-ca7e-11f1-9674-bdf05e7b2f82
# ╟─d609859a-ca7e-11f1-ab80-e7bc0b9dd447
# ╠═d609859a-ca7e-11f1-8750-1f67c3c055c7
# ╠═d60985a6-ca7e-11f1-bddc-63d975b81747
# ╟─d60985ae-ca7e-11f1-96b8-39e55239dfa1
# ╠═d60985b8-ca7e-11f1-9af4-1f1fb3b7b55d
# ╟─d60985b8-ca7e-11f1-9446-df311c3ecc53
# ╟─d60985c2-ca7e-11f1-a2cf-c7e6844615e2
# ╟─d60985cc-ca7e-11f1-8c05-19ac41debd09
# ╠═d60985cc-ca7e-11f1-9991-7131fc208b43
# ╟─d60985d8-ca7e-11f1-bfa7-bdfbf8ffbcd9
# ╟─d60985e0-ca7e-11f1-84eb-d71476a373ab
# ╟─d60985e0-ca7e-11f1-8d26-3bc153968f4c
# ╟─d60985ea-ca7e-11f1-bc1f-c901504332f0
# ╟─d60985f4-ca7e-11f1-a16e-238a90e7c64b
# ╠═d60985f4-ca7e-11f1-a89d-d7a32949ce4b
# ╟─d60985f4-ca7e-11f1-9210-559d60c0761e
# ╟─d60985fe-ca7e-11f1-a94e-df54157eac9a
# ╠═d60985fe-ca7e-11f1-a3ce-a1f2d6ae027c
# ╟─d609860a-ca7e-11f1-9abe-3dfb2536a760
# ╟─d609860a-ca7e-11f1-b29b-45854126a07d
# ╟─d6098612-ca7e-11f1-bdde-a97d8c9b068f
# ╠═d6098612-ca7e-11f1-8eb0-8107b7571c11
# ╠═d609861c-ca7e-11f1-9c55-83195afc1dc7
# ╠═d6098626-ca7e-11f1-a304-bf7f062c07cb
# ╠═d6098626-ca7e-11f1-9ed0-395070c035ec
# ╠═d6098630-ca7e-11f1-9ff3-dd2755fbacdb
# ╠═d6098630-ca7e-11f1-9d99-e3cf169dab27
