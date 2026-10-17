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

# ╔═╡ dcbc834c-ca7e-11f1-8de0-131a67af677b
begin
    # Usa el entorno del curso (Code/Project.toml) en lugar del gestor de
    # paquetes propio de Pluto: un solo entorno para todo el repositorio.
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, Statistics, Random, Graphs, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

# ╔═╡ dcbc8252-ca7e-11f1-b5e3-3f6355d9213c
md"""
# Sesión 18: Cooperación en redes

**Objetivos de la sesión**
1. Entender por qué la cooperación desaparece en el **dilema del prisionero**
   cuando la población está bien mezclada.
2. Simular el dilema del prisionero en una **retícula 2D** (Nowak y May) y
   observar cómo sobreviven clústeres de cooperadores por **reciprocidad
   espacial**.
3. Comparar la retícula con una red libre de escala de **Barabási-Albert** y
   explicar el papel de los concentradores (*hubs*) en el sostén de la
   cooperación.
4. **[Avanzado]** Verificar por simulación la regla ``b/c > k`` de
   **reciprocidad de red** (Ohtsuki et al.).

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s18_cooperacion_redes.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): L.B4 (matriz de adyacencia y
grado) y C.B4 (función logística).
"""

# ╔═╡ dcbc82b6-ca7e-11f1-9a76-dd876da9a38b
TableOfContents(title = "Contenido")

# ╔═╡ dcbc82b6-ca7e-11f1-8791-9bbd6d2ef91d
md"""
## 0. El fenómeno: cooperación entre egoístas cuando traicionar paga más

En el **dilema del prisionero**, traicionar siempre deja un pago individual
mayor que cooperar, sin importar qué haga el otro. Si todos interactúan al azar
con todos (como en la ecuación del replicador de la sesión 17), los traidores
explotan a los cooperadores hasta extinguirlos. Pero si los agentes interactúan
únicamente con sus vecinos en el espacio o en una red (Axelrod, 1984; Nowak y
May, 1992), los cooperadores forman **agrupaciones compactas** (*clústeres*)
que se apoyan entre sí y resisten la invasión.
"""

# ╔═╡ dcbc82ca-ca7e-11f1-8327-1ddaf57efead
md"""
Las dos poblaciones empezaron con el mismo 50% de cooperadores (azul) y 50% de
traidores (naranja), y con la misma tentación a traicionar (``b = 1.35``). Al
emparejarse al azar, la cooperación desaparece; en la retícula 2D, sobrevive en
islas compactas de cooperadores.
"""

# ╔═╡ dcbc82d4-ca7e-11f1-aade-af686ad6820d
md"""
## 1. El dilema del prisionero: mezcla global frente a red local

En el **dilema del prisionero débil** de Nowak y May (1992), cada encuentro
entre dos agentes da los siguientes pagos:
- Si ambos **cooperan** (``C``), cada uno recibe la recompensa ``R = 1``.
- Si un **traidor** (``D``) se enfrenta a un cooperador, el traidor recibe la
  **tentación** ``T = b > 1`` y el cooperador recibe ``S = 0``.
- Si ambos **traicionan**, cada uno recibe ``P = 0``.

| Jugador fila \ columna | Cooperar (``C``) | Traicionar (``D``) |
|---|---|---|
| **Cooperar (``C``)** | ``1`` | ``0`` |
| **Traicionar (``D``)** | ``b > 1`` | ``0`` |

**¿Por qué desaparece la cooperación en una población bien mezclada?** Si todos
se emparejan al azar con todos y una fracción ``x`` coopera, tanto cooperadores
como traidores encuentran la misma fracción ``x`` de cooperadores. Sus pagos
esperados son ``\pi_C(x) = x`` y ``\pi_D(x) = b x``. Como ``b > 1``, en la
ecuación del replicador de la sesión 17:

```math
\dot{x} = x(1 - x)\bigl[\pi_C(x) - \pi_D(x)\bigr] = (1 - b)x^2(1 - x) < 0 ,
```

así que ``x(t) \to 0`` siempre.

**¿Qué cambia en una red?** Cuando las interacciones ocurren sobre un grafo con
matriz de adyacencia ``A_{ij}`` (repaso L.B4):
1. En cada ronda, cada nodo ``i`` juega con todos sus vecinos en la red y suma
   sus pagos ``\Pi_i``. Si tiene ``k_i^C`` vecinos cooperadores, un cooperador
   gana ``\Pi_i = k_i^C`` y un traidor gana ``\Pi_i = b\,k_i^C``.
2. Cada nodo imita a vecinos con mayor pago:
   - En la **retícula 2D** usamos la regla determinista de Nowak y May (1992):
     adoptar la estrategia del mayor pago en el vecindario (incluido el propio).
   - En la **red de Barabási-Albert** usamos la regla de imitación proporcional
     de Santos y Pacheco (2005): elegir un vecino ``j`` al azar y, si
     ``\Pi_j > \Pi_i``, copiar su estrategia con probabilidad
     ``(\Pi_j - \Pi_i) / (b \max(k_i, k_j))`` (otra alternativa muy usada en
     física estadística es copiar con la función logística
     ``\sigma(\beta(\Pi_j - \Pi_i))``, repaso C.B4).
"""

# ╔═╡ dcbc82de-ca7e-11f1-af4c-4794b9a96e8a
md"""
## 2. Laboratorio: clústeres de cooperadores en la retícula 2D

Mueve la tentación a traicionar ``b > 1`` y observa cómo evoluciona una
retícula 2D de ``36 \times 36`` agentes (8 vecinos de Moore) frente a la misma
población emparejada al azar:

Tentación de traicionar ``b``: $(@bind b_pd Slider([1.05, 1.15, 1.25, 1.35, 1.55, 1.75, 1.95], default = 1.35, show_value = true))
"""

# ╔═╡ dcbc82f2-ca7e-11f1-96ce-8ba617e259b2
md"""
Al inicio, los traidores invaden a los cooperadores aislados y la cooperación
cae en las primeras generaciones. Pero en la retícula 2D, los cooperadores que
quedan juntos forman **bloques compactos**: un cooperador en el interior o en la
frontera plana de un bloque recibe el apoyo de muchos vecinos cooperadores,
mientras que un traidor en la orilla sólo puede explotar a los pocos
cooperadores que toca.
"""

# ╔═╡ dcbc82fc-ca7e-11f1-837c-f180a15b3d7f
md"""
## 3. Laboratorio: retícula 2D frente a concentradores en Barabási-Albert

¿Cómo depende la cooperación de la estructura de la red en todo el rango de
``b``? Comparamos tres topologías con ``N = 400`` agentes y el mismo grado
medio (``\langle k \rangle = 8``). La dinámica en la red de Barabási-Albert es
lenta, así que la dejamos correr hasta 5,000 generaciones (o hasta que toda la
red coopere o traicione):
1. **Mezcla al azar** (en cada ronda se sortean 8 rivales al azar);
2. **Retícula 2D** (8 vecinos de Moore fijos);
3. **Red de Barabási-Albert** (libre de escala, ``m = 4 \Rightarrow \langle k \rangle \approx 8``, Santos y Pacheco, 2005).
"""

# ╔═╡ dcbc8306-ca7e-11f1-9857-6dc33d0d7bbd
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

# ╔═╡ dcbc8310-ca7e-11f1-a37d-21c0070ff505
md"""
## 4. [Avanzado] La condición ``b/c > k`` de reciprocidad de red

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

La gráfica muestra ``N \cdot \rho_C`` simulado (30,000 ensayos por punto, con
``N = 36`` y ``w = 0.03``) en un anillo (``k = 2``, umbral finito
``34/16 \approx 2.1``) y en una red regular aleatoria de grado ``k = 4``
(umbral finito ``34/7 \approx 4.9``). Cada curva cruza la deriva neutra
``N \rho_C = 1`` cerca de su umbral finito, un poco a la derecha de
``b/c = k``.
"""

# ╔═╡ dcbc831a-ca7e-11f1-9eed-c593dd7daf91
md"""
## Resumen

| Topología / Mecanismo | Regla microscópica | Comportamiento colectivo |
|---|---|---|
| **Mezcla global** | Emparejamiento al azar (``R=1, S=0, T=b>1, P=0``) | ``\pi_D = bx > x = \pi_C``; los cooperadores se extinguen siempre |
| **Retícula 2D (Nowak y May)** | Imitar al vecino con mayor pago en vecindad de Moore | **Reciprocidad espacial**: sobreviven bloques compactos de cooperadores con escalones en ``b = m/n`` |
| **Barabási-Albert (Santos y Pacheco)** | Imitación proporcional en red libre de escala | Los **concentradores (*hubs*)** cooperadores sostienen la cooperación hasta valores altos de ``b`` |
| **Regla ``b/c > k`` [Avanzado]** | Selección débil (muerte-nacimiento) en grafo de grado ``k`` | Un mutante cooperador supera la deriva neutra si ``b/c > k`` (``N \gg k``) |
"""

# ╔═╡ dcbc8324-ca7e-11f1-9e0f-7985bda764fa
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `pago_nowak_may(estrategia, vecinos_coop, b)`, que devuelva el
pago total de un agente con `estrategia` (`1` si coopera, `0` si traiciona) que
tiene `vecinos_coop` vecinos cooperadores en el dilema del prisionero débil con
tentación `b`.
"""

# ╔═╡ dcbc832e-ca7e-11f1-b6e3-c7b3068e2638
function pago_nowak_may(estrategia, vecinos_coop, b)
    missing  # reemplaza esta línea
end

# ╔═╡ dcbc833a-ca7e-11f1-a128-6fadcb1cb432
md"""
**2.** Escribe `coop_favorecida(b, c, k)`, que devuelva `true` si la razón
beneficio-costo `b / c` supera el grado `k` de la red (regla ``b/c > k`` de
Ohtsuki et al. para redes grandes), y `false` en caso contrario.
"""

# ╔═╡ dcbc833a-ca7e-11f1-b8e8-053ce457fa12
function coop_favorecida(b, c, k)
    missing  # reemplaza esta línea
end

# ╔═╡ dcbc8342-ca7e-11f1-bc64-d1dbd83f71f5
md"""
### Ruta básica: interpretación

3. ¿Por qué un cooperador aislado rodeado de traidores en la retícula 2D
   desaparece en una sola generación, mientras que un bloque cuadrado de
   ``3 \times 3`` cooperadores logra sobrevivir cuando ``b = 1.35``?
4. Explica por qué un concentrador (*hub*) traidor en una red de
   Barabási-Albert termina perdiendo su alto pago inicial cuando sus vecinos lo
   imitan, mientras que un *hub* cooperador refuerza su propio pago.
5. Según la condición ``b/c > k`` de Ohtsuki et al. (2006), ¿en cuál red es más
   fácil que prospere un mutante cooperador bajo selección débil: en un anillo
   (``k = 2``), en una cuadrícula (``k = 4``) o en una red aleatoria densa con
   ``k = 20``? Explica intuitivamente por qué tener demasiados vecinos dificulta
   la reciprocidad local.

### Ruta avanzada

6. En una cuadrícula con vecindad de Moore (8 vecinos), calcula para qué valor
   crítico de ``b`` un traidor en la orilla plana de un bloque de ``3 \times 3``
   cooperadores (que toca a 3 cooperadores) obtiene un pago mayor que un
   cooperador del borde del bloque (que toca a 5 cooperadores). Compáralo con
   la caída de la curva azul en la sección 3.
7. En la sección 3, modifica la regla de actualización de la retícula para usar
   la regla logística de Fermi ``\sigma(\beta(\Pi_j - \Pi_i))`` con un vecino al
   azar. ¿Cómo cambia la fracción de cooperadores cuando sube el ruido
   ``T = 1/\beta``?

## Para la próxima sesión

- Lectura: Nowak y May (1992). Complementaria: Slanina (2014), sección 8.1.
  Para la ruta avanzada, Santos y Pacheco (2005); Ohtsuki et al. (2006); de
  consulta, Schecter y Gintis (2016).
- En la **sesión 19** (`Code/s19_juego_minoria.jl`) estudiaremos qué pasa cuando
  ganar depende de estar en la **minoría**: el problema del bar El Farol y el
  **juego de la minoría** como modelo de eficiencia y volatilidad en un mercado
  financiero.
"""

# ╔═╡ dcbc834c-ca7e-11f1-9161-1b6e57d8c4b9
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ dcbc8356-ca7e-11f1-b959-45afa541efbc
begin
    const AZUL       = "#2a78d6"   # serie principal, cooperadores
    const NARANJA    = "#eb6834"   # segunda serie o teoría, traidores
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const BLANCO     = "#ffffff"
end;

# ╔═╡ dcbc8356-ca7e-11f1-90d9-6549ba4857a2
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

# ╔═╡ dcbc82e8-ca7e-11f1-a2cb-3b4b694fd5fc
sim_pd_ret = simular_dilema_reticula(36, b_pd; pasos = 50, rng = Xoshiro(12));

# ╔═╡ dcbc82e8-ca7e-11f1-a26f-7f79b41d1e44
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

# ╔═╡ dcbc8360-ca7e-11f1-9c9b-f1602074c756
"Ejemplo fijo de la sección 0: comparación entre interacción al azar y retícula 2D con b = 1.35."
ejemplo_dilema_espacial = let
    s = simular_dilema_reticula(36, 1.35; pasos = 50, rng = Xoshiro(12))
    (mezcla = (tablero = s.tablero_mezcla, coop = s.historia_mezcla[end]),
     reticula = (tablero = s.tablero, coop = s.historia[end]))
end;

# ╔═╡ dcbc82c0-ca7e-11f1-a4b1-07be0e04cfe4
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

# ╔═╡ dcbc836c-ca7e-11f1-9245-334a6d7b2c4a
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

# ╔═╡ dcbc836c-ca7e-11f1-b6f4-39b822db9fee
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

# ╔═╡ dcbc82fc-ca7e-11f1-a858-771cdd88e157
barrido_pd = comparar_topologias_pd(400; bs = 1.05:0.10:1.95, M = 14);

# ╔═╡ dcbc8306-ca7e-11f1-a6a0-41f86f268f56
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

# ╔═╡ dcbc8374-ca7e-11f1-957c-251d9235841b
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

# ╔═╡ dcbc837e-ca7e-11f1-b6d1-f564bbdcd35a
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

# ╔═╡ dcbc831a-ca7e-11f1-8900-713465af8fa1
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

# ╔═╡ dcbc8388-ca7e-11f1-83dc-43881757d79a
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

# ╔═╡ dcbc832e-ca7e-11f1-9a28-3ff0a0494e9f
verificar(() -> [pago_nowak_may(1, 5, 1.4), pago_nowak_may(0, 3, 1.4), pago_nowak_may(0, 0, 1.8)],
          [5.0, 3 * 1.4, 0.0];
          pista = md"Si `estrategia == 1` recibe `1.0 * vecinos_coop`; si es `0` recibe `b * vecinos_coop`.")

# ╔═╡ dcbc833a-ca7e-11f1-8adf-adebe1f9289a
verificar(() -> [coop_favorecida(5.0, 1.0, 4), coop_favorecida(3.0, 1.0, 4), coop_favorecida(6.0, 2.0, 2)],
          [true, false, true];
          pista = md"Compara `b / c > k`.")

# ╔═╡ Cell order:
# ╟─dcbc8252-ca7e-11f1-b5e3-3f6355d9213c
# ╟─dcbc82b6-ca7e-11f1-9a76-dd876da9a38b
# ╟─dcbc82b6-ca7e-11f1-8791-9bbd6d2ef91d
# ╟─dcbc82c0-ca7e-11f1-a4b1-07be0e04cfe4
# ╟─dcbc82ca-ca7e-11f1-8327-1ddaf57efead
# ╟─dcbc82d4-ca7e-11f1-aade-af686ad6820d
# ╟─dcbc82de-ca7e-11f1-af4c-4794b9a96e8a
# ╠═dcbc82e8-ca7e-11f1-a2cb-3b4b694fd5fc
# ╟─dcbc82e8-ca7e-11f1-a26f-7f79b41d1e44
# ╟─dcbc82f2-ca7e-11f1-96ce-8ba617e259b2
# ╟─dcbc82fc-ca7e-11f1-837c-f180a15b3d7f
# ╠═dcbc82fc-ca7e-11f1-a858-771cdd88e157
# ╟─dcbc8306-ca7e-11f1-a6a0-41f86f268f56
# ╟─dcbc8306-ca7e-11f1-9857-6dc33d0d7bbd
# ╟─dcbc8310-ca7e-11f1-a37d-21c0070ff505
# ╟─dcbc831a-ca7e-11f1-8900-713465af8fa1
# ╟─dcbc831a-ca7e-11f1-9eed-c593dd7daf91
# ╟─dcbc8324-ca7e-11f1-9e0f-7985bda764fa
# ╠═dcbc832e-ca7e-11f1-b6e3-c7b3068e2638
# ╟─dcbc832e-ca7e-11f1-9a28-3ff0a0494e9f
# ╟─dcbc833a-ca7e-11f1-a128-6fadcb1cb432
# ╠═dcbc833a-ca7e-11f1-b8e8-053ce457fa12
# ╟─dcbc833a-ca7e-11f1-8adf-adebe1f9289a
# ╟─dcbc8342-ca7e-11f1-bc64-d1dbd83f71f5
# ╟─dcbc834c-ca7e-11f1-9161-1b6e57d8c4b9
# ╠═dcbc834c-ca7e-11f1-8de0-131a67af677b
# ╠═dcbc8356-ca7e-11f1-b959-45afa541efbc
# ╠═dcbc8356-ca7e-11f1-90d9-6549ba4857a2
# ╠═dcbc8360-ca7e-11f1-9c9b-f1602074c756
# ╠═dcbc836c-ca7e-11f1-9245-334a6d7b2c4a
# ╠═dcbc836c-ca7e-11f1-b6f4-39b822db9fee
# ╠═dcbc8374-ca7e-11f1-957c-251d9235841b
# ╠═dcbc837e-ca7e-11f1-b6d1-f564bbdcd35a
# ╠═dcbc8388-ca7e-11f1-83dc-43881757d79a
