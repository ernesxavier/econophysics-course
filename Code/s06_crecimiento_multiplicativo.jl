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

# ╔═╡ 7ad21efc-ca89-11f1-a183-4f3484351f7e
begin
    # Usa el entorno del curso (Code/Project.toml) en lugar del gestor de
    # paquetes propio de Pluto: un solo entorno para todo el repositorio.
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, CSV, DataFrames, Statistics, StatsBase, Distributions,
          Random, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

# ╔═╡ 7ad21df8-ca89-11f1-97ba-cdac1f0b5427
md"""
# Sesión 6: Crecimiento multiplicativo

**Objetivos de la sesión**
1. Probar con datos de México la **ley de Gibrat**: ¿crecen igual, en
   proporción, los municipios grandes y los pequeños?
2. Entender por qué el crecimiento multiplicativo produce distribuciones
   **lognormales**.
3. Descubrir un mecanismo genérico que produce **leyes de potencias**: el
   crecimiento multiplicativo con una barrera inferior (proceso de Kesten).
4. Medir en el laboratorio cómo cambia el exponente con los parámetros, y
   entender por qué las ciudades siguen la ley de Zipf.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s06_crecimiento_multiplicativo.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): C.B1 (logaritmos), P.I3
(teorema del límite central), P.I6 (la lognormal), E.I6 (Gibrat y Kesten),
P.I5 (leyes de potencias). Ruta avanzada: E.A4 y E.A5.
"""

# ╔═╡ 7ad21e2a-ca89-11f1-8f14-0fb897c70423
TableOfContents(title = "Contenido")

# ╔═╡ 7ad21e34-ca89-11f1-8dff-b3dc92cb2354
md"""
## 0. El fenómeno: ¿crecen igual los grandes y los pequeños?

En la sesión 3 vimos que el tamaño de las localidades de México sigue la ley
de Zipf. Hoy buscamos un **mecanismo** que la produzca. La hipótesis más
simple es la de Robert Gibrat (1931): cada unidad (una ciudad, una empresa)
crece en una **proporción aleatoria** que no depende de su tamaño. Una ciudad
de un millón y un pueblo de mil habitantes tendrían la misma probabilidad de
crecer 10%.

Podemos ponerla a prueba con los censos de 2010 y 2020 (INEGI). Los datos
están en `Data/municipios_2010_2020.csv` y se regeneran con
`Code/descargar_inegi.jl`.
"""

# ╔═╡ 7ad21e3e-ca89-11f1-80aa-3f29a0e6eecc
municipios = CSV.read(joinpath(@__DIR__, "..", "Data", "municipios_2010_2020.csv"), DataFrame)

# ╔═╡ 7ad21e3e-ca89-11f1-bfe2-41c88266bb20
begin
    p10 = Float64.(municipios.poblacion_2010)
    p20 = Float64.(municipios.poblacion_2020)
    g = log.(p20 ./ p10)          # tasa de crecimiento logarítmica 2010–2020
end;

# ╔═╡ 7ad21e5c-ca89-11f1-b890-57b1e87d312b
gibrat_reg = let x = log.(p10)
    b = cov(x, g) / var(x)
    (b = b, R2 = cor(x, g)^2)
end;

# ╔═╡ 7ad21e54-ca89-11f1-9fcb-e733dbcb5b8f
md"""
La ley de Gibrat dice dos cosas, y los datos responden distinto a cada una:

- **La dispersión del crecimiento no depende del tamaño.** Se cumple bastante
  bien: la desviación estándar por decil se mueve en un rango estrecho, sin
  una tendencia clara.
- **El crecimiento promedio no depende del tamaño.** Se cumple sólo en
  parte: la pendiente de la regresión del crecimiento contra el logaritmo de
  la población es $(round(gibrat_reg.b, digits = 3)) (con ``R^2 = `` $(round(gibrat_reg.R2, digits = 3))).
  Los municipios más pequeños crecieron menos, o perdieron población, que
  los grandes: es la **urbanización**, la migración del campo a las ciudades.

Un cuidado: 13 municipios se crearon entre 2010 y 2020 separándose de otros,
que perdieron población "en el papel". Son pocos y no cambian la imagen.

Gibrat es una buena **primera aproximación**. Veamos qué distribución de
tamaños produce.
"""

# ╔═╡ 7ad21e66-ca89-11f1-bd6f-cbeddd99a19c
md"""
## 1. Repaso: el logaritmo de un producto es una suma

Si una cantidad empieza en ``x_0`` y en cada periodo se multiplica por un
factor aleatorio ``a_t``, después de ``n`` periodos vale

```math
x_n = x_0\, a_1 a_2 \cdots a_n
\qquad\Longrightarrow\qquad
\ln x_n = \ln x_0 + \sum_{t=1}^{n} \ln a_t .
```

El logaritmo convierte el producto en una **suma** de términos independientes.
Por el teorema del límite central, ``\ln x_n`` es aproximadamente normal, así
que ``x_n`` es aproximadamente **lognormal** (repaso P.I6). Mueve ``n``: el
histograma de ``x_n`` se vuelve cada vez más asimétrico, mientras el de
``\ln x_n`` se parece cada vez más a una campana.

Número de periodos ``n``: $(@bind n_factores Slider([1, 2, 5, 10, 25, 50, 100], default = 10, show_value = true))
"""

# ╔═╡ 7ad21e70-ca89-11f1-8522-b18628f9aa50
md"""
## 2. La ley de Gibrat no produce una ley de potencias

Simulemos 10 000 "empresas" que empiezan con tamaño 1 y en cada periodo se
multiplican por ``a_t = e^{\sigma Z_t}``, con ``Z_t`` normal estándar. La
gráfica muestra la CCDF de los tamaños en tres momentos.

Volatilidad del crecimiento ``\sigma``: $(@bind σ_gibrat Slider(0.05:0.05:0.5, default = 0.2, show_value = true))
"""

# ╔═╡ 7ad21e86-ca89-11f1-808e-31a402d2428f
md"""
Dos problemas para explicar la ley de Zipf:

1. En log-log, cada CCDF se **curva**: es una lognormal, no una ley de
   potencias.
2. La distribución **nunca se estabiliza**: la varianza de ``\ln x`` crece como
   ``\sigma^2 t``. En cambio, la distribución de tamaños de las ciudades es
   muy parecida de un censo a otro.

Falta algo que impida que los tamaños se dispersen para siempre.
"""

# ╔═╡ 7ad21e86-ca89-11f1-9940-531e9bcc71a3
md"""
## 3. Laboratorio: Gibrat con una barrera inferior

Agreguemos dos ingredientes a Gibrat:

- En promedio, las unidades **se encogen** un poco: ``\ln a_t`` tiene media
  ``m < 0``. (Si se mide el tamaño en relación con el total, que crece, la
  unidad típica se encoge.)
- Hay un **tamaño mínimo** ``x_{\min}``: una empresa no puede tener menos de un
  empleado; un pueblo que se encoge demasiado desaparece y es sustituido.

```math
x_{t+1} = \max\left(a_t\, x_t,\; x_{\min}\right), \qquad \ln a_t \sim \mathcal{N}(m, s^2).
```

El resultado, en el largo plazo, es una **ley de potencias** con exponente

```math
\alpha = -\frac{2m}{s^2}.
```

El encogimiento empuja hacia la barrera y el ruido aleja de ella; el
exponente mide el balance entre ambos.

Media del crecimiento logarítmico ``m``: $(@bind m_barrera Slider(-0.10:0.005:-0.005, default = -0.02, show_value = true))
  Volatilidad ``s``: $(@bind s_barrera Slider(0.10:0.02:0.40, default = 0.20, show_value = true))
  Mecanismo: $(@bind mecanismo Select(["Barrera inferior", "Kesten: x ← a·x + 1"]))
"""

# ╔═╡ 7ad21ea2-ca89-11f1-b04e-5f18bc5dfe00
md"""
**¿Cómo cambia el exponente con los parámetros?** La gráfica compara el
exponente estimado con el teórico para varias combinaciones de ``m`` y ``s``,
con los dos mecanismos.
"""

# ╔═╡ 7ad21eb4-ca89-11f1-b709-97418b51e0ba
md"""
Para la barrera, los puntos caen sobre la diagonal: ``\alpha = -2m/s^2``. Más
ruido o menos encogimiento dan colas **más pesadas** (``\alpha`` menor). Kesten
sigue la misma regla cuando ``\alpha`` es pequeño; con ``\alpha`` grande su cola
de potencias tarda en aparecer y la estimación con muestras finitas se queda
corta. Es una lección general: un exponente medido en datos finitos puede
diferir del asintótico.
"""

# ╔═╡ 7ad21eb4-ca89-11f1-8860-eb6bde617372
md"""
## 4. ¿Por qué ``\alpha \approx 1`` en las ciudades?

El mecanismo de la sección 3 produce cualquier exponente. ¿Por qué las
ciudades eligen justamente ``\alpha = 1``? Gabaix (1999) observó que si
medimos cada ciudad como **fracción de la población total**, el tamaño
promedio es fijo porque las fracciones suman 1. Con crecimiento de Gibrat y
una barrera ``x_{\min}`` muy pequeña frente al promedio ``\bar x``, el
exponente estacionario es

```math
\alpha = \frac{1}{1 - x_{\min}/\bar x} \;\longrightarrow\; 1
\qquad \text{cuando } x_{\min}/\bar x \to 0 .
```

La ley de Zipf surge cuando las unidades crecen en proporción a su tamaño y
el piso inferior es diminuto frente al promedio (en México,
``\hat\alpha = 0.99`` en la sesión 3). Pero que un mecanismo **pueda**
producir lo observado no prueba que sea **el** mecanismo: en la sección 0 vimos
que el crecimiento de los municipios no es exactamente proporcional, y otros
modelos también dan Zipf.
"""

# ╔═╡ 7ad21ec0-ca89-11f1-9d96-b34a1203700e
md"""
## 5. [Avanzado] GARCH es un proceso de Kesten

En la sesión 5 la varianza de GARCH(1,1) se actualizaba como
``\sigma_t^2 = \omega + a r_{t-1}^2 + b\sigma_{t-1}^2``. Como
``r_{t-1} = \sigma_{t-1}\varepsilon_{t-1}``,

```math
\sigma_t^2 = \underbrace{(a\varepsilon_{t-1}^2 + b)}_{\text{factor aleatorio}}\,\sigma_{t-1}^2 + \omega ,
```

que es un proceso de Kesten. Su cola es una ley de potencias con exponente
``\kappa`` que resuelve ``\mathbb{E}\big[(a\varepsilon^2 + b)^{\kappa}\big] = 1``, y la
de los retornos tiene exponente ``2\kappa`` (repaso E.A5).

Reacción ``a``: $(@bind a_k Slider(0.02:0.005:0.25, default = 0.075, show_value = true))
  Memoria ``b``: $(@bind b_k Slider(0.70:0.0025:0.97, default = 0.9125, show_value = true))
"""

# ╔═╡ 7ad21eca-ca89-11f1-b5cd-6db6bd2c5f54
md"""
## Resumen

| Modelo | Regla | Distribución de largo plazo |
|---|---|---|
| Gibrat puro | ``x_{t+1} = a_t x_t`` | Lognormal que se ensancha sin límite |
| Gibrat con barrera | ``x_{t+1} = \max(a_t x_t, x_{\min})`` | Ley de potencias, ``\alpha = -2m/s^2`` |
| Kesten | ``x_{t+1} = a_t x_t + b_t`` | Ley de potencias, ``\mathbb{E}[a^\alpha] = 1`` |
| Gibrat normalizado (Gabaix) | Fracciones del total con barrera pequeña | ``\alpha \to 1``: ley de Zipf |

Y en los datos: el crecimiento de los municipios de México cumple la parte de
Gibrat sobre la dispersión, pero no del todo la del promedio (urbanización).
"""

# ╔═╡ 7ad21ed4-ca89-11f1-8de3-ffe434df3403
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `paso_barrera(x, m, s, xmin, rng)`, que reciba un vector de
tamaños y devuelva los tamaños después de un periodo del modelo de la
sección 3. Genera los factores con `exp.(m .+ s .* randn(rng, length(x)))`.
"""

# ╔═╡ 7ad21ed4-ca89-11f1-82ab-a775181f6c26
function paso_barrera(x, m, s, xmin, rng)
    missing  # reemplaza esta línea
end

# ╔═╡ 7ad21ede-ca89-11f1-a446-1f469abb1237
md"""
**2.** Escribe `alpha_barrera(m, s)`, que devuelva el exponente teórico de la
sección 3.
"""

# ╔═╡ 7ad21ede-ca89-11f1-a933-03005a8e123a
function alpha_barrera(m, s)
    missing  # reemplaza esta línea
end

# ╔═╡ 7ad21ee6-ca89-11f1-9660-83e297c1fda8
md"""
### Ruta básica: interpretación

3. En la gráfica de la sección 0, ¿qué parte de la ley de Gibrat se cumple
   mejor? Propón una explicación económica para la que no se cumple.
4. En el modelo con barrera, ¿qué le pasa a la cola si la volatilidad ``s``
   se duplica? ¿Y si el encogimiento promedio ``|m|`` se duplica?
5. Explica con tus palabras por qué el crecimiento multiplicativo puro
   produce una lognormal cuya varianza crece sin límite, y por qué basta un
   piso inferior o un término aditivo (Kesten) para estabilizar una cola de
   ley de potencias.

### Ruta avanzada

6. Deriva ``\alpha = -2m/s^2`` para el modelo con barrera: escribe la
   ecuación de Fokker-Planck de ``y = \ln x`` (una difusión con deriva ``m``,
   varianza ``s^2`` y barrera reflejante), encuentra su solución estacionaria
   y transfórmala a ``x``.
7. Para Kesten con ``\ln a \sim \mathcal{N}(m, s^2)``, verifica que la
   condición ``\mathbb{E}[a^\alpha] = 1`` da exactamente ``\alpha = -2m/s^2``.
8. Simula el modelo de Gabaix: crecimiento con ``\mathbb{E}[a_t] = 1`` y
   barrera igual a una fracción ``c`` del promedio. Compara el exponente
   estimado con ``1/(1-c)`` para varios ``c``.

## Para la próxima sesión

- Lectura: Slanina (2014), sección 7.1. Complementaria: Slanina (2014),
  sección 7.3.
- En la sesión 7 pasamos de una serie a **muchas**: correlaciones entre
  activos, matrices aleatorias y la estructura de sectores del mercado.
"""

# ╔═╡ 7ad21ef2-ca89-11f1-9827-070b062f638e
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 7ad21efc-ca89-11f1-a298-33cd54b9945e
begin
    const AZUL    = "#2a78d6"   # serie principal
    const NARANJA = "#eb6834"   # segunda serie o teoría
    const GRIS    = "#52514e"   # referencias
    const RAMPA   = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # azul secuencial
end;

# ╔═╡ 7ad21e66-ca89-11f1-8657-11e4a01f59a0
let
    rng = Xoshiro(11)
    x = [prod(rand(rng, Uniform(0.7, 1.4), n_factores)) for _ in 1:10_000]
    p1 = histogram(x; normalize = :pdf, bins = 60, color = AZUL, lw = 0, alpha = 0.7,
                   xlabel = "xₙ", ylabel = "densidad", title = "Producto de $n_factores factores")
    lx = log.(x)
    p2 = histogram(lx; normalize = :pdf, bins = 60, color = AZUL, lw = 0, alpha = 0.7,
                   xlabel = "ln xₙ", title = "Su logaritmo")
    u = range(minimum(lx), maximum(lx), length = 200)
    plot!(p2, u, pdf.(Normal(mean(lx), std(lx)), u); color = NARANJA, lw = 2)
    plot(p1, p2; layout = (1, 2), size = (720, 300))
end

# ╔═╡ 7ad21efc-ca89-11f1-8f93-5b8e968462c1
"Número entero con comas como separador de miles."
miles(x) = replace(string(round(Int, x)), r"(?<=\d)(?=(\d{3})+$)" => ",")

# ╔═╡ 7ad21f06-ca89-11f1-9639-237eb0b020af
"Datos ordenados de mayor a menor y su CCDF empírica, P(X ≥ x)."
function ccdf_empirica(x)
    s = sort(x; rev = true)
    s, (1:length(s)) ./ length(s)
end

# ╔═╡ 7ad21f06-ca89-11f1-8cdd-7dfe5f302bc2
"Media y desviación estándar de `g` por decil de `x` (centro geométrico de cada decil)."
function por_deciles(x, g)
    cortes = quantile(x, 0:0.1:1)
    grupo = [clamp(searchsortedlast(cortes, xi), 1, 10) for xi in x]
    centro = [exp(mean(log.(x[grupo .== k]))) for k in 1:10]
    (centro = centro, media = [mean(g[grupo .== k]) for k in 1:10],
     sd = [std(g[grupo .== k]) for k in 1:10])
end

# ╔═╡ 7ad21e54-ca89-11f1-8538-47c8a0a74756
let
    d = por_deciles(p10, g)
    p1 = scatter(p10, 100 .* g; xscale = :log10, ms = 2, msw = 0, color = AZUL, alpha = 0.35,
                 xlabel = "población en 2010", ylabel = "crecimiento 2010–2020 (%)",
                 title = "Crecimiento de los $(miles(length(g))) municipios de México",
                 ylims = (-60, 100), label = "Municipios", legend = :topleft,
                 xticks = (10.0 .^ (2:7), ["100", "1 mil", "10 mil", "100 mil", "1 millón", "10 millones"]))
    plot!(p1, d.centro, 100 .* d.media; color = NARANJA, lw = 2, marker = :circle, ms = 4,
          msw = 0, label = "Promedio por decil de tamaño")
    hline!(p1, [0]; color = :black, lw = 0.5, label = "")
    p2 = plot(d.centro, 100 .* d.sd; xscale = :log10, color = NARANJA, lw = 2, marker = :circle,
              ms = 4, msw = 0, ylims = (0, 25), xlabel = "población en 2010",
              ylabel = "desviación estándar (%)", title = "Dispersión del crecimiento por decil",
              xticks = (10.0 .^ (2:7), ["100", "1 mil", "10 mil", "100 mil", "1 millón", "10 millones"]))
    plot(p1, p2; layout = (2, 1), size = (720, 560))
end

# ╔═╡ 7ad21f10-ca89-11f1-8728-fd2db4cb9ae7
"Gibrat puro: guarda los tamaños en los `tiempos` pedidos y la varianza de ln x en cada periodo."
function simular_gibrat(σ; N = 10_000, tiempos = [10, 100, 1000], rng = Xoshiro(12))
    lx = zeros(N)
    fotos = Vector{Vector{Float64}}()
    var_log = zeros(maximum(tiempos))
    for t in 1:maximum(tiempos)
        lx .+= σ .* randn(rng, N)
        var_log[t] = var(lx)
        t in tiempos && push!(fotos, exp.(lx))
    end
    (tiempos = tiempos, x = fotos, var_log = var_log)
end

# ╔═╡ 7ad21e7a-ca89-11f1-a182-9bd1c5cae975
gibrat_sim = simular_gibrat(σ_gibrat; N = 10_000, tiempos = [10, 100, 1000]);

# ╔═╡ 7ad21e7a-ca89-11f1-90ee-3b6584d59af9
let
    p1 = plot(xscale = :log10, yscale = :log10, ylims = (1e-4, 1.2), legend = :bottomleft,
              xlabel = "tamaño x", ylabel = "P(X ≥ x)", title = "Gibrat puro: CCDF en tres momentos")
    for (i, t) in enumerate(gibrat_sim.tiempos)
        s, P = ccdf_empirica(gibrat_sim.x[i])
        plot!(p1, s, P; color = RAMPA[i+1], lw = 2, label = "t = $t")
    end
    t = 1:1000
    p2 = plot(t, gibrat_sim.var_log; color = AZUL, lw = 2, label = "Simulación",
              xlabel = "periodo t", ylabel = "Var(ln x)", legend = :topleft,
              title = "La dispersión crece sin límite")
    plot!(p2, t, σ_gibrat^2 .* t; color = NARANJA, lw = 1.5, ls = :dash, label = "Teoría: σ² t")
    plot(p1, p2; layout = (1, 2), size = (900, 340))
end

# ╔═╡ 7ad21f10-ca89-11f1-844a-a3dadced708c
"Un periodo del modelo de Gibrat con barrera inferior."
paso_barrera_ref(x, m, s, xmin, rng) = max.(x .* exp.(m .+ s .* randn(rng, length(x))), xmin)

# ╔═╡ 7ad21f18-ca89-11f1-a03f-ad363bd04c71
"Gibrat con barrera inferior en x = 1, después de T periodos."
function simular_barrera(m, s; N = 10_000, T = 2_000, rng = Xoshiro(13))
    x = ones(N)
    for _ in 1:T
        x = paso_barrera_ref(x, m, s, 1.0, rng)
    end
    x
end

# ╔═╡ 7ad21f18-ca89-11f1-af90-774821049a4e
"Proceso de Kesten x ← a·x + 1, con ln a normal, después de T periodos."
function simular_kesten(m, s; N = 10_000, T = 2_000, rng = Xoshiro(14))
    x = ones(N)
    for _ in 1:T
        x .= x .* exp.(m .+ s .* randn(rng, N)) .+ 1
    end
    x
end

# ╔═╡ 7ad21e8e-ca89-11f1-a2ae-29bfe177d4da
sim_barrera = mecanismo == "Barrera inferior" ?
    simular_barrera(m_barrera, s_barrera; N = 10_000, T = 2_000) :
    simular_kesten(m_barrera, s_barrera; N = 10_000, T = 2_000);

# ╔═╡ 7ad21e98-ca89-11f1-940b-01fd9fe272f7
let
    α_teo = -2m_barrera / s_barrera^2
    s, P = ccdf_empirica(sim_barrera)
    p = scatter(s, P; xscale = :log10, yscale = :log10, ms = 1.5, msw = 0, color = AZUL,
                label = "Simulación: 10 000 unidades, 2000 periodos", legend = :bottomleft,
                xlabel = "tamaño x", ylabel = "P(X ≥ x)", ylims = (1e-4, 1.2),
                title = "$(mecanismo): distribución de largo plazo")
    x0 = quantile(sim_barrera, 0.9)            # la recta pasa por P = 0.1
    xr = exp10.(range(log10(x0), log10(maximum(s)), length = 100))
    plot!(p, xr, 0.1 .* (xr ./ x0) .^ (-α_teo); color = NARANJA, lw = 2, ls = :dash,
          label = "Teoría: pendiente −α = $(round(-α_teo, digits = 2))")
end

# ╔═╡ 7ad21f24-ca89-11f1-900e-e1f83e82d8a0
"""
Método de Clauset, Shalizi y Newman (2009): elige `xmin` minimizando la
distancia de Kolmogorov-Smirnov y estima α por máxima verosimilitud.
"""
function clauset(x; candidatos = quantile(x, 0.5:0.005:0.995), nmin = 50)
    s = sort(x; rev = true)
    L = cumsum(log.(s))
    res = (xmin = Float64[], α = Float64[], D = Float64[], n = Int[])
    for xm in unique(candidatos)
        n = searchsortedlast(s, xm; rev = true)
        n < nmin && continue
        α = n / (L[n] - n * log(xm))
        cola = @view s[n:-1:1]
        D = maximum(abs.((1:n) ./ n .- (1 .- (xm ./ cola) .^ α)))
        push!(res.xmin, xm); push!(res.α, α); push!(res.D, D); push!(res.n, n)
    end
    i = argmin(res.D)
    (; res..., mejor = (xmin = res.xmin[i], α = res.α[i], D = res.D[i], n = res.n[i],
                        se = res.α[i] / sqrt(res.n[i])))
end

# ╔═╡ 7ad21e8e-ca89-11f1-9363-9bd08e3437d7
ajuste_barrera = clauset(sim_barrera).mejor;

# ╔═╡ 7ad21ea2-ca89-11f1-8373-491b23bd139f
md"""
Con ``m = `` $(m_barrera) y ``s = `` $(s_barrera), la teoría predice
``\alpha = `` $(round(-2m_barrera / s_barrera^2, digits = 2)). Estimado con el
método de la sesión 3: ``\hat\alpha = `` $(round(ajuste_barrera.α, digits = 2)) ± $(round(1.96ajuste_barrera.se, digits = 2)).

$(-2m_barrera / s_barrera^2 > 5 ? md"Con ``\alpha`` tan grande la cola es muy delgada y hay pocos datos para estimarla: prueba con ``m`` más cercano a cero o con ``s`` mayor." : md"")

Cambia al mecanismo de **Kesten**, ``x_{t+1} = a_t x_t + 1``: en lugar de una
barrera, un pequeño ingreso fijo impide que ``x`` llegue a cero. El exponente
teórico es el mismo. La ley de potencias sólo aparece en la cola, lejos de 1,
así que con ``\alpha`` grande la estimación suele quedar por debajo de la
teoría.
"""

# ╔═╡ 7ad21f24-ca89-11f1-9e03-e1de5b14d3ea
"Exponente estimado contra teórico para varias combinaciones de (m, s), con ambos mecanismos."
function barrido_exponentes()
    combos = [(-0.01, 0.2), (-0.02, 0.2), (-0.03, 0.2), (-0.04, 0.2), (-0.03, 0.3),
              (-0.06, 0.3), (-0.09, 0.3), (-0.05, 0.25), (-0.06, 0.2)]
    teo = [-2m / s^2 for (m, s) in combos]
    barrera = [clauset(simular_barrera(m, s; N = 10_000, T = 2_000)).mejor.α for (m, s) in combos]
    kesten = [clauset(simular_kesten(m, s; N = 10_000, T = 2_000)).mejor.α for (m, s) in combos]
    (teo = teo, barrera = barrera, kesten = kesten)
end

# ╔═╡ 7ad21eac-ca89-11f1-a54e-e750fbc099fe
barrido = barrido_exponentes();

# ╔═╡ 7ad21eac-ca89-11f1-8f85-c150cc8d253d
let
    p = plot([0, 3.5], [0, 3.5]; color = GRIS, ls = :dash, lw = 1.5, label = "Estimado = teórico",
             xlabel = "α teórico = −2m/s²", ylabel = "α estimado", legend = :topleft,
             title = "El exponente lo fija el balance entre deriva y ruido",
             xlims = (0, 3.5), ylims = (0, 3.5), aspect_ratio = :equal, size = (520, 480))
    scatter!(p, barrido.teo, barrido.barrera; color = AZUL, ms = 6, msw = 0, label = "Barrera inferior")
    scatter!(p, barrido.teo, barrido.kesten; color = NARANJA, ms = 6, msw = 0, marker = :diamond,
             label = "Kesten")
end

# ╔═╡ 7ad21f2e-ca89-11f1-b3df-efdb84ea03f9
"""
Exponente κ de la cola de σ² en GARCH(1,1): resuelve E[(aZ² + b)^κ] = 1 con Z
normal estándar, por bisección y cuadratura en una malla fina.
"""
function kappa_garch(a, b)
    a + b >= 1 && return NaN
    z = range(-12, 12, length = 40_001)
    w = pdf.(Normal(), z) .* step(z)
    E(κ) = sum(w .* (a .* z .^ 2 .+ b) .^ κ)
    lo, hi = 1.0, 60.0
    for _ in 1:60
        mid = (lo + hi) / 2
        E(mid) > 1 ? (hi = mid) : (lo = mid)
    end
    (lo + hi) / 2
end

# ╔═╡ 7ad21eca-ca89-11f1-b64d-6bd76e3162db
md"""
$(a_k + b_k >= 1 ? md"Con ``a + b \ge 1`` el proceso no es estacionario." : md"Con ``a = `` $(a_k) y ``b = `` $(b_k), el exponente de la cola de los retornos es ``2\kappa = `` $(round(2kappa_garch(a_k, b_k), digits = 2)).")

Con los parámetros estimados en la sesión 5:

| Activo | ``a`` | ``b`` | ``\alpha`` implícito en GARCH | ``\alpha`` de los datos (sesión 3) |
|---|---|---|---|---|
| IPC (BMV) | 0.075 | 0.9125 | $(round(2kappa_garch(0.075, 0.9125), digits = 1)) | 2.8 |
| USD/MXN | 0.12 | 0.8575 | $(round(2kappa_garch(0.12, 0.8575), digits = 1)) | 3.1 |
| S&P 500 | 0.125 | 0.8575 | $(round(2kappa_garch(0.125, 0.8575), digits = 1)) | 3.1 |
| Bitcoin | 0.11 | 0.8475 | $(round(2kappa_garch(0.11, 0.8475), digits = 1)) | 4.0 |

GARCH tiene colas de potencia, pero su exponente asintótico es mayor que el
de los datos: sus colas son más ligeras. Es la explicación teórica de lo que
vimos en el laboratorio de la sesión 5.
"""

# ╔═╡ 7ad21f38-ca89-11f1-ba6a-17af166762fd
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

# ╔═╡ 7ad21ede-ca89-11f1-b1ca-799cba1a7629
verificar(() -> paso_barrera(fill(2.0, 1000), -0.05, 0.3, 1.0, Xoshiro(1)),
          paso_barrera_ref(fill(2.0, 1000), -0.05, 0.3, 1.0, Xoshiro(1));
          pista = md"Multiplica `x` por los factores elemento a elemento (`.*`) y aplica `max.(…, xmin)` para que nadie quede por debajo de la barrera.")

# ╔═╡ 7ad21ee6-ca89-11f1-a573-6bdb0c31d93c
verificar(() -> alpha_barrera(-0.03, 0.25), 0.96;
          pista = md"La fórmula es ``\alpha = -2m/s^2``.")

# ╔═╡ Cell order:
# ╟─7ad21df8-ca89-11f1-97ba-cdac1f0b5427
# ╟─7ad21e2a-ca89-11f1-8f14-0fb897c70423
# ╟─7ad21e34-ca89-11f1-8dff-b3dc92cb2354
# ╠═7ad21e3e-ca89-11f1-80aa-3f29a0e6eecc
# ╠═7ad21e3e-ca89-11f1-bfe2-41c88266bb20
# ╟─7ad21e54-ca89-11f1-8538-47c8a0a74756
# ╟─7ad21e54-ca89-11f1-9fcb-e733dbcb5b8f
# ╟─7ad21e5c-ca89-11f1-b890-57b1e87d312b
# ╟─7ad21e66-ca89-11f1-bd6f-cbeddd99a19c
# ╟─7ad21e66-ca89-11f1-8657-11e4a01f59a0
# ╟─7ad21e70-ca89-11f1-8522-b18628f9aa50
# ╟─7ad21e7a-ca89-11f1-a182-9bd1c5cae975
# ╟─7ad21e7a-ca89-11f1-90ee-3b6584d59af9
# ╟─7ad21e86-ca89-11f1-808e-31a402d2428f
# ╟─7ad21e86-ca89-11f1-9940-531e9bcc71a3
# ╟─7ad21e8e-ca89-11f1-a2ae-29bfe177d4da
# ╟─7ad21e8e-ca89-11f1-9363-9bd08e3437d7
# ╟─7ad21e98-ca89-11f1-940b-01fd9fe272f7
# ╟─7ad21ea2-ca89-11f1-8373-491b23bd139f
# ╟─7ad21ea2-ca89-11f1-b04e-5f18bc5dfe00
# ╟─7ad21eac-ca89-11f1-a54e-e750fbc099fe
# ╟─7ad21eac-ca89-11f1-8f85-c150cc8d253d
# ╟─7ad21eb4-ca89-11f1-b709-97418b51e0ba
# ╟─7ad21eb4-ca89-11f1-8860-eb6bde617372
# ╟─7ad21ec0-ca89-11f1-9d96-b34a1203700e
# ╟─7ad21eca-ca89-11f1-b64d-6bd76e3162db
# ╟─7ad21eca-ca89-11f1-b5cd-6db6bd2c5f54
# ╟─7ad21ed4-ca89-11f1-8de3-ffe434df3403
# ╠═7ad21ed4-ca89-11f1-82ab-a775181f6c26
# ╟─7ad21ede-ca89-11f1-b1ca-799cba1a7629
# ╟─7ad21ede-ca89-11f1-a446-1f469abb1237
# ╠═7ad21ede-ca89-11f1-a933-03005a8e123a
# ╟─7ad21ee6-ca89-11f1-a573-6bdb0c31d93c
# ╟─7ad21ee6-ca89-11f1-9660-83e297c1fda8
# ╟─7ad21ef2-ca89-11f1-9827-070b062f638e
# ╠═7ad21efc-ca89-11f1-a183-4f3484351f7e
# ╠═7ad21efc-ca89-11f1-a298-33cd54b9945e
# ╠═7ad21efc-ca89-11f1-8f93-5b8e968462c1
# ╠═7ad21f06-ca89-11f1-9639-237eb0b020af
# ╠═7ad21f06-ca89-11f1-8cdd-7dfe5f302bc2
# ╠═7ad21f10-ca89-11f1-8728-fd2db4cb9ae7
# ╠═7ad21f10-ca89-11f1-844a-a3dadced708c
# ╠═7ad21f18-ca89-11f1-a03f-ad363bd04c71
# ╠═7ad21f18-ca89-11f1-af90-774821049a4e
# ╠═7ad21f24-ca89-11f1-900e-e1f83e82d8a0
# ╠═7ad21f24-ca89-11f1-9e03-e1de5b14d3ea
# ╠═7ad21f2e-ca89-11f1-b3df-efdb84ea03f9
# ╠═7ad21f38-ca89-11f1-ba6a-17af166762fd
