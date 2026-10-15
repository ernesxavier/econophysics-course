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

# ╔═╡ 0edc3240-c8d2-11f1-8d61-95b560a095b4
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

# ╔═╡ 0edc313c-c8d2-11f1-aa2b-212823c7b0ba
md"""
# Sesión 4: Crecimiento multiplicativo

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
`Code/s04_crecimiento_multiplicativo.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): C.B1 (logaritmos), P.I3
(teorema del límite central), P.I6 (la lognormal), E.I6 (Gibrat y Kesten),
P.I5 (leyes de potencias). Ruta avanzada: E.A4 y E.A5.
"""

# ╔═╡ 0edc316e-c8d2-11f1-b655-71d6be4d9591
TableOfContents(title = "Contenido")

# ╔═╡ 0edc316e-c8d2-11f1-8250-df4a23271648
md"""
## 0. El fenómeno: ¿crecen igual los grandes y los pequeños?

En la sesión 2 vimos que el tamaño de las localidades de México sigue la ley
de Zipf. Hoy buscamos un **mecanismo** que la produzca. La hipótesis más
simple es la de Robert Gibrat (1931): cada unidad (una ciudad, una empresa)
crece en una **proporción aleatoria** que no depende de su tamaño. Una ciudad
de un millón y un pueblo de mil habitantes tendrían la misma probabilidad de
crecer 10%.

Podemos ponerla a prueba con los censos de 2010 y 2020 (INEGI). Los datos
están en `Data/municipios_2010_2020.csv` y se regeneran con
`Code/descargar_inegi.jl`.
"""

# ╔═╡ 0edc3178-c8d2-11f1-ada8-0f7a664e1761
municipios = CSV.read(joinpath(@__DIR__, "..", "Data", "municipios_2010_2020.csv"), DataFrame)

# ╔═╡ 0edc3182-c8d2-11f1-b1e3-d714185a46f1
begin
    p10 = Float64.(municipios.poblacion_2010)
    p20 = Float64.(municipios.poblacion_2020)
    g = log.(p20 ./ p10)          # tasa de crecimiento logarítmica 2010–2020
end;

# ╔═╡ 0edc3196-c8d2-11f1-9ec6-afa935fbfeae
gibrat_reg = let x = log.(p10)
    b = cov(x, g) / var(x)
    (b = b, R2 = cor(x, g)^2)
end;

# ╔═╡ 0edc318e-c8d2-11f1-9ffd-9b2a60fa3a18
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

# ╔═╡ 0edc3196-c8d2-11f1-b86a-3b7deeb2e886
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

# ╔═╡ 0edc31aa-c8d2-11f1-954f-0b4ee40379b9
md"""
## 2. La ley de Gibrat no produce una ley de potencias

Simulemos 10 000 "empresas" que empiezan con tamaño 1 y en cada periodo se
multiplican por ``a_t = e^{\sigma Z_t}``, con ``Z_t`` normal estándar. La
gráfica muestra la CCDF de los tamaños en tres momentos.

Volatilidad del crecimiento ``\sigma``: $(@bind σ_gibrat Slider(0.05:0.05:0.5, default = 0.2, show_value = true))
"""

# ╔═╡ 0edc31bc-c8d2-11f1-bffc-2df7c1419ecd
md"""
Dos problemas para explicar la ley de Zipf:

1. En log-log, cada CCDF se **curva**: es una lognormal, no una ley de
   potencias.
2. La distribución **nunca se estabiliza**: la varianza de ``\ln x`` crece como
   ``\sigma^2 t``. En cambio, la distribución de tamaños de las ciudades es
   muy parecida de un censo a otro.

Falta algo que impida que los tamaños se dispersen para siempre.
"""

# ╔═╡ 0edc31c8-c8d2-11f1-82db-177a41670bbd
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

# ╔═╡ 0edc31e6-c8d2-11f1-8b17-dd1ae445cc4f
md"""
**¿Cómo cambia el exponente con los parámetros?** La gráfica compara el
exponente estimado con el teórico para varias combinaciones de ``m`` y ``s``,
con los dos mecanismos.
"""

# ╔═╡ 0edc31ee-c8d2-11f1-9700-0101f9136fef
md"""
Para la barrera, los puntos caen sobre la diagonal: ``\alpha = -2m/s^2``. Más
ruido o menos encogimiento dan colas **más pesadas** (``\alpha`` menor). Kesten
sigue la misma regla cuando ``\alpha`` es pequeño; con ``\alpha`` grande su cola
de potencias tarda en aparecer y la estimación con muestras finitas se queda
corta. Es una lección general: un exponente medido en datos finitos puede
diferir del asintótico.
"""

# ╔═╡ 0edc31fa-c8d2-11f1-b817-672aeda9b2d4
md"""
## 4. ¿Por qué ``\alpha \approx 1`` en las ciudades?

El mecanismo de la sección 3 produce cualquier exponente. ¿Por qué las
ciudades eligen justamente ``\alpha = 1``? Gabaix (1999) observó que hay una
restricción adicional: si medimos cada ciudad como **fracción de la población
total**, el tamaño promedio no puede crecer ni encogerse, porque las
fracciones suman 1. Con crecimiento de Gibrat, esa restricción y una barrera
``x_{\min}`` muy por debajo del tamaño promedio ``\bar x`` dan

```math
\alpha = \frac{1}{1 - x_{\min}/\bar x} \;\longrightarrow\; 1
\qquad \text{cuando } x_{\min}/\bar x \to 0 .
```

La ley de Zipf sería, entonces, la huella de dos hechos: las ciudades crecen
aproximadamente en proporción a su tamaño (sección 0) y la ciudad más pequeña
es diminuta comparada con la promedio.

Las localidades de México dieron ``\hat\alpha = 0.99`` (sesión 2), consistente
con este argumento. Pero ojo: que un mecanismo **pueda** producir lo observado
no prueba que sea **el** mecanismo. En la sección 0 vimos que el crecimiento
de los municipios no es exactamente proporcional, y hay otros modelos que
también dan Zipf.
"""

# ╔═╡ 0edc3204-c8d2-11f1-bf38-49d45f9f1864
md"""
## 5. [Avanzado] GARCH es un proceso de Kesten

En la sesión 3 la varianza de GARCH(1,1) se actualizaba como
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

# ╔═╡ 0edc320e-c8d2-11f1-9058-4baa57cfb683
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

# ╔═╡ 0edc320e-c8d2-11f1-b99c-63a280fd3757
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `paso_barrera(x, m, s, xmin, rng)`, que reciba un vector de
tamaños y devuelva los tamaños después de un periodo del modelo de la
sección 3. Genera los factores con `exp.(m .+ s .* randn(rng, length(x)))`.
"""

# ╔═╡ 0edc3218-c8d2-11f1-863c-c1f7746bf304
function paso_barrera(x, m, s, xmin, rng)
    missing  # reemplaza esta línea
end

# ╔═╡ 0edc3220-c8d2-11f1-8661-17e027c5cc96
md"""
**2.** Escribe `alpha_barrera(m, s)`, que devuelva el exponente teórico de la
sección 3.
"""

# ╔═╡ 0edc3220-c8d2-11f1-8f79-9fb6b2eb5be8
function alpha_barrera(m, s)
    missing  # reemplaza esta línea
end

# ╔═╡ 0edc322c-c8d2-11f1-8827-d9f5d33506c4
md"""
### Ruta básica: interpretación

3. En la gráfica de la sección 0, ¿qué parte de la ley de Gibrat se cumple
   mejor? Propón una explicación económica para la que no se cumple.
4. Si una variable es el producto de muchos factores aleatorios
   independientes, ¿qué distribución esperas? ¿Por qué no una normal?
5. En el modelo con barrera, ¿qué le pasa a la cola si la volatilidad ``s``
   se duplica? ¿Y si el encogimiento promedio ``|m|`` se duplica?
6. Explica con tus palabras por qué hace falta una barrera (o un ingreso
   fijo) para obtener una ley de potencias estable.
7. Según el argumento de Gabaix, ¿qué esperarías del exponente de Zipf en un
   país donde la ciudad más pequeña considerada tiene un tamaño cercano al
   promedio?

### Ruta avanzada

8. Deriva ``\alpha = -2m/s^2`` para el modelo con barrera: escribe la
   ecuación de Fokker-Planck de ``y = \ln x`` (una difusión con deriva ``m``,
   varianza ``s^2`` y barrera reflejante), encuentra su solución estacionaria
   y transfórmala a ``x``.
9. Para Kesten con ``\ln a \sim \mathcal{N}(m, s^2)``, verifica que la
   condición ``\mathbb{E}[a^\alpha] = 1`` da exactamente ``\alpha = -2m/s^2``.
10. Simula el modelo de Gabaix: crecimiento con ``\mathbb{E}[a_t] = 1`` y
    barrera igual a una fracción ``c`` del promedio. Compara el exponente
    estimado con ``1/(1-c)`` para varios ``c``. ¿Qué pasa cuando ``c`` es muy
    pequeño y por qué influye el número de ciudades?
11. Regresa el crecimiento de los municipios sobre el logaritmo de su tamaño
    agregando como control el estado. ¿Cambia la pendiente?

## Para la próxima sesión

- Lectura: Slanina (2014), secciones 7.1 y 7.3; Gabaix (2009), sección 3.
- En la semana 5 pasamos de una serie a **muchas**: correlaciones entre
  activos, matrices aleatorias y la estructura de sectores del mercado.
"""

# ╔═╡ 0edc3236-c8d2-11f1-a5d9-bb243ebf47cf
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 0edc3240-c8d2-11f1-981f-eb7def6ab448
begin
    const AZUL    = "#2a78d6"   # serie principal
    const NARANJA = "#eb6834"   # segunda serie o teoría
    const GRIS    = "#52514e"   # referencias
    const RAMPA   = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # azul secuencial
end;

# ╔═╡ 0edc31aa-c8d2-11f1-af14-25eb7a494bec
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

# ╔═╡ 0edc324a-c8d2-11f1-95f7-37196b4ff44a
"Número entero con comas como separador de miles."
miles(x) = replace(string(round(Int, x)), r"(?<=\d)(?=(\d{3})+$)" => ",")

# ╔═╡ 0edc3254-c8d2-11f1-8ae7-0f82db3d01eb
"Datos ordenados de mayor a menor y su CCDF empírica, P(X ≥ x)."
function ccdf_empirica(x)
    s = sort(x; rev = true)
    s, (1:length(s)) ./ length(s)
end

# ╔═╡ 0edc3254-c8d2-11f1-9f68-291d303cdce2
"Media y desviación estándar de `g` por decil de `x` (centro geométrico de cada decil)."
function por_deciles(x, g)
    cortes = quantile(x, 0:0.1:1)
    grupo = [clamp(searchsortedlast(cortes, xi), 1, 10) for xi in x]
    centro = [exp(mean(log.(x[grupo .== k]))) for k in 1:10]
    (centro = centro, media = [mean(g[grupo .== k]) for k in 1:10],
     sd = [std(g[grupo .== k]) for k in 1:10])
end

# ╔═╡ 0edc3182-c8d2-11f1-9608-659e13a209aa
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

# ╔═╡ 0edc325e-c8d2-11f1-a3cf-cf1964f7e6af
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

# ╔═╡ 0edc31b4-c8d2-11f1-b273-57506488d67e
gibrat_sim = simular_gibrat(σ_gibrat; N = 10_000, tiempos = [10, 100, 1000]);

# ╔═╡ 0edc31bc-c8d2-11f1-bbe8-33a1e680f1e2
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

# ╔═╡ 0edc325e-c8d2-11f1-9eb1-13290c975554
"Un periodo del modelo de Gibrat con barrera inferior."
paso_barrera_ref(x, m, s, xmin, rng) = max.(x .* exp.(m .+ s .* randn(rng, length(x))), xmin)

# ╔═╡ 0edc3268-c8d2-11f1-8e2a-69d76b3da82d
"Gibrat con barrera inferior en x = 1, después de T periodos."
function simular_barrera(m, s; N = 10_000, T = 2_000, rng = Xoshiro(13))
    x = ones(N)
    for _ in 1:T
        x = paso_barrera_ref(x, m, s, 1.0, rng)
    end
    x
end

# ╔═╡ 0edc3268-c8d2-11f1-8cd3-3d4d9cc7ccdd
"Proceso de Kesten x ← a·x + 1, con ln a normal, después de T periodos."
function simular_kesten(m, s; N = 10_000, T = 2_000, rng = Xoshiro(14))
    x = ones(N)
    for _ in 1:T
        x .= x .* exp.(m .+ s .* randn(rng, N)) .+ 1
    end
    x
end

# ╔═╡ 0edc31d2-c8d2-11f1-840c-65ed1f2798fa
sim_barrera = mecanismo == "Barrera inferior" ?
    simular_barrera(m_barrera, s_barrera; N = 10_000, T = 2_000) :
    simular_kesten(m_barrera, s_barrera; N = 10_000, T = 2_000);

# ╔═╡ 0edc31dc-c8d2-11f1-9c08-7fa50eae8867
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

# ╔═╡ 0edc3272-c8d2-11f1-b839-63354214c1b7
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

# ╔═╡ 0edc31d2-c8d2-11f1-9c02-019bec2ee41a
ajuste_barrera = clauset(sim_barrera).mejor;

# ╔═╡ 0edc31dc-c8d2-11f1-b5ae-dbec35d35acf
md"""
Con ``m = `` $(m_barrera) y ``s = `` $(s_barrera), la teoría predice
``\alpha = `` $(round(-2m_barrera / s_barrera^2, digits = 2)). Estimado con el
método de la sesión 2: ``\hat\alpha = `` $(round(ajuste_barrera.α, digits = 2)) ± $(round(1.96ajuste_barrera.se, digits = 2)).

$(-2m_barrera / s_barrera^2 > 5 ? md"Con ``\alpha`` tan grande la cola es muy delgada y hay pocos datos para estimarla: prueba con ``m`` más cercano a cero o con ``s`` mayor." : md"")

Cambia al mecanismo de **Kesten**, ``x_{t+1} = a_t x_t + 1``: en lugar de una
barrera, un pequeño ingreso fijo impide que ``x`` llegue a cero. El exponente
teórico es el mismo. La ley de potencias sólo aparece en la cola, lejos de 1,
así que con ``\alpha`` grande la estimación suele quedar por debajo de la
teoría.
"""

# ╔═╡ 0edc3272-c8d2-11f1-b1e1-312ac2a1caf1
"Exponente estimado contra teórico para varias combinaciones de (m, s), con ambos mecanismos."
function barrido_exponentes()
    combos = [(-0.01, 0.2), (-0.02, 0.2), (-0.03, 0.2), (-0.04, 0.2), (-0.03, 0.3),
              (-0.06, 0.3), (-0.09, 0.3), (-0.05, 0.25), (-0.06, 0.2)]
    teo = [-2m / s^2 for (m, s) in combos]
    barrera = [clauset(simular_barrera(m, s; N = 10_000, T = 2_000)).mejor.α for (m, s) in combos]
    kesten = [clauset(simular_kesten(m, s; N = 10_000, T = 2_000)).mejor.α for (m, s) in combos]
    (teo = teo, barrera = barrera, kesten = kesten)
end

# ╔═╡ 0edc31e6-c8d2-11f1-b7a6-5d2662a140ab
barrido = barrido_exponentes();

# ╔═╡ 0edc31ee-c8d2-11f1-8c54-af67065d9748
let
    p = plot([0, 3.5], [0, 3.5]; color = GRIS, ls = :dash, lw = 1.5, label = "Estimado = teórico",
             xlabel = "α teórico = −2m/s²", ylabel = "α estimado", legend = :topleft,
             title = "El exponente lo fija el balance entre deriva y ruido",
             xlims = (0, 3.5), ylims = (0, 3.5), aspect_ratio = :equal, size = (520, 480))
    scatter!(p, barrido.teo, barrido.barrera; color = AZUL, ms = 6, msw = 0, label = "Barrera inferior")
    scatter!(p, barrido.teo, barrido.kesten; color = NARANJA, ms = 6, msw = 0, marker = :diamond,
             label = "Kesten")
end

# ╔═╡ 0edc327c-c8d2-11f1-9936-27888f32c63d
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

# ╔═╡ 0edc3204-c8d2-11f1-8ac4-ff18b28b1c01
md"""
$(a_k + b_k >= 1 ? md"Con ``a + b \ge 1`` el proceso no es estacionario." : md"Con ``a = `` $(a_k) y ``b = `` $(b_k), el exponente de la cola de los retornos es ``2\kappa = `` $(round(2kappa_garch(a_k, b_k), digits = 2)).")

Con los parámetros estimados en la sesión 3:

| Activo | ``a`` | ``b`` | ``\alpha`` implícito en GARCH | ``\alpha`` de los datos (sesión 2) |
|---|---|---|---|---|
| IPC (BMV) | 0.075 | 0.9125 | $(round(2kappa_garch(0.075, 0.9125), digits = 1)) | 2.8 |
| USD/MXN | 0.12 | 0.8575 | $(round(2kappa_garch(0.12, 0.8575), digits = 1)) | 3.1 |
| S&P 500 | 0.125 | 0.8575 | $(round(2kappa_garch(0.125, 0.8575), digits = 1)) | 3.1 |
| Bitcoin | 0.11 | 0.8475 | $(round(2kappa_garch(0.11, 0.8475), digits = 1)) | 4.0 |

GARCH tiene colas de potencia, pero su exponente asintótico es mayor que el
de los datos: sus colas son más ligeras. Es la explicación teórica de lo que
vimos en el laboratorio de la sesión 3.
"""

# ╔═╡ 0edc3286-c8d2-11f1-bb54-53c3e97f1a7a
"Muestra si la respuesta de un ejercicio es correcta."
function verificar(f, esperado; rtol = 1e-8, pista = md"")
    caja(tipo, titulo, texto) = Markdown.MD(Markdown.Admonition(tipo, titulo, [texto]))
    resultado = try
        f()
    catch e
        return caja("danger", "Tu código produce un error", md"$(sprint(showerror, e))")
    end
    if resultado === missing
        caja("warning", "Pendiente", md"Reemplaza `missing` por tu código.")
    elseif resultado isa typeof(esperado) && size(resultado) == size(esperado) &&
           all(isapprox.(resultado, esperado; rtol, atol = 1e-12))
        caja("correct", "¡Correcto!", md"Bien hecho.")
    else
        caja("danger", "Todavía no", pista)
    end
end

# ╔═╡ 0edc3218-c8d2-11f1-b932-0511bac219f2
verificar(() -> paso_barrera(fill(2.0, 1000), -0.05, 0.3, 1.0, Xoshiro(1)),
          paso_barrera_ref(fill(2.0, 1000), -0.05, 0.3, 1.0, Xoshiro(1));
          pista = md"Multiplica `x` por los factores elemento a elemento (`.*`) y aplica `max.(…, xmin)` para que nadie quede por debajo de la barrera.")

# ╔═╡ 0edc322c-c8d2-11f1-91e4-2f14de0b56e1
verificar(() -> alpha_barrera(-0.03, 0.25), 0.96;
          pista = md"La fórmula es ``\alpha = -2m/s^2``.")

# ╔═╡ Cell order:
# ╟─0edc313c-c8d2-11f1-aa2b-212823c7b0ba
# ╟─0edc316e-c8d2-11f1-b655-71d6be4d9591
# ╟─0edc316e-c8d2-11f1-8250-df4a23271648
# ╠═0edc3178-c8d2-11f1-ada8-0f7a664e1761
# ╠═0edc3182-c8d2-11f1-b1e3-d714185a46f1
# ╟─0edc3182-c8d2-11f1-9608-659e13a209aa
# ╟─0edc318e-c8d2-11f1-9ffd-9b2a60fa3a18
# ╟─0edc3196-c8d2-11f1-9ec6-afa935fbfeae
# ╟─0edc3196-c8d2-11f1-b86a-3b7deeb2e886
# ╟─0edc31aa-c8d2-11f1-af14-25eb7a494bec
# ╟─0edc31aa-c8d2-11f1-954f-0b4ee40379b9
# ╟─0edc31b4-c8d2-11f1-b273-57506488d67e
# ╟─0edc31bc-c8d2-11f1-bbe8-33a1e680f1e2
# ╟─0edc31bc-c8d2-11f1-bffc-2df7c1419ecd
# ╟─0edc31c8-c8d2-11f1-82db-177a41670bbd
# ╟─0edc31d2-c8d2-11f1-840c-65ed1f2798fa
# ╟─0edc31d2-c8d2-11f1-9c02-019bec2ee41a
# ╟─0edc31dc-c8d2-11f1-9c08-7fa50eae8867
# ╟─0edc31dc-c8d2-11f1-b5ae-dbec35d35acf
# ╟─0edc31e6-c8d2-11f1-8b17-dd1ae445cc4f
# ╟─0edc31e6-c8d2-11f1-b7a6-5d2662a140ab
# ╟─0edc31ee-c8d2-11f1-8c54-af67065d9748
# ╟─0edc31ee-c8d2-11f1-9700-0101f9136fef
# ╟─0edc31fa-c8d2-11f1-b817-672aeda9b2d4
# ╟─0edc3204-c8d2-11f1-bf38-49d45f9f1864
# ╟─0edc3204-c8d2-11f1-8ac4-ff18b28b1c01
# ╟─0edc320e-c8d2-11f1-9058-4baa57cfb683
# ╟─0edc320e-c8d2-11f1-b99c-63a280fd3757
# ╠═0edc3218-c8d2-11f1-863c-c1f7746bf304
# ╟─0edc3218-c8d2-11f1-b932-0511bac219f2
# ╟─0edc3220-c8d2-11f1-8661-17e027c5cc96
# ╠═0edc3220-c8d2-11f1-8f79-9fb6b2eb5be8
# ╟─0edc322c-c8d2-11f1-91e4-2f14de0b56e1
# ╟─0edc322c-c8d2-11f1-8827-d9f5d33506c4
# ╟─0edc3236-c8d2-11f1-a5d9-bb243ebf47cf
# ╠═0edc3240-c8d2-11f1-8d61-95b560a095b4
# ╠═0edc3240-c8d2-11f1-981f-eb7def6ab448
# ╠═0edc324a-c8d2-11f1-95f7-37196b4ff44a
# ╠═0edc3254-c8d2-11f1-8ae7-0f82db3d01eb
# ╠═0edc3254-c8d2-11f1-9f68-291d303cdce2
# ╠═0edc325e-c8d2-11f1-a3cf-cf1964f7e6af
# ╠═0edc325e-c8d2-11f1-9eb1-13290c975554
# ╠═0edc3268-c8d2-11f1-8e2a-69d76b3da82d
# ╠═0edc3268-c8d2-11f1-8cd3-3d4d9cc7ccdd
# ╠═0edc3272-c8d2-11f1-b839-63354214c1b7
# ╠═0edc3272-c8d2-11f1-b1e1-312ac2a1caf1
# ╠═0edc327c-c8d2-11f1-9936-27888f32c63d
# ╠═0edc3286-c8d2-11f1-bb54-53c3e97f1a7a
