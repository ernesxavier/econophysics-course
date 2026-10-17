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

# ╔═╡ 16c83e7e-ca79-11f1-9ffb-1f0c5e61f024
begin
    # Usa el entorno del curso (Code/Project.toml) en lugar del gestor de
    # paquetes propio de Pluto: un solo entorno para todo el repositorio.
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, CSV, DataFrames, Dates, Statistics, StatsBase, Distributions,
          Random, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

# ╔═╡ 16c83d70-ca79-11f1-b8ac-d135e84b438d
md"""
# Sesión 2: Volatilidad y modelos nulos

**Objetivos de la sesión**
1. Distinguir **incorrelación** de **independencia**.
2. Medir los otros dos **hechos estilizados**: **agrupamiento de volatilidad**
   y **normalidad agregacional**.
3. Usar **modelos nulos** (datos barajados y caminata gaussiana) para saber
   dónde vive cada hecho estilizado.
4. Comprobar que los cuatro hechos aparecen en mercados muy distintos
   (**universalidad**).

**Cómo usar este cuaderno.** Pluto es *reactivo*: al mover un deslizador o
cambiar un selector, todas las celdas que dependen de ese valor se recalculan
solas. Las secciones marcadas **[Avanzado]** son opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s02_volatilidad.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): P.B3–B4 (histogramas y la
normal), S.B5 (permutaciones). Para las notas formales: P.I1 (independencia y
correlación), P.I3 (teorema del límite central), E.I1 (autocorrelación), S.I4
(datos dependientes). Ruta avanzada: P.A1 (cumulantes), F.A2
(universalidad), S.A4 (memoria larga).
"""

# ╔═╡ 16c83dbe-ca79-11f1-8ba9-65c659edc421
TableOfContents(title = "Contenido")

# ╔═╡ 16c83dca-ca79-11f1-b855-75e62f6ac58a
md"""
## 0. El fenómeno: los días malos llegan juntos

En la sesión 1 vimos que la **dirección** del mercado casi no tiene memoria:
saber si hoy subió no ayuda a predecir si mañana sube. Pero mira el **tamaño**
de los movimientos del S&P 500:
"""

# ╔═╡ 16c83de8-ca79-11f1-bbf5-0f12ea464605
md"""
## 1. Repaso: incorrelación no es independencia

Dos variables son **independientes** si conocer una no dice nada de la otra.
Son **incorreladas** si no tienen relación *lineal*: ``\mathrm{corr}(X, Y) = 0``.
Independencia implica incorrelación, pero no al revés. El ejemplo clásico:
``X`` normal estándar e ``Y = X^2``.
"""

# ╔═╡ 16c83df2-ca79-11f1-aa42-4176e97203ad
md"""
La correlación es casi cero porque los valores negativos y positivos de ``X``
dan el mismo ``Y``: no hay una recta que describa la relación. Pero ``Y`` está
completamente determinada por ``X``. Con los retornos pasa algo parecido: el
signo de ``r_t`` no predice el de ``r_{t+1}``, pero el **tamaño** ``|r_t|`` sí
dice algo sobre ``|r_{t+1}|``.

## 2. Los datos

Los mismos cuatro activos de la sesión 1 (`Data/precios.csv`).
"""

# ╔═╡ 16c83dfc-ca79-11f1-a3e7-11b82d4a537e
datos = CSV.read(joinpath(@__DIR__, "..", "Data", "precios.csv"), DataFrame)

# ╔═╡ 16c83dfc-ca79-11f1-8f87-4392ed895cbd
md"""
**Activo:** $(@bind activo Select(unique(datos.activo)))
"""

# ╔═╡ 16c83dfc-ca79-11f1-9681-eb413b011b73
begin
    serie = sort(filter(:activo => ==(activo), datos), :fecha)
    fechas = serie.fecha
    precios = serie.cierre
    r = diff(log.(precios))      # retornos logarítmicos diarios
    fechas_r = fechas[2:end]
    μ, σ = mean(r), std(r)
end;

# ╔═╡ 16c83e06-ca79-11f1-86bf-0dd372c0a900
md"""
## 3. Hecho 3: agrupamiento de volatilidad

Ahora miramos el **tamaño** del movimiento, ``|r_t|``, sin importar si fue
subida o bajada.
"""

# ╔═╡ 16c83e10-ca79-11f1-92c6-893211b466d3
md"""
Los días grandes vienen **en racimos**: después de un día turbulento suele
venir otro día turbulento. Medimos esto con la autocorrelación de ``|r|``.

Rezagos: $(@bind L2 Slider(10:10:250, default = 100, show_value = true))
"""

# ╔═╡ 16c83e1a-ca79-11f1-8f8a-cbf396c2c11c
md"""
La autocorrelación de ``r`` (azul) muere de inmediato. La de ``|r|`` (naranja)
es positiva y **decae muy lentamente**: sigue siendo significativa después de
semanas o meses. La dirección del mercado no tiene memoria; **la volatilidad
sí**. Esto es lo que modelan los procesos GARCH (sesión 5).

**[Avanzado]** Grafica la autocorrelación de ``|r|`` en escala log-log. ¿Decae
como una exponencial ``e^{-\tau/\tau_0}`` o como una ley de potencias
``\tau^{-\beta}``? Estima ``\beta``.
"""

# ╔═╡ 16c83e1a-ca79-11f1-8f20-fff0365acc63
md"""
## 4. Hecho 4: normalidad agregacional

Los retornos de ``\Delta t`` días son sumas de ``\Delta t`` retornos diarios.
Si el teorema del límite central aplica, al agregar la distribución debería
acercarse a la normal. ¿Pasa eso?

Horizonte ``\Delta t`` (días): $(@bind Δt Slider([1, 2, 5, 10, 20, 60], default = 1, show_value = true))
"""

# ╔═╡ 16c83e24-ca79-11f1-b872-afecedd39a49
md"""
Al agregar, la curtosis baja y la distribución se acerca a la normal: el
teorema del límite central termina ganando. Pero **gana más lento** de lo que
predice la línea punteada, que supone días independientes.

**[Avanzado]** Demuestra que si ``r_1, \dots, r_n`` son independientes e
idénticamente distribuidos con curtosis en exceso ``\kappa_1``, su suma tiene
curtosis ``\kappa_1/n``. ¿Qué hecho estilizado explica que los datos se alejen
de esa predicción?
"""

# ╔═╡ 16c83e30-ca79-11f1-a6bf-efb0e3bbdaea
md"""
## 5. ¿Qué tan especial es esto? Modelos nulos

Para saber si un patrón es real, lo comparamos con un **modelo nulo**: una
versión de los datos donde el patrón no debería existir.

- **Caminata gaussiana**: retornos normales independientes, con la misma
  ``\mu`` y ``\sigma`` que los datos. Es el modelo de libro de texto.
- **Datos barajados**: los *mismos* retornos reales, pero en orden aleatorio.
  Conserva exactamente la distribución y destruye cualquier estructura
  temporal.

Modelo: $(@bind modelo Select(["Datos reales", "Datos barajados", "Caminata gaussiana"]))
  Semilla: $(@bind semilla NumberField(1:1000, default = 1))
"""

# ╔═╡ 16c83e38-ca79-11f1-8360-8d5f3df50348
r_modelo = if modelo == "Datos reales"
    r
elseif modelo == "Datos barajados"
    shuffle(Xoshiro(semilla), r)
else
    μ .+ σ .* randn(Xoshiro(semilla), length(r))
end;

# ╔═╡ 16c83e42-ca79-11f1-8381-553049270d77
md"""
Compara los tres modelos y responde:
- ¿Cuál conserva las colas pesadas? ¿Cuál conserva el agrupamiento de
  volatilidad?
- Al barajar, la curtosis no cambia pero la autocorrelación de ``|r|``
  desaparece. ¿Qué nos dice eso sobre **dónde** vive cada hecho estilizado?

La función `panel_hechos` resume los cuatro hechos de las sesiones 1 y 2 en una sola figura. La
usaremos en la sesión 20 para evaluar mercados artificiales.
"""

# ╔═╡ 16c83e42-ca79-11f1-a2f4-75c1bb661f54
md"""
## 6. Universalidad

Hasta ahora miramos un activo a la vez. La tabla resume los mismos indicadores
para los cuatro.
"""

# ╔═╡ 16c83e4c-ca79-11f1-a703-b3dc7bfb058d
md"""
Los cuatro mercados son muy distintos: un índice de acciones mexicanas, una
divisa, el mercado más grande del mundo y una criptomoneda sin flujos de
efectivo. Sus volatilidades anuales van de 12% a 66%. Y sin embargo
**los cuatro muestran los mismos hechos estilizados**: curtosis muy por encima
de 0, autocorrelación de ``r`` casi nula, autocorrelación de ``|r|`` positiva
y persistente, y muchos más eventos extremos de los que predice la normal.

En física, cuando sistemas microscópicamente distintos comparten el mismo
comportamiento macroscópico, se habla de **universalidad**. Es una pista de
que el patrón no depende de los detalles de cada mercado, sino de algo más
general: la forma en que interactúan muchos agentes. Esa es la hipótesis de
trabajo de la econofísica.
"""

# ╔═╡ 16c83e56-ca79-11f1-90e2-e79230332530
md"""
## Resumen

| Hecho estilizado | Cómo se mide | Lo que vimos |
|---|---|---|
| Colas pesadas (sesión 1) | Densidad en escala log; curtosis | Eventos de 5σ o más son comunes |
| Sin autocorrelación lineal (sesión 1) | ACF de ``r`` | Indistinguible de cero después del rezago 1 |
| Agrupamiento de volatilidad | ACF de ``\lvert r\rvert`` | Positiva y persistente por meses |
| Normalidad agregacional | Curtosis vs. ``\Delta t`` | Baja al agregar, pero más lento que con días independientes |

| Modelo nulo | Conserva | Destruye |
|---|---|---|
| Datos barajados | La distribución: colas pesadas | El orden temporal: agrupamiento de volatilidad |
| Caminata gaussiana | Media y varianza | Colas pesadas y agrupamiento |
"""

# ╔═╡ 16c83e56-ca79-11f1-b6c3-7f8d3f68df7f
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `agregar_retornos(r, Δt)`, que sume los retornos en bloques
consecutivos de `Δt` días que no se traslapan (si sobran días al final, se
descartan). Por ejemplo, con `Δt = 2` y los retornos `[1, 2, 3, 4, 5]` debe
devolver `[3, 7]`.
"""

# ╔═╡ 16c83e62-ca79-11f1-b0c0-2dfe5e922e71
function agregar_retornos(r, Δt)
    missing  # reemplaza esta línea
end

# ╔═╡ 16c83e62-ca79-11f1-9015-95e0f65c909f
md"""
**2.** Escribe `prob_racha_mia(z, k)`: entre los días con ``|z_t| > k`` (que tengan
un día siguiente), la fracción en que también ``|z_{t+1}| > k``. Es el número
de la sección 0.
"""

# ╔═╡ 16c83e6a-ca79-11f1-a606-25011a70ffa0
function prob_racha_mia(z, k)
    missing  # reemplaza esta línea
end

# ╔═╡ 16c83e74-ca79-11f1-83a7-cfc2693f7eac
md"""
### Ruta básica: interpretación

3. Si los retornos no tienen autocorrelación, ¿por qué ``|r|`` sí la tiene?
   ¿Qué implica para alguien que administra el riesgo de un portafolio?
4. Con el selector de modelos nulos: ¿qué hecho estilizado sobrevive al
   barajado y cuál no? Explica por qué.
5. Al agregar los retornos, la curtosis baja, pero más lento que la línea
   punteada. ¿Qué hecho estilizado explica esa lentitud?
6. Bitcoin y el tipo de cambio USD/MXN son mercados muy distintos. ¿Qué tienen
   en común según la tabla de universalidad? ¿Por qué eso sugiere que el
   mecanismo no depende de los detalles de cada mercado?

### Ruta avanzada

7. Resuelve los dos ejercicios **[Avanzado]** de las secciones 3 y 4.
8. **Efecto apalancamiento.** Calcula ``\mathrm{corr}(r_t, |r_{t+\tau}|)`` para
   ``\tau = 1, \dots, 20`` en el S&P 500. ¿Es simétrica la relación entre
   caídas y volatilidad futura?

## Para la próxima sesión

- Lectura: Cont (2001), *Empirical properties of asset returns: stylized facts
  and statistical issues*, secciones 1 a 3. Complementaria: Mantegna y Stanley
  (1999), capítulos 5 y 7; Requião da Cunha (2021), sección 2.4.
- En la sesión 3 mediremos **qué tan pesadas** son las colas: ¿siguen una ley
  de potencias?
"""

# ╔═╡ 16c83e74-ca79-11f1-9b3c-85b432c0fe6c
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 16c83e88-ca79-11f1-9f2a-f9fe2441e3a7
begin
    const AZUL    = "#2a78d6"   # serie principal
    const NARANJA = "#eb6834"   # segunda serie
    const GRIS    = "#52514e"   # referencias (normal, bandas, eventos)
end;

# ╔═╡ 16c83de8-ca79-11f1-9e50-d3758e5a755a
let
    x = randn(Xoshiro(3), 2_000)
    y = x .^ 2
    scatter(x, y, color = AZUL, ms = 2.5, msw = 0, alpha = 0.5,
            xlabel = "X", ylabel = "Y = X²",
            title = "corr(X, Y) = $(round(cor(x, y), digits = 3)), pero Y depende por completo de X")
end

# ╔═╡ 16c83e10-ca79-11f1-b4e4-5d70d591f2bf
let
    N = length(r)
    banda = 1.96 / sqrt(N)
    τ = 1:L2
    p = plot(τ, fill(banda, L2), fillrange = fill(-banda, L2), color = GRIS,
             fillalpha = 0.15, lw = 0, label = "",
             xlabel = "Rezago τ (días)", ylabel = "ρ(τ)", legend = :topright,
             title = "$activo: autocorrelación de r y de |r|")
    plot!(p, τ, autocor(r, τ), color = AZUL, lw = 2, label = "r")
    plot!(p, τ, autocor(abs.(r), τ), color = NARANJA, lw = 2, label = "|r|")
    hline!(p, [0], color = :black, lw = 0.5, label = "")
end

# ╔═╡ 16c83e88-ca79-11f1-9e16-9d4154f213b5
const EVENTOS = [Date(2008, 10, 10), Date(2016, 11, 9), Date(2020, 3, 16)];

# ╔═╡ 16c83e88-ca79-11f1-aee8-1d2e78570fa9
"Dibuja líneas verticales en los eventos que caen dentro del rango de `fechas`."
function marcar_eventos!(p, fechas)
    dentro = filter(e -> first(fechas) <= e <= last(fechas), EVENTOS)
    vline!(p, dentro, color = GRIS, ls = :dot, lw = 1, label = "")
end

# ╔═╡ 16c83e94-ca79-11f1-8654-7be9bbac4a05
"Retornos agregados en bloques no traslapados de `Δt` días."
agregar(r, Δt) = [sum(@view r[i:i+Δt-1]) for i in 1:Δt:length(r)-Δt+1]

# ╔═╡ 16c83e24-ca79-11f1-a541-35dfdcbf305e
let
    hs = [1, 2, 3, 5, 10, 20, 40, 60]
    κ = [kurtosis(agregar(r, h)) for h in hs]
    p = plot(hs, κ[1] ./ hs, color = GRIS, ls = :dash, lw = 1.5,
             label = "Si los días fueran independientes: κ₁/Δt",
             xscale = :log10, xticks = (hs, string.(hs)), legend = :topright,
             xlabel = "Horizonte Δt (días)", ylabel = "Curtosis en exceso",
             title = "$activo: las colas se adelgazan al agregar")
    plot!(p, hs, κ, color = AZUL, lw = 2, marker = :circle, ms = 5,
          markerstrokecolor = :white, markerstrokewidth = 2, label = "Datos")
    hline!(p, [0], color = :black, lw = 0.5, label = "")
end

# ╔═╡ 16c83e94-ca79-11f1-913c-572c2020d06f
"Densidad empírica en bins de ancho `ancho`; descarta bins vacíos (para escala log)."
function densidad_empirica(z; ancho = 0.25)
    h = fit(Histogram, z, -16:ancho:16)
    centros = collect(h.edges[1][1:end-1]) .+ ancho / 2
    dens = h.weights ./ (length(z) * ancho)
    keep = dens .> 0
    centros[keep], dens[keep]
end

# ╔═╡ 16c83e9c-ca79-11f1-ac25-516a3c7aabb8
dias_por_año(activo) = activo == "Bitcoin" ? 365 : 252

# ╔═╡ 16c83e9c-ca79-11f1-b2fd-fb4394fd12b2
"Límite simétrico del eje z: al menos ±8 y siempre mayor que el dato más extremo."
rango_z(z) = max(8, ceil(Int, maximum(abs, z)) + 1)

# ╔═╡ 16c83e1a-ca79-11f1-bc84-5bd6cb4e08d3
let
    rΔ = agregar(r, Δt)
    zΔ = (rΔ .- mean(rΔ)) ./ std(rΔ)
    x, dens = densidad_empirica(zΔ)
    zmax = rango_z(zΔ)
    p = plot(-zmax:0.05:zmax, pdf.(Normal(), -zmax:0.05:zmax), color = GRIS, ls = :dash,
             lw = 1.5, label = "Normal estándar", yscale = :log10,
             ylims = (1e-4, 1), xlims = (-zmax, zmax), legend = :topright,
             xlabel = "z", ylabel = "Densidad",
             title = "Retornos de $Δt día(s): N = $(length(rΔ)), " *
                     "curtosis en exceso = $(round(kurtosis(rΔ), digits = 1))")
    scatter!(p, x, dens, color = AZUL, ms = 4, markerstrokewidth = 0,
             label = "Datos")
end

# ╔═╡ 16c83e9c-ca79-11f1-b627-3d7c65673134
"Marcas del eje x: un año de cada `paso`, etiquetado sólo con el año."
function ticks_años(fechas; paso = 4)
    años = year(first(fechas)):paso:year(last(fechas))
    marcas = Date.(años, 1, 1)
    (Dates.value.(marcas), string.(años))
end

# ╔═╡ 16c83e06-ca79-11f1-bdd6-b11c234e6334
let
    p = plot(fechas_r, 100abs.(r), color = AZUL, lw = 0.6,
             ylabel = "|r| (%)", title = "$activo: tamaño de los movimientos diarios",
             xticks = ticks_años(fechas))
    marcar_eventos!(p, fechas)
end

# ╔═╡ 16c83ea6-ca79-11f1-a7b9-879a06e80075
function grafica_acf(x, L; titulo = "", color = AZUL)
    banda = 1.96 / sqrt(length(x))
    p = plot(1:L, fill(banda, L), fillrange = fill(-banda, L), color = GRIS,
             fillalpha = 0.15, lw = 0, xlabel = "Rezago τ (días)",
             ylabel = "ρ(τ)", title = titulo, ylims = (-0.15, 0.4))
    bar!(p, 1:L, autocor(x, 1:L), color = color, lw = 0, bar_width = 0.6)
    hline!(p, [0], color = :black, lw = 0.5)
end

# ╔═╡ 16c83ea6-ca79-11f1-9ba9-792259d5cbcd
"Resumen visual de los cuatro hechos estilizados para una serie de retornos."
function panel_hechos(x; titulo = "", L = 100)
    z = (x .- mean(x)) ./ std(x)
    p1 = plot(100x, color = AZUL, lw = 0.5, ylabel = "r (%)", xlabel = "Día",
              title = "Retornos")
    xs, dens = densidad_empirica(z)
    zmax = rango_z(z)
    p2 = plot(-zmax:0.05:zmax, pdf.(Normal(), -zmax:0.05:zmax), color = GRIS, ls = :dash,
              lw = 1.5, yscale = :log10, ylims = (1e-5, 1), xlims = (-zmax, zmax),
              xlabel = "z", title = "Densidad (κ = $(round(kurtosis(x), digits = 1)))")
    scatter!(p2, xs, dens, color = AZUL, ms = 3, markerstrokewidth = 0)
    p3 = grafica_acf(x, 30; titulo = "Autocorrelación de r")
    banda = 1.96 / sqrt(length(x))
    acf_abs = autocor(abs.(x), 1:L)
    p4 = plot(1:L, fill(banda, L), fillrange = fill(-banda, L), color = GRIS,
              fillalpha = 0.15, lw = 0, ylims = (-0.1, max(0.4, maximum(acf_abs) + 0.05)),
              xlabel = "Rezago τ (días)", title = "Autocorrelación de |r|")
    plot!(p4, 1:L, acf_abs, color = NARANJA, lw = 2)
    hline!(p4, [0], color = :black, lw = 0.5)
    plot(p1, p2, p3, p4, layout = (2, 2), size = (720, 520),
         plot_title = titulo, plot_titlefontsize = 12)
end

# ╔═╡ 16c83e38-ca79-11f1-9414-65592db3aaec
panel_hechos(r_modelo; titulo = "$activo: $(lowercase(modelo))")

# ╔═╡ 16c83eb0-ca79-11f1-b4d1-51cd7ca9b10e
function tabla_universalidad(datos)
    grupos = collect(groupby(datos, :activo))
    nombres = [g.activo[1] for g in grupos]
    rs = [diff(log.(sort(g, :fecha).cierre)) for g in grupos]
    zs = [(r .- mean(r)) ./ std(r) for r in rs]
    DataFrame(
        "Activo" => nombres,
        "Días" => length.(rs),
        "Volatilidad anual (%)" => [round(100std(r) * sqrt(dias_por_año(n)), digits = 1) for (r, n) in zip(rs, nombres)],
        "Curtosis en exceso" => [round(kurtosis(r), digits = 1) for r in rs],
        "ACF de r (τ = 1)" => [round(autocor(r, [1])[1], digits = 3) for r in rs],
        "ACF de |r| (promedio τ = 1–20)" => [round(mean(autocor(abs.(r), 1:20)), digits = 3) for r in rs],
        "Días con |z| > 5" => [count(>(5), abs.(z)) for z in zs],
        "Esperados (normal)" => [round(length(r) * 2ccdf(Normal(), 5), sigdigits = 2) for r in rs],
    )
end

# ╔═╡ 16c83e42-ca79-11f1-a666-39d9c0871adb
tabla_universalidad(datos)

# ╔═╡ 16c83eb0-ca79-11f1-9fb8-65beb012f48e
"Fracción de días con |z| > k, entre los que siguen a un día con |z| > k."
function prob_racha(z, k)
    dias = findall(t -> abs(z[t]) > k, 1:length(z)-1)
    mean(abs.(z[dias .+ 1]) .> k)
end

# ╔═╡ 16c83dca-ca79-11f1-ae9c-b1aecfd5a9ae
racha_sp = let
    sp = sort(filter(:activo => ==("S&P 500"), datos), :fecha)
    r_sp = diff(log.(sp.cierre))
    z_sp = (r_sp .- mean(r_sp)) ./ std(r_sp)
    (fechas = sp.fecha, r = r_sp, base = round(100mean(abs.(z_sp) .> 2), digits = 1),
     racha = round(100prob_racha(z_sp, 2), digits = 1))
end;

# ╔═╡ 16c83dd4-ca79-11f1-9b9f-83dd4f593196
let s = racha_sp
    p = plot(s.fechas[2:end], 100abs.(s.r), color = AZUL, lw = 0.6, ylabel = "|r| (%)",
             title = "S&P 500: tamaño de los movimientos diarios", xticks = ticks_años(s.fechas))
    marcar_eventos!(p, s.fechas)
end

# ╔═╡ 16c83dde-ca79-11f1-8840-61808a989611
md"""
En un día cualquiera, la probabilidad de un movimiento de más de 2
desviaciones estándar es de $(racha_sp.base)%. Pero si **ayer** hubo uno, la
probabilidad de que **hoy** haya otro sube a $(racha_sp.racha)%. Los días
turbulentos vienen en racimos. ¿Cómo encaja eso con un mercado sin memoria?
"""

# ╔═╡ 16c83eba-ca79-11f1-9a3a-8942ea96b22d
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

# ╔═╡ 16c83e62-ca79-11f1-b89f-2393c12fb1ed
verificar(() -> agregar_retornos([1.0, 2.0, 3.0, 4.0, 5.0], 2), [3.0, 7.0];
          pista = md"Recorre los inicios de bloque `1:Δt:length(r)-Δt+1` y suma `r[i:i+Δt-1]`.")

# ╔═╡ 16c83e6a-ca79-11f1-9bef-273961ceab64
verificar(() -> prob_racha_mia([3.0, 3.0, 0.0, 3.0, 0.0, 0.0, 3.0, 3.0], 2), 0.5;
          pista = md"Encuentra los días `t` con `abs(z[t]) > k` entre `1` y `length(z) - 1`, y promedia `abs(z[t + 1]) > k` sobre ellos.")

# ╔═╡ Cell order:
# ╟─16c83d70-ca79-11f1-b8ac-d135e84b438d
# ╟─16c83dbe-ca79-11f1-8ba9-65c659edc421
# ╟─16c83dca-ca79-11f1-b855-75e62f6ac58a
# ╟─16c83dca-ca79-11f1-ae9c-b1aecfd5a9ae
# ╟─16c83dd4-ca79-11f1-9b9f-83dd4f593196
# ╟─16c83dde-ca79-11f1-8840-61808a989611
# ╟─16c83de8-ca79-11f1-bbf5-0f12ea464605
# ╟─16c83de8-ca79-11f1-9e50-d3758e5a755a
# ╟─16c83df2-ca79-11f1-aa42-4176e97203ad
# ╠═16c83dfc-ca79-11f1-a3e7-11b82d4a537e
# ╟─16c83dfc-ca79-11f1-8f87-4392ed895cbd
# ╠═16c83dfc-ca79-11f1-9681-eb413b011b73
# ╟─16c83e06-ca79-11f1-86bf-0dd372c0a900
# ╟─16c83e06-ca79-11f1-bdd6-b11c234e6334
# ╟─16c83e10-ca79-11f1-92c6-893211b466d3
# ╟─16c83e10-ca79-11f1-b4e4-5d70d591f2bf
# ╟─16c83e1a-ca79-11f1-8f8a-cbf396c2c11c
# ╟─16c83e1a-ca79-11f1-8f20-fff0365acc63
# ╟─16c83e1a-ca79-11f1-bc84-5bd6cb4e08d3
# ╟─16c83e24-ca79-11f1-a541-35dfdcbf305e
# ╟─16c83e24-ca79-11f1-b872-afecedd39a49
# ╟─16c83e30-ca79-11f1-a6bf-efb0e3bbdaea
# ╠═16c83e38-ca79-11f1-8360-8d5f3df50348
# ╟─16c83e38-ca79-11f1-9414-65592db3aaec
# ╟─16c83e42-ca79-11f1-8381-553049270d77
# ╟─16c83e42-ca79-11f1-a2f4-75c1bb661f54
# ╟─16c83e42-ca79-11f1-a666-39d9c0871adb
# ╟─16c83e4c-ca79-11f1-a703-b3dc7bfb058d
# ╟─16c83e56-ca79-11f1-90e2-e79230332530
# ╟─16c83e56-ca79-11f1-b6c3-7f8d3f68df7f
# ╠═16c83e62-ca79-11f1-b0c0-2dfe5e922e71
# ╟─16c83e62-ca79-11f1-b89f-2393c12fb1ed
# ╟─16c83e62-ca79-11f1-9015-95e0f65c909f
# ╠═16c83e6a-ca79-11f1-a606-25011a70ffa0
# ╟─16c83e6a-ca79-11f1-9bef-273961ceab64
# ╟─16c83e74-ca79-11f1-83a7-cfc2693f7eac
# ╟─16c83e74-ca79-11f1-9b3c-85b432c0fe6c
# ╠═16c83e7e-ca79-11f1-9ffb-1f0c5e61f024
# ╠═16c83e88-ca79-11f1-9f2a-f9fe2441e3a7
# ╠═16c83e88-ca79-11f1-9e16-9d4154f213b5
# ╠═16c83e88-ca79-11f1-aee8-1d2e78570fa9
# ╠═16c83e94-ca79-11f1-8654-7be9bbac4a05
# ╠═16c83e94-ca79-11f1-913c-572c2020d06f
# ╠═16c83e9c-ca79-11f1-ac25-516a3c7aabb8
# ╠═16c83e9c-ca79-11f1-b2fd-fb4394fd12b2
# ╠═16c83e9c-ca79-11f1-b627-3d7c65673134
# ╠═16c83ea6-ca79-11f1-a7b9-879a06e80075
# ╠═16c83ea6-ca79-11f1-9ba9-792259d5cbcd
# ╠═16c83eb0-ca79-11f1-b4d1-51cd7ca9b10e
# ╠═16c83eb0-ca79-11f1-9fb8-65beb012f48e
# ╠═16c83eba-ca79-11f1-9a3a-8942ea96b22d
