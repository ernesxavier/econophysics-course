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

# ╔═╡ 68a07cb6-ca92-11f1-956a-0f4d9e8391a4
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

# ╔═╡ 68a07be4-ca92-11f1-8c2d-812bf3cb5c0b
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

**Conceptos transversales.** Se completan la **media, la varianza y la
correlación**. Se retoman los **retornos logarítmicos** (sesión 1).
"""

# ╔═╡ 68a07c0c-ca92-11f1-a6a6-e30f91840b30
TableOfContents(title = "Contenido")

# ╔═╡ 68a07c16-ca92-11f1-8813-a1d1c36f6730
md"""
## 0. El fenómeno: los días malos llegan juntos

En la sesión 1 vimos que la **dirección** del mercado casi no tiene memoria:
saber si hoy subió no ayuda a predecir si mañana sube. Pero mira el **tamaño**
de los movimientos del S&P 500:
"""

# ╔═╡ 68a07c28-ca92-11f1-acd2-091aa1885efe
md"""
## 1. Repaso: incorrelación no es independencia

Dos variables son **independientes** si conocer una no dice nada de la otra.
Son **incorreladas** si no tienen relación *lineal*: ``\mathrm{corr}(X, Y) = 0``.
Independencia implica incorrelación, pero no al revés. El ejemplo clásico:
``X`` normal estándar e ``Y = X^2``.
"""

# ╔═╡ 68a07c34-ca92-11f1-9702-81fb28c7f8fc
md"""
La correlación es casi cero porque los valores negativos y positivos de ``X``
dan el mismo ``Y``: no hay una recta que describa la relación. Pero ``Y`` está
completamente determinada por ``X``. Con los retornos pasa algo parecido: el
signo de ``r_t`` no predice el de ``r_{t+1}``, pero el **tamaño** ``|r_t|`` sí
dice algo sobre ``|r_{t+1}|``.

## 2. Los datos

Los mismos cuatro activos de la sesión 1 (`Data/precios.csv`).
"""

# ╔═╡ 68a07c3e-ca92-11f1-8746-0ffac5cbd5a7
datos = CSV.read(joinpath(@__DIR__, "..", "Data", "precios.csv"), DataFrame)

# ╔═╡ 68a07c3e-ca92-11f1-b431-0746bae0d2d7
md"""
**Activo:** $(@bind activo Select(unique(datos.activo)))
"""

# ╔═╡ 68a07c3e-ca92-11f1-ab3b-31705a74f3a1
begin
    serie = sort(filter(:activo => ==(activo), datos), :fecha)
    fechas = serie.fecha
    precios = serie.cierre
    r = diff(log.(precios))      # retornos logarítmicos diarios
    fechas_r = fechas[2:end]
    μ, σ = mean(r), std(r)
end;

# ╔═╡ 68a07c48-ca92-11f1-8e0a-bd6a17a56da4
md"""
## 3. Hecho 3: agrupamiento de volatilidad

Ahora miramos el **tamaño** del movimiento, ``|r_t|``, sin importar si fue
subida o bajada.
"""

# ╔═╡ 68a07c52-ca92-11f1-a245-3bf92d57192b
md"""
Los días grandes vienen **en racimos**: después de un día turbulento suele
venir otro día turbulento. Medimos esto con la autocorrelación de ``|r|``.

Rezagos: $(@bind L2 Slider(10:10:250, default = 100, show_value = true))
"""

# ╔═╡ 68a07c5c-ca92-11f1-8adc-45faedfc9e8b
md"""
La autocorrelación de ``r`` (azul) muere de inmediato. La de ``|r|`` (naranja)
es positiva y **decae muy lentamente**: sigue siendo significativa después de
semanas o meses. La dirección del mercado no tiene memoria; **la volatilidad
sí**. Esto es lo que modelan los procesos GARCH (sesión 5).

**[Avanzado]** Grafica la autocorrelación de ``|r|`` en escala log-log. ¿Decae
como una exponencial ``e^{-\tau/\tau_0}`` o como una ley de potencias
``\tau^{-\beta}``? Estima ``\beta``.
"""

# ╔═╡ 68a07c5c-ca92-11f1-9776-c9e747877d64
md"""
## 4. Hecho 4: normalidad agregacional

Los retornos de ``\Delta t`` días son sumas de ``\Delta t`` retornos diarios.
Si el teorema del límite central aplica, al agregar la distribución debería
acercarse a la normal. ¿Pasa eso?

Horizonte ``\Delta t`` (días): $(@bind Δt Slider([1, 2, 5, 10, 20, 60], default = 1, show_value = true))
"""

# ╔═╡ 68a07c66-ca92-11f1-b07c-65311a17037c
md"""
Al agregar, la curtosis baja y la distribución se acerca a la normal: el
teorema del límite central termina ganando. Pero **gana más lento** de lo que
predice la línea punteada, que supone días independientes.

**[Avanzado]** Demuestra que si ``r_1, \dots, r_n`` son independientes e
idénticamente distribuidos con curtosis en exceso ``\kappa_1``, su suma tiene
curtosis ``\kappa_1/n``. ¿Qué hecho estilizado explica que los datos se alejen
de esa predicción?
"""

# ╔═╡ 68a07c70-ca92-11f1-b4e9-4993515b5710
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

# ╔═╡ 68a07c70-ca92-11f1-9d98-573ddf062da1
r_modelo = if modelo == "Datos reales"
    r
elseif modelo == "Datos barajados"
    shuffle(Xoshiro(semilla), r)
else
    μ .+ σ .* randn(Xoshiro(semilla), length(r))
end;

# ╔═╡ 68a07c7a-ca92-11f1-8f07-17b10caae732
md"""
Compara los tres modelos y responde:
- ¿Cuál conserva las colas pesadas? ¿Cuál conserva el agrupamiento de
  volatilidad?
- Al barajar, la curtosis no cambia pero la autocorrelación de ``|r|``
  desaparece. ¿Qué nos dice eso sobre **dónde** vive cada hecho estilizado?

La función `panel_hechos` resume los cuatro hechos de las sesiones 1 y 2 en una sola figura. La
usaremos en la sesión 20 para evaluar mercados artificiales.
"""

# ╔═╡ 68a07c7a-ca92-11f1-aff8-65359beb63a7
md"""
## 6. Universalidad

Hasta ahora miramos un activo a la vez. La tabla resume los mismos indicadores
para los cuatro.
"""

# ╔═╡ 68a07c84-ca92-11f1-8a71-b3d76f31d9b8
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

# ╔═╡ 68a07c8e-ca92-11f1-9ef0-3163b0fe8193
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

# ╔═╡ 68a07c8e-ca92-11f1-a036-f1354a46b762
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `agregar_retornos(r, Δt)`, que sume los retornos en bloques
consecutivos de `Δt` días que no se traslapan (si sobran días al final, se
descartan). Por ejemplo, con `Δt = 2` y los retornos `[1, 2, 3, 4, 5]` debe
devolver `[3, 7]`.
"""

# ╔═╡ 68a07c9a-ca92-11f1-b2be-c766ea145013
function agregar_retornos(r, Δt)
    missing  # reemplaza esta línea
end

# ╔═╡ 68a07c9a-ca92-11f1-a092-6b6c5dfb13dd
md"""
**2.** Escribe `prob_racha_mia(z, k)`: entre los días con ``|z_t| > k`` (que tengan
un día siguiente), la fracción en que también ``|z_{t+1}| > k``. Es el número
de la sección 0.
"""

# ╔═╡ 68a07ca2-ca92-11f1-8771-97e9daa7a217
function prob_racha_mia(z, k)
    missing  # reemplaza esta línea
end

# ╔═╡ 68a07cac-ca92-11f1-8cc9-3d6a4708ad92
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


### Conceptos que vuelven

Preguntas cortas para repasar ideas de sesiones anteriores; se responden con
lo que ya sabes y con el laboratorio de hoy.

- **R1.** En la sección 4 sumamos retornos diarios para obtener los de
  ``\Delta t`` días. ¿Por qué esa suma es correcta con retornos logarítmicos y
  no con retornos simples? (Sesión 1.)

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

# ╔═╡ 68a07cac-ca92-11f1-b16a-abc15ddd8336
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 68a07cb6-ca92-11f1-aa4b-edd560938362
begin
    const AZUL    = "#2a78d6"   # serie principal
    const NARANJA = "#eb6834"   # segunda serie
    const GRIS    = "#52514e"   # referencias (normal, bandas, eventos)
end;

# ╔═╡ 68a07c28-ca92-11f1-aecc-f750e0b31dec
let
    x = randn(Xoshiro(3), 2_000)
    y = x .^ 2
    scatter(x, y, color = AZUL, ms = 2.5, msw = 0, alpha = 0.5,
            xlabel = "X", ylabel = "Y = X²",
            title = "corr(X, Y) = $(round(cor(x, y), digits = 3)), pero Y depende por completo de X")
end

# ╔═╡ 68a07c52-ca92-11f1-93b1-7b6cb97b0391
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

# ╔═╡ 68a07cb6-ca92-11f1-a298-0fa856b0aa5b
const EVENTOS = [Date(2008, 10, 10), Date(2016, 11, 9), Date(2020, 3, 16)];

# ╔═╡ 68a07cc0-ca92-11f1-be02-b538452978d1
"Dibuja líneas verticales en los eventos que caen dentro del rango de `fechas`."
function marcar_eventos!(p, fechas)
    dentro = filter(e -> first(fechas) <= e <= last(fechas), EVENTOS)
    vline!(p, dentro, color = GRIS, ls = :dot, lw = 1, label = "")
end

# ╔═╡ 68a07cc0-ca92-11f1-a0cb-e7ed3fde1f8b
"Retornos agregados en bloques no traslapados de `Δt` días."
agregar(r, Δt) = [sum(@view r[i:i+Δt-1]) for i in 1:Δt:length(r)-Δt+1]

# ╔═╡ 68a07c66-ca92-11f1-9a6e-2bb6a7c963e8
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

# ╔═╡ 68a07ccc-ca92-11f1-b35f-61f8f17deb4d
"Densidad empírica en bins de ancho `ancho`; descarta bins vacíos (para escala log)."
function densidad_empirica(z; ancho = 0.25)
    h = fit(Histogram, z, -16:ancho:16)
    centros = collect(h.edges[1][1:end-1]) .+ ancho / 2
    dens = h.weights ./ (length(z) * ancho)
    keep = dens .> 0
    centros[keep], dens[keep]
end

# ╔═╡ 68a07ccc-ca92-11f1-9e8d-49d632a4b558
dias_por_año(activo) = activo == "Bitcoin" ? 365 : 252

# ╔═╡ 68a07cd4-ca92-11f1-8b6b-33959cc9e643
"Límite simétrico del eje z: al menos ±8 y siempre mayor que el dato más extremo."
rango_z(z) = max(8, ceil(Int, maximum(abs, z)) + 1)

# ╔═╡ 68a07c5c-ca92-11f1-bd86-990ba444e080
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

# ╔═╡ 68a07cd4-ca92-11f1-9a69-53fc7f26cac2
"Marcas del eje x: un año de cada `paso`, etiquetado sólo con el año."
function ticks_años(fechas; paso = 4)
    años = year(first(fechas)):paso:year(last(fechas))
    marcas = Date.(años, 1, 1)
    (Dates.value.(marcas), string.(años))
end

# ╔═╡ 68a07c48-ca92-11f1-bc8c-59fc17978061
let
    p = plot(fechas_r, 100abs.(r), color = AZUL, lw = 0.6,
             ylabel = "|r| (%)", title = "$activo: tamaño de los movimientos diarios",
             xticks = ticks_años(fechas))
    marcar_eventos!(p, fechas)
end

# ╔═╡ 68a07cd4-ca92-11f1-b240-f1d2bd627786
function grafica_acf(x, L; titulo = "", color = AZUL)
    banda = 1.96 / sqrt(length(x))
    p = plot(1:L, fill(banda, L), fillrange = fill(-banda, L), color = GRIS,
             fillalpha = 0.15, lw = 0, xlabel = "Rezago τ (días)",
             ylabel = "ρ(τ)", title = titulo, ylims = (-0.15, 0.4))
    bar!(p, 1:L, autocor(x, 1:L), color = color, lw = 0, bar_width = 0.6)
    hline!(p, [0], color = :black, lw = 0.5)
end

# ╔═╡ 68a07cde-ca92-11f1-94e0-f92a7e045cb7
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

# ╔═╡ 68a07c7a-ca92-11f1-95e6-35efdc2ff389
panel_hechos(r_modelo; titulo = "$activo: $(lowercase(modelo))")

# ╔═╡ 68a07cde-ca92-11f1-92a9-b5205136ebca
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

# ╔═╡ 68a07c84-ca92-11f1-a84c-cf9db2a6f274
tabla_universalidad(datos)

# ╔═╡ 68a07ce8-ca92-11f1-a5b6-075ef1131d2a
"Fracción de días con |z| > k, entre los que siguen a un día con |z| > k."
function prob_racha(z, k)
    dias = findall(t -> abs(z[t]) > k, 1:length(z)-1)
    mean(abs.(z[dias .+ 1]) .> k)
end

# ╔═╡ 68a07c16-ca92-11f1-8f01-437ba7403d9b
racha_sp = let
    sp = sort(filter(:activo => ==("S&P 500"), datos), :fecha)
    r_sp = diff(log.(sp.cierre))
    z_sp = (r_sp .- mean(r_sp)) ./ std(r_sp)
    (fechas = sp.fecha, r = r_sp, base = round(100mean(abs.(z_sp) .> 2), digits = 1),
     racha = round(100prob_racha(z_sp, 2), digits = 1))
end;

# ╔═╡ 68a07c20-ca92-11f1-8c8a-21063a6297ea
let s = racha_sp
    p = plot(s.fechas[2:end], 100abs.(s.r), color = AZUL, lw = 0.6, ylabel = "|r| (%)",
             title = "S&P 500: tamaño de los movimientos diarios", xticks = ticks_años(s.fechas))
    marcar_eventos!(p, s.fechas)
end

# ╔═╡ 68a07c28-ca92-11f1-bafd-67406e63ca6b
md"""
En un día cualquiera, la probabilidad de un movimiento de más de 2
desviaciones estándar es de $(racha_sp.base)%. Pero si **ayer** hubo uno, la
probabilidad de que **hoy** haya otro sube a $(racha_sp.racha)%. Los días
turbulentos vienen en racimos. ¿Cómo encaja eso con un mercado sin memoria?
"""

# ╔═╡ 68a07ce8-ca92-11f1-9c60-7dd449f0a9ca
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

# ╔═╡ 68a07c9a-ca92-11f1-9703-5f73d554de69
verificar(() -> agregar_retornos([1.0, 2.0, 3.0, 4.0, 5.0], 2), [3.0, 7.0];
          pista = md"Recorre los inicios de bloque `1:Δt:length(r)-Δt+1` y suma `r[i:i+Δt-1]`.")

# ╔═╡ 68a07ca2-ca92-11f1-9ee8-87125568694c
verificar(() -> prob_racha_mia([3.0, 3.0, 0.0, 3.0, 0.0, 0.0, 3.0, 3.0], 2), 0.5;
          pista = md"Encuentra los días `t` con `abs(z[t]) > k` entre `1` y `length(z) - 1`, y promedia `abs(z[t + 1]) > k` sobre ellos.")

# ╔═╡ Cell order:
# ╟─68a07be4-ca92-11f1-8c2d-812bf3cb5c0b
# ╟─68a07c0c-ca92-11f1-a6a6-e30f91840b30
# ╟─68a07c16-ca92-11f1-8813-a1d1c36f6730
# ╟─68a07c16-ca92-11f1-8f01-437ba7403d9b
# ╟─68a07c20-ca92-11f1-8c8a-21063a6297ea
# ╟─68a07c28-ca92-11f1-bafd-67406e63ca6b
# ╟─68a07c28-ca92-11f1-acd2-091aa1885efe
# ╟─68a07c28-ca92-11f1-aecc-f750e0b31dec
# ╟─68a07c34-ca92-11f1-9702-81fb28c7f8fc
# ╠═68a07c3e-ca92-11f1-8746-0ffac5cbd5a7
# ╟─68a07c3e-ca92-11f1-b431-0746bae0d2d7
# ╠═68a07c3e-ca92-11f1-ab3b-31705a74f3a1
# ╟─68a07c48-ca92-11f1-8e0a-bd6a17a56da4
# ╟─68a07c48-ca92-11f1-bc8c-59fc17978061
# ╟─68a07c52-ca92-11f1-a245-3bf92d57192b
# ╟─68a07c52-ca92-11f1-93b1-7b6cb97b0391
# ╟─68a07c5c-ca92-11f1-8adc-45faedfc9e8b
# ╟─68a07c5c-ca92-11f1-9776-c9e747877d64
# ╟─68a07c5c-ca92-11f1-bd86-990ba444e080
# ╟─68a07c66-ca92-11f1-9a6e-2bb6a7c963e8
# ╟─68a07c66-ca92-11f1-b07c-65311a17037c
# ╟─68a07c70-ca92-11f1-b4e9-4993515b5710
# ╠═68a07c70-ca92-11f1-9d98-573ddf062da1
# ╟─68a07c7a-ca92-11f1-95e6-35efdc2ff389
# ╟─68a07c7a-ca92-11f1-8f07-17b10caae732
# ╟─68a07c7a-ca92-11f1-aff8-65359beb63a7
# ╟─68a07c84-ca92-11f1-a84c-cf9db2a6f274
# ╟─68a07c84-ca92-11f1-8a71-b3d76f31d9b8
# ╟─68a07c8e-ca92-11f1-9ef0-3163b0fe8193
# ╟─68a07c8e-ca92-11f1-a036-f1354a46b762
# ╠═68a07c9a-ca92-11f1-b2be-c766ea145013
# ╟─68a07c9a-ca92-11f1-9703-5f73d554de69
# ╟─68a07c9a-ca92-11f1-a092-6b6c5dfb13dd
# ╠═68a07ca2-ca92-11f1-8771-97e9daa7a217
# ╟─68a07ca2-ca92-11f1-9ee8-87125568694c
# ╟─68a07cac-ca92-11f1-8cc9-3d6a4708ad92
# ╟─68a07cac-ca92-11f1-b16a-abc15ddd8336
# ╠═68a07cb6-ca92-11f1-956a-0f4d9e8391a4
# ╠═68a07cb6-ca92-11f1-aa4b-edd560938362
# ╠═68a07cb6-ca92-11f1-a298-0fa856b0aa5b
# ╠═68a07cc0-ca92-11f1-be02-b538452978d1
# ╠═68a07cc0-ca92-11f1-a0cb-e7ed3fde1f8b
# ╠═68a07ccc-ca92-11f1-b35f-61f8f17deb4d
# ╠═68a07ccc-ca92-11f1-9e8d-49d632a4b558
# ╠═68a07cd4-ca92-11f1-8b6b-33959cc9e643
# ╠═68a07cd4-ca92-11f1-9a69-53fc7f26cac2
# ╠═68a07cd4-ca92-11f1-b240-f1d2bd627786
# ╠═68a07cde-ca92-11f1-94e0-f92a7e045cb7
# ╠═68a07cde-ca92-11f1-92a9-b5205136ebca
# ╠═68a07ce8-ca92-11f1-a5b6-075ef1131d2a
# ╠═68a07ce8-ca92-11f1-9c60-7dd449f0a9ca
