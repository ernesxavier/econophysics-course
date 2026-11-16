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

# ╔═╡ e7853d76-ca94-11f1-806f-978633577a7a
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

# ╔═╡ e7853c22-ca94-11f1-b634-0d9eaca1d1cf
md"""
# Taller del Módulo I: los hechos estilizados de una acción

**Qué se evalúa.** Que sepas aplicar, sin ayuda, el diagnóstico de las
sesiones 1 a 7 a un activo que no vimos en clase, y que sepas **interpretar**
lo que sale. No se piden derivaciones.

**Tu acción.** El profesor te asigna una de las 95 acciones del S&P 100
(`Data/sp100_precios.csv`, de enero de 2015 a septiembre de 2026). Elígela en
el selector de la sección 1 y no la cambies: todas tus respuestas deben
referirse a ella.

**Cómo entregar.**
1. Completa las funciones de los ejercicios de código (reemplaza `missing`).
   La caja verde de `verificar` te dice si están bien.
2. Responde cada pregunta de interpretación en la celda **Tu respuesta**
   que está debajo, en tres a seis líneas. Cita los números de tus gráficas.
3. Entrega el archivo `.jl` del cuaderno y una exportación a HTML o PDF (menú
   de Pluto, "Export"). El trabajo es individual: puedes discutir ideas con tus
   compañeros, pero el código y la redacción deben ser tuyos.

**Cómo se califica** (100 puntos, 10% de la calificación final):

| Parte | Código | Interpretación | Total |
|---|---|---|---|
| 1. Colas y extremos | 8 | 14 | 22 |
| 2. Memoria y volatilidad | 8 | 14 | 22 |
| 3. Qué tan pesada es la cola | 7 | 14 | 21 |
| 4. ¿Qué modelo la describe? | — | 14 | 14 |
| 5. La acción y el mercado | 7 | 14 | 21 |

Una respuesta completa dice **qué** se observa, **con qué número o
gráfica** y **por qué** importa. La ruta avanzada, al final, es opcional y
cuenta como extra.

**En el semestre** hay una sesión de taller en clase: en los primeros 20
minutos se resuelven dudas del módulo, después se trabaja en parejas en las
partes 1 y 2, y al final cada pareja comenta una respuesta con el grupo. En el
trimestre el taller se resuelve en casa.

**Repaso recomendado** (en `Class_notes/repaso/`): P.B2–B4 (media, varianza,
histogramas y la normal), C.B2 (escalas logarítmicas), S.B5 (permutaciones).
Para la ruta avanzada: S.I2 (Ljung-Box) y S.I3 (estimación de colas).
"""

# ╔═╡ e7853c54-ca94-11f1-81df-17461ae29fb5
TableOfContents(title = "Contenido")

# ╔═╡ e7853c5e-ca94-11f1-b0af-f3a5ab2880e0
md"""
## 1. Datos y retornos

Cargamos las 95 acciones y el índice S&P 500 de `Data/precios.csv`, que
usaremos en la parte 5 como referencia del mercado.
"""

# ╔═╡ e7853c68-ca94-11f1-a3d9-a7c5f6965898
begin
    precios = CSV.read(joinpath(@__DIR__, "..", "Data", "sp100_precios.csv"), DataFrame)
    sectores = CSV.read(joinpath(@__DIR__, "..", "Data", "sp100_sectores.csv"), DataFrame)
    indices = CSV.read(joinpath(@__DIR__, "..", "Data", "precios.csv"), DataFrame)
    opciones = [f.ticker => "$(f.ticker): $(f.nombre) ($(f.sector))" for f in eachrow(sectores)]
end;

# ╔═╡ e7853c68-ca94-11f1-a074-e71355470b65
md"""
**Tu acción:** $(@bind ticker Select(opciones))
"""

# ╔═╡ e7853c72-ca94-11f1-b874-b1e9836f9b9d
begin
    serie = sort(filter(:ticker => ==(ticker), precios), :fecha)
    fechas = serie.fecha[2:end]
    r = diff(log.(serie.cierre))       # retornos logarítmicos diarios
    z = (r .- mean(r)) ./ std(r)       # retornos estandarizados
end;

# ╔═╡ e7853c7a-ca94-11f1-a4e7-1bc765bf436c
Markdown.parse("""
**Tu muestra.** $(length(r)) días de operación. Retorno medio diario
$(round(100mean(r), digits = 3))%, volatilidad diaria $(round(100std(r), digits = 2))%
(anualizada: $(round(100std(r) * sqrt(252), digits = 1))%).
""")

# ╔═╡ e7853c86-ca94-11f1-a05a-a33cf3adcd06
md"""
## 2. Parte 1: colas y extremos (hechos 1 y 2)

### Ejercicio 1 (4 puntos)

Escribe `exceso_extremos(z, k)`: el número de días con ``|z_t| > k``
dividido entre el número que esperarías si los retornos fueran normales,
``n \cdot \mathbb P(|Z| > k) = 2n\,\bar\Phi(k)``. En Julia,
``\bar\Phi(k)`` es `ccdf(Normal(), k)`.
"""

# ╔═╡ e7853c86-ca94-11f1-9c6a-0b68fe6903c2
function exceso_extremos(z, k)
    missing  # reemplaza esta línea
end

# ╔═╡ e7853c9a-ca94-11f1-a0ab-2d02c770eff5
Markdown.parse("""
**Resumen.** Curtosis en exceso: **$(round(kurtosis(r), digits = 1))** (una
normal tiene 0). Días con ``|z| > 4``: $(count(x -> abs(x) > 4, z)); una normal
predice $(round(length(z) * 2ccdf(Normal(), 4), sigdigits = 2)). El peor día fue
el $(fechas[argmin(r)]), con ``r = `` $(round(100minimum(r), digits = 1))%
(``z = `` $(round(minimum(z), digits = 1))).
""")

# ╔═╡ e7853c9a-ca94-11f1-965f-3d14ca109da5
md"""
### Ejercicio 2 (4 puntos)

Escribe `autocorrelacion(x, τ)`, la autocorrelación muestral con rezago
``\tau``:

```math
\hat\rho(\tau) = \frac{\sum_{t=1}^{n-\tau} (x_t - \bar x)(x_{t+\tau} - \bar x)}{\sum_{t=1}^{n} (x_t - \bar x)^2} .
```
"""

# ╔═╡ e7853ca4-ca94-11f1-a5d2-bf1acdd332ec
function autocorrelacion(x, τ)
    missing  # reemplaza esta línea
end

# ╔═╡ e7853cac-ca94-11f1-bc46-8394d0b7dac0
md"""
**Pregunta 1** (7 puntos). ¿Las colas de tu acción son más pesadas que las
de una normal? Usa la curtosis y la gráfica de observados entre esperados.
¿Qué tan seguido habría ocurrido tu peor día si los retornos fueran normales?
"""

# ╔═╡ e7853cb8-ca94-11f1-ac84-c1cd2da7790e
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ e7853cb8-ca94-11f1-be9b-15a60fc9c5c3
md"""
**Pregunta 2** (7 puntos). ¿Los retornos de tu acción tienen autocorrelación?
¿Cuántas barras salen de la banda gris, y cuántas esperarías que salieran por
puro azar si la autocorrelación verdadera fuera cero? ¿Qué te dice eso sobre
la posibilidad de predecir el retorno de mañana con el de hoy?
"""

# ╔═╡ e7853cb8-ca94-11f1-9514-3776181b7b9c
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ e7853cc2-ca94-11f1-9cbb-a71cbaa66298
md"""
## 3. Parte 2: memoria de la volatilidad y agregación (hechos 3 y 4)

### Ejercicio 3 (8 puntos)

Escribe `curtosis_agregada(r, Δts)`: para cada horizonte ``\Delta t`` del
vector `Δts`, suma los retornos en bloques de ``\Delta t`` días que no se
traslapan (los días que sobran al final se descartan) y calcula la curtosis en
exceso de esos retornos agregados con `kurtosis`. Devuelve un vector con un
número por horizonte.
"""

# ╔═╡ e7853cc2-ca94-11f1-af1c-a1e0df7b806e
function curtosis_agregada(r, Δts)
    missing  # reemplaza esta línea
end

# ╔═╡ e7853ccc-ca94-11f1-9336-2349d764710c
md"""
**Modelo nulo:** $(@bind nulo_t1 Select(["Datos originales", "Datos barajados"]))

Barajar los días conserva la distribución de los retornos, pero destruye
cualquier dependencia en el tiempo (sesión 2).
"""

# ╔═╡ e7853cd6-ca94-11f1-839c-af4f97ddc658
md"""
**Pregunta 3** (7 puntos). Compara la autocorrelación de ``r`` con la de
``|r|``. ¿Qué hecho estilizado ves y cuántos días dura? ¿Qué pasa con las dos
curvas al barajar los datos, y qué te dice eso sobre **dónde** vive el
agrupamiento de volatilidad: en la distribución o en el orden de los días?
"""

# ╔═╡ e7853cd6-ca94-11f1-a19a-3781958cf69c
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ e7853ce0-ca94-11f1-8ea9-07f06bb7482b
md"""
**Pregunta 4** (7 puntos). ¿Qué le pasa a la curtosis al agregar los
retornos a 20 o 40 días? ¿Baja más rápido o más lento que la línea punteada?
Con los datos barajados, ¿cambia la respuesta? Explica la diferencia usando la
pregunta 3.
"""

# ╔═╡ e7853ce0-ca94-11f1-9bad-454a648b73bd
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ e7853cea-ca94-11f1-a735-a36ecf99948a
md"""
## 4. Parte 3: ¿qué tan pesada es la cola?

Separamos las **pérdidas** (``-z`` en los días con ``z < 0``) de las
**ganancias** (``z`` en los días con ``z > 0``) y dibujamos la CCDF de cada
una en escala log-log (sesión 3).

### Ejercicio 4 (7 puntos)

Escribe `hill_alpha(x, k)`, el estimador de Hill del exponente de la cola con
los ``k`` datos más grandes de `x`:

```math
\hat\alpha = \frac{k}{\sum_{i=1}^{k} \ln\big(x_{(i)}/x_{(k+1)}\big)} ,
```

donde ``x_{(1)} \ge x_{(2)} \ge \dots`` son los datos ordenados de mayor a
menor.
"""

# ╔═╡ e7853cea-ca94-11f1-84f6-1103f9e01b47
function hill_alpha(x, k)
    missing  # reemplaza esta línea
end

# ╔═╡ e7853cf4-ca94-11f1-9aaf-419c1e766391
md"""
Número de datos de la cola ``k``: $(@bind k_hill Slider(20:10:400, default = 120, show_value = true))
"""

# ╔═╡ e7853cfe-ca94-11f1-8d4b-6302417c207d
md"""
**Pregunta 5** (7 puntos). Mueve ``k`` y elige un valor razonable usando la
curva de Hill. ¿Qué criterio usaste? Reporta ``\hat\alpha`` de las pérdidas
con su error estándar. ¿Cambia mucho la estimación si duplicas ``k``? ¿Por qué
no conviene usar ni un ``k`` muy pequeño ni uno muy grande?
"""

# ╔═╡ e7853d08-ca94-11f1-a4ff-1dff5d2821ca
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ e7853d08-ca94-11f1-96b1-db4fc47b8fa5
md"""
**Pregunta 6** (7 puntos). Con tu ``\hat\alpha``, ¿la varianza de los
retornos es finita? ¿Y el cuarto momento, que es el que mide la curtosis?
(Una ley de potencias con exponente ``\alpha`` tiene finitos sólo los momentos
de orden menor que ``\alpha``, sesión 3.) ¿Qué implica eso para la confianza
que puedes tener en la curtosis de la parte 1? ¿Son las pérdidas más pesadas
que las ganancias?
"""

# ╔═╡ e7853d12-ca94-11f1-ad1f-ffeaf7f917d4
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ e7853d12-ca94-11f1-b0e1-4591a8d673ad
md"""
## 5. Parte 4: ¿qué modelo describe a tu acción?

Comparamos tus datos con dos modelos de la sesión 5, simulados con el mismo
número de días:

- **Caminata aleatoria gaussiana:** retornos normales e independientes con la
  misma volatilidad.
- **GARCH(1,1)** ajustado a tu acción por máxima verosimilitud: la varianza de
  mañana depende del retorno y de la varianza de hoy.
"""

# ╔═╡ e7853d26-ca94-11f1-acc5-3d38bc1b80c6
md"""
**Pregunta 7** (7 puntos). Para cada hecho estilizado (colas pesadas,
ausencia de autocorrelación de ``r``, agrupamiento de volatilidad), di si lo
reproduce la caminata gaussiana y si lo reproduce el GARCH. ¿Hay alguno que el
GARCH reproduzca sólo en parte? Usa los números de la tabla.
"""

# ╔═╡ e7853d30-ca94-11f1-8a1b-637066f729ed
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ e7853d30-ca94-11f1-92e6-454bccea1085
md"""
**Pregunta 8** (7 puntos). La fila del GARCH muestra un rango sobre 10
simulaciones. ¿Por qué no basta con una sola simulación para comparar el
modelo con los datos (sesión 4)? ¿El valor de tu acción cae dentro del rango
del GARCH en todas las columnas? Si en alguna no cae, ¿qué te dice eso del
modelo?
"""

# ╔═╡ e7853d3a-ca94-11f1-92e0-095739837d0d
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ e7853d3a-ca94-11f1-8ff0-5110b44fb1c4
md"""
## 6. Parte 5: la acción y el mercado

Ponemos lado a lado los retornos de tu acción y los del índice S&P 500, en
los días en que hay datos de los dos (sesión 7).

### Ejercicio 5 (7 puntos)

Escribe `beta_mercado(ra, rm)`, la pendiente de la recta de mínimos cuadrados
de los retornos de la acción `ra` contra los del mercado `rm`:
``\hat\beta = \widehat{\mathrm{cov}}(r_a, r_m)/\widehat{\mathrm{var}}(r_m)``.
"""

# ╔═╡ e7853d3a-ca94-11f1-bb8f-f7e1348d9a84
function beta_mercado(ra, rm)
    missing  # reemplaza esta línea
end

# ╔═╡ e7853d44-ca94-11f1-b635-b7a0b913c159
conjunto = let
    sp = sort(filter(:activo => ==("S&P 500"), indices), :fecha)
    a = DataFrame(fecha = serie.fecha, pa = serie.cierre)
    m = DataFrame(fecha = sp.fecha, pm = sp.cierre)
    j = sort(innerjoin(a, m; on = :fecha), :fecha)
    (fecha = j.fecha[2:end], ra = diff(log.(j.pa)), rm = diff(log.(j.pm)))
end;

# ╔═╡ e7853d50-ca94-11f1-9599-a54504cc3b04
mercado_t1 = let c = conjunto
    zm = (c.rm .- mean(c.rm)) ./ std(c.rm)
    agitado = abs.(zm) .> 2
    (β = cov(c.ra, c.rm) / var(c.rm), agitado = agitado, ρ = cor(c.ra, c.rm),
     ρ_calma = cor(c.ra[.!agitado], c.rm[.!agitado]), ρ_agitado = cor(c.ra[agitado], c.rm[agitado]))
end;

# ╔═╡ e7853d58-ca94-11f1-9027-ffd10be68603
Markdown.parse("""
Correlación en los días normales del índice: **$(round(mercado_t1.ρ_calma, digits = 2))**;
en los días agitados (``|z| > 2``, $(count(mercado_t1.agitado)) días):
**$(round(mercado_t1.ρ_agitado, digits = 2))**.
""")

# ╔═╡ e7853d58-ca94-11f1-9785-3b2afac4e1a1
md"""
**Pregunta 9** (7 puntos). Interpreta ``\hat\beta`` y la correlación de tu
acción con el mercado. Si el S&P 500 cae 3% en un día, ¿cuánto esperas que
caiga tu acción? ¿Qué fracción de la varianza de tu acción explica el mercado
(el cuadrado de la correlación)? ¿Qué parte del riesgo se podría eliminar
diversificando?
"""

# ╔═╡ e7853d62-ca94-11f1-863e-a113bc5e17c8
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ e7853d62-ca94-11f1-b0e3-9504f93446e6
md"""
**Pregunta 10** (7 puntos). Compara la correlación en días normales con la de
días agitados. ¿Qué implica para un inversionista que cuenta con la
diversificación para protegerse justo cuando más la necesita? Relaciónalo con
el modo de mercado de la sesión 7.
"""

# ╔═╡ e7853d62-ca94-11f1-9e57-613ec124b35e
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ e7853d6c-ca94-11f1-a531-35c1611ffbd0
md"""
## 7. Ruta avanzada (extra, hasta 10 puntos)

Haz **dos** de los cuatro. Agrega las celdas que necesites al final del
cuaderno.

- **A1. El ``x_{\min}`` de Clauset.** Elige el inicio de la cola de las
  pérdidas con el método de Clauset, Shalizi y Newman (2009) de la sesión 3
  (la función `clauset` está en el apéndice). Compara su ``\hat\alpha`` con tu
  elección de la pregunta 5.
- **A2. Ljung-Box.** Aplica la prueba de Ljung-Box (repaso S.I2) a ``r`` y a
  ``|r|`` con 20 rezagos. ¿Confirma tus respuestas a las preguntas 2 y 3?
- **A3. Efecto apalancamiento.** Calcula
  ``\mathrm{corr}(r_t, |r_{t+\tau}|)`` para ``\tau = 1, \dots, 20``. ¿Las
  caídas anuncian más volatilidad que las subidas?
- **A4. ¿Correlación o volatilidad?** Simula dos series normales con
  correlación constante 0.5 y calcula su correlación en los días en que
  ``|z|`` de la segunda es mayor que 2. ¿Sale 0.5? ¿Qué cambia esto en tu
  respuesta a la pregunta 10? (Forbes y Rigobon, 2002.)
"""

# ╔═╡ e7853d6c-ca94-11f1-8fc9-db9ca03a005d
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para resolver el taller.
"""

# ╔═╡ e7853d76-ca94-11f1-821c-7b6b39e509ee
begin
    const AZUL    = "#2a78d6"   # serie principal
    const NARANJA = "#eb6834"   # segunda serie
    const GRIS    = "#52514e"   # referencias (normal, bandas)
end;

# ╔═╡ e7853d50-ca94-11f1-a662-3565a78e0d1e
let
    c, m = conjunto, mercado_t1
    a = m.agitado
    p = scatter(100c.rm[.!a], 100c.ra[.!a]; color = AZUL, ms = 2, msw = 0, alpha = 0.4,
                label = "Días normales del índice", xlabel = "retorno del S&P 500 (%)",
                ylabel = "retorno de $ticker (%)", legend = :topleft, size = (620, 420),
                title = "β = $(round(m.β, digits = 2)); correlación = $(round(m.ρ, digits = 2))")
    scatter!(p, 100c.rm[a], 100c.ra[a]; color = NARANJA, ms = 3, msw = 0, alpha = 0.7,
             label = "Días con |z| del índice > 2")
    xs = collect(extrema(100c.rm))
    plot!(p, xs, m.β .* xs; color = :black, lw = 1.5, label = "Recta con pendiente β")
end

# ╔═╡ e7853d76-ca94-11f1-a909-b12bcf3a01c6
nombre_de(t) = sectores.nombre[findfirst(==(t), sectores.ticker)]

# ╔═╡ e7853c7a-ca94-11f1-9e40-e74a0b7b4b5c
let
    nombre = nombre_de(ticker)
    p1 = plot(serie.fecha, serie.cierre; color = AZUL, lw = 1, yscale = :log10,
              ylabel = "precio (escala log)", title = "$nombre: precio de cierre")
    p2 = plot(fechas, 100r; color = AZUL, lw = 0.5, ylabel = "r (%)",
              title = "Retornos logarítmicos diarios")
    plot(p1, p2; layout = (2, 1), size = (720, 460))
end

# ╔═╡ e7853d82-ca94-11f1-a08e-1b70ab63c8c8
"Retornos agregados en bloques no traslapados de `Δt` días."
agregar(r, Δt) = [sum(@view r[i:i+Δt-1]) for i in 1:Δt:length(r)-Δt+1]

# ╔═╡ e7853ccc-ca94-11f1-93df-d94292462c3d
let
    x = nulo_t1 == "Datos originales" ? r : shuffle(Xoshiro(1), r)
    L = 60
    banda = 1.96 / sqrt(length(x))
    p1 = plot(1:L, fill(banda, L); fillrange = fill(-banda, L), color = GRIS, fillalpha = 0.15,
              lw = 0, label = "", xlabel = "rezago τ (días)", ylabel = "ρ(τ)",
              ylims = (-0.1, 0.45), legend = :topright,
              title = "Autocorrelación ($(lowercase(nulo_t1)))")
    plot!(p1, 1:L, autocor(x, 1:L); color = AZUL, lw = 1.5, label = "r")
    plot!(p1, 1:L, autocor(abs.(x), 1:L); color = NARANJA, lw = 2, label = "|r|")
    Δts = [1, 2, 5, 10, 20, 40]
    κ = [kurtosis(agregar(x, Δt)) for Δt in Δts]
    p2 = plot(Δts, κ; color = AZUL, lw = 2, marker = :circle, ms = 4, msw = 0, xscale = :log10,
              xticks = (Δts, string.(Δts)), label = ticker, xlabel = "horizonte Δt (días)",
              ylabel = "curtosis en exceso", legend = :topright, title = "Normalidad agregacional")
    plot!(p2, Δts, κ[1] ./ Δts; color = GRIS, ls = :dash, label = "Si los días fueran independientes")
    hline!(p2, [0]; color = :black, lw = 0.5, label = "")
    plot(p1, p2; layout = (1, 2), size = (780, 340))
end

# ╔═╡ e7853d82-ca94-11f1-b283-11dcd210eb1b
"Densidad empírica en bins de ancho `ancho`; descarta bins vacíos (para escala log)."
function densidad_empirica(z; ancho = 0.25)
    h = fit(Histogram, z, -30:ancho:30)
    centros = collect(h.edges[1][1:end-1]) .+ ancho / 2
    dens = h.weights ./ (length(z) * ancho)
    keep = dens .> 0
    centros[keep], dens[keep]
end

# ╔═╡ e7853d8a-ca94-11f1-978f-c944f038e06d
"Límite simétrico del eje z: al menos ±8 y siempre mayor que el dato más extremo."
rango_z(z) = max(8, ceil(Int, maximum(abs, z)) + 1)

# ╔═╡ e7853d8a-ca94-11f1-bb04-f9f9b9557c5c
function grafica_acf(x, L; titulo = "", color = AZUL)
    banda = 1.96 / sqrt(length(x))
    p = plot(1:L, fill(banda, L), fillrange = fill(-banda, L), color = GRIS,
             fillalpha = 0.15, lw = 0, xlabel = "Rezago τ (días)",
             ylabel = "ρ(τ)", title = titulo, ylims = (-0.2, 0.2))
    bar!(p, 1:L, autocor(x, 1:L), color = color, lw = 0, bar_width = 0.6)
    hline!(p, [0], color = :black, lw = 0.5)
end

# ╔═╡ e7853cac-ca94-11f1-916c-25e5a11e8e86
grafica_acf(r, 30; titulo = "Autocorrelación de los retornos de $ticker")

# ╔═╡ e7853d94-ca94-11f1-9b5e-e1a63d331131
"Días con |z| > k entre los esperados con una normal (versión de referencia)."
exceso_extremos_ref(z, k) = count(x -> abs(x) > k, z) / (2 * length(z) * ccdf(Normal(), k))

# ╔═╡ e7853c90-ca94-11f1-ba49-b58c15b9f313
let
    xs, dens = densidad_empirica(z)
    zmax = rango_z(z)
    p1 = plot(-zmax:0.05:zmax, pdf.(Normal(), -zmax:0.05:zmax); color = GRIS, ls = :dash, lw = 1.5,
              yscale = :log10, ylims = (1e-5, 1), xlims = (-zmax, zmax), label = "Normal",
              xlabel = "z", ylabel = "densidad (escala log)", legend = :topright,
              title = "Densidad de los retornos estandarizados")
    scatter!(p1, xs, dens; color = AZUL, ms = 3, msw = 0, label = ticker)
    ks = 2:0.5:5
    ex = [exceso_extremos_ref(z, k) for k in ks]
    p2 = bar(ks, ex; color = NARANJA, lw = 0, bar_width = 0.35, yscale = :log10,
             xlabel = "umbral k", ylabel = "observados / esperados",
             title = "Días con |z| > k frente a la normal")
    hline!(p2, [1]; color = GRIS, ls = :dash)
    plot(p1, p2; layout = (1, 2), size = (760, 330))
end

# ╔═╡ e7853d94-ca94-11f1-ac41-997e8df45089
"Datos ordenados de mayor a menor y su CCDF empírica, P(X ≥ x)."
function ccdf_empirica(x)
    s = sort(x; rev = true)
    s, (1:length(s)) ./ length(s)
end

# ╔═╡ e7853d94-ca94-11f1-984d-bb70dd2ec622
"Estimador de Hill con los k datos más grandes, y su error estándar."
function hill_ref(x, k)
    s = sort(x; rev = true)
    α = k / sum(log.(s[1:k] ./ s[k+1]))
    (α = α, se = α / sqrt(k))
end

# ╔═╡ e7853cfe-ca94-11f1-a875-8ff8b6262357
let
    perdidas, ganancias = -z[z .< 0], z[z .> 0]
    αp = hill_ref(perdidas, k_hill)
    αg = hill_ref(ganancias, k_hill)
    Markdown.parse("""
    Con ``k = `` $k_hill: ``\\hat\\alpha`` de las pérdidas
    ``= `` $(round(αp.α, digits = 2)) ± $(round(αp.se, digits = 2)); de las
    ganancias ``= `` $(round(αg.α, digits = 2)) ± $(round(αg.se, digits = 2)).
    """)
end

# ╔═╡ e7853d9e-ca94-11f1-8f9c-f1e474b9a96c
"Estimador de Hill para cada k en `ks`, ordenando los datos una sola vez."
function curva_hill(x, ks)
    s = sort(x; rev = true)
    L = cumsum(log.(s))
    [k / (L[k] - k * log(s[k+1])) for k in ks]
end

# ╔═╡ e7853cf4-ca94-11f1-bec1-6d350f063b2b
let
    perdidas, ganancias = -z[z .< 0], z[z .> 0]
    xp, Pp = ccdf_empirica(perdidas)
    xg, Pg = ccdf_empirica(ganancias)
    xs = 10 .^ range(-1, log10(maximum(abs, z)), length = 100)
    p1 = plot(xp, Pp; color = NARANJA, lw = 2, xscale = :log10, yscale = :log10,
              label = "Pérdidas", xlabel = "x (desviaciones estándar)", ylabel = "P(X > x)",
              legend = :bottomleft, title = "CCDF de las pérdidas y las ganancias",
              xlims = (0.1, maximum(abs, z) * 1.2), ylims = (1e-4, 1.1))
    plot!(p1, xg, Pg; color = AZUL, lw = 2, label = "Ganancias")
    plot!(p1, xs, ccdf.(Normal(), xs) ./ ccdf(Normal(), 0); color = GRIS, ls = :dash, label = "Normal")
    ks = 20:5:400
    p2 = plot(ks, curva_hill(perdidas, ks); color = NARANJA, lw = 2, label = "Pérdidas",
              xlabel = "k (datos de la cola)", ylabel = "α de Hill", legend = :topright,
              ylims = (0, 8), title = "Curva de Hill")
    plot!(p2, ks, curva_hill(ganancias, ks); color = AZUL, lw = 2, label = "Ganancias")
    vline!(p2, [k_hill]; color = GRIS, ls = :dot, label = "k elegido")
    plot(p1, p2; layout = (1, 2), size = (780, 340))
end

# ╔═╡ e7853d9e-ca94-11f1-8083-8fbe346fcb4e
"""
Método de Clauset, Shalizi y Newman (2009). Para cada candidato a `xmin` estima
α por máxima verosimilitud y calcula la distancia de Kolmogorov-Smirnov entre
la cola empírica y la ley de potencias ajustada. Devuelve todas las curvas y,
en `mejor`, el `xmin` de menor distancia.
"""
function clauset(x; candidatos = quantile(x, 0.5:0.005:0.995), nmin = 50)
    s = sort(x; rev = true)
    L = cumsum(log.(s))
    res = (xmin = Float64[], α = Float64[], D = Float64[], n = Int[])
    for xm in unique(candidatos)
        n = searchsortedlast(s, xm; rev = true)        # datos mayores o iguales a xm
        n < nmin && continue
        α = n / (L[n] - n * log(xm))
        cola = @view s[n:-1:1]                         # cola en orden creciente
        D = maximum(abs.((1:n) ./ n .- (1 .- (xm ./ cola) .^ α)))
        push!(res.xmin, xm); push!(res.α, α); push!(res.D, D); push!(res.n, n)
    end
    i = argmin(res.D)
    (; res..., mejor = (xmin = res.xmin[i], α = res.α[i], D = res.D[i], n = res.n[i],
                        se = res.α[i] / sqrt(res.n[i])))
end

# ╔═╡ e7853da8-ca94-11f1-95dd-93122e2a5b2c
"""
Simula un GARCH(1,1) con varianza de largo plazo `s2`, es decir,
ω = s2 (1 − a − b). Descarta 500 pasos iniciales.
"""
function simular_garch(n, a, b, s2; rng = Xoshiro(10))
    ω = s2 * (1 - a - b)
    σ2 = s2
    for _ in 1:500
        σ2 = ω + (a * randn(rng)^2 + b) * σ2
    end
    r, v = zeros(n), zeros(n)
    for t in 1:n
        v[t] = σ2
        r[t] = sqrt(σ2) * randn(rng)
        σ2 = ω + a * r[t]^2 + b * σ2
    end
    (r = r, σ2 = v)
end

# ╔═╡ e7853da8-ca94-11f1-892d-6f9a5eefd873
"Log-verosimilitud normal del GARCH(1,1) con varianza de largo plazo s2."
function garch_loglik(x, a, b, s2)
    ω = s2 * (1 - a - b)
    σ2 = s2
    ℓ = 0.0
    for t in eachindex(x)
        ℓ -= 0.5 * (log(σ2) + x[t]^2 / σ2)
        σ2 = ω + a * x[t]^2 + b * σ2
    end
    ℓ
end

# ╔═╡ e7853db0-ca94-11f1-8fd2-b5f85c48dea7
"""
Ajuste del GARCH(1,1) por máxima verosimilitud en una cuadrícula de (a, b),
con la varianza de largo plazo fijada en la varianza muestral.
"""
function ajustar_garch(r; as = 0.005:0.005:0.4, bs = 0.5:0.0025:0.9975)
    x = r .- mean(r)
    s2 = var(x)
    mejor = (ℓ = -Inf, a = NaN, b = NaN, s2 = s2)
    for a in as, b in bs
        a + b < 0.9995 || continue
        ℓ = garch_loglik(x, a, b, s2)
        ℓ > mejor.ℓ && (mejor = (ℓ = ℓ, a = a, b = b, s2 = s2))
    end
    mejor
end

# ╔═╡ e7853d1c-ca94-11f1-9f39-318644ef10cb
ajuste = ajustar_garch(r);

# ╔═╡ e7853d1c-ca94-11f1-8566-6746659badac
Markdown.parse("""
**GARCH ajustado:** ``a = `` $(ajuste.a), ``b = `` $(ajuste.b),
``a + b = `` $(round(ajuste.a + ajuste.b, digits = 4)). Una perturbación en la
volatilidad tarda unos $(round(Int, log(0.5) / log(ajuste.a + ajuste.b))) días
en reducirse a la mitad.
""")

# ╔═╡ e7853d1c-ca94-11f1-9358-e998172e60db
modelos_t1 = let
    n = length(r)
    gauss = std(r) .* randn(Xoshiro(3), n)
    garch = [simular_garch(n, ajuste.a, ajuste.b, ajuste.s2; rng = Xoshiro(s)).r for s in 1:10]
    (gauss = gauss, garch = garch)
end;

# ╔═╡ e7853d26-ca94-11f1-9ad9-bf754e913830
let
    fila(nombre, x) = (nombre, kurtosis(x), autocor(x, [1])[1], mean(autocor(abs.(x), 1:20)),
                       hill_ref(-filter(<(0), (x .- mean(x)) ./ std(x)), 120).α)
    filas = [fila("**$ticker** (datos)", r), fila("Caminata gaussiana", modelos_t1.gauss)]
    g = [fila("", x) for x in modelos_t1.garch]
    f(v) = round(v, digits = 2)
    rango(i) = "$(f(mean(x[i] for x in g))) ($(f(minimum(x[i] for x in g))) a $(f(maximum(x[i] for x in g))))"
    lineas = ["| $(x[1]) | $(f(x[2])) | $(round(x[3], digits = 3)) | $(round(x[4], digits = 3)) | $(f(x[5])) |" for x in filas]
    push!(lineas, "| GARCH(1,1), 10 simulaciones: media (mín. a máx.) | $(rango(2)) | $(rango(3)) | $(rango(4)) | $(rango(5)) |")
    Markdown.parse("""
    | Serie | Curtosis en exceso | ACF de r, rezago 1 | ACF de ∣r∣, promedio rezagos 1–20 | α de Hill de las pérdidas (k = 120) |
    |---|---:|---:|---:|---:|
    $(join(lineas, "\n"))
    """)
end

# ╔═╡ e7853db0-ca94-11f1-9e59-15194a5ee906
const RETORNOS_PRUEBA = 0.01 .* randn(Xoshiro(2026), 500);

# ╔═╡ e7853db0-ca94-11f1-be32-85e5f3c5f86a
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

# ╔═╡ e7853c90-ca94-11f1-8cb5-716a8e916dfd
verificar(() -> exceso_extremos([0.0, 3.0, -3.0, 1.0, 0.5], 2), 2 / (5 * 2ccdf(Normal(), 2));
          pista = md"Cuenta con `count(x -> abs(x) > k, z)` y divide entre `2 * length(z) * ccdf(Normal(), k)`.")

# ╔═╡ e7853ca4-ca94-11f1-a530-0f6b51ca0822
verificar(() -> autocorrelacion([1.0, 2.0, 3.0, 4.0, 5.0], 1), 0.4;
          pista = md"Con `d = x .- mean(x)`, calcula `sum(d[1:end-τ] .* d[1+τ:end]) / sum(d .^ 2)`.")

# ╔═╡ e7853ccc-ca94-11f1-8b0b-11c1f3972a01
verificar(() -> curtosis_agregada(RETORNOS_PRUEBA, [1, 5]),
          [kurtosis(RETORNOS_PRUEBA), kurtosis(agregar(RETORNOS_PRUEBA, 5))];
          pista = md"Para cada `Δt`, los bloques son `r[i:i+Δt-1]` con `i` en `1:Δt:length(r)-Δt+1`.")

# ╔═╡ e7853cf4-ca94-11f1-bcb9-8161ccd29a46
verificar(() -> hill_alpha([1.0, 8.0, 2.0, 4.0], 2), 2 / (3log(2));
          pista = md"Ordena con `s = sort(x; rev = true)` y calcula `k / sum(log.(s[1:k] ./ s[k+1]))`.")

# ╔═╡ e7853d44-ca94-11f1-9b05-fdf2641169a5
verificar(() -> beta_mercado([0.0, 3.0, 4.0, 9.0], [0.0, 1.0, 2.0, 3.0]), 2.8;
          pista = md"Usa `cov(ra, rm) / var(rm)`.")

# ╔═╡ Cell order:
# ╟─e7853c22-ca94-11f1-b634-0d9eaca1d1cf
# ╟─e7853c54-ca94-11f1-81df-17461ae29fb5
# ╟─e7853c5e-ca94-11f1-b0af-f3a5ab2880e0
# ╠═e7853c68-ca94-11f1-a3d9-a7c5f6965898
# ╟─e7853c68-ca94-11f1-a074-e71355470b65
# ╠═e7853c72-ca94-11f1-b874-b1e9836f9b9d
# ╟─e7853c7a-ca94-11f1-9e40-e74a0b7b4b5c
# ╟─e7853c7a-ca94-11f1-a4e7-1bc765bf436c
# ╟─e7853c86-ca94-11f1-a05a-a33cf3adcd06
# ╠═e7853c86-ca94-11f1-9c6a-0b68fe6903c2
# ╟─e7853c90-ca94-11f1-8cb5-716a8e916dfd
# ╟─e7853c90-ca94-11f1-ba49-b58c15b9f313
# ╟─e7853c9a-ca94-11f1-a0ab-2d02c770eff5
# ╟─e7853c9a-ca94-11f1-965f-3d14ca109da5
# ╠═e7853ca4-ca94-11f1-a5d2-bf1acdd332ec
# ╟─e7853ca4-ca94-11f1-a530-0f6b51ca0822
# ╟─e7853cac-ca94-11f1-916c-25e5a11e8e86
# ╟─e7853cac-ca94-11f1-bc46-8394d0b7dac0
# ╠═e7853cb8-ca94-11f1-ac84-c1cd2da7790e
# ╟─e7853cb8-ca94-11f1-be9b-15a60fc9c5c3
# ╠═e7853cb8-ca94-11f1-9514-3776181b7b9c
# ╟─e7853cc2-ca94-11f1-9cbb-a71cbaa66298
# ╠═e7853cc2-ca94-11f1-af1c-a1e0df7b806e
# ╟─e7853ccc-ca94-11f1-8b0b-11c1f3972a01
# ╟─e7853ccc-ca94-11f1-9336-2349d764710c
# ╟─e7853ccc-ca94-11f1-93df-d94292462c3d
# ╟─e7853cd6-ca94-11f1-839c-af4f97ddc658
# ╠═e7853cd6-ca94-11f1-a19a-3781958cf69c
# ╟─e7853ce0-ca94-11f1-8ea9-07f06bb7482b
# ╠═e7853ce0-ca94-11f1-9bad-454a648b73bd
# ╟─e7853cea-ca94-11f1-a735-a36ecf99948a
# ╠═e7853cea-ca94-11f1-84f6-1103f9e01b47
# ╟─e7853cf4-ca94-11f1-bcb9-8161ccd29a46
# ╟─e7853cf4-ca94-11f1-9aaf-419c1e766391
# ╟─e7853cf4-ca94-11f1-bec1-6d350f063b2b
# ╟─e7853cfe-ca94-11f1-a875-8ff8b6262357
# ╟─e7853cfe-ca94-11f1-8d4b-6302417c207d
# ╠═e7853d08-ca94-11f1-a4ff-1dff5d2821ca
# ╟─e7853d08-ca94-11f1-96b1-db4fc47b8fa5
# ╠═e7853d12-ca94-11f1-ad1f-ffeaf7f917d4
# ╟─e7853d12-ca94-11f1-b0e1-4591a8d673ad
# ╠═e7853d1c-ca94-11f1-9f39-318644ef10cb
# ╟─e7853d1c-ca94-11f1-8566-6746659badac
# ╠═e7853d1c-ca94-11f1-9358-e998172e60db
# ╟─e7853d26-ca94-11f1-9ad9-bf754e913830
# ╟─e7853d26-ca94-11f1-acc5-3d38bc1b80c6
# ╠═e7853d30-ca94-11f1-8a1b-637066f729ed
# ╟─e7853d30-ca94-11f1-92e6-454bccea1085
# ╠═e7853d3a-ca94-11f1-92e0-095739837d0d
# ╟─e7853d3a-ca94-11f1-8ff0-5110b44fb1c4
# ╠═e7853d3a-ca94-11f1-bb8f-f7e1348d9a84
# ╟─e7853d44-ca94-11f1-9b05-fdf2641169a5
# ╠═e7853d44-ca94-11f1-b635-b7a0b913c159
# ╠═e7853d50-ca94-11f1-9599-a54504cc3b04
# ╟─e7853d50-ca94-11f1-a662-3565a78e0d1e
# ╟─e7853d58-ca94-11f1-9027-ffd10be68603
# ╟─e7853d58-ca94-11f1-9785-3b2afac4e1a1
# ╠═e7853d62-ca94-11f1-863e-a113bc5e17c8
# ╟─e7853d62-ca94-11f1-b0e3-9504f93446e6
# ╠═e7853d62-ca94-11f1-9e57-613ec124b35e
# ╟─e7853d6c-ca94-11f1-a531-35c1611ffbd0
# ╟─e7853d6c-ca94-11f1-8fc9-db9ca03a005d
# ╠═e7853d76-ca94-11f1-806f-978633577a7a
# ╠═e7853d76-ca94-11f1-821c-7b6b39e509ee
# ╠═e7853d76-ca94-11f1-a909-b12bcf3a01c6
# ╠═e7853d82-ca94-11f1-a08e-1b70ab63c8c8
# ╠═e7853d82-ca94-11f1-b283-11dcd210eb1b
# ╠═e7853d8a-ca94-11f1-978f-c944f038e06d
# ╠═e7853d8a-ca94-11f1-bb04-f9f9b9557c5c
# ╠═e7853d94-ca94-11f1-9b5e-e1a63d331131
# ╠═e7853d94-ca94-11f1-ac41-997e8df45089
# ╠═e7853d94-ca94-11f1-984d-bb70dd2ec622
# ╠═e7853d9e-ca94-11f1-8f9c-f1e474b9a96c
# ╠═e7853d9e-ca94-11f1-8083-8fbe346fcb4e
# ╠═e7853da8-ca94-11f1-95dd-93122e2a5b2c
# ╠═e7853da8-ca94-11f1-892d-6f9a5eefd873
# ╠═e7853db0-ca94-11f1-8fd2-b5f85c48dea7
# ╠═e7853db0-ca94-11f1-9e59-15194a5ee906
# ╠═e7853db0-ca94-11f1-be32-85e5f3c5f86a
