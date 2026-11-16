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

# ╔═╡ 5c27066c-ca92-11f1-b8a1-f7014b533d85
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

# ╔═╡ 5c27057c-ca92-11f1-ac41-1d45ed46beac
md"""
# Sesión 1: Retornos y colas pesadas

**Objetivos de la sesión**
1. Entender qué busca la econofísica y por qué empieza por los datos.
2. Repasar retornos logarítmicos, media, varianza y curtosis.
3. Medir en datos reales los dos primeros **hechos estilizados** de los
   mercados financieros: **colas pesadas** y **ausencia de autocorrelación
   lineal**.

**Cómo usar este cuaderno.** Pluto es *reactivo*: al mover un deslizador o
cambiar un selector, todas las celdas que dependen de ese valor se recalculan
solas. No hace falta ejecutar nada en orden. Las secciones marcadas
**[Avanzado]** son opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s01_retornos.jl`. La primera vez, Julia instala los paquetes
del curso.

**Repaso recomendado** (en `Class_notes/repaso/`): C.B1 (logaritmos), P.B1–B5
(variables aleatorias, media, varianza, normal), E.B1 (procesos
estocásticos), S.B1–B3 (estimación y pruebas). Para las notas formales:
P.I1–I2, E.I1–I2 (autocorrelación y martingalas), C.I1–I2 y S.I2. Ruta
avanzada: P.A5.

**Diagnóstico.** Antes del laboratorio, resuelve las diez preguntas de
`Class_notes/repaso/diagnostico_inicial.md` (25 minutos, no cuenta para la
calificación). Al final te dice qué claves de repaso te conviene leer.

**Conceptos transversales.** Hoy se presentan los **retornos logarítmicos**,
las **distribuciones** (histograma y escala logarítmica; la CCDF llega en la
sesión 3) y la **media, la varianza y la correlación** (que se completan en la
sesión 2).
"""

# ╔═╡ 5c2705a4-ca92-11f1-8ec0-e1340af79203
TableOfContents(title = "Contenido")

# ╔═╡ 5c2705b0-ca92-11f1-bb7d-c3c171fe80c7
md"""
## 0. El fenómeno: ¿el mercado es un volado?

Un modelo muy usado dice que el precio de una acción sube o baja cada día
como un volado, con movimientos de tamaño típico parecido. Si fuera así, los
días con caídas enormes serían rarísimos. Veamos el peor día del S&P 500 en
los datos del curso:
"""

# ╔═╡ 5c2705cc-ca92-11f1-861c-01b564f4c0f1
md"""
**¿Qué es la econofísica?** La economía tradicional suele partir de un agente
representativo que optimiza. La econofísica parte de otra pregunta: **¿qué
patrones colectivos surgen cuando muchos agentes heterogéneos interactúan?**
Es la misma pregunta que hace la física estadística con átomos y moléculas:
nadie predice la trayectoria de una molécula, pero la presión de un gas sí se
predice.

Antes de construir modelos, un físico **mide**. Kepler describió las órbitas
antes de que Newton las explicara. En finanzas, las regularidades empíricas que
se repiten en casi todos los mercados y periodos se llaman **hechos
estilizados** (Cont, 2001). Hoy medimos los dos primeros y en la sesión 2 los
otros dos. Durante el curso construiremos modelos que intenten explicarlos, y
en la sesión 20 volveremos a estas mediciones para ver cuáles reproducen.
"""

# ╔═╡ 5c2705cc-ca92-11f1-9515-752c4bcba805
md"""
## 1. Repaso: retornos

Si ``P_t`` es el precio de un activo en el día ``t``, hay dos formas de medir
cuánto cambió:

```math
R_t = \frac{P_t - P_{t-1}}{P_{t-1}} \quad \text{(retorno simple)}, \qquad
r_t = \ln\frac{P_t}{P_{t-1}} = \ln(1 + R_t) \quad \text{(retorno logarítmico)}.
```

Para cambios pequeños son casi iguales, porque ``\ln(1+x) \approx x``. Mueve el
deslizador para ver cuándo dejan de parecerse.
"""

# ╔═╡ 5c2705d6-ca92-11f1-a03e-6bd3ab7c7af8
md"""
Retorno simple ``R``: $(@bind R_simple Slider(-0.5:0.01:0.5, default = 0.1, show_value = true))
"""

# ╔═╡ 5c2705e2-ca92-11f1-9a75-9dbee15bd055
md"""
Con ``R = `` $(round(100R_simple, digits = 1))%, el retorno logarítmico es
``r = `` $(round(100log1p(R_simple), digits = 2))%. La diferencia es
$(round(100abs(R_simple - log1p(R_simple)), digits = 2)) puntos porcentuales.

**¿Por qué usamos ``r``?** Porque **se suma en el tiempo**: el retorno
logarítmico de varios días es la suma de los retornos diarios,

```math
\ln\frac{P_{t+k}}{P_t} = r_{t+1} + r_{t+2} + \dots + r_{t+k}.
```

Eso nos permite usar todo lo que sabemos sobre sumas de variables aleatorias,
como el teorema del límite central.

**Tres estadísticas que usaremos:**

| Medida | Fórmula | Qué mide |
|---|---|---|
| Media | ``\mu = \frac{1}{N}\sum r_t`` | Tendencia promedio |
| Desviación estándar | ``\sigma = \sqrt{\frac{1}{N-1}\sum (r_t - \mu)^2}`` | **Volatilidad**: el tamaño típico de un movimiento |
| Curtosis en exceso | ``\kappa = \frac{1}{N}\sum \left(\frac{r_t-\mu}{\sigma}\right)^4 - 3`` | Peso de las colas: la normal tiene ``\kappa = 0`` |
"""

# ╔═╡ 5c2705e2-ca92-11f1-9003-d5c34654de92
md"""
## 2. Los datos

Usamos precios de cierre diarios de cuatro activos muy distintos: un índice bursátil mexicano, un tipo de cambio, un índice
estadounidense y una criptomoneda. Los datos se descargaron de Yahoo Finance
con `Code/descargar_datos.jl` y están en `Data/precios.csv`.
"""

# ╔═╡ 5c2705ea-ca92-11f1-8a0f-9b6239542270
datos = CSV.read(joinpath(@__DIR__, "..", "Data", "precios.csv"), DataFrame)

# ╔═╡ 5c2705ea-ca92-11f1-9651-3b9586606eb9
md"""
**Activo:** $(@bind activo Select(unique(datos.activo)))
"""

# ╔═╡ 5c2705f4-ca92-11f1-bd2f-55511de8a0a6
begin
    serie = sort(filter(:activo => ==(activo), datos), :fecha)
    fechas = serie.fecha
    precios = serie.cierre
    r = diff(log.(precios))      # retornos logarítmicos diarios
    fechas_r = fechas[2:end]
    μ, σ = mean(r), std(r)
end;

# ╔═╡ 5c2705f4-ca92-11f1-b5cc-6190bddb1971
md"""
$(@bind escala_log CheckBox(default = true)) Precio en escala logarítmica
"""

# ╔═╡ 5c2705fe-ca92-11f1-ba41-e1cfe06c2810
md"""
Las líneas punteadas marcan tres episodios: la **crisis financiera** (octubre
de 2008), la **elección presidencial en EE. UU.** (9 de noviembre de 2016) y el
**inicio de la pandemia de COVID-19** (marzo de 2020).

Observa la gráfica de retornos: fluctúan alrededor de cero, pero **no con la
misma intensidad todo el tiempo**. Hay periodos de calma y periodos de
turbulencia. Esa observación visual es el punto de partida de todo lo que
sigue.
"""

# ╔═╡ 5c270608-ca92-11f1-be55-a181a1eff8e3
md"""
## 3. Hecho 1: colas pesadas

Si los retornos fueran normales, casi todos los días caerían a menos de 3
desviaciones estándar de la media. Para comparar activos con volatilidades
distintas, **estandarizamos**: ``z_t = (r_t - \mu)/\sigma``. Así, ``z = 4``
significa "un movimiento de 4 desviaciones estándar".

La gráfica compara la densidad empírica de ``z`` con la normal estándar.
**Marca la casilla del eje logarítmico**: en escala lineal las colas son
invisibles; en escala logarítmica saltan a la vista.
"""

# ╔═╡ 5c270614-ca92-11f1-abd1-99c6e4c25a38
z = (r .- μ) ./ σ;

# ╔═╡ 5c270614-ca92-11f1-8920-63527588838a
md"""
$(@bind colas_log CheckBox(default = false)) Eje vertical logarítmico
"""

# ╔═╡ 5c27061c-ca92-11f1-934d-31c806ef768f
md"""
**Curtosis en exceso de $activo:** $(round(kurtosis(r), digits = 1)) (la normal tiene 0).

¿Qué tan raros son los eventos extremos? La tabla compara cuántos días con
``|z| > k`` hubo en los datos contra cuántos predice la normal.
"""

# ╔═╡ 5c270630-ca92-11f1-94a1-97a20aa902d6
md"""
## 4. Hecho 2: ausencia de autocorrelación lineal

La **autocorrelación** al rezago ``\tau`` mide si el retorno de hoy se parece
al de hace ``\tau`` días:

```math
\rho(\tau) = \mathrm{corr}(r_t,\, r_{t+\tau}).
```

Si ``\rho(\tau)`` fuera claramente positiva, podríamos predecir el mercado:
"si hoy subió, mañana subirá". La franja gris marca el intervalo
``\pm 1.96/\sqrt{N}``: las barras dentro de ella son indistinguibles de cero.
"""

# ╔═╡ 5c270630-ca92-11f1-acf1-9105dba75bcb
md"""
Rezagos: $(@bind L Slider(5:5:100, default = 30, show_value = true))
"""

# ╔═╡ 5c27063a-ca92-11f1-a125-0f03d417ac62
md"""
Después del primer rezago, casi todas las barras quedan dentro de la franja.
Esto es consistente con la **hipótesis de mercado eficiente**: si los retornos
pasados predijeran los futuros, alguien explotaría esa regularidad y, al
hacerlo, la borraría.

En el rezago 1 suele aparecer una correlación pequeña pero significativa:
positiva en el IPC (≈ 0.08) y negativa en el S&P 500 (≈ −0.10). Una
correlación de 0.1 explica apenas el 1% de la varianza del día siguiente, y
los costos de transacción suelen impedir aprovecharla. *Pregunta:* ¿por qué un
índice con muchas acciones poco negociadas, como el IPC, podría mostrar
correlación positiva en el rezago 1?

Pero cuidado: *ausencia de correlación lineal no significa independencia*.
En la sesión 2 veremos que el **tamaño** de los movimientos sí tiene memoria.
"""

# ╔═╡ 5c270642-ca92-11f1-9418-099053f88683
md"""
## Resumen

| Hecho estilizado | Cómo se mide | Lo que vimos |
|---|---|---|
| Colas pesadas | Densidad en escala log; curtosis; conteo de eventos extremos | Eventos de 5σ o más son comunes |
| Sin autocorrelación lineal | ACF de ``r`` | Pequeña en el rezago 1 e indistinguible de cero después |
"""

# ╔═╡ 5c270642-ca92-11f1-a669-0163a14e6219
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe la función `mis_retornos(p)` que reciba un vector de precios y
devuelva el vector de retornos logarítmicos. Reemplaza `missing` por tu código.
"""

# ╔═╡ 5c27064e-ca92-11f1-8709-c599751aacae
function mis_retornos(p)
    missing  # reemplaza esta línea
end

# ╔═╡ 5c27064e-ca92-11f1-be9c-076a076f3a17
md"""
**2.** Escribe la función `curtosis_exceso(x)` usando la fórmula de la sección 1.
"""

# ╔═╡ 5c270658-ca92-11f1-9a43-fba29ab635f1
function curtosis_exceso(x)
    missing  # reemplaza esta línea
end

# ╔═╡ 5c270662-ca92-11f1-89a4-cd0b5a485563
md"""
### Ruta básica: interpretación

3. Ordena los cuatro activos de colas más pesadas a menos pesadas usando la
   curtosis de cada uno. ¿Coincide con tu intuición sobre qué tan riesgoso es
   cada uno?
4. Según la tabla de eventos extremos, ¿cada cuánto debería ocurrir el peor día
   del S&P 500 si los retornos fueran normales? ¿Cuántos días con ``|z| > 5``
   hubo en realidad?
5. En el rezago 1, el IPC muestra una autocorrelación positiva pequeña. ¿Por
   qué un índice con muchas acciones poco negociadas podría mostrarla? ¿Se
   podría ganar dinero con ella?
6. Un compañero dice: "los retornos tienen media casi cero y el histograma
   parece una campana, así que son normales". ¿Qué gráfica de este cuaderno
   usarías para convencerlo de lo contrario y qué le pedirías que mirara?

### Ruta avanzada

7. Demuestra que si el exceso de retorno es una diferencia de martingala, su
   autocorrelación es cero en todos los rezagos (notas formales, sección 4).
8. La distribución ``t`` de Student con ``\nu > 4`` grados de libertad tiene
   curtosis en exceso ``6/(\nu - 4)``. ¿Qué ``\nu`` daría la curtosis que mediste
   en cada activo? ¿Qué problema aparece si la curtosis es muy grande?

## Para la próxima sesión

- Lectura: Mantegna y Stanley (1999), capítulos 1 y 2. Complementaria: Slanina
  (2014), secciones 1.1 y 1.2.
- En la sesión 2 miramos el **tamaño** de los movimientos: ¿un día turbulento
  anuncia otro? Y aprenderemos a comparar los datos contra **modelos nulos**.
"""

# ╔═╡ 5c27066c-ca92-11f1-bae4-d1fe25583445
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 5c270674-ca92-11f1-9bcc-2536a770240e
begin
    const AZUL    = "#2a78d6"   # serie principal
    const NARANJA = "#eb6834"   # segunda serie
    const GRIS    = "#52514e"   # referencias (normal, bandas, eventos)
end;

# ╔═╡ 5c2705d6-ca92-11f1-9b75-a5b99207b5a4
let
    x = -0.5:0.01:0.5
    p = plot(x, x, color = GRIS, ls = :dash, lw = 1.5,
             xlabel = "Retorno simple R", ylabel = "Retorno logarítmico r",
             title = "r = ln(1 + R)  vs.  r = R")
    plot!(p, x, log1p.(x), color = AZUL, lw = 2)
    scatter!(p, [R_simple], [log1p(R_simple)], color = AZUL, ms = 6,
             markerstrokecolor = :white, markerstrokewidth = 2)
    annotate!(p, -0.45, 0.35, text("línea punteada: r = R", 8, :left, GRIS))
    p
end

# ╔═╡ 5c270674-ca92-11f1-a637-05b607e11d11
const EVENTOS = [Date(2008, 10, 10), Date(2016, 11, 9), Date(2020, 3, 16)];

# ╔═╡ 5c270680-ca92-11f1-bb14-c5033fabf7cb
"Dibuja líneas verticales en los eventos que caen dentro del rango de `fechas`."
function marcar_eventos!(p, fechas)
    dentro = filter(e -> first(fechas) <= e <= last(fechas), EVENTOS)
    vline!(p, dentro, color = GRIS, ls = :dot, lw = 1, label = "")
end

# ╔═╡ 5c270680-ca92-11f1-af61-1f5c9df7b5f5
"Densidad empírica en bins de ancho `ancho`; descarta bins vacíos (para escala log)."
function densidad_empirica(z; ancho = 0.25)
    h = fit(Histogram, z, -16:ancho:16)
    centros = collect(h.edges[1][1:end-1]) .+ ancho / 2
    dens = h.weights ./ (length(z) * ancho)
    keep = dens .> 0
    centros[keep], dens[keep]
end

# ╔═╡ 5c270680-ca92-11f1-b13e-35a72351a63c
dias_por_año(activo) = activo == "Bitcoin" ? 365 : 252

# ╔═╡ 5c27068a-ca92-11f1-b1e8-8b03f59a0c3c
"Límite simétrico del eje z: al menos ±8 y siempre mayor que el dato más extremo."
rango_z(z) = max(8, ceil(Int, maximum(abs, z)) + 1)

# ╔═╡ 5c270614-ca92-11f1-ab67-35fedfdac8fd
let
    x, dens = densidad_empirica(z)
    zmax = rango_z(z)
    p = plot(-zmax:0.05:zmax, pdf.(Normal(), -zmax:0.05:zmax), color = GRIS, ls = :dash,
             lw = 1.5, label = "Normal estándar",
             yscale = colas_log ? :log10 : :identity,
             ylims = colas_log ? (1e-5, 1) : (0, 0.75),
             xlims = colas_log ? (-zmax, zmax) : (-8, 8), legend = :topright,
             xlabel = "z (desviaciones estándar)", ylabel = "Densidad",
             title = "$activo: distribución de retornos estandarizados")
    scatter!(p, x, dens, color = AZUL, ms = 4, markerstrokewidth = 0,
             label = "Datos")
    p
end

# ╔═╡ 5c27068a-ca92-11f1-ba2e-c5d693c455fd
"Marcas del eje x: un año de cada `paso`, etiquetado sólo con el año."
function ticks_años(fechas; paso = 4)
    años = year(first(fechas)):paso:year(last(fechas))
    marcas = Date.(años, 1, 1)
    (Dates.value.(marcas), string.(años))
end

# ╔═╡ 5c270694-ca92-11f1-8e34-7d7851886bb8
"Marcas del eje de precios en escala log: 1, 2 y 5 por cada potencia de 10, con separador de miles."
function ticks_precio(precios)
    lo, hi = extrema(precios)
    marcas = [m * 10.0^k for k in floor(Int, log10(lo)):ceil(Int, log10(hi)) for m in (1, 2, 5)]
    marcas = filter(x -> lo <= x <= hi, marcas)
    etiqueta(x) = x < 1 ? string(round(x, sigdigits = 2)) :
                  replace(string(round(Int, x)), r"(?<=\d)(?=(\d{3})+$)" => ",")
    (marcas, etiqueta.(marcas))
end

# ╔═╡ 5c2705fe-ca92-11f1-b67d-9d86bb77da16
let
    p1 = plot(fechas, precios, color = AZUL, lw = 1.2,
              yscale = escala_log ? :log10 : :identity,
              yticks = escala_log ? ticks_precio(precios) : :auto,
              ylabel = "Precio", title = "$activo: precio de cierre",
              xticks = ticks_años(fechas))
    p2 = plot(fechas_r, 100r, color = AZUL, lw = 0.6,
              ylabel = "r (%)", title = "Retornos logarítmicos diarios",
              xticks = ticks_años(fechas))
    for p in (p1, p2)
        marcar_eventos!(p, fechas)
    end
    plot(p1, p2, layout = (2, 1), size = (720, 520))
end

# ╔═╡ 5c270694-ca92-11f1-b530-c3a355976e8b
"Frecuencia de un evento de `k` desviaciones estándar si los retornos fueran normales."
function frecuencia_normal(k, dpa)
    p = 2ccdf(Normal(), k)          # probabilidad diaria
    por_año = p * dpa
    por_año >= 1 && return "unas $(round(Int, por_año)) veces al año"
    años = 1 / por_año
    años < 1.5  && return "aproximadamente una vez al año"
    años < 1e3  && return "una vez cada $(round(Int, años)) años"
    años < 1e6  && return "una vez cada $(round(Int, años / 1e3)) mil años"
    años < 1e9  && return "una vez cada $(round(Int, años / 1e6)) millones de años"
    años < 1.38e10 && return "una vez cada $(round(años / 1e9, digits = 1)) mil millones de años"
    "menos de una vez en toda la edad del universo"
end

# ╔═╡ 5c2705b8-ca92-11f1-9ef4-e1905175cd72
peor_sp = let
    sp = sort(filter(:activo => ==("S&P 500"), datos), :fecha)
    r_sp = diff(log.(sp.cierre))
    z_sp = (r_sp .- mean(r_sp)) ./ std(r_sp)
    i = argmin(r_sp)
    (fechas = sp.fecha, r = r_sp, i = i, fecha = sp.fecha[i + 1],
     caida = round(-100expm1(r_sp[i]), digits = 1), z = round(abs(z_sp[i]), digits = 1),
     frec = frecuencia_normal(abs(z_sp[i]), 252))
end;

# ╔═╡ 5c2705b8-ca92-11f1-8a55-73837b625ee8
let s = peor_sp
    p = plot(s.fechas[2:end], 100s.r, color = AZUL, lw = 0.6, ylabel = "r (%)",
             title = "S&P 500: retornos diarios", xticks = ticks_años(s.fechas))
    scatter!(p, [s.fechas[s.i + 1]], [100s.r[s.i]], color = NARANJA, ms = 6, msw = 0)
end

# ╔═╡ 5c2705c2-ca92-11f1-8772-3dde7d247402
md"""
El día $(peor_sp.fecha) el S&P 500 cayó **un $(peor_sp.caida)%** en una sola
sesión (punto naranja): un movimiento de $(peor_sp.z) desviaciones estándar.
Si los retornos fueran normales, un día así ocurriría $(peor_sp.frec). Y sin
embargo ocurrió, junto con muchos otros días casi igual de extremos.
"""

# ╔═╡ 5c270626-ca92-11f1-bdd1-8b8ab2b098df
extremo = let i = argmax(abs.(z))
    (fecha = fechas_r[i], r = round(100r[i], digits = 1),
     z = round(z[i], digits = 1),
     frecuencia = frecuencia_normal(abs(z[i]), dias_por_año(activo)))
end;

# ╔═╡ 5c270626-ca92-11f1-b76d-d992af0ae9f8
md"""
El movimiento más extremo de $activo fue el **día $(extremo.fecha)**:
``r = `` $(extremo.r)%, es decir, ``z = `` $(extremo.z).
Según la distribución normal, **un movimiento así ocurriría
$(extremo.frecuencia)**.
"""

# ╔═╡ 5c27069e-ca92-11f1-9fe4-7b26d095223e
function tabla_extremos(z, activo)
    N = length(z)
    dpa = dias_por_año(activo)
    DataFrame(
        "Umbral |z| >" => 3:8,
        "Días observados" => [count(>(k), abs.(z)) for k in 3:8],
        "Días esperados (normal)" => [round(N * 2ccdf(Normal(), k), sigdigits = 2) for k in 3:8],
        "Frecuencia según la normal" => [frecuencia_normal(k, dpa) for k in 3:8],
    )
end

# ╔═╡ 5c27061c-ca92-11f1-b0bb-2d35c3c78901
tabla_extremos(z, activo)

# ╔═╡ 5c27069e-ca92-11f1-b741-3b8d87e42afa
function grafica_acf(x, L; titulo = "", color = AZUL)
    banda = 1.96 / sqrt(length(x))
    p = plot(1:L, fill(banda, L), fillrange = fill(-banda, L), color = GRIS,
             fillalpha = 0.15, lw = 0, xlabel = "Rezago τ (días)",
             ylabel = "ρ(τ)", title = titulo, ylims = (-0.15, 0.4))
    bar!(p, 1:L, autocor(x, 1:L), color = color, lw = 0, bar_width = 0.6)
    hline!(p, [0], color = :black, lw = 0.5)
end

# ╔═╡ 5c27063a-ca92-11f1-bab5-bf53d531654c
grafica_acf(r, L; titulo = "$activo: autocorrelación de los retornos")

# ╔═╡ 5c2706a6-ca92-11f1-8257-6359fd5b9dc8
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

# ╔═╡ 5c27064e-ca92-11f1-a619-7dbf3ebdfce0
verificar(() -> mis_retornos([100.0, 110.0, 99.0, 99.0]),
          [log(1.1), log(0.9), 0.0];
          pista = md"Usa `diff` y `log.`. Recuerda que ``r_t = \ln P_t - \ln P_{t-1}``.")

# ╔═╡ 5c2706b2-ca92-11f1-8151-63ee24de2a45
const VECTOR_PRUEBA = rand(Xoshiro(2026), TDist(5), 100_000);

# ╔═╡ 5c270658-ca92-11f1-bf2f-d770430f0f3a
verificar(() -> curtosis_exceso(VECTOR_PRUEBA), kurtosis(VECTOR_PRUEBA); rtol = 1e-2,
          pista = md"Estandariza con `(x .- mean(x)) ./ std(x)`, eleva a la cuarta potencia con `.^4`, promedia con `mean` y resta 3.")

# ╔═╡ Cell order:
# ╟─5c27057c-ca92-11f1-ac41-1d45ed46beac
# ╟─5c2705a4-ca92-11f1-8ec0-e1340af79203
# ╟─5c2705b0-ca92-11f1-bb7d-c3c171fe80c7
# ╟─5c2705b8-ca92-11f1-9ef4-e1905175cd72
# ╟─5c2705b8-ca92-11f1-8a55-73837b625ee8
# ╟─5c2705c2-ca92-11f1-8772-3dde7d247402
# ╟─5c2705cc-ca92-11f1-861c-01b564f4c0f1
# ╟─5c2705cc-ca92-11f1-9515-752c4bcba805
# ╟─5c2705d6-ca92-11f1-a03e-6bd3ab7c7af8
# ╟─5c2705d6-ca92-11f1-9b75-a5b99207b5a4
# ╟─5c2705e2-ca92-11f1-9a75-9dbee15bd055
# ╟─5c2705e2-ca92-11f1-9003-d5c34654de92
# ╠═5c2705ea-ca92-11f1-8a0f-9b6239542270
# ╟─5c2705ea-ca92-11f1-9651-3b9586606eb9
# ╠═5c2705f4-ca92-11f1-bd2f-55511de8a0a6
# ╟─5c2705f4-ca92-11f1-b5cc-6190bddb1971
# ╟─5c2705fe-ca92-11f1-b67d-9d86bb77da16
# ╟─5c2705fe-ca92-11f1-ba41-e1cfe06c2810
# ╟─5c270608-ca92-11f1-be55-a181a1eff8e3
# ╠═5c270614-ca92-11f1-abd1-99c6e4c25a38
# ╟─5c270614-ca92-11f1-8920-63527588838a
# ╟─5c270614-ca92-11f1-ab67-35fedfdac8fd
# ╟─5c27061c-ca92-11f1-934d-31c806ef768f
# ╟─5c27061c-ca92-11f1-b0bb-2d35c3c78901
# ╟─5c270626-ca92-11f1-bdd1-8b8ab2b098df
# ╟─5c270626-ca92-11f1-b76d-d992af0ae9f8
# ╟─5c270630-ca92-11f1-94a1-97a20aa902d6
# ╟─5c270630-ca92-11f1-acf1-9105dba75bcb
# ╟─5c27063a-ca92-11f1-bab5-bf53d531654c
# ╟─5c27063a-ca92-11f1-a125-0f03d417ac62
# ╟─5c270642-ca92-11f1-9418-099053f88683
# ╟─5c270642-ca92-11f1-a669-0163a14e6219
# ╠═5c27064e-ca92-11f1-8709-c599751aacae
# ╟─5c27064e-ca92-11f1-a619-7dbf3ebdfce0
# ╟─5c27064e-ca92-11f1-be9c-076a076f3a17
# ╠═5c270658-ca92-11f1-9a43-fba29ab635f1
# ╟─5c270658-ca92-11f1-bf2f-d770430f0f3a
# ╟─5c270662-ca92-11f1-89a4-cd0b5a485563
# ╟─5c27066c-ca92-11f1-bae4-d1fe25583445
# ╠═5c27066c-ca92-11f1-b8a1-f7014b533d85
# ╠═5c270674-ca92-11f1-9bcc-2536a770240e
# ╠═5c270674-ca92-11f1-a637-05b607e11d11
# ╠═5c270680-ca92-11f1-bb14-c5033fabf7cb
# ╠═5c270680-ca92-11f1-af61-1f5c9df7b5f5
# ╠═5c270680-ca92-11f1-b13e-35a72351a63c
# ╠═5c27068a-ca92-11f1-b1e8-8b03f59a0c3c
# ╠═5c27068a-ca92-11f1-ba2e-c5d693c455fd
# ╠═5c270694-ca92-11f1-8e34-7d7851886bb8
# ╠═5c270694-ca92-11f1-b530-c3a355976e8b
# ╠═5c27069e-ca92-11f1-9fe4-7b26d095223e
# ╠═5c27069e-ca92-11f1-b741-3b8d87e42afa
# ╠═5c2706a6-ca92-11f1-8257-6359fd5b9dc8
# ╠═5c2706b2-ca92-11f1-8151-63ee24de2a45
