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

# ╔═╡ 7dfd07e6-ca74-11f1-adea-b7ebfbd6c03e
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

# ╔═╡ 7dfd06ce-ca74-11f1-8c99-cb32088a2de7
md"""
# Sesiones 1 y 2: Econofísica y hechos estilizados

**Objetivos de la sesión**
1. Entender qué busca la econofísica y por qué empieza por los datos.
2. Repasar retornos logarítmicos, media, varianza y curtosis.
3. Medir en datos reales los cuatro **hechos estilizados** de los mercados
   financieros: colas pesadas, ausencia de autocorrelación lineal,
   agrupamiento de volatilidad y normalidad agregacional.
4. Comprobar que estos hechos aparecen en mercados muy distintos
   (universalidad).

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
estocásticos), S.B1–B3 y S.B5 (estimación, pruebas y permutaciones). Para
las notas formales: P.I1–I3, E.I1–I2 (autocorrelación y martingalas), C.I1 y
S.I4. Ruta avanzada: P.A1 (cumulantes), F.A2 (universalidad), S.A4.
"""

# ╔═╡ 7dfd0700-ca74-11f1-96fc-531560c2a69b
TableOfContents(title = "Contenido")

# ╔═╡ 7dfd0700-ca74-11f1-b10a-fba3278f137f
md"""
## 0. ¿Qué es la econofísica?

La economía tradicional suele partir de un agente representativo que optimiza.
La econofísica parte de otra pregunta: **¿qué patrones colectivos surgen cuando
muchos agentes heterogéneos interactúan?** Es la misma pregunta que hace la
física estadística con átomos y moléculas: nadie predice la trayectoria de una
molécula, pero la presión de un gas sí se predice.

Antes de construir modelos, un físico **mide**. Kepler describió las órbitas
antes de que Newton las explicara. En finanzas, las regularidades empíricas que
se repiten en casi todos los mercados y periodos se llaman **hechos
estilizados** (Cont, 2001). Hoy los vamos a medir. Durante el semestre
construiremos modelos que intenten explicarlos, y en la sesión 20 volveremos a
este cuaderno para ver cuáles reproducen.
"""

# ╔═╡ 7dfd070a-ca74-11f1-8b0f-0f4175a51ef0
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

# ╔═╡ 7dfd070a-ca74-11f1-ae46-3bb3ca927084
md"""
Retorno simple ``R``: $(@bind R_simple Slider(-0.5:0.01:0.5, default = 0.1, show_value = true))
"""

# ╔═╡ 7dfd071c-ca74-11f1-9661-c3d780a02c9b
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

# ╔═╡ 7dfd071c-ca74-11f1-b30e-adfe8f7e624e
md"""
## 2. Los datos

Usamos precios de cierre diarios de cuatro activos muy distintos: un índice bursátil mexicano, un tipo de cambio, un índice
estadounidense y una criptomoneda. Los datos se descargaron de Yahoo Finance
con `Code/descargar_datos.jl` y están en `Data/precios.csv`.
"""

# ╔═╡ 7dfd0728-ca74-11f1-bf9f-0541dc08ea37
datos = CSV.read(joinpath(@__DIR__, "..", "Data", "precios.csv"), DataFrame)

# ╔═╡ 7dfd0732-ca74-11f1-963e-3fe3d315cba6
md"""
**Activo:** $(@bind activo Select(unique(datos.activo)))
"""

# ╔═╡ 7dfd0732-ca74-11f1-8dff-dba393327886
begin
    serie = sort(filter(:activo => ==(activo), datos), :fecha)
    fechas = serie.fecha
    precios = serie.cierre
    r = diff(log.(precios))      # retornos logarítmicos diarios
    fechas_r = fechas[2:end]
    μ, σ = mean(r), std(r)
end;

# ╔═╡ 7dfd0732-ca74-11f1-8d6a-0919d9fc932c
md"""
$(@bind escala_log CheckBox(default = true)) Precio en escala logarítmica
"""

# ╔═╡ 7dfd0746-ca74-11f1-88cb-fde22e465535
md"""
Las líneas punteadas marcan tres episodios: la **crisis financiera** (octubre
de 2008), la **elección presidencial en EE. UU.** (9 de noviembre de 2016) y el
**inicio de la pandemia de COVID-19** (marzo de 2020).

Observa la gráfica de retornos: fluctúan alrededor de cero, pero **no con la
misma intensidad todo el tiempo**. Hay periodos de calma y periodos de
turbulencia. Esa observación visual es el punto de partida de todo lo que
sigue.
"""

# ╔═╡ 7dfd0746-ca74-11f1-9dce-13dd62311741
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

# ╔═╡ 7dfd074e-ca74-11f1-9fed-334ab50f3ded
z = (r .- μ) ./ σ;

# ╔═╡ 7dfd074e-ca74-11f1-9f05-8d3cd880408b
md"""
$(@bind colas_log CheckBox(default = false)) Eje vertical logarítmico
"""

# ╔═╡ 7dfd075a-ca74-11f1-ba1b-bfa7061d565e
md"""
**Curtosis en exceso de $activo:** $(round(kurtosis(r), digits = 1)) (la normal tiene 0).

¿Qué tan raros son los eventos extremos? La tabla compara cuántos días con
``|z| > k`` hubo en los datos contra cuántos predice la normal.
"""

# ╔═╡ 7dfd076e-ca74-11f1-b22f-f933425bdf48
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

# ╔═╡ 7dfd076e-ca74-11f1-80e7-57e251f15703
md"""
Rezagos: $(@bind L Slider(5:5:100, default = 30, show_value = true))
"""

# ╔═╡ 7dfd0778-ca74-11f1-a6b1-f934b0f8c420
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
"""

# ╔═╡ 7dfd0778-ca74-11f1-bdab-77d0f5e5c032
md"""
## 5. Hecho 3: agrupamiento de volatilidad

Ahora miramos el **tamaño** del movimiento, ``|r_t|``, sin importar si fue
subida o bajada.
"""

# ╔═╡ 7dfd0782-ca74-11f1-8009-e5cfa7d9a716
md"""
Los días grandes vienen **en racimos**: después de un día turbulento suele
venir otro día turbulento. Medimos esto con la autocorrelación de ``|r|``.

Rezagos: $(@bind L2 Slider(10:10:250, default = 100, show_value = true))
"""

# ╔═╡ 7dfd078c-ca74-11f1-8a06-25008333e294
md"""
La autocorrelación de ``r`` (azul) muere de inmediato. La de ``|r|`` (naranja)
es positiva y **decae muy lentamente**: sigue siendo significativa después de
semanas o meses. La dirección del mercado no tiene memoria; **la volatilidad
sí**. Esto es lo que modelan los procesos GARCH (sesión 5).

**[Avanzado]** Grafica la autocorrelación de ``|r|`` en escala log-log. ¿Decae
como una exponencial ``e^{-\tau/\tau_0}`` o como una ley de potencias
``\tau^{-\beta}``? Estima ``\beta``.
"""

# ╔═╡ 7dfd0796-ca74-11f1-9f47-4fa62b66ea9d
md"""
## 6. Hecho 4: normalidad agregacional

Los retornos de ``\Delta t`` días son sumas de ``\Delta t`` retornos diarios.
Si el teorema del límite central aplica, al agregar la distribución debería
acercarse a la normal. ¿Pasa eso?

Horizonte ``\Delta t`` (días): $(@bind Δt Slider([1, 2, 5, 10, 20, 60], default = 1, show_value = true))
"""

# ╔═╡ 7dfd07a0-ca74-11f1-bdd7-1f7e6d1b5813
md"""
Al agregar, la curtosis baja y la distribución se acerca a la normal: el
teorema del límite central termina ganando. Pero **gana más lento** de lo que
predice la línea punteada, que supone días independientes.

**[Avanzado]** Demuestra que si ``r_1, \dots, r_n`` son independientes e
idénticamente distribuidos con curtosis en exceso ``\kappa_1``, su suma tiene
curtosis ``\kappa_1/n``. ¿Qué hecho estilizado explica que los datos se alejen
de esa predicción?
"""

# ╔═╡ 7dfd07a0-ca74-11f1-a5a6-4bc71bb09f8d
md"""
## 7. ¿Qué tan especial es esto? Modelos nulos

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

# ╔═╡ 7dfd07aa-ca74-11f1-b139-f73e015e36be
r_modelo = if modelo == "Datos reales"
    r
elseif modelo == "Datos barajados"
    shuffle(Xoshiro(semilla), r)
else
    μ .+ σ .* randn(Xoshiro(semilla), length(r))
end;

# ╔═╡ 7dfd07b4-ca74-11f1-a8d6-bd85a73a5fe3
md"""
Compara los tres modelos y responde:
- ¿Cuál conserva las colas pesadas? ¿Cuál conserva el agrupamiento de
  volatilidad?
- Al barajar, la curtosis no cambia pero la autocorrelación de ``|r|``
  desaparece. ¿Qué nos dice eso sobre **dónde** vive cada hecho estilizado?

La función `panel_hechos` resume los cuatro hechos en una sola figura. La
usaremos en la sesión 20 para evaluar mercados artificiales.
"""

# ╔═╡ 7dfd07b4-ca74-11f1-a7e1-4971813eb64f
md"""
## 8. Universalidad

Hasta ahora miramos un activo a la vez. La tabla resume los mismos indicadores
para los cuatro.
"""

# ╔═╡ 7dfd07be-ca74-11f1-bf50-756340bea10f
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

# ╔═╡ 7dfd07be-ca74-11f1-869d-e9c28365a099
md"""
## Resumen

| Hecho estilizado | Cómo se mide | Lo que vimos |
|---|---|---|
| Colas pesadas | Densidad en escala log; curtosis; conteo de eventos extremos | Eventos de 5σ o más son comunes |
| Sin autocorrelación lineal | ACF de ``r`` | Pequeña en el rezago 1 e indistinguible de cero después |
| Agrupamiento de volatilidad | ACF de ``\lvert r\rvert`` | Positiva y persistente por meses |
| Normalidad agregacional | Curtosis vs. ``\Delta t`` | Baja al agregar, pero más lento que con días independientes |
"""

# ╔═╡ 7dfd07c8-ca74-11f1-b409-ebb2b700380c
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe la función `mis_retornos(p)` que reciba un vector de precios y
devuelva el vector de retornos logarítmicos. Reemplaza `missing` por tu código.
"""

# ╔═╡ 7dfd07c8-ca74-11f1-b032-cb289cacd26e
function mis_retornos(p)
    missing  # reemplaza esta línea
end

# ╔═╡ 7dfd07d2-ca74-11f1-83bd-b58d669c0597
md"""
**2.** Escribe la función `curtosis_exceso(x)` usando la fórmula de la sección 1.
"""

# ╔═╡ 7dfd07dc-ca74-11f1-a3f0-ebd826124a98
function curtosis_exceso(x)
    missing  # reemplaza esta línea
end

# ╔═╡ 7dfd07e6-ca74-11f1-8163-7fd13fc958aa
md"""
### Ruta básica: interpretación

3. Ordena los cuatro activos de colas más pesadas a menos pesadas. ¿Coincide
   con tu intuición sobre qué tan riesgoso es cada uno?
4. Según la tabla de eventos extremos, ¿cada cuánto debería ocurrir el peor día
   del S&P 500 si los retornos fueran normales? ¿Cuántos días con ``|z| > 5``
   hubo en realidad?
5. Si los retornos no tienen autocorrelación, ¿por qué ``|r|`` sí la tiene?
   ¿Qué implica para alguien que administra el riesgo de un portafolio?
6. Con el selector de modelos nulos: ¿qué hecho estilizado sobrevive al
   barajado y cuál no? Explica por qué.

### Ruta avanzada

7. Resuelve los dos ejercicios **[Avanzado]** de las secciones 5 y 6.
8. **Efecto apalancamiento.** Calcula ``\mathrm{corr}(r_t, |r_{t+\tau}|)`` para
   ``\tau = 1, \dots, 20`` en el S&P 500. ¿Es simétrica la relación entre
   caídas y volatilidad futura?

## Para la próxima sesión

- Lectura: Cont (2001), *Empirical properties of asset returns: stylized facts
  and statistical issues*, secciones 1–3.
- Lectura complementaria: Mantegna & Stanley, capítulos 1–2.
- En la sesión 3 mediremos **qué tan pesadas** son las colas: ¿siguen una ley
  de potencias?
"""

# ╔═╡ 7dfd07e6-ca74-11f1-8d59-eb2ed5ce2656
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 7dfd07f2-ca74-11f1-933f-7d8226cfc3a1
begin
    const AZUL    = "#2a78d6"   # serie principal
    const NARANJA = "#eb6834"   # segunda serie
    const GRIS    = "#52514e"   # referencias (normal, bandas, eventos)
end;

# ╔═╡ 7dfd071c-ca74-11f1-a0bf-6534c69c1f59
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

# ╔═╡ 7dfd078c-ca74-11f1-82dc-73224a064db2
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

# ╔═╡ 7dfd07f2-ca74-11f1-9395-1dcaf9958b68
const EVENTOS = [Date(2008, 10, 10), Date(2016, 11, 9), Date(2020, 3, 16)];

# ╔═╡ 7dfd07fa-ca74-11f1-b350-ade65ad9e00e
"Dibuja líneas verticales en los eventos que caen dentro del rango de `fechas`."
function marcar_eventos!(p, fechas)
    dentro = filter(e -> first(fechas) <= e <= last(fechas), EVENTOS)
    vline!(p, dentro, color = GRIS, ls = :dot, lw = 1, label = "")
end

# ╔═╡ 7dfd07fa-ca74-11f1-bed9-2708c8a279d5
"Retornos agregados en bloques no traslapados de `Δt` días."
agregar(r, Δt) = [sum(@view r[i:i+Δt-1]) for i in 1:Δt:length(r)-Δt+1]

# ╔═╡ 7dfd0796-ca74-11f1-bfc5-5382610f85b5
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

# ╔═╡ 7dfd082c-ca74-11f1-a8d6-61bb0de87719
"Densidad empírica en bins de ancho `ancho`; descarta bins vacíos (para escala log)."
function densidad_empirica(z; ancho = 0.25)
    h = fit(Histogram, z, -16:ancho:16)
    centros = collect(h.edges[1][1:end-1]) .+ ancho / 2
    dens = h.weights ./ (length(z) * ancho)
    keep = dens .> 0
    centros[keep], dens[keep]
end

# ╔═╡ 7dfd082c-ca74-11f1-a6f0-e5b86b3e89d1
dias_por_año(activo) = activo == "Bitcoin" ? 365 : 252

# ╔═╡ 7dfd082c-ca74-11f1-90e2-2f90cda53be6
"Límite simétrico del eje z: al menos ±8 y siempre mayor que el dato más extremo."
rango_z(z) = max(8, ceil(Int, maximum(abs, z)) + 1)

# ╔═╡ 7dfd075a-ca74-11f1-854b-c16a160905a5
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

# ╔═╡ 7dfd0796-ca74-11f1-b869-3fa3c409c3a2
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

# ╔═╡ 7dfd0836-ca74-11f1-ae10-931f55affce2
"Marcas del eje x: un año de cada `paso`, etiquetado sólo con el año."
function ticks_años(fechas; paso = 4)
    años = year(first(fechas)):paso:year(last(fechas))
    marcas = Date.(años, 1, 1)
    (Dates.value.(marcas), string.(años))
end

# ╔═╡ 7dfd0782-ca74-11f1-92ec-4f78dc335ff6
let
    p = plot(fechas_r, 100abs.(r), color = AZUL, lw = 0.6,
             ylabel = "|r| (%)", title = "$activo: tamaño de los movimientos diarios",
             xticks = ticks_años(fechas))
    marcar_eventos!(p, fechas)
end

# ╔═╡ 7dfd0836-ca74-11f1-b33a-f3682ca99ec5
"Marcas del eje de precios en escala log: 1, 2 y 5 por cada potencia de 10, con separador de miles."
function ticks_precio(precios)
    lo, hi = extrema(precios)
    marcas = [m * 10.0^k for k in floor(Int, log10(lo)):ceil(Int, log10(hi)) for m in (1, 2, 5)]
    marcas = filter(x -> lo <= x <= hi, marcas)
    etiqueta(x) = x < 1 ? string(round(x, sigdigits = 2)) :
                  replace(string(round(Int, x)), r"(?<=\d)(?=(\d{3})+$)" => ",")
    (marcas, etiqueta.(marcas))
end

# ╔═╡ 7dfd073c-ca74-11f1-aa3c-a3bc6adc73b6
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

# ╔═╡ 7dfd0840-ca74-11f1-8612-4deb867fe441
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

# ╔═╡ 7dfd0764-ca74-11f1-8549-d3162c44b971
extremo = let i = argmax(abs.(z))
    (fecha = fechas_r[i], r = round(100r[i], digits = 1),
     z = round(z[i], digits = 1),
     frecuencia = frecuencia_normal(abs(z[i]), dias_por_año(activo)))
end;

# ╔═╡ 7dfd0764-ca74-11f1-b921-2549e2ad7f8d
md"""
El movimiento más extremo de $activo fue el **día $(extremo.fecha)**:
``r = `` $(extremo.r)%, es decir, ``z = `` $(extremo.z).
Según la distribución normal, **un movimiento así ocurriría
$(extremo.frecuencia)**.
"""

# ╔═╡ 7dfd0840-ca74-11f1-83c9-b9e8a54d5c84
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

# ╔═╡ 7dfd075a-ca74-11f1-9824-737307a2c2b3
tabla_extremos(z, activo)

# ╔═╡ 7dfd084a-ca74-11f1-b262-fbb327a22395
function grafica_acf(x, L; titulo = "", color = AZUL)
    banda = 1.96 / sqrt(length(x))
    p = plot(1:L, fill(banda, L), fillrange = fill(-banda, L), color = GRIS,
             fillalpha = 0.15, lw = 0, xlabel = "Rezago τ (días)",
             ylabel = "ρ(τ)", title = titulo, ylims = (-0.15, 0.4))
    bar!(p, 1:L, autocor(x, 1:L), color = color, lw = 0, bar_width = 0.6)
    hline!(p, [0], color = :black, lw = 0.5)
end

# ╔═╡ 7dfd0778-ca74-11f1-b5e9-29971f771671
grafica_acf(r, L; titulo = "$activo: autocorrelación de los retornos")

# ╔═╡ 7dfd084a-ca74-11f1-a018-83ed79a78baa
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

# ╔═╡ 7dfd07aa-ca74-11f1-a177-81e9eacb557b
panel_hechos(r_modelo; titulo = "$activo: $(lowercase(modelo))")

# ╔═╡ 7dfd0852-ca74-11f1-9a87-df66de7fcca2
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

# ╔═╡ 7dfd07be-ca74-11f1-aa16-c50724b463d8
tabla_universalidad(datos)

# ╔═╡ 7dfd0852-ca74-11f1-b751-91ab79ec3332
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

# ╔═╡ 7dfd07d2-ca74-11f1-bbe4-2f8426204f91
verificar(() -> mis_retornos([100.0, 110.0, 99.0, 99.0]),
          [log(1.1), log(0.9), 0.0];
          pista = md"Usa `diff` y `log.`. Recuerda que ``r_t = \ln P_t - \ln P_{t-1}``.")

# ╔═╡ 7dfd085e-ca74-11f1-a422-01b498704a66
const VECTOR_PRUEBA = rand(Xoshiro(2026), TDist(5), 100_000);

# ╔═╡ 7dfd07dc-ca74-11f1-8c8d-39ee91579e97
verificar(() -> curtosis_exceso(VECTOR_PRUEBA), kurtosis(VECTOR_PRUEBA); rtol = 1e-2,
          pista = md"Estandariza con `(x .- mean(x)) ./ std(x)`, eleva a la cuarta potencia con `.^4`, promedia con `mean` y resta 3.")

# ╔═╡ Cell order:
# ╟─7dfd06ce-ca74-11f1-8c99-cb32088a2de7
# ╟─7dfd0700-ca74-11f1-96fc-531560c2a69b
# ╟─7dfd0700-ca74-11f1-b10a-fba3278f137f
# ╟─7dfd070a-ca74-11f1-8b0f-0f4175a51ef0
# ╟─7dfd070a-ca74-11f1-ae46-3bb3ca927084
# ╟─7dfd071c-ca74-11f1-a0bf-6534c69c1f59
# ╟─7dfd071c-ca74-11f1-9661-c3d780a02c9b
# ╟─7dfd071c-ca74-11f1-b30e-adfe8f7e624e
# ╠═7dfd0728-ca74-11f1-bf9f-0541dc08ea37
# ╟─7dfd0732-ca74-11f1-963e-3fe3d315cba6
# ╠═7dfd0732-ca74-11f1-8dff-dba393327886
# ╟─7dfd0732-ca74-11f1-8d6a-0919d9fc932c
# ╟─7dfd073c-ca74-11f1-aa3c-a3bc6adc73b6
# ╟─7dfd0746-ca74-11f1-88cb-fde22e465535
# ╟─7dfd0746-ca74-11f1-9dce-13dd62311741
# ╠═7dfd074e-ca74-11f1-9fed-334ab50f3ded
# ╟─7dfd074e-ca74-11f1-9f05-8d3cd880408b
# ╟─7dfd075a-ca74-11f1-854b-c16a160905a5
# ╟─7dfd075a-ca74-11f1-ba1b-bfa7061d565e
# ╟─7dfd075a-ca74-11f1-9824-737307a2c2b3
# ╟─7dfd0764-ca74-11f1-8549-d3162c44b971
# ╟─7dfd0764-ca74-11f1-b921-2549e2ad7f8d
# ╟─7dfd076e-ca74-11f1-b22f-f933425bdf48
# ╟─7dfd076e-ca74-11f1-80e7-57e251f15703
# ╟─7dfd0778-ca74-11f1-b5e9-29971f771671
# ╟─7dfd0778-ca74-11f1-a6b1-f934b0f8c420
# ╟─7dfd0778-ca74-11f1-bdab-77d0f5e5c032
# ╟─7dfd0782-ca74-11f1-92ec-4f78dc335ff6
# ╟─7dfd0782-ca74-11f1-8009-e5cfa7d9a716
# ╟─7dfd078c-ca74-11f1-82dc-73224a064db2
# ╟─7dfd078c-ca74-11f1-8a06-25008333e294
# ╟─7dfd0796-ca74-11f1-9f47-4fa62b66ea9d
# ╟─7dfd0796-ca74-11f1-b869-3fa3c409c3a2
# ╟─7dfd0796-ca74-11f1-bfc5-5382610f85b5
# ╟─7dfd07a0-ca74-11f1-bdd7-1f7e6d1b5813
# ╟─7dfd07a0-ca74-11f1-a5a6-4bc71bb09f8d
# ╠═7dfd07aa-ca74-11f1-b139-f73e015e36be
# ╟─7dfd07aa-ca74-11f1-a177-81e9eacb557b
# ╟─7dfd07b4-ca74-11f1-a8d6-bd85a73a5fe3
# ╟─7dfd07b4-ca74-11f1-a7e1-4971813eb64f
# ╟─7dfd07be-ca74-11f1-aa16-c50724b463d8
# ╟─7dfd07be-ca74-11f1-bf50-756340bea10f
# ╟─7dfd07be-ca74-11f1-869d-e9c28365a099
# ╟─7dfd07c8-ca74-11f1-b409-ebb2b700380c
# ╠═7dfd07c8-ca74-11f1-b032-cb289cacd26e
# ╟─7dfd07d2-ca74-11f1-bbe4-2f8426204f91
# ╟─7dfd07d2-ca74-11f1-83bd-b58d669c0597
# ╠═7dfd07dc-ca74-11f1-a3f0-ebd826124a98
# ╟─7dfd07dc-ca74-11f1-8c8d-39ee91579e97
# ╟─7dfd07e6-ca74-11f1-8163-7fd13fc958aa
# ╟─7dfd07e6-ca74-11f1-8d59-eb2ed5ce2656
# ╠═7dfd07e6-ca74-11f1-adea-b7ebfbd6c03e
# ╠═7dfd07f2-ca74-11f1-933f-7d8226cfc3a1
# ╠═7dfd07f2-ca74-11f1-9395-1dcaf9958b68
# ╠═7dfd07fa-ca74-11f1-b350-ade65ad9e00e
# ╠═7dfd07fa-ca74-11f1-bed9-2708c8a279d5
# ╠═7dfd082c-ca74-11f1-a8d6-61bb0de87719
# ╠═7dfd082c-ca74-11f1-a6f0-e5b86b3e89d1
# ╠═7dfd082c-ca74-11f1-90e2-2f90cda53be6
# ╠═7dfd0836-ca74-11f1-ae10-931f55affce2
# ╠═7dfd0836-ca74-11f1-b33a-f3682ca99ec5
# ╠═7dfd0840-ca74-11f1-8612-4deb867fe441
# ╠═7dfd0840-ca74-11f1-83c9-b9e8a54d5c84
# ╠═7dfd084a-ca74-11f1-b262-fbb327a22395
# ╠═7dfd084a-ca74-11f1-a018-83ed79a78baa
# ╠═7dfd0852-ca74-11f1-9a87-df66de7fcca2
# ╠═7dfd0852-ca74-11f1-b751-91ab79ec3332
# ╠═7dfd085e-ca74-11f1-a422-01b498704a66
