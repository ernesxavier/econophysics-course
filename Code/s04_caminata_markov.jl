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

# ╔═╡ f36fadf8-ca79-11f1-9a33-eda1611a8cef
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

# ╔═╡ f36fad12-ca79-11f1-9f66-cb7b4e1cda68
md"""
# Sesión 4: Caminata aleatoria y cadenas de Markov

**Objetivos de la sesión**
1. Simular una **caminata aleatoria** y comprobar que su dispersión crece como
   ``\sqrt{n}``: la firma de la difusión.
2. Construir un mercado que alterna entre **calma y turbulencia** con una
   **cadena de Markov** y medir cuánto duran sus crisis.
3. Entender por qué **una sola simulación no basta**: hay que promediar sobre
   muchas corridas.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s04_caminata_markov.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): E.B1–B3 (procesos
estocásticos, caminata aleatoria y promediar simulaciones), P.B2 (varianza de
una suma). Para las notas formales: E.I3 (cadenas de Markov), L.I1
(eigenvalores), S.I4 (datos dependientes). Ruta avanzada: E.A1 (movimiento
browniano), E.A3 (ecuación maestra), C.A3 (difusión), L.A5 (Perron-Frobenius).
"""

# ╔═╡ f36fad4c-ca79-11f1-831d-edcd1973fec5
TableOfContents(title = "Contenido")

# ╔═╡ f36fad62-ca79-11f1-a8bd-f33bd95b70a2
datos = CSV.read(joinpath(@__DIR__, "..", "Data", "precios.csv"), DataFrame);

# ╔═╡ f36fad58-ca79-11f1-9b1d-2d4004eba53f
md"""
## 0. El fenómeno: ¿qué proceso genera estas series?

En las sesiones 1 a 3 medimos tres propiedades de los retornos: casi no
tienen autocorrelación, su volatilidad se agrupa en el tiempo y sus colas
siguen una ley de potencias con ``\alpha \approx 3``. En esta sesión y la
siguiente cambiamos la pregunta: **¿qué mecanismo aleatorio produce series
así?**

Compara los retornos reales con retornos normales independientes de la misma
volatilidad. Los dos oscilan alrededor de cero, pero sólo uno tiene periodos
de calma y periodos de turbulencia: un mercado que alterna entre calma y crisis.

**Activo:** $(@bind activo Select(unique(datos.activo)))
"""

# ╔═╡ f36fad62-ca79-11f1-b4b7-8fcbedd3fad5
begin
    serie = sort(filter(:activo => ==(activo), datos), :fecha)
    r = diff(log.(serie.cierre))           # retornos logarítmicos diarios
    σ_r = std(r)
end;

# ╔═╡ f36fad7e-ca79-11f1-a1a3-19f15914136c
md"""
La estrategia de estas dos sesiones es la de un físico: proponer el modelo más
simple, simularlo, compararlo con los datos y, cuando falla, agregar sólo lo
necesario. Empezamos por el más simple de todos.

## 1. Repaso: la varianza de una suma

Si ``\xi_1, \dots, \xi_n`` son independientes,

```math
\mathbb{E}\Big[\sum_i \xi_i\Big] = \sum_i \mathbb{E}[\xi_i], \qquad
\mathrm{Var}\Big[\sum_i \xi_i\Big] = \sum_i \mathrm{Var}[\xi_i].
```

La media de una suma siempre es la suma de las medias; para la varianza hace
falta independencia (o al menos incorrelación). Con esta regla se entienden
la caminata aleatoria y, en la sección 4, el error de una simulación.
"""

# ╔═╡ f36fad7e-ca79-11f1-a924-4dbf60a0979a
md"""
## 2. La caminata aleatoria

En cada paso, un caminante sube 1 con probabilidad ``p`` o baja 1 con
probabilidad ``1 - p``. Su posición después de ``n`` pasos es
``S_n = \xi_1 + \dots + \xi_n``. Cada paso tiene media ``2p - 1`` y varianza
``4p(1-p)``, así que, por la sección 1,

```math
\mathbb{E}[S_n] = n(2p - 1), \qquad \mathrm{Var}[S_n] = 4np(1-p).
```

La media crece como ``n``, pero la dispersión crece como ``\sqrt{n}``: es la
firma de la **difusión**. Con ``p = 1/2`` es el modelo de Bachelier (1900)
para los precios.

Pasos: $(@bind n_pasos Slider(100:100:5000, default = 1000, show_value = true))
  Probabilidad de subir ``p``: $(@bind p_subir Slider(0.40:0.01:0.60, default = 0.50, show_value = true))
"""

# ╔═╡ f36fad94-ca79-11f1-9faa-29d4b523ace5
var_caminata = let rng = Xoshiro(4)
    finales = [sum(2 .* (rand(rng, n_pasos) .< p_subir) .- 1) for _ in 1:200]
    (simulada = var(finales), teorica = 4n_pasos * p_subir * (1 - p_subir))
end;

# ╔═╡ f36fad94-ca79-11f1-bc58-a1c3ab34096f
md"""
En el paso final, la varianza simulada de las 200 caminatas es
$(round(var_caminata.simulada, digits = 0)) y la teórica es
$(round(var_caminata.teorica, digits = 0)). Duplica el número de pasos: la
dispersión no se duplica, crece en un factor ``\sqrt{2} \approx 1.41``.
"""

# ╔═╡ f36fad9e-ca79-11f1-a3c9-67802fd7b237
md"""
## 3. Cadenas de Markov: un mercado con dos estados

Una **cadena de Markov** salta entre estados con probabilidades que sólo
dependen del estado actual. Imagina un mercado que puede estar en **calma**
(volatilidad baja) o en **turbulencia** (volatilidad tres veces mayor). Cada
día, desde la calma se pasa a turbulencia con probabilidad ``p``, y desde la
turbulencia se vuelve a la calma con probabilidad ``q``:

```math
P = \begin{pmatrix} 1-p & p \\ q & 1-q \end{pmatrix}.
```

A la larga, la fracción de días en turbulencia es ``\pi_T = p/(p+q)`` y un
episodio turbulento dura en promedio ``1/q`` días.

Probabilidad de entrar en turbulencia ``p``: $(@bind p_markov Slider(0.002:0.002:0.05, default = 0.01, show_value = true))
  de volver a la calma ``q``: $(@bind q_markov Slider(0.01:0.01:0.2, default = 0.05, show_value = true))
"""

# ╔═╡ f36fadb0-ca79-11f1-a0f6-11c0934039ae
md"""
## 4. Promediar simulaciones: una corrida no basta

La teoría dice que el mercado pasa una fracción ``\pi_T = p/(p+q)`` de los días
en turbulencia. Pero si sólo simulamos **una** vez, ¿qué tan cerca queda la
fracción medida? Repetimos la simulación de 3000 días 300 veces, con los
mismos ``p`` y ``q`` de la sección 3, y miramos cuánto varía el resultado.
"""

# ╔═╡ f36fadc6-ca79-11f1-8007-f5a28d236162
md"""
## Resumen

| Proceso | Regla | Qué reproduce | Qué no |
|---|---|---|---|
| Caminata aleatoria | ``S_{n+1} = S_n \pm 1`` | Difusión: dispersión ``\propto \sqrt{n}`` | Colas, agrupamiento |
| Cadena de Markov (dos regímenes) | Volatilidad salta entre dos niveles | Agrupamiento de corto plazo, colas algo pesadas | Memoria larga |

**Lección de método:** una simulación es una muestra; para estimar algo hay
que promediar muchas corridas.
"""

# ╔═╡ f36fadc6-ca79-11f1-b906-4bc46e7e2fc5
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `var_caminata_teo(n, p)`, la varianza teórica de la posición de una
caminata aleatoria después de `n` pasos con probabilidad `p` de subir.
"""

# ╔═╡ f36fadd0-ca79-11f1-84da-e73a0f79c279
function var_caminata_teo(n, p)
    missing  # reemplaza esta línea
end

# ╔═╡ f36fadda-ca79-11f1-a6f5-ff638356c1a4
md"""
**2.** Escribe `estacionaria_dos_estados(p, q)`, que devuelva el vector
``(\pi_C, \pi_T)`` con la fracción de días en calma y en turbulencia a la
larga.
"""

# ╔═╡ f36fadda-ca79-11f1-9554-b131270117e5
function estacionaria_dos_estados(p, q)
    missing  # reemplaza esta línea
end

# ╔═╡ f36fade4-ca79-11f1-8e30-65cc6da7e602
md"""
### Ruta básica: interpretación

3. En la caminata aleatoria con ``p = 1/2``, ¿cuánto hay que aumentar el
   número de pasos para que la dispersión típica se triplique?
4. En el modelo de dos regímenes con ``p = 0.02`` y ``q = 0.1``, ¿qué fracción
   del tiempo pasa el mercado en turbulencia y cuánto dura un episodio
   turbulento? ¿Por qué este modelo produce agrupamiento de volatilidad?
5. En la sección 4, ¿qué le pasa al ancho del histograma cuando las crisis son
   más largas (``q`` más pequeño)? ¿Por qué 3000 días de un mercado persistente
   "valen" menos que 3000 días independientes?
6. El modelo de dos regímenes reproduce el agrupamiento de volatilidad, pero su
   autocorrelación de ``|r|`` decae exponencialmente. ¿Qué hecho de la sesión 2
   no alcanza a reproducir?

### Ruta avanzada

7. Demuestra que ``\pi = (q/(p+q),\ p/(p+q))`` es un eigenvector izquierdo de
   la matriz ``P`` con eigenvalor 1 y que el otro eigenvalor es ``1 - p - q``.
8. Para la fracción de días turbulentos medida en una corrida de ``n`` días,
   las notas formales dan una varianza
   ``\pi_C\pi_T(1+\lambda)/\big((1-\lambda)n\big)``. Compárala con la
   varianza de las 300 corridas de la sección 4.

## Para la próxima sesión

- Lectura: Slanina (2014), sección 2.1. Complementaria: Requião da Cunha
  (2021), sección 3.1.
- En la sesión 5 agregamos una **tendencia** al ruido: ecuación de Langevin,
  movimiento browniano geométrico y GARCH, un modelo en el que la volatilidad
  tiene memoria.
"""

# ╔═╡ f36fadee-ca79-11f1-a6ad-353ff9cee68c
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ f36fadf8-ca79-11f1-908b-07feb7e53ea8
begin
    const AZUL    = "#2a78d6"   # serie principal
    const NARANJA = "#eb6834"   # segunda serie o modelo
    const GRIS    = "#52514e"   # referencias
end;

# ╔═╡ f36fad76-ca79-11f1-9845-cb449be8418e
let
    n = min(length(r), 2000)
    real = 100 .* r[end-n+1:end]
    iid = 100σ_r .* randn(Xoshiro(3), n)
    lim = maximum(abs, real) * 1.05
    p1 = plot(real; color = AZUL, lw = 0.6, ylims = (-lim, lim), ylabel = "r (%)",
              title = "$activo: últimos $n días")
    p2 = plot(iid; color = GRIS, lw = 0.6, ylims = (-lim, lim), ylabel = "r (%)",
              xlabel = "día", title = "Normales independientes con la misma volatilidad")
    plot(p1, p2; layout = (2, 1), size = (720, 460))
end

# ╔═╡ f36fad8a-ca79-11f1-84cc-43f1d64519fe
let
    rng = Xoshiro(4)
    S = [cumsum(rand(rng, n_pasos) .< p_subir) .* 2 .- (1:n_pasos) for _ in 1:200]
    t = 0:n_pasos
    μ = t .* (2p_subir - 1); σ = sqrt.(4 .* t .* p_subir .* (1 - p_subir))
    p = plot(xlabel = "paso n", ylabel = "posición Sₙ", legend = :topleft,
             title = "200 caminatas aleatorias con p = $(p_subir)")
    for (i, s) in enumerate(S)
        plot!(p, t, [0; s]; color = AZUL, alpha = 0.12, lw = 0.8, label = i == 1 ? "Trayectorias" : "")
    end
    plot!(p, t, μ; color = NARANJA, lw = 2, label = "Media teórica")
    plot!(p, t, μ .+ 2σ; color = NARANJA, lw = 1.5, ls = :dash, label = "Media ± 2 desviaciones")
    plot!(p, t, μ .- 2σ; color = NARANJA, lw = 1.5, ls = :dash, label = "")
end

# ╔═╡ f36fae02-ca79-11f1-81e0-0786cfb62b91
"Mercado con dos regímenes de volatilidad que siguen una cadena de Markov."
function simular_regimenes(p, q; n = 3000, σ_calma = 0.008, factor = 3, rng = Xoshiro(8))
    estado = ones(Int, n)
    for t in 2:n
        u = rand(rng)
        estado[t] = estado[t-1] == 1 ? (u < p ? 2 : 1) : (u < q ? 1 : 2)
    end
    σ = ifelse.(estado .== 1, σ_calma, factor * σ_calma)
    (estado = estado, r = σ .* randn(rng, n))
end

# ╔═╡ f36fad9e-ca79-11f1-b0e2-1b93f7ebceb3
mercado = simular_regimenes(p_markov, q_markov; n = 3000, σ_calma = 0.008, factor = 3);

# ╔═╡ f36fada8-ca79-11f1-99d0-b5c9660a285f
md"""
Fracción de días en turbulencia: $(round(100 * mean(mercado.estado .== 2), digits = 1))%
en la simulación contra $(round(100p_markov / (p_markov + q_markov), digits = 1))% de la
teoría (zonas sombreadas).

Este modelo tan simple ya reproduce dos hechos estilizados: los retornos no
tienen autocorrelación, pero ``|r|`` sí, porque un día turbulento suele seguir
a otro. Y la mezcla de dos normales con varianzas distintas tiene colas más
pesadas que una normal (curtosis en exceso de
$(round(kurtosis(mercado.r), digits = 1))). Pero la autocorrelación de ``|r|``
decae como ``(1 - p - q)^\tau``: exponencialmente, más rápido que en los datos.
"""

# ╔═╡ f36fadb0-ca79-11f1-861c-cdf2cb362127
corridas = let
    fr = [mean(simular_regimenes(p_markov, q_markov; n = 3000, σ_calma = 0.008, factor = 3,
                                 rng = Xoshiro(1000 + k)).estado .== 2) for k in 1:300]
    (fracciones = fr, teoria = p_markov / (p_markov + q_markov))
end;

# ╔═╡ f36fadbc-ca79-11f1-b9d4-350e74469178
let c = corridas
    p = histogram(100 .* c.fracciones; bins = 30, color = AZUL, alpha = 0.6, lw = 0,
                  xlabel = "días en turbulencia (%) en una corrida de 3000 días",
                  ylabel = "número de corridas", label = "300 corridas", legend = :topright,
                  title = "Lo que mide una sola simulación")
    vline!(p, [100c.teoria]; color = NARANJA, lw = 2.5, label = "Teoría p/(p+q)")
    vline!(p, [100mean(c.fracciones)]; color = GRIS, lw = 2, ls = :dash, label = "Promedio de las 300")
end

# ╔═╡ f36fadbc-ca79-11f1-a4e9-1bfe9e9187fe
md"""
Una corrida cualquiera puede quedar lejos de la teoría: en estas 300, la
fracción de días turbulentos va de
$(round(100minimum(corridas.fracciones), digits = 1))% a
$(round(100maximum(corridas.fracciones), digits = 1))%, contra
$(round(100corridas.teoria, digits = 1))% de la teoría. El **promedio** de
las corridas, en cambio, queda muy cerca: $(round(100mean(corridas.fracciones), digits = 1))%.

La razón es la persistencia: como un día turbulento suele seguir a otro, 3000
días no son 3000 observaciones independientes, sino muchas menos. Haz ``q``
más pequeño (crisis más largas) y verás que el histograma se ensancha. Por eso,
en todo el curso, cuando un resultado salga de una simulación, lo
**promediaremos sobre muchas corridas** y reportaremos cuánto varía.
"""

# ╔═╡ f36fae02-ca79-11f1-825c-b3ff8809f0f2
"Agrupa índices consecutivos en pares (inicio, fin) para sombrear intervalos."
function rangos(idx)
    out = Float64[]
    inicio = idx[1]
    for i in 2:length(idx)
        if idx[i] != idx[i-1] + 1
            append!(out, [inicio - 0.5, idx[i-1] + 0.5]); inicio = idx[i]
        end
    end
    append!(out, [inicio - 0.5, idx[end] + 0.5])
end

# ╔═╡ f36fada8-ca79-11f1-b0b8-cbc0c04c78de
let
    m = mercado
    p1 = plot(100 .* m.r; color = AZUL, lw = 0.6, ylabel = "r (%)",
              title = "Retornos simulados con dos regímenes")
    turb = findall(m.estado .== 2)
    isempty(turb) || vspan!(p1, rangos(turb); color = NARANJA, alpha = 0.15, lw = 0)
    p2 = plot(1:60, autocor(abs.(m.r), 1:60); color = AZUL, lw = 2, xlabel = "rezago τ (días)",
              ylabel = "ρ(τ)", title = "Autocorrelación de |r|", ylims = (-0.05, 0.5))
    plot!(p2, 1:60, autocor(m.r, 1:60); color = GRIS, lw = 1.5)
    hline!(p2, [0]; color = :black, lw = 0.5)
    annotate!(p2, 30, 0.42, text("azul: |r| · gris: r", 8, GRIS))
    plot(p1, p2; layout = (2, 1), size = (720, 480))
end

# ╔═╡ f36fae0c-ca79-11f1-8434-e7cf3bbe014d
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

# ╔═╡ f36fadd0-ca79-11f1-aece-4b27575f96c9
verificar(() -> [var_caminata_teo(100, 0.5), var_caminata_teo(1000, 0.6)], [100.0, 960.0];
          pista = md"Cada paso tiene varianza ``4p(1-p)`` y los pasos son independientes.")

# ╔═╡ f36fade4-ca79-11f1-8708-d58e22d50319
verificar(() -> estacionaria_dos_estados(0.01, 0.05), [5 / 6, 1 / 6];
          pista = md"En el estado estacionario entra a la turbulencia tanta probabilidad como sale: ``\pi_C p = \pi_T q``, con ``\pi_C + \pi_T = 1``.")

# ╔═╡ Cell order:
# ╟─f36fad12-ca79-11f1-9f66-cb7b4e1cda68
# ╟─f36fad4c-ca79-11f1-831d-edcd1973fec5
# ╟─f36fad58-ca79-11f1-9b1d-2d4004eba53f
# ╠═f36fad62-ca79-11f1-a8bd-f33bd95b70a2
# ╠═f36fad62-ca79-11f1-b4b7-8fcbedd3fad5
# ╟─f36fad76-ca79-11f1-9845-cb449be8418e
# ╟─f36fad7e-ca79-11f1-a1a3-19f15914136c
# ╟─f36fad7e-ca79-11f1-a924-4dbf60a0979a
# ╟─f36fad8a-ca79-11f1-84cc-43f1d64519fe
# ╟─f36fad94-ca79-11f1-bc58-a1c3ab34096f
# ╟─f36fad94-ca79-11f1-9faa-29d4b523ace5
# ╟─f36fad9e-ca79-11f1-a3c9-67802fd7b237
# ╟─f36fad9e-ca79-11f1-b0e2-1b93f7ebceb3
# ╟─f36fada8-ca79-11f1-b0b8-cbc0c04c78de
# ╟─f36fada8-ca79-11f1-99d0-b5c9660a285f
# ╟─f36fadb0-ca79-11f1-a0f6-11c0934039ae
# ╟─f36fadb0-ca79-11f1-861c-cdf2cb362127
# ╟─f36fadbc-ca79-11f1-b9d4-350e74469178
# ╟─f36fadbc-ca79-11f1-a4e9-1bfe9e9187fe
# ╟─f36fadc6-ca79-11f1-8007-f5a28d236162
# ╟─f36fadc6-ca79-11f1-b906-4bc46e7e2fc5
# ╠═f36fadd0-ca79-11f1-84da-e73a0f79c279
# ╟─f36fadd0-ca79-11f1-aece-4b27575f96c9
# ╟─f36fadda-ca79-11f1-a6f5-ff638356c1a4
# ╠═f36fadda-ca79-11f1-9554-b131270117e5
# ╟─f36fade4-ca79-11f1-8708-d58e22d50319
# ╟─f36fade4-ca79-11f1-8e30-65cc6da7e602
# ╟─f36fadee-ca79-11f1-a6ad-353ff9cee68c
# ╠═f36fadf8-ca79-11f1-9a33-eda1611a8cef
# ╠═f36fadf8-ca79-11f1-908b-07feb7e53ea8
# ╠═f36fae02-ca79-11f1-81e0-0786cfb62b91
# ╠═f36fae02-ca79-11f1-825c-b3ff8809f0f2
# ╠═f36fae0c-ca79-11f1-8434-e7cf3bbe014d
