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

# ╔═╡ 098ea8ce-ca7e-11f1-9fdd-e97268d7a5ad
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

# ╔═╡ 098ea810-ca7e-11f1-8ee4-11801deff225
md"""
# Sesión 16: Schelling: segregación

**Objetivos de la sesión**
1. Simular el **modelo de segregación de Schelling** desde cero en Julia y
   comprobar cómo preferencias individuales leves producen segregación
   residencial fuerte.
2. Entender el mecanismo de **externalidad de vecindario** y por qué el umbral
   crítico separa tres fases: mezcla, segregación y frustración.
3. Analizar la **analogía con el modelo de Ising** a temperatura baja con
   magnetización conservada y sus límites.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s16_schelling.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): F.B4 (transiciones de fase)
y F.I3 (modelo de Ising).
"""

# ╔═╡ 098ea842-ca7e-11f1-98cc-45514a1bfd32
TableOfContents(title = "Contenido")

# ╔═╡ 098ea84c-ca7e-11f1-a5ae-6f90cd7e2ac6
md"""
## 0. El fenómeno: barrios segregados sin intolerancia

En una ciudad con dos grupos (azul y naranja), ningún habitante pide vivir
rodeado únicamente de su propio grupo: cada uno se conforma con que al menos
**tres de sus ocho vecinos** (37.5%) sean como él. Es decir, todos aceptan ser
minoría en su cuadra (hasta 62.5% del otro grupo). Sin embargo, cuando los
insatisfechos se mudan a viviendas vacías, la ciudad termina partida en barrios
casi homogéneos (Schelling, 1971).
"""

# ╔═╡ 098ea860-ca7e-11f1-9b64-d9cf51210681
md"""
## 1. El modelo de segregación de Schelling

En los modelos de votante y de Kirman (sesión 15) los agentes están fijos en
una red y cambian de opinión. En el modelo de **Schelling** (1971) ocurre lo
contrario: cada agente conserva su identidad (azul o naranja) pero **cambia de
lugar** si no está conforme con su vecindario.

**Reglas.**
1. En una cuadrícula de ``L \times L`` con fronteras periódicas, una fracción
   ``\rho_0`` de las celdas está vacía (blanco) y el resto se reparte en partes
   iguales entre dos grupos (azul y naranja).
2. Cada agente mira sus 8 celdas adyacentes (**vecindad de Moore**). Si tiene
   ``n_i`` vecinos ocupados y ``n_i^{=}`` son de su mismo grupo, su
   **similitud local** es
   ```math
   s_i = \frac{n_i^{=}}{n_i} .
   ```
   Si ``s_i \ge \tau`` (el **umbral de tolerancia**), el agente está
   **satisfecho**; si ``s_i < \tau``, está **insatisfecho** (si ``n_i = 0``, se
   considera satisfecho).
3. En cada ronda, los agentes insatisfechos se mudan, en orden aleatorio, a
   celdas vacías elegidas al azar. El proceso se detiene cuando todos están
   satisfechos o al llegar al máximo de rondas.

**¿Por qué una mudanza desencadena otras?** Cuando un agente azul se muda
porque tenía pocos vecinos azules, la fracción de azules **baja** para los
vecinos azules que se quedaron atrás y **sube** en el barrio al que llega. Un
vecino que antes estaba justo en el umbral (por ejemplo, con 3 de 8) pasa a 2
de 7 (``28.6\% < 37.5\%``), se vuelve insatisfecho y también se muda: una
**cascada de mudanzas**.
"""

# ╔═╡ 098ea86a-ca7e-11f1-a122-cb6b7cf72c41
md"""
## 2. Laboratorio: de las preferencias individuales a la ciudad segregada

Mueve el umbral ``\tau`` —la fracción mínima de vecinos propios que pide cada
agente para no mudarse— y la fracción ``\rho_0`` de viviendas vacías.

Umbral requerido ``\tau``: $(@bind τ_sch Slider([0.125, 0.25, 0.375, 0.50, 0.625, 0.75, 0.85], default = 0.375, show_value = true))
  Fracción de celdas vacías ``\rho_0``: $(@bind vac_sch Slider([0.05, 0.15, 0.30], default = 0.15, show_value = true))
"""

# ╔═╡ 098ea874-ca7e-11f1-bc3b-43f7ce624683
md"""
### 2.1 Barrido del umbral de tolerancia: las tres fases

Para ver la transición completa (repaso F.B4), medimos la **similitud promedio
final** y la **fracción de agentes satisfechos** al variar el umbral ``\tau``
de 10% a 85%:
"""

# ╔═╡ 098ea892-ca7e-11f1-8112-fbd603c1be12
md"""
## 3. Analogía con el modelo de Ising y sus límites

Si asignamos ``s_i = +1`` a un grupo, ``s_i = -1`` al otro y ``s_i = 0`` a las
celdas vacías, la energía del modelo de Ising (sesión 13, repaso F.I3)

```math
H = -J \sum_{\langle i,j\rangle} s_i s_j
```

baja cuando hay más pares de vecinos del mismo grupo y menos **fronteras**
(enlaces entre grupos distintos). En el modelo de Schelling el número total de
agentes de cada grupo no cambia —como en un Ising con **magnetización
conservada**— y las mudanzas tienden a acortar el perímetro de las fronteras
entre los dos grupos, igual que las gotas de aceite en agua minimizan su
tensión superficial (Vinković y Kirman, 2006).

Comprobémoslo midiendo en cada ronda de la simulación la fracción de enlaces
entre vecinos ocupados que conectan a dos grupos distintos (**frontera
interfacial**):
"""

# ╔═╡ 098ea89c-ca7e-11f1-8569-95ef9e0fa2af
md"""
**¿Dónde termina la analogía?** En el algoritmo de Metropolis (sesión 13), cada
cambio se acepta según el cambio en la energía **global** ``\Delta H``. En
Schelling, en cambio:

- A cada agente sólo le importa pasar su propio umbral ``s_i \ge \tau``: una
  vez que tiene 3 de 8 vecinos iguales (con ``\tau = 37.5\%``), le da igual
  tener 3 que 8.
- Cuando un agente insatisfecho se muda al azar a una celda vacía, su llegada
  puede incomodar a más vecinos nuevos de los que dejó atrás y **subir**
  momentáneamente la frontera total (como se ve cuando ``\tau`` es muy alto y
  el sistema no logra ordenarse).
"""

# ╔═╡ 098ea8a6-ca7e-11f1-84d7-73f78f0de57b
md"""
## Resumen

| Concepto | Idea |
|---|---|
| Similitud local | Fracción ``s_i = n_i^{=}/n_i`` de vecinos ocupados que pertenecen al mismo grupo |
| Regla de Schelling | Si ``s_i < \tau``, el agente se muda a una celda vacía al azar |
| Externalidad de vecindario | Mudarse reduce ``s_j`` para los antiguos vecinos del mismo grupo y desencadena cascadas |
| Tres fases | ``\tau \le 25\%``: mezcla; ``30\% \le \tau \le 70\%``: segregación fuerte (``\approx 80\%``); ``\tau \ge 75\%``: frustración |
| Analogía con Ising | Equivale aproximadamente a enfriar un Ising con número fijo de espines ``\pm 1``, reduciendo las fronteras entre dominios |
"""

# ╔═╡ 098ea8a6-ca7e-11f1-8442-27a78f706b0f
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `similitud_local(iguales, totales)`, que devuelva la fracción
`iguales / totales` cuando `totales > 0`, y `1.0` cuando `totales == 0`.
"""

# ╔═╡ 098ea8b0-ca7e-11f1-beb8-b9027551565e
function similitud_local(iguales, totales)
    missing  # reemplaza esta línea
end

# ╔═╡ 098ea8ba-ca7e-11f1-8911-f90128439f0e
md"""
**2.** Escribe `satisfecho(iguales, totales, τ)`, que devuelva `true` si un
agente de Schelling con `iguales` vecinos de su mismo grupo, de un total de
`totales` vecinos ocupados, cumple el umbral `τ` (si `totales == 0`, devuelve
`true`).
"""

# ╔═╡ 098ea8ba-ca7e-11f1-b041-09b6566e3592
function satisfecho(iguales, totales, τ)
    missing  # reemplaza esta línea
end

# ╔═╡ 098ea8c4-ca7e-11f1-94bc-33cf6a71df93
md"""
### Ruta básica: interpretación

3. En el modelo de Schelling con ``\tau = 0.375``, ningún agente rechaza ser
   minoría en su barrio (acepta hasta 62.5% de vecinos del otro grupo). Explica
   paso a paso por qué la ciudad no se queda con barrios mezclados al 50%.
4. ¿Qué pasa en el modelo de Schelling cuando el umbral ``\tau`` sube por
   encima de 0.75? ¿Por qué baja la fracción de agentes satisfechos?
5. ¿En qué se parece y en qué se distingue la regla de mudanza de Schelling del
   algoritmo de Metropolis para el modelo de Ising (sesión 13)?

### Ruta avanzada

6. En una cuadrícula sin celdas vacías donde dos agentes de grupos distintos
   intercambian sus casas sólo si ambos aumentan estrictamente su número de
   vecinos iguales, demuestra que el número total de enlaces entre grupos
   distintos disminuye en cada intercambio y que la dinámica no puede ciclar.
7. Modifica `simular_schelling` para que los agentes insatisfechos sólo se
   muden a celdas vacías donde sí queden satisfechos (si existe alguna). ¿Qué
   cambia en el régimen ``\tau \ge 0.75``?

## Para la próxima sesión

- Lectura: Kirman (2011), capítulo 6. Complementaria: Sen y Chakrabarti (2014),
  capítulo 4; Schelling (1971); de consulta, Wilensky y Rand (2015). Para la
  ruta avanzada, Slanina (2014), sección 8.2; Sen y Chakrabarti (2014),
  capítulo 6; Vinković y Kirman (2006).
- En la **sesión 17** (`Code/s17_replicador_goodwin.jl`) los agentes elegirán
  estrategias según sus pagos: veremos la **dinámica del replicador**, el juego
  **Halcón-Paloma** y el **ciclo macroeconómico de Goodwin**.
"""

# ╔═╡ 098ea8ce-ca7e-11f1-a706-a716ce70bf9a
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 098ea8d8-ca7e-11f1-adb2-b75c2a497a01
begin
    const AZUL       = "#2a78d6"   # serie principal, grupo 1
    const NARANJA    = "#eb6834"   # segunda serie o teoría, grupo 2
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const BLANCO     = "#ffffff"
end;

# ╔═╡ 098ea8d8-ca7e-11f1-9de0-35f2992510ec
"Porcentaje con a lo más un decimal: 0.375 → \"37.5%\", 0.25 → \"25%\"."
porcentaje(x) = replace(string(round(100x, digits = 1)), r"\.0$" => "") * "%"

# ╔═╡ 098ea8e2-ca7e-11f1-b06a-49b2737fbc46
"Cuenta vecinos del mismo tipo y vecinos ocupados totales en la vecindad de Moore (8 vecinos, periódica)."
function vecindad_schelling(tablero, r, c)
    L = size(tablero, 1)
    tipo = tablero[r, c]
    iguales, totales = 0, 0
    for dr in -1:1, dc in -1:1
        (dr == 0 && dc == 0) && continue
        v = tablero[mod1(r + dr, L), mod1(c + dc, L)]
        if v != 0
            totales += 1
            v == tipo && (iguales += 1)
        end
    end
    iguales, totales
end

# ╔═╡ 098ea8e2-ca7e-11f1-82a4-837d0b1404aa
"Mide la similitud promedio, la fracción de agentes satisfechos y la fracción de enlaces entre grupos distintos."
function metricas_schelling(tablero, τ)
    L = size(tablero, 1)
    suma_sim, con_vecinos, felices, agentes = 0.0, 0, 0, 0
    enlaces_distintos, enlaces_ocupados = 0, 0
    for r in 1:L, c in 1:L
        tablero[r, c] == 0 && continue
        agentes += 1
        ig, tot = vecindad_schelling(tablero, r, c)
        enlaces_ocupados += tot
        enlaces_distintos += (tot - ig)
        if tot > 0
            suma_sim += ig / tot
            con_vecinos += 1
            ig / tot >= τ && (felices += 1)
        else
            felices += 1
        end
    end
    (similitud = con_vecinos > 0 ? suma_sim / con_vecinos : 1.0,
     felices = felices / max(agentes, 1),
     frontera = enlaces_ocupados > 0 ? enlaces_distintos / enlaces_ocupados : 0.0)
end

# ╔═╡ 098ea8ec-ca7e-11f1-8d00-d5ae5a76c5eb
"""
Simula el modelo de segregación de Schelling en una cuadrícula L × L con
fracción `p_vacio` de celdas vacías y umbral de tolerancia `τ`.
"""
function simular_schelling(L, p_vacio, τ; pasos_max = 120, rng = Xoshiro(21))
    n_celdas = L * L
    n_vacios = round(Int, n_celdas * p_vacio)
    n_agentes = n_celdas - n_vacios
    n_uno = n_agentes ÷ 2
    celdas = vcat(zeros(Int, n_vacios), ones(Int, n_uno), fill(2, n_agentes - n_uno))
    shuffle!(rng, celdas)
    tablero = reshape(celdas, L, L)

    m0 = metricas_schelling(tablero, τ)
    sim_hist = [m0.similitud]
    fel_hist = [m0.felices]
    fro_hist = [m0.frontera]

    for _ in 1:pasos_max
        infelices = Tuple{Int,Int}[]
        vacias = Tuple{Int,Int}[]
        for r in 1:L, c in 1:L
            if tablero[r, c] == 0
                push!(vacias, (r, c))
            else
                ig, tot = vecindad_schelling(tablero, r, c)
                (tot > 0 && ig / tot < τ) && push!(infelices, (r, c))
            end
        end
        isempty(infelices) && break
        shuffle!(rng, infelices)
        shuffle!(rng, vacias)
        for (r, c) in infelices
            isempty(vacias) && break
            vr, vc = popfirst!(vacias)
            tablero[vr, vc] = tablero[r, c]
            tablero[r, c] = 0
            push!(vacias, (r, c))
        end
        mt = metricas_schelling(tablero, τ)
        push!(sim_hist, mt.similitud)
        push!(fel_hist, mt.felices)
        push!(fro_hist, mt.frontera)
    end
    (tablero = tablero, similitud = sim_hist, felices = fel_hist, frontera = fro_hist)
end

# ╔═╡ 098ea874-ca7e-11f1-9b0b-05215f7d8e17
sim_schelling = simular_schelling(40, vac_sch, τ_sch; pasos_max = 120, rng = Xoshiro(21));

# ╔═╡ 098ea874-ca7e-11f1-a9c1-1b5ec4c6e520
let
    s = sim_schelling
    paleta = cgrad([BLANCO, AZUL, NARANJA])
    p1 = heatmap(s.tablero; color = paleta, clims = (0, 2), aspect_ratio = :equal,
                 colorbar = false, axis = false, ticks = false, grid = false,
                 title = "Ciudad final (τ = $(porcentaje(τ_sch))): similitud = $(round(Int, 100s.similitud[end]))%")
    rondas = 0:length(s.similitud)-1
    p2 = plot(rondas, 100 .* s.similitud; color = AZUL, lw = 2.2, ylims = (0, 105),
              xlabel = "ronda", ylabel = "porcentaje (%)", legend = :bottomright,
              label = "Similitud promedio de los vecindarios",
              title = "Evolución en $(length(rondas)-1) rondas")
    plot!(p2, rondas, 100 .* s.felices; color = NARANJA, lw = 2, ls = :dash,
          label = "Agentes satisfechos")
    hline!(p2, [100τ_sch]; color = GRIS, ls = :dot, label = "Umbral pedido τ = $(porcentaje(τ_sch))")
    plot(p1, p2; layout = Plots.grid(1, 2; widths = [0.44, 0.56]), size = (900, 380))
end

# ╔═╡ 098ea892-ca7e-11f1-be80-f7a9bc70dcc6
let
    s = sim_schelling
    rondas = 0:length(s.frontera)-1
    p = plot(rondas, 100 .* s.frontera; color = NARANJA, lw = 2.2,
             xlabel = "ronda", ylabel = "enlaces entre grupos distintos (%)",
             ylims = (0, 60), legend = :topright,
             label = "Frontera entre grupos (∝ energía de interfaz)",
             title = "Las mudanzas reducen el perímetro entre grupos (τ = $(porcentaje(τ_sch)))")
    hline!(p, [50]; color = GRIS_CLARO, ls = :dot, label = "Mezcla al azar (50%)")
end

# ╔═╡ 098ea8f8-ca7e-11f1-8bac-23e0b8c490e3
"Ejemplo fijo de Schelling para la sección 0 (antes y después con τ = 0.375)."
ejemplo_schelling = let
    ini = simular_schelling(40, 0.15, 0.375; pasos_max = 0, rng = Xoshiro(21))
    fin = simular_schelling(40, 0.15, 0.375; pasos_max = 120, rng = Xoshiro(21))
    (ini = (tablero = ini.tablero, similitud = ini.similitud[1], felices = ini.felices[1]),
     fin = (tablero = fin.tablero, similitud = fin.similitud[end], felices = fin.felices[end]))
end;

# ╔═╡ 098ea854-ca7e-11f1-9322-97fd00adc5e0
let
    ini, fin = ejemplo_schelling
    paleta = cgrad([BLANCO, AZUL, NARANJA])
    p1 = heatmap(ini.tablero; color = paleta, clims = (0, 2), aspect_ratio = :equal,
                 colorbar = false, axis = false, ticks = false, grid = false,
                 title = "Inicio al azar (similitud = $(round(Int, 100ini.similitud))%)")
    p2 = heatmap(fin.tablero; color = paleta, clims = (0, 2), aspect_ratio = :equal,
                 colorbar = false, axis = false, ticks = false, grid = false,
                 title = "Equilibrio con τ = 37.5% (similitud = $(round(Int, 100fin.similitud))%)")
    plot(p1, p2; layout = (1, 2), size = (820, 390))
end

# ╔═╡ 098ea854-ca7e-11f1-8096-ab63317295a8
md"""
Nadie buscaba una ciudad con $(round(Int, 100ejemplo_schelling.fin.similitud))%
de vecinos del mismo grupo: a cada agente le bastaba con 37.5%. El exceso de
segregación es una propiedad **emergente**: las intenciones individuales no se
suman linealmente cuando las decisiones de unos alteran el entorno de los demás.
"""

# ╔═╡ 098ea900-ca7e-11f1-98a1-a12ddf8e6f94
"Barrido de la similitud final y fracción de satisfechos en Schelling para varios umbrales τ."
function barrido_schelling(L, p_vacio; τs = 0.10:0.05:0.85)
    sim, fel = Float64[], Float64[]
    for (i, τ) in enumerate(τs)
        res = simular_schelling(L, p_vacio, τ; pasos_max = 100, rng = Xoshiro(30 + i))
        push!(sim, res.similitud[end])
        push!(fel, res.felices[end])
    end
    (τ = collect(τs), similitud = sim, felices = fel)
end

# ╔═╡ 098ea87e-ca7e-11f1-81e7-a34a017a4c91
barrido_sch = barrido_schelling(35, 0.15; τs = 0.10:0.05:0.85);

# ╔═╡ 098ea886-ca7e-11f1-ad45-f5ed08d5ae8b
let
    b = barrido_sch
    p = plot(100 .* b.τ, 100 .* b.similitud; color = AZUL, lw = 2.5, marker = :circle, ms = 4, msw = 0,
             label = "Similitud promedio final", xlabel = "umbral requerido τ (%)",
             ylabel = "porcentaje (%)", ylims = (0, 105), legend = :bottomleft,
             title = "Segregación emergente frente al umbral individual")
    plot!(p, 100 .* b.τ, 100 .* b.felices; color = NARANJA, lw = 2, ls = :dash, marker = :square, ms = 3.5, msw = 0,
          label = "Agentes satisfechos al final")
    plot!(p, 100 .* b.τ, 100 .* b.τ; color = GRIS, lw = 1, ls = :dot, label = "Similitud = umbral τ")
    hline!(p, [50]; color = GRIS_CLARO, lw = 1, label = "Mezcla al azar (50%)")
end

# ╔═╡ 098ea8f8-ca7e-11f1-84a7-0796ff3e285a
"Rango de la similitud final del barrido de Schelling para τ entre `lo` y `hi`, como texto."
function rango_sch(lo, hi)
    k = (barrido_sch.τ .>= lo - 1e-9) .& (barrido_sch.τ .<= hi + 1e-9)
    a, b = extrema(barrido_sch.similitud[k])
    "$(round(Int, 100a))% y $(round(Int, 100b))%"
end

# ╔═╡ 098ea886-ca7e-11f1-afbd-0b54faabf8d3
md"""
La curva muestra tres regímenes muy claros:

1. **``\tau \le 25\%`` (mezcla estable):** Casi todos los agentes ya cumplen su
   umbral en la configuración aleatoria inicial; casi nadie se muda y la
   similitud queda entre $(rango_sch(0.0, 0.25)), cerca de la mezcla al azar.
2. **``30\% \le \tau \le 70\%`` (segregación fuerte con preferencias leves):**
   En cuanto ``\tau`` pasa de 25% a 30%, las mudanzas de los pocos
   insatisfechos iniciales desatan un efecto dominó que sólo se detiene cuando
   la similitud llega a valores entre $(rango_sch(0.3, 0.7)), muy por encima de
   la línea punteada ``\tau`` que cada agente pedía.
3. **``\tau \ge 75\%`` (frustración / congelamiento):** Los agentes exigen
   vecindarios tan puros que no hay suficientes lugares satisfactorios; las
   mudanzas no logran acomodar a todos y la fracción de satisfechos cae.
"""

# ╔═╡ 098ea900-ca7e-11f1-965b-0775624a4b54
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

# ╔═╡ 098ea8b0-ca7e-11f1-b39d-5930fe3027f6
verificar(() -> [similitud_local(3, 8), similitud_local(0, 5), similitud_local(0, 0)],
          [0.375, 0.0, 1.0];
          pista = md"Usa un condicional: `totales == 0 ? 1.0 : iguales / totales`.")

# ╔═╡ 098ea8ba-ca7e-11f1-912b-0db665a8dad0
verificar(() -> [satisfecho(3, 8, 0.375), satisfecho(2, 8, 0.375), satisfecho(0, 0, 0.5)],
          [true, false, true];
          pista = md"Puedes usar `similitud_local(iguales, totales) >= τ` o `totales == 0 || iguales / totales >= τ`.")

# ╔═╡ Cell order:
# ╟─098ea810-ca7e-11f1-8ee4-11801deff225
# ╟─098ea842-ca7e-11f1-98cc-45514a1bfd32
# ╟─098ea84c-ca7e-11f1-a5ae-6f90cd7e2ac6
# ╟─098ea854-ca7e-11f1-9322-97fd00adc5e0
# ╟─098ea854-ca7e-11f1-8096-ab63317295a8
# ╟─098ea860-ca7e-11f1-9b64-d9cf51210681
# ╟─098ea86a-ca7e-11f1-a122-cb6b7cf72c41
# ╠═098ea874-ca7e-11f1-9b0b-05215f7d8e17
# ╟─098ea874-ca7e-11f1-a9c1-1b5ec4c6e520
# ╟─098ea874-ca7e-11f1-bc3b-43f7ce624683
# ╠═098ea87e-ca7e-11f1-81e7-a34a017a4c91
# ╟─098ea886-ca7e-11f1-ad45-f5ed08d5ae8b
# ╟─098ea886-ca7e-11f1-afbd-0b54faabf8d3
# ╟─098ea892-ca7e-11f1-8112-fbd603c1be12
# ╟─098ea892-ca7e-11f1-be80-f7a9bc70dcc6
# ╟─098ea89c-ca7e-11f1-8569-95ef9e0fa2af
# ╟─098ea8a6-ca7e-11f1-84d7-73f78f0de57b
# ╟─098ea8a6-ca7e-11f1-8442-27a78f706b0f
# ╠═098ea8b0-ca7e-11f1-beb8-b9027551565e
# ╟─098ea8b0-ca7e-11f1-b39d-5930fe3027f6
# ╟─098ea8ba-ca7e-11f1-8911-f90128439f0e
# ╠═098ea8ba-ca7e-11f1-b041-09b6566e3592
# ╟─098ea8ba-ca7e-11f1-912b-0db665a8dad0
# ╟─098ea8c4-ca7e-11f1-94bc-33cf6a71df93
# ╟─098ea8ce-ca7e-11f1-a706-a716ce70bf9a
# ╠═098ea8ce-ca7e-11f1-9fdd-e97268d7a5ad
# ╠═098ea8d8-ca7e-11f1-adb2-b75c2a497a01
# ╠═098ea8d8-ca7e-11f1-9de0-35f2992510ec
# ╠═098ea8e2-ca7e-11f1-b06a-49b2737fbc46
# ╠═098ea8e2-ca7e-11f1-82a4-837d0b1404aa
# ╠═098ea8ec-ca7e-11f1-8d00-d5ae5a76c5eb
# ╠═098ea8f8-ca7e-11f1-8bac-23e0b8c490e3
# ╠═098ea8f8-ca7e-11f1-84a7-0796ff3e285a
# ╠═098ea900-ca7e-11f1-98a1-a12ddf8e6f94
# ╠═098ea900-ca7e-11f1-965b-0775624a4b54
