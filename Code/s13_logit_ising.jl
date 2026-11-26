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

# ╔═╡ e17a225a-ca9b-11f1-8c9c-415708d55293
begin
    # Usa el entorno del curso (Code/Project.toml) en lugar del gestor de
    # paquetes propio de Pluto: un solo entorno para todo el repositorio.
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, Statistics, Random, Distributions, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

# ╔═╡ e17a2192-ca9b-11f1-8147-ad971ab2c48a
md"""
# Sesión 13: Logit, Boltzmann e Ising

**Objetivos de la sesión**
1. Repasar la **función logística** y el modelo **logit** de elección
   discreta.
2. Ver que el logit es la **distribución de Boltzmann**: el modelo que usan
   los economistas para estimar demanda es física estadística.
3. Agregar **interacciones sociales** a la elección (Brock y Durlauf) y
   obtener el **modelo de Ising**.
4. Simular el Ising con el algoritmo de **Metropolis** y encontrar su
   transición de fase.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s13_logit_ising.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): C.B4 (función logística),
F.B3 (temperatura), F.B4 (transiciones de fase), F.I1 (Boltzmann), F.I2
(logit = Boltzmann), F.I3 (Ising), F.I5 (Metropolis). Ruta avanzada: E.A3
(ecuación maestra).

**Conceptos transversales.** Hoy se presenta el **logit** (elección discreta),
que vuelve en las sesiones 18 y 20. Se retoman las **cadenas de Markov**
(sesión 4).
"""

# ╔═╡ e17a21ba-ca9b-11f1-8f79-fbaff49e9fec
TableOfContents(title = "Contenido")

# ╔═╡ e17a21c4-ca9b-11f1-b864-45e68b346d49
md"""
## 0. El fenómeno: mismos incentivos, resultados distintos

Glaeser, Sacerdote y Scheinkman (1996) notaron que las tasas de
criminalidad varían entre ciudades de Estados Unidos mucho más de lo que
explican sus diferencias en ingreso, educación o policía. Algo parecido pasa
con la deserción escolar, la fertilidad o la adopción de tecnologías: grupos
con incentivos casi iguales terminan comportándose de forma muy distinta.
Su explicación: las **interacciones sociales**. Cuando lo que hacen mis
vecinos cambia lo que me conviene hacer a mí, pequeñas diferencias
iniciales se amplifican.

Abajo hay dos "ciudades" de 60 × 60 habitantes con **exactamente las mismas
reglas**: cada uno elige entre dos opciones (azul o naranja), prefiere
coincidir con sus cuatro vecinos y tiene algo de ruido en sus decisiones.
Ninguna opción es mejor que la otra. Sólo cambia el azar.
"""

# ╔═╡ e17a21ce-ca9b-11f1-ad6e-cfcb4c3e69c3
md"""
Las dos ciudades empezaron con elecciones al azar, mitad y mitad. Al final,
en cada una domina una opción distinta, y el **consenso** no lo impuso nadie:
emergió de la imitación local. Hoy vamos a construir este modelo desde la
teoría de la elección que usan los economistas.
"""

# ╔═╡ e17a21e2-ca9b-11f1-85ed-49bc9346826e
md"""
## 1. Repaso: la función logística y el modelo logit

La **función logística** convierte cualquier número en una probabilidad
(repaso C.B4):

```math
\sigma(x) = \frac{1}{1 + e^{-x}} .
```

En el modelo **logit**, un agente elige entre A y B, cuyas utilidades
difieren en ``\Delta u = u_A - u_B``. La probabilidad de elegir A es

```math
\mathbb{P}(A) = \sigma(\beta\,\Delta u) = \frac{e^{\beta u_A}}{e^{\beta u_A} + e^{\beta u_B}} .
```

El parámetro ``\beta`` es la **intensidad de elección**. Con ``\beta`` grande,
el agente casi siempre elige la mejor opción; con ``\beta = 0``, elige al azar.

Intensidad de elección ``\beta``: $(@bind β_logit Slider([0.25, 0.5, 1.0, 2.0, 5.0, 20.0], default = 1.0, show_value = true))
"""

# ╔═╡ e17a21ec-ca9b-11f1-9d04-19bb8b8c356e
md"""
## 2. Logit = Boltzmann

¿De dónde sale la fórmula logit? McFadden la derivó suponiendo que la
utilidad tiene una parte que vemos, ``u_i``, y un **choque aleatorio**
``\varepsilon_i`` que no vemos (gustos, humor, información privada):

```math
U_i = u_i + \varepsilon_i .
```

El agente elige la opción con mayor ``U_i``. Si los ``\varepsilon_i`` son
independientes y tienen distribución de **Gumbel** con escala ``1/\beta``,
la probabilidad de que gane la opción ``i`` es exactamente

```math
\mathbb{P}(i) = \frac{e^{\beta u_i}}{\sum_j e^{\beta u_j}} .
```

Ahora compárala con la **distribución de Boltzmann** de la sesión 8
(repaso F.I1): un sistema a temperatura ``T = 1/\beta`` ocupa el estado ``i``,
con energía ``E_i``, con probabilidad ``e^{-\beta E_i}/Z``. Son **la misma
fórmula** con ``E_i = -u_i``: la utilidad es menos la energía, y la
"temperatura" ``1/\beta`` es el tamaño del ruido en las decisiones (repaso
F.I2 y F.B3).

Comprobémoslo simulando: tres opciones con utilidades 1, 0.5 y 0, y 100,000
agentes con choques de Gumbel.
"""

# ╔═╡ e17a21f6-ca9b-11f1-bd25-0186b4dab455
utilidades = [1.0, 0.5, 0.0];

# ╔═╡ e17a21f6-ca9b-11f1-96ca-5586e11f1ccd
elecciones = let
    rng = Xoshiro(1)
    β = 2.0
    [argmax(utilidades .+ rand(rng, Gumbel(0, 1 / β), 3)) for _ in 1:100_000]
end;

# ╔═╡ e17a2200-ca9b-11f1-9783-61aa9da1eb51
md"""
Las barras coinciden. Lo que en física se llama **agitación térmica**, en
economía se llama **heterogeneidad no observada**.
"""

# ╔═╡ e17a220a-ca9b-11f1-9ef5-b71a4d1dfc71
md"""
## 3. Elección con interacciones sociales: el modelo de Ising

Brock y Durlauf (2001) agregaron al logit un ingrediente: a cada agente le
importa **lo que hacen los demás**. Codifiquemos la elección del agente ``i``
como ``s_i = +1`` (adopta) o ``s_i = -1`` (no adopta). Su utilidad es

```math
U_i(s_i) = \underbrace{h\, s_i}_{\text{incentivo privado}}
+ \underbrace{J\, s_i \sum_{j \in \text{vecinos}} s_j}_{\text{conformidad}}
+ \varepsilon_i(s_i) .
```

- ``h`` es el incentivo común a adoptar (un subsidio, un precio).
- ``J > 0`` es el premio por **coincidir** con los vecinos.

Lo que decide cada agente es la diferencia de utilidad entre adoptar y no
adoptar, ``2\big(h + J\sum_{j} s_j\big)``. Es exactamente menos el cambio de
energía al cambiar ``s_i`` en el **modelo de Ising** (repaso F.I3):

```math
H(s) = -J\sum_{\langle i,j\rangle} s_i s_j - h\sum_i s_i ,
```

donde ``\langle i,j\rangle`` recorre cada par de vecinos una sola vez.
(Sumar las utilidades de todos los agentes contaría dos veces cada par: la
equivalencia es entre las decisiones, no entre las sumas.)

Los imanes se ordenan porque cada espín "prefiere" alinearse con sus
vecinos; las poblaciones llegan a consensos porque cada agente prefiere
coincidir con los suyos. El parámetro de orden es la **magnetización**,
que en lectura social es la **adopción neta**:

```math
m = \frac{1}{N}\sum_i s_i \in [-1, 1] .
```

Con ``m = 0`` la población está dividida a la mitad; con ``m = 1`` todos
adoptan.
"""

# ╔═╡ e17a220a-ca9b-11f1-bc35-ab9d1afc35e8
md"""
## 4. Laboratorio: el algoritmo de Metropolis

Para simular el modelo, cada agente revisa su elección de vez en cuando. El
algoritmo de **Metropolis** (repaso F.I5) hace esto:

1. Elige un agente ``i`` al azar.
2. Calcula cuánto cambiaría la energía si cambiara de opción:
   ``\Delta E = 2 s_i \left(J\sum_j s_j + h\right)``.
3. Si ``\Delta E \le 0`` (le conviene), cambia. Si no, cambia con
   probabilidad ``e^{-\beta\,\Delta E}``: a veces se equivoca, por el ruido.

``N`` revisiones forman un **barrido**. Tras muchos barridos, las
configuraciones aparecen con la probabilidad de Boltzmann, ``e^{-\beta H}/Z``
(por balance detallado; repaso E.A3). La regla logit de la sección 1, aplicada
a cada revisión, llega a la misma distribución; se llama **dinámica de baño
térmico**.

Ponemos a los agentes en una cuadrícula de 60 × 60, cada uno con 4 vecinos,
con ``J = 1`` y ``h = 0``. El único parámetro es la temperatura ``T = 1/\beta``:
el ruido.

Temperatura ``T``: $(@bind T_red Slider([1.0, 1.5, 2.0, 2.2, 2.3, 2.5, 3.0, 4.0], default = 2.0, show_value = true))
"""

# ╔═╡ e17a221e-ca9b-11f1-ac69-0939b07cb693
md"""
- Con ``T`` **baja** (poco ruido), la población llega a un consenso: casi todos
  eligen lo mismo. Cuál de las dos opciones gana lo decide el azar.
- Con ``T`` **alta** (mucho ruido), las elecciones se mezclan y ``m`` fluctúa
  alrededor de 0.
- Cerca de ``T \approx 2.27`` aparecen grupos de todos los tamaños, y ``m``
  da grandes saltos. Es la **temperatura crítica**.

Para ver la transición completa, medimos ``|m|`` promediado en la segunda
mitad de la simulación para muchas temperaturas y dos tamaños de red, y lo
comparamos con la solución exacta para una red infinita (Onsager, 1944;
Yang, 1952).
"""

# ╔═╡ e17a2228-ca9b-11f1-abf3-cbaa4d3f82b0
md"""
Debajo de ``T_c`` la simulación sigue a la curva exacta. Por encima, la red
finita conserva una magnetización residual que se reduce al crecer la red:
con 48 × 48 está más cerca de cero que con 16 × 16. Como en la sesión 9, la
transición sólo es abrupta en el límite de muchos agentes (repaso F.B4).
"""

# ╔═╡ e17a2232-ca9b-11f1-92c7-adfe9ce17a8d
md"""
## Resumen

| Concepto | Idea |
|---|---|
| Logit | ``\mathbb P(A) = \sigma(\beta\,\Delta u)``; ``\beta`` es la intensidad de elección |
| Logit = Boltzmann | Choques de Gumbel ↔ agitación térmica; ``E = -u``, ``T = 1/\beta`` |
| Brock-Durlauf | Logit + conformidad = modelo de Ising |
| Metropolis | Revisiones con ruido que muestrean ``e^{-\beta H}/Z`` |
| Transición de fase | Red cuadrada: consenso espontáneo debajo de ``T_c \approx 2.269\,J`` |
"""

# ╔═╡ e17a2232-ca9b-11f1-b136-a30dd93bebd6
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `prob_logit(Δu, β)`, la probabilidad de elegir A cuando su
ventaja en utilidad es `Δu`.
"""

# ╔═╡ e17a223c-ca9b-11f1-a0df-d77d798db4c9
function prob_logit(Δu, β)
    missing  # reemplaza esta línea
end

# ╔═╡ e17a2248-ca9b-11f1-9706-71d981da5b8b
md"""
**2.** *Si alcanza (margen de la clase).* Escribe `cambio_energia(s, suma_vecinos, J, h)`, el cambio de energía
``\Delta E`` si el agente con elección `s` cambia de opción, cuando la suma de
las elecciones de sus vecinos es `suma_vecinos`.
"""

# ╔═╡ e17a2248-ca9b-11f1-8cab-178bb1995d88
function cambio_energia(s, suma_vecinos, J, h)
    missing  # reemplaza esta línea
end

# ╔═╡ e17a2250-ca9b-11f1-a548-615f7fee17e8
md"""
### Ruta básica: interpretación

3. En el modelo logit, ¿qué representa ``\beta``? ¿Qué pasa con las elecciones
   cuando ``\beta \to 0`` y cuando ``\beta \to \infty``?
4. Explica con tus palabras por qué la utilidad de Brock y Durlauf lleva al
   modelo de Ising. ¿Qué papel juega cada término?
5. En la sección 0, dos ciudades con las mismas reglas terminan en consensos
   opuestos. ¿Qué implica esto para un estudio empírico que compara
   ciudades?


### Conceptos que vuelven

Preguntas cortas para repasar ideas de sesiones anteriores; se responden con
lo que ya sabes y con el laboratorio de hoy.

- **R1.** El algoritmo de Metropolis es una cadena de Markov sobre las
  configuraciones de la red. ¿Qué significa que llegue a una distribución
  estacionaria? ¿Por qué, en el barrido de temperatura, promediamos sólo los
  últimos 1500 de los 3000 barridos? (Sesión 4.)

### Ruta avanzada

6. Demuestra que si ``\varepsilon_A`` y ``\varepsilon_B`` son Gumbel
   independientes con escala ``1/\beta``, la diferencia
   ``\varepsilon_A - \varepsilon_B`` tiene distribución logística, y deduce la
   fórmula logit.
7. Repite el barrido de temperatura de la sección 4 con una red de
   Erdős-Rényi con grado promedio 4 en lugar de la cuadrícula. ¿Dónde está
   ahora la transición?

## Para la próxima sesión

- Lectura: Sen y Chakrabarti (2014), capítulo 2. Para la ruta avanzada,
  Brock y Durlauf (2001).
- En la sesión 14 resolveremos el modelo cuando cada agente mira el promedio
  de toda la población (**campo medio**) y veremos cómo cambios graduales en
  los incentivos producen saltos irreversibles (**histéresis**) y avalanchas.
"""

# ╔═╡ e17a225a-ca9b-11f1-9554-dd6ed105775e
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ e17a2264-ca9b-11f1-8ef7-ab926fe34aea
begin
    const AZUL       = "#2a78d6"   # serie principal, opción +1
    const NARANJA    = "#eb6834"   # segunda serie o teoría, opción -1
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const BLANCO     = "#ffffff"
    const RAMPA      = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # azul secuencial
end;

# ╔═╡ e17a2200-ca9b-11f1-8bb6-291db255f714
let
    β = 2.0
    teoria = exp.(β .* utilidades) ./ sum(exp.(β .* utilidades))
    sim = [mean(elecciones .== i) for i in 1:3]
    x = 1:3
    p = bar(x .- 0.18, sim; bar_width = 0.34, color = AZUL, lw = 0, label = "Simulación: elige la mayor U = u + ε",
            xticks = (x, ["Opción 1 (u = 1)", "Opción 2 (u = 0.5)", "Opción 3 (u = 0)"]),
            ylabel = "fracción que la elige", ylims = (0, 0.75), legend = :topright,
            title = "Choques de Gumbel con β = 2 producen el logit")
    bar!(p, x .+ 0.18, teoria; bar_width = 0.34, color = NARANJA, lw = 0, label = "Boltzmann: exp(βu) / Z")
end

# ╔═╡ e17a2264-ca9b-11f1-a790-5d307dd4a945
"Función logística."
σ(x) = 1 / (1 + exp(-x))

# ╔═╡ e17a21ec-ca9b-11f1-94ec-6bb252d6272f
let
    x = range(-5, 5, length = 300)
    p = plot(x, σ.(β_logit .* x); color = AZUL, lw = 2.5, ylims = (0, 1),
             xlabel = "ventaja de A en utilidad, Δu = u_A − u_B", ylabel = "probabilidad de elegir A",
             label = "β = $β_logit", legend = :bottomright, title = "Modelo logit")
    for b in (0.25, 1.0, 5.0)
        b == β_logit && continue
        plot!(p, x, σ.(b .* x); color = GRIS_CLARO, lw = 1.5, label = "β = $b")
    end
    hline!(p, [0.5]; color = GRIS, ls = :dot, label = "")
end

# ╔═╡ e17a226e-ca9b-11f1-b6d8-d5083c405b91
"""
Modelo de Ising en una cuadrícula L × L con fronteras periódicas, J = 1 y
h = 0, simulado con Metropolis. `inicio` puede ser `:azar` u `:ordenado`
(todos en +1). Devuelve la configuración final y la magnetización después
de cada barrido.
"""
function ising_red(L, T; barridos = 1_000, rng = Xoshiro(1), inicio = :azar)
    s = inicio == :azar ? rand(rng, [-1, 1], L, L) : ones(Int, L, L)
    β = 1 / T
    acepta = [exp(-β * ΔE) for ΔE in 1:8]    # probabilidades para ΔE > 0
    m = Float64[]
    for _ in 1:barridos
        for _ in 1:L^2
            i, j = rand(rng, 1:L), rand(rng, 1:L)
            vecinos = s[mod1(i + 1, L), j] + s[mod1(i - 1, L), j] + s[i, mod1(j + 1, L)] + s[i, mod1(j - 1, L)]
            ΔE = 2 * s[i, j] * vecinos
            if ΔE <= 0 || rand(rng) < acepta[ΔE]
                s[i, j] = -s[i, j]
            end
        end
        push!(m, mean(s))
    end
    (s = s, m = m)
end

# ╔═╡ e17a21ce-ca9b-11f1-b736-a1f5290ba6e0
let
    ps = map((3, 8)) do semilla
        s = ising_red(60, 2.0; barridos = 400, rng = Xoshiro(semilla), inicio = :azar).s
        heatmap(s; color = cgrad([NARANJA, AZUL]), clims = (-1, 1), aspect_ratio = :equal,
                colorbar = false, axis = false, ticks = false, grid = false,
                title = "Ciudad $(semilla == 3 ? "A" : "B"): $(round(Int, 100 * mean(s .== 1)))% elige azul")
    end
    plot(ps...; layout = (1, 2), size = (760, 400))
end

# ╔═╡ e17a2216-ca9b-11f1-8f73-b51b12506cfb
simulacion = ising_red(60, T_red; barridos = 1_000, rng = Xoshiro(2));

# ╔═╡ e17a2216-ca9b-11f1-9828-e5f32cabe084
let
    s = simulacion.s
    p1 = heatmap(s; color = cgrad([NARANJA, AZUL]), clims = (-1, 1), aspect_ratio = :equal,
                 colorbar = false, axis = false, ticks = false, grid = false,
                 title = "Configuración final, T = $T_red")
    t = 1:length(simulacion.m)
    p2 = plot(t, simulacion.m; color = AZUL, lw = 1.5, ylims = (-1.05, 1.05), legend = false,
              xlabel = "barrido", ylabel = "adopción neta m",
              title = "Magnetización a lo largo de la simulación")
    hline!(p2, [0]; color = GRIS, ls = :dot)
    plot(p1, p2; layout = grid(1, 2; widths = [0.42, 0.58]), size = (900, 380))
end

# ╔═╡ e17a221e-ca9b-11f1-a082-470a75e4b629
barrido_T = let
    Ts = 1.0:0.1:4.0
    [(L = L, T = Ts, m = [mean(abs.(ising_red(L, T; barridos = 3_000, rng = Xoshiro(round(Int, 100T)),
                                               inicio = :ordenado).m[1501:end])) for T in Ts])
     for L in (16, 48)]
end;

# ╔═╡ e17a227a-ca9b-11f1-ac11-4f38796e27b8
"Magnetización espontánea exacta de la red cuadrada (Onsager y Yang), con J = 1."
onsager(T) = T < 2 / log(1 + sqrt(2)) ? (1 - sinh(2 / T)^(-4))^(1 / 8) : 0.0

# ╔═╡ e17a2228-ca9b-11f1-9f22-f51eaf77a7d9
let
    Tc = 2 / log(1 + sqrt(2))
    x = range(1.0, 4.0, length = 400)
    p = plot(x, onsager.(x); color = NARANJA, lw = 2.5, label = "Solución exacta (red infinita)",
             xlabel = "temperatura T", ylabel = "|m| promedio", ylims = (0, 1.05), legend = :topright,
             title = "Transición de fase del modelo de Ising en la red cuadrada")
    for (d, c, mk) in zip(barrido_T, (RAMPA[1], RAMPA[3]), (:diamond, :circle))
        scatter!(p, d.T, d.m; color = c, ms = 4, msw = 0, marker = mk, label = "Metropolis, $(d.L) × $(d.L)")
    end
    vline!(p, [Tc]; color = GRIS, ls = :dash, label = "T_c ≈ 2.269")
end

# ╔═╡ e17a227a-ca9b-11f1-852a-a7915bd93a13
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

# ╔═╡ e17a223c-ca9b-11f1-84ff-39692be8420a
verificar(() -> [prob_logit(0.0, 3.0), prob_logit(1.0, 2.0), prob_logit(-1.0, 0.5)],
          [0.5, 1 / (1 + exp(-2.0)), 1 / (1 + exp(0.5))];
          pista = md"Es la logística evaluada en ``\beta\,\Delta u``: `1 / (1 + exp(-β * Δu))`.")

# ╔═╡ e17a2250-ca9b-11f1-b0d2-715169c25939
verificar(() -> [cambio_energia(1, 4, 1.0, 0.0), cambio_energia(-1, 2, 1.0, 0.5), cambio_energia(1, -2, 0.5, 0.0)],
          [8.0, -5.0, -2.0];
          pista = md"La fórmula es ``\Delta E = 2 s_i (J \sum_j s_j + h)``.")

# ╔═╡ Cell order:
# ╟─e17a2192-ca9b-11f1-8147-ad971ab2c48a
# ╟─e17a21ba-ca9b-11f1-8f79-fbaff49e9fec
# ╟─e17a21c4-ca9b-11f1-b864-45e68b346d49
# ╟─e17a21ce-ca9b-11f1-b736-a1f5290ba6e0
# ╟─e17a21ce-ca9b-11f1-ad6e-cfcb4c3e69c3
# ╟─e17a21e2-ca9b-11f1-85ed-49bc9346826e
# ╟─e17a21ec-ca9b-11f1-94ec-6bb252d6272f
# ╟─e17a21ec-ca9b-11f1-9d04-19bb8b8c356e
# ╠═e17a21f6-ca9b-11f1-bd25-0186b4dab455
# ╠═e17a21f6-ca9b-11f1-96ca-5586e11f1ccd
# ╟─e17a2200-ca9b-11f1-8bb6-291db255f714
# ╟─e17a2200-ca9b-11f1-9783-61aa9da1eb51
# ╟─e17a220a-ca9b-11f1-9ef5-b71a4d1dfc71
# ╟─e17a220a-ca9b-11f1-bc35-ab9d1afc35e8
# ╠═e17a2216-ca9b-11f1-8f73-b51b12506cfb
# ╟─e17a2216-ca9b-11f1-9828-e5f32cabe084
# ╟─e17a221e-ca9b-11f1-ac69-0939b07cb693
# ╠═e17a221e-ca9b-11f1-a082-470a75e4b629
# ╟─e17a2228-ca9b-11f1-9f22-f51eaf77a7d9
# ╟─e17a2228-ca9b-11f1-abf3-cbaa4d3f82b0
# ╟─e17a2232-ca9b-11f1-92c7-adfe9ce17a8d
# ╟─e17a2232-ca9b-11f1-b136-a30dd93bebd6
# ╠═e17a223c-ca9b-11f1-a0df-d77d798db4c9
# ╟─e17a223c-ca9b-11f1-84ff-39692be8420a
# ╟─e17a2248-ca9b-11f1-9706-71d981da5b8b
# ╠═e17a2248-ca9b-11f1-8cab-178bb1995d88
# ╟─e17a2250-ca9b-11f1-b0d2-715169c25939
# ╟─e17a2250-ca9b-11f1-a548-615f7fee17e8
# ╟─e17a225a-ca9b-11f1-9554-dd6ed105775e
# ╠═e17a225a-ca9b-11f1-8c9c-415708d55293
# ╠═e17a2264-ca9b-11f1-8ef7-ab926fe34aea
# ╠═e17a2264-ca9b-11f1-a790-5d307dd4a945
# ╠═e17a226e-ca9b-11f1-b6d8-d5083c405b91
# ╠═e17a227a-ca9b-11f1-ac11-4f38796e27b8
# ╠═e17a227a-ca9b-11f1-852a-a7915bd93a13
