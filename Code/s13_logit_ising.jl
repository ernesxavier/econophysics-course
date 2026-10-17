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

# ╔═╡ 8bc764f8-ca82-11f1-b48c-d1fac4639544
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

# ╔═╡ 8bc76412-ca82-11f1-bd86-95fcd9064a5e
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
"""

# ╔═╡ 8bc7644e-ca82-11f1-b7fb-4d6754ba2a46
TableOfContents(title = "Contenido")

# ╔═╡ 8bc7644e-ca82-11f1-b560-81661b82a26d
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

# ╔═╡ 8bc76462-ca82-11f1-adbc-6b25970ae0e7
md"""
Las dos ciudades empezaron con elecciones al azar, mitad y mitad. Al final,
en cada una domina una opción distinta, y el **consenso** no lo impuso nadie:
emergió de la imitación local. Hoy vamos a construir este modelo desde la
teoría de la elección que usan los economistas.
"""

# ╔═╡ 8bc7646c-ca82-11f1-bec8-d7ebd834aea4
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

# ╔═╡ 8bc76476-ca82-11f1-b805-73abcb89a177
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

# ╔═╡ 8bc76480-ca82-11f1-9f24-051825b1bac0
utilidades = [1.0, 0.5, 0.0];

# ╔═╡ 8bc7648a-ca82-11f1-9c78-51696d03abc8
elecciones = let
    rng = Xoshiro(1)
    β = 2.0
    [argmax(utilidades .+ rand(rng, Gumbel(0, 1 / β), 3)) for _ in 1:100_000]
end;

# ╔═╡ 8bc7648a-ca82-11f1-8f10-3bb13b231a92
md"""
Las barras coinciden. Lo que en física se llama **agitación térmica**, en
economía se llama **heterogeneidad no observada**.
"""

# ╔═╡ 8bc76494-ca82-11f1-ac8e-115ee83ad93f
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

La suma de las dos primeras partes, sobre todos los agentes, es menos la
energía del **modelo de Ising** (repaso F.I3):

```math
H(s) = -J\sum_{\langle i,j\rangle} s_i s_j - h\sum_i s_i .
```

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

# ╔═╡ 8bc7649e-ca82-11f1-a1fc-156ca27ffdac
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

# ╔═╡ 8bc764a8-ca82-11f1-bf3e-d5707c1b667e
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

# ╔═╡ 8bc764bc-ca82-11f1-915d-b3a39f983380
md"""
Debajo de ``T_c`` la simulación sigue a la curva exacta. Por encima, la red
finita conserva una magnetización residual que se reduce al crecer la red:
con 48 × 48 está más cerca de cero que con 16 × 16. Como en la sesión 9, la
transición sólo es abrupta en el límite de muchos agentes (repaso F.B4).
"""

# ╔═╡ 8bc764bc-ca82-11f1-9a0d-67b1c3fc5743
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

# ╔═╡ 8bc764c6-ca82-11f1-bb07-81215699405a
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `prob_logit(Δu, β)`, la probabilidad de elegir A cuando su
ventaja en utilidad es `Δu`.
"""

# ╔═╡ 8bc764d0-ca82-11f1-8c38-c15321083f90
function prob_logit(Δu, β)
    missing  # reemplaza esta línea
end

# ╔═╡ 8bc764dc-ca82-11f1-ac77-07dd60b45ee8
md"""
**2.** Escribe `cambio_energia(s, suma_vecinos, J, h)`, el cambio de energía
``\Delta E`` si el agente con elección `s` cambia de opción, cuando la suma de
las elecciones de sus vecinos es `suma_vecinos`.
"""

# ╔═╡ 8bc764dc-ca82-11f1-9c63-0b1c3088a1b3
function cambio_energia(s, suma_vecinos, J, h)
    missing  # reemplaza esta línea
end

# ╔═╡ 8bc764e4-ca82-11f1-93f9-0d827c5ed0f7
md"""
### Ruta básica: interpretación

3. En el modelo logit, ¿qué representa ``\beta``? ¿Qué pasa con las elecciones
   cuando ``\beta \to 0`` y cuando ``\beta \to \infty``?
4. Explica con tus palabras por qué la utilidad de Brock y Durlauf lleva al
   modelo de Ising. ¿Qué papel juega cada término?
5. En la sección 0, dos ciudades con las mismas reglas terminan en consensos
   opuestos. ¿Qué implica esto para un estudio empírico que compara
   ciudades?

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

# ╔═╡ 8bc764ee-ca82-11f1-8324-7f59b36397ad
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 8bc764f8-ca82-11f1-b557-6923b433821d
begin
    const AZUL       = "#2a78d6"   # serie principal, opción +1
    const NARANJA    = "#eb6834"   # segunda serie o teoría, opción -1
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const BLANCO     = "#ffffff"
    const RAMPA      = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # azul secuencial
end;

# ╔═╡ 8bc7648a-ca82-11f1-8495-cd69d98c55f4
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

# ╔═╡ 8bc76502-ca82-11f1-b829-0b606b9d44b6
"Función logística."
σ(x) = 1 / (1 + exp(-x))

# ╔═╡ 8bc76476-ca82-11f1-84ca-5125f23115c9
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

# ╔═╡ 8bc76502-ca82-11f1-8187-8f2e6ca90a49
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

# ╔═╡ 8bc76458-ca82-11f1-99d4-db296032e8eb
let
    ps = map((3, 8)) do semilla
        s = ising_red(60, 2.0; barridos = 400, rng = Xoshiro(semilla), inicio = :azar).s
        heatmap(s; color = cgrad([NARANJA, AZUL]), clims = (-1, 1), aspect_ratio = :equal,
                colorbar = false, axis = false, ticks = false, grid = false,
                title = "Ciudad $(semilla == 3 ? "A" : "B"): $(round(Int, 100 * mean(s .== 1)))% elige azul")
    end
    plot(ps...; layout = (1, 2), size = (760, 400))
end

# ╔═╡ 8bc7649e-ca82-11f1-9e8e-03dc487f945c
simulacion = ising_red(60, T_red; barridos = 1_000, rng = Xoshiro(2));

# ╔═╡ 8bc764a8-ca82-11f1-95f0-e7ee76cbd481
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

# ╔═╡ 8bc764b2-ca82-11f1-991b-5776a30a442a
barrido_T = let
    Ts = 1.0:0.1:4.0
    [(L = L, T = Ts, m = [mean(abs.(ising_red(L, T; barridos = 3_000, rng = Xoshiro(round(Int, 100T)),
                                               inicio = :ordenado).m[1501:end])) for T in Ts])
     for L in (16, 48)]
end;

# ╔═╡ 8bc7650e-ca82-11f1-a822-0b0d8b043c30
"Magnetización espontánea exacta de la red cuadrada (Onsager y Yang), con J = 1."
onsager(T) = T < 2 / log(1 + sqrt(2)) ? (1 - sinh(2 / T)^(-4))^(1 / 8) : 0.0

# ╔═╡ 8bc764b2-ca82-11f1-818f-d39b10eaa846
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

# ╔═╡ 8bc7650e-ca82-11f1-ba04-572fa70013c9
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

# ╔═╡ 8bc764d0-ca82-11f1-9c39-0f9536856af7
verificar(() -> [prob_logit(0.0, 3.0), prob_logit(1.0, 2.0), prob_logit(-1.0, 0.5)],
          [0.5, 1 / (1 + exp(-2.0)), 1 / (1 + exp(0.5))];
          pista = md"Es la logística evaluada en ``\beta\,\Delta u``: `1 / (1 + exp(-β * Δu))`.")

# ╔═╡ 8bc764e4-ca82-11f1-8821-316af33bf012
verificar(() -> [cambio_energia(1, 4, 1.0, 0.0), cambio_energia(-1, 2, 1.0, 0.5), cambio_energia(1, -2, 0.5, 0.0)],
          [8.0, -5.0, -2.0];
          pista = md"La fórmula es ``\Delta E = 2 s_i (J \sum_j s_j + h)``.")

# ╔═╡ Cell order:
# ╟─8bc76412-ca82-11f1-bd86-95fcd9064a5e
# ╟─8bc7644e-ca82-11f1-b7fb-4d6754ba2a46
# ╟─8bc7644e-ca82-11f1-b560-81661b82a26d
# ╟─8bc76458-ca82-11f1-99d4-db296032e8eb
# ╟─8bc76462-ca82-11f1-adbc-6b25970ae0e7
# ╟─8bc7646c-ca82-11f1-bec8-d7ebd834aea4
# ╟─8bc76476-ca82-11f1-84ca-5125f23115c9
# ╟─8bc76476-ca82-11f1-b805-73abcb89a177
# ╠═8bc76480-ca82-11f1-9f24-051825b1bac0
# ╠═8bc7648a-ca82-11f1-9c78-51696d03abc8
# ╟─8bc7648a-ca82-11f1-8495-cd69d98c55f4
# ╟─8bc7648a-ca82-11f1-8f10-3bb13b231a92
# ╟─8bc76494-ca82-11f1-ac8e-115ee83ad93f
# ╟─8bc7649e-ca82-11f1-a1fc-156ca27ffdac
# ╠═8bc7649e-ca82-11f1-9e8e-03dc487f945c
# ╟─8bc764a8-ca82-11f1-95f0-e7ee76cbd481
# ╟─8bc764a8-ca82-11f1-bf3e-d5707c1b667e
# ╠═8bc764b2-ca82-11f1-991b-5776a30a442a
# ╟─8bc764b2-ca82-11f1-818f-d39b10eaa846
# ╟─8bc764bc-ca82-11f1-915d-b3a39f983380
# ╟─8bc764bc-ca82-11f1-9a0d-67b1c3fc5743
# ╟─8bc764c6-ca82-11f1-bb07-81215699405a
# ╠═8bc764d0-ca82-11f1-8c38-c15321083f90
# ╟─8bc764d0-ca82-11f1-9c39-0f9536856af7
# ╟─8bc764dc-ca82-11f1-ac77-07dd60b45ee8
# ╠═8bc764dc-ca82-11f1-9c63-0b1c3088a1b3
# ╟─8bc764e4-ca82-11f1-8821-316af33bf012
# ╟─8bc764e4-ca82-11f1-93f9-0d827c5ed0f7
# ╟─8bc764ee-ca82-11f1-8324-7f59b36397ad
# ╠═8bc764f8-ca82-11f1-b48c-d1fac4639544
# ╠═8bc764f8-ca82-11f1-b557-6923b433821d
# ╠═8bc76502-ca82-11f1-b829-0b606b9d44b6
# ╠═8bc76502-ca82-11f1-8187-8f2e6ca90a49
# ╠═8bc7650e-ca82-11f1-a822-0b0d8b043c30
# ╠═8bc7650e-ca82-11f1-ba04-572fa70013c9
