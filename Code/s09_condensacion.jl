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

# ╔═╡ 064ec8e0-ca80-11f1-ba30-1b07f57af930
begin
    # Usa el entorno del curso (Code/Project.toml) en lugar del gestor de
    # paquetes propio de Pluto: un solo entorno para todo el repositorio.
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, CSV, DataFrames, Statistics, Distributions, Random, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

# ╔═╡ 064ec7e4-ca80-11f1-8589-cd7f3d5d737c
md"""
# Sesión 9: Condensación de riqueza y redistribución

**Objetivos de la sesión**
1. Combinar las dos ideas de las sesiones anteriores: **crecimiento
   multiplicativo** (sesión 6) e **intercambio** (sesión 8), en el modelo de
   Bouchaud y Mézard.
2. Obtener una cola de **Pareto** cuyo exponente depende de cuánto
   intercambian los agentes.
3. Encontrar una **transición de fase**: si el intercambio es débil, la
   riqueza se **condensa** en unos pocos agentes.
4. Usar el modelo como laboratorio para **impuestos y redistribución**.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales. Algunas celdas simulan miles de agentes y tardan unos segundos.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s09_condensacion.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): E.B3 (promediar
simulaciones), F.B4 (transiciones de fase), E.I6 (Gibrat y Kesten), F.I6
(modelos de intercambio). Ruta avanzada: E.A4 (Fokker-Planck).
"""

# ╔═╡ 064ec816-ca80-11f1-bc22-3df113ff94d2
TableOfContents(title = "Contenido")

# ╔═╡ 064ec816-ca80-11f1-84d5-df39c2467675
md"""
## 0. El problema: inversión sin intercambio

En la sesión 8 el dinero sólo cambiaba de manos; nadie lo hacía crecer. En
la realidad, la riqueza también **se invierte**: cada agente obtiene un
rendimiento aleatorio, proporcional a lo que tiene. Es la ley de Gibrat de la
sesión 6:

```math
dW_i = W_i\, (\text{rendimiento aleatorio}) .
```

¿Qué pasa con la **fracción** de la riqueza total que tiene cada agente si
sólo hay inversión? Simulemos ``N = 1000`` agentes que empiezan iguales. La
gráfica muestra la fracción que tiene el agente más rico en cada momento.
"""

# ╔═╡ 064ec82c-ca80-11f1-b491-f14cf2752876
md"""
Sin intercambio, la riqueza se **concentra** sin límite: en la sesión 6 vimos
que la varianza de ``\ln W`` crece como ``\sigma^2 t``, así que con el tiempo un
solo agente, el más afortunado, acaba con casi todo. Con intercambio, la
fracción del más rico se estabiliza. La pregunta de hoy es **cuánto**
intercambio hace falta.
"""

# ╔═╡ 064ec836-ca80-11f1-978d-d395939ffd06
md"""
## 1. El modelo de Bouchaud y Mézard

Bouchaud y Mézard (2000) combinaron las dos fuerzas:

```math
dW_i = \underbrace{\sigma\, W_i\, dB_i}_{\text{inversión}}
\;+\; \underbrace{J\,\big(\bar W - W_i\big)\, dt}_{\text{intercambio}} ,
```

donde ``dB_i`` son ruidos independientes (repaso E.I6) y ``\bar W`` es la
riqueza promedio. El término de intercambio dice que cada agente recibe en
promedio ``J \bar W\, dt`` de los demás (compras, salarios, servicios) y
gasta ``J W_i\, dt``: los ricos gastan más, en proporción a su riqueza.

Como sólo nos interesa cómo se **reparte** la riqueza, trabajamos con la
riqueza relativa ``w_i = W_i/\bar W``, cuyo promedio es 1. El único
parámetro que importa es el cociente ``J/\sigma^2``: aquí fijamos ``\sigma = 1``.

**Resultado** (sección 5): la distribución estacionaria de ``w`` es

```math
p(w) = \frac{(\mu-1)^\mu}{\Gamma(\mu)}\; \frac{e^{-(\mu-1)/w}}{w^{1+\mu}},
\qquad \mu = 1 + \frac{2J}{\sigma^2} ,
```

una distribución **gamma inversa** con cola de Pareto de exponente ``\mu``.
Más intercambio, cola más delgada.

Intensidad de intercambio ``J``: $(@bind J_dist Slider([0.1, 0.25, 0.5, 1.0, 2.0], default = 0.5, show_value = true))
"""

# ╔═╡ 064ec854-ca80-11f1-87e0-9771df17749d
md"""
## 2. ¿Qué pasa cuando ``\mu`` se acerca a 1?

Una distribución de Pareto con ``\mu \le 1`` no tiene media finita. En una
población de ``N`` agentes eso significa que la riqueza **total** está
dominada por los más ricos: unos pocos tienen una fracción finita de todo,
por grande que sea ``N``. A esto se le llama **condensación**, por analogía
con la condensación de Bose-Einstein, donde una fracción finita de las
partículas ocupa un solo estado.

Para medirla usamos la **razón de participación inversa**

```math
Y_2 = \sum_{i=1}^{N} \left(\frac{W_i}{\sum_j W_j}\right)^2 .
```

Si todos tienen lo mismo, ``Y_2 = 1/N``; si uno tiene todo, ``Y_2 = 1``. En
general ``1/Y_2`` es el "número efectivo" de agentes que tienen la riqueza.
La prueba de la condensación es el **tamaño**: si al aumentar ``N`` el valor
de ``Y_2`` tiende a cero, la riqueza está repartida; si se queda en un valor
positivo, está condensada. ``Y_2`` es el **parámetro de orden** de la
transición (repaso F.B4).

En el modelo de la sección 1, ``\mu = 1 + 2J/\sigma^2 > 1`` para todo ``J > 0``:
nunca hay condensación, aunque con ``J`` pequeño se acerca. Pero ese modelo
supone que cada agente intercambia con **todos** los demás (campo medio). En
una economía real cada quien comercia con unos pocos.
"""

# ╔═╡ 064ec85e-ca80-11f1-a519-dd4ddc68e300
md"""
## 3. Laboratorio: la transición de condensación

Pongamos a los agentes en una **red aleatoria** donde cada uno tiene en
promedio 4 socios comerciales, y sólo intercambia con ellos:

```math
dW_i = \sigma\, W_i\, dB_i + \frac{J}{c_i}\sum_{j \in \text{socios de } i} (W_j - W_i)\, dt ,
```

con ``c_i`` el número de socios de ``i``. La gráfica muestra ``Y_2`` en el
estado estacionario para varios valores de ``J`` y tres tamaños de población,
de 1,000 a 20,000 agentes (mediana de 4 simulaciones). Compara cada modelo
consigo mismo al crecer ``N``.

Estas simulaciones tardan varios minutos, así que se corrieron una vez con
`Code/simular_condensacion.jl` y aquí sólo se leen sus resultados, en
`Data/condensacion_bm.csv`.
"""

# ╔═╡ 064ec85e-ca80-11f1-a178-8d908092580f
barrido = let
    d = CSV.read(joinpath(@__DIR__, "..", "Data", "condensacion_bm.csv"), DataFrame)
    combine(groupby(d, [:modelo, :J, :N]), :Y2 => median => :Y2)
end

# ╔═╡ 064ec872-ca80-11f1-8614-014a54cb5fcc
md"""
- En **campo medio** (izquierda), ``Y_2`` baja al aumentar ``N`` para todo
  ``J``: no hay condensación, aunque con ``J`` pequeño la bajada es lenta.
- En la **red** (derecha), con ``J \le 0.1`` las tres curvas coinciden:
  ``Y_2`` **no baja** al pasar de 1,000 a 20,000 agentes, y unos cuantos
  tienen una fracción finita de toda la riqueza. Con ``J \ge 0.4``, ``Y_2``
  vuelve a bajar con ``N``. Entre ambos valores hay una **transición de
  fase** en un ``J_c`` finito.

Con estos tamaños no podemos ubicar ``J_c`` con precisión: entre 0.15 y 0.3
las curvas se cruzan, y cerca de una transición las simulaciones finitas
son muy ruidosas. Bouchaud y Mézard, midiendo el exponente de la cola,
estimaron ``J/\sigma_{BM}^2 \approx 0.3`` para esta red; con su convención,
``\sigma_{BM}^2 = \sigma^2/2``, eso es ``J_c \approx 0.15`` en nuestras
unidades, dentro de la zona de transición que vemos.

**¿Por qué la red favorece la condensación?** Un agente afortunado sólo
"derrama" su riqueza a sus pocos socios, que a su vez la pasan lentamente
al resto. La riqueza fluye peor que en campo medio, y la inversión gana la
carrera.

**Visualmente.** Ordena a los agentes de más rico a más pobre y acumula su
riqueza. Con ``J`` pequeño en la red aparece una "escalera": unos pocos
escalones grandes concentran casi todo.

``J`` en la red: $(@bind J_esc Slider([0.05, 0.1, 0.2, 0.4, 1.0], default = 0.05, show_value = true))
"""

# ╔═╡ 064ec87c-ca80-11f1-8198-3f8aa828e2ef
md"""
## 4. Laboratorio: impuestos y redistribución

Un **impuesto a la riqueza** con tasa ``\varphi`` que se reparte por igual a
todos agrega un término

```math
\varphi\,\big(\bar W - W_i\big)\, dt ,
```

que tiene **exactamente la misma forma** que el intercambio en campo medio.
Por eso, en campo medio, el impuesto equivale a aumentar ``J``:

```math
\mu = 1 + \frac{2(J + \varphi)}{\sigma^2} .
```

Fijemos un intercambio débil, ``J = 0.1``, y variemos el impuesto.

Tasa del impuesto ``\varphi``: $(@bind φ_imp Slider([0.0, 0.025, 0.05, 0.1, 0.2, 0.4], default = 0.1, show_value = true))
"""

# ╔═╡ 064ec890-ca80-11f1-82ba-d9ba11b78736
barrido_impuesto = let
    d = CSV.read(joinpath(@__DIR__, "..", "Data", "condensacion_impuesto.csv"), DataFrame)
    combine(groupby(d, [:φ, :N]), :Y2 => median => :Y2)
end

# ╔═╡ 064ec89a-ca80-11f1-a0bf-e5aa1e6e4444
md"""
- **Sin impuesto**, ``Y_2`` no cambia al pasar de 1,000 a 20,000 agentes: la
  economía está condensada.
- Con ``\varphi = 0.05``, ``Y_2`` ya es la cuarta parte que sin impuesto con
  20,000 agentes, aunque baja despacio al crecer ``N``.
- Con ``\varphi = 0.1`` o más, ``Y_2`` baja claramente al crecer ``N``: la
  economía sale de la fase condensada.

El impuesto **redistribuye a todos**, no sólo a los socios, y reconecta a la
economía como si fuera de campo medio. Un impuesto del mismo tamaño que el
intercambio, ``J = 0.1``, es mucho más eficaz porque llega a todos.

**Pero el diseño importa.** Bouchaud y Mézard muestran que un impuesto a la
riqueza que **no** se redistribuye, combinado con un impuesto al ingreso,
puede **aumentar** la desigualdad: lo que reduce la desigualdad no es
cobrar, sino repartir.
"""

# ╔═╡ 064ec8a4-ca80-11f1-a789-a349ee856539
md"""
## 5. [Avanzado] La solución estacionaria de Fokker-Planck

Para ``N`` grande, el promedio ``\bar W`` crece de forma determinista y la
riqueza relativa ``w = W/\bar W`` sigue (en la convención de Itô, con
``\sigma`` y ``J`` como arriba)

```math
dw = J\,(1 - w)\, dt + \sigma\, w\, dB .
```

Es una difusión con deriva ``A(w) = J(1-w)``, que empuja hacia el promedio,
y ruido ``B(w) = \sigma w``, proporcional a la riqueza. Su densidad cumple la
ecuación de Fokker-Planck (repaso E.A4)

```math
\frac{\partial p}{\partial t} = -\frac{\partial}{\partial w}\big[J(1-w)\,p\big]
+ \frac{\sigma^2}{2}\frac{\partial^2}{\partial w^2}\big[w^2 p\big] .
```

La solución estacionaria con flujo nulo es

```math
p(w) \propto \frac{1}{\sigma^2 w^2}\exp\!\left(\int^w \frac{2J(1-y)}{\sigma^2 y^2}\, dy\right)
= \frac{C}{w^2}\exp\!\left(-\frac{2J}{\sigma^2}\Big(\frac{1}{w} + \ln w\Big)\right)
= C\, \frac{e^{-(\mu-1)/w}}{w^{1+\mu}},
```

con ``\mu = 1 + 2J/\sigma^2``. Para ``w`` grande el factor exponencial tiende a
1 y queda la ley de potencias ``w^{-1-\mu}``. Para ``w`` cercano a cero el factor
``e^{-(\mu-1)/w}`` aplasta la densidad: nadie queda en la ruina, porque el
intercambio siempre trae algo.

Es el mismo mecanismo que el proceso de Kesten de la sesión 6: un término
multiplicativo (la inversión) y uno aditivo (lo que se recibe de los demás).
"""

# ╔═╡ 064ec8ae-ca80-11f1-a0e9-59bfe2882fce
md"""
## Resumen

| Situación | Distribución de la riqueza relativa | ¿Condensación? |
|---|---|---|
| Sólo inversión (``J = 0``) | Lognormal que se ensancha sin límite | Sí: un agente acaba con casi todo |
| Inversión + intercambio, campo medio | Gamma inversa, Pareto con ``\mu = 1 + 2J/\sigma^2`` | No para ``J > 0`` |
| Inversión + intercambio en red (4 socios) | Pareto en la cola | Sí para ``J < J_c``, con ``J_c`` entre 0.1 y 0.4 (``\sigma = 1``) |
| Con impuesto redistribuido ``\varphi`` | En campo medio, ``\mu = 1 + 2(J + \varphi)/\sigma^2`` | Un impuesto pequeño la elimina |
"""

# ╔═╡ 064ec8ae-ca80-11f1-b425-0b8ab50dd798
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `participacion(w)`, que calcule la razón de participación
inversa ``Y_2`` de un vector de riquezas (no necesariamente normalizado).
"""

# ╔═╡ 064ec8ba-ca80-11f1-bd5d-3ba0b1ff6df2
function participacion(w)
    missing  # reemplaza esta línea
end

# ╔═╡ 064ec8c2-ca80-11f1-b0ff-4bfb46c063f1
md"""
**2.** Escribe `exponente_pareto(J, φ, σ)`, que devuelva el exponente ``\mu``
de la cola en campo medio con intercambio ``J``, impuesto redistribuido
``\varphi`` y volatilidad ``\sigma``.
"""

# ╔═╡ 064ec8cc-ca80-11f1-8a71-e95e82753218
function exponente_pareto(J, φ, σ)
    missing  # reemplaza esta línea
end

# ╔═╡ 064ec8d6-ca80-11f1-9a89-9d864cda096e
md"""
### Ruta básica: interpretación

3. ¿Por qué, sin intercambio, la riqueza acaba concentrada en un solo
   agente aunque todos tengan las mismas oportunidades?
4. En el modelo, ¿qué representa ``J`` en una economía real? Da dos ejemplos
   de cosas que lo aumentarían.
5. Explica con tus palabras por qué la condensación aparece en la red y no
   en campo medio.
6. Si ``Y_2 = 0.05``, ¿cuántos agentes "efectivos" tienen la riqueza?
7. Según el modelo, ¿qué es más eficaz contra la condensación: más
   intercambio entre socios o un impuesto pequeño redistribuido a todos?
   ¿Qué supuestos del modelo podrían cambiar esta conclusión?

### Ruta avanzada

8. Verifica que la densidad gamma inversa de la sección 1 está normalizada
   y tiene media 1 (repaso C.A4).
9. Demuestra que si ``W`` tiene cola de Pareto con ``\mu < 1``, la suma de ``N``
   copias independientes está dominada por la más grande: compara
   ``\mathbb{E}[\max_i W_i]`` con ``\mathbb{E}[\sum_i W_i]`` para ``N`` grande.
10. Repite el laboratorio de la sección 3 con redes de 2 y de 8 socios en
    promedio. ¿Cómo cambia ``J_c``?
11. Modifica la simulación para que el impuesto se cobre pero **no** se
    redistribuya (la riqueza recaudada desaparece). ¿Cambia la distribución
    de la riqueza **relativa**? Explica por qué.

## Para la próxima sesión

- Lectura: Chakrabarti et al. (2013), sección 4.5; para la ruta avanzada,
  Bouchaud y Mézard (2000).
- En las sesiones 10 y 11 estudiamos las **redes** como objeto en sí mismo:
  grados, centralidad y modelos de redes aleatorias. Hoy vimos que la forma de
  la red cambia el resultado; allá aprenderemos a medirla, primero con la red
  de matrimonios de las familias de Florencia (sesión 10) y después con una red
  construida con las correlaciones de la sesión 7 (sesión 11).
"""

# ╔═╡ 064ec8d6-ca80-11f1-88bc-cdb38180b014
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 064ec8e0-ca80-11f1-a37c-1f91c81ff2e2
begin
    const AZUL    = "#2a78d6"   # serie principal
    const NARANJA = "#eb6834"   # segunda serie o teoría
    const GRIS    = "#52514e"   # referencias y texto
    const RAMPA   = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # azul secuencial
end;

# ╔═╡ 064ec89a-ca80-11f1-9e92-9beb1eaa61e6
let
    marcas = [1e-3, 1e-2, 1e-1]
    p = plot(; xscale = :log10, yscale = :log10, xlabel = "número de agentes N", ylabel = "Y₂ (escala log)",
             xticks = ([1000, 5000, 20000], ["1,000", "5,000", "20,000"]),
             yticks = (marcas, string.(marcas)), ylims = (1e-3, 0.2), legend = :outerright,
             title = "Red con J = 0.1: Y₂ al crecer la población", size = (760, 360))
    φs = sort(unique(barrido_impuesto.φ))
    for (k, φ) in enumerate(φs)
        d = sort(barrido_impuesto[barrido_impuesto.φ .== φ, :], :N)
        plot!(p, d.N, d.Y2; color = φ == 0 ? NARANJA : RAMPA[min(k - 1, 4)], lw = 2, marker = :circle,
              ms = 4, msw = 0, label = φ == 0 ? "Sin impuesto" : "φ = $φ")
    end
    p
end

# ╔═╡ 064ec8ec-ca80-11f1-a7bb-b33654da21eb
"Número entero con comas como separador de miles."
miles(x) = replace(string(round(Int, x)), r"(?<=\d)(?=(\d{3})+$)" => ",")

# ╔═╡ 064ec868-ca80-11f1-a64c-332d6c4c2bad
let
    marcas = [1e-4, 1e-3, 1e-2, 1e-1]
    ejes = (xlabel = "intensidad de intercambio J", yscale = :log10, ylims = (8e-5, 0.4),
            xlims = (0, 0.85), yticks = (marcas, string.(marcas)))
    p1 = plot(; ejes..., ylabel = "Y₂ (escala log)", legend = false,
              title = "Campo medio: intercambio con todos")
    p2 = plot(; ejes..., legend = :bottomleft, title = "Red: 4 socios en promedio")
    vspan!(p2, [0, 0.1]; color = NARANJA, alpha = 0.12, label = "Y₂ no cambia con N")
    for (k, N) in enumerate(sort(unique(barrido.N)))
        for (p, modelo) in ((p1, "campo medio"), (p2, "red"))
            d = sort(barrido[(barrido.modelo .== modelo) .& (barrido.N .== N), :], :J)
            plot!(p, d.J, d.Y2; color = RAMPA[k+1], lw = 2, marker = :circle, ms = 4, msw = 0,
                  label = "N = $(miles(N))")
        end
    end
    plot(p1, p2; layout = (1, 2), size = (920, 380))
end

# ╔═╡ 064ec8ec-ca80-11f1-bbd3-af00aed7d190
"Datos ordenados de mayor a menor y su CCDF empírica, P(X ≥ x)."
function ccdf_empirica(x)
    s = sort(x; rev = true)
    s, (1:length(s)) ./ length(s)
end

# ╔═╡ 064ec8f4-ca80-11f1-b0bb-3bf550cbc3b8
"Índice de Gini de una muestra."
function gini(x)
    s = sort(x)
    n = length(s)
    sum((2 .* (1:n) .- n .- 1) .* s) / (n * sum(s))
end

# ╔═╡ 064ec8f4-ca80-11f1-b7bf-e3116e69b230
"""
Gini de la gamma inversa con exponente μ: es igual al de una gamma con forma
μ - 1 (ver las notas formales de la sesión 9).
"""
gini_bm(μ) = exp(Distributions.SpecialFunctions.loggamma(μ - 0.5) -
                 Distributions.SpecialFunctions.loggamma(μ - 1)) / ((μ - 1) * sqrt(π))

# ╔═╡ 064ec8fe-ca80-11f1-97b5-93d243b871e6
"Red aleatoria con grado promedio c: se emparejan al azar c 'medias aristas' por agente."
function red_aleatoria(N, c; rng = Xoshiro(7))
    vecinos = [Int[] for _ in 1:N]
    extremos = shuffle(rng, repeat(1:N, c))
    for k in 1:2:length(extremos)-1
        a, b = extremos[k], extremos[k+1]
        a != b && (push!(vecinos[a], b); push!(vecinos[b], a))
    end
    for i in 1:N
        isempty(vecinos[i]) && push!(vecinos[i], mod1(i + 1, N))   # nadie queda aislado
    end
    vecinos
end

# ╔═╡ 064ec8fe-ca80-11f1-874b-b3b9bcaf824c
"""
Modelo de Bouchaud-Mézard con σ = 1 y paso dt. En cada paso: rendimiento
aleatorio exacto (sin cambiar la media), intercambio con todos (`vecinos =
nothing`) o con los socios de la red, impuesto `φ` redistribuido a todos, y
normalización para que la riqueza promedio sea 1. Devuelve la riqueza
relativa final; con `Y2 = true` devuelve también el promedio de Y₂ en la
segunda mitad de la simulación.
"""
function simular_bm(N, J; vecinos = nothing, φ = 0.0, dt = 0.01, T = 200.0,
                    rng = Xoshiro(1), Y2 = false)
    w = ones(N)
    nuevo = similar(w)
    pasos = round(Int, T / dt)
    suma_Y2, cuenta = 0.0, 0
    for k in 1:pasos
        w .*= exp.(sqrt(dt) .* randn(rng, N) .- dt / 2)
        if vecinos === nothing
            w .+= J * dt .* (mean(w) .- w)
        else
            for i in 1:N
                s = 0.0
                for j in vecinos[i]
                    s += w[j] - w[i]
                end
                nuevo[i] = w[i] + J * dt * s / length(vecinos[i])
            end
            w .= nuevo
        end
        φ > 0 && (w .+= φ * dt .* (mean(w) .- w))
        w ./= mean(w)
        if Y2 && k > pasos ÷ 2 && k % 10 == 0
            suma_Y2 += sum(w .^ 2) / N^2
            cuenta += 1
        end
    end
    Y2 ? (w, suma_Y2 / cuenta) : w
end

# ╔═╡ 064ec840-ca80-11f1-9b5c-6763d2fbe0a8
riqueza = simular_bm(10_000, J_dist; T = 150.0);

# ╔═╡ 064ec848-ca80-11f1-80da-bfa442083edd
let
    μ = 1 + 2J_dist
    s, P = ccdf_empirica(riqueza)
    p = scatter(s, P; xscale = :log10, yscale = :log10, ms = 1.5, msw = 0, color = AZUL,
                label = "Simulación: 10 000 agentes", legend = :bottomleft, ylims = (1e-4, 1.2),
                xlabel = "riqueza relativa w", ylabel = "P(W ≥ w)",
                title = "Distribución estacionaria con J = $J_dist")
    x = exp10.(range(log10(minimum(s)), log10(maximum(s)), length = 300))
    plot!(p, x, ccdf.(InverseGamma(μ, μ - 1), x); color = NARANJA, lw = 2,
          label = "Teoría: gamma inversa, μ = $(round(μ, digits = 2))")
end

# ╔═╡ 064ec854-ca80-11f1-8921-1979ee9cbb55
md"""
La teoría coincide con la simulación en todo el rango. El Gini de esta
distribución es $(round(gini(riqueza), digits = 3)) en la simulación y
$(round(gini_bm(1 + 2J_dist), digits = 3)) en la teoría. Con ``J = 0.1`` la
cola tiene ``\mu = 1.2``, muy cerca de 1: la riqueza está tan concentrada que
el promedio apenas existe.
"""

# ╔═╡ 064ec872-ca80-11f1-9dc5-1d73eb58303d
let
    w = simular_bm(1000, J_esc; vecinos = red_aleatoria(1000, 4), T = 300.0)
    s = cumsum(sort(w ./ sum(w); rev = true))
    p = plot(1:50, s[1:50]; color = AZUL, lw = 2, seriestype = :steppost,
             xlabel = "número de agentes, del más rico al más pobre",
             ylabel = "fracción acumulada de la riqueza", ylims = (0, 1), xlims = (0, 50),
             label = "Red, J = $J_esc", legend = :bottomright,
             title = "Los 50 agentes más ricos de 1000")
    plot!(p, [0, 50], [0, 0.05]; color = GRIS, ls = :dash, label = "Todos iguales")
    annotate!(p, 25, 0.3, text("Y₂ = $(round(sum((w ./ sum(w)) .^ 2), digits = 3)); " *
                                "el más rico tiene el $(round(Int, 100s[1]))%", 9, GRIS))
end

# ╔═╡ 064ec886-ca80-11f1-bfc8-f944c5234308
con_impuesto = (campo_medio = simular_bm(5_000, 0.1; φ = φ_imp, T = 200.0),
                red = simular_bm(5_000, 0.1; φ = φ_imp, vecinos = red_aleatoria(5_000, 4), T = 300.0));

# ╔═╡ 064ec886-ca80-11f1-bfac-19edbd3efa33
let
    φs = range(0, 0.4, length = 100)
    p1 = plot(φs, gini_bm.(1 .+ 2 .* (0.1 .+ φs)); color = NARANJA, lw = 2, ylims = (0, 1),
              xlabel = "tasa del impuesto φ", ylabel = "índice de Gini", legend = :topright,
              label = "Teoría (campo medio)", title = "Gini con J = 0.1")
    scatter!(p1, [φ_imp], [gini(con_impuesto.campo_medio)]; color = AZUL, ms = 7, msw = 0,
             label = "Simulación, campo medio")
    scatter!(p1, [φ_imp], [gini(con_impuesto.red)]; color = AZUL, ms = 7, msw = 0,
             marker = :diamond, label = "Simulación, red")
    w = con_impuesto.red ./ sum(con_impuesto.red)
    s = cumsum(sort(w; rev = true))
    p2 = plot(1:100, s[1:100]; color = AZUL, lw = 2, seriestype = :steppost, ylims = (0, 1),
              xlabel = "agentes, del más rico al más pobre", ylabel = "fracción acumulada",
              label = "Red con φ = $φ_imp", legend = :bottomright,
              title = "Los 100 más ricos de 5000 (red)")
    plot!(p2, [0, 100], [0, 0.02]; color = GRIS, ls = :dash, label = "Todos iguales")
    plot(p1, p2; layout = (1, 2), size = (920, 360))
end

# ╔═╡ 064ec890-ca80-11f1-b6c0-974cc15c088b
md"""
Con el impuesto $(φ_imp), el Gini en la red es
$(round(gini(con_impuesto.red), digits = 2)) y los 10 agentes más ricos tienen el
$(round(Int, 100 * sum(sort(con_impuesto.red ./ sum(con_impuesto.red); rev = true)[1:10])))%
de la riqueza.

**¿Saca el impuesto a la economía de la fase condensada?** Como en la
sección 3, la prueba es el tamaño. Con `Code/simular_condensacion.jl` se
corrieron simulaciones de la red con ``J = 0.1`` y varios impuestos, de
1,000 a 20,000 agentes (mediana de 4 simulaciones, en
`Data/condensacion_impuesto.csv`):
"""

# ╔═╡ 064ec908-ca80-11f1-a0fe-7b24838a7365
"Fracción de la riqueza del agente más rico a lo largo del tiempo (campo medio)."
function historia_maximo(N, J; T = 100.0, dt = 0.01, rng = Xoshiro(2))
    w = ones(N)
    t, maximo = Float64[], Float64[]
    for k in 1:round(Int, T / dt)
        w .*= exp.(sqrt(dt) .* randn(rng, N) .- dt / 2)
        w .+= J * dt .* (mean(w) .- w)
        w ./= mean(w)
        if k % 20 == 0
            push!(t, k * dt)
            push!(maximo, maximum(w) / N)
        end
    end
    (t = t, maximo = maximo)
end

# ╔═╡ 064ec822-ca80-11f1-9b20-9ff14e211294
let
    p = plot(; yscale = :log10, xlabel = "tiempo", ylabel = "fracción de la riqueza total",
             ylims = (5e-5, 1.5), legend = :bottomright,
             title = "Fracción del agente más rico (N = 1000)")
    for (J, c, nombre) in ((0.0, NARANJA, "Sólo inversión (J = 0)"), (0.5, AZUL, "Inversión e intercambio (J = 0.5)"))
        h = historia_maximo(1000, J; T = 100.0)
        plot!(p, h.t, h.maximo; color = c, lw = 2, label = nombre)
    end
    hline!(p, [1 / 1000]; color = GRIS, ls = :dash, lw = 1, label = "Todos iguales: 1/N")
end

# ╔═╡ 064ec908-ca80-11f1-a77e-b7865037c210
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

# ╔═╡ 064ec8c2-ca80-11f1-a46e-b33b82d16a63
verificar(() -> participacion([1.0, 1.0, 2.0]), 0.375;
          pista = md"Primero divide entre el total: `f = w ./ sum(w)`. Luego suma los cuadrados: `sum(f .^ 2)`.")

# ╔═╡ 064ec8cc-ca80-11f1-b918-bf712a3cdf4d
verificar(() -> exponente_pareto(0.3, 0.1, 2.0), 1.2;
          pista = md"La fórmula es ``\mu = 1 + 2(J + \varphi)/\sigma^2``.")

# ╔═╡ Cell order:
# ╟─064ec7e4-ca80-11f1-8589-cd7f3d5d737c
# ╟─064ec816-ca80-11f1-bc22-3df113ff94d2
# ╟─064ec816-ca80-11f1-84d5-df39c2467675
# ╟─064ec822-ca80-11f1-9b20-9ff14e211294
# ╟─064ec82c-ca80-11f1-b491-f14cf2752876
# ╟─064ec836-ca80-11f1-978d-d395939ffd06
# ╠═064ec840-ca80-11f1-9b5c-6763d2fbe0a8
# ╟─064ec848-ca80-11f1-80da-bfa442083edd
# ╟─064ec854-ca80-11f1-8921-1979ee9cbb55
# ╟─064ec854-ca80-11f1-87e0-9771df17749d
# ╟─064ec85e-ca80-11f1-a519-dd4ddc68e300
# ╠═064ec85e-ca80-11f1-a178-8d908092580f
# ╟─064ec868-ca80-11f1-a64c-332d6c4c2bad
# ╟─064ec872-ca80-11f1-8614-014a54cb5fcc
# ╟─064ec872-ca80-11f1-9dc5-1d73eb58303d
# ╟─064ec87c-ca80-11f1-8198-3f8aa828e2ef
# ╠═064ec886-ca80-11f1-bfc8-f944c5234308
# ╟─064ec886-ca80-11f1-bfac-19edbd3efa33
# ╟─064ec890-ca80-11f1-b6c0-974cc15c088b
# ╠═064ec890-ca80-11f1-82ba-d9ba11b78736
# ╟─064ec89a-ca80-11f1-9e92-9beb1eaa61e6
# ╟─064ec89a-ca80-11f1-a0bf-e5aa1e6e4444
# ╟─064ec8a4-ca80-11f1-a789-a349ee856539
# ╟─064ec8ae-ca80-11f1-a0e9-59bfe2882fce
# ╟─064ec8ae-ca80-11f1-b425-0b8ab50dd798
# ╠═064ec8ba-ca80-11f1-bd5d-3ba0b1ff6df2
# ╟─064ec8c2-ca80-11f1-a46e-b33b82d16a63
# ╟─064ec8c2-ca80-11f1-b0ff-4bfb46c063f1
# ╠═064ec8cc-ca80-11f1-8a71-e95e82753218
# ╟─064ec8cc-ca80-11f1-b918-bf712a3cdf4d
# ╟─064ec8d6-ca80-11f1-9a89-9d864cda096e
# ╟─064ec8d6-ca80-11f1-88bc-cdb38180b014
# ╠═064ec8e0-ca80-11f1-ba30-1b07f57af930
# ╠═064ec8e0-ca80-11f1-a37c-1f91c81ff2e2
# ╠═064ec8ec-ca80-11f1-a7bb-b33654da21eb
# ╠═064ec8ec-ca80-11f1-bbd3-af00aed7d190
# ╠═064ec8f4-ca80-11f1-b0bb-3bf550cbc3b8
# ╠═064ec8f4-ca80-11f1-b7bf-e3116e69b230
# ╠═064ec8fe-ca80-11f1-97b5-93d243b871e6
# ╠═064ec8fe-ca80-11f1-874b-b3b9bcaf824c
# ╠═064ec908-ca80-11f1-a0fe-7b24838a7365
# ╠═064ec908-ca80-11f1-a77e-b7865037c210
