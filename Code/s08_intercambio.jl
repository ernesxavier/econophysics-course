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

# ╔═╡ bb4677f4-ca92-11f1-92f7-6911629eb229
begin
    # Usa el entorno del curso (Code/Project.toml) en lugar del gestor de
    # paquetes propio de Pluto: un solo entorno para todo el repositorio.
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, CSV, DataFrames, Statistics, StatsBase, Distributions, Random, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

# ╔═╡ bb467696-ca92-11f1-8683-d943e5aa14ec
md"""
# Sesión 8: Desigualdad y modelos de intercambio

**Objetivos de la sesión**
1. Medir la desigualdad del ingreso en México con la **curva de Lorenz** y el
   **índice de Gini**, usando la ENIGH 2024.
2. Entender por qué, si el dinero se conserva y se intercambia al azar, su
   distribución es **exponencial**, igual que la energía en un gas.
3. Ver cómo el **ahorro** cambia la forma de la distribución: ahorro
   uniforme (gamma) y ahorro heterogéneo (cola de Pareto).
4. Comparar los tres modelos con los datos de México en el laboratorio.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s08_intercambio.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): P.B3 (histogramas), C.B1
(exponenciales), C.B6 (la integral como área), F.B1–B3 (de lo micro a lo
macro, entropía, temperatura), F.I1 (distribución de Boltzmann), F.I6
(modelos de intercambio). Ruta avanzada: C.A1 (multiplicadores de Lagrange),
C.A4 (función gamma).

**Conceptos transversales.** Se retoman las **distribuciones, la CCDF y las
escalas logarítmicas** (sesiones 1 y 3).
"""

# ╔═╡ bb4676dc-ca92-11f1-bc48-11e123ee30ca
TableOfContents(title = "Contenido")

# ╔═╡ bb4676e8-ca92-11f1-857b-bf3d27d5247c
md"""
## 0. El fenómeno: ¿cómo se reparte el ingreso en México?

La Encuesta Nacional de Ingresos y Gastos de los Hogares (ENIGH) 2024 de
INEGI entrevistó a 91,414 hogares. Cada uno tiene un **factor de expansión**:
el número de hogares del país que representa. Todas las estadísticas de esta
sesión están ponderadas con él. Los datos están en
`Data/enigh2024_hogares.csv` y se regeneran con `Code/descargar_enigh.jl`.

Usamos el **ingreso corriente mensual**: trabajo, negocios, rentas,
transferencias y el valor estimado de vivir en casa propia, en pesos de 2024.

Unidad de análisis: $(@bind unidad Select(["Hogar", "Persona (ingreso per cápita)"]))
"""

# ╔═╡ bb4676f0-ca92-11f1-8281-85e4839bce42
enigh = CSV.read(joinpath(@__DIR__, "..", "Data", "enigh2024_hogares.csv"), DataFrame)

# ╔═╡ bb4676f0-ca92-11f1-94f5-2764f2d1b659
ingreso, peso = if unidad == "Hogar"
    enigh.ingreso, Float64.(enigh.factor)
else
    # cada hogar cuenta tantas veces como integrantes tiene
    enigh.ingreso ./ enigh.integrantes, Float64.(enigh.factor .* enigh.integrantes)
end;

# ╔═╡ bb46771a-ca92-11f1-af8e-8b073d1a8b73
md"""
## 1. Curva de Lorenz e índice de Gini

Ordena a todos de menor a mayor ingreso. La **curva de Lorenz** ``L(p)`` es la
fracción del ingreso total que recibe la fracción ``p`` más pobre. Si todos
ganaran lo mismo, ``L(p) = p``: la diagonal. Mientras más se aleja la curva
de la diagonal, más desigual es la distribución.

El **índice de Gini** es el doble del área entre la diagonal y la curva
(repaso C.B6):

```math
G = 1 - 2\int_0^1 L(p)\, dp .
```

Vale 0 con igualdad perfecta y se acerca a 1 si una sola unidad recibe todo.
"""

# ╔═╡ bb467722-ca92-11f1-b4cf-1b64cde4898d
md"""
Con hogares como unidad, el Gini de México es de alrededor de 0.40; con
personas, un poco mayor, porque los hogares pobres tienen más integrantes.
Para comparar con los modelos usaremos ingresos divididos entre su media:
el Gini no cambia si se multiplica a todos por la misma cantidad.
"""

# ╔═╡ bb46772c-ca92-11f1-9c98-5d0792fb3542
md"""
## 2. Dinero conservado: la distribución exponencial

**El argumento verbal.** Imagina ``N`` agentes que empiezan con el mismo
dinero. En cada transacción, dos agentes al azar juntan su dinero y se lo
reparten con una proporción aleatoria. El dinero total **se conserva**, como
la energía en los choques entre moléculas de un gas. Para un gas sabemos la
respuesta (repaso F.I1): la energía se reparte según la **distribución de
Boltzmann**, ``p(E) \propto e^{-E/k_B T}``. Por la misma razón, el dinero
debería terminar con distribución **exponencial**,

```math
p(m) = \frac{1}{T}\, e^{-m/T}, \qquad T = \langle m \rangle = \frac{M}{N},
```

donde la "temperatura" ``T`` es el dinero promedio (Dragulescu y Yakovenko,
2000). La razón profunda es la **entropía**: de todas las formas de
repartir ``M`` pesos entre ``N`` agentes, la exponencial es la que se
obtiene de más maneras (sección 6).

**Comprobémoslo numéricamente** con 5,000 agentes que empiezan con 1 peso
cada uno.

Regla de intercambio: $(@bind regla Select(["Reparto aleatorio", "Monto fijo"]))
  Transacciones por agente: $(@bind t_dy Select([0, 1, 2, 5, 20, 100, 500], default = 100))
"""

# ╔═╡ bb467740-ca92-11f1-8ad3-a5037d90e95f
md"""
## 3. Ahorro uniforme: la distribución gamma

En el modelo anterior un agente puede perder casi todo en una sola
transacción. Chakraborti y Chakrabarti (2000) agregaron **ahorro**: cada
agente guarda una fracción ``\lambda`` de su dinero y sólo intercambia el
resto,

```math
\begin{aligned}
m_i' &= \lambda\, m_i + \varepsilon\,(1-\lambda)(m_i + m_j), \\
m_j' &= \lambda\, m_j + (1-\varepsilon)(1-\lambda)(m_i + m_j),
\end{aligned}
```

con ``\varepsilon`` uniforme en ``[0, 1]``. El dinero total se sigue
conservando. Con ``\lambda = 0`` se recupera el modelo de la sección 2.

Fracción ahorrada ``\lambda``: $(@bind λ_cc Slider(0:0.05:0.9, default = 0.5, show_value = true))
"""

# ╔═╡ bb46775e-ca92-11f1-bd6e-f7666d30201d
md"""
## 4. Ahorro heterogéneo: la cola de Pareto

¿Y si cada agente ahorra una fracción **distinta** ``\lambda_i``? Chatterjee,
Chakrabarti y Manna (2004) mostraron que entonces aparece una **cola de
Pareto**: los agentes que ahorran casi todo (``\lambda_i`` cercano a 1)
acumulan fortunas muy grandes, porque casi no arriesgan nada en cada
transacción.

Asignamos los ahorros con ``\lambda_i = 1 - U_i^{1/a}``, con ``U_i`` uniforme. El
parámetro ``a`` controla cuántos agentes ahorran casi todo: con ``a = 1`` los
``\lambda_i`` son uniformes en ``[0, 1]``; con ``a`` mayor, hay menos agentes
cerca de 1.

Parámetro ``a``: $(@bind a_het Slider([1.0, 1.5, 2.0, 2.5, 3.0], default = 1.0, show_value = true))
"""

# ╔═╡ bb467790-ca92-11f1-91f0-cd384a2ed288
md"""
**Resumen numérico.** Cada columna se compara con la de la ENIGH (todas con
media 1):
"""

# ╔═╡ bb4677a4-ca92-11f1-8a56-635b4cf1f73d
md"""
## 6. [Avanzado] Máxima entropía

**El problema.** Buscamos la densidad ``p(m)`` en ``m \ge 0`` con **máxima
entropía** (repaso F.B2),

```math
S[p] = -\int_0^\infty p(m) \ln p(m)\, dm ,
```

sujeta a dos restricciones: ``p`` está normalizada, ``\int_0^\infty p\, dm = 1``,
y el dinero promedio es ``T``, ``\int_0^\infty m\, p\, dm = T``.

**Multiplicadores de Lagrange** (repaso C.A1). Igualamos a cero la derivada
de ``S`` menos ``\alpha`` por la primera restricción y ``\beta`` por la
segunda, con respecto al valor de ``p`` en cada punto:

```math
-\ln p(m) - 1 - \alpha - \beta m = 0
\quad\Longrightarrow\quad
p(m) = e^{-1-\alpha}\, e^{-\beta m}.
```

Las restricciones fijan ``e^{-1-\alpha} = \beta`` y ``\beta = 1/T``: es la
exponencial, con entropía máxima ``S = 1 + \ln T``.

**La simulación sube la entropía.** La gráfica muestra la entropía de la
distribución simulada (estimada con un histograma) a lo largo del tiempo, en
el modelo sin ahorro con ``T = 1``. Parte de un valor muy bajo (todos con lo
mismo) y se acerca a la máxima posible, ``S = 1``: es la **segunda ley** en
versión económica.
"""

# ╔═╡ bb4677b8-ca92-11f1-8cf5-05507f1e0f03
md"""
**[Avanzado] La aproximación gamma.** La fórmula ``n = 1 + 3\lambda/(1-\lambda)``
de la sección 3 sale de igualar la **varianza** de la distribución simulada
con la de una gamma. Para comprobarla, ajustamos una gamma por máxima
verosimilitud a simulaciones con varios ``\lambda``:
"""

# ╔═╡ bb4677cc-ca92-11f1-982d-3fb1c3f4dc7b
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `gini_simple(x)`, que calcule el índice de Gini de un vector
de ingresos (sin pesos) con la fórmula

```math
G = \frac{\sum_{i=1}^{n}\sum_{j=1}^{n} \lvert x_i - x_j\rvert}{2 n^2 \bar x} .
```
"""

# ╔═╡ bb4677cc-ca92-11f1-966c-d76caf51e546
function gini_simple(x)
    missing  # reemplaza esta línea
end

# ╔═╡ bb4677d6-ca92-11f1-bad9-2d58708083bf
md"""
**2.** *Si alcanza (margen de la clase).* Escribe `intercambio(mi, mj, λ, ε)`, que devuelva el dinero de los dos
agentes después de una transacción con ahorro ``\lambda`` (sección 3), como
un vector `[mi_nuevo, mj_nuevo]`.
"""

# ╔═╡ bb4677e0-ca92-11f1-9779-e139ee1b4028
function intercambio(mi, mj, λ, ε)
    missing  # reemplaza esta línea
end

# ╔═╡ bb4677ea-ca92-11f1-a585-6980e252bdd8
md"""
### Ruta básica: interpretación

3. En el modelo sin ahorro nadie es más hábil que otro. ¿Qué nos dice el
   modelo sobre la desigualdad que observamos? ¿Qué **no** nos dice?
4. Explica con tus palabras por qué ahorrar reduce la desigualdad en el
   modelo de la sección 3.
5. ¿Por qué las encuestas de hogares como la ENIGH subestiman la cola de la
   distribución? ¿Qué otra fuente de datos usarías para medirla?


### Conceptos que vuelven

Preguntas cortas para repasar ideas de sesiones anteriores; se responden con
lo que ya sabes y con el laboratorio de hoy.

- **R1.** En la CCDF del ingreso de la sección 5, ¿cómo distingues una cola
  exponencial de una cola de Pareto? ¿En qué escala se ve recta cada una?
  (Sesiones 1 y 3.)

### Ruta avanzada

6. Demuestra que el Gini de una exponencial es 1/2 y que el de una gamma
   con forma ``n`` es ``\Gamma(n + 1/2) / \big(n\, \Gamma(n) \sqrt{\pi}\big)``.
7. Calcula la varianza estacionaria del dinero en el modelo con ahorro
   uniforme: escribe ``\mathbb{E}[m_i'^2]`` en términos de ``\mathbb{E}[m^2]`` y
   ``\mathbb{E}[m]^2``, impón estacionariedad y deduce ``n = 1 + 3\lambda/(1-\lambda)``.
8. Calcula el Gini de la ENIGH con intervalos de confianza bootstrap
   remuestreando hogares (repaso S.I1). ¿Es estadísticamente distinto el
   Gini de las localidades rurales (`tam_loc == 4`) del de las ciudades
   grandes (`tam_loc == 1`)?

## Para la próxima sesión

- Lectura: Yakovenko y Rosser (2009). Complementaria: Requião da Cunha
  (2021), sección 8.3; Slanina (2014), sección 7.2. Para la ruta avanzada,
  Chakrabarti et al. (2013), capítulos 2 a 5; de consulta, Garibaldi y Scalas
  (2010).
- En la sesión 9 agregamos crecimiento multiplicativo al intercambio (modelo
  de Bouchaud-Mézard): veremos que la riqueza puede **condensarse** en unos
  pocos agentes, y probaremos impuestos y redistribución.
"""

# ╔═╡ bb4677f4-ca92-11f1-b6c6-f148f1458e57
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ bb4677fe-ca92-11f1-91e9-3b4dc8897e8c
begin
    const AZUL    = "#2a78d6"   # serie principal
    const NARANJA = "#eb6834"   # segunda serie o teoría
    const AQUA    = "#1baf7a"   # tercera serie
    const GRIS    = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # datos de fondo
end;

# ╔═╡ bb467808-ca92-11f1-933e-0f634c54bc95
"Número entero con comas como separador de miles."
miles(x) = replace(string(round(Int, x)), r"(?<=\d)(?=(\d{3})+$)" => ",")

# ╔═╡ bb467704-ca92-11f1-9c4f-57c88d0d9a1f
let
    m = mean(ingreso, Weights(peso))
    q = quantile(ingreso, Weights(peso), [0.5, 0.9, 0.99])
    lim = quantile(ingreso, Weights(peso), 0.98)
    bordes = range(0, lim, length = 81)
    p = histogram(ingreso; weights = peso, bins = bordes, normalize = :probability,
                  color = AZUL, lw = 0, alpha = 0.8, xlims = (0, lim),
                  xlabel = "ingreso corriente mensual (pesos)",
                  ylabel = "fracción por intervalo de $(miles(step(bordes))) pesos",
                  xformatter = x -> miles(x), legend = :topright,
                  title = "Ingreso mensual por $(lowercase(first(unidad, 7))): ENIGH 2024 (98% más bajo)",
                  label = "")
    vline!(p, [q[1]]; color = GRIS, lw = 2, ls = :dash, label = "Mediana: $(miles(q[1]))")
    vline!(p, [m]; color = NARANJA, lw = 2, label = "Media: $(miles(m))")
end

# ╔═╡ bb46770e-ca92-11f1-a95b-e93110c00b28
md"""
La distribución es muy asimétrica: la media queda a la derecha de la mediana
porque unos pocos ingresos muy altos la jalan. El 10% de mayores ingresos
empieza en $(miles(quantile(ingreso, Weights(peso), 0.9))) pesos al mes y el
1% más alto en $(miles(quantile(ingreso, Weights(peso), 0.99))). Hoy
buscamos un mecanismo que explique esta forma.
"""

# ╔═╡ bb467808-ca92-11f1-929e-45b9f436fb30
"Curva de Lorenz ponderada: fracción acumulada de la población y del ingreso."
function lorenz(x, w)
    o = sortperm(x)
    W, S = cumsum(w[o]), cumsum(x[o] .* w[o])
    [0; W ./ W[end]], [0; S ./ S[end]]
end

# ╔═╡ bb467808-ca92-11f1-9cb1-19a1bff46546
"Índice de Gini con pesos (por defecto, todos iguales): 1 menos el doble del área bajo la curva de Lorenz."
function gini(x, w = ones(length(x)))
    p, L = lorenz(x, w)
    1 - sum(diff(p) .* (L[2:end] .+ L[1:end-1]))
end

# ╔═╡ bb46771a-ca92-11f1-84f2-cd2a708351d7
let
    p, L = lorenz(ingreso, peso)
    G = gini(ingreso, peso)
    fig = plot([0, 1], [0, 1]; color = GRIS, ls = :dash, lw = 1.5, label = "Igualdad perfecta",
               xlabel = "fracción de la población, de menor a mayor ingreso",
               ylabel = "fracción del ingreso total", aspect_ratio = :equal,
               xlims = (0, 1), ylims = (0, 1), legend = :topleft, size = (520, 500),
               title = "Curva de Lorenz, ENIGH 2024: Gini = $(round(G, digits = 3))")
    plot!(fig, p, L; fillrange = p, fillcolor = AZUL, fillalpha = 0.15, color = AZUL, lw = 2,
          label = "México")
    i = findfirst(>=(0.5), p)
    scatter!(fig, [p[i]], [L[i]]; color = AZUL, ms = 5, msw = 0, label = "")
    annotate!(fig, p[i] + 0.02, L[i] - 0.03,
              text("el 50% más pobre recibe\nel $(round(Int, 100L[i]))% del ingreso", 8, :left, GRIS))
end

# ╔═╡ bb467812-ca92-11f1-92b1-6b062133d363
"Gini de una gamma con forma n."
gini_gamma(n) = exp(Distributions.SpecialFunctions.loggamma(n + 0.5) -
                    Distributions.SpecialFunctions.loggamma(n)) / (n * sqrt(π))

# ╔═╡ bb46781e-ca92-11f1-bcae-77d2e4c8be32
"Datos ordenados de mayor a menor y su CCDF empírica, P(X ≥ x)."
function ccdf_empirica(x)
    s = sort(x; rev = true)
    s, (1:length(s)) ./ length(s)
end

# ╔═╡ bb46781e-ca92-11f1-aaa5-9dcecf67f626
"CCDF con pesos: P(X ≥ x) usando el factor de expansión."
function ccdf_ponderada(x, w)
    o = sortperm(x; rev = true)
    x[o], cumsum(w[o]) ./ sum(w)
end

# ╔═╡ bb467826-ca92-11f1-ba5c-6b8b85a0b295
"Estimador de Hill del exponente de la cola, con la fracción `q` de datos más grandes."
function hill_cola(x, q)
    s = sort(x; rev = true)
    n = round(Int, q * length(s))
    n / sum(log.(s[1:n] ./ s[n]))
end

# ╔═╡ bb467826-ca92-11f1-9ece-23a228f16f2d
"""
Modelo de Dragulescu-Yakovenko sin ahorro: `t` transacciones por agente.
"Reparto aleatorio": los dos juntan su dinero y lo reparten con proporción
uniforme. "Monto fijo": uno le da 0.1 al otro, si lo tiene.
"""
function simular_dy(regla, t; N = 5_000, rng = Xoshiro(1))
    m = ones(N)
    for _ in 1:round(Int, t * N)
        i, j = rand(rng, 1:N), rand(rng, 1:N)
        i == j && continue
        if regla == "Reparto aleatorio"
            s = m[i] + m[j]
            ε = rand(rng)
            m[i], m[j] = ε * s, (1 - ε) * s
        elseif m[i] >= 0.1
            m[i] -= 0.1
            m[j] += 0.1
        end
    end
    m
end

# ╔═╡ bb467736-ca92-11f1-ab0b-edb805ba317f
dinero_dy = simular_dy(regla, t_dy; N = 5_000);

# ╔═╡ bb467736-ca92-11f1-a840-67f771eb6d17
let
    m = dinero_dy
    bordes = 0:0.1:6
    x = range(0, 6, length = 200)
    p1 = histogram(m; bins = bordes, normalize = :pdf, color = AZUL, lw = 0, alpha = 0.7,
                   xlabel = "dinero m (promedio = 1)", ylabel = "densidad", xlims = (0, 6),
                   ylims = (0, 1.3), label = "Simulación", legend = :topright,
                   title = "$t_dy transacciones por agente")
    plot!(p1, x, exp.(-x); color = NARANJA, lw = 2, label = "Exponencial: e⁻ᵐ")
    p2 = plot(; yscale = :log10, xlabel = "dinero m", ylabel = "densidad (escala log)",
              xlims = (0, 6), ylims = (1e-3, 2), title = "Escala semilog: una recta")
    if t_dy > 0
        h = fit(Histogram, m, bordes)
        d = h.weights ./ (length(m) * step(bordes))
        k = d .> 0
        scatter!(p2, (bordes[1:end-1] .+ 0.05)[k], d[k]; color = AZUL, ms = 4, msw = 0)
    end
    plot!(p2, x, exp.(-x); color = NARANJA, lw = 2)
    plot(p1, p2; layout = (1, 2), size = (900, 340))
end

# ╔═╡ bb467740-ca92-11f1-9236-0357511d771c
md"""
Al principio todos tienen 1 peso (Gini = 0). Después de unas pocas
transacciones por agente la distribución ya es exponencial, con
Gini = $(round(gini(dinero_dy), digits = 3)) en esta simulación; el de una
exponencial es exactamente ``1/2``. Con la regla de **monto fijo** (en cada
encuentro uno le da 0.1 pesos al otro, si los tiene) se llega a la misma
distribución, aunque más despacio: prueba con 500 transacciones por agente.
El detalle microscópico no importa, sólo la conservación y el azar.

Ojo con la interpretación: nadie en el modelo es más hábil ni trabaja más.
La desigualdad aparece **sólo por azar**, y es considerable.
"""

# ╔═╡ bb467830-ca92-11f1-afb2-b9ec1b88696f
"Modelo con ahorro: el agente i guarda la fracción λ[i] de su dinero. `pasos` transacciones por agente."
function simular_ahorro(λ; pasos = 500, rng = Xoshiro(2))
    N = length(λ)
    m = ones(N)
    for _ in 1:pasos*N
        i, j = rand(rng, 1:N), rand(rng, 1:N)
        i == j && continue
        ε = rand(rng)
        s = (1 - λ[i]) * m[i] + (1 - λ[j]) * m[j]
        m[i], m[j] = λ[i] * m[i] + ε * s, λ[j] * m[j] + (1 - ε) * s
    end
    m
end

# ╔═╡ bb467748-ca92-11f1-926c-7dabaeb29911
dinero_cc = simular_ahorro(fill(λ_cc, 5_000); pasos = 500);

# ╔═╡ bb467754-ca92-11f1-a352-a5d84930e0c0
let
    n = 1 + 3λ_cc / (1 - λ_cc)
    x = range(0.001, 4, length = 300)
    p1 = histogram(dinero_cc; bins = 0:0.05:4, normalize = :pdf, color = AZUL, lw = 0, alpha = 0.7,
                   xlabel = "dinero m (promedio = 1)", ylabel = "densidad", xlims = (0, 4),
                   label = "Simulación", legend = :topright, title = "Ahorro λ = $λ_cc")
    plot!(p1, x, pdf.(Gamma(n, 1 / n), x); color = NARANJA, lw = 2,
          label = "Gamma, n = $(round(n, digits = 2))")
    λs = 0:0.05:0.9
    p2 = plot(λs, gini_gamma.(1 .+ 3λs ./ (1 .- λs)); color = NARANJA, lw = 2, xlims = (0, 0.9),
              ylims = (0, 0.55), xlabel = "fracción ahorrada λ", ylabel = "índice de Gini",
              label = "Gamma de Patriarca", legend = :topright, title = "El ahorro reduce la desigualdad")
    scatter!(p2, [λ_cc], [gini(dinero_cc)]; color = AZUL, ms = 7, msw = 0, label = "Esta simulación")
    plot(p1, p2; layout = (1, 2), size = (900, 340))
end

# ╔═╡ bb467754-ca92-11f1-b4f8-d71db90c4419
md"""
Con ahorro, casi nadie se queda con muy poco: la densidad se anula en
``m = 0`` y tiene un máximo en un valor positivo. Patriarca, Chakraborti y
Kaski (2004) encontraron que la distribución es muy cercana a una **gamma**,

```math
p(m) = \frac{n^n}{\Gamma(n)}\, m^{n-1} e^{-n m},
\qquad n = 1 + \frac{3\lambda}{1-\lambda},
```

(con dinero promedio 1). Para ``\lambda = 0``, ``n = 1`` y es la exponencial. Con
``\lambda`` = $(λ_cc), el Gini baja a $(round(gini(dinero_cc), digits = 3)).
Una interpretación física: ``n`` es como el número de "grados de libertad"
entre los que cada agente reparte su dinero.
"""

# ╔═╡ bb46783a-ca92-11f1-ae92-e3bfbda89f9c
"Ahorros heterogéneos λᵢ = 1 - Uᵢ^(1/a), con Uᵢ uniforme."
ahorros_heterogeneos(a; N = 5_000, rng = Xoshiro(3)) = 1 .- rand(rng, N) .^ (1 / a)

# ╔═╡ bb467768-ca92-11f1-8f60-a3e97e934a28
dinero_het = simular_ahorro(ahorros_heterogeneos(a_het; N = 5_000); pasos = 5_000);

# ╔═╡ bb467768-ca92-11f1-96f7-8fc37d61185e
let
    s, P = ccdf_empirica(dinero_het)
    p = plot(s, P; xscale = :log10, yscale = :log10, color = AZUL, lw = 2,
             label = "Ahorro heterogéneo (a = $a_het)", xlabel = "dinero m (promedio = 1)",
             ylabel = "P(M ≥ m)", ylims = (1 / 5_000, 1.2), xlims = (1e-2, maximum(s) * 1.5),
             legend = :bottomleft, title = "Distribución acumulada complementaria (log-log)")
    x = exp10.(range(-2, log10(maximum(s)), length = 200))
    plot!(p, x, exp.(-x); color = NARANJA, lw = 2, ls = :dash, label = "Exponencial (sin ahorro)")
    α = hill_cola(dinero_het, 0.05)
    xm = s[round(Int, 0.05 * length(s))]
    xr = exp10.(range(log10(xm), log10(maximum(s)), length = 50))
    plot!(p, xr, 0.05 .* (xr ./ xm) .^ (-α); color = GRIS, lw = 1.5, ls = :dot,
          label = "Pareto con α = $(round(α, digits = 2)) (5% más rico)")
end

# ╔═╡ bb467772-ca92-11f1-911f-6b73490de096
md"""
La cola es casi una recta en log-log: una ley de potencias. Con ``a = 1`` su
exponente es cercano a 1, la ley de Zipf. Al aumentar ``a`` la cola se vuelve
más delgada y el Gini baja: $(round(gini(dinero_het), digits = 3)) en esta
simulación. Fíjate en el rango del eje horizontal: con ``a = 1`` el agente
más rico tiene cientos de veces el dinero promedio.
"""

# ╔═╡ bb46783a-ca92-11f1-b2ec-9fa6fdbb7d28
"Ajustes a la ENIGH (ingreso / media): gamma por máxima verosimilitud y distancias KS."
ajuste_enigh = let
    k = ingreso .> 0
    x, w = ingreso[k], peso[k]
    u = x ./ mean(x, Weights(w))
    g = fit_mle(Gamma, u, w)
    n = shape(g)
    ln = LogNormal(mean(log.(u), Weights(w)), std(log.(u), Weights(w)))
    o = sortperm(u)
    F = cumsum(w[o]) ./ sum(w)
    ks(d) = maximum(abs.(F .- cdf.(d, u[o])))
    (n = n, λ = (n - 1) / (n + 2), ks_gamma = ks(Gamma(n, 1 / n)),
     ks_exp = ks(Exponential(1)), ks_lognormal = ks(ln))
end

# ╔═╡ bb46777a-ca92-11f1-af21-4b0501fcb4c9
md"""
## 5. Laboratorio: los modelos frente a la ENIGH

Ahora comparamos los tres modelos con el ingreso de México, dividido entre
su media para que todos tengan promedio 1. Para el ahorro uniforme usamos el
``\lambda`` que mejor ajusta los datos: la gamma de máxima verosimilitud para
la ENIGH tiene forma ``n = `` $(round(ajuste_enigh.n, digits = 2)), que corresponde
a ``\lambda = (n-1)/(n+2) = `` $(round(ajuste_enigh.λ, digits = 2)).

Parámetro ``a`` del ahorro heterogéneo: $(@bind a_lab Slider([1.0, 1.5, 2.0, 2.5, 3.0], default = 2.5, show_value = true))
"""

# ╔═╡ bb46777a-ca92-11f1-8a25-f39b9bacdfcc
modelos = let
    λ = ajuste_enigh.λ
    (dy = simular_dy("Reparto aleatorio", 500; N = 5_000),
     cc = simular_ahorro(fill(λ, 5_000); pasos = 500),
     het = simular_ahorro(ahorros_heterogeneos(a_lab; N = 5_000); pasos = 5_000))
end;

# ╔═╡ bb467786-ca92-11f1-a92c-69c23328eee8
let
    x = ingreso ./ mean(ingreso, Weights(peso))
    series = [(modelos.dy, "Sin ahorro (exponencial)", AZUL),
              (modelos.cc, "Ahorro uniforme, λ = $(round(ajuste_enigh.λ, digits = 2))", NARANJA),
              (modelos.het, "Ahorro heterogéneo, a = $a_lab", AQUA)]
    # Izquierda: densidad en escala lineal, donde se ve la parte baja
    bordes = 0:0.1:4
    densidad(v, w) = fit(Histogram, v, Weights(w), bordes).weights ./ (sum(w) * step(bordes))
    p1 = bar(bordes[1:end-1] .+ 0.05, densidad(x, peso); bar_width = 0.1, color = GRIS_CLARO, lw = 0,
             label = "ENIGH 2024", xlims = (0, 4), legend = false, xlabel = "ingreso / media",
             ylabel = "densidad", title = "El cuerpo (escala lineal)")
    for (m, nombre, c) in series
        plot!(p1, bordes, [densidad(m, ones(length(m))); 0]; seriestype = :steppost, color = c, lw = 2)
    end
    # Derecha: CCDF en log-log, donde se ve la cola
    s, P = ccdf_ponderada(x, peso)
    k = 1:max(1, length(s) ÷ 3000):length(s)            # menos puntos para graficar
    p2 = plot(; xscale = :log10, yscale = :log10, xlims = (1e-2, 100), ylims = (1e-4, 1.2),
              legend = :bottomleft, xlabel = "ingreso / media", ylabel = "P(X ≥ x)",
              title = "La cola (log-log)")
    for (m, nombre, c) in series
        sm, Pm = ccdf_empirica(m)
        plot!(p2, sm, Pm; color = c, lw = 2, label = nombre)
    end
    plot!(p2, s[k], P[k]; color = :black, lw = 2.5, ls = :dash, label = "ENIGH 2024")
    plot(p1, p2; layout = (1, 2), size = (920, 380))
end

# ╔═╡ bb46779a-ca92-11f1-9f66-bd7c6f7e6b91
md"""
Lo que dicen los datos:

- El modelo **sin ahorro** exagera la desigualdad en la parte baja: predice
  demasiados hogares con ingreso cercano a cero y un Gini de 0.5.
- El **ahorro uniforme** reproduce bien el cuerpo de la distribución y el
  Gini, pero su cola cae demasiado rápido: casi no hay ingresos de 5 veces la
  media.
- Con hogares como unidad, el **ahorro heterogéneo** con ``a = 2.5`` reproduce la cola de la encuesta,
  pero deja a los más pobres con menos de lo que tienen en los datos (mira
  el percentil 10). Con ``a = 1`` la cola es mucho más pesada.

Ninguno de los tres ajusta todo a la vez: un solo mecanismo no basta para
explicar la distribución completa.

Dos advertencias. Primero, estos modelos describen **riqueza** que se
intercambia, y la ENIGH mide **ingreso**: la analogía es útil, pero no es la
misma variable. Segundo, las encuestas captan mal los ingresos más altos (los
hogares más ricos casi no aparecen o declaran de menos), así que la cola real
de México es más pesada que la de la ENIGH.

**[Avanzado]** Una lognormal, el resultado del crecimiento multiplicativo de
la sesión 6, ajusta el cuerpo de la ENIGH aún mejor: la distancia de
Kolmogorov-Smirnov es $(round(ajuste_enigh.ks_lognormal, digits = 3)), contra
$(round(ajuste_enigh.ks_gamma, digits = 3)) de la gamma y
$(round(ajuste_enigh.ks_exp, digits = 3)) de la exponencial (repaso S.I3). Que
un modelo ajuste no prueba que su mecanismo sea el correcto.
"""

# ╔═╡ bb4677c2-ca92-11f1-966b-63191c8c5ae4
md"""
## Resumen

| Modelo | Regla | Distribución estacionaria | Gini |
|---|---|---|---|
| Sin ahorro (Dragulescu-Yakovenko) | Juntar y repartir al azar | Exponencial: máxima entropía | 0.5 |
| Ahorro uniforme (Chakraborti-Chakrabarti) | Cada quien guarda ``\lambda`` | Gamma, ``n = 1 + 3\lambda/(1-\lambda)`` | Baja con ``\lambda`` |
| Ahorro heterogéneo (Chatterjee-Chakrabarti-Manna) | Cada quien guarda su ``\lambda_i`` | Cuerpo gamma y cola de Pareto | Depende de la distribución de ``\lambda_i`` |
| **ENIGH 2024** | | Cuerpo cercano a una gamma con ``\lambda \approx`` $(round(ajuste_enigh.λ, digits = 2)), cola más pesada | $(round(gini(enigh.ingreso, Float64.(enigh.factor)), digits = 2)) (hogares) |
"""

# ╔═╡ bb467844-ca92-11f1-a478-a5998847542d
"Tabla con estadísticas de la ENIGH y de los tres modelos, todos con media 1."
function tabla_comparacion(ingreso, peso, modelos, λ, a)
    W = Weights(peso)
    x = ingreso ./ mean(ingreso, W)
    fila(nombre, v, w) = (nombre, gini(v, w), quantile(v, Weights(w), 0.5),
                          quantile(v, Weights(w), 0.1), quantile(v, Weights(w), 0.99))
    unos(v) = ones(length(v))
    filas = [fila("ENIGH 2024", x, peso),
             fila("Sin ahorro", modelos.dy, unos(modelos.dy)),
             fila("Ahorro uniforme (λ = $(round(λ, digits = 2)))", modelos.cc, unos(modelos.cc)),
             fila("Ahorro heterogéneo (a = $a)", modelos.het, unos(modelos.het))]
    r(v) = round(v, digits = 2)
    filas_md = join(["| $(f[1]) | $(r(f[2])) | $(r(f[3])) | $(r(f[4])) | $(r(f[5])) |" for f in filas], "\n")
    Markdown.parse("""
    | | Gini | Mediana | Percentil 10 | Percentil 99 |
    |---|---|---|---|---|
    $filas_md
    """)
end

# ╔═╡ bb467790-ca92-11f1-9d8f-fb9c5c235eb6
tabla_comparacion(ingreso, peso, modelos, ajuste_enigh.λ, a_lab)

# ╔═╡ bb467850-ca92-11f1-8a1d-57e947932f33
"Entropía diferencial estimada con un histograma, en varios momentos del modelo sin ahorro."
function entropia_en_el_tiempo(; N = 5_000, tiempos = [1, 10, 100], rng = Xoshiro(4))
    m = ones(N)
    hecho = 0.0
    S = Float64[]
    Δ = 0.05
    for t in tiempos
        for _ in 1:round(Int, (t - hecho) * N)
            i, j = rand(rng, 1:N), rand(rng, 1:N)
            i == j && continue
            s = m[i] + m[j]; ε = rand(rng)
            m[i], m[j] = ε * s, (1 - ε) * s
        end
        hecho = t
        h = fit(Histogram, m, 0:Δ:maximum(m)+Δ)
        p = h.weights ./ (N * Δ)
        push!(S, -sum(pk * log(pk) * Δ for pk in p if pk > 0))
    end
    (t = tiempos, S = S)
end

# ╔═╡ bb4677ac-ca92-11f1-99e2-973f73dead30
let
    e = entropia_en_el_tiempo(; N = 5_000, tiempos = [0.25, 0.5, 1, 2, 3, 5, 10, 20, 50])
    plot(e.t, e.S; xscale = :log10, color = AZUL, lw = 2, marker = :circle, ms = 5, msw = 0,
         label = "Entropía de la simulación", xlabel = "transacciones por agente",
         ylabel = "entropía S", legend = :bottomright, ylims = (-1, 1.2),
         title = "El intercambio aleatorio maximiza la entropía")
    hline!([1]; color = NARANJA, lw = 2, ls = :dash, label = "Máximo: exponencial, S = 1")
end

# ╔═╡ bb467858-ca92-11f1-9f32-7508b2229db7
"Forma de la gamma ajustada por máxima verosimilitud para varios valores de λ."
function barrido_forma_gamma()
    λs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    (λ = λs, n = [shape(fit_mle(Gamma, simular_ahorro(fill(λ, 5_000); pasos = 500))) for λ in λs])
end

# ╔═╡ bb4677b8-ca92-11f1-be9d-0b771bc1eef4
let
    b = barrido_forma_gamma()
    λs = range(0, 0.9, length = 200)
    p = plot(λs, 1 .+ 3λs ./ (1 .- λs); color = NARANJA, lw = 2, yscale = :log10,
             label = "n = 1 + 3λ/(1-λ)", xlabel = "fracción ahorrada λ", ylabel = "forma n de la gamma",
             legend = :topleft, title = "Forma de la distribución con ahorro uniforme")
    scatter!(p, b.λ, b.n; color = AZUL, ms = 6, msw = 0, label = "Ajuste por máxima verosimilitud")
end

# ╔═╡ bb467862-ca92-11f1-ad6a-6bc15f3710b4
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

# ╔═╡ bb4677d6-ca92-11f1-af6d-61f4809f1072
verificar(() -> gini_simple([1.0, 2.0, 3.0, 4.0, 10.0]), 0.4;
          pista = md"Puedes usar una comprensión doble: `sum(abs(a - b) for a in x, b in x)`, y dividir entre `2 * length(x)^2 * mean(x)`.")

# ╔═╡ bb4677e0-ca92-11f1-a428-716f601978f4
verificar(() -> intercambio(3.0, 1.0, 0.5, 0.1), [1.7, 2.3];
          pista = md"Calcula primero lo que se pone en juego, `s = (1 - λ) * (mi + mj)`; luego cada agente recibe lo que guardó (`λ * mi`) más su parte de `s`.")

# ╔═╡ Cell order:
# ╟─bb467696-ca92-11f1-8683-d943e5aa14ec
# ╟─bb4676dc-ca92-11f1-bc48-11e123ee30ca
# ╟─bb4676e8-ca92-11f1-857b-bf3d27d5247c
# ╠═bb4676f0-ca92-11f1-8281-85e4839bce42
# ╠═bb4676f0-ca92-11f1-94f5-2764f2d1b659
# ╟─bb467704-ca92-11f1-9c4f-57c88d0d9a1f
# ╟─bb46770e-ca92-11f1-a95b-e93110c00b28
# ╟─bb46771a-ca92-11f1-af8e-8b073d1a8b73
# ╟─bb46771a-ca92-11f1-84f2-cd2a708351d7
# ╟─bb467722-ca92-11f1-b4cf-1b64cde4898d
# ╟─bb46772c-ca92-11f1-9c98-5d0792fb3542
# ╠═bb467736-ca92-11f1-ab0b-edb805ba317f
# ╟─bb467736-ca92-11f1-a840-67f771eb6d17
# ╟─bb467740-ca92-11f1-9236-0357511d771c
# ╟─bb467740-ca92-11f1-8ad3-a5037d90e95f
# ╠═bb467748-ca92-11f1-926c-7dabaeb29911
# ╟─bb467754-ca92-11f1-a352-a5d84930e0c0
# ╟─bb467754-ca92-11f1-b4f8-d71db90c4419
# ╟─bb46775e-ca92-11f1-bd6e-f7666d30201d
# ╠═bb467768-ca92-11f1-8f60-a3e97e934a28
# ╟─bb467768-ca92-11f1-96f7-8fc37d61185e
# ╟─bb467772-ca92-11f1-911f-6b73490de096
# ╟─bb46777a-ca92-11f1-af21-4b0501fcb4c9
# ╠═bb46777a-ca92-11f1-8a25-f39b9bacdfcc
# ╟─bb467786-ca92-11f1-a92c-69c23328eee8
# ╟─bb467790-ca92-11f1-91f0-cd384a2ed288
# ╟─bb467790-ca92-11f1-9d8f-fb9c5c235eb6
# ╟─bb46779a-ca92-11f1-9f66-bd7c6f7e6b91
# ╟─bb4677a4-ca92-11f1-8a56-635b4cf1f73d
# ╟─bb4677ac-ca92-11f1-99e2-973f73dead30
# ╟─bb4677b8-ca92-11f1-8cf5-05507f1e0f03
# ╟─bb4677b8-ca92-11f1-be9d-0b771bc1eef4
# ╟─bb4677c2-ca92-11f1-966b-63191c8c5ae4
# ╟─bb4677cc-ca92-11f1-982d-3fb1c3f4dc7b
# ╠═bb4677cc-ca92-11f1-966c-d76caf51e546
# ╟─bb4677d6-ca92-11f1-af6d-61f4809f1072
# ╟─bb4677d6-ca92-11f1-bad9-2d58708083bf
# ╠═bb4677e0-ca92-11f1-9779-e139ee1b4028
# ╟─bb4677e0-ca92-11f1-a428-716f601978f4
# ╟─bb4677ea-ca92-11f1-a585-6980e252bdd8
# ╟─bb4677f4-ca92-11f1-b6c6-f148f1458e57
# ╠═bb4677f4-ca92-11f1-92f7-6911629eb229
# ╠═bb4677fe-ca92-11f1-91e9-3b4dc8897e8c
# ╠═bb467808-ca92-11f1-933e-0f634c54bc95
# ╠═bb467808-ca92-11f1-929e-45b9f436fb30
# ╠═bb467808-ca92-11f1-9cb1-19a1bff46546
# ╠═bb467812-ca92-11f1-92b1-6b062133d363
# ╠═bb46781e-ca92-11f1-bcae-77d2e4c8be32
# ╠═bb46781e-ca92-11f1-aaa5-9dcecf67f626
# ╠═bb467826-ca92-11f1-ba5c-6b8b85a0b295
# ╠═bb467826-ca92-11f1-9ece-23a228f16f2d
# ╠═bb467830-ca92-11f1-afb2-b9ec1b88696f
# ╠═bb46783a-ca92-11f1-ae92-e3bfbda89f9c
# ╠═bb46783a-ca92-11f1-b2ec-9fa6fdbb7d28
# ╠═bb467844-ca92-11f1-a478-a5998847542d
# ╠═bb467850-ca92-11f1-8a1d-57e947932f33
# ╠═bb467858-ca92-11f1-9f32-7508b2229db7
# ╠═bb467862-ca92-11f1-ad6a-6bc15f3710b4
