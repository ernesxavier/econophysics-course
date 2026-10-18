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

# ╔═╡ 7db70cda-ca99-11f1-8c6d-1f44cb6bfb1a
begin
    # Usa el entorno del curso (Code/Project.toml) en lugar del gestor de
    # paquetes propio de Pluto: un solo entorno para todo el repositorio.
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, CSV, DataFrames, Statistics, StatsBase, Distributions,
          Random, Graphs, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

# ╔═╡ 7db70b92-ca99-11f1-be11-598adb5d183f
md"""
# Taller del Módulo II: desigualdad y contagio

**Qué se evalúa.** Que sepas usar las herramientas de las sesiones 8 a 14
(curva de Lorenz, Gini, CCDF, modelos de intercambio, redes y cascadas) en un
caso nuevo, y que sepas **interpretar** lo que sale. No se piden
derivaciones.

**Tu caso.** El profesor te asigna una entidad federativa. En la parte 1
estudias la distribución del ingreso de sus hogares con la ENIGH 2024
(`Data/enigh2024_hogares.csv`), y en la parte 2 la comparas con un modelo de
intercambio. La parte 3 es común a todos: cascadas en redes simuladas, con
una pregunta que no vimos en clase, **a quién conviene proteger**.

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
| 1. La desigualdad en tu entidad | 10 | 20 | 30 |
| 2. ¿La explica un modelo de intercambio? | 10 | 20 | 30 |
| 3. ¿A quién proteger del contagio? | 10 | 30 | 40 |

Una respuesta completa dice **qué** se observa, **con qué número o
gráfica** y **por qué** importa. La ruta avanzada, al final, es opcional y
cuenta como extra.

**En el semestre** hay una sesión de taller en clase: en los primeros 20
minutos se resuelven dudas del módulo, después se trabaja en parejas en las
partes 1 y 3, y al final cada pareja comenta una respuesta con el grupo. En el
trimestre el taller se resuelve en casa.

**Repaso recomendado** (en `Class_notes/repaso/`): P.B1 y P.B3
(distribuciones e histogramas), C.B2 (escalas logarítmicas), F.B1–B3 (de lo
micro a lo macro, entropía y temperatura), L.B4 (matriz de adyacencia), E.B3
(promediar simulaciones). Para la ruta avanzada: F.I3–I4 (Ising y campo
medio) y C.I5 (puntos fijos).
"""

# ╔═╡ 7db70be0-ca99-11f1-8548-4d2f73814b79
TableOfContents(title = "Contenido")

# ╔═╡ 7db70bea-ca99-11f1-9549-47938cb3db28
md"""
## 1. Parte 1: la desigualdad en tu entidad

Cada hogar de la ENIGH representa a muchos hogares del país: el **factor de
expansión** dice a cuántos. Todas las medidas de esta parte usan ese factor
como peso (sesión 8).
"""

# ╔═╡ 7db70bf2-ca99-11f1-af27-43263d7baa23
enigh = CSV.read(joinpath(@__DIR__, "..", "Data", "enigh2024_hogares.csv"), DataFrame);

# ╔═╡ 7db70c08-ca99-11f1-af3c-db9bac6e2676
md"""
### Ejercicio 1 (10 puntos)

Escribe `participacion_top(x, w, q)`: la fracción del ingreso total que
reciben los hogares más ricos que juntos son la fracción `q` de la población,
contando cada hogar con su peso `w`. Por ejemplo, con `q = 0.1` es la
participación del 10% más rico.

Para simplificar, ordena los hogares de mayor a menor ingreso y toma los
primeros hasta que sus pesos acumulados alcancen `q` del peso total (incluye
el hogar con el que se alcanza o rebasa `q`).
"""

# ╔═╡ 7db70c08-ca99-11f1-a06e-cb71e41dc855
function participacion_top(x, w, q)
    missing  # reemplaza esta línea
end

# ╔═╡ 7db70c24-ca99-11f1-bef4-fdfcc6468358
md"""
**Pregunta 1** (10 puntos). ¿Qué tan desigual es el ingreso en tu entidad
comparado con el país? Usa el Gini y la participación del 10% y del 1% más
ricos. ¿Dice lo mismo la curva de Lorenz que el Gini? ¿Por qué dos entidades
con el mismo Gini podrían tener curvas de Lorenz distintas?
"""

# ╔═╡ 7db70c30-ca99-11f1-b6a8-f95402e32f24
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ 7db70c3a-ca99-11f1-8901-2bf3f09c3053
md"""
Ahora la **forma** de la distribución. Dividimos el ingreso entre su media
ponderada y dibujamos la CCDF en dos escalas: en **semilog** una exponencial
se ve recta, y en **log-log** se ve recta una ley de Pareto (sesiones 3 y 8).
"""

# ╔═╡ 7db70c44-ca99-11f1-b537-190a6584f26b
md"""
**Pregunta 2** (10 puntos). ¿En qué parte de la distribución se parece tu
entidad a una exponencial, y en qué parte se separa? ¿Hacia dónde se separa:
hay más o menos hogares muy ricos que los que predice la exponencial? ¿En
qué escala se ve mejor cada cosa y por qué?
"""

# ╔═╡ 7db70c44-ca99-11f1-9670-71c4849ac597
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ 7db70c4e-ca99-11f1-87d8-8940e485d01f
md"""
## 2. Parte 2: ¿la explica un modelo de intercambio?

En el modelo de Chakraborti y Chakrabarti (sesión 8), cada agente guarda la
fracción ``\lambda`` de su dinero y pone el resto en un intercambio al azar.
La distribución estacionaria es casi una gamma con forma

```math
n = \frac{1 + 2\lambda}{1 - \lambda} ,
```

cuyo Gini depende sólo de ``n``. Con ``\lambda = 0`` es la exponencial
(``n = 1``, Gini = 0.5); al ahorrar más, ``n`` crece y el Gini baja.

### Ejercicio 2 (10 puntos)

Escribe `lambda_para_gini(G)`: el ahorro ``\lambda`` entre 0 y 0.99 con el
que el modelo tiene Gini igual a `G`. Usa **bisección**: el Gini del modelo,
`gini_gamma((1 + 2λ) / (1 - λ))`, baja cuando ``\lambda`` sube. Repite 60
veces: toma el punto medio del intervalo y quédate con la mitad donde está la
solución.
"""

# ╔═╡ 7db70c4e-ca99-11f1-86bc-15e3c2e108d6
function lambda_para_gini(G)
    missing  # reemplaza esta línea
end

# ╔═╡ 7db70c6c-ca99-11f1-8040-6b2521b2ec04
md"""
**Pregunta 3** (10 puntos). Calibramos el modelo para que tenga el mismo
Gini que tu entidad. ¿Reproduce también el cuerpo de la distribución (la
mediana y el percentil 10)? ¿Y la cola (el percentil 99 y el máximo)? Si
coincide el Gini pero no la cola, ¿qué te dice eso sobre resumir la
desigualdad con un solo número?
"""

# ╔═╡ 7db70c6c-ca99-11f1-b7a4-45e63be46bb9
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ 7db70c76-ca99-11f1-9ce4-5f64b0eff02d
md"""
**Pregunta 4** (10 puntos). ¿Qué significa en tu entidad el valor de
``\lambda`` que encontraste? Menciona dos razones por las que el modelo de
intercambio no puede ser la explicación completa de la distribución del
ingreso (piensa en qué mide la ENIGH, en el trabajo y en la herencia), y una
razón por la que de todos modos es útil (sesiones 8 y 9).
"""

# ╔═╡ 7db70c76-ca99-11f1-9f99-cb9524652177
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ 7db70c80-ca99-11f1-ad86-8bcbb50b0075
md"""
## 3. Parte 3: ¿a quién proteger del contagio?

Usamos el modelo de umbral de Watts de la sesión 12: un nodo cae cuando la
fracción de sus vecinos caídos es al menos ``\varphi``. En una red de bancos,
``\varphi = K/\theta`` es el capital dividido entre los préstamos
interbancarios.

Un regulador puede **proteger** a algunos bancos (exigirles tanto capital que
no caen nunca). ¿A quién conviene proteger: a los más conectados o a bancos
al azar? Y si cae un banco, ¿importa cuál?

### Ejercicio 3 (5 puntos)

Escribe `mas_conectados(grados, m)`: los índices de los `m` nodos con mayor
grado (en cualquier orden).
"""

# ╔═╡ 7db70c80-ca99-11f1-bdc6-1b7821a8ebf8
function mas_conectados(grados, m)
    missing  # reemplaza esta línea
end

# ╔═╡ 7db70c8a-ca99-11f1-ae73-f11b1e95c48f
md"""
### Ejercicio 4 (5 puntos)

Escribe `fraccion_vulnerables(grados, φ)`: la fracción de nodos que caen con
**un solo** vecino caído, es decir, con grado ``k \ge 1`` y
``1/k \ge \varphi``.
"""

# ╔═╡ 7db70c96-ca99-11f1-8dc9-7fed5f7ad10c
function fraccion_vulnerables(grados, φ)
    missing  # reemplaza esta línea
end

# ╔═╡ 7db70c9e-ca99-11f1-a2c1-2b912b2da96a
md"""
**El experimento.** Para cada configuración simulamos 100 redes de 2,000
nodos con umbral ``\varphi = 0.18``. En cada red comparamos:

- **¿Quién cae primero?** Un nodo al azar contra el nodo de mayor grado.
- **¿A quién se protege?** Una fracción de los nodos, elegidos al azar o entre
  los más conectados. Después cae un nodo al azar (no protegido).

Tipo de red: $(@bind red_t2 Select(["Erdős-Rényi", "Barabási-Albert"]))

Grado promedio ``z``: $(@bind z_t2 Select([2, 4, 6, 8], default = 4))
"""

# ╔═╡ 7db70cb2-ca99-11f1-bb49-9987b621a2b5
md"""
**Pregunta 5** (10 puntos). ¿Importa quién cae primero? Compara un nodo al
azar con el de mayor grado en las dos redes con ``z = 4``. ¿En cuál de las
dos redes la diferencia es más grande y por qué? Explica qué les pasa a los
vecinos del nodo de mayor grado cuando éste cae.
"""

# ╔═╡ 7db70cbc-ca99-11f1-9999-13e23f32a422
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ 7db70cbc-ca99-11f1-af47-c3145c896e45
md"""
**Pregunta 6** (10 puntos). En una epidemia, vacunar a los más conectados es
mucho más eficaz que vacunar al azar. ¿Pasa lo mismo aquí? Compara las dos
estrategias de protección en Erdős-Rényi con ``z = 4`` y con ``z = 2``, y
después en Barabási-Albert. Explica en qué casos conviene proteger a los más
conectados usando la idea de nodo **vulnerable** (ejercicio 4): ¿son
vulnerables los nodos más conectados? ¿Qué papel tienen en la cascada si nadie
los protege, cuando la red está muy conectada y cuando apenas lo está?
"""

# ╔═╡ 7db70cc8-ca99-11f1-8b21-87fe79a2cf49
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ 7db70cc8-ca99-11f1-820d-e3e7e5672402
md"""
**Pregunta 7** (10 puntos). Escribe una recomendación de cinco a ocho líneas
para un regulador bancario con base en las preguntas 5 y 6. ¿Dónde pondrías
la vigilancia y dónde los requisitos de capital? Menciona una limitación del
modelo que podría cambiar tu recomendación (por ejemplo, que en la red real
los bancos grandes también son los que más prestan, sesión 12).
"""

# ╔═╡ 7db70cd0-ca99-11f1-9bc3-694f119f9aac
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ 7db70cd0-ca99-11f1-847a-0da4c3783373
md"""
## 4. Ruta avanzada (extra, hasta 10 puntos)

Haz **dos** de los cuatro. Agrega las celdas que necesites al final del
cuaderno.

- **A1. La cola de tu entidad.** Estima el exponente de Pareto de la cola con
  el estimador de Hill (`hill_cola` está en el apéndice) usando el 1%, 2% y 5%
  más rico. ¿Es estable? La ENIGH subestima los ingresos más altos (los
  hogares muy ricos casi no responden): ¿hacia dónde sesga eso tu estimación?
- **A2. Ahorro heterogéneo.** Repite la parte 2 con ahorros distintos para
  cada agente, ``\lambda_i = 1 - U_i^{1/a}`` (`ahorros_heterogeneos`, sesión
  8). ¿Aparece una cola de Pareto? ¿Se parece más a tu entidad?
- **A3. La ventana de cascadas.** Para Barabási-Albert, barre ``z`` de 2 a 12
  con un nodo al azar como semilla y compara la ventana con la de
  Erdős-Rényi. Explica la diferencia con la fracción de vulnerables.
- **A4. Histéresis en el campo medio.** Resuelve ``m = \tanh(\beta(Jm + h))``
  iterando desde un valor inicial, con ``\beta J = 1.5``. Sube ``h`` de -1 a 1
  usando como punto de partida la solución anterior, y después bájalo de 1 a
  -1. Dibuja las dos curvas. ¿Qué tienen en común con una cascada global: un
  cambio pequeño que produce un cambio grande (sesiones 13 y 14)?
"""

# ╔═╡ 7db70cda-ca99-11f1-bb32-dd7a6e47e176
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para resolver el taller.
"""

# ╔═╡ 7db70ce4-ca99-11f1-8e5a-5bc92af01953
begin
    const AZUL    = "#2a78d6"   # serie principal
    const NARANJA = "#eb6834"   # segunda serie
    const GRIS    = "#52514e"   # referencias
end;

# ╔═╡ 7db70cee-ca99-11f1-a5b4-73a597a989a9
"Nombres de las entidades federativas, en el orden de su clave del INEGI."
const ENTIDADES = ["Aguascalientes", "Baja California", "Baja California Sur", "Campeche",
                   "Coahuila", "Colima", "Chiapas", "Chihuahua", "Ciudad de México", "Durango",
                   "Guanajuato", "Guerrero", "Hidalgo", "Jalisco", "Estado de México",
                   "Michoacán", "Morelos", "Nayarit", "Nuevo León", "Oaxaca", "Puebla",
                   "Querétaro", "Quintana Roo", "San Luis Potosí", "Sinaloa", "Sonora",
                   "Tabasco", "Tamaulipas", "Tlaxcala", "Veracruz", "Yucatán", "Zacatecas"];

# ╔═╡ 7db70bf2-ca99-11f1-a127-351330b93036
md"""
**Tu entidad:** $(@bind entidad Select(collect(1:32) .=> ENTIDADES, default = 1))
"""

# ╔═╡ 7db70bfe-ca99-11f1-8b73-a77797fc6a7a
begin
    hogares = filter(:entidad => ==(entidad), enigh)
    x_ent, w_ent = hogares.ingreso, Float64.(hogares.factor)
    x_nac, w_nac = enigh.ingreso, Float64.(enigh.factor)
end;

# ╔═╡ 7db70cee-ca99-11f1-8cc6-2be3f0d3e03a
miles(x) = replace(string(round(Int, x)), r"(?<=\d)(?=(\d{3})+$)" => ",")

# ╔═╡ 7db70cfa-ca99-11f1-bbe7-15d3fa44abba
"Curva de Lorenz ponderada: fracción acumulada de la población y del ingreso."
function lorenz(x, w)
    o = sortperm(x)
    W, S = cumsum(w[o]), cumsum(x[o] .* w[o])
    [0; W ./ W[end]], [0; S ./ S[end]]
end

# ╔═╡ 7db70c1c-ca99-11f1-8c11-ad0361929a12
let
    pe, Le = lorenz(x_ent, w_ent)
    pn, Ln = lorenz(x_nac, w_nac)
    p = plot([0, 1], [0, 1]; color = GRIS, ls = :dash, label = "Igualdad perfecta",
             xlabel = "fracción de los hogares (de menor a mayor ingreso)",
             ylabel = "fracción del ingreso", legend = :topleft, size = (480, 420),
             aspect_ratio = 1, xlims = (0, 1), ylims = (0, 1), title = "Curva de Lorenz")
    plot!(p, pn, Ln; color = GRIS, lw = 1.5, label = "Nacional")
    plot!(p, pe, Le; color = AZUL, lw = 2.5, label = ENTIDADES[entidad])
end

# ╔═╡ 7db70cfa-ca99-11f1-8b29-d3eacdf2e5a1
"Índice de Gini con pesos (por defecto, todos iguales): 1 menos el doble del área bajo la curva de Lorenz."
function gini(x, w = ones(length(x)))
    p, L = lorenz(x, w)
    1 - sum(diff(p) .* (L[2:end] .+ L[1:end-1]))
end

# ╔═╡ 7db70d02-ca99-11f1-b167-9597daac4068
"Participación en el ingreso de la fracción `q` más rica (versión de referencia)."
function participacion_top_ref(x, w, q)
    o = sortperm(x; rev = true)
    n = findfirst(>=(q), cumsum(w[o]) ./ sum(w))
    sum(x[o][1:n] .* w[o][1:n]) / sum(x .* w)
end

# ╔═╡ 7db70c24-ca99-11f1-98dc-2d961441b910
let
    fila(nombre, x, w) = "| $nombre | $(length(x)) | $(round(gini(x, w), digits = 3)) | " *
        "$(miles(quantile(x, Weights(w), 0.5))) | $(round(100participacion_top_ref(x, w, 0.1), digits = 1))% | " *
        "$(round(100participacion_top_ref(x, w, 0.01), digits = 1))% |"
    ginis = [gini(g.ingreso, Float64.(g.factor)) for g in groupby(enigh, :entidad; sort = true)]
    lugar = count(>(gini(x_ent, w_ent)), ginis) + 1
    Markdown.parse("""
    | | Hogares en la muestra | Gini | Mediana mensual (pesos) | Ingreso del 10% más rico | Ingreso del 1% más rico |
    |---|---:|---:|---:|---:|---:|
    $(fila("**" * ENTIDADES[entidad] * "**", x_ent, w_ent))
    $(fila("Nacional", x_nac, w_nac))

    Entre las 32 entidades, la tuya es la número **$lugar** de mayor a menor
    Gini (van de $(round(minimum(ginis), digits = 3)) a $(round(maximum(ginis), digits = 3))).
    """)
end

# ╔═╡ 7db70d02-ca99-11f1-a66d-7b47cb5bf6f8
"Gini de una gamma con forma n."
gini_gamma(n) = exp(Distributions.SpecialFunctions.loggamma(n + 0.5) -
                    Distributions.SpecialFunctions.loggamma(n)) / (n * sqrt(π))

# ╔═╡ 7db70d0c-ca99-11f1-8d93-b98ed2ef8f95
"Ahorro λ con el que el modelo de intercambio tiene Gini `G` (versión de referencia)."
function lambda_para_gini_ref(G)
    a, b = 0.0, 0.99
    for _ in 1:60
        c = (a + b) / 2
        gini_gamma((1 + 2c) / (1 - c)) > G ? (a = c) : (b = c)
    end
    (a + b) / 2
end

# ╔═╡ 7db70d0c-ca99-11f1-aeec-7fb4d776be67
"Datos ordenados de mayor a menor y su CCDF empírica, P(X ≥ x)."
function ccdf_empirica(x)
    s = sort(x; rev = true)
    s, (1:length(s)) ./ length(s)
end

# ╔═╡ 7db70d0c-ca99-11f1-919a-3788476e4ff6
"CCDF con pesos: P(X ≥ x) usando el factor de expansión."
function ccdf_ponderada(x, w)
    o = sortperm(x; rev = true)
    x[o], cumsum(w[o]) ./ sum(w)
end

# ╔═╡ 7db70c3a-ca99-11f1-842c-23a21e9bda10
let
    m = mean(x_ent, Weights(w_ent))
    k = x_ent .> 0
    xs, P = ccdf_ponderada(x_ent[k] ./ m, w_ent[k])
    u = range(0, 6, length = 100)
    p1 = plot(xs, P; color = AZUL, lw = 2, yscale = :log10, xlims = (0, 6), ylims = (1e-3, 1.05),
              label = ENTIDADES[entidad], xlabel = "ingreso / media", ylabel = "P(X ≥ x)",
              legend = :topright, title = "Semilog: ¿cuerpo exponencial?")
    plot!(p1, u, exp.(-u); color = GRIS, ls = :dash, label = "Exponencial con media 1")
    p2 = plot(xs, P; color = AZUL, lw = 2, xscale = :log10, yscale = :log10,
              xlims = (0.05, maximum(xs) * 1.2), ylims = (minimum(P) / 2, 1.05),
              label = ENTIDADES[entidad], xlabel = "ingreso / media", ylabel = "P(X ≥ x)",
              legend = :bottomleft, title = "Log-log: ¿cola de Pareto?")
    plot!(p2, u[2:end], exp.(-u[2:end]); color = GRIS, ls = :dash, label = "Exponencial con media 1")
    plot(p1, p2; layout = (1, 2), size = (780, 340))
end

# ╔═╡ 7db70d16-ca99-11f1-a12a-d99f99c579aa
"Estimador de Hill del exponente de la cola, con la fracción `q` de datos más grandes."
function hill_cola(x, q)
    s = sort(x; rev = true)
    n = round(Int, q * length(s))
    n / sum(log.(s[1:n] ./ s[n]))
end

# ╔═╡ 7db70d16-ca99-11f1-a035-2fa2f243110b
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

# ╔═╡ 7db70c58-ca99-11f1-8f9a-55b7597be0ad
modelo_ent = let
    λ = lambda_para_gini_ref(gini(x_ent, w_ent))
    (λ = λ, m = simular_ahorro(fill(λ, 5_000)))
end;

# ╔═╡ 7db70c62-ca99-11f1-8996-6d837d64f631
let
    m = mean(x_ent, Weights(w_ent))
    k = x_ent .> 0
    xs, P = ccdf_ponderada(x_ent[k] ./ m, w_ent[k])
    xm, Pm = ccdf_empirica(modelo_ent.m ./ mean(modelo_ent.m))
    etiqueta = "Modelo, λ = $(round(modelo_ent.λ, digits = 2))"
    p1 = plot(xs, P; color = AZUL, lw = 2, yscale = :log10, xlims = (0, 6), ylims = (1e-3, 1.05),
              label = ENTIDADES[entidad], xlabel = "ingreso / media", ylabel = "P(X ≥ x)",
              legend = :topright, title = "Semilog")
    plot!(p1, xm, Pm; color = NARANJA, lw = 2, ls = :dash, label = etiqueta)
    p2 = plot(xs, P; color = AZUL, lw = 2, xscale = :log10, yscale = :log10,
              xlims = (0.05, maximum(xs) * 1.2), ylims = (minimum(P) / 2, 1.05),
              label = ENTIDADES[entidad], xlabel = "ingreso / media", ylabel = "P(X ≥ x)",
              legend = :bottomleft, title = "Log-log")
    plot!(p2, xm, Pm; color = NARANJA, lw = 2, ls = :dash, label = etiqueta)
    plot(p1, p2; layout = (1, 2), size = (780, 340))
end

# ╔═╡ 7db70c62-ca99-11f1-b50d-69bf1edaeb6f
let
    W = Weights(w_ent)
    u = x_ent ./ mean(x_ent, W)
    mm = modelo_ent.m ./ mean(modelo_ent.m)
    f(v) = round(v, digits = 2)
    Markdown.parse("""
    Con ``\\lambda = `` $(f(modelo_ent.λ)), el modelo tiene el mismo Gini que tu
    entidad. Comparación (ingreso entre la media):

    | | Gini | Percentil 10 | Mediana | Percentil 99 | Máximo |
    |---|---:|---:|---:|---:|---:|
    | **$(ENTIDADES[entidad])** | $(f(gini(x_ent, w_ent))) | $(f(quantile(u, W, 0.1))) | $(f(quantile(u, W, 0.5))) | $(f(quantile(u, W, 0.99))) | $(f(maximum(u))) |
    | Modelo con ahorro | $(f(gini(mm))) | $(f(quantile(mm, 0.1))) | $(f(quantile(mm, 0.5))) | $(f(quantile(mm, 0.99))) | $(f(maximum(mm))) |
    """)
end

# ╔═╡ 7db70d20-ca99-11f1-8455-b519d54c2b90
"Ahorros heterogéneos λᵢ = 1 - Uᵢ^(1/a), con Uᵢ uniforme."
ahorros_heterogeneos(a; N = 5_000, rng = Xoshiro(3)) = 1 .- rand(rng, N) .^ (1 / a)

# ╔═╡ 7db70d20-ca99-11f1-8f23-d3d354127634
"Fracción de nodos vulnerables (versión de referencia)."
fraccion_vulnerables_ref(grados, φ) = count(k -> k >= 1 && 1 / k >= φ, grados) / length(grados)

# ╔═╡ 7db70d28-ca99-11f1-8b27-912053711850
"""
Cascada de umbral a partir de `semilla`: un nodo cae cuando la fracción de
sus vecinos caídos es al menos `φ`. Los nodos `protegidos` no caen nunca.
Devuelve la ronda en que cayó cada nodo (-1 si no cayó).
"""
function rondas_cascada(g, φ, semilla; protegidos = falses(nv(g)))
    ronda = fill(-1, nv(g))
    ronda[semilla] = 0
    caidos = zeros(Int, nv(g))
    frente = [semilla]
    r = 0
    while !isempty(frente)
        r += 1
        nuevo = Int[]
        for v in frente, u in neighbors(g, v)
            (ronda[u] >= 0 || protegidos[u]) && continue
            caidos[u] += 1
            if caidos[u] >= φ * degree(g, u)
                ronda[u] = r
                push!(nuevo, u)
            end
        end
        frente = nuevo
    end
    ronda
end

# ╔═╡ 7db70d34-ca99-11f1-8a2f-0f88a57f3c92
"""
Experimento de la parte 3: en `corridas` redes de `n` nodos, la probabilidad
de una cascada global (más del 5%) y su tamaño promedio, según quién cae
primero y a quién se protege.
"""
function experimento_proteccion(red, z; n = 2000, φ = 0.18, corridas = 100,
                                fracciones = [0.0, 0.025, 0.05, 0.1, 0.15, 0.2])
    generar(s) = red == "Erdős-Rényi" ? erdos_renyi(n, round(Int, n * z / 2); seed = s) :
                                         barabasi_albert(n, z ÷ 2; seed = s)
    resumen(t) = (p = mean(t .> 0.05), t = mean(t))
    tam(ronda) = count(>=(0), ronda) / n
    azar, hub = Float64[], Float64[]
    top = [Float64[] for _ in fracciones]
    alea = [Float64[] for _ in fracciones]
    k_ejemplo = Int[]
    for s in 1:corridas
        g = generar(s)
        k = degree(g)
        s == 1 && (k_ejemplo = k)
        rng = Xoshiro(s)
        push!(azar, tam(rondas_cascada(g, φ, rand(rng, 1:n))))
        push!(hub, tam(rondas_cascada(g, φ, argmax(k))))
        orden = sortperm(k; rev = true)
        mezcla = randperm(rng, n)
        for (i, f) in enumerate(fracciones)
            m = round(Int, f * n)
            for (lista, elegidos) in ((top[i], orden[1:m]), (alea[i], mezcla[1:m]))
                prot = falses(n)
                prot[elegidos] .= true
                semilla = rand(rng, findall(!, prot))
                push!(lista, tam(rondas_cascada(g, φ, semilla; protegidos = prot)))
            end
        end
    end
    (azar = resumen(azar), hub = resumen(hub), fracciones = fracciones,
     prot_top = (p = [resumen(t).p for t in top], t = [resumen(t).t for t in top]),
     prot_azar = (p = [resumen(t).p for t in alea], t = [resumen(t).t for t in alea]),
     grados_ejemplo = k_ejemplo)
end

# ╔═╡ 7db70c9e-ca99-11f1-8839-297befae210f
experimento_t2 = experimento_proteccion(red_t2, z_t2);

# ╔═╡ 7db70ca8-ca99-11f1-8f0a-3fb7367ebee3
let
    e = experimento_t2
    k = e.grados_ejemplo
    Markdown.parse("""
    **$red_t2 con ``z = `` $z_t2.** En una de las redes, el nodo de mayor grado
    tiene $(maximum(k)) vecinos y el $(round(100fraccion_vulnerables_ref(k, 0.18), digits = 1))%
    de los nodos es vulnerable.

    | Primero cae... | Probabilidad de cascada global (más del 5%) | Tamaño promedio |
    |---|---:|---:|
    | Un nodo al azar | $(round(e.azar.p, digits = 2)) | $(round(100e.azar.t, digits = 1))% |
    | El nodo de mayor grado | $(round(e.hub.p, digits = 2)) | $(round(100e.hub.t, digits = 1))% |
    """)
end

# ╔═╡ 7db70ca8-ca99-11f1-9283-9b056f24e1ee
let
    e = experimento_t2
    fr = 100 .* e.fracciones
    p1 = plot(fr, e.prot_azar.p; color = AZUL, lw = 2, marker = :circle, ms = 4, msw = 0,
              label = "Protegidos al azar", xlabel = "nodos protegidos (%)",
              ylabel = "probabilidad de cascada global", ylims = (0, 1.05), legend = :topright,
              title = "Probabilidad")
    plot!(p1, fr, e.prot_top.p; color = NARANJA, lw = 2, marker = :diamond, ms = 4, msw = 0,
          label = "Protegidos los más conectados")
    p2 = plot(fr, 100 .* e.prot_azar.t; color = AZUL, lw = 2, marker = :circle, ms = 4, msw = 0,
              label = "Protegidos al azar", xlabel = "nodos protegidos (%)",
              ylabel = "tamaño promedio (% de la red)", ylims = (0, 105), legend = :topright,
              title = "Tamaño")
    plot!(p2, fr, 100 .* e.prot_top.t; color = NARANJA, lw = 2, marker = :diamond, ms = 4, msw = 0,
          label = "Protegidos los más conectados")
    plot(p1, p2; layout = (1, 2), size = (780, 330),
         plot_title = "$red_t2, z = $z_t2, φ = 0.18", plot_titlefontsize = 11)
end

# ╔═╡ 7db70d3e-ca99-11f1-837e-55657c87fbdd
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

# ╔═╡ 7db70c12-ca99-11f1-8c41-735ad79d3288
verificar(() -> [participacion_top([1.0, 2.0, 3.0, 4.0], [1.0, 1.0, 1.0, 1.0], 0.25),
                 participacion_top([1.0, 2.0, 3.0, 4.0], [1.0, 1.0, 1.0, 1.0], 0.5),
                 participacion_top([10.0, 1.0], [1.0, 9.0], 0.1)],
          [0.4, 0.7, 10 / 19];
          pista = md"Con `o = sortperm(x; rev = true)`, acumula `cumsum(w[o]) ./ sum(w)`; `n = findfirst(>=(q), acumulado)`; la respuesta es `sum(x[o][1:n] .* w[o][1:n]) / sum(x .* w)`.")

# ╔═╡ 7db70c58-ca99-11f1-b35e-4f505a8c89c4
verificar(() -> [lambda_para_gini(gini_gamma(1.0)), lambda_para_gini(gini_gamma(4.0))], [0.0, 0.5];
          rtol = 1e-6,
          pista = md"Empieza con `a, b = 0.0, 0.99`. En cada paso, `c = (a + b) / 2`; si `gini_gamma((1 + 2c) / (1 - c)) > G`, la solución está a la derecha (`a = c`); si no, `b = c`.")

# ╔═╡ 7db70c8a-ca99-11f1-b7ed-c178031138d9
verificar(() -> (r = mas_conectados([3, 9, 1, 7, 2], 2); r === missing ? missing : sort(r)), [2, 4];
          pista = md"`sortperm(grados; rev = true)[1:m]` o `partialsortperm(grados, 1:m; rev = true)`.")

# ╔═╡ 7db70c96-ca99-11f1-868c-4f192c9a4764
verificar(() -> fraccion_vulnerables([0, 1, 2, 5, 6, 10], 0.18), 0.5;
          pista = md"Cuenta los `k` con `k >= 1 && 1 / k >= φ` y divide entre `length(grados)`.")

# ╔═╡ Cell order:
# ╟─7db70b92-ca99-11f1-be11-598adb5d183f
# ╟─7db70be0-ca99-11f1-8548-4d2f73814b79
# ╟─7db70bea-ca99-11f1-9549-47938cb3db28
# ╠═7db70bf2-ca99-11f1-af27-43263d7baa23
# ╟─7db70bf2-ca99-11f1-a127-351330b93036
# ╠═7db70bfe-ca99-11f1-8b73-a77797fc6a7a
# ╟─7db70c08-ca99-11f1-af3c-db9bac6e2676
# ╠═7db70c08-ca99-11f1-a06e-cb71e41dc855
# ╟─7db70c12-ca99-11f1-8c41-735ad79d3288
# ╟─7db70c1c-ca99-11f1-8c11-ad0361929a12
# ╟─7db70c24-ca99-11f1-98dc-2d961441b910
# ╟─7db70c24-ca99-11f1-bef4-fdfcc6468358
# ╠═7db70c30-ca99-11f1-b6a8-f95402e32f24
# ╟─7db70c3a-ca99-11f1-8901-2bf3f09c3053
# ╟─7db70c3a-ca99-11f1-842c-23a21e9bda10
# ╟─7db70c44-ca99-11f1-b537-190a6584f26b
# ╠═7db70c44-ca99-11f1-9670-71c4849ac597
# ╟─7db70c4e-ca99-11f1-87d8-8940e485d01f
# ╠═7db70c4e-ca99-11f1-86bc-15e3c2e108d6
# ╟─7db70c58-ca99-11f1-b35e-4f505a8c89c4
# ╠═7db70c58-ca99-11f1-8f9a-55b7597be0ad
# ╟─7db70c62-ca99-11f1-8996-6d837d64f631
# ╟─7db70c62-ca99-11f1-b50d-69bf1edaeb6f
# ╟─7db70c6c-ca99-11f1-8040-6b2521b2ec04
# ╠═7db70c6c-ca99-11f1-b7a4-45e63be46bb9
# ╟─7db70c76-ca99-11f1-9ce4-5f64b0eff02d
# ╠═7db70c76-ca99-11f1-9f99-cb9524652177
# ╟─7db70c80-ca99-11f1-ad86-8bcbb50b0075
# ╠═7db70c80-ca99-11f1-bdc6-1b7821a8ebf8
# ╟─7db70c8a-ca99-11f1-b7ed-c178031138d9
# ╟─7db70c8a-ca99-11f1-ae73-f11b1e95c48f
# ╠═7db70c96-ca99-11f1-8dc9-7fed5f7ad10c
# ╟─7db70c96-ca99-11f1-868c-4f192c9a4764
# ╟─7db70c9e-ca99-11f1-a2c1-2b912b2da96a
# ╠═7db70c9e-ca99-11f1-8839-297befae210f
# ╟─7db70ca8-ca99-11f1-8f0a-3fb7367ebee3
# ╟─7db70ca8-ca99-11f1-9283-9b056f24e1ee
# ╟─7db70cb2-ca99-11f1-bb49-9987b621a2b5
# ╠═7db70cbc-ca99-11f1-9999-13e23f32a422
# ╟─7db70cbc-ca99-11f1-af47-c3145c896e45
# ╠═7db70cc8-ca99-11f1-8b21-87fe79a2cf49
# ╟─7db70cc8-ca99-11f1-820d-e3e7e5672402
# ╠═7db70cd0-ca99-11f1-9bc3-694f119f9aac
# ╟─7db70cd0-ca99-11f1-847a-0da4c3783373
# ╟─7db70cda-ca99-11f1-bb32-dd7a6e47e176
# ╠═7db70cda-ca99-11f1-8c6d-1f44cb6bfb1a
# ╠═7db70ce4-ca99-11f1-8e5a-5bc92af01953
# ╠═7db70cee-ca99-11f1-a5b4-73a597a989a9
# ╠═7db70cee-ca99-11f1-8cc6-2be3f0d3e03a
# ╠═7db70cfa-ca99-11f1-bbe7-15d3fa44abba
# ╠═7db70cfa-ca99-11f1-8b29-d3eacdf2e5a1
# ╠═7db70d02-ca99-11f1-b167-9597daac4068
# ╠═7db70d02-ca99-11f1-a66d-7b47cb5bf6f8
# ╠═7db70d0c-ca99-11f1-8d93-b98ed2ef8f95
# ╠═7db70d0c-ca99-11f1-aeec-7fb4d776be67
# ╠═7db70d0c-ca99-11f1-919a-3788476e4ff6
# ╠═7db70d16-ca99-11f1-a12a-d99f99c579aa
# ╠═7db70d16-ca99-11f1-a035-2fa2f243110b
# ╠═7db70d20-ca99-11f1-8455-b519d54c2b90
# ╠═7db70d20-ca99-11f1-8f23-d3d354127634
# ╠═7db70d28-ca99-11f1-8b27-912053711850
# ╠═7db70d34-ca99-11f1-8a2f-0f88a57f3c92
# ╠═7db70d3e-ca99-11f1-837e-55657c87fbdd
