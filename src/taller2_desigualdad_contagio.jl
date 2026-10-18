#%% hide
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

#%% hide
TableOfContents(title = "Contenido")

#%% hide
md"""
## 1. Parte 1: la desigualdad en tu entidad

Cada hogar de la ENIGH representa a muchos hogares del país: el **factor de
expansión** dice a cuántos. Todas las medidas de esta parte usan ese factor
como peso (sesión 8).
"""

#%%
enigh = CSV.read(joinpath(@__DIR__, "..", "Data", "enigh2024_hogares.csv"), DataFrame);

#%% hide
md"""
**Tu entidad:** $(@bind entidad Select(collect(1:32) .=> ENTIDADES, default = 1))
"""

#%%
begin
    hogares = filter(:entidad => ==(entidad), enigh)
    x_ent, w_ent = hogares.ingreso, Float64.(hogares.factor)
    x_nac, w_nac = enigh.ingreso, Float64.(enigh.factor)
end;

#%% hide
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

#%%
function participacion_top(x, w, q)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> [participacion_top([1.0, 2.0, 3.0, 4.0], [1.0, 1.0, 1.0, 1.0], 0.25),
                 participacion_top([1.0, 2.0, 3.0, 4.0], [1.0, 1.0, 1.0, 1.0], 0.5),
                 participacion_top([10.0, 1.0], [1.0, 9.0], 0.1)],
          [0.4, 0.7, 10 / 19];
          pista = md"Con `o = sortperm(x; rev = true)`, acumula `cumsum(w[o]) ./ sum(w)`; `n = findfirst(>=(q), acumulado)`; la respuesta es `sum(x[o][1:n] .* w[o][1:n]) / sum(x .* w)`.")

#%% hide
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

#%% hide
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

#%% hide
md"""
**Pregunta 1** (10 puntos). ¿Qué tan desigual es el ingreso en tu entidad
comparado con el país? Usa el Gini y la participación del 10% y del 1% más
ricos. ¿Dice lo mismo la curva de Lorenz que el Gini? ¿Por qué dos entidades
con el mismo Gini podrían tener curvas de Lorenz distintas?
"""

#%%
md"""
**Tu respuesta.** Escribe aquí.
"""

#%% hide
md"""
Ahora la **forma** de la distribución. Dividimos el ingreso entre su media
ponderada y dibujamos la CCDF en dos escalas: en **semilog** una exponencial
se ve recta, y en **log-log** se ve recta una ley de Pareto (sesiones 3 y 8).
"""

#%% hide
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

#%% hide
md"""
**Pregunta 2** (10 puntos). ¿En qué parte de la distribución se parece tu
entidad a una exponencial, y en qué parte se separa? ¿Hacia dónde se separa:
hay más o menos hogares muy ricos que los que predice la exponencial? ¿En
qué escala se ve mejor cada cosa y por qué?
"""

#%%
md"""
**Tu respuesta.** Escribe aquí.
"""

#%% hide
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

#%%
function lambda_para_gini(G)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> [lambda_para_gini(gini_gamma(1.0)), lambda_para_gini(gini_gamma(4.0))], [0.0, 0.5];
          rtol = 1e-6,
          pista = md"Empieza con `a, b = 0.0, 0.99`. En cada paso, `c = (a + b) / 2`; si `gini_gamma((1 + 2c) / (1 - c)) > G`, la solución está a la derecha (`a = c`); si no, `b = c`.")

#%%
modelo_ent = let
    λ = lambda_para_gini_ref(gini(x_ent, w_ent))
    (λ = λ, m = simular_ahorro(fill(λ, 5_000)))
end;

#%% hide
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

#%% hide
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

#%% hide
md"""
**Pregunta 3** (10 puntos). Calibramos el modelo para que tenga el mismo
Gini que tu entidad. ¿Reproduce también el cuerpo de la distribución (la
mediana y el percentil 10)? ¿Y la cola (el percentil 99 y el máximo)? Si
coincide el Gini pero no la cola, ¿qué te dice eso sobre resumir la
desigualdad con un solo número?
"""

#%%
md"""
**Tu respuesta.** Escribe aquí.
"""

#%% hide
md"""
**Pregunta 4** (10 puntos). ¿Qué significa en tu entidad el valor de
``\lambda`` que encontraste? Menciona dos razones por las que el modelo de
intercambio no puede ser la explicación completa de la distribución del
ingreso (piensa en qué mide la ENIGH, en el trabajo y en la herencia), y una
razón por la que de todos modos es útil (sesiones 8 y 9).
"""

#%%
md"""
**Tu respuesta.** Escribe aquí.
"""

#%% hide
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

#%%
function mas_conectados(grados, m)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> (r = mas_conectados([3, 9, 1, 7, 2], 2); r === missing ? missing : sort(r)), [2, 4];
          pista = md"`sortperm(grados; rev = true)[1:m]` o `partialsortperm(grados, 1:m; rev = true)`.")

#%% hide
md"""
### Ejercicio 4 (5 puntos)

Escribe `fraccion_vulnerables(grados, φ)`: la fracción de nodos que caen con
**un solo** vecino caído, es decir, con grado ``k \ge 1`` y
``1/k \ge \varphi``.
"""

#%%
function fraccion_vulnerables(grados, φ)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> fraccion_vulnerables([0, 1, 2, 5, 6, 10], 0.18), 0.5;
          pista = md"Cuenta los `k` con `k >= 1 && 1 / k >= φ` y divide entre `length(grados)`.")

#%% hide
md"""
**El experimento.** Para cada configuración simulamos 100 redes de 2,000
nodos con umbral ``\varphi = 0.18``. En cada red comparamos:

- **¿Quién cae primero?** Un nodo al azar contra el nodo de mayor grado.
- **¿A quién se protege?** Una fracción de los nodos, elegidos al azar o entre
  los más conectados. Después cae un nodo al azar (no protegido).

Tipo de red: $(@bind red_t2 Select(["Erdős-Rényi", "Barabási-Albert"]))

Grado promedio ``z``: $(@bind z_t2 Select([2, 4, 6, 8], default = 4))
"""

#%%
experimento_t2 = experimento_proteccion(red_t2, z_t2);

#%% hide
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

#%% hide
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

#%% hide
md"""
**Pregunta 5** (10 puntos). ¿Importa quién cae primero? Compara un nodo al
azar con el de mayor grado en las dos redes con ``z = 4``. ¿En cuál de las
dos redes la diferencia es más grande y por qué? Explica qué les pasa a los
vecinos del nodo de mayor grado cuando éste cae.
"""

#%%
md"""
**Tu respuesta.** Escribe aquí.
"""

#%% hide
md"""
**Pregunta 6** (10 puntos). En una epidemia, vacunar a los más conectados es
mucho más eficaz que vacunar al azar. ¿Pasa lo mismo aquí? Compara las dos
estrategias de protección en Erdős-Rényi con ``z = 4`` y con ``z = 2``, y
después en Barabási-Albert. Explica en qué casos conviene proteger a los más
conectados usando la idea de nodo **vulnerable** (ejercicio 4): ¿son
vulnerables los nodos más conectados? ¿Qué papel tienen en la cascada si nadie
los protege, cuando la red está muy conectada y cuando apenas lo está?
"""

#%%
md"""
**Tu respuesta.** Escribe aquí.
"""

#%% hide
md"""
**Pregunta 7** (10 puntos). Escribe una recomendación de cinco a ocho líneas
para un regulador bancario con base en las preguntas 5 y 6. ¿Dónde pondrías
la vigilancia y dónde los requisitos de capital? Menciona una limitación del
modelo que podría cambiar tu recomendación (por ejemplo, que en la red real
los bancos grandes también son los que más prestan, sesión 12).
"""

#%%
md"""
**Tu respuesta.** Escribe aquí.
"""

#%% hide
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

#%% hide
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para resolver el taller.
"""

#%%
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

#%%
begin
    const AZUL    = "#2a78d6"   # serie principal
    const NARANJA = "#eb6834"   # segunda serie
    const GRIS    = "#52514e"   # referencias
end;

#%%
"Nombres de las entidades federativas, en el orden de su clave del INEGI."
const ENTIDADES = ["Aguascalientes", "Baja California", "Baja California Sur", "Campeche",
                   "Coahuila", "Colima", "Chiapas", "Chihuahua", "Ciudad de México", "Durango",
                   "Guanajuato", "Guerrero", "Hidalgo", "Jalisco", "Estado de México",
                   "Michoacán", "Morelos", "Nayarit", "Nuevo León", "Oaxaca", "Puebla",
                   "Querétaro", "Quintana Roo", "San Luis Potosí", "Sinaloa", "Sonora",
                   "Tabasco", "Tamaulipas", "Tlaxcala", "Veracruz", "Yucatán", "Zacatecas"];

#%%
miles(x) = replace(string(round(Int, x)), r"(?<=\d)(?=(\d{3})+$)" => ",")

#%%
"Curva de Lorenz ponderada: fracción acumulada de la población y del ingreso."
function lorenz(x, w)
    o = sortperm(x)
    W, S = cumsum(w[o]), cumsum(x[o] .* w[o])
    [0; W ./ W[end]], [0; S ./ S[end]]
end

#%%
"Índice de Gini con pesos (por defecto, todos iguales): 1 menos el doble del área bajo la curva de Lorenz."
function gini(x, w = ones(length(x)))
    p, L = lorenz(x, w)
    1 - sum(diff(p) .* (L[2:end] .+ L[1:end-1]))
end

#%%
"Participación en el ingreso de la fracción `q` más rica (versión de referencia)."
function participacion_top_ref(x, w, q)
    o = sortperm(x; rev = true)
    n = findfirst(>=(q), cumsum(w[o]) ./ sum(w))
    sum(x[o][1:n] .* w[o][1:n]) / sum(x .* w)
end

#%%
"Gini de una gamma con forma n."
gini_gamma(n) = exp(Distributions.SpecialFunctions.loggamma(n + 0.5) -
                    Distributions.SpecialFunctions.loggamma(n)) / (n * sqrt(π))

#%%
"Ahorro λ con el que el modelo de intercambio tiene Gini `G` (versión de referencia)."
function lambda_para_gini_ref(G)
    a, b = 0.0, 0.99
    for _ in 1:60
        c = (a + b) / 2
        gini_gamma((1 + 2c) / (1 - c)) > G ? (a = c) : (b = c)
    end
    (a + b) / 2
end

#%%
"Datos ordenados de mayor a menor y su CCDF empírica, P(X ≥ x)."
function ccdf_empirica(x)
    s = sort(x; rev = true)
    s, (1:length(s)) ./ length(s)
end

#%%
"CCDF con pesos: P(X ≥ x) usando el factor de expansión."
function ccdf_ponderada(x, w)
    o = sortperm(x; rev = true)
    x[o], cumsum(w[o]) ./ sum(w)
end

#%%
"Estimador de Hill del exponente de la cola, con la fracción `q` de datos más grandes."
function hill_cola(x, q)
    s = sort(x; rev = true)
    n = round(Int, q * length(s))
    n / sum(log.(s[1:n] ./ s[n]))
end

#%%
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

#%%
"Ahorros heterogéneos λᵢ = 1 - Uᵢ^(1/a), con Uᵢ uniforme."
ahorros_heterogeneos(a; N = 5_000, rng = Xoshiro(3)) = 1 .- rand(rng, N) .^ (1 / a)

#%%
"Fracción de nodos vulnerables (versión de referencia)."
fraccion_vulnerables_ref(grados, φ) = count(k -> k >= 1 && 1 / k >= φ, grados) / length(grados)

#%%
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

#%%
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

#%%
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
