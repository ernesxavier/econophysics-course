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

# ╔═╡ b6167956-ca96-11f1-a3bd-4fcefb85f14e
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

# ╔═╡ b6167848-ca96-11f1-981a-cb3852251aa1
md"""
# Taller del Módulo III: cierre del ciclo

**Qué se evalúa.** El curso empezó midiendo hechos estilizados en datos
reales (sesiones 1 y 2) y terminó construyendo mercados artificiales con
agentes (sesión 20). En este taller cierras el ciclo: aplicas el mismo
diagnóstico a series simuladas, averiguas **qué mecanismo** del modelo produce
cada hecho y tratas de calibrarlo. Se evalúa la **interpretación**; no se
piden derivaciones.

**Tu serie.** El profesor te asigna un activo: uno de los cuatro de las
sesiones 1 y 2 (`Data/precios.csv`) o una de las acciones del S&P 100
(`Data/sp100_precios.csv`). Elígelo en el selector de la sección 1 y no lo
cambies.

**Cómo entregar.**
1. Completa las funciones de los ejercicios de código (reemplaza `missing`).
   La caja verde de `verificar` te dice si están bien.
2. Responde cada pregunta de interpretación en la celda **Tu respuesta**
   que está debajo, en tres a seis líneas. Cita los números de tus tablas.
3. Entrega el archivo `.jl` del cuaderno y una exportación a HTML o PDF (menú
   de Pluto, "Export"). El trabajo es individual: puedes discutir ideas con tus
   compañeros, pero el código y la redacción deben ser tuyos.

**Cómo se califica** (100 puntos, 10% de la calificación final):

| Parte | Código | Interpretación | Total |
|---|---|---|---|
| 1. Tres modelos frente a tus datos | 10 | 20 | 30 |
| 2. ¿Qué mecanismo produce cada hecho? | 10 | 30 | 40 |
| 3. Calibrar a mano | 10 | 20 | 30 |

Una respuesta completa dice **qué** se observa, **con qué número o
gráfica** y **por qué** importa. La ruta avanzada, al final, es opcional y
cuenta como extra.

**En el semestre** hay una sesión de taller en clase: en los primeros 20
minutos se resuelven dudas del módulo, después se trabaja en parejas en las
partes 1 y 2, y al final cada pareja comenta una respuesta con el grupo. En el
trimestre el taller se resuelve en casa.

**Repaso recomendado** (en `Class_notes/repaso/`): E.B3 (promediar
simulaciones), S.B5 (pruebas por simulación), C.B4 (la función logística),
F.B4 (transiciones de fase). Para la ruta avanzada: S.I4 (datos dependientes)
y S.A3 (calibración por momentos simulados).
"""

# ╔═╡ b616787a-ca96-11f1-9997-43beb8c952c9
TableOfContents(title = "Contenido")

# ╔═╡ b616787a-ca96-11f1-b90b-21115ee8cb53
md"""
## 1. Parte 1: tres modelos frente a tus datos

Comparamos tu serie con tres modelos, cada uno simulado con el mismo número
de días y con **10 semillas** distintas:

- **Caminata gaussiana:** retornos normales e independientes (sesión 5).
- **GARCH(1,1)** ajustado a tus datos por máxima verosimilitud (sesión 5).
- **Franke-Westerhoff** con los parámetros que los autores estimaron para el
  S&P 500 (sesión 20): fundamentalistas y chartistas que cambian de estrategia
  por contagio.
"""

# ╔═╡ b6167884-ca96-11f1-baf8-e97912478392
begin
    indices = CSV.read(joinpath(@__DIR__, "..", "Data", "precios.csv"), DataFrame)
    acciones = CSV.read(joinpath(@__DIR__, "..", "Data", "sp100_precios.csv"), DataFrame)
    sectores = CSV.read(joinpath(@__DIR__, "..", "Data", "sp100_sectores.csv"), DataFrame)
    opciones = [[a => a for a in unique(indices.activo)];
                [f.ticker => "$(f.ticker): $(f.nombre)" for f in eachrow(sectores)]]
end;

# ╔═╡ b616788e-ca96-11f1-8797-5fb338c8b969
md"""
**Tu serie:** $(@bind activo Select(opciones, default = "S&P 500"))
"""

# ╔═╡ b6167898-ca96-11f1-b4f6-715dc22a62cb
r = let
    d = activo in indices.activo ? filter(:activo => ==(activo), indices) :
                                   rename(filter(:ticker => ==(activo), acciones), :ticker => :activo)
    diff(log.(sort(d, :fecha).cierre))
end;

# ╔═╡ b6167898-ca96-11f1-b70f-13738b505be3
md"""
### Ejercicio 1 (10 puntos)

Escribe `resumen_hechos(r)`. Estandariza los retornos,
``z_t = (r_t - \bar r)/s``, y devuelve un vector con cuatro números:

1. la curtosis en exceso de ``z`` (`kurtosis`);
2. la autocorrelación de ``z`` con rezago 1;
3. la autocorrelación de ``|z|`` con rezago 5;
4. la fracción de días con ``|z| > 3``.

Para las autocorrelaciones usa `autocor(x, [τ])[1]`.
"""

# ╔═╡ b61678a4-ca96-11f1-bac3-59ed3ec9ce19
function resumen_hechos(r)
    missing  # reemplaza esta línea
end

# ╔═╡ b61678c0-ca96-11f1-8fc6-a56c14ff1e4a
md"""
**Pregunta 1** (10 puntos). ¿Qué hechos estilizados reproduce cada uno de los
tres modelos? ¿Cuál se acerca más a tu serie en la curtosis y cuál en la
memoria de la volatilidad? ¿Cae el valor de tu serie dentro del rango de las
10 semillas de algún modelo?
"""

# ╔═╡ b61678c0-ca96-11f1-91a1-afaa112842f0
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ b61678ca-ca96-11f1-aed0-49778a7cfeda
md"""
**Pregunta 2** (10 puntos). El GARCH y el modelo de Franke-Westerhoff
producen agrupamiento de volatilidad, pero de manera muy distinta. ¿Qué
**supone** el GARCH para obtenerlo y qué **explica** el modelo de agentes?
Si fueras un regulador que quiere saber qué pasaría con una nueva regla del
mercado, ¿cuál de los dos te serviría más y por qué?
"""

# ╔═╡ b61678ca-ca96-11f1-a04b-852589f67ecc
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ b61678d2-ca96-11f1-a604-85db8fb1b157
md"""
## 2. Parte 2: ¿qué mecanismo produce cada hecho?

El modelo de Franke-Westerhoff tiene varios ingredientes:

- **Extrapolación** ``\chi``: los chartistas compran si el precio subió.
- **Ruido distinto:** la demanda de los chartistas es mucho más ruidosa que la
  de los fundamentalistas (``\sigma_c = 2.087`` contra ``\sigma_f = 0.758``).
- **Cambio de estrategia** con una probabilidad logit cuyo atractivo tiene
  tres partes: el **contagio** ``\alpha_n`` (conviene hacer lo que hace la
  mayoría), la **predisposición** ``\alpha_0`` y el **temor al
  desalineamiento** ``\alpha_p`` (si el precio se aleja mucho del valor
  fundamental, conviene ser fundamentalista).

Para saber qué ingrediente produce cada hecho, **apagamos uno a la vez** y
repetimos el diagnóstico con 10 semillas.
"""

# ╔═╡ b61678e8-ca96-11f1-a9e6-cf5e05ccd7f1
md"""
### Ejercicio 2 (10 puntos)

Escribe `volatilidad_por_regimen(r, z; umbral = 0.5)`. Aquí `z` es la
fracción de chartistas en cada día. Para los días ``t = 2, \dots, n``, separa
los retornos según si la fracción de chartistas **del día anterior**,
``z_{t-1}``, es mayor que `umbral` o no, y devuelve `[a, b]`: el promedio de
``|r_t|`` en el primer grupo (dominan los chartistas) y en el segundo.
"""

# ╔═╡ b61678e8-ca96-11f1-99de-4b8b8ea42f79
function volatilidad_por_regimen(r, z; umbral = 0.5)
    missing  # reemplaza esta línea
end

# ╔═╡ b61678fc-ca96-11f1-82db-0bb1ba9d2d85
md"""
**Pregunta 3** (10 puntos). Según la tabla de la ablación, ¿qué ingredientes
son **indispensables** para que aparezcan las colas pesadas y el agrupamiento
de volatilidad, y cuáles casi no importan? ¿Te sorprende alguno? (Piensa en
la extrapolación, que es lo que define a un chartista.)
"""

# ╔═╡ b6167904-ca96-11f1-82ee-557f4ac59214
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ b6167904-ca96-11f1-8041-dd0b74f74f92
md"""
**Pregunta 4** (10 puntos). Con la gráfica y los números del ejercicio 2,
explica en tus palabras **cómo** el modelo produce agrupamiento de
volatilidad. ¿Por qué hace falta que los dos grupos tengan ruidos distintos?
¿Qué papel juega el contagio en que los periodos agitados **duren**?
"""

# ╔═╡ b6167910-ca96-11f1-bc8b-23e368db32cf
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ b6167910-ca96-11f1-a1cf-6b7482c34fb4
md"""
**Pregunta 5** (10 puntos). Que un mecanismo **pueda** producir lo observado
no prueba que **sea** el mecanismo: el GARCH, por ejemplo, también produce
agrupamiento. ¿Qué evidencia adicional, fuera de los precios, buscarías para
saber si en los mercados reales hay cambios entre fundamentalistas y
chartistas? Da al menos dos ejemplos concretos.
"""

# ╔═╡ b6167910-ca96-11f1-b4e8-f5c23bd64a81
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ b616791a-ca96-11f1-8685-71dfef940fcb
md"""
## 3. Parte 3: calibrar a mano

Mueve el **contagio** ``\alpha_n`` y el **ruido de los chartistas**
``\sigma_c`` para que el modelo se parezca a tu serie en la curtosis y en la
memoria de la volatilidad. Cada punto usa 10 semillas.

- Contagio ``\alpha_n``: $(@bind αn_t3 Slider(0.0:0.1:3.2, default = 1.8, show_value = true))
- Ruido de los chartistas ``\sigma_c``: $(@bind σc_t3 Slider(1.0:0.25:5.0, default = 2.0, show_value = true))

### Ejercicio 3 (10 puntos)

Escribe `distancia_momentos(m, d)`: la suma de los cuadrados de las
diferencias **relativas** entre los momentos del modelo `m` y los de los
datos `d` (dos vectores del mismo tamaño):
``\sum_i \big((m_i - d_i)/d_i\big)^2``. Es la medida que usaremos para
decir qué tan cerca estás.
"""

# ╔═╡ b616791a-ca96-11f1-b570-d195f6460d2c
function distancia_momentos(m, d)
    missing  # reemplaza esta línea
end

# ╔═╡ b6167936-ca96-11f1-847e-2d234c371af4
md"""
**Pregunta 6** (10 puntos). Reporta los mejores valores de ``\alpha_n`` y
``\sigma_c`` que encontraste y su distancia. ¿Qué momento fue más fácil de
igualar? Con esos parámetros, ¿qué tanto varía la curtosis entre semillas?
¿Por qué compararías la media de 10 semillas con tus datos y no una sola
simulación? ¿Te bastarían 10 si quisieras distinguir dos parámetros
cercanos?
"""

# ╔═╡ b6167942-ca96-11f1-bf97-2937b97c1039
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ b6167942-ca96-11f1-b7ca-2b71aff44c9d
md"""
**Pregunta 7** (10 puntos). Sube el contagio a ``\alpha_n = 3`` o más. ¿Qué
les pasa a la curtosis y al agrupamiento? Mira la gráfica de la fracción de
chartistas y explica por qué. Relaciónalo con el campo medio de las sesiones
13 y 14: ¿qué pasa con un sistema de espines cuando la interacción es mucho
mayor que el ruido?
"""

# ╔═╡ b6167942-ca96-11f1-8ccf-73a0d5374769
md"""
**Tu respuesta.** Escribe aquí.
"""

# ╔═╡ b616794c-ca96-11f1-ac94-f1444fc53a27
md"""
## 4. Ruta avanzada (extra, hasta 10 puntos)

Haz **dos** de los cuatro. Agrega las celdas que necesites al final del
cuaderno.

- **A1. Momentos simulados.** Barre una cuadrícula de ``(\alpha_n, \sigma_c)``
  y elige el par que minimiza la distancia, ahora con tres momentos: agrega la
  ACF de ``|r|`` con rezago 25. ¿Se parece a tu calibración a mano? (Repaso
  S.A3; Franke y Westerhoff, 2012.)
- **A2. Brock-Hommes con ruido.** Agrega a la ecuación de precios de
  Brock-Hommes (sesión 20, `simular_brock_hommes`) un choque normal pequeño
  cada periodo y aplica el diagnóstico de la parte 1 a ``x_t - x_{t-1}``.
  ¿Qué hechos reproduce?
- **A3. Ljung-Box.** Aplica la prueba de Ljung-Box (repaso S.I2) con 20
  rezagos a ``r`` y a ``|r|`` de tu serie y de cada modelo. ¿Distingue la
  prueba a los modelos entre sí?
- **A4. Ventanas de tiempo.** Divide tu serie en dos mitades y repite la
  calibración en cada una. ¿Cambian los parámetros? ¿Qué te dice eso sobre la
  idea de que los parámetros del modelo son constantes?
"""

# ╔═╡ b616794c-ca96-11f1-802e-a3300e661634
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para resolver el taller.
"""

# ╔═╡ b6167956-ca96-11f1-a936-8df2b346e4b1
begin
    const AZUL    = "#2a78d6"   # serie principal
    const NARANJA = "#eb6834"   # segunda serie
    const GRIS    = "#52514e"   # referencias
end;

# ╔═╡ b6167960-ca96-11f1-9acf-31dd9c82d894
"Curtosis, ACF de r (rezago 1), ACF de |r| (rezago 5) y fracción de días con |z| > 3."
function resumen_hechos_ref(r)
    z = (r .- mean(r)) ./ std(r)
    [kurtosis(z), autocor(z, [1])[1], autocor(abs.(z), [5])[1], mean(abs.(z) .> 3)]
end

# ╔═╡ b6167960-ca96-11f1-ae5a-81b4a4ea72e5
"Retorno absoluto medio cuando dominaban los chartistas el día anterior y cuando no."
function volatilidad_por_regimen_ref(r, z; umbral = 0.5)
    alto = z[1:end-1] .> umbral
    x = abs.(r[2:end])
    [mean(x[alto]), mean(x[.!alto])]
end

# ╔═╡ b616796a-ca96-11f1-aa5c-bbedc915d2c8
distancia_momentos_ref(m, d) = sum(((m .- d) ./ d) .^ 2)

# ╔═╡ b616796a-ca96-11f1-9418-3713fca7660b
"""
Simula un GARCH(1,1) con varianza de largo plazo `s2`, es decir,
ω = s2 (1 − a − b). Descarta 500 pasos iniciales.
"""
function simular_garch(n, a, b, s2; rng = Xoshiro(10))
    ω = s2 * (1 - a - b)
    σ2 = s2
    for _ in 1:500
        σ2 = ω + (a * randn(rng)^2 + b) * σ2
    end
    r, v = zeros(n), zeros(n)
    for t in 1:n
        v[t] = σ2
        r[t] = sqrt(σ2) * randn(rng)
        σ2 = ω + a * r[t]^2 + b * σ2
    end
    (r = r, σ2 = v)
end

# ╔═╡ b6167974-ca96-11f1-8890-95878475aae7
"Log-verosimilitud normal del GARCH(1,1) con varianza de largo plazo s2."
function garch_loglik(x, a, b, s2)
    ω = s2 * (1 - a - b)
    σ2 = s2
    ℓ = 0.0
    for t in eachindex(x)
        ℓ -= 0.5 * (log(σ2) + x[t]^2 / σ2)
        σ2 = ω + a * x[t]^2 + b * σ2
    end
    ℓ
end

# ╔═╡ b6167974-ca96-11f1-a1ea-5719f3d778d1
"""
Ajuste del GARCH(1,1) por máxima verosimilitud en una cuadrícula de (a, b),
con la varianza de largo plazo fijada en la varianza muestral.
"""
function ajustar_garch(r; as = 0.005:0.005:0.4, bs = 0.5:0.0025:0.9975)
    x = r .- mean(r)
    s2 = var(x)
    mejor = (ℓ = -Inf, a = NaN, b = NaN, s2 = s2)
    for a in as, b in bs
        a + b < 0.9995 || continue
        ℓ = garch_loglik(x, a, b, s2)
        ℓ > mejor.ℓ && (mejor = (ℓ = ℓ, a = a, b = b, s2 = s2))
    end
    mejor
end

# ╔═╡ b616797e-ca96-11f1-912d-237c96559d9e
"""
Simula el modelo DCA-HPM de Franke y Westerhoff (2012): fundamentalistas y
chartistas con contagio (herding), predisposición y temor al desalineamiento.
Los valores por omisión son los que los autores estimaron para el S&P 500.
Devuelve el log-precio `p`, los retornos `r` y la fracción de chartistas `z`.
"""
function simular_franke_westerhoff(; pasos = 2500, calent = 1000, beta = 1.0, mu = 0.01,
                                   phi = 0.12, chi = 1.5, sigma_f = 0.758, sigma_c = 2.087,
                                   alpha_0 = -0.327, alpha_n = 1.79, alpha_p = 18.43,
                                   p_fund = 0.0, rng = Xoshiro(42))
    T = calent + pasos + 1
    p = fill(p_fund, T)
    nf = fill(0.5, T)
    for t in 3:T
        d_f = phi * (p_fund - p[t - 1]) + sigma_f * randn(rng)
        d_c = chi * (p[t - 1] - p[t - 2]) + sigma_c * randn(rng)
        p[t] = p[t - 1] + mu * (nf[t - 1] * d_f + (1 - nf[t - 1]) * d_c)
        # Atractivo relativo del fundamentalismo en t - 1: contagio + predisposición + desalineamiento
        a = alpha_n * (2 * nf[t - 1] - 1) + alpha_0 + alpha_p * (p[t - 1] - p_fund)^2
        nf[t] = 1 / (1 + exp(-beta * a))
    end
    idx = (calent + 2):T
    (p = p[idx], r = p[idx] .- p[idx .- 1], z = 1 .- nf[idx])
end

# ╔═╡ b61678ac-ca96-11f1-ad54-6b2f544cf8c2
modelos_t3 = let
    n = length(r)
    aj = ajustar_garch(r)
    semillas = 1:10
    (garch = aj,
     gauss = [randn(Xoshiro(s), n) for s in semillas],
     garch_r = [simular_garch(n, aj.a, aj.b, aj.s2; rng = Xoshiro(s)).r for s in semillas],
     fw = [simular_franke_westerhoff(; pasos = n, rng = Xoshiro(s)).r for s in semillas])
end;

# ╔═╡ b61678b6-ca96-11f1-944a-af98ed8336d6
let
    m = modelos_t3
    lags = 1:40
    zr = (r .- mean(r)) ./ std(r)
    xs = range(0.2, 8, length = 150)
    p1 = plot(xs, 2ccdf.(Normal(), xs); color = GRIS, ls = :dash, lw = 1.5, yscale = :log10,
              ylims = (1e-4, 1), xlims = (0, 8), label = "Normal", legend = :topright,
              xlabel = "x (desviaciones estándar)", ylabel = "P(|z| > x)", title = "Colas")
    lineas = [("Tus datos", r, :black, :solid), ("GARCH", m.garch_r[1], NARANJA, :dot),
              ("Franke-Westerhoff", m.fw[1], AZUL, :solid)]
    for (nombre, x, c, ls) in lineas
        z = sort(abs.((x .- mean(x)) ./ std(x)); rev = true)
        plot!(p1, z, (1:length(z)) ./ length(z); color = c, ls = ls, lw = 2, label = nombre)
    end
    p2 = plot(; xlabel = "rezago τ (días)", ylabel = "ACF de |r|", ylims = (-0.05, 0.45),
              legend = :topright, title = "Memoria de la volatilidad")
    for (nombre, x, c, ls) in [lineas; [("Caminata gaussiana", m.gauss[1], GRIS, :dash)]]
        plot!(p2, lags, autocor(abs.(x), lags); color = c, ls = ls, lw = 2, label = nombre)
    end
    plot(p1, p2; layout = (1, 2), size = (780, 340))
end

# ╔═╡ b61678b6-ca96-11f1-a062-67f1b5ac19ab
let
    m = modelos_t3
    f(v) = round(v, digits = 3)
    celda(sims, i) = (v = [resumen_hechos_ref(x)[i] for x in sims];
                      "$(f(mean(v))) ($(f(minimum(v))) a $(f(maximum(v))))")
    d = resumen_hechos_ref(r)
    fila(nombre, sims) = "| $nombre | " * join([celda(sims, i) for i in 1:4], " | ") * " |"
    Markdown.parse("""
    Simulaciones: media de 10 semillas (mínimo a máximo). En una normal, el
    0.27% de los días tiene ``|z| > 3``.

    | Serie | Curtosis en exceso | ACF de r, rezago 1 | ACF de ∣r∣, rezago 5 | Días con ∣z∣ > 3 |
    |---|---:|---:|---:|---:|
    | **$activo** | $(join(f.(d), " | ")) |
    $(fila("Caminata gaussiana", m.gauss))
    $(fila("GARCH(1,1), a = $(m.garch.a), b = $(m.garch.b)", m.garch_r))
    $(fila("Franke-Westerhoff", m.fw))
    """)
end

# ╔═╡ b61678d2-ca96-11f1-b47f-433e0c474e22
ablacion = let
    n = length(r)
    casos = [("Modelo completo", (;)),
             ("Sin extrapolación (χ = 0)", (chi = 0.0,)),
             ("Mismo ruido en los dos grupos", (sigma_c = 0.758,)),
             ("Sin contagio (α_n = 0)", (alpha_n = 0.0,)),
             ("Sin predisposición (α_0 = 0)", (alpha_0 = 0.0,)),
             ("Sin temor al desalineamiento (α_p = 0)", (alpha_p = 0.0,))]
    [(nombre = nombre,
      res = [resumen_hechos_ref(simular_franke_westerhoff(; pasos = n, rng = Xoshiro(s), kw...).r) for s in 1:10])
     for (nombre, kw) in casos]
end;

# ╔═╡ b61678de-ca96-11f1-a7ae-31857df7dd1d
let
    f(v) = round(v, digits = 3)
    celda(res, i) = (v = getindex.(res, i); "$(f(mean(v))) ($(f(minimum(v))) a $(f(maximum(v))))")
    filas = ["| $(a.nombre) | " * join([celda(a.res, i) for i in 1:4], " | ") * " |" for a in ablacion]
    Markdown.parse("""
    | Variante | Curtosis en exceso | ACF de r, rezago 1 | ACF de ∣r∣, rezago 5 | Días con ∣z∣ > 3 |
    |---|---:|---:|---:|---:|
    $(join(filas, "\n"))
    """)
end

# ╔═╡ b61678f2-ca96-11f1-9ff0-532529c0616e
regimenes = let
    sim = simular_franke_westerhoff(; pasos = length(r), rng = Xoshiro(1))
    (sim = sim, v = volatilidad_por_regimen_ref(sim.r, sim.z))
end;

# ╔═╡ b61678f2-ca96-11f1-b93f-e7db22f2f35d
let
    sim = regimenes.sim
    t = 1:min(1500, length(sim.r))
    p1 = plot(t, sim.r[t]; color = AZUL, lw = 0.6, ylabel = "retorno", title = "Franke-Westerhoff, semilla 1")
    p2 = plot(t, sim.z[t]; color = NARANJA, lw = 1.2, ylabel = "fracción de chartistas",
              xlabel = "día", ylims = (0, 1))
    hline!(p2, [0.5]; color = GRIS, ls = :dash)
    plot(p1, p2; layout = (2, 1), size = (720, 420))
end

# ╔═╡ b61678fc-ca96-11f1-8cf3-15ad6ab68aac
let
    sim, v = regimenes.sim, regimenes.v
    Markdown.parse("""
    En esta simulación, los chartistas dominan el $(round(100mean(sim.z .> 0.5), digits = 1))%
    de los días. El retorno absoluto medio es $(round(v[1], digits = 3)) cuando
    dominaban el día anterior y $(round(v[2], digits = 3)) cuando no:
    **$(round(v[1] / v[2], digits = 1)) veces** más.
    """)
end

# ╔═╡ b6167924-ca96-11f1-a53f-f575f2928aa6
calibracion = let
    n = length(r)
    sims = [simular_franke_westerhoff(; pasos = n, alpha_n = αn_t3, sigma_c = σc_t3, rng = Xoshiro(s))
            for s in 1:10]
    res = [resumen_hechos_ref(s.r) for s in sims]
    momentos(v) = [v[1], v[3]]                     # curtosis y ACF de |r| con rezago 5
    dm = momentos(resumen_hechos_ref(r))
    mm = mean(momentos.(res))
    (res = res, datos = dm, modelo = mm, distancia = distancia_momentos_ref(mm, dm),
     z = sims[1].z, chartistas = mean(mean(s.z) for s in sims))
end;

# ╔═╡ b616792e-ca96-11f1-aebb-1d9b6dd5a56c
let
    c = calibracion
    f(v) = round(v, digits = 3)
    rango(i) = (v = getindex.(c.res, i); "$(f(minimum(v))) a $(f(maximum(v)))")
    Markdown.parse("""
    | | Curtosis en exceso | ACF de ∣r∣, rezago 5 |
    |---|---:|---:|
    | **$activo** | $(f(c.datos[1])) | $(f(c.datos[2])) |
    | Modelo, media de 10 semillas | $(f(c.modelo[1])) | $(f(c.modelo[2])) |
    | Modelo, rango de las semillas | $(rango(1)) | $(rango(3)) |

    **Distancia:** $(round(c.distancia, digits = 3)). Fracción media de
    chartistas: $(round(c.chartistas, digits = 2)).
    """)
end

# ╔═╡ b616792e-ca96-11f1-a9cb-e5eee0348df6
let
    c = calibracion
    t = 1:min(1500, length(c.z))
    p = plot(t, c.z[t]; color = NARANJA, lw = 1, ylims = (0, 1), xlabel = "día",
             ylabel = "fracción de chartistas", size = (720, 220),
             title = "α_n = $αn_t3, σ_c = $σc_t3 (semilla 1)")
    hline!(p, [0.5]; color = GRIS, ls = :dash)
end

# ╔═╡ b6167988-ca96-11f1-8a8a-cdd4f9ed6fd3
"""
Modelo de Brock y Hommes (1998) con fundamentalistas y seguidores de
tendencia (sesión 20). Devuelve la desviación `x` y la fracción `n2` de
seguidores de tendencia.
"""
function simular_brock_hommes(beta, g; R = 1.10, C = 1.0, pasos = 200, calent = 500, x0 = 0.25)
    T = calent + pasos
    x = zeros(Float64, T + 3)
    n2 = fill(0.5, T + 3)
    x[1:3] .= (x0, x0 * 1.02, x0 * 0.98)

    for t in 4:(T + 3)
        ret_real = x[t - 1] - R * x[t - 2]
        pi1 = -R * x[t - 2] * ret_real - C
        pi2 = (g * x[t - 3] - R * x[t - 2]) * ret_real
        diff_pi = clamp(beta * (pi1 - pi2), -40.0, 40.0)
        n2[t] = 1.0 / (1.0 + exp(diff_pi))
        x[t] = (n2[t] * g * x[t - 1]) / R
    end

    idx = (calent + 4):(T + 3)
    (x = x[idx], n2 = n2[idx])
end

# ╔═╡ b6167988-ca96-11f1-a5f3-fb77ec69b2d0
const RETORNOS_PRUEBA = 0.01 .* randn(Xoshiro(2026), 500);

# ╔═╡ b6167988-ca96-11f1-b721-dde064ea24e4
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

# ╔═╡ b61678ac-ca96-11f1-a663-2f2e219d910b
verificar(() -> resumen_hechos(RETORNOS_PRUEBA), resumen_hechos_ref(RETORNOS_PRUEBA);
          pista = md"Con `z = (r .- mean(r)) ./ std(r)`, devuelve `[kurtosis(z), autocor(z, [1])[1], autocor(abs.(z), [5])[1], mean(abs.(z) .> 3)]`.")

# ╔═╡ b61678f2-ca96-11f1-85d8-9765c04f0251
verificar(() -> volatilidad_por_regimen([0.1, -0.2, 0.05, 0.3, -0.1], [0.8, 0.2, 0.9, 0.1, 0.5]),
          [0.25, 0.075];
          pista = md"Con `alto = z[1:end-1] .> umbral` y `x = abs.(r[2:end])`, devuelve `[mean(x[alto]), mean(x[.!alto])]`.")

# ╔═╡ b6167924-ca96-11f1-bf2d-db786205e345
verificar(() -> distancia_momentos([2.0, 0.1], [4.0, 0.2]), 0.5;
          pista = md"Usa `sum(((m .- d) ./ d) .^ 2)`.")

# ╔═╡ Cell order:
# ╟─b6167848-ca96-11f1-981a-cb3852251aa1
# ╟─b616787a-ca96-11f1-9997-43beb8c952c9
# ╟─b616787a-ca96-11f1-b90b-21115ee8cb53
# ╠═b6167884-ca96-11f1-baf8-e97912478392
# ╟─b616788e-ca96-11f1-8797-5fb338c8b969
# ╠═b6167898-ca96-11f1-b4f6-715dc22a62cb
# ╟─b6167898-ca96-11f1-b70f-13738b505be3
# ╠═b61678a4-ca96-11f1-bac3-59ed3ec9ce19
# ╟─b61678ac-ca96-11f1-a663-2f2e219d910b
# ╠═b61678ac-ca96-11f1-ad54-6b2f544cf8c2
# ╟─b61678b6-ca96-11f1-944a-af98ed8336d6
# ╟─b61678b6-ca96-11f1-a062-67f1b5ac19ab
# ╟─b61678c0-ca96-11f1-8fc6-a56c14ff1e4a
# ╠═b61678c0-ca96-11f1-91a1-afaa112842f0
# ╟─b61678ca-ca96-11f1-aed0-49778a7cfeda
# ╠═b61678ca-ca96-11f1-a04b-852589f67ecc
# ╟─b61678d2-ca96-11f1-a604-85db8fb1b157
# ╠═b61678d2-ca96-11f1-b47f-433e0c474e22
# ╟─b61678de-ca96-11f1-a7ae-31857df7dd1d
# ╟─b61678e8-ca96-11f1-a9e6-cf5e05ccd7f1
# ╠═b61678e8-ca96-11f1-99de-4b8b8ea42f79
# ╟─b61678f2-ca96-11f1-85d8-9765c04f0251
# ╠═b61678f2-ca96-11f1-9ff0-532529c0616e
# ╟─b61678f2-ca96-11f1-b93f-e7db22f2f35d
# ╟─b61678fc-ca96-11f1-8cf3-15ad6ab68aac
# ╟─b61678fc-ca96-11f1-82db-0bb1ba9d2d85
# ╠═b6167904-ca96-11f1-82ee-557f4ac59214
# ╟─b6167904-ca96-11f1-8041-dd0b74f74f92
# ╠═b6167910-ca96-11f1-bc8b-23e368db32cf
# ╟─b6167910-ca96-11f1-a1cf-6b7482c34fb4
# ╠═b6167910-ca96-11f1-b4e8-f5c23bd64a81
# ╟─b616791a-ca96-11f1-8685-71dfef940fcb
# ╠═b616791a-ca96-11f1-b570-d195f6460d2c
# ╟─b6167924-ca96-11f1-bf2d-db786205e345
# ╠═b6167924-ca96-11f1-a53f-f575f2928aa6
# ╟─b616792e-ca96-11f1-aebb-1d9b6dd5a56c
# ╟─b616792e-ca96-11f1-a9cb-e5eee0348df6
# ╟─b6167936-ca96-11f1-847e-2d234c371af4
# ╠═b6167942-ca96-11f1-bf97-2937b97c1039
# ╟─b6167942-ca96-11f1-b7ca-2b71aff44c9d
# ╠═b6167942-ca96-11f1-8ccf-73a0d5374769
# ╟─b616794c-ca96-11f1-ac94-f1444fc53a27
# ╟─b616794c-ca96-11f1-802e-a3300e661634
# ╠═b6167956-ca96-11f1-a3bd-4fcefb85f14e
# ╠═b6167956-ca96-11f1-a936-8df2b346e4b1
# ╠═b6167960-ca96-11f1-9acf-31dd9c82d894
# ╠═b6167960-ca96-11f1-ae5a-81b4a4ea72e5
# ╠═b616796a-ca96-11f1-aa5c-bbedc915d2c8
# ╠═b616796a-ca96-11f1-9418-3713fca7660b
# ╠═b6167974-ca96-11f1-8890-95878475aae7
# ╠═b6167974-ca96-11f1-a1ea-5719f3d778d1
# ╠═b616797e-ca96-11f1-912d-237c96559d9e
# ╠═b6167988-ca96-11f1-8a8a-cdd4f9ed6fd3
# ╠═b6167988-ca96-11f1-a5f3-fb77ec69b2d0
# ╠═b6167988-ca96-11f1-b721-dde064ea24e4
