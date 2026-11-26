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

# ╔═╡ 76345140-ca92-11f1-920d-09ea8533e420
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

# ╔═╡ 76344ff6-ca92-11f1-a45e-07935e3e3590
md"""
# Sesión 3: Colas pesadas y leyes de potencias

**Objetivos de la sesión**
1. Reconocer una ley de potencias y entender por qué es una recta en escala
   log-log.
2. Medir la ley de Zipf en las localidades de México y la ley cúbica inversa
   en los retornos financieros.
3. Distinguir una ley de potencias de una lognormal, y entender por qué no
   basta con "ver una recta".
4. Estimar el exponente de una cola con una receta confiable, el estimador de
   Hill. **[Avanzado]** Elegir ``x_{\min}`` con el método de Clauset, Shalizi y
   Newman.

**Cómo usar este cuaderno.** Como en la sesión 1, Pluto es reactivo: al mover
un deslizador o cambiar un selector, todo lo que depende de ese valor se
recalcula. Las secciones marcadas **[Avanzado]** son opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s03_colas_pesadas.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): P.I5 (leyes de potencias),
C.B2 (escalas logarítmicas), S.B4 (mínimos cuadrados) y, para la ruta
avanzada, P.A2–P.A4.

**Conceptos transversales.** Se completan las **distribuciones, la CCDF y las
escalas logarítmicas**, que reaparecen en las sesiones 6, 8, 9 y 10.
"""

# ╔═╡ 7634502a-ca92-11f1-a8d4-57b0220dbc2a
TableOfContents(title = "Contenido")

# ╔═╡ 76345032-ca92-11f1-94a4-f738e4f94c34
md"""
## 0. El fenómeno: un mundo sin tamaño típico

Piensa en la estatura de los adultos. La persona más alta que conoces mide
quizá 1.2 veces el promedio, y nadie mide diez veces más que otra persona.
Las estaturas tienen un **tamaño típico**: casi todas están cerca del
promedio.

Ahora mira el tamaño de las localidades de México según el Censo de
Población y Vivienda 2020 (INEGI). Los datos están en
`Data/localidades_2020.csv` y se regeneran con `Code/descargar_inegi.jl`.
"""

# ╔═╡ 76345032-ca92-11f1-8043-5ba190907dd5
localidades = CSV.read(joinpath(@__DIR__, "..", "Data", "localidades_2020.csv"), DataFrame)

# ╔═╡ 7634503c-ca92-11f1-91fa-c3ecb13dedbd
pob = Float64.(localidades.poblacion);    # habitantes por localidad (al menos 100)

# ╔═╡ 76345046-ca92-11f1-b85e-034fb8247066
md"""
Las diez localidades más pobladas. Ojo: en la Ciudad de México cada alcaldía
cuenta como una localidad distinta, así que la capital aparece partida.
"""

# ╔═╡ 76345050-ca92-11f1-ae71-73097f4e539a
md"""
## 1. Repaso: por qué una ley de potencias es una recta en log-log

Una **ley de potencias** dice que la probabilidad de superar un valor ``x``
cae como una potencia de ``x``:

```math
P(X > x) = \left(\frac{x}{x_{\min}}\right)^{-\alpha}, \qquad x \ge x_{\min}.
```

Al tomar logaritmos de ambos lados,

```math
\ln P(X > x) = \alpha \ln x_{\min} - \alpha \ln x,
```

que es una **recta** en las variables ``\ln x`` y ``\ln P`` con pendiente
``-\alpha``. En cambio, una exponencial ``P(X > x) = e^{-\lambda x}`` da
``\ln P = -\lambda x``: es recta en escala semilog, no en log-log.

La función ``P(X > x)`` se llama **CCDF** (función de distribución acumulada
complementaria). Para estudiar colas es mejor que el histograma: no depende
de cómo se elijan los intervalos y aprovecha cada dato de la cola.

Exponente ``\alpha``: $(@bind α_demo Slider(0.5:0.1:3.0, default = 1.5, show_value = true))
  Escala: $(@bind escala_demo Select(["log-log", "semilog", "lineal"]))
"""

# ╔═╡ 76345064-ca92-11f1-980a-eb9cf4bb3778
md"""
## 2. Ley de Zipf: las localidades de México

Ordena las localidades de mayor a menor y llama **rango** ``R`` a su lugar en
la lista: la más grande tiene rango 1, la segunda rango 2, etcétera. La
localidad de rango ``R`` tiene exactamente ``R`` localidades de su tamaño o
mayores, así que

```math
R(x) = N \cdot P(X \ge x).
```

Graficar el rango contra el tamaño en log-log es graficar la CCDF, desplazada
por ``\ln N``. Si la cola es una ley de potencias, ``R \propto x^{-\alpha}``,
o bien ``x \propto R^{-1/\alpha}``. La **ley de Zipf** es el caso
``\alpha = 1``: la segunda localidad tendría la mitad de habitantes que la
primera, la tercera un tercio, y así sucesivamente.
"""

# ╔═╡ 76345082-ca92-11f1-a0aa-09ac83e16860
md"""
## 3. ¿Ley de potencias o lognormal?

Una variable es **lognormal** si su logaritmo es normal. Aparece cuando algo
crece multiplicándose por muchos factores aleatorios (sesión 6). Su CCDF en
log-log **no** es una recta: se curva hacia abajo. Pero si la dispersión
``\sigma`` del logaritmo es grande, la curvatura es tan suave que en dos o
tres órdenes de magnitud parece recta.

Mueve ``\sigma`` y compara la CCDF de 10 000 datos lognormales con la recta que
mejor ajusta al 5% más grande de ellos. El eje horizontal abarca siempre los
mismos nueve órdenes de magnitud.

Dispersión ``\sigma``: $(@bind σ_ln Slider(0.5:0.25:3.0, default = 2.0, show_value = true))
  Semilla: $(@bind semilla_ln NumberField(1:1000, default = 1))
"""

# ╔═╡ 7634508e-ca92-11f1-b3ed-b58ed15f9d34
muestra_lognormal = rand(Xoshiro(semilla_ln), LogNormal(0, σ_ln), 10_000);

# ╔═╡ 763450a0-ca92-11f1-afbf-53ab6428bc65
datos = CSV.read(joinpath(@__DIR__, "..", "Data", "precios.csv"), DataFrame);

# ╔═╡ 76345096-ca92-11f1-936d-53b6f8e07494
md"""
## 4. La ley cúbica inversa de los retornos

En la sesión 1 vimos que los retornos tienen colas más pesadas que la normal.
Gopikrishnan, Meyer, Amaral y Stanley (1998) midieron qué tan pesadas, con
datos de miles de acciones, y encontraron

```math
P(|z| > x) \sim x^{-\alpha}, \qquad \alpha \approx 3,
```

donde ``z`` es el retorno estandarizado. Se conoce como la **ley cúbica
inversa**. Veamos si se cumple en nuestros cuatro activos.

**Activo:** $(@bind activo Select(unique(datos.activo)))
  $(@bind separar CheckBox(default = false)) Separar caídas y subidas
"""

# ╔═╡ 763450a0-ca92-11f1-bcd5-b767cbada7a8
begin
    serie = sort(filter(:activo => ==(activo), datos), :fecha)
    r = diff(log.(serie.cierre))           # retornos logarítmicos diarios
    z = (r .- mean(r)) ./ std(r)           # retornos estandarizados
end;

# ╔═╡ 763450aa-ca92-11f1-84c2-53e9b0d7c5ba
md"""
Lejos del centro, los datos caen paralelos a la línea de pendiente ``-3`` y
muy por encima de la normal. Que ``\alpha \approx 3`` significa que al
**duplicar** el umbral, la frecuencia de los eventos se divide entre
``2^3 = 8``. Con la normal se dividiría entre un número astronómico.

Marca la casilla para separar caídas y subidas. ¿Son simétricas las colas en
los cuatro activos?
"""

# ╔═╡ 763450b4-ca92-11f1-bb43-5350b6cb2bac
md"""
## 5. ¿Cómo se mide α?

### 5.1 La tentación de ajustar una recta frente a la receta de Hill

Si la CCDF es una recta en log-log, ¿por qué no ajustar una recta por mínimos
cuadrados y leer la pendiente? Hagamos un experimento: generamos 1000
muestras de 200 datos de una Pareto con ``\alpha = 2`` (sabemos la respuesta)
y comparamos la recta por mínimos cuadrados con la **receta de Hill**, que
promedia los logaritmos de los datos respecto al umbral ``x_{\min}``:

```math
\hat\alpha = \frac{n}{\sum_{i=1}^{n} \ln(x_i / x_{\min})} .
```
"""

# ╔═╡ 763450d2-ca92-11f1-af0f-f12f95d4597e
md"""
### 5.2 El estimador de Hill en la cola

En la práctica, sólo la **cola** sigue una ley de potencias. El **estimador de
Hill** aplica la receta anterior a los ``k`` datos más grandes, usando como
``x_{\min}`` el dato número ``k + 1``:

```math
\hat\alpha_k = \frac{k}{\sum_{i=1}^{k} \ln\left(x_{(i)} / x_{(k+1)}\right)},
```

donde ``x_{(1)} \ge x_{(2)} \ge \dots`` son los datos ordenados de mayor a
menor. El problema es elegir ``k``. La **gráfica de Hill** muestra
``\hat\alpha_k`` para cada ``k``: se busca una zona donde sea estable.

Fracción de datos en la cola: $(@bind frac_k Slider(0.005:0.005:0.20, default = 0.05, show_value = true))
"""

# ╔═╡ 763450e6-ca92-11f1-9670-9ff9b4b0d9de
md"""
### 5.3 [Avanzado] Elegir ``x_{\min}`` con el método de Clauset, Shalizi y Newman

Para elegir ``x_{\min}`` sin hacerlo a ojo, Clauset, Shalizi y Newman (2009)
proponen:

1. Para cada candidato a ``x_{\min}``, estimar ``\alpha`` por máxima
   verosimilitud con los datos mayores o iguales a ``x_{\min}``.
2. Medir qué tan lejos está la cola empírica de la ley de potencias ajustada
   con la **distancia de Kolmogorov-Smirnov**: la mayor diferencia entre las
   dos funciones de distribución (repaso P.A4).
3. Quedarse con el ``x_{\min}`` de menor distancia.
"""

# ╔═╡ 763450ee-ca92-11f1-abb5-177e37336108
md"""
## 6. Universalidad de las colas

*Si alcanza (margen de la clase).* La tabla resume los exponentes estimados para los cuatro activos (por separado
para caídas, subidas y valor absoluto) y para las localidades de México:
"""

# ╔═╡ 763450fa-ca92-11f1-8081-33245d0ba927
md"""
Los exponentes de los índices y del tipo de cambio rondan 3, como predice la
ley cúbica inversa, aunque con intervalos amplios: hay pocos días en la cola.
Bitcoin es el caso raro, con colas muy asimétricas. Las localidades, en
cambio, tienen ``\alpha \approx 1``: una cola mucho más pesada.

Gabaix (2009) reúne muchas leyes de potencias en economía: tamaños de
ciudades y empresas (``\alpha \approx 1``), ingresos y riqueza de los más ricos
(``\alpha`` entre 1.5 y 3), retornos y volumen de los mercados
(``\alpha \approx 3`` y 1.5). Que fenómenos tan distintos compartan
exponentes es una pista de mecanismos comunes, que buscaremos a lo largo del
curso.
"""

# ╔═╡ 76345104-ca92-11f1-a88a-799c2bfb0f85
md"""
## 7. [Avanzado] Cuando el teorema del límite central no aplica

El TLC dice que la suma de muchas variables independientes con **varianza
finita** tiende a la normal. Una ley de potencias con ``\alpha < 2`` tiene
varianza infinita, y entonces la suma tiende a otra distribución: una
**distribución estable** de índice ``\alpha``, con la misma cola
``x^{-\alpha}``. Además, la suma crece como ``n^{1/\alpha}``, más rápido que
``\sqrt{n}`` (repaso P.A2).

La gráfica suma ``n`` variables de Pareto simétricas y normaliza: con
``n^{1/\alpha}`` si ``\alpha < 2``, y con ``\sigma\sqrt{n}`` si ``\alpha > 2``.

Exponente de cada sumando: $(@bind α_est Slider(1.1:0.1:3.5, default = 1.5, show_value = true))
"""

# ╔═╡ 7634510e-ca92-11f1-ba90-a3809feaf26f
md"""
Con ``\alpha < 2``, las tres curvas casi coinciden: sumar no cambia la forma,
sólo la escala (``n^{1/\alpha}``). Es la propiedad que define a una
distribución estable. Con ``\alpha > 2``, el centro se acerca a la normal al
crecer ``n``, pero la cola de ley de potencias sigue ahí, cada vez más lejos.
Los retornos, con ``\alpha \approx 3``, están en el segundo caso: por eso la
normalidad agregacional de la sesión 2 llega, pero lentamente.
"""

# ╔═╡ 76345118-ca92-11f1-9cce-3d90dc3e5338
md"""
## Resumen

| Concepto | Qué dice | Lo que vimos |
|---|---|---|
| Ley de potencias | ``P(X > x) \propto x^{-\alpha}``: recta de pendiente ``-\alpha`` en log-log | Los momentos de orden ``n \ge \alpha`` no existen |
| Ley de Zipf | ``\alpha = 1`` en tamaños | Localidades de México: ``\hat\alpha \approx 1`` |
| Ley cúbica inversa | ``\alpha \approx 3`` en retornos | IPC, USD/MXN y S&P 500 cerca de 3; Bitcoin asimétrico |
| Lognormal | Curva en log-log, pero puede parecer recta | Ver una recta no basta |
| Estimación | Receta de Hill en la cola; elección de ``x_{\min}`` por KS [Avanzado] | La recta por mínimos cuadrados es menos precisa |
"""

# ╔═╡ 76345118-ca92-11f1-966d-0d862228a8a7
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe la función `alpha_hill(x, k)` que devuelva el estimador de Hill
con los ``k`` datos más grandes de `x`, usando la fórmula de la sección 5.2.
Reemplaza `missing` por tu código.
"""

# ╔═╡ 76345120-ca92-11f1-af05-4123b2a59636
function alpha_hill(x, k)
    missing  # reemplaza esta línea
end

# ╔═╡ 7634512c-ca92-11f1-b79f-9d78ddc934fa
md"""
**2.** *Si alcanza (margen de la clase).* Escribe la función `pendiente_loglog(x, y)` que devuelva la pendiente
de la recta de mínimos cuadrados de ``\ln y`` contra ``\ln x``
(repaso S.B4).
"""

# ╔═╡ 7634512c-ca92-11f1-8c69-576e4dab8e0a
function pendiente_loglog(x, y)
    missing  # reemplaza esta línea
end

# ╔═╡ 76345136-ca92-11f1-ba29-cb7e975a61b1
md"""
### Ruta básica: interpretación

3. Si las localidades siguieran exactamente la ley de Zipf (``\alpha = 1``),
   ¿cuántos habitantes tendría la décima localidad comparada con la primera?
   ¿Se cumple en la tabla de la sección 0?
4. Si los días con ``|z| > 3`` ocurren 12 veces al año en un mercado con
   ``\alpha = 3``, ¿cuántas veces al año esperarías días con ``|z| > 6``?
   ¿Cuántas si los retornos fueran normales?
5. En la gráfica de Hill de la sección 5.2, ¿por qué el intervalo de confianza
   se ensancha a la izquierda (``k`` pequeño) y el estimador se aleja de 3 a la
   derecha (``k`` grande)?

### Ruta avanzada

6. Demuestra que si ``P(X > x) = (x/x_{\min})^{-\alpha}``, entonces
   ``\mathbb{E}[X^n]`` es finito si y sólo si ``n < \alpha``, y deriva el
   estimador de máxima verosimilitud de ``\alpha``.
7. **Ley de potencias contra lognormal.** Para las localidades por encima de
   ``x_{\min}``, ajusta también una lognormal truncada por máxima verosimilitud
   y compara ambas con la prueba de Vuong (repaso S.A2). ¿Pueden los datos
   distinguirlas?
8. Con la simulación de la sección 7, mide la mediana de ``|S_n|`` para varios
   ``n`` y verifica que crece como ``n^{1/\alpha}`` cuando ``\alpha < 2``.

## Para la próxima sesión

- Lectura: Requião da Cunha (2021), sección 2.2.7. Complementaria: Gabaix
  (2009), *Power laws in economics and finance*; Thurner, Klimek y Hanel
  (2018), capítulo 3.
- Ruta avanzada: Mantegna y Stanley (1999), capítulos 3 y 4; Requião da Cunha
  (2021), secciones 2.2.5 y 2.2.6; Clauset, Shalizi y Newman (2009),
  *Power-law distributions in empirical data*, secciones 1 a 4.
- En la **sesión 4** pasamos de describir distribuciones a construir **procesos
  en el tiempo**: la caminata aleatoria y las cadenas de Markov.
"""

# ╔═╡ 76345140-ca92-11f1-a717-d1fde816066c
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 7634514a-ca92-11f1-8ebd-f18e9573468e
begin
    const AZUL    = "#2a78d6"   # serie principal
    const NARANJA = "#eb6834"   # segunda serie
    const GRIS    = "#52514e"   # referencias (normal, bandas)
    const RAMPA   = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # azul secuencial
end;

# ╔═╡ 7634505c-ca92-11f1-b631-b98f440718ae
let
    x = exp10.(range(0, 3, length = 400))
    pareto = x .^ (-α_demo)
    λ = log(2) / (2^(1 / α_demo) - 1)          # misma mediana que la Pareto
    expo = exp.(-λ .* (x .- 1))
    xs = escala_demo == "log-log" ? :log10 : :identity
    ys = escala_demo == "lineal" ? :identity : :log10
    xmax = escala_demo == "log-log" ? 1000 : 30
    p = plot(x, pareto; color = AZUL, lw = 2, label = "Ley de potencias, α = $(α_demo)",
             xscale = xs, yscale = ys, xlims = (1, xmax),
             ylims = ys == :log10 ? (1e-8, 1.2) : (0, 1.05), legend = :topright,
             yticks = ys == :log10 ? 10.0 .^ (0:-2:-8) : :auto,
             xlabel = "x", ylabel = "P(X > x)",
             title = "CCDF en escala $(escala_demo) (ambas con la misma mediana)")
    plot!(p, x, expo; color = NARANJA, lw = 2, label = "Exponencial")
end

# ╔═╡ 7634514a-ca92-11f1-808d-d97052afa054
"Número entero con comas como separador de miles."
miles(x) = replace(string(round(Int, x)), r"(?<=\d)(?=(\d{3})+$)" => ",")

# ╔═╡ 7634514a-ca92-11f1-87d2-794c67be8ec8
"Frase sobre qué momentos son finitos para una cola con exponente α."
function momentos(α)
    α <= 1 && return "la media, la varianza y la curtosis son infinitas"
    α <= 2 && return "la media es finita, pero la varianza y la curtosis son infinitas"
    α <= 4 && return "la media y la varianza son finitas, pero la curtosis es infinita"
    "la media, la varianza y la curtosis son finitas"
end

# ╔═╡ 76345064-ca92-11f1-a016-03d6aaed8412
md"""
Con ``\alpha = `` $(α_demo), $(momentos(α_demo)). La razón: el momento
``\mathbb{E}[X^n]`` es finito sólo si ``n < \alpha``.

| Exponente | Media | Varianza | Curtosis | Dónde aparece |
|---|---|---|---|---|
| ``\alpha \le 1`` | infinita | infinita | infinita | tamaño de ciudades (ley de Zipf) |
| ``1 < \alpha \le 2`` | finita | infinita | infinita | riqueza de los más ricos |
| ``2 < \alpha \le 4`` | finita | finita | infinita | retornos financieros (``\alpha \approx 3``) |
| ``\alpha > 4`` | finita | finita | finita | |

La tercera fila explica algo que vimos en la sesión 1: si ``\alpha \approx 3``,
la curtosis de los retornos **no existe**, y por eso su estimación cambia
tanto de una muestra a otra.
"""

# ╔═╡ 76345154-ca92-11f1-9535-4da6fc71626e
"Datos ordenados de mayor a menor y su CCDF empírica, P(X ≥ x)."
function ccdf_empirica(x)
    s = sort(x; rev = true)
    s, (1:length(s)) ./ length(s)
end

# ╔═╡ 763450aa-ca92-11f1-ab3e-1dead9ae70bc
let
    p = plot(xscale = :log10, yscale = :log10, xlims = (0.5, 30), ylims = (1e-4, 1.2),
             xticks = ([0.5, 1, 2, 5, 10, 20], ["0.5", "1", "2", "5", "10", "20"]),
             yticks = 10.0 .^ (0:-1:-4),
             xlabel = "x (desviaciones estándar)", ylabel = "P(|z| ≥ x)", legend = :bottomleft,
             title = "$activo: cola de los retornos estandarizados")
    xr = exp10.(range(log10(0.5), log10(30), length = 200))
    plot!(p, xr, 2 .* ccdf.(Normal(), xr); color = GRIS, ls = :dash, lw = 1.5, label = "Normal")
    if separar
        for (v, c, nombre) in ((-z[z .< 0], NARANJA, "Caídas"), (z[z .> 0], AZUL, "Subidas"))
            s, P = ccdf_empirica(v)
            scatter!(p, s, P; ms = 2.5, msw = 0, color = c, label = nombre)
        end
    else
        s, P = ccdf_empirica(abs.(z))
        scatter!(p, s, P; ms = 2.5, msw = 0, color = AZUL, label = "Datos: |z|")
    end
    P2 = mean(abs.(z) .>= 2)
    xl = exp10.(range(log10(1.5), log10(30), length = 50))
    plot!(p, xl, P2 .* (xl ./ 2) .^ (-3); color = GRIS, lw = 1.5, ls = :dot,
          label = "Pendiente −3 (ley cúbica)")
end

# ╔═╡ 76345104-ca92-11f1-8326-0d148cd4172f
let
    rng = Xoshiro(2026)
    M = 20_000
    p = plot(xscale = :log10, yscale = :log10, xlims = (0.1, 1e3), ylims = (1e-4, 1.2),
             xlabel = "x", ylabel = "P(|suma normalizada| ≥ x)", legend = :bottomleft,
             title = "Sumas de n variables con cola α = $(α_est)")
    for (i, n) in enumerate((1, 10, 100))
        S = [sum(rand(rng, (-1, 1)) * rand(rng)^(-1 / α_est) for _ in 1:n) for _ in 1:M]
        escala = α_est < 2 ? n^(1 / α_est) : sqrt(α_est / (α_est - 2)) * sqrt(n)
        s, P = ccdf_empirica(abs.(S) ./ escala)
        plot!(p, s, P; color = RAMPA[i+1], lw = 2, label = "n = $n")
    end
    if α_est > 2
        xr = exp10.(range(-1, 1, length = 100))
        plot!(p, xr, 2 .* ccdf.(Normal(), xr); color = GRIS, ls = :dash, lw = 1.5, label = "Normal")
    end
    p
end

# ╔═╡ 76345154-ca92-11f1-9496-0fe193f793d5
"Pendiente de la recta de mínimos cuadrados de v contra u."
pendiente(u, v) = cov(u, v) / var(u)

# ╔═╡ 7634508e-ca92-11f1-9a72-037f6353b27a
pendiente_lognormal = let
    s, P = ccdf_empirica(muestra_lognormal)
    cola = s .>= quantile(muestra_lognormal, 0.95)
    pendiente(log.(s[cola]), log.(P[cola]))
end;

# ╔═╡ 76345082-ca92-11f1-b16a-cfdd66b28bdb
let
    x = muestra_lognormal
    s, P = ccdf_empirica(x)
    cola = s .>= quantile(x, 0.95)
    b = pendiente_lognormal
    a = mean(log.(P[cola])) - b * mean(log.(s[cola]))
    p = scatter(s, P; xscale = :log10, yscale = :log10, ms = 2, msw = 0, color = AZUL,
                label = "Datos lognormales (σ = $(σ_ln))", legend = :bottomleft,
                xlabel = "x", ylabel = "P(X ≥ x)", ylims = (5e-5, 1.2), xlims = (1e-4, 1e5),
                xticks = 10.0 .^ (-4:1:5),
                title = "Una lognormal puede parecer ley de potencias")
    xr = exp10.(range(-4, 5, length = 100))
    plot!(p, xr, min.(exp.(a .+ b .* log.(xr)), 10); color = NARANJA, lw = 2, ls = :dash,
          label = "Recta ajustada al 5% superior: pendiente $(round(b, digits = 2))")
end

# ╔═╡ 76345096-ca92-11f1-909b-2ddf89fd4019
md"""
Con ``\sigma`` pequeña, toda la caída ocurre en menos de un orden de magnitud:
nadie la confundiría con una ley de potencias. Con ``\sigma \ge 2``, los datos
se extienden por varios órdenes de magnitud y la parte alta parece una recta.
Con ``\sigma = `` $(σ_ln), la recta ajustada al 5% superior da
``\alpha \approx`` $(round(-pendiente_lognormal, digits = 1)), y alguien podría
reportar "una ley de potencias" aunque los datos son lognormales. Por eso,
antes de afirmar que algo sigue una ley de potencias, hace falta:

1. Que la recta abarque **varios órdenes de magnitud** en ``x``.
2. Estimar ``\alpha`` con un **método correcto** (sección 5), no con una regla.
3. **Comparar** contra alternativas como la lognormal (ruta avanzada).

Clauset, Shalizi y Newman (2009) revisaron muchas "leyes de potencias"
publicadas y encontraron que varias no resisten estas pruebas.
"""

# ╔═╡ 7634515e-ca92-11f1-aa30-3d4401305c32
"Máxima verosimilitud del exponente α de la CCDF, con los datos ≥ xmin."
function alpha_mle(x, xmin)
    cola = filter(>=(xmin), x)
    length(cola) / sum(log.(cola ./ xmin))
end

# ╔═╡ 7634515e-ca92-11f1-8b37-4589605f8179
"Estimador de Hill con los k datos más grandes, y su error estándar."
function hill(x, k)
    s = sort(x; rev = true)
    α = k / sum(log.(s[1:k] ./ s[k+1]))
    (α = α, se = α / sqrt(k))
end

# ╔═╡ 763450dc-ca92-11f1-8e97-25092c91b175
hill_actual = (k = round(Int, frac_k * length(z)), hill(abs.(z), round(Int, frac_k * length(z)))...);

# ╔═╡ 763450dc-ca92-11f1-92da-4f87ba585da3
md"""
Con los $(miles(hill_actual.k)) datos más grandes (``k/N = `` $(frac_k)),
el estimador de Hill da ``\hat\alpha = `` $(round(hill_actual.α, digits = 2)) ± $(round(1.96hill_actual.se, digits = 2)).

Con ``k`` muy pequeño hay pocos datos y mucho ruido; con ``k`` muy grande se
incluye el centro de la distribución, que no es ley de potencias, y el
estimador se sesga.
"""

# ╔═╡ 76345168-ca92-11f1-a303-4bffde54f400
"Estimador de Hill para cada k en `ks`, ordenando los datos una sola vez."
function curva_hill(x, ks)
    s = sort(x; rev = true)
    L = cumsum(log.(s))
    [k / (L[k] - k * log(s[k+1])) for k in ks]
end

# ╔═╡ 763450d2-ca92-11f1-85da-c176444fe982
let
    a = abs.(z)
    N = length(a)
    ks = 10:max(11, round(Int, 0.25N))
    αs = curva_hill(a, ks)
    k = round(Int, frac_k * N)
    p = plot(ks ./ N, αs .+ 1.96 .* αs ./ sqrt.(ks); fillrange = αs .- 1.96 .* αs ./ sqrt.(ks),
             color = AZUL, fillalpha = 0.15, lw = 0, label = "Intervalo de 95%",
             xlabel = "fracción de datos en la cola (k/N)", ylabel = "α estimado",
             title = "$activo: gráfica de Hill para |z|", ylims = (0, 8), legend = :topright)
    plot!(p, ks ./ N, αs; color = AZUL, lw = 2, label = "Hill")
    hline!(p, [3]; color = GRIS, ls = :dot, lw = 1, label = "α = 3")
    vline!(p, [k / N]; color = NARANJA, lw = 1.5, label = "k elegido")
end

# ╔═╡ 76345168-ca92-11f1-88aa-31c98b536685
"""
Método de Clauset, Shalizi y Newman (2009). Para cada candidato a `xmin` estima
α por máxima verosimilitud y calcula la distancia de Kolmogorov-Smirnov entre
la cola empírica y la ley de potencias ajustada. Devuelve todas las curvas y,
en `mejor`, el `xmin` de menor distancia.
"""
function clauset(x; candidatos = quantile(x, 0.5:0.005:0.995), nmin = 50)
    s = sort(x; rev = true)
    L = cumsum(log.(s))
    res = (xmin = Float64[], α = Float64[], D = Float64[], n = Int[])
    for xm in unique(candidatos)
        n = searchsortedlast(s, xm; rev = true)        # datos mayores o iguales a xm
        n < nmin && continue
        α = n / (L[n] - n * log(xm))
        cola = @view s[n:-1:1]                         # cola en orden creciente
        D = maximum(abs.((1:n) ./ n .- (1 .- (xm ./ cola) .^ α)))
        push!(res.xmin, xm); push!(res.α, α); push!(res.D, D); push!(res.n, n)
    end
    i = argmin(res.D)
    (; res..., mejor = (xmin = res.xmin[i], α = res.α[i], D = res.D[i], n = res.n[i],
                        se = res.α[i] / sqrt(res.n[i])))
end

# ╔═╡ 7634506e-ca92-11f1-9b99-9195d39f1ec9
ajuste_loc = clauset(pob; candidatos = quantile(pob, 0.5:0.0025:0.999)).mejor;

# ╔═╡ 7634506e-ca92-11f1-92e3-4be9682a92a8
let
    s = sort(pob; rev = true)
    R = 1:length(s)
    p = plot(s, R; xscale = :log10, yscale = :log10, color = AZUL, lw = 2,
             xlabel = "habitantes x", ylabel = "rango R", yticks = 10.0 .^ (0:5),
             title = "Rango contra tamaño de las localidades de México (Censo 2020)",
             xlims = (100, 3e6), ylims = (0.8, 1e5), label = "Datos", legend = :topright)
    xm, α = ajuste_loc.xmin, ajuste_loc.α
    nxm = ajuste_loc.n
    xr = exp10.(range(log10(xm), log10(2e6), length = 50))
    plot!(p, xr, nxm .* (xr ./ xm) .^ (-α); color = NARANJA, lw = 2, ls = :dash,
          label = "Ley de potencias ajustada, α = $(round(α, digits = 2))")
    vspan!(p, [100, xm]; color = GRIS, alpha = 0.08, label = "Por debajo de x_min")
    plot!(p, xr, nxm .* (xr ./ xm) .^ (-1.0); color = GRIS, lw = 1, ls = :dot,
          label = "Zipf exacta (α = 1)")
end

# ╔═╡ 76345078-ca92-11f1-bba7-c702dcf8ecf5
md"""
Por encima de ``x_{\min} \approx`` $(miles(round(Int, ajuste_loc.xmin))) habitantes
(el método para elegirlo está en la sección 5), el exponente estimado es
``\hat\alpha = `` $(round(ajuste_loc.α, digits = 2)) ± $(round(1.96ajuste_loc.se, digits = 2)),
con $(miles(ajuste_loc.n)) localidades en la cola. Es la ley de Zipf casi exacta.

Dos cuidados:
- Una **localidad** no es una ciudad: la Ciudad de México está partida en
  alcaldías, y Monterrey o Guadalajara en varios municipios. Con zonas
  metropolitanas la parte alta de la gráfica cambia. (Es una buena idea de
  proyecto final.)
- La gráfica dice **qué** pasa, no **por qué**. En la sesión 6 veremos un
  mecanismo que genera leyes de potencias: el crecimiento multiplicativo
  (Gibrat y Kesten).
"""

# ╔═╡ 763450e6-ca92-11f1-b42a-fdbe6aab74f8
ajuste_activo = clauset(abs.(z));

# ╔═╡ 763450e6-ca92-11f1-87ec-51c079380f5c
let
    m = ajuste_activo.mejor
    p1 = plot(ajuste_activo.xmin, ajuste_activo.D; color = AZUL, lw = 2,
              ylabel = "distancia KS", title = "$activo: elección de x_min para |z|",
              xlabel = "")
    vline!(p1, [m.xmin]; color = NARANJA, lw = 1.5)
    p2 = plot(ajuste_activo.xmin, ajuste_activo.α; color = AZUL, lw = 2,
              ylabel = "α estimado", xlabel = "x_min (desviaciones estándar)")
    vline!(p2, [m.xmin]; color = NARANJA, lw = 1.5)
    hline!(p2, [3]; color = GRIS, ls = :dot, lw = 1)
    plot(p1, p2; layout = (2, 1), size = (720, 480))
end

# ╔═╡ 763450ee-ca92-11f1-86cf-916b4a072df9
md"""
La distancia es mínima en ``x_{\min} = `` $(round(ajuste_activo.mejor.xmin, digits = 2))
desviaciones estándar. Ahí, ``\hat\alpha = `` $(round(ajuste_activo.mejor.α, digits = 2)) ± $(round(1.96ajuste_activo.mejor.se, digits = 2)),
con $(ajuste_activo.mejor.n) días en la cola.
"""

# ╔═╡ 76345172-ca92-11f1-8bb9-15f15c256caa
"Resumen del tamaño de las localidades para el texto de la sección 0."
resumen_loc = let s = sort(pob; rev = true), i = argmax(localidades.poblacion)
    (n = length(s), mediana = median(s), maximo = s[1],
     mayor = "$(localidades.nombre_localidad[i]) ($(localidades.nombre_entidad[i]))",
     razon = s[1] / median(s),
     top1 = round(100 * sum(s[1:ceil(Int, 0.01length(s))]) / sum(s), digits = 1))
end;

# ╔═╡ 76345046-ca92-11f1-834e-e954e652d852
md"""
- Hay $(miles(resumen_loc.n)) localidades con al menos 100 habitantes.
- La mediana tiene $(miles(resumen_loc.mediana)) habitantes. La más grande,
  $(resumen_loc.mayor), tiene $(miles(resumen_loc.maximo)): unas
  $(miles(resumen_loc.razon)) veces la mediana.
- El 1% de localidades más grandes concentra el $(resumen_loc.top1)% de la
  población de todas ellas.

No hay un "tamaño típico" de localidad: hay muchísimas pequeñas y unas pocas
enormes, y las enormes concentran a buena parte de la población. Ese es el
sello de una **cola pesada**. En la sesión 1 vimos que los retornos
financieros también la tienen. Hoy vamos a medir **qué tan pesada** es.
"""

# ╔═╡ 76345172-ca92-11f1-9df6-5d0c6e90a7d1
function tabla_mayores(loc; n = 10)
    t = first(sort(loc, :poblacion; rev = true), n)
    DataFrame("Rango" => 1:n, "Localidad" => t.nombre_localidad,
              "Estado" => t.nombre_entidad, "Habitantes" => miles.(t.poblacion))
end

# ╔═╡ 76345050-ca92-11f1-94b1-0d61e3f53c47
tabla_mayores(localidades)

# ╔═╡ 76345172-ca92-11f1-a876-61ef528d8e58
"Compara la recta por mínimos cuadrados y la máxima verosimilitud en muestras de Pareto."
function experimento_estimadores(; α = 2.0, n = 200, M = 1000, rng = Xoshiro(7))
    ols, mle = Float64[], Float64[]
    for _ in 1:M
        x = rand(rng, Pareto(α, 1.0), n)
        s, P = ccdf_empirica(x)
        push!(ols, -pendiente(log.(s), log.(P)))
        push!(mle, alpha_mle(x, 1.0))
    end
    (; ols, mle)
end

# ╔═╡ 763450bc-ca92-11f1-8544-e743c64b8f30
experimento = experimento_estimadores(α = 2.0, n = 200, M = 1000);

# ╔═╡ 763450bc-ca92-11f1-9206-bb07d4f07ce0
let
    bordes = range(1.0, 3.2, length = 45)
    p = histogram(experimento.ols; bins = bordes, color = NARANJA, alpha = 0.55, lw = 0,
                  label = "Recta por mínimos cuadrados", legend = :topright,
                  xlabel = "α estimado", ylabel = "número de muestras",
                  title = "1000 muestras de una Pareto con α = 2 (200 datos cada una)")
    histogram!(p, experimento.mle; bins = bordes, color = AZUL, alpha = 0.55, lw = 0,
               label = "Receta de Hill (máxima verosimilitud)")
    vline!(p, [2.0]; color = GRIS, lw = 1.5, ls = :dash, label = "Valor verdadero")
end

# ╔═╡ 763450c8-ca92-11f1-9749-71b280f16d8f
md"""
| Método | Promedio de ``\hat\alpha`` | Desviación estándar |
|---|---|---|
| Recta por mínimos cuadrados | $(round(mean(experimento.ols), digits = 2)) | $(round(std(experimento.ols), digits = 2)) |
| Receta de Hill (máxima verosimilitud) | $(round(mean(experimento.mle), digits = 2)) | $(round(std(experimento.mle), digits = 2)) |

La recta es más imprecisa: los puntos de la CCDF no son independientes (cada
uno incluye a todos los que están a su derecha) y los últimos, que tienen muy
pocos datos, pesan lo mismo que los demás. La receta de Hill (que en la ruta
avanzada se obtiene por máxima verosimilitud, repaso P.A3) es más precisa y su
error estándar es ``\hat\alpha/\sqrt{n}``.
"""

# ╔═╡ 7634517c-ca92-11f1-bd6a-2de581945e90
function tabla_colas(datos)
    filas = map(collect(groupby(datos, :activo))) do g
        r = diff(log.(sort(g, :fecha).cierre))
        z = (r .- mean(r)) ./ std(r)
        caidas, subidas, todo = (clauset(v).mejor for v in (-z[z .< 0], z[z .> 0], abs.(z)))
        (g.activo[1], round(caidas.α, digits = 2), round(subidas.α, digits = 2),
         round(todo.α, digits = 2), round(1.96todo.se, digits = 2),
         round(todo.xmin, digits = 2), todo.n)
    end
    DataFrame([getindex.(filas, j) for j in 1:7],
              ["Activo", "α caídas", "α subidas", "α de |z|", "± (95%)",
               "x_min (σ)", "Días en la cola"])
end

# ╔═╡ 763450fa-ca92-11f1-b6db-9161f01da4db
tabla_colas(datos)

# ╔═╡ 7634517c-ca92-11f1-a4b5-03ace847a396
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

# ╔═╡ 76345136-ca92-11f1-9338-11b663a0df6c
verificar(() -> pendiente_loglog(collect(1.0:100.0), 3 .* collect(1.0:100.0) .^ (-1.5)), -1.5;
          pista = md"La pendiente de mínimos cuadrados es `cov(u, v) / var(u)`, con `u = log.(x)` y `v = log.(y)`.")

# ╔═╡ 76345186-ca92-11f1-89a1-8bcdfec7bb0e
const VECTOR_PARETO = rand(Xoshiro(2026), Pareto(2.5, 1.0), 10_000);

# ╔═╡ 76345120-ca92-11f1-8a77-7bfad20ed01f
verificar(() -> alpha_hill(VECTOR_PARETO, 500), hill(VECTOR_PARETO, 500).α;
          pista = md"Ordena con `sort(x; rev = true)`. Los ``k`` mayores son `s[1:k]` y ``x_{(k+1)}`` es `s[k+1]`. Después aplica la fórmula con `sum` y `log.`.")

# ╔═╡ Cell order:
# ╟─76344ff6-ca92-11f1-a45e-07935e3e3590
# ╟─7634502a-ca92-11f1-a8d4-57b0220dbc2a
# ╟─76345032-ca92-11f1-94a4-f738e4f94c34
# ╠═76345032-ca92-11f1-8043-5ba190907dd5
# ╠═7634503c-ca92-11f1-91fa-c3ecb13dedbd
# ╟─76345046-ca92-11f1-834e-e954e652d852
# ╟─76345046-ca92-11f1-b85e-034fb8247066
# ╟─76345050-ca92-11f1-94b1-0d61e3f53c47
# ╟─76345050-ca92-11f1-ae71-73097f4e539a
# ╟─7634505c-ca92-11f1-b631-b98f440718ae
# ╟─76345064-ca92-11f1-a016-03d6aaed8412
# ╟─76345064-ca92-11f1-980a-eb9cf4bb3778
# ╟─7634506e-ca92-11f1-9b99-9195d39f1ec9
# ╟─7634506e-ca92-11f1-92e3-4be9682a92a8
# ╟─76345078-ca92-11f1-bba7-c702dcf8ecf5
# ╟─76345082-ca92-11f1-a0aa-09ac83e16860
# ╟─76345082-ca92-11f1-b16a-cfdd66b28bdb
# ╟─7634508e-ca92-11f1-b3ed-b58ed15f9d34
# ╟─7634508e-ca92-11f1-9a72-037f6353b27a
# ╟─76345096-ca92-11f1-909b-2ddf89fd4019
# ╟─76345096-ca92-11f1-936d-53b6f8e07494
# ╠═763450a0-ca92-11f1-afbf-53ab6428bc65
# ╠═763450a0-ca92-11f1-bcd5-b767cbada7a8
# ╟─763450aa-ca92-11f1-ab3e-1dead9ae70bc
# ╟─763450aa-ca92-11f1-84c2-53e9b0d7c5ba
# ╟─763450b4-ca92-11f1-bb43-5350b6cb2bac
# ╟─763450bc-ca92-11f1-8544-e743c64b8f30
# ╟─763450bc-ca92-11f1-9206-bb07d4f07ce0
# ╟─763450c8-ca92-11f1-9749-71b280f16d8f
# ╟─763450d2-ca92-11f1-af0f-f12f95d4597e
# ╟─763450d2-ca92-11f1-85da-c176444fe982
# ╟─763450dc-ca92-11f1-92da-4f87ba585da3
# ╟─763450dc-ca92-11f1-8e97-25092c91b175
# ╟─763450e6-ca92-11f1-9670-9ff9b4b0d9de
# ╟─763450e6-ca92-11f1-b42a-fdbe6aab74f8
# ╟─763450e6-ca92-11f1-87ec-51c079380f5c
# ╟─763450ee-ca92-11f1-86cf-916b4a072df9
# ╟─763450ee-ca92-11f1-abb5-177e37336108
# ╟─763450fa-ca92-11f1-b6db-9161f01da4db
# ╟─763450fa-ca92-11f1-8081-33245d0ba927
# ╟─76345104-ca92-11f1-a88a-799c2bfb0f85
# ╟─76345104-ca92-11f1-8326-0d148cd4172f
# ╟─7634510e-ca92-11f1-ba90-a3809feaf26f
# ╟─76345118-ca92-11f1-9cce-3d90dc3e5338
# ╟─76345118-ca92-11f1-966d-0d862228a8a7
# ╠═76345120-ca92-11f1-af05-4123b2a59636
# ╟─76345120-ca92-11f1-8a77-7bfad20ed01f
# ╟─7634512c-ca92-11f1-b79f-9d78ddc934fa
# ╠═7634512c-ca92-11f1-8c69-576e4dab8e0a
# ╟─76345136-ca92-11f1-9338-11b663a0df6c
# ╟─76345136-ca92-11f1-ba29-cb7e975a61b1
# ╟─76345140-ca92-11f1-a717-d1fde816066c
# ╠═76345140-ca92-11f1-920d-09ea8533e420
# ╠═7634514a-ca92-11f1-8ebd-f18e9573468e
# ╠═7634514a-ca92-11f1-808d-d97052afa054
# ╠═7634514a-ca92-11f1-87d2-794c67be8ec8
# ╠═76345154-ca92-11f1-9535-4da6fc71626e
# ╠═76345154-ca92-11f1-9496-0fe193f793d5
# ╠═7634515e-ca92-11f1-aa30-3d4401305c32
# ╠═7634515e-ca92-11f1-8b37-4589605f8179
# ╠═76345168-ca92-11f1-a303-4bffde54f400
# ╠═76345168-ca92-11f1-88aa-31c98b536685
# ╠═76345172-ca92-11f1-8bb9-15f15c256caa
# ╠═76345172-ca92-11f1-9df6-5d0c6e90a7d1
# ╠═76345172-ca92-11f1-a876-61ef528d8e58
# ╠═7634517c-ca92-11f1-bd6a-2de581945e90
# ╠═7634517c-ca92-11f1-a4b5-03ace847a396
# ╠═76345186-ca92-11f1-89a1-8bcdfec7bb0e
