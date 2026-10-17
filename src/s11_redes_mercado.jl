#%% hide
md"""
# Sesión 11: Redes: el mercado y la robustez

**Objetivos de la sesión**
1. Construir la **red de umbral** de correlaciones del S&P 100 (datos de la
   sesión 7) y compararla con los modelos de Erdős-Rényi y Barabási-Albert.
2. Medir la **robustez** de una red: qué pasa con el componente gigante cuando
   fallan nodos al azar frente a cuando se ataca a los más conectados.
3. **[Avanzado]** Analizar el **Laplaciano** (vector de Fiedler), el **modelo
   de configuración** con la **asortatividad**, y el **umbral de percolación**
   de Molloy-Reed y Cohen et al.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s11_redes_mercado.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): F.B4 (transiciones de fase).
Para las notas formales: L.I2 (matrices simétricas y teorema espectral), L.I5
(distancia entre activos). Ruta avanzada: L.A1 (el Laplaciano de una red),
F.A5 (percolación).
"""

#%% hide
TableOfContents(title = "Contenido")

#%% hide
md"""
## 0. El fenómeno: ¿qué pasa si cae el nodo más conectado?

En la sesión 10 vimos que la conexión preferencial (Barabási-Albert) crea
**concentradores** (*hubs*): unos pocos nodos conectados con buena parte de la
red, como los grandes bancos en el mercado interbancario o los aeropuertos
principales en una red aérea.

Mira esta red de 45 nodos con concentradores. A la izquierda quitamos **6 nodos
al azar**; a la derecha quitamos los **6 nodos más conectados**. En cada panel,
los nodos en color forman el componente gigante que queda, y los grises son
fragmentos que se separaron:
"""

#%% hide
let
    e = ejemplo_ataque_inicial
    p1 = dibujar_red(e.g_azar, e.X_azar; tam = 4.5, colores = e.col_azar, size = (440, 360),
                     title = "Quitar 6 nodos al azar: gigante = $(e.gig_azar)%")
    p2 = dibujar_red(e.g_ataque, e.X_ataque; tam = 4.5, colores = e.col_ataque, size = (440, 360),
                     title = "Quitar los 6 mayores hubs: gigante = $(e.gig_ataque)%")
    plot(p1, p2; layout = (1, 2), size = (880, 360))
end

#%% hide
md"""
Quitar 6 nodos al azar apenas toca la red: el componente gigante conserva el
$(ejemplo_ataque_inicial.gig_azar)% de los nodos. Quitar los 6 concentradores la
rompe en $(ejemplo_ataque_inicial.piezas_ataque) pedazos, y el más grande
conserva sólo el $(ejemplo_ataque_inicial.gig_ataque)%. Hoy mediremos primero la estructura de una red
financiera real y luego pondremos a prueba la **robustez** de las redes ante
fallas aleatorias y ataques dirigidos.
"""

#%% hide
md"""
## 1. Laboratorio: la red del mercado

En la sesión 7 construimos el árbol de expansión mínima de 95 acciones del
S&P 100. Un árbol no tiene triángulos por construcción, así que su
agrupamiento es cero: no sirve para compararlo con los modelos de la sesión 10.
Usemos otra red, más simple: conectamos dos acciones si su correlación supera un
**umbral** ``u``,

```math
A_{ij} = 1 \quad\text{si}\quad \rho_{ij} > u \quad (i \ne j) .
```
"""

#%%
begin
    precios = CSV.read(joinpath(@__DIR__, "..", "Data", "sp100_precios.csv"), DataFrame)
    sectores = CSV.read(joinpath(@__DIR__, "..", "Data", "sp100_sectores.csv"), DataFrame)
    ancho = unstack(precios, :fecha, :ticker, :cierre)    # una columna por acción
    tickers = names(ancho)[2:end]
    R = diff(log.(Matrix{Float64}(ancho[:, 2:end])); dims = 1)   # rendimientos: T × N
    C = cor(R)
    N = size(C, 1)
end;

#%% hide
md"""
Umbral de correlación ``u``: $(@bind u_mercado Slider([0.4, 0.45, 0.5, 0.55, 0.6], default = 0.5, show_value = true))
"""

#%%
mercado = SimpleGraph((C .> u_mercado) .& .!I(N));

#%% hide
let
    sector = [sectores.sector[findfirst(==(t), sectores.ticker)] for t in tickers]
    destacados = ["Finanzas", "Tecnología", "Servicios públicos"]
    paleta = Dict("Finanzas" => AZUL, "Tecnología" => NARANJA, "Servicios públicos" => VERDE)
    colores = [get(paleta, s, GRIS_CLARO) for s in sector]
    kk = degree(mercado)
    p = dibujar_red(mercado, posiciones_mercado; tam = 3 .+ 0.35 .* kk, colores,
                    size = (720, 520), alpha_enlaces = 0.12,
                    title = "Acciones con correlación mayor que $u_mercado: $(ne(mercado)) enlaces")
    for s in destacados
        scatter!(p, [NaN], [NaN]; color = paleta[s], ms = 6, msw = 0, label = s)
    end
    scatter!(p, [NaN], [NaN]; color = GRIS_CLARO, ms = 6, msw = 0, label = "Otros sectores")
    for i in sortperm(kk; rev = true)[1:5]
        annotate!(p, posiciones_mercado[i, 1], posiciones_mercado[i, 2] + 0.06, text(tickers[i], 8, GRIS))
    end
    plot!(p; legend = :bottomleft)
end

#%% hide
md"""
Las posiciones vienen de las distancias ``d_{ij} = \sqrt{2(1-\rho_{ij})}`` de la
sesión 7 (repaso L.I5): acciones cercanas en la gráfica se mueven juntas. El
tamaño de cada punto es su grado; las cinco acciones con más enlaces tienen su
ticker. Los bancos y gestores de inversión forman el núcleo más denso.

Ahora la comparación: dos redes aleatorias con el **mismo número de nodos y
de enlaces** que la red del mercado. (En Barabási-Albert el número de enlaces
es aproximado: cada nodo nuevo llega con un número entero de enlaces.)
"""

#%%
comparacion = let
    m = ne(mercado)
    er_m = erdos_renyi(N, m; seed = 1)
    ba_m = barabasi_albert(N, max(1, round(Int, m / N)); seed = 1)
    [("Mercado", mercado), ("Erdős-Rényi", er_m), ("Barabási-Albert", ba_m)]
end;

#%% hide
let
    filas = map(comparacion) do (nombre, g)
        kk = degree(g)
        "| $nombre | $(ne(g)) | $(round(mean(kk), digits = 1)) | $(maximum(kk)) | " *
        "$(round(global_clustering_coefficient(g), digits = 2)) | " *
        "$(round(distancia_promedio(g), digits = 2)) | $(round(Int, 100 * componente_gigante(g)))% |"
    end
    Markdown.parse("""
    | Red | Enlaces | Grado promedio | Grado máximo | Agrupamiento | Distancia promedio | Componente gigante |
    |---|---|---|---|---|---|---|
    $(join(filas, "\n"))
    """)
end

#%% hide
let
    p = plot(; xlabel = "grado k", ylabel = "P(K ≥ k)", legend = :topright, ylims = (0, 1.02),
             title = "Distribución de grados: mercado contra los modelos")
    for ((nombre, g), color, m) in zip(comparacion, (AZUL, NARANJA, GRIS), (:circle, :diamond, :square))
        s, P = ccdf_empirica(degree(g))
        plot!(p, s, P; seriestype = :steppost, color, lw = 2, label = nombre)
    end
    p
end

#%% hide
md"""
Con ``u = 0.5``:

- El **agrupamiento** del mercado es varias veces el de las dos redes
  aleatorias. Si JPM está correlacionada con BAC y con C, BAC y C también lo
  están: los sectores forman grupos cerrados. Ni Erdős-Rényi ni
  Barabási-Albert lo reproducen.
- Los **grados** son mucho más desiguales que en Erdős-Rényi: hay acciones
  con más de 30 enlaces y varias aisladas. Pero tampoco es una cola de
  Pareto limpia como la de Barabási-Albert: con 95 nodos no hay espacio para
  una cola larga.
- Una parte de las acciones queda **fuera** del componente gigante: casi
  todas son defensivas (farmacéuticas, telefónicas, supermercados) o tienen
  una dinámica muy propia (Tesla, Netflix), poco correlacionada con el
  resto.

Mueve el umbral. Con ``u`` bajo casi todo está conectado con todo; con ``u``
alto la red se rompe en islas por sector. El agrupamiento alto se mantiene
en todo el rango: es una propiedad de los datos, no del umbral.
"""

#%% hide
md"""
## 2. Laboratorio: robustez ante fallas y ataques

¿Qué tan frágil es una red? Quitemos nodos y midamos qué fracción de los
``N`` nodos originales sigue en el componente gigante (repaso F.B4). Dos
escenarios:

- **Fallas aleatorias**: los nodos se quitan en orden al azar (una empresa
  quiebra por razones propias).
- **Ataque dirigido**: se quitan primero los de mayor grado (un regulador o
  un choque que pega justo a los más conectados).

Usamos redes de Erdős-Rényi y de Barabási-Albert con 10,000 nodos y
``\langle k\rangle = 4``.
"""

#%%
begin
    ba = barabasi_albert(10_000, 2; seed = 1)
    er = erdos_renyi(10_000, ne(ba); seed = 1)
    robustez = Dict(
        (nombre, escenario) => curva_robustez(g, escenario == "Fallas" ? randperm(Xoshiro(3), nv(g)) :
                                                  sortperm(degree(g); rev = true))
        for (nombre, g) in (("Erdős-Rényi", er), ("Barabási-Albert", ba)), escenario in ("Fallas", "Ataque"))
end;

#%% hide
let
    f = (0:nv(er)) ./ nv(er)
    p = plot(; xlabel = "fracción de nodos eliminados f", ylabel = "fracción en el componente gigante",
             xlims = (0, 1), ylims = (0, 1.02), legend = :topright,
             title = "Robustez de dos redes con ⟨k⟩ = 4 (N = 10,000)")
    for (nombre, color) in (("Erdős-Rényi", NARANJA), ("Barabási-Albert", AZUL))
        plot!(p, f, robustez[(nombre, "Fallas")]; color, lw = 2, label = "$nombre, fallas al azar")
        plot!(p, f, robustez[(nombre, "Ataque")]; color, lw = 2, ls = :dash, label = "$nombre, ataque")
    end
    p
end

#%% hide
md"""
- Ante **fallas al azar**, las dos redes se degradan poco a poco y casi
  igual. La de Barabási-Albert conserva un componente gigante un poco más
  de tiempo, hasta que se quita más del 80% de los nodos: al azar casi
  siempre se quita un nodo de grado bajo, y los concentradores mantienen
  unido lo que queda. En una red libre de escala infinita no habría umbral
  (sección 3.3).
- Ante un **ataque**, la diferencia es enorme. La red de Barabási-Albert se
  rompe con apenas el 15% de sus nodos: al quitar los concentradores
  desaparece justo lo que la mantenía conectada. La de Erdős-Rényi resiste
  hasta cerca del 37%, porque no tiene nodos especiales.

Albert, Jeong y Barabási (2000) mostraron este contraste en redes reales,
como internet: las redes libres de escala son **robustas ante fallas pero
frágiles ante ataques**. En la sesión 12 veremos la versión
financiera: una red interbancaria con bancos muy conectados resiste bien
los choques pequeños y dispersos, y es vulnerable cuando cae un banco
central.

**Ahora con el mercado.** ¿Qué pasa si quitamos primero las acciones más
conectadas de la red de correlaciones?
"""

#%% hide
let
    f = (0:N) ./ N
    azar = mean(curva_robustez(mercado, randperm(Xoshiro(s), N)) for s in 1:200)
    ataque = curva_robustez(mercado, sortperm(degree(mercado); rev = true))
    plot(f, azar; color = AZUL, lw = 2, label = "Fallas al azar (promedio de 200)",
         xlabel = "fracción de acciones eliminadas", ylabel = "fracción en el componente gigante",
         ylims = (0, 1.02), xlims = (0, 1), legend = :topright,
         title = "Red del mercado con u = $u_mercado")
    plot!(f, ataque; color = NARANJA, lw = 2, ls = :dash, label = "Ataque a las más conectadas")
end

#%% hide
md"""
Quitar las acciones más conectadas (casi todas financieras) parte la red
mucho antes que quitarlas al azar. En la red de correlaciones, "quitar" una
acción no es una quiebra: significa que, sin los bancos y los gestores de
inversión, el resto del mercado se ve como grupos sectoriales casi
independientes. El sector financiero es el que **liga** a los demás.
"""

#%% hide
md"""
## 3. [Avanzado] Espectro, modelo de configuración y percolación

### 3.1 El Laplaciano

El **Laplaciano** de la red es ``L = D - A``, con ``D`` la matriz diagonal
de grados (repaso L.I2 y L.A1). Sus eigenvalores son reales y no negativos, y:

- el número de eigenvalores **iguales a cero** es el número de componentes
  conexas;
- el segundo eigenvalor más pequeño, ``\lambda_2``, mide qué tan difícil es
  partir la red en dos; su eigenvector, el **vector de Fiedler**, sugiere la
  partición: los nodos con componente positiva de un lado, los negativos del
  otro.
"""

#%%
espectro = let
    λ = eigvals(Symmetric(Float64.(Matrix(laplacian_matrix(mercado)))))
    gigante = connected_components(mercado)[argmax(length.(connected_components(mercado)))]
    sub, _ = induced_subgraph(mercado, gigante)
    F = eigen(Symmetric(Float64.(Matrix(laplacian_matrix(sub)))))
    (λ = λ, ceros = count(<(1e-9), λ), componentes = length(connected_components(mercado)),
     gigante = gigante, fiedler = F.vectors[:, 2])
end;

#%% hide
let
    sector = [sectores.sector[findfirst(==(tickers[i]), sectores.ticker)] for i in espectro.gigante]
    nombres = sort(unique(sector))
    pos = [count(sector .== s .&& espectro.fiedler .> 0) for s in nombres]
    neg = [count(sector .== s .&& espectro.fiedler .<= 0) for s in nombres]
    orden = sortperm(pos .- neg)
    p = bar(1:length(nombres), pos[orden]; color = AZUL, lw = 0, label = "Componente positiva",
            xticks = (1:length(nombres), nombres[orden]), xrotation = 35, ylabel = "acciones",
            legend = :topleft, bottom_margin = 10Plots.mm, size = (760, 360),
            ylims = (-maximum(neg) - 1.5, maximum(pos) + 1.5),
            title = "Vector de Fiedler del componente gigante (u = $u_mercado)")
    bar!(p, 1:length(nombres), -neg[orden]; color = NARANJA, lw = 0, label = "Componente negativa")
    hline!(p, [0]; color = GRIS, lw = 1, label = "")
end

#%% hide
md"""
La red tiene $(espectro.componentes) componentes conexas (contando las
acciones aisladas) y el Laplaciano tiene exactamente $(espectro.ceros)
eigenvalores iguales a cero. El vector de Fiedler del componente gigante
separa a las acciones **defensivas** (servicios públicos, consumo básico,
farmacéuticas) de las **cíclicas** (tecnología, finanzas, industria): la
misma división que el segundo eigenvector de la matriz de correlación en la
sesión 7, encontrada ahora sólo con la red.

### 3.2 Modelo de configuración y asortatividad

¿El agrupamiento alto del mercado se debe a que algunas acciones tienen
muchos enlaces? Para saberlo, comparamos con un **modelo de configuración**:
una red al azar con **exactamente los mismos grados** que la real. Se
obtiene recableando la red real muchas veces: se toman dos enlaces
``a{-}b`` y ``c{-}d`` y se cambian por ``a{-}d`` y ``c{-}b``, lo que no
cambia el grado de nadie.

La **asortatividad** es la correlación entre los grados de los dos extremos
de cada enlace. Si es positiva, los nodos muy conectados se conectan entre
sí; si es negativa, los concentradores se conectan con nodos de grado bajo.
"""

#%%
recableada = recablear(mercado, 20ne(mercado));

#%% hide
let
    ag(g) = round(global_clustering_coefficient(g), digits = 2)
    as(g) = round(assortativity(g), digits = 2)
    ba_m = comparacion[3][2]
    Markdown.parse("""
    | Red | Agrupamiento | Asortatividad |
    |---|---|---|
    | Mercado | $(ag(mercado)) | $(as(mercado)) |
    | Configuración (mismos grados) | $(ag(recableada)) | $(as(recableada)) |
    | Barabási-Albert | $(ag(ba_m)) | $(as(ba_m)) |
    """)
end

#%% hide
md"""
Con los mismos grados, el agrupamiento baja de forma notable: la
heterogeneidad de los grados explica una parte, y el resto es estructura de
sectores. La asortatividad del mercado es positiva, como en las redes
sociales: las acciones muy conectadas se conectan entre sí (el núcleo
financiero). En la red recableada y en Barabási-Albert es negativa: con
grados tan desiguales y una sola arista por par, los concentradores no
tienen suficientes concentradores con quién conectarse.

### 3.3 El umbral de percolación

Quitar nodos al azar con probabilidad ``f`` es un problema de **percolación**
(repaso F.A5). Para una red con distribución de grados dada y sin otra
estructura, el componente gigante existe si un vecino típico tiene al menos
otro vecino más (criterio de Molloy y Reed):

```math
\kappa = \frac{\langle k^2\rangle}{\langle k\rangle} > 2 .
```

Al quitar una fracción ``f`` de los nodos al azar, ``\kappa`` cambia y el
componente gigante desaparece en

```math
f_c = 1 - \frac{1}{\kappa - 1} .
```

Para Erdős-Rényi, ``\kappa = c + 1`` y ``f_c = 1 - 1/c``. Para una red libre
de escala con ``P(K \ge k) \sim k^{-2}``, ``\langle k^2\rangle`` diverge al
crecer ``N`` y ``f_c \to 1``: **no hay umbral**. La derivación está en las
notas formales.
"""

#%% hide
let
    filas = map((("Erdős-Rényi", er), ("Barabási-Albert", ba))) do (nombre, g)
        kk = degree(g)
        κ = mean(kk .^ 2) / mean(kk)
        S = robustez[(nombre, "Fallas")]
        f_sim = (findfirst(<(0.01), S) - 1) / nv(g)
        "| $nombre | $(round(κ, digits = 1)) | $(round(1 - 1 / (κ - 1), digits = 2)) | $(round(f_sim, digits = 2)) |"
    end
    Markdown.parse("""
    | Red (N = 10,000) | ``\\kappa`` | ``f_c`` teórico | ``f`` donde el gigante baja de 1% |
    |---|---|---|---|
    $(join(filas, "\n"))
    """)
end

#%% hide
md"""
La teoría acierta en Erdős-Rényi. En Barabási-Albert, ``f_c`` es cercano a 1
y crece (lentamente, como ``\ln N``) con el tamaño de la red; cerca del umbral
el componente gigante es tan pequeño que la simulación lo pierde antes.
"""

#%% hide
md"""
## Resumen

| Concepto | Idea |
|---|---|
| Red del mercado (``\rho_{ij} > u``) | Agrupamiento alto por sectores; núcleo financiero que liga a los demás |
| Robustez ante fallas vs. ataques | Libres de escala: aguantan fallas al azar, caen ante ataques a los concentradores |
| Criterio de Molloy-Reed | Hay componente gigante si ``\kappa = \langle k^2\rangle/\langle k\rangle > 2``; umbral ``f_c = 1 - 1/(\kappa - 1)`` |
| Laplaciano y Fiedler | Número de ceros = componentes conexas; vector de Fiedler separa grupos |
| Modelo de configuración | Red aleatoria con los mismos grados; separa efecto del grado del efecto sectorial |
"""

#%% hide
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `red_umbral(C, u)`, que reciba una matriz de correlación `C` de
``N \times N`` y un umbral `u`, y devuelva la matriz de adyacencia booleana o de
enteros `0/1` con ``A_{ij} = 1`` si ``C_{ij} > u`` e ``i \ne j``, y `0` en la
diagonal.
"""

#%%
function red_umbral(C, u)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> red_umbral([1.0 0.7 0.2; 0.7 1.0 0.6; 0.2 0.6 1.0], 0.5),
          [0 1 0; 1 0 1; 0 1 0];
          pista = md"Usa `(C .> u) .& .!I(size(C, 1))`.")

#%% hide
md"""
**2.** Escribe `umbral_fallas(k)`, que reciba el vector de grados `k` de una
red y devuelva el umbral teórico de fallas aleatorias
``f_c = 1 - \frac{1}{\kappa - 1}``, donde ``\kappa = \frac{\langle k^2\rangle}{\langle k\rangle}``.
"""

#%%
function umbral_fallas(k)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> umbral_fallas([2, 4, 4, 6]),
          1.0 - 1.0 / (mean([4, 16, 16, 36]) / mean([2, 4, 4, 6]) - 1.0);
          pista = md"Calcula `κ = mean(k .^ 2) / mean(k)` y devuelve `1 - 1 / (κ - 1)`.")

#%% hide
md"""
### Ruta básica: interpretación

3. En la red del mercado, ¿qué representa un enlace? ¿Qué significa que una
   acción tenga grado alto?
4. Las cadenas de suministro tienen algunos proveedores de los que dependen
   muchas empresas (por ejemplo, fabricantes de semiconductores). Según la
   sección 2, ¿ante qué tipo de choque son vulnerables?
5. ¿Por qué no usamos el árbol de expansión mínima de la sesión 7 para
   comparar el agrupamiento con los modelos?

### Ruta avanzada

6. Repite la tabla de la sección 3.3 con redes de Barabási-Albert de 1,000,
   10,000 y 100,000 nodos. ¿Cómo cambia ``\kappa``? ¿Y ``f_c``?
7. Modifica `curva_robustez` para un ataque **adaptativo**: después de
   cada eliminación, recalcula los grados y quita el de mayor grado en la
   red que queda. ¿Es más dañino que el ataque de la sección 2?

## Para la próxima sesión

- Lectura: Slanina (2014), capítulo 6. Para la ruta avanzada, Newman (2018),
  capítulos 12, 13 y 15; de consulta, Jackson (2008), capítulos 2 a 5.
- En la **sesión 12** (`Code/s12_contagio.jl`) la red deja de ser el objeto de
  estudio y se vuelve el **medio** por donde se propaga algo: cascadas de
  quiebras y contagio de pérdidas en una red interbancaria.
"""

#%% hide
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

#%%
begin
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, CSV, DataFrames, Statistics, LinearAlgebra, Random, Graphs, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

#%%
begin
    const AZUL       = "#2a78d6"   # serie principal
    const NARANJA    = "#eb6834"   # segunda serie o teoría
    const VERDE      = "#2e9e6b"   # tercera serie
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const BLANCO     = "#ffffff"
end;

#%%
"Datos ordenados de mayor a menor y su CCDF empírica, P(X ≥ x)."
function ccdf_empirica(x)
    s = sort(x; rev = true)
    s, (1:length(s)) ./ length(s)
end

#%%
"Fracción de los nodos que está en el componente conexo más grande."
componente_gigante(g) = maximum(length.(connected_components(g))) / nv(g)

#%%
"Distancia promedio entre pares de nodos del componente gigante."
function distancia_promedio(g)
    cc = connected_components(g)
    sub, _ = induced_subgraph(g, cc[argmax(length.(cc))])
    n = nv(sub)
    n < 2 && return 0.0
    sum(sum(gdistances(sub, v)) for v in vertices(sub)) / (n * (n - 1))
end

#%%
"""
Posiciones 2D para dibujar una red: escalamiento multidimensional clásico de
las distancias en la red, refinado con mayorización de estrés.
"""
function disposicion(g; iter = 200)
    n = nv(g)
    δ = Float64.(reduce(hcat, [gdistances(g, s) for s in 1:n]))
    lejos = maximum(δ[δ .< typemax(Int) / 2]) + 1
    δ[δ .> typemax(Int) / 2] .= lejos
    Jc = I - fill(1 / n, n, n)
    F = eigen(Symmetric(-Jc * (δ .^ 2) * Jc / 2))
    X = F.vectors[:, end-1:end] .* sqrt.(max.(F.values[end-1:end], 0))'
    W = [i == j ? 0.0 : 1 / δ[i, j]^2 for i in 1:n, j in 1:n]
    Vp = pinv(Diagonal(vec(sum(W; dims = 2))) - W)
    for _ in 1:iter
        B = [i == j ? 0.0 : -W[i, j] * δ[i, j] / max(norm(X[i, :] - X[j, :]), 1e-9)
             for i in 1:n, j in 1:n]
        B[diagind(B)] .= -vec(sum(B; dims = 2))
        X = Vp * B * X
    end
    X
end

#%%
"Posiciones de las acciones a partir de la distancia de correlación (escalamiento multidimensional clásico)."
posiciones_mercado = let
    D2 = 2 .* (1 .- C)
    Jc = I - fill(1 / N, N, N)
    F = eigen(Symmetric(-Jc * D2 * Jc / 2))
    F.vectors[:, end-1:end] .* sqrt.(F.values[end-1:end])'
end;

#%%
"Dibuja una red con posiciones `X` (n × 2). Los enlaces en `resaltar` van en naranja."
function dibujar_red(g, X; tam = 4, colores = AZUL, alpha_enlaces = 0.5, resaltar = Tuple{Int,Int}[],
                     kwargs...)
    xs, ys = Float64[], Float64[]
    for e in edges(g)
        (src(e), dst(e)) in resaltar && continue
        append!(xs, [X[src(e), 1], X[dst(e), 1], NaN])
        append!(ys, [X[src(e), 2], X[dst(e), 2], NaN])
    end
    p = plot(xs, ys; color = GRIS, alpha = alpha_enlaces, lw = 1, label = "", axis = false,
             grid = false, ticks = false, aspect_ratio = :equal, kwargs...)
    for (a, b) in resaltar
        plot!(p, X[[a, b], 1], X[[a, b], 2]; color = NARANJA, lw = 1.5, label = "")
    end
    scatter!(p, X[:, 1], X[:, 2]; ms = tam, color = colores, msw = 0.5, msc = BLANCO, label = "")
end

#%%
"Ejemplo visual inicial: quitar 6 nodos al azar vs. los 6 mayores hubs en una red Barabási-Albert."
ejemplo_ataque_inicial = let
    g0 = barabasi_albert(45, 2; seed = 2)
    X0 = disposicion(g0)
    n0 = nv(g0)

    quitar_azar = randperm(Xoshiro(11), n0)[1:6]
    keep_azar = setdiff(1:n0, quitar_azar)
    g_azar, _ = induced_subgraph(g0, keep_azar)
    cc_a = connected_components(g_azar)
    gig_a = cc_a[argmax(length.(cc_a))]
    col_azar = [v in gig_a ? AZUL : GRIS_CLARO for v in 1:nv(g_azar)]

    quitar_hubs = sortperm(degree(g0); rev = true)[1:6]
    keep_ataque = setdiff(1:n0, quitar_hubs)
    g_ataque, _ = induced_subgraph(g0, keep_ataque)
    cc_t = connected_components(g_ataque)
    gig_t = cc_t[argmax(length.(cc_t))]
    col_ataque = [v in gig_t ? NARANJA : GRIS_CLARO for v in 1:nv(g_ataque)]

    (g_azar = g_azar, X_azar = X0[keep_azar, :], col_azar = col_azar,
     gig_azar = round(Int, 100 * length(gig_a) / n0),
     g_ataque = g_ataque, X_ataque = X0[keep_ataque, :], col_ataque = col_ataque,
     gig_ataque = round(Int, 100 * length(gig_t) / n0), piezas_ataque = length(cc_t))
end;

#%%
"""
Fracción de los nodos originales en el componente gigante al quitar los
nodos en el `orden` dado.
"""
function curva_robustez(g, orden)
    n = nv(g)
    jefe, tam = collect(1:n), ones(Int, n)
    activo = falses(n)
    raiz(i) = jefe[i] == i ? i : (jefe[i] = raiz(jefe[i]))
    S = zeros(n + 1)
    mayor = 0
    for (t, v) in enumerate(reverse(orden))
        activo[v] = true
        for u in neighbors(g, v)
            activo[u] || continue
            a, b = raiz(u), raiz(v)
            a == b && continue
            tam[a] < tam[b] && ((a, b) = (b, a))
            jefe[b] = a
            tam[a] += tam[b]
        end
        mayor = max(mayor, tam[raiz(v)])
        S[n-t+1] = mayor / n
    end
    S
end

#%%
"Recablea `veces` pares de enlaces (a-b, c-d) → (a-d, c-b) sin cambiar los grados."
function recablear(g, veces; rng = Xoshiro(3))
    h = copy(g)
    E = [(src(e), dst(e)) for e in edges(h)]
    hechos = 0
    while hechos < veces
        i, j = rand(rng, 1:length(E), 2)
        i == j && continue
        (a, b), (c, d) = E[i], E[j]
        rand(rng, Bool) && ((c, d) = (d, c))
        (a == d || c == b || has_edge(h, a, d) || has_edge(h, c, b)) && continue
        rem_edge!(h, a, b); rem_edge!(h, c, d)
        add_edge!(h, a, d); add_edge!(h, c, b)
        E[i], E[j] = (a, d), (c, b)
        hechos += 1
    end
    h
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
