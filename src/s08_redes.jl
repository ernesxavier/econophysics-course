#%% hide
md"""
# Sesión 8: Redes: estructura y modelos

**Objetivos de la sesión**
1. Representar una red con su **matriz de adyacencia** y calcular con ella
   grados, caminos y triángulos.
2. Medir la importancia de un nodo con distintas **centralidades** y la
   estructura local con el **agrupamiento**.
3. Conocer tres modelos de redes aleatorias: **Erdős-Rényi**,
   **Watts-Strogatz** y **Barabási-Albert**, y qué propiedad de las redes
   reales explica cada uno.
4. Comparar los modelos con una red real: la red de correlaciones del
   S&P 100 de la sesión 5.
5. Medir la **robustez**: qué pasa cuando fallan nodos al azar y cuando se
   ataca a los más conectados.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s08_redes.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): L.B4 (matriz de
adyacencia), F.B4 (transiciones de fase), P.I5 (leyes de potencias), L.I4
(potencias de matrices y centralidad). Ruta avanzada: L.A1 (Laplaciano),
L.A5 (Perron-Frobenius), F.A5 (percolación).
"""

#%% hide
TableOfContents(title = "Contenido")

#%% hide
md"""
## 0. El fenómeno: ¿por qué ganaron los Médici?

En 1434, los **Médici** tomaron el control de Florencia y, con algunas
interrupciones, lo conservaron por tres siglos. No eran la familia más rica ni la que tenía más
asientos en el gobierno: los **Strozzi** los superaban en ambas cosas.
Padgett y Ansell (1993) propusieron que la diferencia estaba en la **red**:
con quién se casaba cada familia.

Abajo está la red de matrimonios entre 15 familias de la élite florentina.
Cada línea es al menos un matrimonio entre las dos familias.
"""

#%% hide
let
    X = disposicion(red_medici)
    tam = 4 .+ 3 .* degree(red_medici)
    colores = [f == "Médici" ? NARANJA : f == "Strozzi" ? AZUL : GRIS_CLARO for f in familias]
    p = dibujar_red(red_medici, X; tam, colores, size = (620, 460),
                    title = "Matrimonios entre familias de Florencia, c. 1430")
    annotate!(p, [(X[i, 1], X[i, 2] - 0.32, text(familias[i], 8, GRIS)) for i in 1:nv(red_medici)])
    p
end

#%% hide
md"""
Los Médici están en el centro: **conectan** a familias que no tienen otra
manera de comunicarse entre sí. Los Strozzi tienen casi tantos enlaces, pero
todos dentro del mismo grupo. Hoy vamos a convertir esta intuición en
números, y después a preguntarnos qué **mecanismos** generan las redes que
vemos en la economía.
"""

#%% hide
md"""
## 1. Repaso: la matriz de adyacencia

Una red con ``n`` nodos se guarda en una matriz ``A`` de ``n \times n``:
``A_{ij} = 1`` si hay un enlace entre ``i`` y ``j``, y 0 si no (repaso L.B4).
En una red **no dirigida**, como la de matrimonios, ``A`` es simétrica.

La gráfica muestra la matriz de las familias florentinas: cada cuadro
oscuro es un 1.
"""

#%%
A = Matrix(adjacency_matrix(red_medici));

#%% hide
heatmap(familias, familias, A; color = cgrad([BLANCO, AZUL]), yflip = true,
        aspect_ratio = :equal, colorbar = false, xrotation = 60, size = (560, 520),
        title = "Matriz de adyacencia (oscuro = 1)", tickfontsize = 7)

#%% hide
md"""
Casi todo lo que vamos a medir sale de operaciones con ``A``:

- El **grado** ``k_i`` (número de vecinos) es la suma de la fila ``i``:
  ``k = A\,\mathbf{1}``.
- El número de enlaces es ``m = \tfrac{1}{2}\sum_i k_i``.
- ``(A^2)_{ij}`` cuenta los **caminos de longitud 2** de ``i`` a ``j``, es
  decir, los vecinos que ``i`` y ``j`` tienen en común.
- ``(A^3)_{ii}`` cuenta los caminos de longitud 3 que salen de ``i`` y
  regresan a ``i``: cada **triángulo** se cuenta 6 veces en total (3 vértices
  por 2 sentidos), así que hay ``\mathrm{tr}(A^3)/6`` triángulos (repaso L.I4).
"""

#%%
k = A * ones(Int, size(A, 1))   # grados

#%%
sum(k) ÷ 2, tr(A^3) ÷ 6        # enlaces y triángulos

#%% hide
md"""
Hay $(sum(k) ÷ 2) enlaces y sólo $(tr(A^3) ÷ 6) triángulos. Los Médici tienen
grado $(k[findfirst(==("Médici"), familias)]), el mayor de la red, pero los
Strozzi y los Guadagni tienen 4: el grado por sí solo no explica la
diferencia.
"""

#%% hide
md"""
## 2. Medidas: centralidad y agrupamiento

Hay varias maneras de decir que un nodo es "importante". Cada una responde
a una pregunta distinta.

| Medida | Pregunta | Definición |
|---|---|---|
| **Grado** | ¿Con cuántos se relaciona? | ``k_i = \sum_j A_{ij}`` |
| **Cercanía** | ¿Qué tan rápido llega a todos? | ``1/`` (distancia promedio de ``i`` a los demás) |
| **Intermediación** | ¿Por cuántos caminos cortos pasa? | Fracción de los caminos más cortos entre otros dos nodos que pasan por ``i`` |
| **Eigenvector** | ¿Está conectado con nodos importantes? | ``x_i \propto \sum_j A_{ij} x_j``, el primer eigenvector de ``A`` |

La **distancia** entre dos nodos es el número de enlaces del camino más
corto que los une. La intermediación es la medida del "intermediario": si
la familia ``i`` está en el camino entre otras dos, puede controlar lo que
fluye entre ellas (información, crédito, favores políticos).

El **coeficiente de agrupamiento** de un nodo mide si sus vecinos se
conocen entre sí:

```math
C_i = \frac{\text{enlaces entre los vecinos de } i}{k_i(k_i - 1)/2} .
```

Vale 1 si todos los vecinos están conectados entre sí y 0 si ninguno lo
está. En redes sociales suele ser alto: los amigos de mis amigos tienden a
ser mis amigos.
"""

#%%
medidas = DataFrame(familia = familias,
                    grado = degree(red_medici),
                    cercania = closeness_centrality(red_medici),
                    intermediacion = betweenness_centrality(red_medici),
                    eigenvector = eigenvector_centrality(red_medici),
                    agrupamiento = local_clustering_coefficient(red_medici));

#%% hide
let
    d = sort(medidas, :intermediacion; rev = true)
    redondear(x) = x isa AbstractFloat ? round(x, digits = 2) : x
    d = mapcols(c -> redondear.(c), d)
    rename!(d, ["Familia", "Grado", "Cercanía", "Intermediación", "Eigenvector", "Agrupamiento"])
end

#%% hide
let
    d = sort(medidas, :intermediacion; rev = true)
    f(c) = d[!, c] ./ maximum(d[!, c])
    x = 1:nrow(d)
    p = plot(; xticks = (x, d.familia), xrotation = 45, ylabel = "relativa al máximo",
             legend = :topright, size = (760, 340), bottom_margin = 8Plots.mm,
             title = "Cuatro centralidades, ordenadas por intermediación")
    for (c, color, nombre, m) in ((:grado, GRIS, "Grado", :square), (:cercania, RAMPA[2], "Cercanía", :diamond),
                                  (:eigenvector, RAMPA[4], "Eigenvector", :utriangle),
                                  (:intermediacion, NARANJA, "Intermediación", :circle))
        plot!(p, x, f(c); color, lw = 1.5, marker = m, ms = 4, msw = 0, label = nombre)
    end
    p
end

#%% hide
md"""
Las cuatro medidas ponen a los Médici en primer lugar, pero la
**intermediación** es la que los separa del resto: es el doble que la de los
Guadagni, la segunda familia. Más de la mitad de los caminos más cortos entre
otras familias pasa por ellos. Padgett y Ansell atribuyen a esta posición
buena parte del ascenso de los Médici: eran el único puente entre grupos
que no se hablaban, y sus aliados dependían de ellos, no unos de otros.

El agrupamiento de los Médici es casi cero: sus vecinos **no** están
conectados entre sí. Ser el centro de una estrella, y no parte de un grupo
cerrado, es precisamente lo que les daba poder.
"""

#%% hide
md"""
## 3. Tres modelos de redes aleatorias

Para entender una red real hay que compararla con algo. Un **modelo nulo**
es una red aleatoria que comparte algunas propiedades con la real (número de
nodos, número de enlaces) y en lo demás es puro azar. Si la red real se
parece al modelo, la propiedad se explica por azar; si no, hay un mecanismo
que buscar.

### 3.1 Erdős-Rényi: enlaces al azar

**Reglas.** Hay ``N`` nodos. Cada uno de los ``N(N-1)/2`` pares se conecta
con probabilidad ``p``, independientemente de los demás. El grado promedio
es ``c = p\,(N-1)``.

El grado de cada nodo es la suma de ``N - 1`` monedas con probabilidad ``p``:
una binomial, que para ``N`` grande y ``c`` fijo es una **Poisson** con media
``c`` (repaso P.B1). Todos los nodos tienen más o menos el mismo grado.

Mueve el grado promedio y observa el **componente gigante** (en naranja): el
grupo conectado más grande.

Grado promedio ``c``: $(@bind c_er Slider([0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0], default = 1.2, show_value = true))
"""

#%% hide
let
    g = erdos_renyi(300, c_er / 299; seed = 4)
    X = disposicion(g)
    cc = connected_components(g)
    gigante = cc[argmax(length.(cc))]
    colores = [v in gigante ? NARANJA : GRIS_CLARO for v in vertices(g)]
    p1 = dibujar_red(g, X; tam = 3, colores, size = (460, 420),
                     title = "Erdős-Rényi, N = 300, c = $c_er")
    S = [mean(componente_gigante(erdos_renyi(2000, c / 1999; seed = s)) for s in 1:5) for c in grid_c]
    p2 = plot(grid_c, S_teoria.(grid_c); color = NARANJA, lw = 2, label = "Teoría, N → ∞",
              xlabel = "grado promedio c", ylabel = "fracción en el componente gigante",
              legend = :bottomright, title = "Aparece en c = 1")
    scatter!(p2, grid_c, S; color = AZUL, ms = 3, msw = 0, label = "Simulación, N = 2,000")
    vline!(p2, [c_er]; color = GRIS, ls = :dash, label = "")
    plot(p1, p2; layout = (1, 2), size = (920, 400))
end

#%% hide
md"""
Con ``c < 1`` la red está hecha de islas pequeñas. En ``c = 1``, en promedio
**un vecino por nodo**, aparece de golpe un componente gigante que contiene
una fracción finita de todos los nodos. Es una **transición de fase**
(repaso F.B4), como la condensación de la sesión 7: la fracción en el
componente gigante es el parámetro de orden y ``c`` el parámetro de control.

**Lo que no explica.** En una red de Erdős-Rényi grande, la probabilidad de
que dos vecinos de un nodo estén conectados es la misma ``p`` que para
cualquier par: el agrupamiento es ``C = p = c/(N-1)``, casi cero. Las redes
sociales y económicas tienen agrupamiento alto.
"""

#%% hide
md"""
### 3.2 Watts-Strogatz: mundo pequeño

**Reglas.** Los ``N`` nodos se acomodan en un círculo y cada uno se conecta
con sus ``K`` vecinos más cercanos (``K/2`` de cada lado): mucho
agrupamiento, como un barrio. Luego, cada enlace se **recablea** con
probabilidad ``p``: uno de sus extremos se mueve a un nodo al azar.

Probabilidad de recableado ``p``: $(@bind p_ws Slider([0.0, 0.02, 0.05, 0.1, 0.3, 1.0], default = 0.05, show_value = true))
"""

#%% hide
let
    n = 40
    g = watts_strogatz(n, 4, p_ws; seed = 2)
    θ = 2π .* (0:n-1) ./ n
    X = [cos.(θ) sin.(θ)]
    largos = [(src(e), dst(e)) for e in edges(g) if min(abs(src(e) - dst(e)), n - abs(src(e) - dst(e))) > 2]
    p1 = dibujar_red(g, X; tam = 4, colores = fill(AZUL, n), size = (420, 420),
                     title = "N = 40, K = 4, p = $p_ws", resaltar = largos)
    p2 = plot(barrido_ws.p, barrido_ws.C ./ barrido_ws.C[1]; xscale = :log10, color = AZUL, lw = 2,
              marker = :circle, ms = 3, msw = 0, label = "Agrupamiento C(p)/C(0)",
              xlabel = "probabilidad de recableado p", ylims = (0, 1.05), legend = :bottomleft,
              title = "N = 1,000, K = 10")
    plot!(p2, barrido_ws.p, barrido_ws.L ./ barrido_ws.L[1]; color = NARANJA, lw = 2,
          marker = :square, ms = 3, msw = 0, label = "Distancia promedio L(p)/L(0)")
    p_ws > 0 && vline!(p2, [p_ws]; color = GRIS, ls = :dash, label = "")
    plot(p1, p2; layout = (1, 2), size = (920, 400))
end

#%% hide
md"""
Los enlaces recableados (en naranja) son **atajos**. Basta con que el 1% de
los enlaces sea un atajo para que la distancia promedio caiga a una fracción
de su valor original, mientras el agrupamiento casi no cambia. Hay un rango
amplio de ``p``, entre 0.01 y 0.1, donde la red tiene **las dos cosas**:
grupos cerrados como un barrio y distancias cortas como una red al azar.
Watts y Strogatz (1998) lo llamaron **mundo pequeño**, por los "seis grados
de separación" de Milgram.

**Lo que no explica.** Todos los nodos siguen teniendo un grado cercano a
``K``. No hay nodos con cientos de enlaces.
"""

#%% hide
md"""
### 3.3 Barabási-Albert: conexión preferencial

**Reglas.** La red **crece**. Cada nuevo nodo llega con ``m`` enlaces y se
conecta a nodos existentes con probabilidad **proporcional a su grado**:

```math
\Pi(i) = \frac{k_i}{\sum_j k_j} .
```

Los que ya tienen muchos enlaces atraen más: "el rico se hace más rico".
Es el crecimiento proporcional de la sesión 4 (Gibrat), ahora aplicado a
enlaces: el grado esperado de un nodo crece en proporción a su grado. En
economía el mismo mecanismo se conoce como proceso de Yule o de Simon. Barabási y Albert (1999) mostraron que
genera una distribución de grados con **cola de Pareto**:

```math
P(K \ge k) \approx \frac{m(m+1)}{k(k+1)} \sim k^{-2} .
```

Comparemos 10,000 nodos de Barabási-Albert con una red de Erdős-Rényi con el
mismo grado promedio, ``\langle k\rangle = 4``.
"""

#%%
ba = barabasi_albert(10_000, 2; seed = 1);

#%%
er = erdos_renyi(10_000, ne(ba); seed = 1);

#%% hide
let
    p = plot(; xscale = :log10, yscale = :log10, xlabel = "grado k", ylabel = "P(K ≥ k)",
             legend = :topright, ylims = (5e-5, 1.5), title = "Distribución de grados, ⟨k⟩ = 4")
    for (g, color, nombre) in ((er, NARANJA, "Erdős-Rényi"), (ba, AZUL, "Barabási-Albert"))
        s, P = ccdf_grados(degree(g))
        scatter!(p, s, P; color, ms = 3, msw = 0, label = "$nombre: grado máximo $(maximum(s))")
    end
    x = exp10.(range(log10(2), log10(maximum(degree(ba))), length = 100))
    plot!(p, x, 6 ./ (x .* (x .+ 1)); color = GRIS, ls = :dash, lw = 1.5,
          label = "m(m+1)/(k(k+1)), pendiente -2")
end

#%% hide
md"""
En Erdős-Rényi nadie tiene más de $(maximum(degree(er))) enlaces: la cola de la
Poisson cae más rápido que cualquier potencia. En Barabási-Albert hay
**concentradores** (*hubs*) con cientos de enlaces, y la CCDF es una recta
de pendiente ``-2`` en escala log-log (repaso P.I5). A una red con esta cola
se le llama **libre de escala**.

**Lo que no explica.** El agrupamiento de Barabási-Albert también tiende a
cero al crecer ``N``: $(round(global_clustering_coefficient(ba), digits = 3))
en esta red. Ningún modelo simple tiene todo; cada uno aísla un mecanismo.

| Modelo | Mecanismo | Distribución de grados | Agrupamiento | Distancia típica |
|---|---|---|---|---|
| Erdős-Rényi | Enlaces al azar | Poisson | ``c/N``, casi 0 | Corta: ``\ln N/\ln c`` |
| Watts-Strogatz | Barrio + atajos | Concentrada en ``K`` | Alto | Corta si ``p \gtrsim 0.01`` |
| Barabási-Albert | Crecimiento + conexión preferencial | Pareto, ``P(K \ge k) \sim k^{-2}`` | Casi 0 | Muy corta |
"""

#%% hide
md"""
## 4. Laboratorio: la red del mercado

En la sesión 5 construimos el árbol de expansión mínima de 95 acciones del
S&P 100. Un árbol no tiene triángulos por construcción, así que su
agrupamiento es cero: no sirve para compararlo con los modelos. Usemos otra
red, más simple: conectamos dos acciones si su correlación supera un
**umbral** ``u``,

```math
A_{ij} = 1 \quad\text{si}\quad \rho_{ij} > u .
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
sesión 5: acciones cercanas en la gráfica se mueven juntas. El tamaño de
cada punto es su grado; las cinco acciones con más enlaces tienen su
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
## 5. Laboratorio: robustez ante fallas y ataques

¿Qué tan frágil es una red? Quitemos nodos y midamos qué fracción de los
``N`` nodos originales sigue en el componente gigante. Dos escenarios:

- **Fallas aleatorias**: los nodos se quitan en orden al azar (una empresa
  quiebra por razones propias).
- **Ataque dirigido**: se quitan primero los de mayor grado (un regulador o
  un choque que pega justo a los más conectados).

Usamos redes de Erdős-Rényi y de Barabási-Albert con 10,000 nodos y
``\langle k\rangle = 4``, las de la sección 3.3.
"""

#%%
robustez = Dict(
    (nombre, escenario) => curva_robustez(g, escenario == "Fallas" ? randperm(Xoshiro(3), nv(g)) :
                                              sortperm(degree(g); rev = true))
    for (nombre, g) in (("Erdős-Rényi", er), ("Barabási-Albert", ba)), escenario in ("Fallas", "Ataque"));

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
  (sección 6.3).
- Ante un **ataque**, la diferencia es enorme. La red de Barabási-Albert se
  rompe con apenas el 15% de sus nodos: al quitar los concentradores
  desaparece justo lo que la mantenía conectada. La de Erdős-Rényi resiste
  hasta cerca del 37%, porque no tiene nodos especiales.

Albert, Jeong y Barabási (2000) mostraron este contraste en redes reales,
como internet: las redes libres de escala son **robustas ante fallas pero
frágiles ante ataques**. En la semana 9 veremos la versión
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
## 6. [Avanzado] Espectro, modelo de configuración y percolación

### 6.1 El Laplaciano

El **Laplaciano** de la red es ``L = D - A``, con ``D`` la matriz diagonal
de grados (repaso L.A1). Sus eigenvalores son reales y no negativos, y:

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
sesión 5, encontrada ahora sólo con la red.

### 6.2 Modelo de configuración y asortatividad

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

### 6.3 El umbral de percolación

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
| Matriz de adyacencia ``A`` | Grados ``k = A\mathbf 1``; triángulos ``\mathrm{tr}(A^3)/6`` |
| Centralidades | Grado, cercanía, intermediación, eigenvector: cada una mide un tipo de importancia |
| Agrupamiento | ¿Mis vecinos se conocen? Alto en redes sociales y en el mercado |
| Erdős-Rényi | Grados Poisson; componente gigante en ``c = 1``; agrupamiento casi 0 |
| Watts-Strogatz | Pocos atajos bastan para distancias cortas con agrupamiento alto |
| Barabási-Albert | Conexión preferencial → ``P(K \ge k) \sim k^{-2}``, concentradores |
| Robustez | Libres de escala: aguantan fallas al azar, caen ante ataques a los concentradores |
| Red del mercado | Agrupamiento alto por sectores; núcleo financiero que liga a los demás |
"""

#%% hide
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `grados(A)`, que devuelva el vector de grados de una red a
partir de su matriz de adyacencia.
"""

#%%
function grados(A)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> grados([0 1 1 0; 1 0 1 0; 1 1 0 1; 0 0 1 0]), [2, 2, 3, 1];
          pista = md"El grado es la suma de cada fila: `vec(sum(A; dims = 2))`, o bien `A * ones(Int, size(A, 1))`.")

#%% hide
md"""
**2.** Escribe `agrupamiento(A, i)`, el coeficiente de agrupamiento del nodo
``i``. Pista: ``(A^3)_{ii}`` es el doble del número de enlaces entre los
vecinos de ``i``. Si el grado es menor que 2, devuelve `0.0`.
"""

#%%
function agrupamiento(A, i)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> [agrupamiento([0 1 1 0; 1 0 1 0; 1 1 0 1; 0 0 1 0], i) for i in 1:4],
          [1.0, 1.0, 1 / 3, 0.0];
          pista = md"Con `k = sum(A[i, :])`, el resultado es `(A^3)[i, i] / (k * (k - 1))` si `k ≥ 2`.")

#%% hide
md"""
### Ruta básica: interpretación

3. Da un ejemplo económico de una red donde la **intermediación** importe
   más que el grado, y explica por qué.
4. ¿Por qué el agrupamiento de una red de Erdős-Rényi grande es casi cero?
   ¿Qué propiedad de las redes sociales viola ese modelo?
5. En la red del mercado, ¿qué representa un enlace? ¿Qué significa que una
   acción tenga grado alto?
6. Las cadenas de suministro tienen algunos proveedores de los que dependen
   muchas empresas (por ejemplo, fabricantes de semiconductores). Según la
   sección 5, ¿ante qué tipo de choque son vulnerables?
7. ¿Por qué no usamos el árbol de expansión mínima de la sesión 5 para
   comparar el agrupamiento con los modelos?

### Ruta avanzada

8. Demuestra que en Erdős-Rényi el grado tiende a una Poisson con media
   ``c`` cuando ``N \to \infty`` con ``c = p(N-1)`` fijo.
9. Calcula el vector de Fiedler de la red de los Médici. ¿Qué dos grupos
   separa? ¿Dónde quedan los Médici y los Strozzi?
10. Repite la tabla de la sección 6.3 con redes de Barabási-Albert de 1,000,
    10,000 y 100,000 nodos. ¿Cómo cambia ``\kappa``? ¿Y ``f_c``?
11. Modifica `curva_robustez` para un ataque **adaptativo**: después de
    cada eliminación, recalcula los grados y quita el de mayor grado en la
    red que queda. ¿Es más dañino que el ataque de la sección 5?

## Para la próxima sesión

- Lectura: Newman (2018), capítulos 6 y 7 (medidas) y 10 y 11 (redes
  aleatorias); para la ruta avanzada, capítulos 12, 13 y 15.
- En la semana 9 la red deja de ser el objeto de estudio y se vuelve el
  **medio** por donde se propaga algo: cascadas de quiebras en una red
  interbancaria.
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
    # Usa el entorno del curso (Code/Project.toml) en lugar del gestor de
    # paquetes propio de Pluto: un solo entorno para todo el repositorio.
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
    const RAMPA      = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # azul secuencial
end;

#%%
# Red de matrimonios entre familias florentinas (Padgett y Ansell, 1993), sin
# los Pucci, que no tienen enlaces.
familias, red_medici = let
    nombres = ["Acciaiuoli", "Albizzi", "Barbadori", "Bischeri", "Castellani", "Ginori",
               "Guadagni", "Lamberteschi", "Médici", "Pazzi", "Peruzzi", "Ridolfi",
               "Salviati", "Strozzi", "Tornabuoni"]
    enlaces = [("Médici", "Acciaiuoli"), ("Médici", "Albizzi"), ("Médici", "Barbadori"),
               ("Médici", "Ridolfi"), ("Médici", "Salviati"), ("Médici", "Tornabuoni"),
               ("Castellani", "Barbadori"), ("Castellani", "Peruzzi"), ("Castellani", "Strozzi"),
               ("Strozzi", "Peruzzi"), ("Strozzi", "Ridolfi"), ("Strozzi", "Bischeri"),
               ("Peruzzi", "Bischeri"), ("Guadagni", "Bischeri"), ("Guadagni", "Albizzi"),
               ("Guadagni", "Tornabuoni"), ("Guadagni", "Lamberteschi"), ("Albizzi", "Ginori"),
               ("Ridolfi", "Tornabuoni"), ("Salviati", "Pazzi")]
    g = SimpleGraph(length(nombres))
    for (a, b) in enlaces
        add_edge!(g, findfirst(==(a), nombres), findfirst(==(b), nombres))
    end
    nombres, g
end;

#%%
"Datos ordenados de mayor a menor y su CCDF empírica, P(X ≥ x)."
function ccdf_empirica(x)
    s = sort(x; rev = true)
    s, (1:length(s)) ./ length(s)
end

#%%
"Grados distintos (mayores que cero) y la fracción de nodos con grado mayor o igual, P(K ≥ k)."
function ccdf_grados(k)
    s = sort(unique(filter(>(0), k)))
    s, [mean(k .>= x) for x in s]
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
"Solución de S = 1 - exp(-cS), la fracción en el componente gigante de Erdős-Rényi."
function S_teoria(c)
    S = 1.0
    for _ in 1:2000
        S = 1 - exp(-c * S)
    end
    S < 1e-6 ? 0.0 : S
end

#%%
grid_c = 0.2:0.2:4.0;

#%%
"Agrupamiento y distancia promedio de Watts-Strogatz (N = 1000, K = 10) en función de p."
barrido_ws = let
    ps = [0.0; exp10.(range(-4, 0, length = 13))]
    res = map(ps) do p
        g = watts_strogatz(1000, 10, p; seed = 1)
        Lm = mean(mean(gdistances(g, s)) * 1000 / 999 for s in 1:20:1000)
        (global_clustering_coefficient(g), Lm)
    end
    (p = [1e-5; ps[2:end]], C = first.(res), L = last.(res))   # p = 0 se dibuja en 10⁻⁵
end;

#%%
"""
Posiciones 2D para dibujar una red: escalamiento multidimensional clásico de
las distancias en la red, refinado con mayorización de estrés. Los pares sin
camino se tratan como si estuvieran a la distancia máxima más uno.
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
"""
Fracción de los nodos originales en el componente gigante al quitar los
nodos en el `orden` dado. El resultado `S[q + 1]` corresponde a `q` nodos
quitados. Se calcula al revés: se agregan los nodos del último al primero y
se unen sus grupos con una estructura de conjuntos disjuntos.
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
