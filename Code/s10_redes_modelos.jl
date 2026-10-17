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

# ╔═╡ 3f168f76-ca82-11f1-b535-53b8acc1e32f
begin
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, DataFrames, Statistics, LinearAlgebra, Random, Graphs, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

# ╔═╡ 3f168e60-ca82-11f1-b2a1-275515b25d6d
md"""
# Sesión 10: Redes: medidas y modelos

**Objetivos de la sesión**
1. Representar una red con su **matriz de adyacencia** y calcular con ella
   grados, caminos y triángulos.
2. Medir la importancia de un nodo con distintas **centralidades** y la
   estructura local con el **agrupamiento**.
3. Conocer tres modelos de redes aleatorias: **Erdős-Rényi**,
   **Watts-Strogatz** y **Barabási-Albert**, y qué propiedad de las redes
   reales explica cada uno.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s10_redes_modelos.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): P.B1 (distribuciones),
L.B4 (matriz de adyacencia), C.B2 (escalas log-log), F.B4 (transiciones de
fase). Para las notas formales: P.I5 (leyes de potencias), L.I4 (potencias
de matrices y centralidad), S.I3 (comparar datos con una distribución). Ruta
avanzada: L.A5 (Perron-Frobenius), F.A5 (percolación).
"""

# ╔═╡ 3f168e9a-ca82-11f1-a205-d5479b2f705b
TableOfContents(title = "Contenido")

# ╔═╡ 3f168ea4-ca82-11f1-b1a6-4fcaa351fdc5
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

# ╔═╡ 3f168eae-ca82-11f1-b67c-010594f09892
md"""
Los Médici están en el centro: **conectan** a familias que no tienen otra
manera de comunicarse entre sí. Los Strozzi tienen casi tantos enlaces, pero
todos dentro del mismo grupo. Hoy vamos a convertir esta intuición en
números, y después a preguntarnos qué **mecanismos** generan las redes que
vemos en la economía.
"""

# ╔═╡ 3f168eb8-ca82-11f1-9a27-5177aab28cf3
md"""
## 1. Repaso: la matriz de adyacencia

Una red con ``n`` nodos se guarda en una matriz ``A`` de ``n \times n``:
``A_{ij} = 1`` si hay un enlace entre ``i`` y ``j``, y 0 si no (repaso L.B4).
En una red **no dirigida**, como la de matrimonios, ``A`` es simétrica.

La gráfica muestra la matriz de las familias florentinas: cada cuadro
oscuro es un 1.
"""

# ╔═╡ 3f168ecc-ca82-11f1-ad98-053ded8a0aae
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

# ╔═╡ 3f168ee0-ca82-11f1-a269-0dcfb46d62d5
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

# ╔═╡ 3f168ef2-ca82-11f1-9898-45cd56b102d3
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

# ╔═╡ 3f168efe-ca82-11f1-8262-03f3ae9376ce
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

# ╔═╡ 3f168f12-ca82-11f1-815e-c9d1b9e36d39
md"""
Con ``c < 1`` la red está hecha de islas pequeñas. En ``c = 1``, en promedio
**un vecino por nodo**, aparece de golpe un componente gigante que contiene
una fracción finita de todos los nodos. Es una **transición de fase**
(repaso F.B4), como la condensación de la sesión 9: la fracción en el
componente gigante es el parámetro de orden y ``c`` el parámetro de control.

**Lo que no explica.** En una red de Erdős-Rényi grande, la probabilidad de
que dos vecinos de un nodo estén conectados es la misma ``p`` que para
cualquier par: el agrupamiento es ``C = p = c/(N-1)``, casi cero. Las redes
sociales y económicas tienen agrupamiento alto.
"""

# ╔═╡ 3f168f12-ca82-11f1-b278-97a8f3ea4964
md"""
### 3.2 Watts-Strogatz: mundo pequeño

**Reglas.** Los ``N`` nodos se acomodan en un círculo y cada uno se conecta
con sus ``K`` vecinos más cercanos (``K/2`` de cada lado): mucho
agrupamiento, como un barrio. Luego, cada enlace se **recablea** con
probabilidad ``p``: uno de sus extremos se mueve a un nodo al azar.

Probabilidad de recableado ``p``: $(@bind p_ws Slider([0.0, 0.02, 0.05, 0.1, 0.3, 1.0], default = 0.05, show_value = true))
"""

# ╔═╡ 3f168f1c-ca82-11f1-8e1b-c54925c08ef1
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

# ╔═╡ 3f168f24-ca82-11f1-a35e-1d3548b26407
md"""
### 3.3 Barabási-Albert: conexión preferencial

**Reglas.** La red **crece**. Cada nuevo nodo llega con ``m`` enlaces y se
conecta a nodos existentes con probabilidad **proporcional a su grado**:

```math
\Pi(i) = \frac{k_i}{\sum_j k_j} .
```

Los que ya tienen muchos enlaces atraen más: "el rico se hace más rico".
Es el crecimiento proporcional de la sesión 6 (Gibrat), ahora aplicado a
enlaces: el grado esperado de un nodo crece en proporción a su grado. En
economía el mismo mecanismo se conoce como proceso de Yule o de Simon. Barabási y Albert (1999) mostraron que
genera una distribución de grados con **cola de Pareto**:

```math
P(K \ge k) \approx \frac{m(m+1)}{k(k+1)} \sim k^{-2} .
```

Comparemos 10,000 nodos de Barabási-Albert con una red de Erdős-Rényi con el
mismo grado promedio, ``\langle k\rangle = 4``.
"""

# ╔═╡ 3f168f30-ca82-11f1-9dbd-f380f9381143
ba = barabasi_albert(10_000, 2; seed = 1);

# ╔═╡ 3f168f30-ca82-11f1-b6f5-a5ce329e4e0e
er = erdos_renyi(10_000, ne(ba); seed = 1);

# ╔═╡ 3f168f3a-ca82-11f1-a0c1-8bcc646864e4
md"""
En Erdős-Rényi nadie tiene más de $(maximum(degree(er))) enlaces: la cola de la
Poisson cae más rápido que cualquier potencia. En Barabási-Albert hay
**concentradores** (*hubs*) con cientos de enlaces, y la CCDF es una recta
de pendiente ``-2`` en escala log-log (repaso C.B2 y P.I5). A una red con esta
cola se le llama **libre de escala**.

**Lo que no explica.** El agrupamiento de Barabási-Albert también tiende a
cero al crecer ``N``: $(round(global_clustering_coefficient(ba), digits = 3))
en esta red. Ningún modelo simple tiene todo; cada uno aísla un mecanismo.

| Modelo | Mecanismo | Distribución de grados | Agrupamiento | Distancia típica |
|---|---|---|---|---|
| Erdős-Rényi | Enlaces al azar | Poisson | ``c/N``, casi 0 | Corta: ``\ln N/\ln c`` |
| Watts-Strogatz | Barrio + atajos | Concentrada en ``K`` | Alto | Corta si ``p \gtrsim 0.01`` |
| Barabási-Albert | Crecimiento + conexión preferencial | Pareto, ``P(K \ge k) \sim k^{-2}`` | Casi 0 | Muy corta |
"""

# ╔═╡ 3f168f44-ca82-11f1-a20c-452131c965df
md"""
## Resumen

| Concepto | Idea |
|---|---|
| Matriz de adyacencia ``A`` | Grados ``k = A\mathbf 1``; triángulos ``\mathrm{tr}(A^3)/6`` |
| Centralidades | Grado, cercanía, intermediación, eigenvector: cada una mide un tipo de importancia |
| Agrupamiento | ¿Mis vecinos se conocen? Alto en redes sociales y económicas |
| Erdős-Rényi | Grados Poisson; componente gigante en ``c = 1``; agrupamiento casi 0 |
| Watts-Strogatz | Pocos atajos bastan para distancias cortas con agrupamiento alto |
| Barabási-Albert | Conexión preferencial → ``P(K \ge k) \sim k^{-2}``, concentradores |
"""

# ╔═╡ 3f168f4e-ca82-11f1-b6de-058fd9cd343b
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `grados(A)`, que devuelva el vector de grados de una red a
partir de su matriz de adyacencia.
"""

# ╔═╡ 3f168f4e-ca82-11f1-a9dd-7f9fb19c2dd3
function grados(A)
    missing  # reemplaza esta línea
end

# ╔═╡ 3f168f58-ca82-11f1-91b0-71335c5a5e6a
md"""
**2.** Escribe `agrupamiento(A, i)`, el coeficiente de agrupamiento del nodo
``i``. Pista: ``(A^3)_{ii}`` es el doble del número de enlaces entre los
vecinos de ``i``. Si el grado es menor que 2, devuelve `0.0`.
"""

# ╔═╡ 3f168f58-ca82-11f1-b994-17d25bbb50e5
function agrupamiento(A, i)
    missing  # reemplaza esta línea
end

# ╔═╡ 3f168f62-ca82-11f1-beb1-3b9d1c028062
md"""
### Ruta básica: interpretación

3. Da un ejemplo económico de una red donde la **intermediación** importe
   más que el grado, y explica por qué.
4. ¿Por qué el agrupamiento de una red de Erdős-Rényi grande es casi cero?
   ¿Qué propiedad de las redes sociales viola ese modelo?
5. ¿Qué ingrediente adicional convierte una red de enlaces aleatorios en una
   red con concentradores (*hubs*) como la de Barabási-Albert?

### Ruta avanzada

6. Demuestra que en Erdős-Rényi el grado tiende a una Poisson con media
   ``c`` cuando ``N \to \infty`` con ``c = p(N-1)`` fijo.
7. Deriva ``C(0) = 3(K-2)/(4(K-1))`` para el anillo de Watts-Strogatz.

## Para la próxima sesión

- Lectura: Newman (2018), capítulo 7. Complementaria: Newman (2018), capítulos
  6 y 11; Requião da Cunha (2021), sección 8.1.
- En la **sesión 11** (`Code/s11_redes_mercado.jl`) pondremos estos tres
  modelos frente a la red de correlaciones del S&P 100 y mediremos qué pasa
  cuando fallan nodos al azar o cuando se ataca a los más conectados.
"""

# ╔═╡ 3f168f6c-ca82-11f1-9762-b1a2e851437f
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 3f168f80-ca82-11f1-a8e6-4d3e9d87ea31
begin
    const AZUL       = "#2a78d6"   # serie principal
    const NARANJA    = "#eb6834"   # segunda serie o teoría
    const VERDE      = "#2e9e6b"   # tercera serie
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const BLANCO     = "#ffffff"
    const RAMPA      = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # azul secuencial
end;

# ╔═╡ 3f168f80-ca82-11f1-8b53-a531693860ee
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

# ╔═╡ 3f168ec0-ca82-11f1-9d29-514c13e72aaa
A = Matrix(adjacency_matrix(red_medici));

# ╔═╡ 3f168ed6-ca82-11f1-8772-afa2996b5d2d
k = A * ones(Int, size(A, 1))   # grados

# ╔═╡ 3f168ed6-ca82-11f1-9174-81bfed4dad63
sum(k) ÷ 2, tr(A^3) ÷ 6        # enlaces y triángulos

# ╔═╡ 3f168ec0-ca82-11f1-a0b9-832aef4cba76
heatmap(familias, familias, A; color = cgrad([BLANCO, AZUL]), yflip = true,
        aspect_ratio = :equal, colorbar = false, xrotation = 60, size = (560, 520),
        title = "Matriz de adyacencia (oscuro = 1)", tickfontsize = 7)

# ╔═╡ 3f168ed6-ca82-11f1-84ea-cd66d600ce9c
md"""
Hay $(sum(k) ÷ 2) enlaces y sólo $(tr(A^3) ÷ 6) triángulos. Los Médici tienen
grado $(k[findfirst(==("Médici"), familias)]), el mayor de la red, pero los
Strozzi y los Guadagni tienen 4: el grado por sí solo no explica la
diferencia.
"""

# ╔═╡ 3f168eea-ca82-11f1-87e5-2b82644e6a30
medidas = DataFrame(familia = familias,
                    grado = degree(red_medici),
                    cercania = closeness_centrality(red_medici),
                    intermediacion = betweenness_centrality(red_medici),
                    eigenvector = eigenvector_centrality(red_medici),
                    agrupamiento = local_clustering_coefficient(red_medici));

# ╔═╡ 3f168eea-ca82-11f1-a8d4-69eea812bf8d
let
    d = sort(medidas, :intermediacion; rev = true)
    redondear(x) = x isa AbstractFloat ? round(x, digits = 2) : x
    d = mapcols(c -> redondear.(c), d)
    rename!(d, ["Familia", "Grado", "Cercanía", "Intermediación", "Eigenvector", "Agrupamiento"])
end

# ╔═╡ 3f168ef2-ca82-11f1-b468-9be276aff310
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

# ╔═╡ 3f168f8a-ca82-11f1-a8b3-11f2c246c5d6
"Grados distintos (mayores que cero) y la fracción de nodos con grado mayor o igual, P(K ≥ k)."
function ccdf_grados(k)
    s = sort(unique(filter(>(0), k)))
    s, [mean(k .>= x) for x in s]
end

# ╔═╡ 3f168f3a-ca82-11f1-9133-99ce75287e12
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

# ╔═╡ 3f168f8a-ca82-11f1-9d23-5b5ee5adf503
"Fracción de los nodos que está en el componente conexo más grande."
componente_gigante(g) = maximum(length.(connected_components(g))) / nv(g)

# ╔═╡ 3f168f96-ca82-11f1-b6d4-971e318ded6d
"Solución de S = 1 - exp(-cS), la fracción en el componente gigante de Erdős-Rényi."
function S_teoria(c)
    S = 1.0
    for _ in 1:2000
        S = 1 - exp(-c * S)
    end
    S < 1e-6 ? 0.0 : S
end

# ╔═╡ 3f168f96-ca82-11f1-8133-837202ef05d7
grid_c = 0.2:0.2:4.0;

# ╔═╡ 3f168f9e-ca82-11f1-bd7e-e30aa4886e4b
"Agrupamiento y distancia promedio de Watts-Strogatz (N = 1000, K = 10) en función de p."
barrido_ws = let
    ps = [0.0; exp10.(range(-4, 0, length = 13))]
    res = map(ps) do p
        g = watts_strogatz(1000, 10, p; seed = 1)
        Lm = mean(mean(gdistances(g, s)) * 1000 / 999 for s in 1:20:1000)
        (global_clustering_coefficient(g), Lm)
    end
    (p = [1e-5; ps[2:end]], C = first.(res), L = last.(res))
end;

# ╔═╡ 3f168f9e-ca82-11f1-a823-911cbe0acd7a
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

# ╔═╡ 3f168fa8-ca82-11f1-b24c-d13687954c93
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

# ╔═╡ 3f168ea4-ca82-11f1-a480-3bbb84dbb72e
let
    X = disposicion(red_medici)
    tam = 4 .+ 3 .* degree(red_medici)
    colores = [f == "Médici" ? NARANJA : f == "Strozzi" ? AZUL : GRIS_CLARO for f in familias]
    p = dibujar_red(red_medici, X; tam, colores, size = (620, 460),
                    title = "Matrimonios entre familias de Florencia, c. 1430")
    annotate!(p, [(X[i, 1], X[i, 2] - 0.32, text(familias[i], 8, GRIS)) for i in 1:nv(red_medici)])
    p
end

# ╔═╡ 3f168f08-ca82-11f1-85d0-dbecb1eb997e
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

# ╔═╡ 3f168f1c-ca82-11f1-a055-31dd3e3cc6a3
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

# ╔═╡ 3f168fb2-ca82-11f1-a577-55735d69e8dc
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

# ╔═╡ 3f168f58-ca82-11f1-82a1-3fe7ec7a7ca7
verificar(() -> grados([0 1 1 0; 1 0 1 0; 1 1 0 1; 0 0 1 0]), [2, 2, 3, 1];
          pista = md"El grado es la suma de cada fila: `vec(sum(A; dims = 2))`, o bien `A * ones(Int, size(A, 1))`.")

# ╔═╡ 3f168f62-ca82-11f1-9ce0-c7aaef4381d7
verificar(() -> [agrupamiento([0 1 1 0; 1 0 1 0; 1 1 0 1; 0 0 1 0], i) for i in 1:4],
          [1.0, 1.0, 1 / 3, 0.0];
          pista = md"Con `k = sum(A[i, :])`, el resultado es `(A^3)[i, i] / (k * (k - 1))` si `k ≥ 2`.")

# ╔═╡ Cell order:
# ╟─3f168e60-ca82-11f1-b2a1-275515b25d6d
# ╟─3f168e9a-ca82-11f1-a205-d5479b2f705b
# ╟─3f168ea4-ca82-11f1-b1a6-4fcaa351fdc5
# ╟─3f168ea4-ca82-11f1-a480-3bbb84dbb72e
# ╟─3f168eae-ca82-11f1-b67c-010594f09892
# ╟─3f168eb8-ca82-11f1-9a27-5177aab28cf3
# ╠═3f168ec0-ca82-11f1-9d29-514c13e72aaa
# ╟─3f168ec0-ca82-11f1-a0b9-832aef4cba76
# ╟─3f168ecc-ca82-11f1-ad98-053ded8a0aae
# ╠═3f168ed6-ca82-11f1-8772-afa2996b5d2d
# ╠═3f168ed6-ca82-11f1-9174-81bfed4dad63
# ╟─3f168ed6-ca82-11f1-84ea-cd66d600ce9c
# ╟─3f168ee0-ca82-11f1-a269-0dcfb46d62d5
# ╠═3f168eea-ca82-11f1-87e5-2b82644e6a30
# ╟─3f168eea-ca82-11f1-a8d4-69eea812bf8d
# ╟─3f168ef2-ca82-11f1-b468-9be276aff310
# ╟─3f168ef2-ca82-11f1-9898-45cd56b102d3
# ╟─3f168efe-ca82-11f1-8262-03f3ae9376ce
# ╟─3f168f08-ca82-11f1-85d0-dbecb1eb997e
# ╟─3f168f12-ca82-11f1-815e-c9d1b9e36d39
# ╟─3f168f12-ca82-11f1-b278-97a8f3ea4964
# ╟─3f168f1c-ca82-11f1-a055-31dd3e3cc6a3
# ╟─3f168f1c-ca82-11f1-8e1b-c54925c08ef1
# ╟─3f168f24-ca82-11f1-a35e-1d3548b26407
# ╠═3f168f30-ca82-11f1-9dbd-f380f9381143
# ╠═3f168f30-ca82-11f1-b6f5-a5ce329e4e0e
# ╟─3f168f3a-ca82-11f1-9133-99ce75287e12
# ╟─3f168f3a-ca82-11f1-a0c1-8bcc646864e4
# ╟─3f168f44-ca82-11f1-a20c-452131c965df
# ╟─3f168f4e-ca82-11f1-b6de-058fd9cd343b
# ╠═3f168f4e-ca82-11f1-a9dd-7f9fb19c2dd3
# ╟─3f168f58-ca82-11f1-82a1-3fe7ec7a7ca7
# ╟─3f168f58-ca82-11f1-91b0-71335c5a5e6a
# ╠═3f168f58-ca82-11f1-b994-17d25bbb50e5
# ╟─3f168f62-ca82-11f1-9ce0-c7aaef4381d7
# ╟─3f168f62-ca82-11f1-beb1-3b9d1c028062
# ╟─3f168f6c-ca82-11f1-9762-b1a2e851437f
# ╠═3f168f76-ca82-11f1-b535-53b8acc1e32f
# ╠═3f168f80-ca82-11f1-a8e6-4d3e9d87ea31
# ╠═3f168f80-ca82-11f1-8b53-a531693860ee
# ╠═3f168f8a-ca82-11f1-a8b3-11f2c246c5d6
# ╠═3f168f8a-ca82-11f1-9d23-5b5ee5adf503
# ╠═3f168f96-ca82-11f1-b6d4-971e318ded6d
# ╠═3f168f96-ca82-11f1-8133-837202ef05d7
# ╠═3f168f9e-ca82-11f1-bd7e-e30aa4886e4b
# ╠═3f168f9e-ca82-11f1-a823-911cbe0acd7a
# ╠═3f168fa8-ca82-11f1-b24c-d13687954c93
# ╠═3f168fb2-ca82-11f1-a577-55735d69e8dc
