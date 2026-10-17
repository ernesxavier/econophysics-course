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

# ╔═╡ 823ee49e-ca80-11f1-acd2-3d414bc95fa4
begin
    # Usa el entorno del curso (Code/Project.toml) en lugar del gestor de
    # paquetes propio de Pluto: un solo entorno para todo el repositorio.
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, CSV, DataFrames, Dates, Statistics, LinearAlgebra, Random, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

# ╔═╡ 823ee34a-ca80-11f1-906d-37eff8d968b3
md"""
# Sesión 7: Correlaciones y matrices aleatorias

**Objetivos de la sesión**
1. Construir la **matriz de correlación** de 95 acciones del S&P 100 y leerla
   a través de sus eigenvalores y eigenvectores.
2. Separar **señal de ruido** comparando el espectro real con el de los datos
   barajados.
3. Identificar el **modo de mercado** y los **sectores** en los eigenvectores.
4. Convertir correlaciones en distancias y construir el **árbol de expansión
   mínima**, el puente hacia las redes del Módulo II.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s07_correlaciones.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): P.B2 (covarianza y
correlación), L.B1–B3 (vectores y matrices), S.B5 (pruebas por permutación),
L.I1–I3 (eigenvalores, matrices simétricas, componentes principales), L.I5
(distancia entre activos). Ruta avanzada: L.A2 (Marchenko-Pastur).
"""

# ╔═╡ 823ee370-ca80-11f1-995e-9b104afdea92
TableOfContents(title = "Contenido")

# ╔═╡ 823ee37c-ca80-11f1-a695-c5b70d702a47
md"""
## 0. El fenómeno: las acciones se mueven juntas

Hasta ahora estudiamos cada serie por separado. Hoy vemos **muchas a la
vez**: los precios diarios de 95 acciones del S&P 100, de enero de 2015 a
septiembre de 2026. Los datos están en `Data/sp100_precios.csv` y
`Data/sp100_sectores.csv`, y se regeneran con `Code/descargar_sp100.jl`.

Una advertencia: usamos las empresas que están **hoy** en el índice. Las que
salieron (por quiebra, compra o pérdida de valor) no aparecen, así que la
muestra tiene **sesgo de supervivencia**.
"""

# ╔═╡ 823ee390-ca80-11f1-aca8-6b8e059ea3a4
begin
    precios = CSV.read(joinpath(@__DIR__, "..", "Data", "sp100_precios.csv"), DataFrame)
    sectores = CSV.read(joinpath(@__DIR__, "..", "Data", "sp100_sectores.csv"), DataFrame)
    ancho = unstack(precios, :fecha, :ticker, :cierre)    # una columna por acción
    tickers = names(ancho)[2:end]
    R = diff(log.(Matrix{Float64}(ancho[:, 2:end])); dims = 1)   # rendimientos: T × N
    fechas = ancho.fecha[2:end]
    T, N = size(R)
end;

# ╔═╡ 823ee390-ca80-11f1-861c-bd9aacb29bff
md"""
Tenemos ``T = `` $(T) días y ``N = `` $(N) acciones. Elige dos acciones y
compara sus rendimientos diarios, día por día:

Acción 1: $(@bind accion1 Select(tickers, default = "JPM"))
  Acción 2: $(@bind accion2 Select(tickers, default = "BAC"))
"""

# ╔═╡ 823ee3a2-ca80-11f1-95aa-537ae946b2c0
md"""
## 1. Repaso: la matriz de correlación

Si ``z_{ti}`` es el rendimiento de la acción ``i`` el día ``t``, estandarizado
(media 0 y desviación estándar 1), la correlación entre ``i`` y ``j`` es

```math
\rho_{ij} = \frac{1}{T-1}\sum_{t=1}^{T} z_{ti}\, z_{tj},
\qquad\text{o en forma matricial}\qquad
C = \frac{1}{T-1}\, Z^\top Z .
```

``C`` es una matriz ``N \times N`` simétrica con unos en la diagonal. Abajo
está completa, con las acciones ordenadas por sector: los bloques oscuros
sobre la diagonal son sectores cuyas acciones se mueven juntas.
"""

# ╔═╡ 823ee3ae-ca80-11f1-a874-ef132643d574
C = cor(R);

# ╔═╡ 823ee3b8-ca80-11f1-a4f2-dbae20b60a02
md"""
**Eigenvalores como direcciones principales** (repaso L.I1–I3). Como ``C`` es
simétrica, tiene ``N`` eigenvalores reales ``\lambda_1 \ge \lambda_2 \ge \dots \ge \lambda_N \ge 0``
con eigenvectores ortogonales ``v_1, \dots, v_N``. Cada eigenvector es una
**cartera**: una combinación de las ``N`` acciones. Su eigenvalor es la
varianza de esa cartera (con rendimientos estandarizados), y los eigenvalores
suman ``N``:

```math
C\, v_k = \lambda_k\, v_k, \qquad \sum_{k=1}^{N} \lambda_k = \operatorname{tr} C = N .
```

Si las acciones fueran independientes, ``C`` sería la identidad y todos los
eigenvalores valdrían 1. Un eigenvalor grande indica una dirección en la que
muchas acciones se mueven juntas.
"""

# ╔═╡ 823ee3c2-ca80-11f1-8a23-2144a189a357
espectro = eigen(Symmetric(C));   # eigenvalores en orden creciente

# ╔═╡ 823ee3c2-ca80-11f1-8a19-bf1625be78ea
λ = reverse(espectro.values);     # λ[1] es el mayor

# ╔═╡ 823ee3cc-ca80-11f1-9b7d-3933bc5ad931
md"""
Los cinco mayores eigenvalores son $(join(round.(λ[1:5], digits = 1), ", ")).
El primero, ``\lambda_1 = `` $(round(λ[1], digits = 1)), concentra el
$(round(Int, 100λ[1] / N))% de la varianza total. Pero ¿cuáles de los demás
son significativos? Con ``T`` finito, incluso acciones independientes darían
eigenvalores distintos de 1.
"""

# ╔═╡ 823ee3cc-ca80-11f1-92f7-291318cf2056
md"""
## 2. Laboratorio: ¿señal o ruido?

La prueba es la de permutaciones de la sesión 2 (repaso S.B5): **barajamos
cada serie por separado**. Cada acción conserva exactamente sus rendimientos
(su media, su varianza y sus colas), pero los días quedan desalineados, así
que las correlaciones entre acciones desaparecen. Lo que sobresale del
espectro barajado es **señal**.

$(@bind barajar Button("Barajar otra vez"))
"""

# ╔═╡ 823ee3d6-ca80-11f1-986f-c577c6ad93d8
λ_barajado = let
    barajar
    rng = Xoshiro()
    Rb = similar(R)
    for j in 1:N
        Rb[:, j] = shuffle(rng, R[:, j])
    end
    sort(eigvals(Symmetric(cor(Rb))); rev = true)
end;

# ╔═╡ 823ee3e0-ca80-11f1-83d2-f1c92b705e6f
md"""
Al barajar, todos los eigenvalores quedan entre
$(round(minimum(λ_barajado), digits = 2)) y $(round(maximum(λ_barajado), digits = 2)):
es el **ruido de estimación** de tener ``T`` finito. En los datos reales,
$(count(>(maximum(λ_barajado)), λ)) eigenvalores quedan por encima de ese
máximo. Esos son la señal; el resto del espectro real se parece al barajado
y no se puede distinguir del ruido. Presiona el botón varias veces: el
espectro barajado casi no cambia.

Nota que el grueso del espectro real queda un poco **a la izquierda** del
barajado: los eigenvalores grandes se llevan parte de la suma ``N``, y al
resto le toca menos.
"""

# ╔═╡ 823ee3ea-ca80-11f1-82e4-43177fc88f50
md"""
## 3. El modo de mercado

El primer eigenvector ``v_1`` asigna a cada acción un peso. La gráfica los
muestra multiplicados por ``\sqrt{N}``, así que una cartera con el mismo peso
en todas las acciones daría barras de altura 1.
"""

# ╔═╡ 823ee3ea-ca80-11f1-81fc-233c2e60b46b
v1 = espectro.vectors[:, end] .* sign(sum(espectro.vectors[:, end]));   # el signo es arbitrario

# ╔═╡ 823ee3f4-ca80-11f1-8bac-89c9e80546f2
md"""
Todas las componentes son positivas y parecidas: ``v_1`` es, casi, comprar
un poco de cada acción. Es el **modo de mercado**: el factor común que mueve
a todas a la vez. Para comprobarlo, comparemos el rendimiento diario de esa
cartera con el del índice S&P 500 completo (de `Data/precios.csv`):
"""

# ╔═╡ 823ee3fe-ca80-11f1-86e9-5d5ebc9bab20
md"""
La correlación es casi 1: el primer eigenvector **es** el mercado, aunque no
le dijimos nada sobre el índice.

**Laboratorio: ¿cambia el modo de mercado con el tiempo?** Calculemos
``\lambda_1/N``, la fracción de la varianza que explica el mercado, en una
ventana móvil.

Longitud de la ventana (días hábiles): $(@bind ventana Select([60 => "60 (3 meses)", 125 => "125 (6 meses)", 250 => "250 (1 año)", 500 => "500 (2 años)"], default = 250))
"""

# ╔═╡ 823ee412-ca80-11f1-bb82-9b9320842c63
md"""
En marzo de 2020 (la caída por la pandemia) el mercado llegó a explicar más
de la mitad de la varianza: **en las crisis todo cae junto**. Es una mala
noticia para la diversificación, que protege justo cuando las correlaciones
son bajas. Con ventanas cortas la serie es más ruidosa, porque cada matriz
se estima con menos días.
"""

# ╔═╡ 823ee412-ca80-11f1-99fb-37271add34f8
md"""
## 4. Sectores en los eigenvectores

Los demás eigenvectores significativos son ortogonales a ``v_1``: describen
cómo se mueven unas acciones **respecto a otras**, una vez descontado el
mercado. Mueve ``k`` para ver los eigenvectores que salieron del ruido. El
signo de cada eigenvector es arbitrario: lo que importa es qué acciones
tienen el mismo signo y cuáles el contrario.

Eigenvector ``k``: $(@bind k_vec Slider(2:8, default = 2, show_value = true))
"""

# ╔═╡ 823ee41c-ca80-11f1-8e2d-45ea67166fa7
md"""
Cada eigenvector cuenta una historia económica:

- ``v_2`` opone las acciones **defensivas** (servicios públicos, consumo
  básico, farmacéuticas, telefónicas) a las de **tecnología**: cuando la
  tecnología sube, las defensivas tienden a quedarse atrás, y al revés.
- ``v_3`` separa a los **bancos** de las empresas de software.
- ``v_4`` separa a los **semiconductores** (AMAT, LRCX, MU, INTC) del
  software, dentro del mismo sector de tecnología.
- ``v_5`` opone las **farmacéuticas** a los servicios públicos y las
  tabacaleras, y ``v_6`` junta a la **energía** con la industria de
  **defensa**, frente a comercio y bancos.

Los sectores oficiales no son la última palabra: los datos encuentran grupos
que cruzan sectores (defensivas) o que los parten (semiconductores contra
software).
"""

# ╔═╡ 823ee426-ca80-11f1-b481-c5c9b9792590
md"""
## 5. Del mercado a una red: el árbol de expansión mínima

Para dibujar el mercado como una red necesitamos una **distancia** entre
acciones (repaso L.I5):

```math
d_{ij} = \sqrt{2\,(1 - \rho_{ij})} .
```

Vale 0 si ``\rho_{ij} = 1`` y 2 si ``\rho_{ij} = -1``, y cumple la desigualdad
del triángulo. Con ``N`` acciones hay ``N(N-1)/2`` distancias, demasiadas
para dibujarlas. Mantegna (1999) propuso quedarse con el **árbol de expansión
mínima**: conectar las ``N`` acciones con ``N - 1`` enlaces, sin ciclos, de
modo que la suma de las distancias sea mínima.

El **algoritmo de Kruskal** lo construye así: ordena todos los pares de
menor a mayor distancia y recórrelos; agrega el enlace si une dos grupos que
todavía no estaban conectados, y descártalo si cerraría un ciclo.
"""

# ╔═╡ 823ee430-ca80-11f1-81d8-77076281b8b8
D = sqrt.(2 .* (1 .- C));

# ╔═╡ 823ee430-ca80-11f1-84e5-bf1bcdcbcdff
md"""
Resalta dos sectores:

Sector en azul: $(@bind sector_a Select(sort(unique(sectores.sector)), default = "Finanzas"))
  Sector en naranja: $(@bind sector_b Select(sort(unique(sectores.sector)), default = "Tecnología"))
"""

# ╔═╡ 823ee44e-ca80-11f1-a5e3-4d471590983b
md"""
## 6. [Avanzado] La ley de Marchenko-Pastur

Barajar da el espectro del ruido numéricamente. Para ``N`` series
independientes y ``N, T \to \infty`` con ``Q = T/N \ge 1`` fijo, la densidad
de eigenvalores tiene forma exacta (repaso L.A2):

```math
\rho(\lambda) = \frac{Q}{2\pi}\, \frac{\sqrt{(\lambda_+ - \lambda)(\lambda - \lambda_-)}}{\lambda},
\qquad \lambda_\pm = \left(1 \pm \frac{1}{\sqrt{Q}}\right)^2 ,
```

para ``\lambda_- \le \lambda \le \lambda_+``, y cero fuera. Con menos días
(``Q`` menor) el intervalo de ruido es más ancho y se tragan más
eigenvalores. Usa sólo los últimos ``T = QN`` días:

``Q = T/N``: $(@bind Q_mp Select([q => string(q) for q in (1.5, 2, 3, 5, 10, 20, round(T / N, digits = 1))], default = 5))
"""

# ╔═╡ 823ee46c-ca80-11f1-8038-5901ebcc0b46
md"""
## Resumen

| Herramienta | Qué revela en el S&P 100 |
|---|---|
| Matriz de correlación ``C`` | Casi todas las correlaciones son positivas; bloques por sector |
| Espectro barajado | El ruido de estimación; sólo unos cuantos eigenvalores son señal |
| ``v_1`` (modo de mercado) | Todas las acciones con el mismo signo; es el índice. Su peso sube en las crisis |
| ``v_2, v_3, \dots`` | Grupos económicos: defensivas vs. tecnología, bancos, semiconductores |
| Árbol de expansión mínima | Los sectores emergen de los datos, sin etiquetarlos |
| Marchenko-Pastur | La forma exacta del ruido para series independientes |
"""

# ╔═╡ 823ee46c-ca80-11f1-81ee-d3882dd60763
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `correlacion(R)`, que reciba una matriz de rendimientos
(``T`` filas, ``N`` columnas) y devuelva su matriz de correlación usando la
fórmula ``C = Z^\top Z/(T-1)`` de la sección 1, sin usar `cor`. Estandariza
cada columna con `mean(R; dims = 1)` y `std(R; dims = 1)`.
"""

# ╔═╡ 823ee478-ca80-11f1-a421-5d3f3b4bab54
function correlacion(R)
    missing  # reemplaza esta línea
end

# ╔═╡ 823ee478-ca80-11f1-b36e-87156de7e8dd
md"""
**2.** Escribe `distancia(ρ)`, que convierta una correlación en la distancia
de la sección 5.
"""

# ╔═╡ 823ee480-ca80-11f1-b5f4-fdb9b79e2106
function distancia(ρ)
    missing  # reemplaza esta línea
end

# ╔═╡ 823ee48a-ca80-11f1-84ae-7faef8d64005
md"""
### Ruta básica: interpretación

3. ¿Por qué barajar cada serie por separado destruye las correlaciones pero
   conserva la distribución de rendimientos de cada acción? ¿Qué pasaría si
   barajáramos los días de todas las acciones con la **misma** permutación?
4. Si el mercado explica el 36% de la varianza, ¿qué significa eso para un
   inversionista que cree estar diversificado porque tiene 95 acciones?
5. En la sección 3, ``\lambda_1/N`` subió en 2020. Explica con tus palabras
   por qué la diversificación falla justo en las crisis.
6. Elige una rama del árbol que no sea de un solo sector e investiga qué
   tienen en común sus empresas.
7. ¿Cómo afecta el sesgo de supervivencia a lo que vimos? Piensa en qué
   tipo de empresas salen del índice.

### Ruta avanzada

8. Para el modelo de un factor, ``C = (1-\rho)I + \rho\,\mathbf{1}\mathbf{1}^\top``,
   demuestra que ``\lambda_1 = 1 + (N-1)\rho`` con ``v_1 = \mathbf{1}/\sqrt{N}``
   y que los demás eigenvalores valen ``1 - \rho``. Usa la correlación
   promedio de la sección 1 para predecir ``\lambda_1`` y compárala con la
   observada.
9. Simula ``N = 100`` series normales independientes con ``T = 400`` y
   verifica la ley de Marchenko-Pastur. Repite con series ``t`` de Student
   con 3 grados de libertad: ¿cambia el borde ``\lambda_+``?
10. Quita el modo de mercado: regresa cada serie sobre la cartera ``v_1`` y
    calcula la matriz de correlación de los residuos. ¿Qué eigenvalores
    quedan fuera del ruido? ¿Cambia el árbol?
11. Construye el árbol con los datos de 2015–2019 y con los de 2022–2026.
    ¿Qué enlaces sobreviven? Propón una medida de cuánto cambió.

## Para la próxima sesión

- Lectura: Mantegna y Stanley (1999), capítulos 12 y 13. Complementaria:
  Slanina (2014), sección 1.4; Mantegna (1999). Para la ruta avanzada,
  Requião da Cunha (2021), sección 5.2; Laloux et al. (1999).
- En la sesión 8 empieza el Módulo II: la **distribución del dinero** como
  un problema de física estadística (máxima entropía y modelos de
  intercambio).
"""

# ╔═╡ 823ee49e-ca80-11f1-be1d-e3b642437755
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 823ee4a6-ca80-11f1-9fcd-cd7e3ef0258a
begin
    const AZUL       = "#2a78d6"   # serie principal o primer sector resaltado
    const NARANJA    = "#eb6834"   # segunda serie, teoría o segundo sector
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const RAMPA_SEC  = ["#f4f8fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
                        "#256abf", "#184f95", "#0d366b"]   # azul secuencial, 0 → 1
end;

# ╔═╡ 823ee3e0-ca80-11f1-a95a-739d80f8b883
let
    bordes = 0:0.1:7.5
    p = histogram(λ_barajado; bins = bordes, color = GRIS, alpha = 0.45, lw = 0,
                  label = "Datos barajados", xlabel = "eigenvalor λ", ylabel = "número de eigenvalores",
                  title = "Espectro de la matriz de correlación", xlims = (0, 7.5), legend = :topright)
    histogram!(p, filter(<(7.5), λ); bins = bordes, color = AZUL, alpha = 0.75, lw = 0,
               label = "Datos reales")
    vline!(p, [maximum(λ_barajado)]; color = GRIS, ls = :dash, lw = 1.5, label = "Máximo barajado")
    annotate!(p, 7.4, 4, text("λ₁ = $(round(λ[1], digits = 1)) →\nfuera de la gráfica", 8, :right, GRIS))
end

# ╔═╡ 823ee4a6-ca80-11f1-b4de-3d4aa73ad187
"Número entero con comas como separador de miles."
miles(x) = replace(string(round(Int, x)), r"(?<=\d)(?=(\d{3})+$)" => ",")

# ╔═╡ 823ee3a2-ca80-11f1-ace0-dfb419e29361
md"""
Con JPMorgan y Bank of America (dos bancos) la nube es una elipse alargada:
cuando uno sube, el otro tiende a subir. Prueba ahora con JPM y DUK (Duke
Energy, una eléctrica): la nube se vuelve más redonda. Y casi cualquier par
tiene correlación **positiva**: hay algo que mueve a todas las acciones a la
vez.

Con ``N = `` $(N) acciones hay ``N(N-1)/2 = `` $(miles(N * (N - 1) ÷ 2)) pares.
No podemos mirarlos uno por uno: necesitamos una manera de resumirlos.
"""

# ╔═╡ 823ee4a6-ca80-11f1-9e4e-5570f2b667a2
"Nombre de la empresa a partir de su ticker."
nombre_de(t) = sectores.nombre[findfirst(==(t), sectores.ticker)]

# ╔═╡ 823ee39a-ca80-11f1-9ac7-5159cd791c15
let
    i, j = findfirst(==(accion1), tickers), findfirst(==(accion2), tickers)
    ρ = cor(R[:, i], R[:, j])
    lim = quantile(abs.(vcat(R[:, i], R[:, j])), 0.999)
    scatter(100 .* R[:, i], 100 .* R[:, j]; ms = 2, msw = 0, color = AZUL, alpha = 0.3,
            xlims = (-100lim, 100lim), ylims = (-100lim, 100lim), aspect_ratio = :equal,
            xlabel = "$accion1: rendimiento diario (%)", ylabel = "$accion2: rendimiento diario (%)",
            title = "$(nombre_de(accion1)) y $(nombre_de(accion2)): correlación ρ = $(round(ρ, digits = 2))",
            size = (520, 500))
end

# ╔═╡ 823ee4b2-ca80-11f1-b7a9-e5cf44dfe042
"Orden de las acciones por sector. Los sectores pequeños van entre los grandes para que sus etiquetas no se encimen."
orden_sector = let
    orden = ["Tecnología", "Energía", "Finanzas", "Servicios públicos", "Salud", "Bienes raíces",
             "Industria", "Materiales", "Comunicaciones", "Consumo básico", "Consumo discrecional"]
    sortperm(collect(zip(sectores.sector, tickers)); by = x -> (findfirst(==(x[1]), orden), x[2]))
end;

# ╔═╡ 823ee4b2-ca80-11f1-ad83-ed1a077a35f7
"Centro y borde de cada sector en el orden `orden_sector`."
grupos = let
    s = sectores.sector[orden_sector]
    inicio = [1; findall(s[2:end] .!= s[1:end-1]) .+ 1]
    fin = [inicio[2:end] .- 1; length(s)]
    (sector = s[inicio], centro = (inicio .+ fin) ./ 2, borde = fin[1:end-1] .+ 0.5)
end;

# ╔═╡ 823ee3b8-ca80-11f1-8cc8-932b94815cdd
let
    o = orden_sector
    p1 = heatmap(C[o, o]; color = cgrad(RAMPA_SEC), clims = (0, 1), yflip = true,
                 aspect_ratio = :equal, xlims = (0.5, N + 0.5), ylims = (0.5, N + 0.5),
                 xticks = (grupos.centro, grupos.sector), xrotation = 60,
                 yticks = (grupos.centro, grupos.sector), tickfontsize = 7,
                 colorbar_title = "ρ", title = "Matriz de correlación (ordenada por sector)",
                 framestyle = :box)
    for b in grupos.borde
        vline!(p1, [b]; color = :white, lw = 1)
        hline!(p1, [b]; color = :white, lw = 1)
    end
    ρs = [C[i, j] for j in 1:N for i in 1:j-1]
    p2 = histogram(ρs; bins = 50, color = AZUL, lw = 0, alpha = 0.8, xlabel = "ρᵢⱼ", label = "",
                   ylabel = "número de pares", title = "Las $(miles(length(ρs))) correlaciones")
    vline!(p2, [mean(ρs)]; color = NARANJA, lw = 2, label = "Promedio = $(round(mean(ρs), digits = 2))",
           legend = :topright)
    plot(p1, p2; layout = grid(1, 2, widths = [0.6, 0.4]), size = (900, 480), bottom_margin = 8Plots.mm)
end

# ╔═╡ 823ee4bc-ca80-11f1-9b68-63f4e9bdcce0
"Barras de un vector de N componentes, agrupadas por sector con un hueco entre sectores."
function barras_por_sector(v; ylabel = "", titulo = "", referencia = 0)
    s = sectores.sector[orden_sector]
    x = cumsum([1; [s[i] == s[i-1] ? 1 : 3 for i in 2:N]])
    inicio = [1; findall(s[2:end] .!= s[1:end-1]) .+ 1]
    fin = [inicio[2:end] .- 1; N]
    centro = (x[inicio] .+ x[fin]) ./ 2
    a, b = extrema([v; referencia])
    p = bar(x, v[orden_sector]; color = AZUL, lw = 0, bar_width = 0.8, ylabel = ylabel,
            ylims = (min(a, 0) - 0.1abs(b - a), max(b, 0) + 0.1abs(b - a)),
            title = titulo, xticks = (centro, grupos.sector), xrotation = 35, size = (760, 380),
            xlims = (0, x[end] + 1), bottom_margin = 10Plots.mm)
    hline!(p, [referencia]; color = GRIS, lw = 1, ls = referencia == 0 ? :solid : :dash)
end

# ╔═╡ 823ee3f4-ca80-11f1-8469-2959bf5f94ab
barras_por_sector(sqrt(N) .* v1; ylabel = "√N · componente de v₁",
                  titulo = "Primer eigenvector: todas las acciones con el mismo signo",
                  referencia = 1)

# ╔═╡ 823ee41c-ca80-11f1-b7f4-d5bd61ac3eee
barras_por_sector(sqrt(N) .* espectro.vectors[:, end-k_vec+1]; ylabel = "√N · componente de v$(k_vec)",
                  titulo = "Eigenvector $k_vec (λ = $(round(λ[k_vec], digits = 2)))", referencia = 0)

# ╔═╡ 823ee4c6-ca80-11f1-92cd-0bd73f2a3d13
"Rendimiento diario del índice S&P 500 y de la cartera del modo de mercado, en las mismas fechas."
modo_vs_indice = let
    p = CSV.read(joinpath(@__DIR__, "..", "Data", "precios.csv"), DataFrame)
    sp = Dict(zip(p.fecha[p.activo .== "S&P 500"], p.cierre[p.activo .== "S&P 500"]))
    Z = (R .- mean(R; dims = 1)) ./ std(R; dims = 1)
    modo = Z * v1
    f0 = ancho.fecha
    k = [haskey(sp, f0[t]) && haskey(sp, f0[t+1]) for t in 1:T]
    (indice = [log(sp[f0[t+1]] / sp[f0[t]]) for t in findall(k)], modo = modo[k])
end;

# ╔═╡ 823ee3fe-ca80-11f1-95c8-fd427e12fee4
let
    m = modo_vs_indice
    lim = quantile(abs.(m.indice), 0.999)
    limy = quantile(abs.(m.modo), 0.999)
    scatter(100 .* m.indice, m.modo; ms = 2, msw = 0, color = AZUL, alpha = 0.3,
            xlims = (-100lim, 100lim), ylims = (-limy, limy), xlabel = "S&P 500: rendimiento diario (%)",
            ylabel = "cartera v₁ (rendimientos estandarizados)",
            title = "Modo de mercado vs. índice: correlación = $(round(cor(m.indice, m.modo), digits = 3))",
            size = (620, 400))
end

# ╔═╡ 823ee4c6-ca80-11f1-8421-9f073e52979d
"λ₁/N en ventanas móviles de `ventana` días, avanzando de 5 en 5."
function lambda1_movil(R, fechas, ventana)
    T, N = size(R)
    inicios = 1:5:T-ventana+1
    (fecha = [fechas[t+ventana-1] for t in inicios],
     fraccion = [eigmax(Symmetric(cor(R[t:t+ventana-1, :]))) / N for t in inicios])
end

# ╔═╡ 823ee408-ca80-11f1-8063-ad26ed159013
movil = lambda1_movil(R, fechas, ventana);

# ╔═╡ 823ee408-ca80-11f1-bc3e-316d4ad67163
let
    años = Date(2016):Year(2):Date(2026)
    p = plot(movil.fecha, movil.fraccion; color = AZUL, lw = 2, ylims = (0, 1),
             xticks = (Dates.value.(años), string.(year.(años))),
             xlabel = "fin de la ventana", ylabel = "λ₁ / N",
             title = "Fracción de la varianza explicada por el mercado (ventana de $ventana días)")
    i = argmax(movil.fraccion)
    scatter!(p, [movil.fecha[i]], [movil.fraccion[i]]; color = NARANJA, ms = 6, msw = 0)
    annotate!(p, movil.fecha[i], movil.fraccion[i] + 0.06,
              text("$(round(movil.fraccion[i], digits = 2)) al $(movil.fecha[i])", 8, GRIS))
end

# ╔═╡ 823ee4d0-ca80-11f1-ae6a-ed5661ccb30f
"""
Árbol de expansión mínima con el algoritmo de Kruskal. Usa una estructura de
conjuntos disjuntos (`jefe`) para saber si dos nodos ya están conectados.
"""
function kruskal(D)
    N = size(D, 1)
    pares = sort([(D[i, j], i, j) for j in 1:N for i in 1:j-1])
    jefe = collect(1:N)
    raiz(i) = jefe[i] == i ? i : (jefe[i] = raiz(jefe[i]))
    enlaces = Tuple{Int,Int}[]
    for (_, i, j) in pares
        a, b = raiz(i), raiz(j)
        a == b && continue          # cerraría un ciclo
        jefe[a] = b
        push!(enlaces, (i, j))
        length(enlaces) == N - 1 && break
    end
    enlaces
end

# ╔═╡ 823ee430-ca80-11f1-aaf1-0f5dcc7c07df
arbol = kruskal(D);   # lista de enlaces (i, j)

# ╔═╡ 823ee43a-ca80-11f1-9793-f52dceb0ec58
md"""
Sin decirle al algoritmo a qué sector pertenece cada acción, el árbol agrupa
a los bancos en una cadena, a la tecnología en ramas propias y a los
servicios públicos y el consumo básico juntos. En el centro aparece BRK-B
(Berkshire Hathaway), un conglomerado que posee partes de muchos sectores:
es casi un pequeño mercado.

¿Podría ser casualidad? Contemos cuántos de los $(length(arbol)) enlaces unen
acciones del **mismo sector** y comparemos con lo que pasaría si las
etiquetas de sector se repartieran al azar entre las acciones (otra prueba
de permutación).
"""

# ╔═╡ 823ee4d0-ca80-11f1-a1eb-8b3ab3186e0d
"""
Posiciones 2D para dibujar un árbol: escalamiento multidimensional clásico de
las distancias en la red, refinado con mayorización de estrés.
"""
function disposicion(enlaces, N; iter = 300)
    vecinos = [Int[] for _ in 1:N]
    for (i, j) in enlaces
        push!(vecinos[i], j); push!(vecinos[j], i)
    end
    δ = zeros(N, N)                     # distancia en número de enlaces
    for s in 1:N
        visto = fill(-1, N); visto[s] = 0; cola = [s]
        while !isempty(cola)
            u = popfirst!(cola)
            for v in vecinos[u]
                visto[v] < 0 && (visto[v] = visto[u] + 1; push!(cola, v))
            end
        end
        δ[s, :] = visto
    end
    J = I - fill(1 / N, N, N)
    F = eigen(Symmetric(-J * (δ .^ 2) * J / 2))
    X = F.vectors[:, end-1:end] .* sqrt.(F.values[end-1:end])'
    W = [i == j ? 0.0 : 1 / δ[i, j]^2 for i in 1:N, j in 1:N]
    Vp = pinv(Diagonal(vec(sum(W; dims = 2))) - W)
    for _ in 1:iter
        B = [i == j ? 0.0 : -W[i, j] * δ[i, j] / max(norm(X[i, :] - X[j, :]), 1e-9)
             for i in 1:N, j in 1:N]
        B[diagind(B)] .= -vec(sum(B; dims = 2))
        X = Vp * B * X
    end
    X
end

# ╔═╡ 823ee4d8-ca80-11f1-9a0c-f5bc0e40c6c5
posiciones_arbol = disposicion(arbol, N);

# ╔═╡ 823ee43a-ca80-11f1-85ca-1523f499cdfd
let
    X = posiciones_arbol
    color = [s == sector_a ? AZUL : s == sector_b ? NARANJA : GRIS_CLARO for s in sectores.sector]
    p = plot(; aspect_ratio = :equal, axis = false, grid = false, size = (820, 760),
             legend = :bottomleft, title = "Árbol de expansión mínima de $N acciones del S&P 100")
    for (i, j) in arbol
        plot!(p, X[[i, j], 1], X[[i, j], 2]; color = GRIS, lw = 0.8, label = "")
    end
    for (s, c) in ((sector_a, AZUL), (sector_b, NARANJA))
        k = sectores.sector .== s
        scatter!(p, X[k, 1], X[k, 2]; color = c, ms = 8, msw = 0, label = s)
    end
    k = .!(in.(sectores.sector, Ref((sector_a, sector_b))))
    scatter!(p, X[k, 1], X[k, 2]; color = GRIS_CLARO, ms = 8, msw = 0, label = "Otros sectores")
    annotate!(p, X[:, 1], X[:, 2] .+ 0.35, text.(tickers, 6, GRIS))
end

# ╔═╡ 823ee4d8-ca80-11f1-b7ee-c5a4bf39d155
"Enlaces del árbol dentro del mismo sector, observados y con etiquetas al azar."
prueba_sectores = let
    rng = Xoshiro(5)
    mismos(s) = count(((i, j),) -> s[i] == s[j], arbol)
    (observado = mismos(sectores.sector),
     azar = [mismos(shuffle(rng, sectores.sector)) for _ in 1:999])
end;

# ╔═╡ 823ee446-ca80-11f1-84e1-851b35df297b
let
    t = prueba_sectores
    p = histogram(t.azar; bins = -0.5:1:maximum(t.azar)+0.5, color = GRIS, alpha = 0.6, lw = 0,
                  label = "Sectores al azar (999 repeticiones)", xlims = (0, length(arbol)),
                  xlabel = "enlaces dentro del mismo sector", ylabel = "repeticiones",
                  title = "El árbol respeta los sectores", legend = :topright, size = (720, 300))
    vline!(p, [t.observado]; color = AZUL, lw = 3, label = "Árbol real: $(t.observado) de $(length(arbol))")
end

# ╔═╡ 823ee44e-ca80-11f1-8726-d358025fe479
md"""
Con sectores al azar, en promedio sólo
$(round(mean(prueba_sectores.azar), digits = 1)) enlaces unirían acciones del
mismo sector; el árbol real tiene $(prueba_sectores.observado). Ninguna de
las 999 permutaciones se acerca.

El árbol es nuestra primera **red**. En las sesiones 10 y 11 estudiaremos las redes
en general: grados, centralidad y robustez. Ahí veremos que, para comparar el
mercado con modelos de redes, conviene otra red construida con las mismas
correlaciones, porque un árbol no tiene triángulos.
"""

# ╔═╡ 823ee4e4-ca80-11f1-ad10-d523cc9f0f08
"Densidad de Marchenko-Pastur para Q = T/N."
function densidad_mp(λ, Q)
    λm, λp = (1 - 1 / sqrt(Q))^2, (1 + 1 / sqrt(Q))^2
    λm <= λ <= λp ? Q / (2π) * sqrt((λp - λ) * (λ - λm)) / λ : 0.0
end

# ╔═╡ 823ee4e4-ca80-11f1-a9b8-03617fd3db38
"Espectros real y barajado usando los últimos T = Q·N días."
function espectros_mp(R, Q; rng = Xoshiro(6))
    N = size(R, 2)
    T = min(round(Int, Q * N), size(R, 1))
    X = R[end-T+1:end, :]
    Xb = reduce(hcat, [shuffle(rng, X[:, j]) for j in 1:N])
    q = T / N
    (T = T, real = eigvals(Symmetric(cor(X))), barajado = eigvals(Symmetric(cor(Xb))),
     λmenos = (1 - 1 / sqrt(q))^2, λmas = (1 + 1 / sqrt(q))^2)
end

# ╔═╡ 823ee458-ca80-11f1-9745-2d10c6b53297
mp = espectros_mp(R, Q_mp);

# ╔═╡ 823ee458-ca80-11f1-93e6-7fcc794ec19a
let
    xmax = max(4, mp.λmas + 1)
    Δ = xmax / 60
    bordes = 0:Δ:xmax
    # densidad: fracción de los N eigenvalores en cada intervalo, entre su ancho
    densidad(λs) = [count(b -> bordes[k] <= b < bordes[k+1], λs) / (N * Δ) for k in 1:length(bordes)-1]
    db, dr = densidad(mp.barajado), densidad(mp.real)
    x = range(mp.λmenos, mp.λmas, length = 300)
    ymax = 1.1 * maximum([db; dr; densidad_mp.(x, Q_mp)])
    p = bar(bordes[1:end-1] .+ Δ / 2, db; bar_width = Δ, color = GRIS, alpha = 0.45, lw = 0,
            label = "Datos barajados", xlabel = "eigenvalor λ", ylabel = "densidad",
            xlims = (0, xmax), ylims = (0, ymax), legend = :topright,
            title = "Últimos $(mp.T) días (Q = $(Q_mp)): ruido teórico vs. datos")
    plot!(p, bordes, [dr; dr[end]]; seriestype = :steppost, color = AZUL, lw = 2, label = "Datos reales")
    plot!(p, x, densidad_mp.(x, Q_mp); color = NARANJA, lw = 2, label = "Marchenko-Pastur")
end

# ╔═╡ 823ee462-ca80-11f1-b626-334257c8b3af
md"""
Con ``Q = `` $(Q_mp), el ruido va de ``\lambda_- = `` $(round(mp.λmenos, digits = 2))
a ``\lambda_+ = `` $(round(mp.λmas, digits = 2)). Los datos barajados siguen la curva teórica, y
$(count(>(mp.λmas), mp.real)) eigenvalores reales quedan por encima de
``\lambda_+`` (algunos fuera de la gráfica). Laloux et al. (1999) usaron esta
comparación para mostrar que la mayor parte de una matriz de correlación
financiera es ruido, lo que importa mucho al construir carteras.

Dos sutilezas: la curva supone rendimientos con varianza finita y sin
dependencia temporal, y los rendimientos reales tienen colas pesadas y
volatilidad agrupada. Barajar conserva las colas; por eso es la comparación
más honesta.
"""

# ╔═╡ 823ee4e4-ca80-11f1-a9b8-4d451b1741f4
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

# ╔═╡ 823ee478-ca80-11f1-bd3b-699c128af538
verificar(() -> correlacion(R[1:300, 1:5]), cor(R[1:300, 1:5]);
          pista = md"Calcula `Z = (R .- mean(R; dims = 1)) ./ std(R; dims = 1)`; luego `Z' * Z` y divide entre `size(R, 1) - 1`.")

# ╔═╡ 823ee480-ca80-11f1-b0e3-65a77006585f
verificar(() -> distancia.([1.0, 0.5, 0.0, -1.0]), [0.0, 1.0, sqrt(2), 2.0];
          pista = md"La fórmula es ``d = \sqrt{2(1-\rho)}``. Para que funcione con `distancia.(…)` basta escribirla para un número.")

# ╔═╡ Cell order:
# ╟─823ee34a-ca80-11f1-906d-37eff8d968b3
# ╟─823ee370-ca80-11f1-995e-9b104afdea92
# ╟─823ee37c-ca80-11f1-a695-c5b70d702a47
# ╠═823ee390-ca80-11f1-aca8-6b8e059ea3a4
# ╟─823ee390-ca80-11f1-861c-bd9aacb29bff
# ╟─823ee39a-ca80-11f1-9ac7-5159cd791c15
# ╟─823ee3a2-ca80-11f1-ace0-dfb419e29361
# ╟─823ee3a2-ca80-11f1-95aa-537ae946b2c0
# ╠═823ee3ae-ca80-11f1-a874-ef132643d574
# ╟─823ee3b8-ca80-11f1-8cc8-932b94815cdd
# ╟─823ee3b8-ca80-11f1-a4f2-dbae20b60a02
# ╠═823ee3c2-ca80-11f1-8a23-2144a189a357
# ╠═823ee3c2-ca80-11f1-8a19-bf1625be78ea
# ╟─823ee3cc-ca80-11f1-9b7d-3933bc5ad931
# ╟─823ee3cc-ca80-11f1-92f7-291318cf2056
# ╠═823ee3d6-ca80-11f1-986f-c577c6ad93d8
# ╟─823ee3e0-ca80-11f1-a95a-739d80f8b883
# ╟─823ee3e0-ca80-11f1-83d2-f1c92b705e6f
# ╟─823ee3ea-ca80-11f1-82e4-43177fc88f50
# ╠═823ee3ea-ca80-11f1-81fc-233c2e60b46b
# ╟─823ee3f4-ca80-11f1-8469-2959bf5f94ab
# ╟─823ee3f4-ca80-11f1-8bac-89c9e80546f2
# ╟─823ee3fe-ca80-11f1-95c8-fd427e12fee4
# ╟─823ee3fe-ca80-11f1-86e9-5d5ebc9bab20
# ╟─823ee408-ca80-11f1-8063-ad26ed159013
# ╟─823ee408-ca80-11f1-bc3e-316d4ad67163
# ╟─823ee412-ca80-11f1-bb82-9b9320842c63
# ╟─823ee412-ca80-11f1-99fb-37271add34f8
# ╟─823ee41c-ca80-11f1-b7f4-d5bd61ac3eee
# ╟─823ee41c-ca80-11f1-8e2d-45ea67166fa7
# ╟─823ee426-ca80-11f1-b481-c5c9b9792590
# ╠═823ee430-ca80-11f1-81d8-77076281b8b8
# ╠═823ee430-ca80-11f1-aaf1-0f5dcc7c07df
# ╟─823ee430-ca80-11f1-84e5-bf1bcdcbcdff
# ╟─823ee43a-ca80-11f1-85ca-1523f499cdfd
# ╟─823ee43a-ca80-11f1-9793-f52dceb0ec58
# ╟─823ee446-ca80-11f1-84e1-851b35df297b
# ╟─823ee44e-ca80-11f1-8726-d358025fe479
# ╟─823ee44e-ca80-11f1-a5e3-4d471590983b
# ╟─823ee458-ca80-11f1-9745-2d10c6b53297
# ╟─823ee458-ca80-11f1-93e6-7fcc794ec19a
# ╟─823ee462-ca80-11f1-b626-334257c8b3af
# ╟─823ee46c-ca80-11f1-8038-5901ebcc0b46
# ╟─823ee46c-ca80-11f1-81ee-d3882dd60763
# ╠═823ee478-ca80-11f1-a421-5d3f3b4bab54
# ╟─823ee478-ca80-11f1-bd3b-699c128af538
# ╟─823ee478-ca80-11f1-b36e-87156de7e8dd
# ╠═823ee480-ca80-11f1-b5f4-fdb9b79e2106
# ╟─823ee480-ca80-11f1-b0e3-65a77006585f
# ╟─823ee48a-ca80-11f1-84ae-7faef8d64005
# ╟─823ee49e-ca80-11f1-be1d-e3b642437755
# ╠═823ee49e-ca80-11f1-acd2-3d414bc95fa4
# ╠═823ee4a6-ca80-11f1-9fcd-cd7e3ef0258a
# ╠═823ee4a6-ca80-11f1-b4de-3d4aa73ad187
# ╠═823ee4a6-ca80-11f1-9e4e-5570f2b667a2
# ╠═823ee4b2-ca80-11f1-b7a9-e5cf44dfe042
# ╠═823ee4b2-ca80-11f1-ad83-ed1a077a35f7
# ╠═823ee4bc-ca80-11f1-9b68-63f4e9bdcce0
# ╠═823ee4c6-ca80-11f1-92cd-0bd73f2a3d13
# ╠═823ee4c6-ca80-11f1-8421-9f073e52979d
# ╠═823ee4d0-ca80-11f1-ae6a-ed5661ccb30f
# ╠═823ee4d0-ca80-11f1-a1eb-8b3ab3186e0d
# ╠═823ee4d8-ca80-11f1-9a0c-f5bc0e40c6c5
# ╠═823ee4d8-ca80-11f1-b7ee-c5a4bf39d155
# ╠═823ee4e4-ca80-11f1-ad10-d523cc9f0f08
# ╠═823ee4e4-ca80-11f1-a9b8-03617fd3db38
# ╠═823ee4e4-ca80-11f1-a9b8-4d451b1741f4
