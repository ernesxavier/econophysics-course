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

# ╔═╡ 28289cae-c8f8-11f1-9965-85d9a28392eb
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

# ╔═╡ 28289a80-c8f8-11f1-b651-5101b43ca65e
md"""
# Sesión 5: Correlaciones y matrices aleatorias

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
`Code/s05_correlaciones.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): P.B2 (covarianza y
correlación), L.B1–B3 (vectores y matrices), S.B5 (pruebas por permutación),
L.I1–I3 (eigenvalores, matrices simétricas, componentes principales), L.I5
(distancia entre activos). Ruta avanzada: L.A2 (Marchenko-Pastur).
"""

# ╔═╡ 28289ad8-c8f8-11f1-a2fc-0f21169a9b0e
TableOfContents(title = "Contenido")

# ╔═╡ 28289ae0-c8f8-11f1-b686-19dca2480b14
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

# ╔═╡ 28289b00-c8f8-11f1-a022-39e539c6d5ff
begin
    precios = CSV.read(joinpath(@__DIR__, "..", "Data", "sp100_precios.csv"), DataFrame)
    sectores = CSV.read(joinpath(@__DIR__, "..", "Data", "sp100_sectores.csv"), DataFrame)
    ancho = unstack(precios, :fecha, :ticker, :cierre)    # una columna por acción
    tickers = names(ancho)[2:end]
    R = diff(log.(Matrix{Float64}(ancho[:, 2:end])); dims = 1)   # rendimientos: T × N
    fechas = ancho.fecha[2:end]
    T, N = size(R)
end;

# ╔═╡ 28289b00-c8f8-11f1-87da-391aa9ee3dd5
md"""
Tenemos ``T = `` $(T) días y ``N = `` $(N) acciones. Elige dos acciones y
compara sus rendimientos diarios, día por día:

Acción 1: $(@bind accion1 Select(tickers, default = "JPM"))
  Acción 2: $(@bind accion2 Select(tickers, default = "BAC"))
"""

# ╔═╡ 28289b28-c8f8-11f1-9245-997857c90de3
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

# ╔═╡ 28289b28-c8f8-11f1-a639-fdf7410167d7
C = cor(R);

# ╔═╡ 28289b3c-c8f8-11f1-b411-cbcaa67f1745
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

# ╔═╡ 28289b50-c8f8-11f1-ac22-038c6b3d1480
espectro = eigen(Symmetric(C));   # eigenvalores en orden creciente

# ╔═╡ 28289b5a-c8f8-11f1-8f42-b3fc734b8d85
λ = reverse(espectro.values);     # λ[1] es el mayor

# ╔═╡ 28289b5a-c8f8-11f1-b704-a5dd3cfbaf3e
md"""
Los cinco mayores eigenvalores son $(join(round.(λ[1:5], digits = 1), ", ")).
El primero, ``\lambda_1 = `` $(round(λ[1], digits = 1)), concentra el
$(round(Int, 100λ[1] / N))% de la varianza total. Pero ¿cuáles de los demás
son significativos? Con ``T`` finito, incluso acciones independientes darían
eigenvalores distintos de 1.
"""

# ╔═╡ 28289b6e-c8f8-11f1-b0cb-03710b5a67d1
md"""
## 2. Laboratorio: ¿señal o ruido?

La prueba es la de permutaciones de la sesión 1 (repaso S.B5): **barajamos
cada serie por separado**. Cada acción conserva exactamente sus rendimientos
(su media, su varianza y sus colas), pero los días quedan desalineados, así
que las correlaciones entre acciones desaparecen. Lo que sobresale del
espectro barajado es **señal**.

$(@bind barajar Button("Barajar otra vez"))
"""

# ╔═╡ 28289b78-c8f8-11f1-96b2-a39bc304cf23
λ_barajado = let
    barajar
    rng = Xoshiro()
    Rb = similar(R)
    for j in 1:N
        Rb[:, j] = shuffle(rng, R[:, j])
    end
    sort(eigvals(Symmetric(cor(Rb))); rev = true)
end;

# ╔═╡ 28289b84-c8f8-11f1-a8fa-23588ade43bf
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

# ╔═╡ 28289b8c-c8f8-11f1-a0d0-01bc01d0c5ce
md"""
## 3. El modo de mercado

El primer eigenvector ``v_1`` asigna a cada acción un peso. La gráfica los
muestra multiplicados por ``\sqrt{N}``, así que una cartera con el mismo peso
en todas las acciones daría barras de altura 1.
"""

# ╔═╡ 28289b96-c8f8-11f1-842e-65a5769e3628
v1 = espectro.vectors[:, end] .* sign(sum(espectro.vectors[:, end]));   # el signo es arbitrario

# ╔═╡ 28289baa-c8f8-11f1-9fa5-01c1b1ee42f5
md"""
Todas las componentes son positivas y parecidas: ``v_1`` es, casi, comprar
un poco de cada acción. Es el **modo de mercado**: el factor común que mueve
a todas a la vez. Para comprobarlo, comparemos el rendimiento diario de esa
cartera con el del índice S&P 500 completo (de `Data/precios.csv`):
"""

# ╔═╡ 28289bb6-c8f8-11f1-bfa7-bfa2238547c7
md"""
La correlación es casi 1: el primer eigenvector **es** el mercado, aunque no
le dijimos nada sobre el índice.

**Laboratorio: ¿cambia el modo de mercado con el tiempo?** Calculemos
``\lambda_1/N``, la fracción de la varianza que explica el mercado, en una
ventana móvil.

Longitud de la ventana (días hábiles): $(@bind ventana Select([60 => "60 (3 meses)", 125 => "125 (6 meses)", 250 => "250 (1 año)", 500 => "500 (2 años)"], default = 250))
"""

# ╔═╡ 28289bd2-c8f8-11f1-8bd3-c71a72b75817
md"""
En marzo de 2020 (la caída por la pandemia) el mercado llegó a explicar más
de la mitad de la varianza: **en las crisis todo cae junto**. Es una mala
noticia para la diversificación, que protege justo cuando las correlaciones
son bajas. Con ventanas cortas la serie es más ruidosa, porque cada matriz
se estima con menos días.
"""

# ╔═╡ 28289bdc-c8f8-11f1-9477-13acc0730dee
md"""
## 4. Sectores en los eigenvectores

Los demás eigenvectores significativos son ortogonales a ``v_1``: describen
cómo se mueven unas acciones **respecto a otras**, una vez descontado el
mercado. Mueve ``k`` para ver los eigenvectores que salieron del ruido. El
signo de cada eigenvector es arbitrario: lo que importa es qué acciones
tienen el mismo signo y cuáles el contrario.

Eigenvector ``k``: $(@bind k_vec Slider(2:8, default = 2, show_value = true))
"""

# ╔═╡ 28289be4-c8f8-11f1-96e3-256a8b9a3ea3
md"""
Cada eigenvector cuenta una historia económica:

- ``v_2`` opone las acciones **defensivas** (servicios públicos, consumo
  básico, farmacéuticas, telefónicas) a las de **tecnología**: cuando la
  tecnología sube, las defensivas tienden a quedarse atrás, y al revés.
- ``v_3`` separa a los **bancos** de las empresas de software.
- ``v_4`` separa a los **semiconductores** (AMAT, LRCX, MU, INTC) del
  software, dentro del mismo sector de tecnología.
- ``v_5`` reúne a las **farmacéuticas** y ``v_6`` a la **energía** con la
  industria de **defensa**.

Los sectores oficiales no son la última palabra: los datos encuentran grupos
que cruzan sectores (defensivas) o que los parten (semiconductores contra
software).
"""

# ╔═╡ 28289bf0-c8f8-11f1-ba1e-718f4af564df
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

# ╔═╡ 28289bfa-c8f8-11f1-a614-e11215f270fc
D = sqrt.(2 .* (1 .- C));

# ╔═╡ 28289c04-c8f8-11f1-a3c8-7fd0f2ae6e8d
md"""
Resalta dos sectores:

Sector en azul: $(@bind sector_a Select(sort(unique(sectores.sector)), default = "Finanzas"))
  Sector en naranja: $(@bind sector_b Select(sort(unique(sectores.sector)), default = "Tecnología"))
"""

# ╔═╡ 28289c36-c8f8-11f1-8d09-638aa0d7e324
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

# ╔═╡ 28289c5e-c8f8-11f1-9f9e-97f0911d4d49
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

# ╔═╡ 28289c68-c8f8-11f1-b55c-d9bb9aa65538
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `correlacion(R)`, que reciba una matriz de rendimientos
(``T`` filas, ``N`` columnas) y devuelva su matriz de correlación usando la
fórmula ``C = Z^\top Z/(T-1)`` de la sección 1, sin usar `cor`. Estandariza
cada columna con `mean(R; dims = 1)` y `std(R; dims = 1)`.
"""

# ╔═╡ 28289c72-c8f8-11f1-ac74-e3a059d2d812
function correlacion(R)
    missing  # reemplaza esta línea
end

# ╔═╡ 28289c7c-c8f8-11f1-b0cf-b7a7c9a5db52
md"""
**2.** Escribe `distancia(ρ)`, que convierta una correlación en la distancia
de la sección 5.
"""

# ╔═╡ 28289c86-c8f8-11f1-9973-fd1f0e2b0cac
function distancia(ρ)
    missing  # reemplaza esta línea
end

# ╔═╡ 28289c9a-c8f8-11f1-a3f8-dbfbc322c4fa
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

- Lectura: Mantegna y Stanley (1999), capítulos 12 y 13; Mantegna (1999).
- **Entrega**: examen para casa del Módulo I.
- En la semana 6 empieza el Módulo II: la **distribución del dinero** como
  un problema de física estadística (máxima entropía y modelos de
  intercambio).
"""

# ╔═╡ 28289ca4-c8f8-11f1-aae6-77fc40cb118d
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ 28289cba-c8f8-11f1-81ab-7bc248c4f5f3
begin
    const AZUL       = "#2a78d6"   # serie principal o primer sector resaltado
    const NARANJA    = "#eb6834"   # segunda serie, teoría o segundo sector
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const RAMPA_SEC  = ["#f4f8fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
                        "#256abf", "#184f95", "#0d366b"]   # azul secuencial, 0 → 1
end;

# ╔═╡ 28289b78-c8f8-11f1-8481-65c0c8c80e36
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

# ╔═╡ 28289cba-c8f8-11f1-a806-275494fe53fc
"Número entero con comas como separador de miles."
miles(x) = replace(string(round(Int, x)), r"(?<=\d)(?=(\d{3})+$)" => ",")

# ╔═╡ 28289b1e-c8f8-11f1-bc89-f19dce6122d0
md"""
Con JPMorgan y Bank of America (dos bancos) la nube es una elipse alargada:
cuando uno sube, el otro tiende a subir. Prueba ahora con JPM y DUK (Duke
Energy, una eléctrica): la nube se vuelve más redonda. Y casi cualquier par
tiene correlación **positiva**: hay algo que mueve a todas las acciones a la
vez.

Con ``N = `` $(N) acciones hay ``N(N-1)/2 = `` $(miles(N * (N - 1) ÷ 2)) pares.
No podemos mirarlos uno por uno: necesitamos una manera de resumirlos.
"""

# ╔═╡ 28289cc2-c8f8-11f1-aba6-ff0b911faf53
"Nombre de la empresa a partir de su ticker."
nombre_de(t) = sectores.nombre[findfirst(==(t), sectores.ticker)]

# ╔═╡ 28289b14-c8f8-11f1-ab22-e73f64b8b0e1
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

# ╔═╡ 28289ccc-c8f8-11f1-aae6-69f834ce23cf
"Orden de las acciones por sector. Los sectores pequeños van entre los grandes para que sus etiquetas no se encimen."
orden_sector = let
    orden = ["Tecnología", "Energía", "Finanzas", "Servicios públicos", "Salud", "Bienes raíces",
             "Industria", "Materiales", "Comunicaciones", "Consumo básico", "Consumo discrecional"]
    sortperm(collect(zip(sectores.sector, tickers)); by = x -> (findfirst(==(x[1]), orden), x[2]))
end;

# ╔═╡ 28289cd6-c8f8-11f1-886c-df97eda44f40
"Centro y borde de cada sector en el orden `orden_sector`."
grupos = let
    s = sectores.sector[orden_sector]
    inicio = [1; findall(s[2:end] .!= s[1:end-1]) .+ 1]
    fin = [inicio[2:end] .- 1; length(s)]
    (sector = s[inicio], centro = (inicio .+ fin) ./ 2, borde = fin[1:end-1] .+ 0.5)
end;

# ╔═╡ 28289b32-c8f8-11f1-8d2e-4dcdce2c6937
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

# ╔═╡ 28289cd6-c8f8-11f1-a755-4912d15936a1
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

# ╔═╡ 28289ba0-c8f8-11f1-9dc2-0f34ac969ffd
barras_por_sector(sqrt(N) .* v1; ylabel = "√N · componente de v₁",
                  titulo = "Primer eigenvector: todas las acciones con el mismo signo",
                  referencia = 1)

# ╔═╡ 28289be4-c8f8-11f1-847e-d53c57f01fa5
barras_por_sector(sqrt(N) .* espectro.vectors[:, end-k_vec+1]; ylabel = "√N · componente de v$(k_vec)",
                  titulo = "Eigenvector $k_vec (λ = $(round(λ[k_vec], digits = 2)))", referencia = 0)

# ╔═╡ 28289ce0-c8f8-11f1-87be-07af1780c37e
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

# ╔═╡ 28289baa-c8f8-11f1-98aa-ed0f3fba8658
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

# ╔═╡ 28289cec-c8f8-11f1-970f-e394d5fe5734
"λ₁/N en ventanas móviles de `ventana` días, avanzando de 5 en 5."
function lambda1_movil(R, fechas, ventana)
    T, N = size(R)
    inicios = 1:5:T-ventana+1
    (fecha = [fechas[t+ventana-1] for t in inicios],
     fraccion = [eigmax(Symmetric(cor(R[t:t+ventana-1, :]))) / N for t in inicios])
end

# ╔═╡ 28289bbe-c8f8-11f1-9028-0f58433fbe4e
movil = lambda1_movil(R, fechas, ventana);

# ╔═╡ 28289bc8-c8f8-11f1-97b2-ef654a8ce200
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

# ╔═╡ 28289cf4-c8f8-11f1-b1ba-e5cfd7688d59
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

# ╔═╡ 28289c04-c8f8-11f1-9cb3-290cbfb8b03d
arbol = kruskal(D);   # lista de enlaces (i, j)

# ╔═╡ 28289c22-c8f8-11f1-8551-2d30c40db5d7
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

# ╔═╡ 28289cfe-c8f8-11f1-8853-8df07d76d644
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

# ╔═╡ 28289d08-c8f8-11f1-a70d-d31d3931b490
posiciones_arbol = disposicion(arbol, N);

# ╔═╡ 28289c0e-c8f8-11f1-8d32-719353243497
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

# ╔═╡ 28289d12-c8f8-11f1-87cc-19a1cc27020c
"Enlaces del árbol dentro del mismo sector, observados y con etiquetas al azar."
prueba_sectores = let
    rng = Xoshiro(5)
    mismos(s) = count(((i, j),) -> s[i] == s[j], arbol)
    (observado = mismos(sectores.sector),
     azar = [mismos(shuffle(rng, sectores.sector)) for _ in 1:999])
end;

# ╔═╡ 28289c22-c8f8-11f1-9d8b-e1e0f88162c7
let
    t = prueba_sectores
    p = histogram(t.azar; bins = -0.5:1:maximum(t.azar)+0.5, color = GRIS, alpha = 0.6, lw = 0,
                  label = "Sectores al azar (999 repeticiones)", xlims = (0, length(arbol)),
                  xlabel = "enlaces dentro del mismo sector", ylabel = "repeticiones",
                  title = "El árbol respeta los sectores", legend = :topright, size = (720, 300))
    vline!(p, [t.observado]; color = AZUL, lw = 3, label = "Árbol real: $(t.observado) de $(length(arbol))")
end

# ╔═╡ 28289c2c-c8f8-11f1-86cd-0d6e718c3227
md"""
Con sectores al azar, en promedio sólo
$(round(mean(prueba_sectores.azar), digits = 1)) enlaces unirían acciones del
mismo sector; el árbol real tiene $(prueba_sectores.observado). Ninguna de
las 999 permutaciones se acerca.

El árbol es nuestra primera **red**. En la semana 8 lo usaremos para medir
grados, centralidad y robustez.
"""

# ╔═╡ 28289d1e-c8f8-11f1-b541-21ae983c55b0
"Densidad de Marchenko-Pastur para Q = T/N."
function densidad_mp(λ, Q)
    λm, λp = (1 - 1 / sqrt(Q))^2, (1 + 1 / sqrt(Q))^2
    λm <= λ <= λp ? Q / (2π) * sqrt((λp - λ) * (λ - λm)) / λ : 0.0
end

# ╔═╡ 28289d1e-c8f8-11f1-9419-f749c2da1bf0
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

# ╔═╡ 28289c40-c8f8-11f1-9d39-3104f263b874
mp = espectros_mp(R, Q_mp);

# ╔═╡ 28289c48-c8f8-11f1-9b83-45b766229077
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

# ╔═╡ 28289c54-c8f8-11f1-adb9-6f46f9a5796b
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

# ╔═╡ 28289d26-c8f8-11f1-a3df-07c80896e5c3
"Muestra si la respuesta de un ejercicio es correcta."
function verificar(f, esperado; rtol = 1e-8, pista = md"")
    caja(tipo, titulo, texto) = Markdown.MD(Markdown.Admonition(tipo, titulo, [texto]))
    resultado = try
        f()
    catch e
        return caja("danger", "Tu código produce un error", md"$(sprint(showerror, e))")
    end
    if resultado === missing
        caja("warning", "Pendiente", md"Reemplaza `missing` por tu código.")
    elseif resultado isa typeof(esperado) && size(resultado) == size(esperado) &&
           all(isapprox.(resultado, esperado; rtol, atol = 1e-12))
        caja("correct", "¡Correcto!", md"Bien hecho.")
    else
        caja("danger", "Todavía no", pista)
    end
end

# ╔═╡ 28289c7c-c8f8-11f1-9e87-9d4678cb071a
verificar(() -> correlacion(R[1:300, 1:5]), cor(R[1:300, 1:5]);
          pista = md"Calcula `Z = (R .- mean(R; dims = 1)) ./ std(R; dims = 1)`; luego `Z' * Z` y divide entre `size(R, 1) - 1`.")

# ╔═╡ 28289c90-c8f8-11f1-a9ac-0b7dfefb490f
verificar(() -> distancia.([1.0, 0.5, 0.0, -1.0]), [0.0, 1.0, sqrt(2), 2.0];
          pista = md"La fórmula es ``d = \sqrt{2(1-\rho)}``. Para que funcione con `distancia.(…)` basta escribirla para un número.")

# ╔═╡ Cell order:
# ╟─28289a80-c8f8-11f1-b651-5101b43ca65e
# ╟─28289ad8-c8f8-11f1-a2fc-0f21169a9b0e
# ╟─28289ae0-c8f8-11f1-b686-19dca2480b14
# ╠═28289b00-c8f8-11f1-a022-39e539c6d5ff
# ╟─28289b00-c8f8-11f1-87da-391aa9ee3dd5
# ╟─28289b14-c8f8-11f1-ab22-e73f64b8b0e1
# ╟─28289b1e-c8f8-11f1-bc89-f19dce6122d0
# ╟─28289b28-c8f8-11f1-9245-997857c90de3
# ╠═28289b28-c8f8-11f1-a639-fdf7410167d7
# ╟─28289b32-c8f8-11f1-8d2e-4dcdce2c6937
# ╟─28289b3c-c8f8-11f1-b411-cbcaa67f1745
# ╠═28289b50-c8f8-11f1-ac22-038c6b3d1480
# ╠═28289b5a-c8f8-11f1-8f42-b3fc734b8d85
# ╟─28289b5a-c8f8-11f1-b704-a5dd3cfbaf3e
# ╟─28289b6e-c8f8-11f1-b0cb-03710b5a67d1
# ╠═28289b78-c8f8-11f1-96b2-a39bc304cf23
# ╟─28289b78-c8f8-11f1-8481-65c0c8c80e36
# ╟─28289b84-c8f8-11f1-a8fa-23588ade43bf
# ╟─28289b8c-c8f8-11f1-a0d0-01bc01d0c5ce
# ╠═28289b96-c8f8-11f1-842e-65a5769e3628
# ╟─28289ba0-c8f8-11f1-9dc2-0f34ac969ffd
# ╟─28289baa-c8f8-11f1-9fa5-01c1b1ee42f5
# ╟─28289baa-c8f8-11f1-98aa-ed0f3fba8658
# ╟─28289bb6-c8f8-11f1-bfa7-bfa2238547c7
# ╟─28289bbe-c8f8-11f1-9028-0f58433fbe4e
# ╟─28289bc8-c8f8-11f1-97b2-ef654a8ce200
# ╟─28289bd2-c8f8-11f1-8bd3-c71a72b75817
# ╟─28289bdc-c8f8-11f1-9477-13acc0730dee
# ╟─28289be4-c8f8-11f1-847e-d53c57f01fa5
# ╟─28289be4-c8f8-11f1-96e3-256a8b9a3ea3
# ╟─28289bf0-c8f8-11f1-ba1e-718f4af564df
# ╠═28289bfa-c8f8-11f1-a614-e11215f270fc
# ╠═28289c04-c8f8-11f1-9cb3-290cbfb8b03d
# ╟─28289c04-c8f8-11f1-a3c8-7fd0f2ae6e8d
# ╟─28289c0e-c8f8-11f1-8d32-719353243497
# ╟─28289c22-c8f8-11f1-8551-2d30c40db5d7
# ╟─28289c22-c8f8-11f1-9d8b-e1e0f88162c7
# ╟─28289c2c-c8f8-11f1-86cd-0d6e718c3227
# ╟─28289c36-c8f8-11f1-8d09-638aa0d7e324
# ╟─28289c40-c8f8-11f1-9d39-3104f263b874
# ╟─28289c48-c8f8-11f1-9b83-45b766229077
# ╟─28289c54-c8f8-11f1-adb9-6f46f9a5796b
# ╟─28289c5e-c8f8-11f1-9f9e-97f0911d4d49
# ╟─28289c68-c8f8-11f1-b55c-d9bb9aa65538
# ╠═28289c72-c8f8-11f1-ac74-e3a059d2d812
# ╟─28289c7c-c8f8-11f1-9e87-9d4678cb071a
# ╟─28289c7c-c8f8-11f1-b0cf-b7a7c9a5db52
# ╠═28289c86-c8f8-11f1-9973-fd1f0e2b0cac
# ╟─28289c90-c8f8-11f1-a9ac-0b7dfefb490f
# ╟─28289c9a-c8f8-11f1-a3f8-dbfbc322c4fa
# ╟─28289ca4-c8f8-11f1-aae6-77fc40cb118d
# ╠═28289cae-c8f8-11f1-9965-85d9a28392eb
# ╠═28289cba-c8f8-11f1-81ab-7bc248c4f5f3
# ╠═28289cba-c8f8-11f1-a806-275494fe53fc
# ╠═28289cc2-c8f8-11f1-aba6-ff0b911faf53
# ╠═28289ccc-c8f8-11f1-aae6-69f834ce23cf
# ╠═28289cd6-c8f8-11f1-886c-df97eda44f40
# ╠═28289cd6-c8f8-11f1-a755-4912d15936a1
# ╠═28289ce0-c8f8-11f1-87be-07af1780c37e
# ╠═28289cec-c8f8-11f1-970f-e394d5fe5734
# ╠═28289cf4-c8f8-11f1-b1ba-e5cfd7688d59
# ╠═28289cfe-c8f8-11f1-8853-8df07d76d644
# ╠═28289d08-c8f8-11f1-a70d-d31d3931b490
# ╠═28289d12-c8f8-11f1-87cc-19a1cc27020c
# ╠═28289d1e-c8f8-11f1-b541-21ae983c55b0
# ╠═28289d1e-c8f8-11f1-9419-f749c2da1bf0
# ╠═28289d26-c8f8-11f1-a3df-07c80896e5c3
