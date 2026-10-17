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

# ╔═╡ ba30de1a-c9eb-11f1-b180-47208819aec9
begin
    # Usa el entorno del curso (Code/Project.toml) en lugar del gestor de
    # paquetes propio de Pluto: un solo entorno para todo el repositorio.
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using PlutoUI, Statistics, Random, Distributions, Graphs, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

# ╔═╡ ba30dcf8-c9eb-11f1-9110-e936cab76302
md"""
# Sesión 11: Opinión, herding y segregación

**Objetivos de la sesión**
1. Construir **modelos basados en agentes (ABM)** en los que reglas locales de
   imitación o mudanza generan patrones macroscópicos sin coordinación central.
2. Simular el **modelo de votante**, medir su **tiempo de consenso** y comparar
   el efecto de distintas topologías de red.
3. Incorporar cambios autónomos de opinión en el **modelo de hormigas de
   Kirman** para explicar el **herding** y la transición hacia una distribución
   **bimodal**.
4. Simular el **modelo de segregación de Schelling** y comprobar cómo
   preferencias individuales leves producen segregación residencial fuerte.
5. **[Avanzado]** Contrastar la simulación con el tiempo exacto de absorción del
   votante y con la solución estacionaria de la ecuación maestra de Kirman.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s11_opinion.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): E.B2–B3 (caminata aleatoria
y promedios), L.B4 (matriz de adyacencia), C.B4 (función logística), F.B3–B4
(temperatura y transiciones de fase), E.I3 (cadenas de Markov), F.I3–I4
(Ising y campo medio). Ruta avanzada: E.A3 (ecuación maestra), E.A6 (tiempos
de primer paso y absorción), L.A1 (Laplaciano).
"""

# ╔═╡ ba30dd2a-c9eb-11f1-a955-f795c27e6017
TableOfContents(title = "Contenido")

# ╔═╡ ba30dd34-c9eb-11f1-9106-7533c1d56ae3
md"""
## 0. El fenómeno: intenciones individuales y resultados colectivos

En el Módulo II estudiamos sistemas en equilibrio estadístico, como el modelo
de Ising. Hoy abrimos el **Módulo III** con modelos donde los agentes siguen
reglas directas de comportamiento —copiar a un vecino, seguir a la mayoría o
mudarse de casa— sin resolver un problema de optimización global.

Dos acertijos clásicos muestran por qué el agregado no es la simple suma de las
preferencias individuales:

1. **Dos restaurantes idénticos** (Becker, 1991; Kirman, 1993). Dos locales con
   la misma comida y los mismos precios están uno frente al otro. En vez de
   repartirse la clientela mitad y mitad, uno se llena y el otro queda vacío;
   tiempo después, sin que cambie la calidad de ninguno, la multitud se pasa al
   otro. Los biólogos observaron exactamente lo mismo con hormigas frente a dos
   fuentes idénticas de alimento.
2. **Barrios segregados sin intolerancia** (Schelling, 1971). En una ciudad con
   dos grupos, ningún habitante pide vivir rodeado únicamente de su propio
   grupo: cada uno se conforma con que al menos **tres de sus ocho vecinos**
   (37.5%) sean como él. Sin embargo, cuando los insatisfechos se mudan, la
   ciudad termina partida en barrios casi homogéneos.
"""

# ╔═╡ ba30dd52-c9eb-11f1-b7ce-61b65d44eeaf
md"""
## 1. Repaso: estados absorbentes frente a ruido

En una cadena de Markov donde los agentes **únicamente copian** las opiniones
existentes, si en algún momento todos coinciden (``k = 0`` o ``k = N`` agentes
con la opinión 1), ya no queda nadie con la opinión contraria a quien copiar.
Esos dos consensos son **estados absorbentes** (repaso E.A6): una vez que el
sistema entra en uno, se queda ahí para siempre, y la pregunta relevante es
**cuánto tarda en absorberse**.

Si además cada agente tiene una pequeña probabilidad de cambiar de opinión por
cuenta propia (**ruido idiosincrático**, repaso E.I3 y F.B3), los extremos
dejan de ser absorbentes: siempre puede reaparecer un disidente. El sistema ya
no se detiene y converge a una **distribución estacionaria** ``P_k``.
"""

# ╔═╡ ba30dd52-c9eb-11f1-bc01-db539696fd6d
# En imitación pura, la probabilidad de cambiar cuando k de N tienen opinión 1
# es proporcional a k(N - k): se anula en k = 0 y en k = N.
tasa_imitacion_pura(k, N) = (k / N) * ((N - k) / (N - 1))

# ╔═╡ ba30dd5c-c9eb-11f1-9249-51a366d88393
tasa_imitacion_pura(0, 100), tasa_imitacion_pura(50, 100), tasa_imitacion_pura(100, 100)

# ╔═╡ ba30dd66-c9eb-11f1-8675-f30b9c80414d
md"""
## 2. El modelo de votante

El **modelo de votante** (Clifford y Sudbury, 1973; Holley y Liggett, 1975) es
el modelo más simple de formación de opinión por imitación pura.

**Reglas.**
1. Hay ``N`` agentes en los nodos de una red conexa; cada uno tiene opinión
   ``s_i \in \{0, 1\}``.
2. En cada paso elemental se elige un agente ``i`` al azar y uno de sus vecinos
   ``j`` al azar; ``i`` copia la opinión de ``j`` (``s_i \leftarrow s_j``).
3. ``N`` pasos elementales forman un **barrido** (una unidad de tiempo, en la
   que cada agente actualiza su opinión una vez en promedio).

En una red donde todos están conectados con todos (red completa), si hay ``k``
agentes con opinión 1, en un paso elemental ``k`` sube en 1 si se elige un
agente en 0 y copia a uno en 1, y baja en 1 en el caso simétrico:

```math
\mathbb{P}(k \to k+1) = \frac{N-k}{N}\frac{k}{N-1}, \qquad
\mathbb{P}(k \to k-1) = \frac{k}{N}\frac{N-k}{N-1}.
```

Como las dos probabilidades son **iguales**, el número ``k_t`` de agentes con
opinión 1 no tiene tendencia: es una **martingala** (repaso E.I2), igual que la
ruina del jugador (repaso E.A6). De ahí se siguen dos consecuencias inmediatas:

- En cualquier red finita conexa, el modelo **siempre llega al consenso**
  (``k = 0`` o ``k = N``).
- Si empieza con una fracción ``x_0 = k_0 / N`` de agentes en 1, la
  probabilidad de terminar en el consenso ``k = N`` es exactamente ``x_0``.

Fracción inicial ``x_0``: $(@bind x0_vot Slider(0.1:0.1:0.9, default = 0.3, show_value = true))
  Número de agentes ``N``: $(@bind N_vot Slider([40, 80, 150], default = 80, show_value = true))
"""

# ╔═╡ ba30dd7a-c9eb-11f1-9d97-59567d04119d
md"""
### 2.1 Laboratorio: ¿cómo afecta la red al tiempo de consenso?

En la red completa, el tiempo promedio hasta el consenso desde ``x_0 = 1/2``
crece linealmente con el número de agentes:
``\mathbb{E}[\tau] \approx (N - 1)\ln 2 \approx 0.693\,N`` barridos. ¿Qué pasa
cuando los agentes sólo ven a sus vecinos en una red?

Comparamos cuatro redes con ``N = 144`` nodos, todas empezando mitad y mitad
(``x_0 = 0.5``):
- **Anillo 1D** (cada nodo tiene 2 vecinos);
- **Cuadrícula 2D** de ``12 \times 12`` (4 vecinos);
- **Erdős-Rényi** conexa (grado promedio ``\approx 4``);
- **Barabási-Albert** (grado promedio ``\approx 4``).
"""

# ╔═╡ ba30dd8e-c9eb-11f1-8b3d-b36fc8f1d91f
md"""
La estructura de la red cambia radicalmente la velocidad del consenso:

- En el **anillo 1D**, se forman dominios de 0 y de 1 cuyas fronteras caminan al
  azar hasta chocar; recorrer una distancia de orden ``N`` por difusión toma del
  orden de ``N^2`` barridos (repaso E.B2). Por eso tarda **muchísimo más** que
  las demás redes.
- En **Erdős-Rényi** y la **cuadrícula 2D**, los atajos o las dos dimensiones
  evitan que las fronteras se estanquen y el tiempo vuelve a ser de orden ``N``
  (en 2D, ``N \ln N``).
- En **Barabási-Albert**, los concentradores (*hubs*) son imitados por muchos
  vecinos a la vez y aceleran todavía más la llegada al consenso.
"""

# ╔═╡ ba30dd8e-c9eb-11f1-8021-b31c7b64e531
md"""
## 3. El modelo de hormigas de Kirman: herding y bimodalidad

Las hormigas de Kirman (1993) —y los inversionistas de un mercado financiero— no
se quedan congelados para siempre en un consenso: siempre hay una pequeña
probabilidad de cambiar de opción de manera autónoma.

**Reglas.** Hay ``N`` agentes con estado ``0`` o ``1`` (por ejemplo, pesimistas
y optimistas, o fundamentalistas y chartistas). Sea ``k`` el número de agentes
en el estado ``1``. En cada paso se elige un agente al azar. Al decidir su
estado:
- Con peso ``N - 1`` imita a uno de los demás agentes elegido al azar (**herding**).
- Con peso ``2a`` actúa por cuenta propia y elige ``0`` o ``1`` con igual
  probabilidad (es decir, peso ``a`` a favor de cada opción).

Si el agente elegido estaba en ``0`` (probabilidad ``(N-k)/N``), pasa a ``1``
con probabilidad ``(a + k)/(N - 1 + 2a)``. Así, las probabilidades de transición
en un paso son:

```math
W(k \to k+1) = \frac{N-k}{N}\,\frac{a + k}{N - 1 + 2a}, \qquad
W(k \to k-1) = \frac{k}{N}\,\frac{a + N - k}{N - 1 + 2a}.
```

El parámetro clave es ``a > 0``, que mide la fuerza de la **decisión autónoma
frente a la imitación** (cuando ``a \to 0`` recuperamos el modelo de votante).

### Laboratorio: de la campana a las olas de opinión

Mueve ``a`` y observa tanto la trayectoria en el tiempo de la fracción
``x_t = k_t / N`` como el histograma de los valores que visita.

Cambio autónomo frente a imitación ``a``: $(@bind a_kirman Slider([0.2, 0.5, 1.0, 2.0, 5.0], default = 0.2, show_value = true))
"""

# ╔═╡ ba30dda2-c9eb-11f1-ba7f-89022dd9cd47
md"""
Observa las tres fases del modelo al variar ``a``:

- **``a < 1`` (herding fuerte, p. ej. ``a = 0.2``):** La distribución es
  **bimodal**, con picos en ``x \approx 0`` y ``x \approx 1``. La mayoría se
  concentra en una opción durante largos periodos y, de pronto, una racha de
  conversiones voltea a toda la población hacia la otra opción. **No hubo
  ninguna noticia externa**: los cambios de régimen son endógenos.
- **``a = 1`` (umbral crítico):** La distribución estacionaria es **exactamente
  uniforme**: todas las fracciones ``k/N`` son igual de probables.
- **``a > 1`` (decisión autónoma fuerte, p. ej. ``a = 5``):** La imitación no
  alcanza a arrastrar al grupo y la fracción fluctúa alrededor de ``1/2``
  (**unimodal**).

**Diferencia con el modelo de Ising (sesión 10).** En el Ising, cada agente
responde de forma no lineal (logística) al promedio; en Kirman, la tasa de
imitación es **lineal** en ``k/N``. Por eso en Kirman la tendencia promedio
siempre apunta hacia ``x = 1/2``: la bimodalidad cuando ``a < 1`` proviene de
que el **ruido de imitación se apaga cerca de los bordes** ``x = 0`` y
``x = 1`` (donde casi todos opinan igual), atrapando temporalmente al sistema
en consensos cuasi-absorbentes.
"""

# ╔═╡ ba30ddac-c9eb-11f1-934f-4f29b55971fc
md"""
## 4. El modelo de segregación de Schelling

En los modelos de votante y de Kirman los agentes están fijos y cambian de
estado. En el modelo de **Schelling** (1971) ocurre lo contrario: cada agente
conserva su identidad (azul o naranja) pero **cambia de lugar** si no le gusta
su vecindario.

**Reglas.**
1. En una cuadrícula de ``L \times L`` con fronteras periódicas, una fracción
   ``\rho_0`` de las celdas está vacía y el resto se reparte en partes iguales
   entre dos grupos (azul y naranja).
2. Cada agente mira sus 8 celdas vecinas (vecindad de Moore). Si entre sus
   vecinos ocupados la fracción de su mismo grupo es al menos el **umbral de
   tolerancia** ``\tau``, el agente está **satisfecho**; si es menor, está
   **insatisfecho** (si no tiene vecinos ocupados, está satisfecho).
3. En cada ronda, los agentes insatisfechos se mudan, en orden aleatorio, a
   celdas vacías elegidas al azar. El proceso se detiene cuando todos están
   satisfechos o al llegar al máximo de rondas.

### Laboratorio: de las preferencias individuales a la ciudad segregada

Mueve el umbral ``\tau``: la fracción mínima de vecinos propios que pide cada
agente para no mudarse.

Umbral requerido ``\tau``: $(@bind τ_sch Slider([0.125, 0.25, 0.375, 0.50, 0.625, 0.75, 0.85], default = 0.375, show_value = true))
  Fracción de celdas vacías ``\rho_0``: $(@bind vac_sch Slider([0.05, 0.15, 0.30], default = 0.15, show_value = true))
"""

# ╔═╡ ba30ddb6-c9eb-11f1-97c8-ff604b3e1731
md"""
Para ver la transición completa, medimos la **similitud promedio final** y la
**fracción de agentes satisfechos** al variar el umbral ``\tau`` de 0 a 0.85:
"""

# ╔═╡ ba30ddd4-c9eb-11f1-963a-552a0eb453a0
md"""
## 5. [Avanzado] Ecuación maestra y límite continuo de Kirman

En el modelo de Kirman, la probabilidad estacionaria ``P_k`` de observar ``k``
agentes en el estado 1 cumple la condición de **balance detallado** (repaso
E.A3):

```math
P_k\,W(k \to k+1) = P_{k+1}\,W(k+1 \to k)
\quad\Longrightarrow\quad
\frac{P_{k+1}}{P_k} = \frac{N - k}{k + 1}\,\frac{k + a}{N - k - 1 + a}.
```

Esta recurrencia tiene solución exacta para cualquier ``N``: la distribución
**Beta-binomial** (o de Pólya),

```math
P_k = \binom{N}{k}\frac{\mathrm{B}(k + a,\, N - k + a)}{\mathrm{B}(a, a)}.
```

Cuando ``N \to \infty`` con ``a`` fijo, y midiendo el tiempo en unidades de
``N(N - 1 + 2a)/2 \approx N^2/2`` pasos (unos ``N/2`` barridos), la variable
continua ``x = k / N \in (0, 1)`` cumple la ecuación de Fokker-Planck (repaso
E.A4) con deriva ``A(x) = a(1/2 - x)`` y coeficiente de difusión
``B(x) = x(1-x)``. Su densidad estacionaria es la **distribución Beta**
simétrica:

```math
p(x) = \frac{\Gamma(2a)}{\Gamma(a)^2}\,x^{a-1}(1 - x)^{a-1}.
```

La gráfica compara la distribución estacionaria exacta ``N \cdot P_k`` con la
densidad continua ``\mathrm{Beta}(a, a)`` para los tres regímenes ``a = 0.4``,
``a = 1`` y ``a = 3``.
"""

# ╔═╡ ba30dde8-c9eb-11f1-a088-2f2baabc45b9
md"""
## Resumen

| Modelo | Regla microscópica | Comportamiento colectivo |
|---|---|---|
| **Votante** | Copiar a un vecino al azar (sin ruido) | Consenso absorbente; ``\mathbb{P}(\text{gana } 1) = x_0``; ``\mathbb{E}[\tau] \propto N^2`` en 1D y ``\propto N`` en red completa |
| **Kirman** | Imitación (herding) + cambio autónomo ``a`` | Sin estados absorbentes; estacionaria ``\mathrm{Beta}(a, a)``: bimodal si ``a < 1``, uniforme si ``a = 1``, unimodal si ``a > 1`` |
| **Schelling** | Mudarse si la fracción de vecinos similares es ``< \tau`` | Preferencias leves (``\tau \approx 3/8``) desencadenan segregación fuerte (similitud de 75 a 80%, contra 50% al azar) |
"""

# ╔═╡ ba30dde8-c9eb-11f1-a562-eb0c02114edc
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `tasa_kirman(k, N, a)`, la probabilidad ``W(k \to k+1)`` de que
en un paso del modelo de Kirman el número de agentes con opinión 1 suba de
`k` a `k + 1`, usando la fórmula de la sección 3.
"""

# ╔═╡ ba30ddf2-c9eb-11f1-b6f4-e1ef068699ca
function tasa_kirman(k, N, a)
    missing  # reemplaza esta línea
end

# ╔═╡ ba30ddfc-c9eb-11f1-a2c9-1d6be8fcbf70
md"""
**2.** Escribe `satisfecho(iguales, totales, τ)`, que devuelva `true` si un
agente de Schelling con `iguales` vecinos de su mismo grupo, de un total de
`totales` vecinos ocupados, cumple el umbral `τ` (si `totales == 0`, devuelve
`true`).
"""

# ╔═╡ ba30ddfc-c9eb-11f1-8381-bdc27d4c436e
function satisfecho(iguales, totales, τ)
    missing  # reemplaza esta línea
end

# ╔═╡ ba30de06-c9eb-11f1-9460-856842467eac
md"""
### Ruta básica: interpretación

3. En el modelo de votante con ``N = 100`` agentes en la red completa, si al
   inicio sólo 15 agentes tienen la opinión 1, ¿cuál es la probabilidad de que
   toda la población termine con la opinión 1? Compara esa respuesta con lo que
   ocurriría en el modelo de Ising a temperatura muy baja (sesión 10).
4. ¿Por qué el consenso tarda mucho más en un anillo que en una red de
   Erdős-Rényi o de Barabási-Albert del mismo tamaño?
5. En el modelo de Kirman con ``a = 0.2``, la población pasa largos periodos en
   un extremo y de pronto cambia al otro. ¿Qué implicación tiene esto para
   alguien que busca una noticia económica detrás de cada cambio brusco del
   mercado?
6. En el modelo de Schelling con ``\tau = 0.375``, ningún agente rechaza ser
   minoría en su barrio (acepta hasta 62.5% de vecinos del otro grupo). Explica
   paso a paso por qué la ciudad no se queda con barrios mezclados al 50%.
7. ¿Qué pasa en el modelo de Schelling cuando el umbral ``\tau`` sube por
   encima de 0.75? ¿Por qué baja la fracción de agentes satisfechos?

### Ruta avanzada

8. Plantea la recurrencia para la probabilidad de victoria ``v_k`` de la
   opinión 1 desde el estado ``k`` en el modelo de votante en la red completa y
   demuestra que ``v_k = k / N``.
9. Verifica que la distribución Beta-binomial de la sección 5 cumple la
   ecuación de balance detallado ``P_k W(k \to k+1) = P_{k+1} W(k+1 \to k)`` y
   comprueba que para ``a = 1`` se reduce a ``P_k = 1 / (N + 1)``.
10. Usa la fórmula estacionaria de la ecuación de Fokker-Planck (repaso E.A4)
    con ``A(x) = a(1/2 - x)`` y ``B(x) = x(1-x)`` para obtener la densidad
    ``p(x) \propto x^{a-1}(1-x)^{a-1}``.
11. En una red cualquiera donde un nodo se elige al azar y copia a un vecino al
    azar, demuestra que ``\sum_i k_i s_i`` (la opinión ponderada por el grado
    ``k_i``) es una martingala. ¿Qué implica eso para la probabilidad de que
    gane una opinión que empieza en un concentrador (*hub*) de una red de
    Barabási-Albert?

## Para la próxima sesión

- Lectura: Sen y Chakrabarti (2014), capítulos 3 y 4; Kirman (2011), capítulos
  2, 3 y 6; Requião da Cunha (2021), sección 8.2. Para la ruta avanzada,
  Slanina (2014), capítulo 3 y sección 8.2; Sen y Chakrabarti (2014),
  capítulo 6; Castellano, Fortunato y Loreto (2009); Kirman (1993). De
  consulta: Wilensky y Rand (2015).
- En la semana 12 los agentes ya no copian a ciegas: eligen estrategias según
  sus pagos (**juegos evolutivos, dinámica del replicador, ciclo de Goodwin y
  cooperación en redes**).
"""

# ╔═╡ ba30de12-c9eb-11f1-90a2-0bf76f839b53
md"""
---
## Apéndice: paquetes y funciones auxiliares

Estas celdas definen herramientas que usa el cuaderno. No es necesario
leerlas para seguir la sesión.
"""

# ╔═╡ ba30de1a-c9eb-11f1-8e7e-a124f3e02890
begin
    const AZUL       = "#2a78d6"   # serie principal, grupo 1
    const NARANJA    = "#eb6834"   # segunda serie o teoría, grupo 2
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const BLANCO     = "#ffffff"
end;

# ╔═╡ ba30de24-c9eb-11f1-8ed9-a710a0ad859a
"Porcentaje con a lo más un decimal: 0.375 → \"37.5%\", 0.25 → \"25%\"."
porcentaje(x) = replace(string(round(100x, digits = 1)), r"\.0$" => "") * "%"

# ╔═╡ ba30de24-c9eb-11f1-a08d-9f54f2b0ce27
"""
Tiempo esperado exacto (en barridos) hasta el consenso en la red completa
desde k₀ agentes en 1. Para N grande se aproxima por
-(N - 1)[x ln x + (1 - x) ln(1 - x)] con x = k₀/N.
"""
function tiempo_teorico_votante(N, k0)
    (k0 == 0 || k0 == N) && return 0.0
    (N - 1) / N * ((N - k0) * sum(1 / (N - j) for j in 1:k0-1; init = 0.0) +
                   k0 * sum(1 / j for j in k0:N-1))
end

# ╔═╡ ba30de2e-c9eb-11f1-91a5-fdf7172ba36b
"""
Simula M trayectorias del modelo de votante en la red completa con N agentes,
empezando con una fracción x0 en 1. Guarda las primeras 12 trayectorias (en
barridos) y devuelve la fracción de victorias de 1 y el tiempo promedio.
"""
function trayectorias_votante(N, x0; M = 200, rng = Xoshiro(7))
    k0 = clamp(round(Int, x0 * N), 1, N - 1)
    curvas = Vector{Float64}[]
    gana1 = 0
    tiempos = Float64[]
    for m in 1:M
        k = k0
        pasos = 0
        tr = m <= 12 ? [k / N] : Float64[]
        while 0 < k < N
            pasos += 1
            # k cambia con probabilidad 2 k (N - k) / (N (N - 1)), subiendo o bajando con prob 1/2
            if rand(rng) < 2 * k * (N - k) / (N * (N - 1))
                k += rand(rng, Bool) ? 1 : -1
            end
            if m <= 12 && (pasos % N == 0 || k == 0 || k == N)
                push!(tr, k / N)
            end
        end
        m <= 12 && push!(curvas, tr)
        gana1 += (k == N)
        push!(tiempos, pasos / N)
    end
    (curvas = curvas, gana1 = gana1 / M, t_medio = mean(tiempos))
end

# ╔═╡ ba30dd66-c9eb-11f1-92ab-11cc0b323895
sim_votante = trayectorias_votante(N_vot, x0_vot; M = 200);

# ╔═╡ ba30dd6e-c9eb-11f1-b5f5-8da193774ba6
let
    d = sim_votante
    p = plot(; xlabel = "tiempo (barridos)", ylabel = "fracción con opinión 1",
             ylims = (-0.03, 1.03), legend = :outerbottom, legend_columns = 3, size = (760, 400),
             title = "Votante en red completa (N = $N_vot, x₀ = $x0_vot): 12 trayectorias")
    for tr in d.curvas
        c = tr[end] == 1.0 ? AZUL : NARANJA
        plot!(p, 0:length(tr)-1, tr; color = c, alpha = 0.45, lw = 1.1, label = "")
    end
    plot!(p, [NaN], [NaN]; color = AZUL, lw = 2, label = "Terminan en consenso 1")
    plot!(p, [NaN], [NaN]; color = NARANJA, lw = 2, label = "Terminan en consenso 0")
    hline!(p, [x0_vot]; color = GRIS, ls = :dash, lw = 1.5, label = "Fracción inicial x₀")
    p
end

# ╔═╡ ba30dd6e-c9eb-11f1-9ce8-f3fab3a8037f
let
    pemp = round(Int, 100 * sim_votante.gana1)
    pteo = round(Int, 100 * x0_vot)
    temp = round(sim_votante.t_medio, digits = 1)
    tteo = round(tiempo_teorico_votante(N_vot, round(Int, x0_vot * N_vot)), digits = 1)
    frase = x0_vot < 0.5 ? "una minoría del $pteo% gana el $pteo% de las veces" :
            x0_vot > 0.5 ? "una mayoría del $pteo% pierde el $(100 - pteo)% de las veces" :
            "con la población dividida a la mitad, cada opinión gana la mitad de las veces"
    Markdown.parse("""
De 200 simulaciones que empezaron en ``x_0 = $x0_vot``, la opinión 1 ganó en el
$pemp% de las corridas (teoría: $pteo%), y el consenso tardó en promedio
$temp barridos (teoría exacta: $tteo barridos).

A diferencia del modelo de Ising a temperatura baja (sesión 10), donde la
mayoría inicial casi siempre se impone, en el modelo de votante $frase:
nadie empuja hacia la mayoría; el consenso llega únicamente porque las
fluctuaciones aleatorias terminan topando con uno de los dos bordes
absorbentes.
""")
end

# ╔═╡ ba30de38-c9eb-11f1-9748-71256e67926f
"Simula el modelo de votante en un grafo `g` desde k0 nodos en 1 al azar; devuelve los barridos hasta el consenso."
function votante_en_red(g, k0, rng; max_barridos = 25_000)
    N = nv(g)
    s = zeros(Int, N)
    s[randperm(rng, N)[1:k0]] .= 1
    suma = k0
    adj = [neighbors(g, i) for i in 1:N]
    pasos_max = max_barridos * N
    pasos = 0
    while 0 < suma < N && pasos < pasos_max
        pasos += 1
        i = rand(rng, 1:N)
        vec = adj[i]
        isempty(vec) && continue
        j = vec[rand(rng, 1:length(vec))]
        if s[i] != s[j]
            suma += s[j] - s[i]
            s[i] = s[j]
        end
    end
    pasos / N
end

# ╔═╡ ba30de38-c9eb-11f1-ab2b-83d7d582b2e8
"Genera una red de Erdős-Rényi conexa con N nodos y grado promedio c."
function er_conexa(N, c, semilla)
    m = round(Int, N * c / 2)
    for s in semilla:semilla+1000
        g = erdos_renyi(N, m; seed = s)
        is_connected(g) && return g
    end
    erdos_renyi(N, m; seed = semilla)
end

# ╔═╡ ba30de44-c9eb-11f1-922e-bd6a199de4e5
"Compara el tiempo al consenso en función de N (red completa) y entre distintas topologías de red."
function comparar_redes_votante(N; M = 120)
    rng = Xoshiro(15)
    Ns = [20, 40, 80, 120, 160, 200]
    t_comp = [trayectorias_votante(n, 0.5; M = 200, rng = Xoshiro(100 + n)).t_medio for n in Ns]

    L = round(Int, sqrt(N))
    anillo = watts_strogatz(N, 2, 0.0; seed = 1)
    cuadricula = Graphs.grid([L, L]; periodic = true)
    er = er_conexa(N, 4, 10)
    ba = barabasi_albert(N, 2; seed = 20)

    redes = map((("Anillo 1D", anillo), ("Cuadrícula 2D", cuadricula),
                 ("Erdős-Rényi", er), ("Barabási-Albert", ba))) do (nombre, g)
        ts = [votante_en_red(g, N ÷ 2, rng) for _ in 1:M]
        (nombre = nombre, media = mean(ts), se = std(ts) / sqrt(M))
    end
    (Ns = Ns, t_completa = t_comp, redes = redes)
end

# ╔═╡ ba30dd84-c9eb-11f1-bbd9-d3efeb45dde2
tiempos_redes = comparar_redes_votante(144; M = 120);

# ╔═╡ ba30dd84-c9eb-11f1-85c5-efef78039671
let
    Ns = [20, 40, 80, 120, 160, 200]
    p1 = plot(Ns, tiempo_teorico_votante.(Ns, Ns .÷ 2); color = NARANJA, lw = 2.2,
              label = "Teoría exacta (≈ (N − 1) ln 2)", xlabel = "número de agentes N",
              ylabel = "tiempo al consenso (barridos)", legend = :topleft,
              title = "Red completa: el tiempo crece como N")
    scatter!(p1, tiempos_redes.Ns, tiempos_redes.t_completa; color = AZUL, ms = 5, msw = 0,
             label = "Simulación (x₀ = 0.5)")

    nombres = [r.nombre for r in tiempos_redes.redes]
    medias = [r.media for r in tiempos_redes.redes]
    errores = [1.96 * r.se for r in tiempos_redes.redes]
    p2 = bar(1:length(nombres), medias; yerror = errores, color = AZUL, lw = 0,
             bar_width = 0.55, xticks = (1:length(nombres), nombres),
             ylabel = "barridos hasta el consenso", legend = false,
             title = "Efecto de la topología (N = 144, x₀ = 0.5)")
    plot(p1, p2; layout = (1, 2), size = (900, 360), bottom_margin = 5Plots.mm)
end

# ╔═╡ ba30de44-c9eb-11f1-b07e-bb8bd54bdb70
"""
Simula el modelo de hormigas de Kirman con N agentes y parámetro de cambio
autónomo `a` durante `barridos` barridos (descarta 200 de calentamiento).
Devuelve el número `k` de agentes en 1 después de cada barrido.
"""
function simular_kirman(N, a; barridos = 4_000, rng = Xoshiro(11))
    k = N ÷ 2
    denom = N * (N - 1 + 2a)
    ks = Int[]
    for b in 1:(barridos + 200)
        for _ in 1:N
            u = rand(rng) * denom
            sube = (N - k) * (a + k)
            baja = k * (a + N - k)
            if u < sube
                k += 1
            elseif u < sube + baja
                k -= 1
            end
        end
        b > 200 && push!(ks, k)
    end
    (k = ks,)
end

# ╔═╡ ba30dd98-c9eb-11f1-b878-bb6a34c88b64
sim_kirman = simular_kirman(100, a_kirman; barridos = 4_000, rng = Xoshiro(11));

# ╔═╡ ba30de4c-c9eb-11f1-bc40-43b8e3c5e1e0
"Distribución estacionaria exacta P_k (k = 0, …, N) del modelo de Kirman: Beta-binomial(N, a, a)."
kirman_estacionaria(N, a) = pdf.(BetaBinomial(N, a, a), 0:N)

# ╔═╡ ba30dda2-c9eb-11f1-8619-1d310df67452
let
    N = 100
    x = sim_kirman.k ./ N
    p1 = plot(1:600, x[1:600]; color = AZUL, lw = 1.2, ylims = (0, 1),
              xlabel = "barrido", ylabel = "fracción x = k / N", legend = false,
              title = "Primeros 600 barridos (N = $N, a = $a_kirman)")
    hline!(p1, [0.5]; color = GRIS, ls = :dot)

    ks = 0:N
    emp = [ count(==(k), sim_kirman.k) / length(sim_kirman.k) for k in ks ]
    teo = kirman_estacionaria(N, a_kirman)
    p2 = bar(ks ./ N, emp .* N; bar_width = 1 / N, color = AZUL, alpha = 0.55, lw = 0,
             label = "Simulación", xlabel = "fracción x = k / N", ylabel = "densidad",
             xlims = (-0.02, 1.02), legend = :top,
             title = a_kirman < 1 ? "a < 1: bimodal (herding)" :
                     a_kirman == 1 ? "a = 1: uniforme (umbral)" : "a > 1: unimodal")
    plot!(p2, ks ./ N, teo .* N; color = NARANJA, lw = 2.5,
          label = "Solución exacta (ecuación maestra)")
    plot(p1, p2; layout = (1, 2), size = (900, 360), bottom_margin = 5Plots.mm)
end

# ╔═╡ ba30dde0-c9eb-11f1-b517-a5b08cdc14b1
let
    N = 80
    ks = 1:N-1
    xs = range(0.01, 0.99, length = 300)
    p = plot(; xlabel = "fracción x = k / N", ylabel = "densidad p(x)",
             ylims = (0, 3.2), xlims = (0, 1), legend = :top,
             title = "Ecuación maestra discreta (puntos) vs. Fokker-Planck Beta(a, a) (curvas)")
    for (a, c, lb) in ((0.4, NARANJA, "a = 0.4 (bimodal)"),
                       (1.0, GRIS, "a = 1.0 (uniforme)"),
                       (3.0, AZUL, "a = 3.0 (unimodal)"))
        plot!(p, xs, pdf.(Beta(a, a), xs); color = c, lw = 2.2, label = lb)
        P = kirman_estacionaria(N, a)
        scatter!(p, ks ./ N, N .* P[2:end-1]; color = c, ms = 3, msw = 0, label = "")
    end
    p
end

# ╔═╡ ba30de4c-c9eb-11f1-aaf4-3fb0f2adde31
"Cuenta vecinos del mismo tipo y vecinos ocupados totales en la vecindad de Moore (8 vecinos, periódica)."
function vecindad_schelling(tablero, r, c)
    L = size(tablero, 1)
    tipo = tablero[r, c]
    iguales, totales = 0, 0
    for dr in -1:1, dc in -1:1
        (dr == 0 && dc == 0) && continue
        v = tablero[mod1(r + dr, L), mod1(c + dc, L)]
        if v != 0
            totales += 1
            v == tipo && (iguales += 1)
        end
    end
    iguales, totales
end

# ╔═╡ ba30de56-c9eb-11f1-aab1-81d30399b603
"Mide la similitud promedio y la fracción de agentes satisfechos en el tablero de Schelling."
function metricas_schelling(tablero, τ)
    L = size(tablero, 1)
    suma_sim, con_vecinos, felices, agentes = 0.0, 0, 0, 0
    for r in 1:L, c in 1:L
        tablero[r, c] == 0 && continue
        agentes += 1
        ig, tot = vecindad_schelling(tablero, r, c)
        if tot > 0
            suma_sim += ig / tot
            con_vecinos += 1
            ig / tot >= τ && (felices += 1)
        else
            felices += 1
        end
    end
    (similitud = con_vecinos > 0 ? suma_sim / con_vecinos : 1.0,
     felices = felices / max(agentes, 1))
end

# ╔═╡ ba30de56-c9eb-11f1-8921-27917a7fd6fc
"""
Simula el modelo de segregación de Schelling en una cuadrícula L × L con
fracción `p_vacio` de celdas vacías y umbral de tolerancia `τ`.
"""
function simular_schelling(L, p_vacio, τ; pasos_max = 120, rng = Xoshiro(21))
    n_celdas = L * L
    n_vacios = round(Int, n_celdas * p_vacio)
    n_agentes = n_celdas - n_vacios
    n_uno = n_agentes ÷ 2
    celdas = vcat(zeros(Int, n_vacios), ones(Int, n_uno), fill(2, n_agentes - n_uno))
    shuffle!(rng, celdas)
    tablero = reshape(celdas, L, L)

    m0 = metricas_schelling(tablero, τ)
    sim_hist = [m0.similitud]
    fel_hist = [m0.felices]

    for _ in 1:pasos_max
        infelices = Tuple{Int,Int}[]
        vacias = Tuple{Int,Int}[]
        for r in 1:L, c in 1:L
            if tablero[r, c] == 0
                push!(vacias, (r, c))
            else
                ig, tot = vecindad_schelling(tablero, r, c)
                (tot > 0 && ig / tot < τ) && push!(infelices, (r, c))
            end
        end
        isempty(infelices) && break
        shuffle!(rng, infelices)
        shuffle!(rng, vacias)
        for (r, c) in infelices
            isempty(vacias) && break
            vr, vc = popfirst!(vacias)
            tablero[vr, vc] = tablero[r, c]
            tablero[r, c] = 0
            push!(vacias, (r, c))
        end
        mt = metricas_schelling(tablero, τ)
        push!(sim_hist, mt.similitud)
        push!(fel_hist, mt.felices)
    end
    (tablero = tablero, similitud = sim_hist, felices = fel_hist)
end

# ╔═╡ ba30ddac-c9eb-11f1-aeb3-77a493e4aa4d
sim_schelling = simular_schelling(40, vac_sch, τ_sch; pasos_max = 120, rng = Xoshiro(21));

# ╔═╡ ba30ddb6-c9eb-11f1-bfe3-c18094e9e328
let
    s = sim_schelling
    paleta = cgrad([BLANCO, AZUL, NARANJA])
    p1 = heatmap(s.tablero; color = paleta, clims = (0, 2), aspect_ratio = :equal,
                 colorbar = false, axis = false, ticks = false, grid = false,
                 title = "Ciudad final (τ = $(porcentaje(τ_sch))): similitud = $(round(Int, 100s.similitud[end]))%")
    rondas = 0:length(s.similitud)-1
    p2 = plot(rondas, 100 .* s.similitud; color = AZUL, lw = 2.2, ylims = (0, 105),
              xlabel = "ronda", ylabel = "porcentaje (%)", legend = :bottomright,
              label = "Similitud promedio de los vecindarios",
              title = "Evolución en $(length(rondas)-1) rondas")
    plot!(p2, rondas, 100 .* s.felices; color = NARANJA, lw = 2, ls = :dash,
          label = "Agentes satisfechos")
    hline!(p2, [100τ_sch]; color = GRIS, ls = :dot, label = "Umbral pedido τ = $(porcentaje(τ_sch))")
    plot(p1, p2; layout = Plots.grid(1, 2; widths = [0.44, 0.56]), size = (900, 380))
end

# ╔═╡ ba30de60-c9eb-11f1-8abb-a5b49378d25d
"Ejemplo fijo de Schelling para la sección 0 (antes y después con τ = 0.375)."
ejemplo_schelling = let
    ini = simular_schelling(40, 0.15, 0.375; pasos_max = 0, rng = Xoshiro(21))
    fin = simular_schelling(40, 0.15, 0.375; pasos_max = 120, rng = Xoshiro(21))
    (ini = (tablero = ini.tablero, similitud = ini.similitud[1], felices = ini.felices[1]),
     fin = (tablero = fin.tablero, similitud = fin.similitud[end], felices = fin.felices[end]))
end;

# ╔═╡ ba30dd3c-c9eb-11f1-b1b0-7b5be74895ab
let
    ini, fin = ejemplo_schelling
    paleta = cgrad([BLANCO, AZUL, NARANJA])
    p1 = heatmap(ini.tablero; color = paleta, clims = (0, 2), aspect_ratio = :equal,
                 colorbar = false, axis = false, ticks = false, grid = false,
                 title = "Inicio al azar (similitud = $(round(Int, 100ini.similitud))%)")
    p2 = heatmap(fin.tablero; color = paleta, clims = (0, 2), aspect_ratio = :equal,
                 colorbar = false, axis = false, ticks = false, grid = false,
                 title = "Equilibrio con τ = 37.5% (similitud = $(round(Int, 100fin.similitud))%)")
    plot(p1, p2; layout = (1, 2), size = (820, 390))
end

# ╔═╡ ba30dd3c-c9eb-11f1-a8e8-7bd7b08b20e2
md"""
Nadie buscaba una ciudad con $(round(Int, 100ejemplo_schelling.fin.similitud))% de vecinos del mismo grupo: bastaba con
37.5%. El exceso de segregación es una propiedad **emergente**. En esta sesión
estudiaremos tres mecanismos fundamentales de dinámica social: el **modelo de
votante**, el **modelo de Kirman** y el **modelo de Schelling**.
"""

# ╔═╡ ba30de6a-c9eb-11f1-ab55-5ddac6518317
"Barrido de la similitud final y fracción de satisfechos en Schelling para varios umbrales τ."
function barrido_schelling(L, p_vacio; τs = 0.10:0.05:0.85)
    sim, fel = Float64[], Float64[]
    for (i, τ) in enumerate(τs)
        res = simular_schelling(L, p_vacio, τ; pasos_max = 100, rng = Xoshiro(30 + i))
        push!(sim, res.similitud[end])
        push!(fel, res.felices[end])
    end
    (τ = collect(τs), similitud = sim, felices = fel)
end

# ╔═╡ ba30ddc0-c9eb-11f1-8aee-8b4ddc4557ad
barrido_sch = barrido_schelling(35, 0.15; τs = 0.10:0.05:0.85);

# ╔═╡ ba30ddc0-c9eb-11f1-bcd1-9365670f99ba
let
    b = barrido_sch
    p = plot(100 .* b.τ, 100 .* b.similitud; color = AZUL, lw = 2.5, marker = :circle, ms = 4, msw = 0,
             label = "Similitud promedio final", xlabel = "umbral requerido τ (%)",
             ylabel = "porcentaje (%)", ylims = (0, 105), legend = :bottomleft,
             title = "Segregación emergente frente al umbral individual")
    plot!(p, 100 .* b.τ, 100 .* b.felices; color = NARANJA, lw = 2, ls = :dash, marker = :square, ms = 3.5, msw = 0,
          label = "Agentes satisfechos al final")
    plot!(p, 100 .* b.τ, 100 .* b.τ; color = GRIS, lw = 1, ls = :dot, label = "Similitud = umbral τ")
    hline!(p, [50]; color = GRIS_CLARO, lw = 1, label = "Mezcla al azar (50%)")
end

# ╔═╡ ba30de6a-c9eb-11f1-8fb9-17405055c55e
"Rango de la similitud final del barrido de Schelling para τ entre `lo` y `hi`, como texto."
function rango_sch(lo, hi)
    k = (barrido_sch.τ .>= lo - 1e-9) .& (barrido_sch.τ .<= hi + 1e-9)
    a, b = extrema(barrido_sch.similitud[k])
    "$(round(Int, 100a))% y $(round(Int, 100b))%"
end

# ╔═╡ ba30ddca-c9eb-11f1-bac4-ebc2f299cc04
md"""
La curva muestra tres regímenes muy claros:

1. **``\tau \le 25\%`` (mezcla estable):** Casi todos los agentes ya cumplen su
   umbral en la configuración aleatoria inicial; casi nadie se muda y la
   similitud queda entre $(rango_sch(0.0, 0.25)), cerca de la mezcla al azar.
2. **``30\% \le \tau \le 65\%`` (segregación fuerte con preferencias leves):**
   En cuanto ``\tau`` cruza ``\approx 1/3``, las mudanzas de los pocos
   insatisfechos iniciales reducen la fracción de su grupo en el barrio que
   dejan y la aumentan en el que llegan. Eso vuelve insatisfechos a otros
   vecinos y desata un efecto dominó que sólo se detiene cuando la similitud
   llega a valores entre $(rango_sch(0.3, 0.65)), muy por encima de la línea
   punteada ``\tau`` que cada agente pedía.
3. **``\tau \ge 75\%`` (frustración / congelamiento):** Los agentes exigen
   vecindarios tan puros que no hay suficientes lugares satisfactorios; las
   mudanzas no logran acomodar a todos y la fracción de satisfechos cae.

**Una analogía con el modelo de Ising.** Si asignamos ``s_i = +1`` a un
grupo, ``s_i = -1`` al otro y ``s_i = 0`` a las celdas vacías, la energía de
Ising ``H = -J \sum_{\langle i,j\rangle} s_i s_j`` baja cuando hay más pares
de vecinos del mismo grupo y menos fronteras entre grupos. Las mudanzas de
Schelling suelen reducir esas fronteras, conservando el número de agentes de
cada grupo, y la ciudad se separa en dominios parecidos a los del Ising debajo
de ``T_c`` (sesión 10). Es una **analogía**, no una equivalencia: un agente se
muda para cumplir su propio umbral, no para bajar ``H``, y una mudanza puede
aumentar la energía total (Vinković y Kirman, 2006).
"""

# ╔═╡ ba30de72-c9eb-11f1-b6a5-97d4fcf62c13
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

# ╔═╡ ba30ddf2-c9eb-11f1-be88-318a1442eaac
verificar(() -> [tasa_kirman(0, 100, 1.0), tasa_kirman(50, 100, 1.0), tasa_kirman(100, 100, 1.0)],
          [1.0 / 101.0, 0.5 * 51.0 / 101.0, 0.0];
          pista = md"La fórmula es ``\frac{N - k}{N}\,\frac{a + k}{N - 1 + 2a}``.")

# ╔═╡ ba30ddfc-c9eb-11f1-a70f-49fb7f691238
verificar(() -> [satisfecho(3, 8, 0.375), satisfecho(2, 8, 0.375), satisfecho(0, 0, 0.5)],
          [true, false, true];
          pista = md"Si `totales == 0` devuelve `true`; en caso contrario compara `iguales / totales >= τ`.")

# ╔═╡ Cell order:
# ╟─ba30dcf8-c9eb-11f1-9110-e936cab76302
# ╟─ba30dd2a-c9eb-11f1-a955-f795c27e6017
# ╟─ba30dd34-c9eb-11f1-9106-7533c1d56ae3
# ╟─ba30dd3c-c9eb-11f1-b1b0-7b5be74895ab
# ╟─ba30dd3c-c9eb-11f1-a8e8-7bd7b08b20e2
# ╟─ba30dd52-c9eb-11f1-b7ce-61b65d44eeaf
# ╠═ba30dd52-c9eb-11f1-bc01-db539696fd6d
# ╠═ba30dd5c-c9eb-11f1-9249-51a366d88393
# ╟─ba30dd66-c9eb-11f1-8675-f30b9c80414d
# ╟─ba30dd66-c9eb-11f1-92ab-11cc0b323895
# ╟─ba30dd6e-c9eb-11f1-b5f5-8da193774ba6
# ╟─ba30dd6e-c9eb-11f1-9ce8-f3fab3a8037f
# ╟─ba30dd7a-c9eb-11f1-9d97-59567d04119d
# ╠═ba30dd84-c9eb-11f1-bbd9-d3efeb45dde2
# ╟─ba30dd84-c9eb-11f1-85c5-efef78039671
# ╟─ba30dd8e-c9eb-11f1-8b3d-b36fc8f1d91f
# ╟─ba30dd8e-c9eb-11f1-8021-b31c7b64e531
# ╠═ba30dd98-c9eb-11f1-b878-bb6a34c88b64
# ╟─ba30dda2-c9eb-11f1-8619-1d310df67452
# ╟─ba30dda2-c9eb-11f1-ba7f-89022dd9cd47
# ╟─ba30ddac-c9eb-11f1-934f-4f29b55971fc
# ╠═ba30ddac-c9eb-11f1-aeb3-77a493e4aa4d
# ╟─ba30ddb6-c9eb-11f1-bfe3-c18094e9e328
# ╟─ba30ddb6-c9eb-11f1-97c8-ff604b3e1731
# ╠═ba30ddc0-c9eb-11f1-8aee-8b4ddc4557ad
# ╟─ba30ddc0-c9eb-11f1-bcd1-9365670f99ba
# ╟─ba30ddca-c9eb-11f1-bac4-ebc2f299cc04
# ╟─ba30ddd4-c9eb-11f1-963a-552a0eb453a0
# ╟─ba30dde0-c9eb-11f1-b517-a5b08cdc14b1
# ╟─ba30dde8-c9eb-11f1-a088-2f2baabc45b9
# ╟─ba30dde8-c9eb-11f1-a562-eb0c02114edc
# ╠═ba30ddf2-c9eb-11f1-b6f4-e1ef068699ca
# ╟─ba30ddf2-c9eb-11f1-be88-318a1442eaac
# ╟─ba30ddfc-c9eb-11f1-a2c9-1d6be8fcbf70
# ╠═ba30ddfc-c9eb-11f1-8381-bdc27d4c436e
# ╟─ba30ddfc-c9eb-11f1-a70f-49fb7f691238
# ╟─ba30de06-c9eb-11f1-9460-856842467eac
# ╟─ba30de12-c9eb-11f1-90a2-0bf76f839b53
# ╠═ba30de1a-c9eb-11f1-b180-47208819aec9
# ╠═ba30de1a-c9eb-11f1-8e7e-a124f3e02890
# ╠═ba30de24-c9eb-11f1-8ed9-a710a0ad859a
# ╠═ba30de24-c9eb-11f1-a08d-9f54f2b0ce27
# ╠═ba30de2e-c9eb-11f1-91a5-fdf7172ba36b
# ╠═ba30de38-c9eb-11f1-9748-71256e67926f
# ╠═ba30de38-c9eb-11f1-ab2b-83d7d582b2e8
# ╠═ba30de44-c9eb-11f1-922e-bd6a199de4e5
# ╠═ba30de44-c9eb-11f1-b07e-bb8bd54bdb70
# ╠═ba30de4c-c9eb-11f1-bc40-43b8e3c5e1e0
# ╠═ba30de4c-c9eb-11f1-aaf4-3fb0f2adde31
# ╠═ba30de56-c9eb-11f1-aab1-81d30399b603
# ╠═ba30de56-c9eb-11f1-8921-27917a7fd6fc
# ╠═ba30de60-c9eb-11f1-8abb-a5b49378d25d
# ╠═ba30de6a-c9eb-11f1-8fb9-17405055c55e
# ╠═ba30de6a-c9eb-11f1-ab55-5ddac6518317
# ╠═ba30de72-c9eb-11f1-b6a5-97d4fcf62c13
