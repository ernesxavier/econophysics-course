#%% hide
md"""
# Sesión 15: Votante y Kirman: imitación y herding

**Objetivos de la sesión**
1. Construir **modelos basados en agentes (ABM)** en los que reglas locales de
   imitación generan patrones macroscópicos sin coordinación central.
2. Distinguir entre cadenas de Markov con **estados absorbentes** (imitación
   pura) y con **distribución estacionaria** (imitación más ruido).
3. Simular el **modelo de votante**, comprobar que la opinión media es una
   martingala y comparar el **tiempo de consenso** en distintas redes.
4. Incorporar cambios autónomos de opinión en el **modelo de hormigas de
   Kirman** para explicar el **herding** y la transición hacia una distribución
   **bimodal**.
5. **[Avanzado]** Contrastar la simulación con el tiempo exacto de absorción del
   votante y con la solución estacionaria de la ecuación maestra de Kirman.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s15_votante_kirman.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): E.B2–B3 (caminata aleatoria
y promedios), L.B4 (matriz de adyacencia), C.B4 (función logística), F.B3–B4
(temperatura y transiciones de fase), E.I2 (martingalas), E.I3 (cadenas de
Markov), F.I3–I4 (Ising y campo medio). Ruta avanzada: E.A2 (lema de Itô),
E.A3 (ecuación maestra), E.A4 (Fokker-Planck), E.A6 (tiempos de primer paso y
absorción).
"""

#%% hide
TableOfContents(title = "Contenido")

#%% hide
md"""
## 0. El fenómeno: dos restaurantes idénticos; uno se llena y el otro no

En el Módulo II estudiamos sistemas en equilibrio estadístico, como el modelo
de Ising. Hoy abrimos el **Módulo III** con modelos donde los agentes siguen
reglas directas de comportamiento —copiar a un vecino o seguir a la multitud—
sin resolver un problema de optimización global.

Becker (1991) planteó un acertijo cotidiano: dos restaurantes con la misma
comida y los mismos precios están uno frente al otro. En vez de repartirse la
clientela mitad y mitad, uno se llena y el otro queda casi vacío; tiempo
después, sin que cambie la calidad de ninguno, la multitud se pasa al otro.
Kirman (1993) observó que los entomólogos habían visto exactamente lo mismo con
hormigas frente a dos fuentes idénticas de alimento: cuando cada individuo
recluta a otros con frecuencia y rara vez explora por su cuenta, la mayoría se
concentra en una sola opción y alterna en olas espontáneas.
"""

#%% hide
let
    N = 100
    c_herd = simular_kirman(N, 0.2; barridos = 800, rng = Xoshiro(11)).k ./ N
    c_ind  = simular_kirman(N, 4.0; barridos = 800, rng = Xoshiro(11)).k ./ N
    p1 = plot(1:800, c_herd; color = AZUL, lw = 1.4, ylims = (0, 1),
              xlabel = "tiempo (barridos)", ylabel = "fracción en el restaurante A",
              label = "Con herding fuerte (a = 0.2)", legend = :outerbottom,
              title = "Dos restaurantes idénticos (100 clientes)")
    plot!(p1, 1:800, c_ind; color = NARANJA, lw = 1.2, alpha = 0.75,
          label = "Casi independientes (a = 4.0)")
    hline!(p1, [0.5]; color = GRIS, ls = :dot, label = "")

    ks = 0:N
    p2 = plot(ks ./ N, kirman_estacionaria(N, 0.2) .* N; color = AZUL, lw = 2.5,
              xlabel = "fracción en el restaurante A", ylabel = "densidad estacionaria",
              label = "a = 0.2: uno lleno y el otro vacío", legend = :top,
              title = "Distribución del reparto de clientes")
    plot!(p2, ks ./ N, kirman_estacionaria(N, 4.0) .* N; color = NARANJA, lw = 2.5,
          label = "a = 4.0: mitad y mitad")
    plot(p1, p2; layout = (1, 2), size = (880, 350), bottom_margin = 4Plots.mm)
end

#%% hide
md"""
Ninguno de los dos restaurantes es mejor, ni ocurrió ningún choque externo
cuando la multitud cambió de restaurante: las olas nacen de la **imitación** combinada con
unos pocos cambios individuales. En esta sesión empezaremos con la imitación
pura (**modelo de votante**) y luego añadiremos esos cambios autónomos
(**modelo de Kirman**).
"""

#%% hide
md"""
## 1. Repaso: estados absorbentes frente a ruido

En una cadena de Markov donde los agentes **únicamente copian** las opiniones
existentes, si en algún momento todos coinciden (``k = 0`` o ``k = N`` agentes
con la opinión 1), ya no queda nadie con la opinión contraria a quien copiar.
Esos dos consensos son **estados absorbentes** (repaso E.I3 y E.A6): una vez que
el sistema entra en uno, se queda ahí para siempre, y la pregunta relevante es
**quién gana** y **cuánto tarda en absorberse**.

Si además cada agente tiene una pequeña probabilidad de cambiar de opinión por
cuenta propia (**ruido idiosincrático**, repaso F.B3), los extremos dejan de
ser absorbentes: siempre puede reaparecer un disidente. El sistema ya no se
detiene y converge a una **distribución estacionaria** ``P_k``.
"""

#%%
# En imitación pura, la probabilidad de que un agente en 0 copie a uno en 1
# cuando k de N tienen opinión 1 es proporcional a k(N - k): se anula en k = 0 y en k = N.
tasa_imitacion_pura(k, N) = (k / N) * ((N - k) / (N - 1))

#%%
tasa_imitacion_pura(0, 100), tasa_imitacion_pura(50, 100), tasa_imitacion_pura(100, 100)

#%% hide
md"""
## 2. El modelo de votante

El **modelo de votante** (Clifford y Sudbury, 1973; Holley y Liggett, 1975) es
el modelo más simple de formación de opinión por imitación pura.

**Reglas.**
1. Hay ``N`` agentes en los nodos de una red conexa (repaso L.B4); cada uno
   tiene opinión ``s_i \in \{0, 1\}``.
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
  probabilidad de terminar en el consenso ``k = N`` es exactamente ``x_0``
  (lo comprobamos promediando sobre muchas corridas; repaso E.B3).

Fracción inicial ``x_0``: $(@bind x0_vot Slider(0.1:0.1:0.9, default = 0.3, show_value = true))
  Número de agentes ``N``: $(@bind N_vot Slider([40, 80, 150], default = 80, show_value = true))
"""

#%% hide
sim_votante = trayectorias_votante(N_vot, x0_vot; M = 200);

#%% hide
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

#%% hide
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

A diferencia del modelo de Ising a temperatura baja (sesión 13, repaso F.I3), donde la
mayoría inicial casi siempre se impone, en el modelo de votante $frase:
nadie empuja hacia la mayoría; el consenso llega únicamente porque las
fluctuaciones aleatorias terminan topando con uno de los dos bordes
absorbentes.
""")
end

#%% hide
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

#%%
tiempos_redes = comparar_redes_votante(144; M = 120);

#%% hide
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
    p2 = bar(1:length(nombres), medias; color = AZUL, lw = 0,
             bar_width = 0.55, xticks = (1:length(nombres), nombres),
             ylabel = "barridos hasta el consenso", legend = false,
             ylims = (0, 1.08 * maximum(medias .+ errores)),
             title = "Efecto de la topología (N = 144, IC 95%)")
    scatter!(p2, 1:length(nombres), medias; yerror = errores, ms = 0, msc = GRIS,
             msw = 1.5, color = GRIS, label = "")
    plot(p1, p2; layout = (1, 2), size = (900, 360), bottom_margin = 5Plots.mm)
end

#%% hide
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

#%% hide
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

#%%
sim_kirman = simular_kirman(100, a_kirman; barridos = 4_000, rng = Xoshiro(11));

#%% hide
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

#%% hide
md"""
Observa las tres fases del modelo al variar ``a`` (repaso F.B4):

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

**Diferencia con el modelo de Ising (sesión 14, repaso F.I4).** En el Ising,
cada agente responde de forma no lineal (logística, repaso C.B4) al promedio;
en Kirman, la tasa de imitación es **lineal** en ``k/N``. Por eso en Kirman la
tendencia promedio siempre apunta hacia ``x = 1/2``: la bimodalidad cuando
``a < 1`` proviene de que el **ruido de imitación se apaga cerca de los bordes**
``x = 0`` y ``x = 1`` (donde casi todos opinan igual), atrapando
temporalmente al sistema en consensos cuasi-absorbentes.
"""

#%% hide
md"""
## 4. [Avanzado] Ecuación maestra y límite continuo de Kirman

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
continua ``x = k / N \in (0, 1)`` cumple la ecuación de Fokker-Planck de
``dx = A(x)\,dt + B(x)\,dW`` (repaso E.A2 y E.A4) con deriva
``A(x) = a(1/2 - x)`` e intensidad de ruido ``B(x)^2 = x(1-x)``. Su densidad
estacionaria es la **distribución Beta** simétrica:

```math
p(x) = \frac{\Gamma(2a)}{\Gamma(a)^2}\,x^{a-1}(1 - x)^{a-1}.
```

La gráfica compara la distribución estacionaria exacta ``N \cdot P_k`` con la
densidad continua ``\mathrm{Beta}(a, a)`` para los tres regímenes ``a = 0.4``,
``a = 1`` y ``a = 3``.
"""

#%% hide
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

#%% hide
md"""
## Resumen

| Modelo | Regla microscópica | Comportamiento colectivo |
|---|---|---|
| **Votante** | Copiar a un vecino al azar (sin ruido) | Consenso absorbente; ``\mathbb{P}(\text{gana } 1) = x_0``; ``\mathbb{E}[\tau] \propto N^2`` en 1D y ``\propto N`` en red completa |
| **Kirman** | Imitación (herding) + cambio autónomo ``a`` | Sin estados absorbentes; estacionaria ``\mathrm{Beta}(a, a)``: bimodal si ``a < 1``, uniforme si ``a = 1``, unimodal si ``a > 1`` |
"""

#%% hide
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `prob_cambio_votante(k, N)`, la probabilidad total
``\mathbb{P}(k \to k+1) + \mathbb{P}(k \to k-1) = \frac{2 k (N - k)}{N (N - 1)}``
de que en un paso elemental del modelo de votante en la red completa cambie el
número de agentes con opinión 1.
"""

#%%
function prob_cambio_votante(k, N)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> [prob_cambio_votante(0, 50), prob_cambio_votante(25, 50), prob_cambio_votante(50, 50)],
          [0.0, 2 * 25 * 25 / (50 * 49), 0.0];
          pista = md"Suma las dos probabilidades de transición: `2 * k * (N - k) / (N * (N - 1))`.")

#%% hide
md"""
**2.** Escribe `tasa_kirman(k, N, a)`, la probabilidad ``W(k \to k+1)`` de que
en un paso del modelo de Kirman el número de agentes con opinión 1 suba de
`k` a `k + 1`, usando la fórmula de la sección 3.
"""

#%%
function tasa_kirman(k, N, a)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> [tasa_kirman(0, 100, 1.0), tasa_kirman(50, 100, 1.0), tasa_kirman(100, 100, 1.0)],
          [1.0 / 101.0, 0.5 * 51.0 / 101.0, 0.0];
          pista = md"La fórmula es ``\frac{N - k}{N}\,\frac{a + k}{N - 1 + 2a}``.")

#%% hide
md"""
### Ruta básica: interpretación

3. En el modelo de votante con ``N = 100`` agentes en la red completa, si al
   inicio sólo 15 agentes tienen la opinión 1, ¿cuál es la probabilidad de que
   toda la población termine con la opinión 1? Compara esa respuesta con lo que
   ocurriría en el modelo de Ising a temperatura muy baja (sesión 13).
4. ¿Por qué el consenso tarda mucho más en un anillo que en una red de
   Erdős-Rényi o de Barabási-Albert del mismo tamaño?
5. En el modelo de Kirman con ``a = 0.2``, la población pasa largos periodos en
   un extremo y de pronto cambia al otro. ¿Qué implicación tiene esto para
   alguien que busca una noticia económica detrás de cada cambio brusco del
   mercado?

### Ruta avanzada

6. Plantea la recurrencia para la probabilidad de victoria ``v_k`` de la
   opinión 1 desde el estado ``k`` en el modelo de votante en la red completa y
   demuestra que ``v_k = k / N``.
7. Verifica que la distribución Beta-binomial de la sección 4 cumple la
   ecuación de balance detallado ``P_k W(k \to k+1) = P_{k+1} W(k+1 \to k)`` y
   comprueba que para ``a = 1`` se reduce a ``P_k = 1 / (N + 1)``.
8. En una red cualquiera donde un nodo se elige al azar y copia a un vecino al
   azar, demuestra que ``\sum_i k_i s_i`` (la opinión ponderada por el grado
   ``k_i``) es una martingala. ¿Qué implica eso para la probabilidad de que
   gane una opinión que empieza en un concentrador (*hub*) de una red de
   Barabási-Albert?

## Para la próxima sesión

- Lectura: Sen y Chakrabarti (2014), capítulo 3. Complementaria: Kirman (2011),
  capítulos 2 y 3. Para la ruta avanzada, Slanina (2014), capítulo 3;
  Castellano, Fortunato y Loreto (2009); Kirman (1993).
- En la **sesión 16** (`Code/s16_schelling.jl`) estudiaremos qué ocurre cuando
  los agentes no cambian de estado sino de **lugar**: el modelo de segregación
  residencial de Schelling.
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
    using PlutoUI, Statistics, Random, Distributions, Graphs, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

#%%
begin
    const AZUL       = "#2a78d6"   # serie principal, grupo 1
    const NARANJA    = "#eb6834"   # segunda serie o teoría, grupo 2
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const BLANCO     = "#ffffff"
end;

#%%
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

#%%
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

#%%
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

#%%
"Genera una red de Erdős-Rényi conexa con N nodos y grado promedio c."
function er_conexa(N, c, semilla)
    m = round(Int, N * c / 2)
    for s in semilla:semilla+1000
        g = erdos_renyi(N, m; seed = s)
        is_connected(g) && return g
    end
    erdos_renyi(N, m; seed = semilla)
end

#%%
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

#%%
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

#%%
"Distribución estacionaria exacta P_k (k = 0, …, N) del modelo de Kirman: Beta-binomial(N, a, a)."
kirman_estacionaria(N, a) = pdf.(BetaBinomial(N, a, a), 0:N)

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
