#%% hide
md"""
# Sesiones 13 y 14: Ising y elección discreta

**Objetivos de la sesión**
1. Repasar la **función logística** y el modelo **logit** de elección
   discreta.
2. Ver que el logit es la **distribución de Boltzmann**: el modelo que usan
   los economistas para estimar demanda es física estadística.
3. Agregar **interacciones sociales** a la elección (Brock y Durlauf) y
   obtener el **modelo de Ising**.
4. Simular el Ising con el algoritmo de **Metropolis** y encontrar su
   transición de fase.
5. Resolver el **campo medio** gráficamente y compararlo con la simulación.
6. Observar la **histéresis**: cambios graduales de incentivos que producen
   saltos de comportamiento que no se revierten.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s13_logit_ising.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): C.B4 (función logística),
F.B3 (temperatura), F.B4 (transiciones de fase), F.I1 (Boltzmann), F.I2
(logit = Boltzmann), F.I3 (Ising), F.I4 (campo medio), F.I5 (Metropolis),
C.I5 (puntos fijos). Ruta avanzada: F.A1 (Landau), F.A3 (campo aleatorio y
avalanchas).
"""

#%% hide
TableOfContents(title = "Contenido")

#%% hide
md"""
## 0. El fenómeno: mismos incentivos, resultados distintos

Glaeser, Sacerdote y Scheinkman (1996) notaron que las tasas de
criminalidad varían entre ciudades de Estados Unidos mucho más de lo que
explican sus diferencias en ingreso, educación o policía. Algo parecido pasa
con la deserción escolar, la fertilidad o la adopción de tecnologías: grupos
con incentivos casi iguales terminan comportándose de forma muy distinta.
Su explicación: las **interacciones sociales**. Cuando lo que hacen mis
vecinos cambia lo que me conviene hacer a mí, pequeñas diferencias
iniciales se amplifican.

Abajo hay dos "ciudades" de 60 × 60 habitantes con **exactamente las mismas
reglas**: cada uno elige entre dos opciones (azul o naranja), prefiere
coincidir con sus cuatro vecinos y tiene algo de ruido en sus decisiones.
Ninguna opción es mejor que la otra. Sólo cambia el azar.
"""

#%% hide
let
    ps = map((3, 8)) do semilla
        s = ising_red(60, 2.0; barridos = 400, rng = Xoshiro(semilla), inicio = :azar).s
        heatmap(s; color = cgrad([NARANJA, AZUL]), clims = (-1, 1), aspect_ratio = :equal,
                colorbar = false, axis = false, ticks = false, grid = false,
                title = "Ciudad $(semilla == 3 ? "A" : "B"): $(round(Int, 100 * mean(s .== 1)))% elige azul")
    end
    plot(ps...; layout = (1, 2), size = (760, 400))
end

#%% hide
md"""
Las dos ciudades empezaron con elecciones al azar, mitad y mitad. Al final,
en cada una domina una opción distinta, y el **consenso** no lo impuso nadie:
emergió de la imitación local. Hoy vamos a construir este modelo desde la
teoría de la elección que usan los economistas.
"""

#%% hide
md"""
## 1. Repaso: la función logística y el modelo logit

La **función logística** convierte cualquier número en una probabilidad
(repaso C.B4):

```math
\sigma(x) = \frac{1}{1 + e^{-x}} .
```

En el modelo **logit**, un agente elige entre A y B, cuyas utilidades
difieren en ``\Delta u = u_A - u_B``. La probabilidad de elegir A es

```math
\mathbb{P}(A) = \sigma(\beta\,\Delta u) = \frac{e^{\beta u_A}}{e^{\beta u_A} + e^{\beta u_B}} .
```

El parámetro ``\beta`` es la **intensidad de elección**. Con ``\beta`` grande,
el agente casi siempre elige la mejor opción; con ``\beta = 0``, elige al azar.

Intensidad de elección ``\beta``: $(@bind β_logit Slider([0.25, 0.5, 1.0, 2.0, 5.0, 20.0], default = 1.0, show_value = true))
"""

#%% hide
let
    x = range(-5, 5, length = 300)
    p = plot(x, σ.(β_logit .* x); color = AZUL, lw = 2.5, ylims = (0, 1),
             xlabel = "ventaja de A en utilidad, Δu = u_A − u_B", ylabel = "probabilidad de elegir A",
             label = "β = $β_logit", legend = :bottomright, title = "Modelo logit")
    for b in (0.25, 1.0, 5.0)
        b == β_logit && continue
        plot!(p, x, σ.(b .* x); color = GRIS_CLARO, lw = 1.5, label = "β = $b")
    end
    hline!(p, [0.5]; color = GRIS, ls = :dot, label = "")
end

#%% hide
md"""
## 2. Logit = Boltzmann

¿De dónde sale la fórmula logit? McFadden la derivó suponiendo que la
utilidad tiene una parte que vemos, ``u_i``, y un **choque aleatorio**
``\varepsilon_i`` que no vemos (gustos, humor, información privada):

```math
U_i = u_i + \varepsilon_i .
```

El agente elige la opción con mayor ``U_i``. Si los ``\varepsilon_i`` son
independientes y tienen distribución de **Gumbel** con escala ``1/\beta``,
la probabilidad de que gane la opción ``i`` es exactamente

```math
\mathbb{P}(i) = \frac{e^{\beta u_i}}{\sum_j e^{\beta u_j}} .
```

Ahora compárala con la **distribución de Boltzmann** de la sesión 8
(repaso F.I1): un sistema a temperatura ``T = 1/\beta`` ocupa el estado ``i``,
con energía ``E_i``, con probabilidad ``e^{-\beta E_i}/Z``. Son **la misma
fórmula** con ``E_i = -u_i``: la utilidad es menos la energía, y la
"temperatura" ``1/\beta`` es el tamaño del ruido en las decisiones (repaso
F.I2 y F.B3).

Comprobémoslo simulando: tres opciones con utilidades 1, 0.5 y 0, y 100,000
agentes con choques de Gumbel.
"""

#%%
utilidades = [1.0, 0.5, 0.0];

#%%
elecciones = let
    rng = Xoshiro(1)
    β = 2.0
    [argmax(utilidades .+ rand(rng, Gumbel(0, 1 / β), 3)) for _ in 1:100_000]
end;

#%% hide
let
    β = 2.0
    teoria = exp.(β .* utilidades) ./ sum(exp.(β .* utilidades))
    sim = [mean(elecciones .== i) for i in 1:3]
    x = 1:3
    p = bar(x .- 0.18, sim; bar_width = 0.34, color = AZUL, lw = 0, label = "Simulación: elige la mayor U = u + ε",
            xticks = (x, ["Opción 1 (u = 1)", "Opción 2 (u = 0.5)", "Opción 3 (u = 0)"]),
            ylabel = "fracción que la elige", ylims = (0, 0.75), legend = :topright,
            title = "Choques de Gumbel con β = 2 producen el logit")
    bar!(p, x .+ 0.18, teoria; bar_width = 0.34, color = NARANJA, lw = 0, label = "Boltzmann: exp(βu) / Z")
end

#%% hide
md"""
Las barras coinciden. Lo que en física se llama **agitación térmica**, en
economía se llama **heterogeneidad no observada**.
"""

#%% hide
md"""
## 3. Elección con interacciones sociales: el modelo de Ising

Brock y Durlauf (2001) agregaron al logit un ingrediente: a cada agente le
importa **lo que hacen los demás**. Codifiquemos la elección del agente ``i``
como ``s_i = +1`` (adopta) o ``s_i = -1`` (no adopta). Su utilidad es

```math
U_i(s_i) = \underbrace{h\, s_i}_{\text{incentivo privado}}
+ \underbrace{J\, s_i \sum_{j \in \text{vecinos}} s_j}_{\text{conformidad}}
+ \varepsilon_i(s_i) .
```

- ``h`` es el incentivo común a adoptar (un subsidio, un precio).
- ``J > 0`` es el premio por **coincidir** con los vecinos.

La suma de las dos primeras partes, sobre todos los agentes, es menos la
energía del **modelo de Ising** (repaso F.I3):

```math
H(s) = -J\sum_{\langle i,j\rangle} s_i s_j - h\sum_i s_i .
```

Los imanes se ordenan porque cada espín "prefiere" alinearse con sus
vecinos; las poblaciones llegan a consensos porque cada agente prefiere
coincidir con los suyos. El parámetro de orden es la **magnetización**,
que en lectura social es la **adopción neta**:

```math
m = \frac{1}{N}\sum_i s_i \in [-1, 1] .
```

Con ``m = 0`` la población está dividida a la mitad; con ``m = 1`` todos
adoptan.
"""

#%% hide
md"""
## 4. Laboratorio: el algoritmo de Metropolis

Para simular el modelo, cada agente revisa su elección de vez en cuando. El
algoritmo de **Metropolis** (repaso F.I5) hace esto:

1. Elige un agente ``i`` al azar.
2. Calcula cuánto cambiaría la energía si cambiara de opción:
   ``\Delta E = 2 s_i \left(J\sum_j s_j + h\right)``.
3. Si ``\Delta E \le 0`` (le conviene), cambia. Si no, cambia con
   probabilidad ``e^{-\beta\,\Delta E}``: a veces se equivoca, por el ruido.

``N`` revisiones forman un **barrido**. Tras muchos barridos, las
configuraciones aparecen con la probabilidad de Boltzmann, ``e^{-\beta H}/Z``.
(La regla logit de la sección 1, aplicada a cada revisión, llega a la misma
distribución; se llama **dinámica de baño térmico**.)

Ponemos a los agentes en una cuadrícula de 60 × 60, cada uno con 4 vecinos,
con ``J = 1`` y ``h = 0``. El único parámetro es la temperatura ``T = 1/\beta``:
el ruido.

Temperatura ``T``: $(@bind T_red Slider([1.0, 1.5, 2.0, 2.2, 2.3, 2.5, 3.0, 4.0], default = 2.0, show_value = true))
"""

#%%
simulacion = ising_red(60, T_red; barridos = 1_000, rng = Xoshiro(2));

#%% hide
let
    s = simulacion.s
    p1 = heatmap(s; color = cgrad([NARANJA, AZUL]), clims = (-1, 1), aspect_ratio = :equal,
                 colorbar = false, axis = false, ticks = false, grid = false,
                 title = "Configuración final, T = $T_red")
    t = 1:length(simulacion.m)
    p2 = plot(t, simulacion.m; color = AZUL, lw = 1.5, ylims = (-1.05, 1.05), legend = false,
              xlabel = "barrido", ylabel = "adopción neta m",
              title = "Magnetización a lo largo de la simulación")
    hline!(p2, [0]; color = GRIS, ls = :dot)
    plot(p1, p2; layout = grid(1, 2; widths = [0.42, 0.58]), size = (900, 380))
end

#%% hide
md"""
- Con ``T`` **baja** (poco ruido), la población llega a un consenso: casi todos
  eligen lo mismo. Cuál de las dos opciones gana lo decide el azar.
- Con ``T`` **alta** (mucho ruido), las elecciones se mezclan y ``m`` fluctúa
  alrededor de 0.
- Cerca de ``T \approx 2.27`` aparecen grupos de todos los tamaños, y ``m``
  da grandes saltos. Es la **temperatura crítica**.

Para ver la transición completa, medimos ``|m|`` promediado en la segunda
mitad de la simulación para muchas temperaturas y dos tamaños de red, y lo
comparamos con la solución exacta para una red infinita (Onsager, 1944;
Yang, 1952).
"""

#%%
barrido_T = let
    Ts = 1.0:0.1:4.0
    [(L = L, T = Ts, m = [mean(abs.(ising_red(L, T; barridos = 3_000, rng = Xoshiro(round(Int, 100T)),
                                               inicio = :ordenado).m[1501:end])) for T in Ts])
     for L in (16, 48)]
end;

#%% hide
let
    Tc = 2 / log(1 + sqrt(2))
    x = range(1.0, 4.0, length = 400)
    p = plot(x, onsager.(x); color = NARANJA, lw = 2.5, label = "Solución exacta (red infinita)",
             xlabel = "temperatura T", ylabel = "|m| promedio", ylims = (0, 1.05), legend = :topright,
             title = "Transición de fase del modelo de Ising en la red cuadrada")
    for (d, c, mk) in zip(barrido_T, (RAMPA[1], RAMPA[3]), (:diamond, :circle))
        scatter!(p, d.T, d.m; color = c, ms = 4, msw = 0, marker = mk, label = "Metropolis, $(d.L) × $(d.L)")
    end
    vline!(p, [Tc]; color = GRIS, ls = :dash, label = "T_c ≈ 2.269")
end

#%% hide
md"""
Debajo de ``T_c`` la simulación sigue a la curva exacta. Por encima, la red
finita conserva una magnetización residual que se reduce al crecer la red:
con 48 × 48 está más cerca de cero que con 16 × 16. Como en la sesión 9, la
transición sólo es abrupta en el límite de muchos agentes (repaso F.B4).
"""

#%% hide
md"""
## 5. Campo medio: la solución gráfica

En la red cada agente sólo ve a 4 vecinos. Brock y Durlauf estudiaron el
caso en que a cada agente le importa el **promedio** de toda la población:
su utilidad social es ``J\, s_i\, m``. En ese caso, el campo que siente el
agente es ``h + Jm``, el mismo para todos, y la regla logit da (repaso F.I1
y F.I4)

```math
\mathbb{P}(s_i = +1) = \sigma\big(2\beta(h + Jm)\big)
\quad\Longrightarrow\quad
\mathbb{E}[s_i] = \tanh\big(\beta(h + Jm)\big) .
```

Si la población es grande, el promedio de las elecciones es ``m`` mismo:

```math
m = \tanh\big(\beta(h + Jm)\big) .
```

Es una ecuación de **punto fijo** (repaso C.I5): ``m`` aparece en los dos
lados. Se resuelve buscando dónde se cruzan la recta ``y = m`` y la curva
``y = \tanh(\beta(h + Jm))``.

Fijamos ``\beta = 1`` (así ``J`` y ``h`` se miden en unidades del ruido).

Intensidad de la interacción ``J``: $(@bind J_cm Slider([0.5, 0.8, 1.0, 1.2, 1.5, 2.0], default = 1.5, show_value = true))

Incentivo común ``h``: $(@bind h_cm Slider([-0.4, -0.2, -0.1, 0.0, 0.1, 0.2, 0.4], default = 0.0, show_value = true))
"""

#%% hide
let
    x = range(-1.2, 1.2, length = 400)
    sol = soluciones_cm(J_cm, h_cm)
    p = plot(x, x; color = GRIS, lw = 1.5, label = "y = m", xlims = (-1.2, 1.2), ylims = (-1.2, 1.2),
             aspect_ratio = :equal, xlabel = "adopción neta m", ylabel = "y", legend = :topleft,
             title = "J = $J_cm, h = $h_cm: $(length(sol) == 1 ? "1 solución" : "$(length(sol)) soluciones")",
             size = (520, 500))
    plot!(p, x, tanh.(h_cm .+ J_cm .* x); color = AZUL, lw = 2.5, label = "y = tanh(h + J m)")
    for m in sol
        estable = J_cm * (1 - tanh(h_cm + J_cm * m)^2) < 1
        scatter!(p, [m], [m]; ms = 7, color = estable ? NARANJA : BLANCO, msc = NARANJA, msw = 2,
                 label = "")
    end
    scatter!(p, [NaN], [NaN]; ms = 7, color = NARANJA, msw = 0, label = "Equilibrio estable")
    scatter!(p, [NaN], [NaN]; ms = 7, color = BLANCO, msc = NARANJA, msw = 2, label = "Equilibrio inestable")
end

#%% hide
md"""
- Con ``J < 1`` y ``h = 0`` hay una sola solución, ``m = 0``: la población se
  divide a la mitad. La conformidad es demasiado débil frente al ruido.
- Con ``J > 1`` y ``h = 0``, la pendiente de ``\tanh`` en el origen es mayor que
  1 y aparecen **dos soluciones más**, ``\pm m^\ast``: dos consensos posibles.
  ``m = 0`` se vuelve **inestable**: cualquier desviación se amplifica.
- La transición está en ``\beta J = 1``. En lenguaje de Brock y Durlauf: hay
  **equilibrios múltiples** cuando la interacción social es más fuerte que el
  ruido.
- Con ``h \ne 0`` y ``J`` grande puede haber tres soluciones, y una de ellas va
  **contra** el incentivo: la población puede estar atrapada en el consenso
  "equivocado".

**¿Funciona?** Simulamos 2,000 agentes que interactúan con el promedio de
todos (con ``h = 0``) y los comparamos con la solución gráfica.
"""

#%%
barrido_cm = let
    Js = 0.4:0.1:2.0
    (J = Js, m = [mean(abs.(campo_medio(2_000, J, 0.0; barridos = 400, rng = Xoshiro(round(Int, 10J))).m[201:end]))
                  for J in Js])
end;

#%% hide
let
    x = range(0.4, 2.0, length = 300)
    p = plot(x, [maximum(abs.(soluciones_cm(J, 0.0))) for J in x]; color = NARANJA, lw = 2.5,
             label = "Campo medio: mayor solución de m = tanh(J m)", xlabel = "interacción βJ",
             ylabel = "|m|", ylims = (0, 1.05), legend = :topleft,
             title = "2,000 agentes que miran el promedio de todos")
    scatter!(p, barrido_cm.J, barrido_cm.m; color = AZUL, ms = 5, msw = 0, label = "Simulación (baño térmico)")
    vline!(p, [1]; color = GRIS, ls = :dash, label = "βJ = 1")
end

#%% hide
md"""
Aquí la simulación sigue a la teoría en todo el rango: cuando todos
interactúan con todos, el campo medio es **exacto** para ``N`` grande. En la
red cuadrada no lo es: con ``z = 4`` vecinos el campo medio predice
``T_c = Jz = 4``, mientras que la solución exacta es 2.269. El campo medio
ignora las fluctuaciones locales y **sobreestima** el orden (repaso F.I4).
"""

#%% hide
md"""
## 6. Laboratorio: histéresis

Ahora fijamos la interacción y movemos **lentamente** el incentivo ``h``:
primero de ``-0.6`` a ``0.6`` (por ejemplo, un subsidio que crece) y luego de
regreso. En cada valor de ``h`` dejamos que la población se ajuste unos
barridos y medimos ``m``.

Interacción ``\beta J``: $(@bind J_hist Slider([0.8, 1.2, 1.5, 2.0], default = 1.5, show_value = true))
"""

#%%
ciclo = histeresis(2_000, J_hist; hs = -0.6:0.02:0.6);

#%% hide
let
    p = plot(ciclo.h_sube, ciclo.m_sube; color = AZUL, lw = 2, marker = :circle, ms = 2.5, msw = 0,
             label = "Simulación: h sube", xlabel = "incentivo común h", ylabel = "adopción neta m",
             ylims = (-1.05, 1.05), legend = :topleft, title = "Ciclo de histéresis con βJ = $J_hist")
    plot!(p, ciclo.h_baja, ciclo.m_baja; color = NARANJA, lw = 2, marker = :square, ms = 2.5, msw = 0,
          label = "Simulación: h baja")
    for (ramas, c) in ((rama_cm(J_hist, -0.6:0.005:0.6, -1.0), AZUL), (rama_cm(J_hist, 0.6:-0.005:-0.6, 1.0), NARANJA))
        plot!(p, ramas.h, ramas.m; color = c, lw = 1, ls = :dash, label = "")
    end
    plot!(p, [NaN], [NaN]; color = GRIS, ls = :dash, label = "Campo medio")
    if J_hist > 1
        hc = salto_cm(J_hist)
        vline!(p, [-hc, hc]; color = GRIS_CLARO, ls = :dot, label = "Saltos teóricos: h = ±$(round(hc, digits = 2))")
    end
    p
end

#%% hide
md"""
Con ``\beta J > 1``:

- Al subir ``h``, la población se queda en la rama de **no adopción** aunque
  el incentivo ya sea positivo. Sigue ahí hasta que esa rama desaparece, y
  entonces **salta** de golpe a casi la adopción total.
- Al bajar ``h``, se queda en la rama de **adopción** hasta un valor negativo
  de ``h``, y entonces salta hacia abajo.
- Entre los dos saltos, el estado de la población depende de **su historia**,
  no sólo del incentivo actual.

**Lecturas económicas.** Un cambio gradual en los incentivos (precios,
regulación, información) puede producir un cambio de comportamiento
**repentino**: adopción de tecnologías, cambios de normas sociales, corridas
bancarias. Y revertir el incentivo no basta para revertir el
comportamiento: hay que ir mucho más allá. Con ``\beta J < 1`` (mueve el
deslizador a 0.8) la respuesta es suave y reversible.

Los saltos de la simulación ocurren un poco **después** que los teóricos.
Cerca del punto donde desaparece una rama, la población se ajusta muy
despacio (en física se llama *frenado crítico*), y 20 barridos por valor de
``h`` no bastan para completar el salto. Con un cambio más lento del
incentivo, los saltos se acercan a los teóricos.
"""

#%% hide
md"""
## 7. [Avanzado] Campo aleatorio, avalanchas y crisis

En la sección 6 todos los agentes eran iguales salvo por el ruido. En la
realidad, cada agente tiene su propio incentivo privado ``h_i``: unos están
casi convencidos de adoptar, otros muy lejos. Este es el **modelo de Ising
con campo aleatorio** (repaso F.A3). A temperatura cero (sin ruido), el
agente ``i`` adopta si

```math
h + h_i + J m > 0 ,
```

con ``h_i`` normal con media 0 y desviación estándar ``R``: la
**heterogeneidad**. Cuando ``h`` sube un poco, adopta el agente más cercano
a su umbral; eso sube ``m``, lo que puede convencer al siguiente, y así: una
**avalancha**.

Heterogeneidad ``R``: $(@bind R_rfim Slider([0.4, 0.6, 0.8, 1.0, 1.5], default = 0.6, show_value = true))
"""

#%%
avalanchas = campo_aleatorio(10_000, 1.0, R_rfim);

#%% hide
let
    d = avalanchas
    p1 = plot(d.h, d.m; color = AZUL, lw = 2, seriestype = :steppost, ylims = (-1.05, 1.05),
              xlabel = "incentivo común h", ylabel = "adopción neta m", legend = false, xlims = (-2, 2),
              title = "J = 1, R = $R_rfim")
    s = sort(unique(d.tamaños))
    P = [mean(d.tamaños .>= x) for x in s]
    p2 = scatter(s, P; xscale = :log10, yscale = :log10, color = AZUL, ms = 3, msw = 0, legend = false,
                 xlabel = "tamaño de la avalancha (agentes)", ylabel = "P(tamaño ≥ s)",
                 title = "Avalancha más grande: $(miles(maximum(d.tamaños))) de 10,000")
    plot(p1, p2; layout = (1, 2), size = (920, 360), bottom_margin = 5Plots.mm)
end

#%% hide
md"""
- Con **poca heterogeneidad** (``R = 0.4`` o ``0.6``), hay una avalancha
  gigante: una fracción grande de la población cambia de golpe. Es una
  **crisis**.
- Con **mucha heterogeneidad** (``R = 1.5``), todas las avalanchas son
  pequeñas y la adopción crece suavemente.
- En medio, cerca de ``R_c = J\sqrt{2/\pi} \approx 0.8``, las avalanchas tienen
  **todos los tamaños**: su distribución es una ley de potencias.

Bouchaud (2013) usa este modelo para explicar cambios abruptos de opinión,
pánicos y caídas del consumo: la diversidad de opiniones protege contra los
saltos colectivos, y cuando la diversidad baja, el sistema se vuelve frágil.
La derivación de ``R_c`` está en las notas formales.
"""

#%% hide
md"""
## Resumen

| Concepto | Idea |
|---|---|
| Logit | ``\mathbb P(A) = \sigma(\beta\,\Delta u)``; ``\beta`` es la intensidad de elección |
| Logit = Boltzmann | Choques de Gumbel ↔ agitación térmica; ``E = -u``, ``T = 1/\beta`` |
| Brock-Durlauf | Logit + conformidad = modelo de Ising |
| Metropolis | Revisiones con ruido que muestrean ``e^{-\beta H}/Z`` |
| Transición de fase | Red cuadrada: consenso espontáneo debajo de ``T_c \approx 2.269\,J`` |
| Campo medio | ``m = \tanh(\beta(h + Jm))``; equilibrios múltiples si ``\beta J > 1`` |
| Histéresis | Saltos al mover ``h``; el estado depende de la historia |
| Campo aleatorio [Avanzado] | Poca heterogeneidad → crisis; mucha → cambio suave |
"""

#%% hide
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `prob_logit(Δu, β)`, la probabilidad de elegir A cuando su
ventaja en utilidad es `Δu`.
"""

#%%
function prob_logit(Δu, β)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> [prob_logit(0.0, 3.0), prob_logit(1.0, 2.0), prob_logit(-1.0, 0.5)],
          [0.5, 1 / (1 + exp(-2.0)), 1 / (1 + exp(0.5))];
          pista = md"Es la logística evaluada en ``\beta\,\Delta u``: `1 / (1 + exp(-β * Δu))`.")

#%% hide
md"""
**2.** Escribe `cambio_energia(s, suma_vecinos, J, h)`, el cambio de energía
``\Delta E`` si el agente con elección `s` cambia de opción, cuando la suma de
las elecciones de sus vecinos es `suma_vecinos`.
"""

#%%
function cambio_energia(s, suma_vecinos, J, h)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> [cambio_energia(1, 4, 1.0, 0.0), cambio_energia(-1, 2, 1.0, 0.5), cambio_energia(1, -2, 0.5, 0.0)],
          [8.0, -5.0, -2.0];
          pista = md"La fórmula es ``\Delta E = 2 s_i (J \sum_j s_j + h)``.")

#%% hide
md"""
### Ruta básica: interpretación

3. En el modelo logit, ¿qué representa ``\beta``? ¿Qué pasa con las elecciones
   cuando ``\beta \to 0`` y cuando ``\beta \to \infty``?
4. Explica con tus palabras por qué la utilidad de Brock y Durlauf lleva al
   modelo de Ising. ¿Qué papel juega cada término?
5. En la sección 0, dos ciudades con las mismas reglas terminan en consensos
   opuestos. ¿Qué implica esto para un estudio empírico que compara
   ciudades?
6. Un gobierno subsidia una tecnología y, cuando la adopción ya es alta,
   retira el subsidio. Según la sección 6, ¿qué puede pasar? ¿De qué
   depende?
7. ¿Por qué el campo medio funciona perfecto en la sección 5 y no en la red
   cuadrada de la sección 4?

### Ruta avanzada

8. Demuestra que si ``\varepsilon_A`` y ``\varepsilon_B`` son Gumbel
   independientes con escala ``1/\beta``, la diferencia
   ``\varepsilon_A - \varepsilon_B`` tiene distribución logística, y deduce la
   fórmula logit.
9. Usa el desarrollo de Taylor de ``\tanh`` para mostrar que, con ``h = 0`` y
   ``\beta J`` un poco mayor que 1, ``m^\ast \approx \sqrt{3(\beta J - 1)}``.
   Compara con la gráfica de la sección 5.
10. Repite el barrido de temperatura de la sección 4 con una red de
    Erdős-Rényi con grado promedio 4 en lugar de la cuadrícula. ¿Dónde está
    ahora la transición? ¿Se acerca más al campo medio?
11. En la sección 7, mide el tamaño de la avalancha más grande para varios
    valores de ``R`` y para ``N = 1{,}000`` y ``10{,}000``. ¿Cómo ubicarías
    ``R_c`` con estos datos?

## Para la próxima sesión

- Lectura: Sen y Chakrabarti (2014), capítulos 2 y 3; Requião da Cunha
  (2021), capítulo 6. Para la ruta avanzada, Brock y Durlauf (2001) y
  Bouchaud, Marsili y Nadal (2023).
- En la sesión 15 empieza el módulo de agentes: modelos de opinión
  (votante, Kirman), herding en mercados y segregación de Schelling. La
  imitación de hoy será el mecanismo central.
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
    using PlutoUI, Statistics, Random, Distributions, Plots, Markdown
    gr()
    default(size = (720, 320), legend = false, grid = :y, gridalpha = 0.15,
            framestyle = :axes, titlefontsize = 11, guidefontsize = 9,
            tickfontsize = 8, legendfontsize = 8, foreground_color_legend = nothing,
            background_color_legend = nothing, margin = 3Plots.mm)
end

#%%
begin
    const AZUL       = "#2a78d6"   # serie principal, opción +1
    const NARANJA    = "#eb6834"   # segunda serie o teoría, opción -1
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const BLANCO     = "#ffffff"
    const RAMPA      = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # azul secuencial
end;

#%%
"Función logística."
σ(x) = 1 / (1 + exp(-x))

#%%
"Número entero con comas como separador de miles."
miles(x) = replace(string(round(Int, x)), r"(?<=\d)(?=(\d{3})+$)" => ",")

#%%
"Datos ordenados de mayor a menor y su CCDF empírica, P(X ≥ x)."
function ccdf_empirica(x)
    s = sort(x; rev = true)
    s, (1:length(s)) ./ length(s)
end

#%%
"""
Modelo de Ising en una cuadrícula L × L con fronteras periódicas, J = 1 y
h = 0, simulado con Metropolis. `inicio` puede ser `:azar` u `:ordenado`
(todos en +1). Devuelve la configuración final y la magnetización después
de cada barrido.
"""
function ising_red(L, T; barridos = 1_000, rng = Xoshiro(1), inicio = :azar)
    s = inicio == :azar ? rand(rng, [-1, 1], L, L) : ones(Int, L, L)
    β = 1 / T
    acepta = [exp(-β * ΔE) for ΔE in 1:8]    # probabilidades para ΔE > 0
    m = Float64[]
    for _ in 1:barridos
        for _ in 1:L^2
            i, j = rand(rng, 1:L), rand(rng, 1:L)
            vecinos = s[mod1(i + 1, L), j] + s[mod1(i - 1, L), j] + s[i, mod1(j + 1, L)] + s[i, mod1(j - 1, L)]
            ΔE = 2 * s[i, j] * vecinos
            if ΔE <= 0 || rand(rng) < acepta[ΔE]
                s[i, j] = -s[i, j]
            end
        end
        push!(m, mean(s))
    end
    (s = s, m = m)
end

#%%
"Magnetización espontánea exacta de la red cuadrada (Onsager y Yang), con J = 1."
onsager(T) = T < 2 / log(1 + sqrt(2)) ? (1 - sinh(2 / T)^(-4))^(1 / 8) : 0.0

#%%
"""
Modelo de campo medio con N agentes y β = 1: cada revisión elige un agente
al azar y le asigna +1 con probabilidad logística σ(2(h + J m)), donde m es
la adopción neta de los demás (baño térmico). Devuelve la última
configuración y m después de cada barrido.
"""
function campo_medio(N, J, h; barridos = 400, rng = Xoshiro(1), s = rand(rng, [-1, 1], N))
    s = copy(s)
    total = sum(s)
    m = Float64[]
    for _ in 1:barridos
        for _ in 1:N
            i = rand(rng, 1:N)
            m_otros = (total - s[i]) / (N - 1)
            nuevo = rand(rng) < σ(2 * (h + J * m_otros)) ? 1 : -1
            total += nuevo - s[i]
            s[i] = nuevo
        end
        push!(m, total / N)
    end
    (s = s, m = m)
end

#%%
"Soluciones de m = tanh(h + J m) en [-1, 1], encontradas por cambio de signo y bisección."
function soluciones_cm(J, h)
    f(m) = tanh(h + J * m) - m
    xs = range(-1.0, 1.0, length = 2001)
    raices = Float64[]
    for k in 1:length(xs)-1
        a, b = xs[k], xs[k+1]
        f(a) == 0 && push!(raices, a)
        if f(a) * f(b) < 0
            for _ in 1:60
                c = (a + b) / 2
                f(a) * f(c) <= 0 ? (b = c) : (a = c)
            end
            push!(raices, (a + b) / 2)
        end
    end
    raices
end

#%%
"Rama de campo medio que se sigue al mover h lentamente: se itera m = tanh(h + J m) desde la solución anterior."
function rama_cm(J, hs, m0)
    m = m0
    ms = map(hs) do h
        for _ in 1:5_000
            m = tanh(h + J * m)
        end
        m
    end
    (h = collect(hs), m = ms)
end

#%%
"Valor de h donde desaparece una rama en el campo medio (J > 1)."
function salto_cm(J)
    m = sqrt(1 - 1 / J)       # donde la pendiente de tanh(h + J m) es 1
    J * m - atanh(m)
end

#%%
"Ciclo de histéresis simulado: h sube y luego baja, con 20 barridos en cada valor."
function histeresis(N, J; hs = -0.6:0.02:0.6, barridos = 20, rng = Xoshiro(4))
    s = fill(-1, N)
    m_sube, m_baja = Float64[], Float64[]
    for h in hs
        r = campo_medio(N, J, h; barridos, rng, s)
        s = r.s
        push!(m_sube, mean(r.m[end-4:end]))
    end
    for h in reverse(hs)
        r = campo_medio(N, J, h; barridos, rng, s)
        s = r.s
        push!(m_baja, mean(r.m[end-4:end]))
    end
    (h_sube = collect(hs), m_sube = m_sube, h_baja = collect(reverse(hs)), m_baja = m_baja)
end

#%%
"""
Ising con campo aleatorio a temperatura cero en campo medio. Los incentivos
privados h_i son normales con desviación R. Se sube h desde -∞: los agentes
adoptan en orden de h_i decreciente, y cada adopción sube el campo de todos
en 2J/N. Devuelve m después de cada avalancha, el h en que ocurrió y el
tamaño de cada avalancha.
"""
function campo_aleatorio(N, J, R; rng = Xoshiro(5))
    hi = sort(R .* randn(rng, N); rev = true)
    arriba = 0                       # agentes que ya adoptaron
    hs, ms, tamaños = Float64[], Float64[], Int[]
    while arriba < N
        m = 2arriba / N - 1
        h = -hi[arriba+1] - J * m    # h justo para que adopte el siguiente
        inicio = arriba
        arriba += 1
        while arriba < N && h + hi[arriba+1] + J * (2arriba / N - 1) > 0
            arriba += 1
        end
        push!(hs, h)
        push!(ms, 2arriba / N - 1)
        push!(tamaños, arriba - inicio)
    end
    (h = [-5.0; hs], m = [-1.0; ms], tamaños = tamaños)
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
