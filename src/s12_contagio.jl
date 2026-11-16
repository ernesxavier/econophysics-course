#%% hide
md"""
# Sesión 12: Contagio y riesgo sistémico

**Objetivos de la sesión**
1. Usar la red como el **medio** por donde se propaga algo: una adopción,
   una quiebra, una pérdida.
2. Simular las **cascadas de umbral** de Watts y encontrar la **ventana de
   cascadas**: el rango de conectividad donde un choque pequeño puede
   arrastrar a toda la red.
3. Aplicar la misma idea a una **red interbancaria**: cuándo la quiebra de un
   banco provoca otras.
4. Medir con **DebtRank** cómo se amplifican las pérdidas antes de que haya
   quiebras, e identificar a los bancos "demasiado centrales para caer".
5. Entender por qué los sistemas financieros son **robustos pero frágiles**.

**Cómo usar este cuaderno.** Pluto es reactivo: al mover un deslizador, todo
lo que depende de él se recalcula. Las secciones marcadas **[Avanzado]** son
opcionales. Algunas celdas simulan miles de cascadas y tardan unos segundos.

**Cómo abrirlo.** Desde la raíz del repositorio:
`julia --project=Code -e 'using Pluto; Pluto.run()'` y abrir
`Code/s12_contagio.jl`.

**Repaso recomendado** (en `Class_notes/repaso/`): L.B4 (matriz de
adyacencia), C.B5 (serie geométrica), L.I4 (potencias de matrices), C.I5
(puntos fijos). Ruta avanzada: L.A4 (series de matrices y puntos fijos),
F.A5 (percolación).

**Conceptos transversales.** Hoy se presentan los **puntos fijos y su
estabilidad** (un choque que se repite atenuándose); se completan en la sesión
14.
"""

#%% hide
TableOfContents(title = "Contenido")

#%% hide
md"""
## 0. El fenómeno: de una quiebra a una crisis

El 15 de septiembre de 2008 quebró Lehman Brothers, un banco de inversión
que no era el más grande de Estados Unidos. En las semanas siguientes, un
fondo de mercado de dinero que tenía deuda de Lehman "rompió el dólar", la
aseguradora AIG tuvo que ser rescatada porque había garantizado deudas de
medio sistema financiero, y los bancos dejaron de prestarse entre sí. Una
falla local se volvió una crisis global.

El **riesgo sistémico** es justamente eso: que el sistema falle por la forma
en que sus partes están conectadas, aunque cada parte, por separado, parezca
sana. Abajo hay una versión de juguete: 200 nodos en una red al azar. Cada
uno cae si al menos el 18% de sus vecinos ha caído. Empezamos tumbando **un
solo** nodo; el color dice en qué ronda cayó cada uno.
"""

#%% hide
let
    g, semilla = ejemplo_cascada
    ronda = rondas_cascada(g, 0.18, semilla)
    X = disposicion(g)
    maxr = maximum(filter(>=(0), ronda))
    paleta = cgrad([NARANJA, "#f6c6ad"])
    colores = [r < 0 ? GRIS_CLARO : paleta[r / max(maxr, 1)] for r in ronda]
    p = dibujar_red(g, X; tam = 4, colores, size = (620, 480),
                    title = "Una cascada: $(count(>=(0), ronda)) de $(nv(g)) nodos caen a partir de uno")
    scatter!(p, [X[semilla, 1]], [X[semilla, 2]]; ms = 9, color = NARANJA, msw = 2, msc = GRIS,
             label = "Nodo inicial")
    scatter!(p, [NaN], [NaN]; color = NARANJA, ms = 6, msw = 0, label = "Cae en las primeras rondas")
    scatter!(p, [NaN], [NaN]; color = "#f6c6ad", ms = 6, msw = 0, label = "Cae en las últimas (hasta la ronda $maxr)")
    scatter!(p, [NaN], [NaN]; color = GRIS_CLARO, ms = 6, msw = 0, label = "No cae")
    plot!(p; legend = :bottomleft)
end

#%% hide
md"""
Con otra red, o con otro nodo inicial, la misma regla puede no producir casi
nada. Hoy queremos saber **cuándo** un choque pequeño se vuelve sistémico.
"""

#%% hide
md"""
## 1. Repaso: un choque que se repite atenuándose

Supón que cada banco, al perder ``1`` peso, hace perder ``\ell`` pesos en
total a sus acreedores, con ``0 < \ell < 1``. Un choque inicial de 1 produce
``\ell`` en la segunda ronda, ``\ell^2`` en la tercera, y así. La pérdida total
es una **serie geométrica** (repaso C.B5):

```math
1 + \ell + \ell^2 + \cdots = \frac{1}{1 - \ell} .
```

Con ``\ell = 0.5`` el choque se duplica; con ``\ell = 0.9`` se multiplica por
10; con ``\ell \ge 1`` la serie no converge: el choque **no se apaga**. En una
red, ``\ell`` se vuelve una matriz ``\Lambda`` y la serie es
``s + \Lambda s + \Lambda^2 s + \cdots`` (repaso L.I4 y L.A4). Lo veremos en la
sección 4.
"""

#%%
amplificacion(ℓ; rondas = 200) = sum(ℓ^t for t in 0:rondas)

#%%
amplificacion(0.5), amplificacion(0.9), amplificacion(1.1)

#%% hide
md"""
## 2. Cascadas de umbral

Watts (2002) propuso el modelo más simple de una cascada en red. Sirve para
adoptar una tecnología, sumarse a una corrida bancaria o quebrar.

**Reglas.**
1. Cada nodo tiene un **umbral** ``\varphi``: cae (o adopta) cuando la
   **fracción** de sus vecinos que ya cayó es al menos ``\varphi``.
2. Al principio cae un solo nodo, elegido al azar.
3. Se revisan los nodos una y otra vez hasta que nadie más cambia.

Un nodo con grado ``k`` cae con **un solo** vecino caído si ``1/k \ge \varphi``,
es decir, si ``k \le 1/\varphi``. A esos nodos les llamamos **vulnerables**.
Con ``\varphi = 0.18``, son vulnerables los nodos con 5 vecinos o menos.

Para que una cascada llegue lejos tienen que pasar dos cosas:

- La red debe estar **conectada**: si el grado promedio ``z`` es menor que 1,
  no hay por dónde propagarse (sesión 10).
- Debe haber suficientes nodos **vulnerables** conectados entre sí: si ``z``
  es grande, casi nadie es vulnerable y la cascada se detiene.

### Laboratorio: la ventana de cascadas

Simulamos 100 cascadas en redes de Erdős-Rényi con 5,000 nodos para cada
grado promedio ``z``, cada una a partir de un nodo al azar. Una cascada es
**global** si alcanza a más del 5% de la red.

Umbral ``\varphi``: $(@bind φ_watts Slider([0.12, 0.15, 0.18, 0.22, 0.25], default = 0.18, show_value = true))
"""

#%%
ventana_er = barrido_cascadas(φ_watts, zs_ventana;
                              generar = (n, z, r) -> erdos_renyi(n, round(Int, n * z / 2); seed = r));

#%% hide
let
    lim = frontera_watts(φ_watts)
    p = plot(; xlabel = "grado promedio z", ylabel = "fracción", ylims = (0, 1.05),
             xlims = (0, maximum(zs_ventana)), legend = :right, size = (760, 340),
             title = "Ventana de cascadas en Erdős-Rényi, φ = $φ_watts")
    length(lim) == 2 && vspan!(p, lim; color = NARANJA, alpha = 0.1,
                               label = "Ventana teórica: $(round(lim[1], digits = 1)) < z < $(round(lim[2], digits = 1))")
    plot!(p, zs_ventana, ventana_er.probabilidad; color = AZUL, lw = 2, marker = :circle, ms = 3, msw = 0,
          label = "Probabilidad de una cascada global")
    plot!(p, zs_ventana, ventana_er.tamaño; color = NARANJA, lw = 2, marker = :square, ms = 3, msw = 0,
          label = "Tamaño promedio de las cascadas globales")
end

#%% hide
md"""
Hay una **ventana**: las cascadas globales sólo ocurren para valores
intermedios de ``z``. Los dos bordes tienen causas distintas.

- **Borde inferior** (``z \approx 1``): es la transición de percolación de la
  sesión 10. Con pocos enlaces, la red está partida en islas.
- **Borde superior**: con muchos enlaces, cada vecino pesa poco (``1/k``) y
  un solo vecino caído no basta. La red se **estabiliza por diversificación**.

Mira las dos curvas cerca del borde superior. La probabilidad de una cascada
global cae a casi cero, pero el **tamaño** de las que sí ocurren sube a casi
toda la red. Un sistema en esa zona parece muy estable: los choques
pequeños casi nunca pasan a mayores, y por eso no hay experiencia que
advierta del peligro. Pero cuando una cascada logra arrancar, se lleva todo.
Es el sentido de **robusto pero frágil**.

Mueve el umbral. Un ``\varphi`` más alto (nodos más resistentes) estrecha la
ventana, y por encima de ``\varphi = 0.25`` la cierra por completo; uno más
bajo la abre hacia ``z`` grandes. La sección 6 explica de
dónde sale la ventana teórica.
"""

#%% hide
md"""
### Laboratorio: Erdős-Rényi contra Barabási-Albert

¿Importa la forma de la red, y no sólo cuántos enlaces tiene? Repitamos el
experimento en redes de Barabási-Albert con el mismo grado promedio, con
``\varphi`` fijo en 0.18.
"""

#%%
comparacion_cascadas = let
    zs = 2:2:12
    er = barrido_cascadas(0.18, zs; generar = (n, z, r) -> erdos_renyi(n, round(Int, n * z / 2); seed = r))
    ba = barrido_cascadas(0.18, zs; generar = (n, z, r) -> barabasi_albert(n, z ÷ 2; seed = r))
    (z = zs, er = er, ba = ba)
end;

#%% hide
let
    d = comparacion_cascadas
    p = plot(; xlabel = "grado promedio z", ylabel = "probabilidad de una cascada global",
             ylims = (0, 1.05), legend = :topright, xticks = d.z,
             title = "Mismo grado promedio, distinta forma (φ = 0.18)")
    plot!(p, d.z, d.er.probabilidad; color = NARANJA, lw = 2, marker = :circle, ms = 4, msw = 0,
          label = "Erdős-Rényi")
    plot!(p, d.z, d.ba.probabilidad; color = AZUL, lw = 2, marker = :square, ms = 4, msw = 0,
          label = "Barabási-Albert")
end

#%% hide
md"""
Con conectividad baja, en Barabási-Albert las cascadas son **mucho más
raras**: los concentradores tienen tantos vecinos que no son vulnerables, y
como casi todos los caminos pasan por ellos, actúan como cortafuegos. (Con
``z = 2`` la red de Barabási-Albert es un árbol en el que casi todos los
caminos cruzan algún concentrador, y la cascada no pasa.)

Cerca del borde superior las dos redes se parecen, y con ``z = 8`` la de
Barabási-Albert todavía tiene alguna cascada global y la de Erdős-Rényi ya
no. En Barabási-Albert todos los nodos tienen al menos ``z/2`` enlaces, y ése
es el grado más frecuente: con ``z = 8``, un tercio de los nodos tiene grado 4
y casi otro quinto grado 5, ambos todavía vulnerables. La forma de la red cambia **dónde** está la ventana,
no sólo su tamaño.
"""

#%% hide
md"""
## 3. Laboratorio: contagio en una red interbancaria

Gai y Kapadia (2010) tradujeron el modelo de Watts a bancos. Cada banco tiene
un **balance**:

| Activos | Pasivos |
|---|---|
| Préstamos a otros bancos (interbancarios): ``\theta`` | Depósitos y otras deudas |
| Otros activos (créditos, bonos): ``1 - \theta`` | **Capital**: ``K`` |

El capital es lo que el banco puede perder sin quebrar. Si un banco presta
sus activos interbancarios en partes iguales a ``k`` deudores y ``f`` de ellos
quiebran sin pagar nada, pierde ``\theta f/k``. Quiebra cuando esa pérdida
supera su capital:

```math
\frac{f}{k} > \frac{K}{\theta} .
```

¡Es el modelo de Watts con umbral ``\varphi = K/\theta``! Con ``\theta = 20\%``
de los activos en préstamos interbancarios y ``K = 4\%`` de capital, el umbral
es ``\varphi = 0.2``. La red ahora es **dirigida**: si ``j`` quiebra, pierden los
bancos que le prestaron, no los que le deben.

Simulamos 200 sistemas de 1,000 bancos; en cada uno quiebra un banco al azar.

Capital ``K`` (fracción de los activos): $(@bind K_bancos Slider([0.02, 0.03, 0.04, 0.05, 0.06], default = 0.04, show_value = true))
"""

#%%
ventana_bancos = barrido_bancos(K_bancos, 1:12);

#%% hide
let
    d = ventana_bancos
    p = plot(; xlabel = "número promedio de deudores por banco", ylabel = "fracción",
             ylims = (0, 1.05), legend = :right, xticks = 1:12, size = (760, 340),
             title = "Contagio interbancario con capital K = $(round(Int, 100K_bancos))% y θ = 20%")
    plot!(p, 1:12, d.probabilidad; color = AZUL, lw = 2, marker = :circle, ms = 4, msw = 0,
          label = "Probabilidad de contagio (más del 5% de los bancos)")
    plot!(p, 1:12, d.tamaño; color = NARANJA, lw = 2, marker = :square, ms = 4, msw = 0,
          label = "Fracción de bancos que quiebran, si hay contagio")
end

#%% hide
md"""
**Diversificación contra contagio.** Prestar a más bancos tiene dos
efectos opuestos:

- Cada préstamo es más pequeño, así que la quiebra de un deudor pesa menos
  (``\theta/k`` baja). Eso **protege**.
- Hay más canales por donde llega el contagio. Eso **expone**.

Con pocos deudores domina el segundo efecto y el contagio aumenta con la
conectividad; con muchos domina el primero y desaparece. Pero cerca del
borde, cuando hay contagio, quiebra prácticamente **todo** el sistema. Antes
de 2008 se pensaba que más conexiones siempre significaban más
diversificación y más estabilidad; este modelo muestra por qué eso no basta.

Sube el capital: la ventana se estrecha. Un capital del 6% con 20% de
préstamos interbancarios da ``\varphi = 0.3``, y ningún banco con más de 3
deudores cae por la quiebra de uno solo.
"""

#%% hide
md"""
## 4. DebtRank: pérdidas antes de la quiebra

En el modelo anterior sólo importan las **quiebras**: un banco que pierde la
mitad de su capital no le hace nada a nadie. Pero si un deudor está en
problemas, el préstamo que le hice vale menos aunque todavía no haya
quebrado. Battiston y sus coautores (2012) propusieron medir cómo se
propaga ese **deterioro**.

Sea ``h_i`` la fracción del capital que el banco ``i`` ha perdido: 0 si está
intacto, 1 si quebró. Si ``A_{ij}`` es lo que ``i`` le prestó a ``j`` y ``E_i``
el capital de ``i``, la regla (en la versión de Bardoscia et al., 2015) es

```math
h_i(t+1) = \min\left\{1,\ h_i(t) + \sum_j \Lambda_{ij}\,\big[h_j(t) - h_j(t-1)\big]\right\},
\qquad \Lambda_{ij} = \frac{A_{ij}}{E_i} .
```

En palabras: si mi deudor ``j`` pierde una fracción adicional de su capital,
el valor de mi préstamo baja en esa misma fracción, y la pérdida, medida
contra mi capital, es ``\Lambda_{ij}`` veces mayor. La matriz ``\Lambda`` es la
matriz de **impacto**.

### Un sistema bancario sintético

Construimos 50 bancos con tamaños de **Pareto** (unos pocos bancos muy
grandes, como en la realidad; sesión 3). Cada banco presta a 3 bancos,
elegidos con probabilidad proporcional a su tamaño, y le presta más a los
más grandes. El parámetro clave es el **apalancamiento interbancario**

```math
\ell = \frac{\text{préstamos interbancarios}}{\text{capital}} ,
```

igual para todos. Cada fila de ``\Lambda`` suma ``\ell``.

Apalancamiento interbancario ``\ell``: $(@bind ℓ_banco Slider([0.25, 0.5, 0.75, 0.9, 1.0, 1.5, 2.0], default = 0.75, show_value = true))
"""

#%%
bancos = sistema_bancario(50; ℓ = ℓ_banco);

#%% hide
let
    ℓs = [0.1, 0.25, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.93, 0.96]
    amp = map(ℓs) do ℓ
        b = sistema_bancario(50; ℓ)
        h0 = fill(0.01, 50)
        perdida(b, debtrank(b, h0)) / perdida(b, h0)
    end
    x = range(0, 0.97, length = 200)
    p1 = plot(x, 1 ./ (1 .- x); color = NARANJA, lw = 2, yscale = :log10, ylims = (1, 100),
              yticks = ([1, 2, 5, 10, 20, 50, 100], ["1", "2", "5", "10", "20", "50", "100"]),
              xlabel = "apalancamiento interbancario ℓ", ylabel = "pérdida final / pérdida inicial",
              label = "Serie geométrica: 1/(1 − ℓ)", legend = :topleft, xlims = (0, 1.05),
              title = "Choque del 1% del capital a todos los bancos")
    scatter!(p1, ℓs, amp; color = AZUL, ms = 4, msw = 0, label = "DebtRank (simulación)")
    vline!(p1, [1]; color = GRIS, ls = :dash, label = "ℓ = 1: el choque no se apaga")
    hist = historia_debtrank(bancos, fill(0.01, 50); rondas = 100)
    p2 = plot(0:100, 100 .* hist; color = AZUL, lw = 2,
              yscale = :log10, ylims = (0.8, 130), xlabel = "ronda de contagio",
              yticks = ([1, 2, 5, 10, 20, 50, 100], ["1", "2", "5", "10", "20", "50", "100"]),
              ylabel = "capital perdido por el sistema (%)", label = "DebtRank, ℓ = $ℓ_banco",
              legend = :bottomright, title = "La pérdida, ronda por ronda")
    ℓ_banco < 1 && hline!(p2, [1 / (1 - ℓ_banco)]; color = NARANJA, ls = :dash, lw = 1.5,
                          label = "Límite: 1/(1 − ℓ) = $(round(1 / (1 - ℓ_banco), digits = 1))%")
    plot(p1, p2; layout = (1, 2), size = (920, 360))
end

#%% hide
md"""
Damos a todos los bancos un choque pequeño: el 1% de su capital. La
gráfica de la derecha sigue la pérdida del sistema ronda por ronda. Como cada fila de ``\Lambda`` suma ``\ell``, cada ronda
de contagio repite la pérdida anterior multiplicada por ``\ell``, y la
pérdida total es exactamente la serie geométrica de la sección 1.

- Con ``\ell < 1`` el choque se **amplifica pero se apaga**: con ``\ell = 0.75``,
  la pérdida final es 4 veces la inicial.
- Con ``\ell \ge 1`` el choque **no se apaga**: por pequeño que sea, crece
  hasta que todos los bancos pierden todo su capital (el 100%). Mueve el
  deslizador a 1.0 o más y mira la gráfica de la derecha.

El umbral ``\ell = 1`` no depende de la forma de la red, sólo de que el
**radio espectral** de ``\Lambda`` (su mayor eigenvalor en valor absoluto)
cruce 1 (repaso L.A4). No es un valor exagerado: un banco con 10% de sus
activos en préstamos a otros bancos y 5% de capital ya tiene ``\ell = 2``.
Bardoscia et al. (2017) estudian cómo se acercó el sistema bancario europeo
a esta zona inestable antes de la crisis.

### Laboratorio: demasiado central para caer

El **DebtRank** de un banco ``j`` es la fracción del capital del **resto** del
sistema que se pierde si ``j`` quiebra (``h_j = 1`` al principio y los demás en
0). ¿Son los bancos más grandes los más peligrosos?
"""

#%%
impacto = [debtrank_de(bancos, j) for j in 1:50];

#%% hide
let
    b = bancos
    acreedores = vec(count(>(0), b.A; dims = 1))
    p = scatter(b.activos, impacto; xscale = :log10, xticks = ([0.3, 1, 3, 10], ["0.3", "1", "3", "10"]), ms = 3 .+ 0.5 .* acreedores, color = AZUL,
                alpha = 0.7, msw = 0, ylims = (0, 1.05), legend = false,
                xlabel = "tamaño del banco (activos, escala log)",
                ylabel = "DebtRank: capital perdido por los demás",
                title = "Impacto de la quiebra de cada banco (ℓ = $ℓ_banco; tamaño del punto: acreedores)")
    for j in sortperm(impacto; rev = true)[1:3]
        annotate!(p, b.activos[j], impacto[j] + 0.05,
                  text("$(acreedores[j]) acreedores", 8, GRIS))
    end
    p
end

#%% hide
md"""
El tamaño importa, pero no lo explica todo. Algunos bancos medianos tienen
un impacto comparable al de los más grandes: les prestan bancos que a su
vez son acreedores de otros, y su caída se propaga en cadena. Battiston y
sus coautores (2012) encontraron lo mismo con los datos de los préstamos de
emergencia de la Reserva Federal en 2008-2010: un grupo de unas 20
instituciones estaba tan interconectado que cada una era **demasiado
central para caer**, no sólo demasiado grande.

Con ``\ell \ge 1`` todos los bancos tienen DebtRank cercano a 1: cualquier
quiebra hunde al sistema entero.
"""

#%% hide
md"""
## 5. [Avanzado] Eisenberg-Noe: compensación como punto fijo

DebtRank supone que una pérdida se transmite completa. Eisenberg y Noe (2001)
plantearon la pregunta contable precisa: si algunos bancos no pueden pagar
todo lo que deben, **¿cuánto paga cada uno?**

Sea ``L_{ij}`` lo que ``i`` le debe a ``j``, ``\bar p_i = \sum_j L_{ij}`` la deuda
total de ``i`` y ``\Pi_{ij} = L_{ij}/\bar p_i`` la fracción que le toca a cada
acreedor. Cada banco tiene activos externos ``e_i`` y una deuda externa ``v_i``
que tiene prioridad. El **vector de pagos** ``p`` cumple

```math
p_i = \max\left\{0,\ \min\left\{\bar p_i,\ e_i - v_i + \sum_j \Pi_{ji}\, p_j\right\}\right\}.
```

Cada banco paga todo lo que debe o, si no le alcanza, todo lo que tiene. Es
una ecuación de **punto fijo** ``p = \Phi(p)`` (repaso C.I5 y L.A4).
``\Phi`` es creciente y está acotada, así que iterar desde ``p = \bar p``
produce una sucesión que baja monótonamente hasta el mayor punto fijo.

**Anillo contra red completa.** Acemoglu, Ozdaglar y Tahbaz-Salehi (2015)
compararon dos estructuras con 10 bancos que deben 1 cada uno a otros bancos
y tienen un colchón de 0.08 (activos externos menos deuda externa):

- **Anillo**: cada banco le debe todo al siguiente.
- **Red completa**: cada banco le debe ``1/9`` a cada uno de los demás.

Un choque reduce los activos externos del banco 1 en ``\varepsilon``.
"""

#%%
acemoglu = let
    n = 10
    anillo = zeros(n, n)
    for i in 1:n
        anillo[i, mod1(i + 1, n)] = 1.0
    end
    completa = fill(1 / (n - 1), n, n)
    completa[diagind(completa)] .= 0
    εs = 0:0.02:1.4
    quiebras(L, ε) = (e = fill(1.08, n); e[1] -= ε; count(compensacion(L, e, fill(1.0, n)).quiebra))
    (ε = εs, anillo = [quiebras(anillo, ε) for ε in εs], completa = [quiebras(completa, ε) for ε in εs])
end;

#%% hide
let
    d = acemoglu
    p = plot(d.ε, d.anillo; color = NARANJA, lw = 2, seriestype = :steppost, label = "Anillo",
             xlabel = "tamaño del choque ε al banco 1", ylabel = "bancos en quiebra",
             ylims = (0, 10.5), yticks = 0:2:10, legend = :topleft,
             title = "Eisenberg-Noe: 10 bancos con colchón de 0.08")
    plot!(p, d.ε, d.completa; color = AZUL, lw = 2, seriestype = :steppost, label = "Red completa")
    vline!(p, [0.8]; color = GRIS, ls = :dash, label = "ε = colchón total del sistema")
end

#%% hide
md"""
- Con **choques pequeños**, la red completa es la más estable: la pérdida del
  banco 1 se reparte entre 9 acreedores, cada uno pierde poco y nadie más
  quiebra. En el anillo, toda la pérdida cae sobre un solo acreedor, que la
  pasa al siguiente, y caen varios en cadena.
- Con **choques grandes**, la red completa pasa de golpe de 1 quiebra a 10.
  En cuanto el choque supera el colchón de **todo** el sistema, repartirlo
  sólo garantiza que todos lo sufran.

Es una **transición de fase** en el tamaño del choque: la misma conexión que
absorbe los choques pequeños propaga los grandes. Acemoglu y sus coautores
demuestran que, para choques grandes, las redes más estables son las
**débilmente conectadas**: grupos densos unidos por pocos enlaces.
"""

#%% hide
md"""
## 6. [Avanzado] De dónde sale la ventana de cascadas

En la sección 2, una cascada global requiere un **grupo conectado de nodos
vulnerables** que crezca sin límite. Es un problema de percolación (sesión 11,
repaso F.A5) restringido a los nodos vulnerables.

Al seguir un enlace se llega a un nodo de grado ``k`` con probabilidad
``k p_k/z``. Ese nodo es vulnerable si ``k \le 1/\varphi``, y en ese caso abre
``k - 1`` caminos nuevos. El número promedio de nodos vulnerables nuevos que
se alcanzan es

```math
\sum_{k \le 1/\varphi} (k - 1)\,\frac{k\, p_k}{z} .
```

La cascada crece si este número es mayor que 1 (Watts, 2002):

```math
\sum_{k=1}^{\lfloor 1/\varphi\rfloor} k(k-1)\, p_k > z .
```

Para Erdős-Rényi, ``p_k = e^{-z} z^k/k!``. Los dos valores de ``z`` donde se
cumple la igualdad son los bordes de la ventana. Con ``\varphi = 0.18`` son
$(join(round.(frontera_watts(0.18), digits = 2), " y ")), que coinciden con la
simulación.
"""

#%% hide
let
    φs = 0.08:0.005:0.26
    bordes = [frontera_watts(φ) for φ in φs]
    con = findall(b -> length(b) == 2, bordes)
    p = plot(first.(bordes[con]), φs[con]; color = AZUL, lw = 2, label = "Borde inferior",
             xlabel = "grado promedio z", ylabel = "umbral φ", xlims = (0, 16),
             legend = :topright, title = "Ventana de cascadas de Erdős-Rényi (teoría)")
    plot!(p, last.(bordes[con]), φs[con]; color = NARANJA, lw = 2, label = "Borde superior")
    annotate!(p, 4, 0.12, text("cascadas globales", 9, GRIS))
    annotate!(p, 12, 0.22, text("sin cascadas", 9, GRIS))
end

#%% hide
md"""
## Resumen

| Concepto | Idea |
|---|---|
| Cascada de umbral (Watts) | Un nodo cae si cae una fracción ``\varphi`` de sus vecinos; vulnerables: ``k \le 1/\varphi`` |
| Ventana de cascadas | Sólo hay cascadas globales con conectividad intermedia |
| Robusto pero frágil | Cerca del borde superior, las cascadas son raras pero abarcan casi todo |
| Gai-Kapadia | Bancos: umbral ``\varphi`` = capital / préstamos interbancarios |
| DebtRank | Las pérdidas se propagan antes de las quiebras; se amplifican por ``1/(1-\ell)`` |
| Radio espectral | Si el de ``\Lambda`` es mayor o igual que 1, un choque pequeño no se apaga |
| Demasiado central para caer | El impacto de una quiebra depende de la red, no sólo del tamaño |
| Eisenberg-Noe [Avanzado] | Pagos como punto fijo; transición de fase entre choques pequeños y grandes |
"""

#%% hide
md"""
## Ejercicios

### Ruta básica: completa el código

**1.** Escribe `cae(caidos, k, φ)`, que diga si un nodo con ``k`` vecinos, de
los cuales `caidos` ya cayeron, cae con umbral ``\varphi``. Un nodo sin
vecinos no cae.
"""

#%%
function cae(caidos, k, φ)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> [cae(1, 5, 0.18), cae(1, 6, 0.18), cae(2, 6, 0.18), cae(0, 0, 0.18)],
          [true, false, true, false];
          pista = md"Cae si `k > 0` y `caidos / k ≥ φ`, es decir, `k > 0 && caidos ≥ φ * k`.")

#%% hide
md"""
**2.** Escribe `matriz_impacto(A, E)`, que devuelva ``\Lambda_{ij} = A_{ij}/E_i``:
cada fila de ``A`` dividida entre el capital del banco de esa fila.
"""

#%%
function matriz_impacto(A, E)
    missing  # reemplaza esta línea
end

#%% hide
verificar(() -> matriz_impacto([0.0 2.0; 1.0 0.0], [4.0, 2.0]), [0.0 0.5; 0.5 0.0];
          pista = md"Con difusión: `A ./ E` divide cada fila ``i`` entre `E[i]` si `E` es un vector columna.")

#%% hide
md"""
### Ruta básica: interpretación

3. Explica con tus palabras los dos bordes de la ventana de cascadas. ¿Cuál
   de los dos tiene que ver con la diversificación?
4. Un regulador observa diez años sin ninguna crisis bancaria y concluye
   que el sistema es estable. Según la sección 3, ¿por qué podría estar
   equivocado?
5. Da un argumento a favor y uno en contra de rescatar a un banco que es
   "demasiado central para caer", usando lo que mide DebtRank.

### Ruta avanzada

6. Modifica `rondas_cascada` para que el umbral de cada nodo sea aleatorio,
   con distribución uniforme entre 0.1 y 0.3. ¿Cómo cambia la ventana?
7. Repite la sección 4 con un apalancamiento distinto para cada banco
   (por ejemplo, ``\ell_i`` uniforme entre 0.3 y 1.2). Calcula el radio
   espectral de ``\Lambda`` y compara con lo que pasa ante un choque del 1%.
8. Demuestra que ``\Phi`` de la sección 5 es creciente y que la iteración
   desde ``p = \bar p`` es una sucesión decreciente. ¿Por qué converge?

## Para la próxima sesión

- Lectura: Battiston et al. (2016). Complementaria: Watts (2002). Para la ruta
  avanzada, Eisenberg y Noe (2001) y Battiston et al. (2012); de consulta,
  Acemoglu et al. (2015) y Glasserman y Young (2016).
- En la sesión 13 los nodos dejan de ser sólo "caídos" o "en pie": cada
  agente **elige** entre dos opciones influido por sus vecinos. Es el modelo
  de Ising. Las cascadas reaparecerán en la ruta avanzada como
  **avalanchas**: cambios de opinión que se arrastran unos a otros.
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
    using PlutoUI, Statistics, LinearAlgebra, Random, Distributions, StatsBase, Graphs, Plots, Markdown
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
    const GRIS       = "#52514e"   # referencias y texto
    const GRIS_CLARO = "#c9c7c0"   # elementos de fondo
    const BLANCO     = "#ffffff"
end;

#%%
"""
Cascada de umbral a partir de `semilla`: un nodo cae cuando la fracción de
sus vecinos caídos es al menos `φ`. Devuelve la ronda en que cayó cada nodo
(-1 si no cayó).
"""
function rondas_cascada(g, φ, semilla)
    ronda = fill(-1, nv(g))
    ronda[semilla] = 0
    caidos = zeros(Int, nv(g))
    frente = [semilla]
    r = 0
    while !isempty(frente)
        r += 1
        nuevo = Int[]
        for v in frente, u in neighbors(g, v)
            ronda[u] >= 0 && continue
            caidos[u] += 1
            if caidos[u] >= φ * degree(g, u)
                ronda[u] = r
                push!(nuevo, u)
            end
        end
        frente = nuevo
    end
    ronda
end

#%%
"Red y nodo inicial del ejemplo de la sección 0 (una cascada que llega lejos)."
ejemplo_cascada = let
    g = erdos_renyi(200, 300; seed = 11)
    tam = [count(>=(0), rondas_cascada(g, 0.18, s)) for s in 1:200]
    g, findfirst(>=(0.5 * maximum(tam)), tam)
end;

#%%
zs_ventana = 0.5:0.5:10.0;

#%%
"""
Para cada grado promedio `z`, `corridas` cascadas en redes nuevas de `n`
nodos, cada una desde un nodo al azar. Devuelve la probabilidad de una
cascada global (más del 5% de la red) y su tamaño promedio.
"""
function barrido_cascadas(φ, zs; generar, n = 5000, corridas = 100)
    res = map(zs) do z
        rng = Xoshiro(1)
        tam = map(1:corridas) do r
            g = generar(n, z, r)
            count(>=(0), rondas_cascada(g, φ, rand(rng, 1:n))) / n
        end
        glob = tam .> 0.05
        (mean(glob), any(glob) ? mean(tam[glob]) : NaN)
    end
    (probabilidad = first.(res), tamaño = last.(res))
end

#%%
"Bordes de la ventana de cascadas de Watts para Erdős-Rényi con umbral φ."
function frontera_watts(φ; zs = 0.01:0.01:20)
    K = floor(Int, 1 / φ)
    f(z) = sum(k * (k - 1) * pdf(Poisson(z), k) for k in 1:K) - z
    s = sign.(f.(zs))
    i = findall(diff(s) .!= 0)
    [(zs[j] + zs[j+1]) / 2 for j in i]
end

#%%
"""
Modelo de Gai y Kapadia: `n` bancos en una red dirigida al azar, cada uno
con `z` deudores en promedio. Quiebra un banco al azar; un banco quiebra si
la fracción de sus deudores quebrados por θ supera su capital K. Devuelve
la fracción de bancos quebrados.
"""
function cascada_bancos(n, z, K; θ = 0.2, rng = Xoshiro(1))
    p = z / (n - 1)
    deudores = [[j for j in 1:n if j != i && rand(rng) < p] for i in 1:n]
    acreedores = [Int[] for _ in 1:n]
    for i in 1:n, j in deudores[i]
        push!(acreedores[j], i)
    end
    s = rand(rng, 1:n)
    quebrado = falses(n)
    quebrado[s] = true
    cola = [s]
    perdidos = zeros(Int, n)
    while !isempty(cola)
        j = popfirst!(cola)
        for i in acreedores[j]
            quebrado[i] && continue
            perdidos[i] += 1
            if θ * perdidos[i] / length(deudores[i]) > K
                quebrado[i] = true
                push!(cola, i)
            end
        end
    end
    mean(quebrado)
end

#%%
"Probabilidad y tamaño del contagio interbancario para cada número promedio de deudores."
function barrido_bancos(K, zs; n = 1000, corridas = 200)
    res = map(zs) do z
        tam = [cascada_bancos(n, z, K; rng = Xoshiro(r)) for r in 1:corridas]
        glob = tam .> 0.05
        (mean(glob), any(glob) ? mean(tam[glob]) : NaN)
    end
    (probabilidad = first.(res), tamaño = last.(res))
end

#%%
"""
Sistema bancario sintético: tamaños de Pareto con media 1; cada banco presta
el 20% de sus activos a `k` bancos elegidos con probabilidad proporcional a
su tamaño, en proporción a ese tamaño. `A[i, j]` es lo que `i` le prestó a
`j`, y el capital es `E = préstamos / ℓ`.
"""
function sistema_bancario(n; ℓ = 0.75, k = 3, θ = 0.2, α = 1.5, rng = Xoshiro(1))
    a = rand(rng, Pareto(α, 1.0), n)
    a .*= n / sum(a)
    A = zeros(n, n)
    for i in 1:n
        otros = setdiff(1:n, i)
        js = sample(rng, otros, Weights(a[otros]), k; replace = false)
        A[i, js] = θ * a[i] .* a[js] ./ sum(a[js])
    end
    (activos = a, A = A, E = θ .* a ./ ℓ)
end

#%%
"DebtRank lineal (Bardoscia et al., 2015) a partir del deterioro inicial `h0`."
function debtrank(b, h0; rondas = 2000)
    Λ = b.A ./ b.E
    h, antes = copy(h0), zeros(length(h0))
    for _ in 1:rondas
        nuevo = min.(1.0, h .+ Λ * (h .- antes))
        antes, h = h, nuevo
        maximum(abs.(h .- antes)) < 1e-12 && break
    end
    h
end

#%%
"Pérdida del sistema (fracción del capital total) en cada ronda de DebtRank."
function historia_debtrank(b, h0; rondas = 40)
    Λ = b.A ./ b.E
    h, antes = copy(h0), zeros(length(h0))
    hist = [perdida(b, h)]
    for _ in 1:rondas
        nuevo = min.(1.0, h .+ Λ * (h .- antes))
        antes, h = h, nuevo
        push!(hist, perdida(b, h))
    end
    hist
end

#%%
"Pérdida de capital del sistema, como fracción del capital total."
perdida(b, h) = sum(b.E .* h) / sum(b.E)

#%%
"DebtRank del banco j: capital perdido por los demás bancos si j quiebra."
function debtrank_de(b, j)
    h0 = zeros(length(b.E))
    h0[j] = 1.0
    h = debtrank(b, h0)
    (sum(b.E .* h) - b.E[j]) / sum(b.E)
end

#%%
"""
Vector de pagos de Eisenberg-Noe con deuda externa prioritaria `v`. `L[i, j]`
es lo que `i` le debe a `j` y `e` son los activos externos. Se itera desde el
pago completo. Un banco quiebra si no cubre su deuda externa más la
interbancaria.
"""
function compensacion(L, e, v; tol = 1e-12)
    p̄ = vec(sum(L; dims = 2))
    Π = L ./ p̄
    p = copy(p̄)
    for _ in 1:10_000
        nuevo = clamp.(e .- v .+ Π' * p, 0.0, p̄)
        maximum(abs.(nuevo .- p)) < tol && (p = nuevo; break)
        p = nuevo
    end
    recursos = e .+ Π' * p
    (pagos = p, quiebra = recursos .< v .+ p̄ .- 1e-9)
end

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
"Dibuja una red con posiciones `X` (n × 2)."
function dibujar_red(g, X; tam = 4, colores = AZUL, alpha_enlaces = 0.4, kwargs...)
    xs, ys = Float64[], Float64[]
    for e in edges(g)
        append!(xs, [X[src(e), 1], X[dst(e), 1], NaN])
        append!(ys, [X[src(e), 2], X[dst(e), 2], NaN])
    end
    p = plot(xs, ys; color = GRIS, alpha = alpha_enlaces, lw = 1, label = "", axis = false,
             grid = false, ticks = false, aspect_ratio = :equal, kwargs...)
    scatter!(p, X[:, 1], X[:, 2]; ms = tam, color = colores, msw = 0.5, msc = BLANCO, label = "")
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
