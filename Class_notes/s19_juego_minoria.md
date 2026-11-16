# Sesión 19: Juego de la minoría. Notas formales

> **Cómo usar estas notas.** La clase combina el cuaderno
> (`Code/s19_juego_minoria.jl`) y el pizarrón. El **plan de clase** de abajo
> dice qué se proyecta, qué secciones de estas notas se escriben en el
> pizarrón, cuándo se trabaja en parejas y qué queda para leer. En el
> pizarrón van las definiciones, las proposiciones y los resultados en
> recuadro; las demostraciones, sólo donde el plan lo indica. El texto que las
> rodea se explica de palabra. Como todo está aquí, en clase no hace falta
> copiar. El cuerpo de las notas es la **ruta básica**. Las demostraciones y
> los temas de la ruta avanzada, opcionales, están al final, en **Para
> profundizar**, y no se ven en clase; donde estaban queda una línea que los
> enlaza.
>
> Repaso recomendado: [P.B2](repaso/probabilidad.md#pb2) (media y varianza de
> una suma de variables independientes),
> [E.B3](repaso/procesos_estocasticos.md#eb3) (estacionariedad y promedios
> temporales), [F.B4](repaso/fisica_estadistica.md#fb4) (transiciones de fase:
> parámetro de control y de orden) y [S.I4](repaso/estadistica.md#si4) (datos
> dependientes y promedios sobre realizaciones); para la ruta avanzada,
> [F.A2](repaso/fisica_estadistica.md#fa2) (exponentes críticos y
> universalidad) y [F.A4](repaso/fisica_estadistica.md#fa4) (vidrios de espín y
> frustración).
>
> Lecturas: Slanina (2014), sec. 5.1. Complementaria: Challet & Zhang (1997);
> Requião da Cunha (2021), sec. 7.1.1; Slanina (2014), sec. 5.3; Arthur
> (1994). Ruta avanzada: Slanina
> (2014), sec. 5.2; Savit, Manuca & Riolo (1999); Challet & Marsili (1999);
> Challet, Marsili & Zecchina (2000); de consulta: Challet, Marsili & Zhang
> (2005).

**Plan de clase** (120 minutos)

| Minutos | Dónde | Qué | Secciones |
|---|---|---|---|
| 0–10 | Proyector | El fenómeno: el bar El Farol y la trampa de pensar igual | Cuaderno, sección 0 |
| 10–18 | Pizarra | El problema de El Farol: por qué no puede haber una regla predictiva compartida; la racionalidad inductiva | 1.1 |
| 18–33 | Pizarra | Las reglas del juego de la minoría: acción agregada y minoría ganadora; historias, $P = 2^M$; estrategias y puntajes virtuales | 1.2 |
| 33–40 | Pizarra | El juego como modelo de mercado: exceso de demanda y por qué gana la minoría | 2.1 |
| 40–50 | Pizarra | Repaso: la varianza de una suma; el número de ganadores y el referente aleatorio $\sigma^2/N = 1$, con la demostración | 2.2 |
| 50–100 | Parejas | Laboratorio: del efecto manada a la coordinación espontánea, y la curva $\sigma^2/N$ contra $\alpha$ para varios $N$; ejercicios de código 1 y 2 | Cuaderno, secciones 2 y 3 y ejercicios 1–2 |
| 100–112 | Pizarra | Lo que mostró el laboratorio: el colapso con $\alpha = 2^M/N$ y las tres fases | 3.1–3.2 |
| 112–120 | Pizarra | Cierre: puesta en común de las preguntas 3 a 5 | Cuaderno, preguntas 3–5 |
| Después | Lectura | El desorden congelado (1.2), el resumen (5) y Para profundizar | |

**Pregunta de la sesión.** Si conviene estar en la minoría, ¿pueden muchos
agentes que aprenden coordinarse mejor que el azar?

---

## 1. Del problema del bar El Farol al juego de la minoría

### 1.1 El problema del bar El Farol y el razonamiento inductivo

En los juegos de coordinación y en el modelo de Ising (sesiones 13, 14 y 17), a
cada agente le convenía hacer lo mismo que la mayoría. Sin embargo, en muchos
contextos económicos —congestión vial, aprovechamiento de recursos escasos o
especulación financiera— **conviene estar en la minoría**: cuando todos quieren
comprar el mismo activo al mismo tiempo, el precio se encarece; cuando todos
acuden al mismo lugar, se satura.

Arthur (1994) planteó esta tensión con el **problema del bar El Farol**:

- Cada jueves por la noche, $N = 100$ personas deciden de manera independiente
  si asisten o no al bar El Farol en Santa Fe.
- El bar sólo resulta agradable si asisten a lo sumo $60$ personas. Si van más
  de $60$, el lugar está demasiado lleno y quienes fueron habrían preferido
  quedarse en casa; si van $60$ o menos, quienes fueron disfrutan la noche y
  quienes se quedaron en casa lamentan no haber ido.

**Imposibilidad de una regla predictiva compartida.** Supongamos que existiera
un modelo deductivo único y público que predijera con exactitud la asistencia
$A_{\text{pred}}(t)$ de esta semana:

- Si el modelo predice $A_{\text{pred}}(t) \le 60$, todos los agentes racionales
  deciden ir, por lo que asisten $100 > 60$ y la predicción se refuta a sí
  misma.
- Si el modelo predice $A_{\text{pred}}(t) > 60$, nadie va, por lo que asisten
  $0 \le 60$ y la predicción vuelve a fallar.

Ante la ausencia de un predictor deductivo común, Arthur propuso que los
agentes emplean **racionalidad acotada e inductiva**: cada persona posee un
pequeño repertorio heterogéneo de reglas heurísticas que mapean la historia
reciente de asistencias en una predicción, lleva la cuenta de qué tan bien ha
funcionado cada regla y actúa según la más exitosa en ese momento. Challet,
Marsili & Ottino (2004) mostraron que con reglas de este tipo (tablas sobre las
últimas $M$ semanas, puntuadas por sus aciertos, como en la sección 1.2) la
asistencia media converge al cupo, que es lo que muestra la sección 0 del
cuaderno.

### 1.2 Formulación del juego de la minoría

Challet & Zhang (1997) destilaron el problema de El Farol en un modelo
binario simétrico que puede analizarse con las herramientas de la mecánica
estadística: el **juego de la minoría** (*Minority Game*, MG).

**Agentes y acción ganadora.** Un número impar $N$ de agentes toma en cada paso
discreto $t = 1, 2, \dots$ una decisión binaria

$$
a_i(t) \in \lbrace -1, +1\rbrace, \qquad i = 1, \dots, N .
$$

La **acción agregada** es

$$
A(t) = \sum_{i=1}^N a_i(t) \in \lbrace -N, -N+2, \dots, N-2, N\rbrace .
$$

Como $N$ es impar, $A(t) \ne 0$ siempre. El grupo minoritario es el que eligió
el signo opuesto a $A(t)$, por lo que el **resultado ganador (el lado de la
minoría)** en el paso $t$ es

$$
a_{\text{min}}(t) = -\operatorname{sgn}\bigl(A(t)\bigr) \in \lbrace -1, +1\rbrace .
$$

El agente $i$ gana en el paso $t$ si y sólo si $a_i(t) = a_{\text{min}}(t)$, es
decir, si $a_i(t) A(t) < 0$.

**Memoria e historias.** Todos los agentes tienen acceso a la misma información
pública: la cadena formada por los últimos $M$ resultados ganadores,

$$
\mu(t) = \bigl(a_{\text{min}}(t - M),\dots, a_{\text{min}}(t - 1)\bigr) .
$$

Como cada entrada toma dos valores posibles, existen exactamente

$$
P = 2^M
$$

historias distintas, que indexamos con un entero
$\mu(t) \in \lbrace 1, 2, \dots, P\rbrace$.

**Estrategias y aprendizaje por refuerzo.** Una **estrategia** $s$ es una
función (tabla de consulta) que asigna a cada una de las $P = 2^M$ historias
posibles una acción recomendada:

$$
s: \lbrace 1, \dots, P\rbrace \longrightarrow \lbrace -1, +1\rbrace,
\qquad \mu \longmapsto s(\mu) .
$$

En total existen $2^P = 2^{2^M}$ estrategias posibles. Antes de comenzar el
juego ($t = 0$), cada agente $i$ recibe $S$ estrategias fijas
$s_{i,1}, \dots, s_{i,S}$, sorteadas al azar de manera independiente y uniforme
entre las $2^P$ tablas posibles (el caso estándar es $S = 2$).

Cada agente mantiene un **puntaje virtual** $U_{i,s}(t)$ para cada una de sus
$S$ estrategias (con $U_{i,s}(0) = 0$). En el paso $t$:

1. Observada la historia $\mu = \mu(t)$, cada agente elige su estrategia activa
   $s_i^\ast(t) = \arg\max_{s \in \lbrace 1,\dots,S\rbrace} U_{i,s}(t)$ (si hay
   empate, elige al azar entre las empatadas) y ejecuta la acción

   $$a_i(t) = s_{i, s_i^\ast(t)}\bigl(\mu(t)\bigr) .$$

2. Se calcula la acción agregada $A(t) = \sum_{i=1}^N a_i(t)$ y todos los
   agentes actualizan los puntajes de **todas** sus $S$ estrategias (tanto la
   que usaron como las que no usaron):

   $$U_{i,s}(t + 1) = U_{i,s}(t) - s_{i,s}\bigl(\mu(t)\bigr)\thinspace\frac{A(t)}{N}$$

   (o de forma equivalente sumando $-s_{i,s}(\mu(t)) A(t)$, ya que multiplicar
   todos los puntajes por $N$ no altera cuál estrategia tiene el máximo). Si la
   estrategia $s$ recomendaba el lado minoritario, $-s(\mu(t))A(t) > 0$ y su
   puntaje sube; si recomendaba la mayoría, su puntaje baja.
3. La historia avanza incorporando
   $a_{\text{min}}(t) = -\operatorname{sgn} A(t)$.

> **Desorden congelado (*quenched disorder*).** Las tablas de estrategias
> $s_{i,s}(\mu)$ quedan fijas desde $t = 0$ para cada agente. Por tanto, para
> caracterizar las propiedades estadísticas del modelo hay que promediar tanto
> en el tiempo (una vez alcanzado el estado estacionario,
> [E.B3](repaso/procesos_estocasticos.md#eb3)) como sobre muchas realizaciones
> independientes del sorteo inicial de estrategias
> ([S.I4](repaso/estadistica.md#si4)).

---

## 2. El juego de la minoría como modelo de mercado

### 2.1 Exceso de demanda, impacto en el precio y ganancia de la minoría

En un mercado financiero estilizado, en cada periodo $t$ el agente $i$ envía
una orden de mercado de tamaño unitario:

- $a_i(t) = +1$: orden de **compra**;
- $a_i(t) = -1$: orden de **venta**.

La suma $A(t) = \sum_{i=1}^N a_i(t) = N_{\text{compra}}(t) - N_{\text{venta}}(t)$
es el **exceso de demanda** (o desbalance de órdenes). Un formador de mercado
(*market maker*) absorbe el desbalance ajustando el precio en proporción a
$A(t)$ dividida por la liquidez $\lambda_M > 0$:

$$
r(t) = \ln p(t + 1) - \ln p(t) = \frac{A(t)}{\lambda_M} ,
$$

donde escribimos el precio $p(t)$ en minúscula para no confundirlo con el
número de historias $P$.

¿Por qué el pago de corto plazo tiene la estructura de minoría $-a_i(t) A(t)$?
Cuando hay exceso de compradores ($A(t) > 0$), la presión compradora empuja el
precio de ejecución hacia arriba: los muchos compradores adquieren el activo
encarecido, mientras que los pocos vendedores (la minoría, $a_i = -1$) logran
vender a un precio elevado. De manera simétrica, cuando casi todos venden
($A(t) < 0$), el precio se desploma y la minoría compradora adquiere el activo
con descuento.

### 2.2 Eficiencia colectiva, volatilidad y el referente aleatorio

Si en el paso $t$ exactamente $k$ agentes eligen $+1$ y $N - k$ eligen $-1$,
entonces $A(t) = 2k - N$ y el número de agentes en el grupo minoritario (los
ganadores del turno) es

$$
N_{\text{gan}}(t) = \min(k, N - k) = \frac{N - \lvert A(t)\rvert}{2} .
$$

El máximo número posible de ganadores en un turno es $(N - 1)/2$, que se
alcanza cuando $\lvert A(t)\rvert = 1$ (la población se divide casi mitad y
mitad). En cambio, si todos hacen lo mismo ($\lvert A(t)\rvert = N$),
$N_{\text{gan}}(t) = 0$: nadie queda en la minoría y el recurso se desperdicia
por completo.

Por tanto, cuanto menores son las fluctuaciones de $A(t)$ alrededor de $0$,
mayor es el número medio de ganadores $\langle N_{\text{gan}}\rangle$ y menor
es la volatilidad del mercado. La medida natural de las fluctuaciones (y de la
**ineficiencia asignativa**) es la **volatilidad reescalada**:

$$
\frac{\sigma^2}{N} = \frac{\bigl\langle A(t)^2\bigr\rangle}{N} .
$$

**Proposición** (referente aleatorio). Si los $N$ agentes eligen
$a_i(t) \in \lbrace -1, +1\rbrace$ al azar de manera independiente y simétrica
($\mathbb P(a_i = +1) = \mathbb P(a_i = -1) = 1/2$), entonces

$$
\langle A\rangle = 0, \qquad
\sigma_{\text{azar}}^2 = \operatorname{Var}(A) = N
\quad\Longrightarrow\quad
\boxed{\ \frac{\sigma_{\text{azar}}^2}{N} = 1\ }
$$

*Demostración.* Cada $a_i$ tiene media $\mathbb E[a_i] = 0$ y varianza
$\operatorname{Var}(a_i) = \mathbb E[a_i^2] = 1$. Por independencia entre los
$N$ agentes ([P.B2](repaso/probabilidad.md#pb2)),
$\operatorname{Var}\bigl(\sum_{i=1}^N a_i\bigr) = \sum_{i=1}^N \operatorname{Var}(a_i) = N$. $\square$

El valor $\sigma^2/N = 1$ divide el comportamiento colectivo en dos regímenes:

- **$\sigma^2 / N > 1$ (peor que el azar):** Los agentes se estorban entre sí;
  la correlación positiva entre sus decisiones genera sobre-reacción y efecto
  manada (*herding*).
- **$\sigma^2 / N < 1$ (mejor que el azar):** Las estrategias de los agentes se
  anticoordinan parcialmente; el sistema reparte el recurso con menor
  desperdicio que si lanzaran monedas.

---

## 3. Escalamiento con $\alpha = 2^M/N$ y las dos fases del modelo

### 3.1 El parámetro de control $\alpha = P/N = 2^M/N$

Savit, Manuca & Riolo (1999) descubrieron que al simular el juego de la minoría
para distintos tamaños de población $N$ y distintas memorias $M$, las curvas de
$\sigma^2/N$ contra $M$ se desplazan hacia la derecha al crecer $N$, pero
**colapsan sobre una única curva universal** cuando se grafican en función del
parámetro adimensional

$$
\boxed{\ \alpha = \frac{P}{N} = \frac{2^M}{N}\ }
$$

El parámetro de control $\alpha$ ([F.B4](repaso/fisica_estadistica.md#fb4))
mide el cociente entre la **complejidad de la información disponible** ($P = 2^M$
historias distinguibles) y el **número de agentes** $N$ que compiten usando esa
información.

### 3.2 Fase hacinada, punto óptimo y fase diluida

Para un número fijo de estrategias por agente ($S = 2$), la curva universal de
$\sigma^2/N$ contra $\alpha$ presenta un mínimo pronunciado en un valor crítico

$$
\alpha_c \approx 0.3374 ,
$$

que separa dos fases con propiedades opuestas:

1. **Fase simétrica o hacinada ($\alpha < \alpha_c$, es decir,
   $2^M < \alpha_c N \approx N/3$):** Hay pocas historias posibles frente al
   número de agentes. En
   consecuencia, el espacio de estrategias relevantes está **sobrepoblado**
   (*crowded*): muchos agentes poseen estrategias idénticas o muy parecidas,
   evalúan los mismos puntajes y cambian al mismo tiempo hacia la misma acción
   cuando una estrategia pasa a ser la mejor. Esta sincronización involuntaria
   provoca grandes avalanchas de órdenes del mismo signo ($A(t) \propto N$).
   Por tanto, $\operatorname{Var}(A) \propto N^2 / P$ y la volatilidad
   reescalada diverge como una ley de potencias al disminuir $\alpha$:

   $$\frac{\sigma^2}{N} \sim C\thinspace\alpha^{-1} \qquad (\alpha \to 0) .$$

   En esta región $\sigma^2/N \gg 1$: el aprendizaje colectivo produce una
   volatilidad mucho mayor que el azar.
2. **Coordinación óptima cerca del punto crítico,
   $\alpha \approx \alpha_c \approx 0.34$:** El número de historias $P = 2^M$ guarda una proporción
   intermedia con $N$. Los agentes logran diversificar sus acciones sin caer en
   el hacinamiento y la volatilidad cae hasta $\sigma^2/N \approx 0.2 < 1$,
   bastante por debajo del límite aleatorio.
3. **Fase asimétrica o diluida ($\alpha > \alpha_c$, es decir,
   $2^M > \alpha_c N$):** El número de historias $P$ es tan grande que las $S = 2$ tablas de
   cada agente son prácticamente ortogonales a las de los demás. Al crecer
   $\alpha$, las correlaciones entre agentes se diluyen y $\sigma^2/N$ sube
   monótonamente desde su mínimo hacia la asíntota aleatoria
   $\lim_{\alpha \to \infty} \sigma^2/N = 1$.

---

## 4. Teoría del punto crítico $\alpha_c$ y vidrios de espín [Avanzado]

El desarrollo de esta sección está en [Para profundizar, P.1](#pp1).


---

## 5. Resumen formal

| Concepto / Régimen | Fórmula / Resultado clave |
|---|---|
| Acción agregada y minoría ganadora | $A(t) = \sum_{i=1}^N a_i(t)$, $\quad a_{\text{min}}(t) = -\operatorname{sgn} A(t)$, $\quad N_{\text{gan}}(t) = (N - \lvert A(t)\rvert)/2$ |
| Actualización de puntajes virtuales | $U_{i,s}(t+1) = U_{i,s}(t) - s_{i,s}(\mu(t))\thinspace A(t)/N$; el agente juega $\arg\max_s U_{i,s}(t)$ |
| Referente de elección aleatoria | Si $a_i = \pm 1$ son i.i.d. con probabilidad $1/2$, entonces $\sigma_{\text{azar}}^2 / N = \langle A^2\rangle / N = 1$ |
| Parámetro de control universal | $\alpha = P / N = 2^M / N$, donde $M$ es la memoria y $P = 2^M$ el número de historias |
| Fase simétrica / hacinada ($\alpha < \alpha_c$) | Predictibilidad $H/N = 0$ (sin arbitraje), pero volatilidad gigante $\sigma^2/N \sim \alpha^{-1} > 1$ por efecto manada |
| Punto crítico ($\alpha_c \approx 0.3374$, $S = 2$) | Mínimo global de volatilidad ($\sigma^2/N \approx 0.2 < 1$); balance de grados de libertad $\alpha_c = 1 - \phi(\alpha_c)$ |
| Fase asimétrica / diluida ($\alpha > \alpha_c$) | Predictibilidad $H/N > 0$ (información sin explotar), fracción $\phi(\alpha) > 0$ de agentes congelados y $\sigma^2/N \to 1$ |

---

## 6. Ejercicios formales

**Ruta básica**

- **F1.** En un juego de la minoría con $N = 101$ agentes, en un turno dado $63$
  agentes eligen $+1$ (comprar) y $38$ eligen $-1$ (vender). Calcula la acción
  agregada $A(t)$, el signo de la acción ganadora $a_{\text{min}}(t)$ y el
  número de agentes ganadores $N_{\text{gan}}(t)$. ¿Qué fracción de los agentes
  ganó en ese turno?
- **F2.** Si $N = 201$ agentes juegan con memoria $M = 6$, ¿cuántas historias
  posibles $P$ existen, cuántas estrategias totales $2^P$ hay en el universo de
  tablas posibles y cuánto vale el parámetro de control $\alpha$? ¿En qué fase
  del modelo se encuentra el sistema?
- **F3.** Demuestra que si en cada turno los $N$ agentes eligen $+1$ con
  probabilidad $p \in (0, 1)$ y $-1$ con probabilidad $1 - p$ de forma
  independiente, entonces $\mathbb E[A] = N(2p - 1)$ y
  $\operatorname{Var}(A) = 4 N p(1 - p)$. Comprueba que la varianza es máxima
  en $p = 1/2$ ($\operatorname{Var}(A)/N = 1$), pero que
  $\langle A^2\rangle/N = 4p(1-p) + N(2p-1)^2$ crece como $\mathcal O(N)$ en cuanto $p \ne 1/2$.
- **F4.** Un agente tiene $S = 2$ estrategias para $M = 2$ ($P = 4$ historias).
  Ante la historia $\mu = 3$, la estrategia $1$ recomienda $s_1(3) = +1$ y la
  estrategia $2$ recomienda $s_2(3) = -1$. Si antes del turno sus puntajes eran
  $U_1 = 1.2$ y $U_2 = 1.0$, y la acción agregada del turno resulta ser
  $A(t) = +15$, ¿qué acción jugó el agente en ese turno, cuáles son los nuevos
  puntajes $U_1, U_2$ (usando la regla $\Delta U_s = -s(\mu)A$) y qué estrategia
  usará la próxima vez que aparezca la historia $\mu = 3$?

**Ruta avanzada**

- **F5.** Demuestra que para cualquier par de estrategias binarias
  $s_{i,1}(\mu), s_{i,2}(\mu) \in \lbrace -1, +1\rbrace$, la descomposición
  $\bar s_i(\mu) = \frac{1}{2}[s_{i,1}(\mu) + s_{i,2}(\mu)]$ y
  $\xi_i(\mu) = \frac{1}{2}[s_{i,1}(\mu) - s_{i,2}(\mu)]$ cumple para toda
  historia $\mu$ que $\bar s_i(\mu)\thinspace\xi_i(\mu) = 0$ y
  $\bar s_i(\mu)^2 + \xi_i(\mu)^2 = 1$. Si las dos tablas se sortean de manera
  independiente y uniforme, ¿con qué probabilidad ocurre $\xi_i(\mu) = 0$ en una
  historia dada $\mu$?
- **F6.** En el extremo $\alpha \ll \alpha_c$ (por ejemplo,
  $M = 1 \Rightarrow P = 2$ historias y $N \gg 1$), existen únicamente $2^P = 4$ tablas posibles:
  dos constantes ($+1,+1$ y $-1,-1$) y dos alternantes ($+1,-1$ y $-1,+1$).
  Explica por qué una fracción finita de los $N$ agentes recibe exactamente el
  mismo par de estrategias y deduce a partir de ello el escalamiento
  $\operatorname{Var}(A) \propto N^2$, es decir,
  $\sigma^2 / N \propto N \propto \alpha^{-1}$.
- **F7.** Usa la descomposición
  $a_i(t) = \bar s_i(\mu(t)) + \xi_i(\mu(t)) \sigma_i(t)$ en el límite $\alpha \to \infty$ ($P \gg N$), donde las
  correlaciones entre agentes desaparecen. Supón que las historias $\mu$ se
  visitan con probabilidad uniforme $1/P$ y que cada $\sigma_i(t) = \pm 1$ no
  depende de las tablas de los demás agentes. Calcula el valor esperado sobre el
  desorden de estrategias de $\langle A^2\rangle / N$ y verifica que coincide
  con el referente aleatorio $1$. (Pista: usa F5 para obtener
  $\mathbb E[a_i^2] = 1$.)
- **F8.** Verifica explícitamente que la derivada parcial de la predictibilidad
  $H(m_1, \dots, m_N) = \frac{1}{P}\sum_{\mu=1}^P \bigl[\Omega(\mu) + \sum_{j=1}^N \xi_j(\mu) m_j\bigr]^2$ respecto a $m_i$ es proporcional a la deriva media del
  puntaje diferencial
  $\langle \Delta q_i\rangle = -\frac{1}{PN}\sum_{\mu=1}^P \xi_i(\mu)\langle A \mid \mu\rangle$, demostrando que $H$ es una función de
  Lyapunov de la dinámica de aprendizaje.

---

## Para profundizar

Material opcional de la ruta avanzada: demostraciones y temas que no se ven
en clase. Cada apartado indica de qué sección de las notas viene.

<a id="pp1"></a>

### P.1 Teoría del punto crítico $\alpha_c$ y vidrios de espín (sección 4)


#### 4.1 Predictibilidad $H$ como parámetro de orden

¿Por qué el mínimo de volatilidad en $\alpha_c \approx 0.3374$ corresponde a
una verdadera **transición de fase** en el sentido de la física estadística?
Para verlo, examinamos si la serie de excesos de demanda $A(t)$ contiene
información predecible condicionada en la historia $\mu$.

Sea $\langle A \mid \mu\rangle$ el promedio estacionario de $A(t)$ en los pasos
en que la historia observada es $\mu \in \lbrace 1, \dots, P\rbrace$. Definimos
la **predictibilidad** (reescalada por $N$) como

$$
\frac{H}{N} = \frac{1}{P N}\sum_{\mu=1}^P \langle A \mid \mu\rangle^2 .
$$

- Si $\langle A \mid \mu\rangle = 0$ para toda historia $\mu$, entonces
  $H/N = 0$: conociendo la historia de longitud $M$, un observador externo no
  puede predecir el signo medio de $A(t)$. El mercado es **informacionalmente
  eficiente (fase simétrica)**.
- Si $\langle A \mid \mu\rangle \ne 0$ para algunas historias $\mu$, entonces
  $H/N > 0$: después de ciertas historias el mercado tiende sistemáticamente a
  comprar más que a vender (o viceversa), dejando **oportunidades de ganancia
  sin explotar (fase asimétrica)**.

El comportamiento numérico y analítico muestra que $H/N$ actúa como un
**parámetro de orden** ([F.B4](repaso/fisica_estadistica.md#fb4),
[F.A2](repaso/fisica_estadistica.md#fa2)):

$$
\frac{H}{N} = 0 \quad \text{para todo } \alpha \le \alpha_c,
\qquad
\frac{H}{N} \sim c_H\thinspace(\alpha - \alpha_c)^2 > 0 \quad \text{para } \alpha > \alpha_c .
$$

#### 4.2 Descomposición de estrategias y agentes congelados

Para $S = 2$, cada agente $i$ tiene dos estrategias
$s_{i,1}(\mu), s_{i,2}(\mu) \in \lbrace -1, +1\rbrace$, que podemos reescribir separando su parte común y
su semidiferencia:

$$
\bar s_i(\mu) = \frac{s_{i,1}(\mu) + s_{i,2}(\mu)}{2} \in \lbrace -1, 0, +1\rbrace,
\qquad
\xi_i(\mu) = \frac{s_{i,1}(\mu) - s_{i,2}(\mu)}{2} \in \lbrace -1, 0, +1\rbrace .
$$

Cuando las dos estrategias de $i$ coinciden ante la historia $\mu$,
$\xi_i(\mu) = 0$ y el agente toma la acción fija $\bar s_i(\mu) = \pm 1$ sin
importar cuál de sus dos estrategias tenga mayor puntaje. Cuando difieren,
$\bar s_i(\mu) = 0$, $\xi_i(\mu) = \pm 1$ y la acción del agente es
$a_i(t) = \xi_i(\mu(t))\thinspace \sigma_i(t)$, donde $\sigma_i(t) = +1$ si el
agente usa su estrategia $1$ y $\sigma_i(t) = -1$ si usa la $2$.

Sea $m_i = \langle \sigma_i(t)\rangle \in [-1, +1]$ la preferencia temporal
media del agente $i$ entre sus dos estrategias. Promediando sobre el tiempo
condicionado en $\mu$,

$$
\langle A \mid \mu\rangle = \Omega(\mu) + \sum_{i=1}^N \xi_i(\mu)\thinspace m_i,
\qquad
\text{donde } \Omega(\mu) = \sum_{i=1}^N \bar s_i(\mu) .
$$

La semidiferencia de puntajes $q_i(t) = \frac{1}{2}[U_{i,1}(t) - U_{i,2}(t)]$
cambia en cada paso en $\Delta q_i(t) = -\xi_i(\mu(t)) A(t)/N$. Si todas las
historias $\mu$ se visitan con frecuencia $\approx 1/P$, la deriva media de
$q_i$ es

$$
\langle \Delta q_i\rangle \approx -\frac{1}{P N}\sum_{\mu=1}^P \xi_i(\mu)\thinspace\langle A \mid \mu\rangle .
$$

De aquí surgen dos tipos de agentes en el estado estacionario:

1. **Agentes congelados (*frozen agents*, $\lvert m_i\rvert = 1$):** Si
   $\langle \Delta q_i\rangle \ne 0$, la diferencia de puntajes $q_i(t)$ crece
   linealmente con el tiempo hacia $+\infty$ o $-\infty$. Una de las dos
   estrategias domina permanentemente a la otra y el agente deja de alternar:
   queda "congelado" en $\sigma_i(t) = +1$ ($m_i = +1$) o en $\sigma_i(t) = -1$
   ($m_i = -1$).
2. **Agentes volátiles o libres ($\lvert m_i\rvert < 1$):** Si
   $\langle \Delta q_i\rangle = 0$, la diferencia de puntajes $q_i(t)$ fluctúa
   alrededor de cero y el agente sigue alternando entre sus dos estrategias.

#### 4.3 Conteo de grados de libertad

Challet, Marsili & Zecchina (2000) mostraron que el estado estacionario del
juego minimiza la predictibilidad

$$
H(m_1, \dots, m_N) = \frac{1}{P}\sum_{\mu=1}^P \left(\Omega(\mu) + \sum_{i=1}^N \xi_i(\mu) m_i\right)^2
$$

sobre el hipercubo $[-1, +1]^N$: los agentes libres anulan su deriva
$\langle \Delta q_i\rangle = -\frac{1}{2N}\frac{\partial H}{\partial m_i}$ y los
congelados quedan en la frontera. Como $H$ es cuadrática en las $m_i$ con
coeficientes aleatorios, es un hamiltoniano de **vidrio de espín**
([F.A4](repaso/fisica_estadistica.md#fa4)); el análisis completo con esas
herramientas está en Challet, Marsili & Zhang (2005). Aquí basta un argumento
de conteo.

**Resultado (conteo de restricciones)** (Challet & Marsili, 1999;
Challet, Marsili & Zecchina, 2000). Sea $\phi(\alpha) \in [0, 1]$ la fracción
de agentes congelados ($\lvert m_i\rvert = 1$). El punto crítico $\alpha_c$
ocurre exactamente cuando la fracción de agentes libres $1 - \phi(\alpha)$
iguala al parámetro $\alpha = P/N$:

$$
1 - \phi(\alpha_c) = \alpha_c \qquad (\alpha_c \approx 0.3374 \text{ para } S = 2) .
$$

*Argumento.* Hacer $H = 0$ equivale a satisfacer las $P$ ecuaciones lineales

$$
\sum_{i=1}^N \xi_i(\mu) m_i = -\Omega(\mu), \qquad \mu = 1, \dots, P ,
$$

sujetas a la restricción de caja $-1 \le m_i \le +1$. Las variables que están
clavadas en la frontera $m_i = \pm 1$ (los $\phi N$ agentes congelados) no
pueden ajustarse libremente; las verdaderas incógnitas continuas son las
preferencias $m_i \in (-1, +1)$ de los $(1 - \phi)N$ agentes libres:

- Cuando $\alpha > \alpha_c$, el número de ecuaciones $P = \alpha N$ es mayor
  que el número de variables libres $(1 - \phi)N$: el sistema lineal está
  **sobredeterminado** y no existe solución con $H = 0$, por lo que el mínimo
  global tiene $H > 0$. Al disminuir $\alpha$ hacia $\alpha_c$, más agentes
  tocan la frontera $m_i = \pm 1$ y la fracción de congelados $\phi(\alpha)$
  aumenta.
- En el instante en que $P = (1 - \phi(\alpha_c))N$, es decir,
  $\alpha_c = 1 - \phi(\alpha_c)$, el número de ecuaciones iguala al número de
  variables libres y el mínimo alcanza $H = 0$.
- Para $\alpha < \alpha_c$, como $H = 0$ implica $\langle A \mid \mu\rangle = 0$
  para todo $\mu$, la deriva $\langle \Delta q_i\rangle$ se anula para todos los
  agentes, la fracción de congelados cae a $\phi = 0$ y el sistema queda
  **subdeterminado** ($P < N$), con un subespacio de dimensión $N - P$ de
  soluciones que dan $H = 0$. Las fluctuaciones a lo largo de esas direcciones
  no restringidas disparan la volatilidad $\sigma^2/N \sim \alpha^{-1}$.
  $\square$

---

## Referencias

- Arthur, W. B. (1994). Inductive reasoning and bounded rationality. *American
  Economic Review*, 84(2), 406–411. https://www.jstor.org/stable/2117868
- Challet, D., & Marsili, M. (1999). Phase transition and symmetry breaking in
  the minority game. *Physical Review E*, 60(6), R6271–R6274.
  https://doi.org/10.1103/PhysRevE.60.R6271
- Challet, D., Marsili, M., & Ottino, G. (2004). Shedding light on El Farol.
  *Physica A*, 332, 469–482. https://doi.org/10.1016/j.physa.2003.06.003
- Challet, D., Marsili, M., & Zhang, Y.-C. (2000). Modeling market mechanism
  with minority game. *Physica A*, 276(1–2), 284–315.
  https://doi.org/10.1016/S0378-4371(99)00446-X
- Challet, D., Marsili, M., & Zecchina, R. (2000). Statistical mechanics of
  systems with heterogeneous agents: Minority games. *Physical Review Letters*,
  84(8), 1824–1827. https://doi.org/10.1103/PhysRevLett.84.1824
- Challet, D., Marsili, M., & Zhang, Y.-C. (2005). *Minority Games: Interacting
  Agents in Financial Markets*. Oxford University Press.
- Challet, D., & Zhang, Y.-C. (1997). Emergence of cooperation and organization
  in an evolutionary game. *Physica A*, 246(3–4), 407–418.
  https://doi.org/10.1016/S0378-4371(97)00419-6
- Coolen, A. C. C. (2005). *The Mathematical Theory of Minority Games:
  Statistical Mechanics of Interacting Agents*. Oxford University Press.
- Moro, E. (2004). The minority game: An introductory guide. En *Advances in
  Condensed Matter and Statistical Mechanics* (pp. 263–286). Nova Science
  Publishers. https://arxiv.org/abs/cond-mat/0402651
- Requião da Cunha, C. (2021). *Introduction to Econophysics: Contemporary
  Approaches with Python Simulations*. CRC Press.
  https://doi.org/10.1201/9781003127956
- Savit, R., Manuca, R., & Riolo, R. (1999). Adaptive competition, market
  efficiency, and phase transitions. *Physical Review Letters*, 82(10),
  2203–2206. https://doi.org/10.1103/PhysRevLett.82.2203
- Slanina, F. (2014). *Essentials of Econophysics Modelling*. Oxford
  University Press.
