# Sesión 16: Schelling: segregación. Notas formales

> **Cómo leer estas notas.** Siguen el orden del programa: primero se observa
> en el cuaderno (`Code/s16_schelling.jl`), después se formaliza aquí. Lo que no
> tiene marca es **ruta básica**. Las secciones y demostraciones marcadas
> **[Avanzado]** son opcionales.
>
> Repaso recomendado: [F.B4](repaso/fisica_estadistica.md#fb4) (transiciones de
> fase) y [F.I3](repaso/fisica_estadistica.md#fi3) (modelo de Ising).
>
> Lecturas: Kirman (2011), cap. 6. Complementaria: Sen & Chakrabarti (2014),
> cap. 4; Schelling (1971); Wilensky & Rand (2015, consulta). Ruta avanzada:
> Slanina (2014), sec. 8.2; Sen & Chakrabarti (2014), cap. 6; Vinković &
> Kirman (2006); Grauwin et al. (2009).

---

## 1. El modelo de segregación de Schelling

### 1.1 Reglas, similitud y satisfacción

Schelling (1971) estudió cómo preferencias individuales moderadas por no ser
una minoría extrema en el propio vecindario generan una fuerte segregación
macroscópica.

**Configuración.** En una cuadrícula de $L \times L$ celdas con condiciones de
frontera periódicas, una fracción $\rho_0 \in (0, 1)$ de las celdas permanece
vacía ($c_i = 0$) y las celdas restantes están ocupadas en partes iguales por
agentes de dos grupos, $c_i \in \lbrace 1, 2\rbrace$.

**Vecindad y similitud local.** Cada agente en la celda $i$ observa su vecindad
de Moore $\mathcal N_i$ (las $8$ celdas adyacentes). Sea $n_i = \sum_{j \in
\mathcal N_i} \mathbf 1(c_j \ne 0)$ el número de vecinos ocupados y $n_i^{=} =
\sum_{j \in \mathcal N_i} \mathbf 1(c_j = c_i)$ el número de vecinos del mismo
grupo. La **similitud local** del agente $i$ es

$$
s_i = \frac{n_i^{=}}{n_i} \qquad (\text{con } s_i = 1 \text{ si } n_i = 0) .
$$

**Regla de mudanza.** Todos los agentes comparten un **umbral de tolerancia**
$\tau \in [0, 1]$:

- El agente $i$ está **satisfecho** si $s_i \ge \tau$.
- Si $s_i < \tau$, está **insatisfecho** y se muda en la ronda siguiente a una
  celda vacía elegida al azar.

La **similitud media** $\bar s = \frac{1}{N_{\text{occ}}}\sum_{c_i \ne 0} s_i$
mide el grado de segregación espacial: en una distribución aleatoria con dos
grupos iguales, $\mathbb E[\bar s] = 1/2$; cuando se forman dominios puros de
cada grupo, $\bar s$ se acerca a $1$ (salvo por los agentes en las fronteras
entre dominios).

---

## 2. Umbral crítico, cascadas de mudanza y frustración

Como una vecindad de Moore tiene a lo sumo $8$ vecinos, el umbral $\tau$
opera en escalones discretos ($1/8 = 0.125$, $2/8 = 0.25$, $3/8 = 0.375$,
etc., más las fracciones con denominadores $n_i < 8$ cuando hay celdas
vacías). Al variar $\tau$ aparecen tres regímenes bien diferenciados
([F.B4](repaso/fisica_estadistica.md#fb4)):

1. **Tolerancia alta ($\tau \le 1/4 = 0.25$):** En la configuración aleatoria
   inicial, la probabilidad binomial de que un agente con $7$ u $8$ vecinos
   tenga menos de $2$ del mismo grupo es pequeña. Pocos agentes se mudan y la
   ciudad permanece integrada ($\bar s$ entre 0.52 y 0.60 en el barrido del
   cuaderno).
2. **Segregación emergente ($0.3 \le \tau \le 0.7$ en el barrido, en
   particular $\tau = 3/8 = 0.375$ y $\tau = 1/2$):** Con celdas vacías, un
   agente puede tener menos de $8$ vecinos ocupados, y fracciones como
   $2/7 \approx 0.29$ ya no cumplen $\tau = 0.3$. Aunque ningún agente desea una
   ciudad segregada —cada uno acepta hasta $5/8$ (62.5%) de vecinos del otro
   grupo cuando $\tau = 3/8$—, cada mudanza disminuye la fracción del propio
   grupo en el vecindario abandonado y la aumenta en el de llegada. Esta
   **externalidad de vecindario** desencadena una cascada de reubicaciones que
   sólo se detiene cuando la ciudad se divide en grandes dominios homogéneos
   con $\bar s$ entre 0.74 y casi 1 (en el barrido del cuaderno) y
   satisfacción final del 100%.
3. **Intolerancia extrema y frustración ($\tau \gtrsim 3/4 = 0.75$):** Casi
   ninguna frontera entre los dos grupos permite que ambos lados tengan
   simultáneamente $s_i \ge \tau$. El sistema no encuentra una configuración
   absorbente y queda congelado o ciclando en un estado desordenado con baja
   satisfacción.

---

## 3. Analogía física: Ising a temperatura cero con magnetización conservada

El modelo de Schelling es el análogo fuera del equilibrio de un **modelo de
Ising ferromagnético con espín conservado** ([F.I3](repaso/fisica_estadistica.md#fi3);
Vinković y Kirman, 2006), conocido en física como **dinámica de Kawasaki**
(dos celdas intercambian sus estados en lugar de voltear un espín individual,
pues el número de agentes de cada grupo es fijo):

- **Energía de interfaz.** En un modelo de Ising con estados $\pm 1$ y vacíos
  $0$, la energía de enlaces entre vecinos distintos es proporcional al
  perímetro total de las fronteras entre el grupo $1$ y el grupo $2$. Mudarse
  hacia vecindarios con $s_i \ge \tau$ suele reducir el número de enlaces entre
  grupos distintos (aunque una mudanza concreta puede aumentarlo), engrosando
  los dominios (*coarsening*) como un ferromagneto enfriado por debajo de su
  temperatura crítica ($T \to 0$).
- **[Avanzado] Falta de función potencial exacta y modelo de bloques.** A
  diferencia del modelo de Ising estándar ([sesión 13](s13_logit_ising.md)), la
  utilidad individual de Schelling es una función escalón
  $u_i = \mathbf 1(s_i \ge \tau)$, que depende del cociente $n_i^{=}/n_i$ y
  cuya ganancia al mudarse no es exactamente igual al cambio de una energía
  global. Grauwin et al. (2009) resuelven exactamente una variante de Schelling
  (donde la utilidad depende de la composición del bloque donde vive el agente)
  en la que un parámetro interpola entre dinámica individual y cooperativa. Con
  una función de potencial que liga el nivel individual con el global muestran
  que la dinámica egoísta lleva a una fase segregada de baja utilidad total, y
  la cooperativa, a una fase mezclada de alta utilidad.
- **Óptimo individual frente a óptimo social.** Con $\tau = 3/8$, una
  configuración integrada en tablero de ajedrez ($s_i = 4/8$ para todos) o en
  franjas de dos celdas de ancho ($s_i = 5/8$) satisface al 100% de los
  agentes con $\bar s$ entre 0.5 y 0.625, pero es dinámicamente inestable ante
  pequeñas fluctuaciones: la dinámica descentralizada converge hacia el estado
  segregado de baja energía interfacial ($\bar s \approx 0.8$).

---

## 4. Resumen formal

| Concepto | Resultado clave |
|---|---|
| Similitud local | $s_i = n_i^{=}/n_i$ en la vecindad de Moore ($8$ celdas) |
| Regla de Schelling | Satisfecho si $s_i \ge \tau$; si $s_i < \tau$, se muda a una celda vacía |
| Tres regímenes | $\tau \le 0.25$: mezcla; $0.3 \le \tau \le 0.7$: segregación fuerte ($\bar s \approx 0.8$); $\tau \ge 0.75$: frustración |
| Externalidad de vecindario | Mudarse reduce la similitud de los antiguos vecinos del mismo grupo y desata cascadas |
| Analogía física | Ising ferromagnético con espín conservado (Kawasaki); minimiza aproximadamente el perímetro interfacial |

---

## 5. Ejercicios formales

**Ruta básica**

- **F1.** Un agente del modelo de Schelling tiene $6$ vecinos ocupados ($2$ de
  su grupo y $4$ del otro) y $2$ celdas vecinas vacías. Si $\tau = 3/8 =
  0.375$, ¿está satisfecho? Si uno de sus vecinos del otro grupo se muda y deja
  su celda vacía, ¿cambia el estado de satisfacción de nuestro agente sin que
  él se haya movido?
- **F2.** En una cuadrícula sin celdas vacías donde los dos grupos se alternan
  como un tablero de ajedrez, calcula la similitud local $s_i$ de cada agente
  en la vecindad de Moore (8 vecinos). ¿Hasta qué umbral $\tau$ estarían todos
  satisfechos en esa configuración perfectamente integrada?
- **F3.** Explica dos diferencias fundamentales entre la dinámica del modelo de
  Schelling y la dinámica de Metropolis del modelo de Ising de la sesión 13.

**Ruta avanzada**

- **F4.** En un grafo con matriz de adyacencia $A_{ij}$ y sin celdas vacías,
  supón que dos agentes de grupos distintos intercambian lugares si y sólo si
  la suma de vecinos iguales de ambos aumenta estrictamente. Demuestra que el
  número total de aristas entre grupos distintos
  $E_{\ne} = \sum_{\langle i,j\rangle} \mathbf 1(c_i \ne c_j)$ decrece
  estrictamente en cada intercambio y concluye que el proceso debe detenerse en
  un número finito de pasos.
- **F5.** En el problema anterior, muestra con un ejemplo concreto de una
  vecindad con celdas vacías que cuando un agente insatisfecho se muda a una
  celda vacía donde sí cumple su umbral $\tau = 1/2$, el número total de
  aristas entre grupos distintos $E_{\ne}$ puede *aumentar*.

---

## Referencias

- Grauwin, S., Bertin, E., Lemoy, R., & Jensen, P. (2009). Competition between
  collective and individual dynamics. *Proceedings of the National Academy of
  Sciences*, 106(49), 20622–20626. https://doi.org/10.1073/pnas.0906263106
- Kirman, A. (2011). *Complex Economics: Individual and Collective
  Rationality*. Routledge.
- Requião da Cunha, B. (2021). *Introduction to Econophysics: Contemporary
  Approaches with Python Simulations*. CRC Press.
  https://doi.org/10.1201/9781003127956
- Schelling, T. C. (1971). Dynamic models of segregation. *Journal of
  Mathematical Sociology*, 1(2), 143–186.
  https://doi.org/10.1080/0022250X.1971.9989794
- Sen, P., & Chakrabarti, B. K. (2014). *Sociophysics: An Introduction*.
  Oxford University Press.
- Slanina, F. (2014). *Essentials of Econophysics Modelling*. Oxford
  University Press.
- Vinković, D., & Kirman, A. (2006). A physical analogue of the Schelling
  model. *Proceedings of the National Academy of Sciences*, 103(51),
  19261–19265. https://doi.org/10.1073/pnas.0609371103
- Wilensky, U., & Rand, W. (2015). *An Introduction to Agent-Based Modeling:
  Modeling Natural, Social, and Engineered Complex Systems with NetLogo*.
  MIT Press.
