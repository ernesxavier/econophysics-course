# Sesión 7: Correlaciones y matrices aleatorias. Notas formales

> **Cómo usar estas notas.** La clase combina el cuaderno
> (`Code/s07_correlaciones.jl`) y el pizarrón. El **plan de clase** de abajo
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
> Repaso recomendado: [P.B2](repaso/probabilidad.md#pb2) (covarianza y
> correlación), [L.B1–B3](repaso/algebra_lineal.md#lb1) (vectores y matrices),
> [S.B5](repaso/estadistica.md#sb5) (pruebas por permutación),
> [L.I1–I3](repaso/algebra_lineal.md#li1) (eigenvalores, teorema espectral,
> componentes principales), [L.I5](repaso/algebra_lineal.md#li5) (distancia
> entre activos); para la ruta avanzada,
> [L.A2](repaso/algebra_lineal.md#la2) (Marchenko-Pastur).
>
> Lecturas: Mantegna & Stanley (1999), caps. 12–13. Complementaria: Slanina
> (2014), sec. 1.4; Mantegna (1999). Ruta avanzada: Requião da Cunha (2021),
> sec. 5.2; Laloux et al. (1999).

**Plan de clase** (120 minutos)

| Minutos | Dónde | Qué | Secciones |
|---|---|---|---|
| 0–10 | Proyector | El fenómeno: las acciones se mueven juntas | Cuaderno, sección 0 |
| 10–22 | Pizarra | Repaso: eigenvalores como direcciones principales; la matriz de correlación y sus propiedades; cada eigenvalor es la varianza de una cartera, y $\lambda_1$ es la máxima, con la demostración | 1.1–1.2 |
| 22–35 | Pizarra | El modelo de un factor: $\rho_{ij} = \beta_i\beta_j$ y $\lambda_1 = 1 + (N-1)\rho$, con la cuenta de $C\mathbf{1}$ (el resto es el ejercicio 6 del cuaderno); $\lambda_1/N$ mide la correlación promedio | 2.1 |
| 35–45 | Pizarra | El ruido: la prueba por permutación y el ancho del espectro, $1/Q$ | 3.1–3.2 |
| 45–50 | Pizarra | La distancia $d_{ij}$ y el algoritmo de Kruskal | 4.1–4.2 |
| 50–100 | Parejas | Laboratorio: ¿señal o ruido?, el modo de mercado, los sectores y el árbol de expansión mínima; ejercicios de código 1 y 2 | Cuaderno, secciones 2 a 5 y ejercicios 1–2 |
| 100–120 | Pizarra | Cierre: puesta en común de las preguntas 3 a 5; ¿respeta el árbol los sectores? | Cuaderno, preguntas 3–5; notas, 4.3 |
| Después | Lectura | Los pesos del modo de mercado (2.2), los sectores en los eigenvectores (2.3), el resumen (5) y Para profundizar | |

**Pregunta de la sesión.** Las acciones se mueven juntas: ¿cuánto de eso es
señal y cuánto es ruido de estimación?

---

## 1. La matriz de correlación

### 1.1 Definición y propiedades

Sean $r_{ti}$ los rendimientos logarítmicos de $N$ activos durante $T$ días
($t = 1, \dots, T$, $i = 1, \dots, N$). Estandarizamos cada columna,

$$
z_{ti} = \frac{r_{ti} - \bar r_i}{s_i},
$$

con $\bar r_i$ y $s_i$ la media y la desviación estándar muestrales del activo
$i$, y formamos la matriz $Z$ de $T \times N$.

**Definición.** La matriz de correlación muestral es

$$
C = \frac{1}{T-1}\thinspace Z^\top Z,
\qquad
\rho_{ij} = C_{ij} = \frac{1}{T-1}\sum_{t=1}^{T} z_{ti}\thinspace z_{tj}.
$$

**Proposición.** $C$ es simétrica, tiene unos en la diagonal, cumple
$\lvert\rho_{ij}\rvert \le 1$ y es **semidefinida positiva**:
$w^\top C w \ge 0$ para todo $w \in \mathbb{R}^N$. Su rango es a lo más
$\min(N, T-1)$.

*Demostración.* La simetría es inmediata. Para todo $w$,

$$
w^\top C w = \frac{1}{T-1}\thinspace (Zw)^\top (Zw) = \frac{\Vert Zw\Vert^2}{T-1} \ge 0 .
$$

La diagonal vale 1 porque cada columna de $Z$ tiene varianza muestral 1, y
$\lvert\rho_{ij}\rvert \le 1$ es la desigualdad de Cauchy-Schwarz entre las
columnas $i$ y $j$ ([L.B1](repaso/algebra_lineal.md#lb1)). Las columnas de
$Z$ tienen media cero, así que viven en un subespacio de dimensión $T-1$; de
ahí la cota del rango. $\square$

**Consecuencia práctica.** Si $T - 1 < N$, $C$ tiene eigenvalores iguales a
cero: hay carteras con varianza muestral nula aunque ningún activo sea
constante. Para estudiar $N$ activos hacen falta bastantes más de $N$ días.
En el cuaderno, $N = 95$ y $T = 2952$.

### 1.2 Descomposición espectral

Por el teorema espectral ([L.I2](repaso/algebra_lineal.md#li2)), $C$ tiene
eigenvalores reales $\lambda_1 \ge \dots \ge \lambda_N \ge 0$ y eigenvectores
ortonormales $v_1, \dots, v_N$, con

$$
C = \sum_{k=1}^{N} \lambda_k\thinspace v_k v_k^\top,
\qquad
\sum_{k=1}^{N} \lambda_k = \operatorname{tr} C = N .
$$

**Interpretación como carteras.** Un vector $w$ con $\Vert w\Vert = 1$ define
una cartera de rendimientos estandarizados, $Zw$, cuya varianza muestral es
$w^\top C w$. Para $w = v_k$ esa varianza vale exactamente $\lambda_k$, y las
carteras $Zv_k$ y $Zv_l$ no están correlacionadas si $k \ne l$.

**Proposición (dirección de máxima varianza).**

$$
\boxed{\ \lambda_1 = \max_{\Vert w\Vert = 1} w^\top C w\ }
$$

y el máximo se alcanza en $w = v_1$.

*Demostración.* Escribe $w = \sum_k c_k v_k$ con $\sum_k c_k^2 = 1$. Entonces
$w^\top C w = \sum_k \lambda_k c_k^2 \le \lambda_1 \sum_k c_k^2 = \lambda_1$,
con igualdad si $w = v_1$. $\square$

El cociente $\lambda_k / N$ es la **fracción de la varianza total** que
explica la dirección $v_k$ ([L.I3](repaso/algebra_lineal.md#li3)). En los
datos del S&P 100, $\lambda_1 = 34.5$: el primer eigenvector explica el 36%
de la varianza.

---

## 2. El modo de mercado

### 2.1 El modelo de un factor

**Modelo.** Cada rendimiento estandarizado es la suma de un factor común
$f_t$ y un ruido propio:

$$
z_{ti} = \beta_i\thinspace f_t + \sqrt{1 - \beta_i^2}\thinspace\varepsilon_{ti},
$$

con $f_t$ y los $\varepsilon_{ti}$ independientes, de media 0 y varianza 1, y
$0 < \beta_i < 1$. Es el **modelo de índice único** de Sharpe (1963): el
factor común es el mercado.

**Proposición.** En el modelo de un factor, $\rho_{ij} = \beta_i \beta_j$ para
$i \ne j$. Si además todos los $\beta_i$ valen $\beta$, con $\rho = \beta^2$,

$$
C = (1-\rho)\thinspace I + \rho\thinspace\mathbf{1}\mathbf{1}^\top,
\qquad
\boxed{\ \lambda_1 = 1 + (N-1)\rho\ }
$$

cuyo mayor eigenvalor es $\lambda_1 = 1 + (N-1)\rho$, con
$v_1 = \mathbf{1}/\sqrt{N}$, y los otros $N-1$ eigenvalores valen $1 - \rho$.

*Demostración.* Por independencia,
$\mathrm{Cov}(z_i, z_j) = \beta_i\beta_j\thinspace\mathrm{Var}(f) = \beta_i\beta_j$.
Para el espectro: $C\mathbf{1} = (1-\rho)\mathbf{1} + \rho N \mathbf{1}$, y si
$u \perp \mathbf{1}$, $Cu = (1-\rho)u$
([L.I3](repaso/algebra_lineal.md#li3)). $\square$

**Comparación con los datos.** La correlación promedio entre los 4,465
pares es $\bar\rho = 0.346$, y el modelo predice
$\lambda_1 \approx 1 + 94 \times 0.346 = 33.5$. El valor observado es 34.5.
Un solo número, la correlación promedio, explica casi todo el primer
eigenvalor.

**Consecuencia.** Para $N$ grande, $\lambda_1 / N \approx \rho$: la fracción
de varianza explicada por el mercado **mide la correlación promedio**. Por eso
la gráfica de $\lambda_1/N$ en ventana móvil del cuaderno sube en las crisis:
tras la caída de marzo de 2020 llegó a 0.64 con ventanas de un año y a 0.78
con ventanas de tres meses.

### 2.2 Los pesos del modo de mercado

Con $\beta_i$ distintos, $v_1$ ya no es exactamente $\mathbf{1}/\sqrt{N}$, pero
sus componentes tienen una interpretación directa.

**Proposición.** Sea $f = Zv_1/\sqrt{\lambda_1}$ la cartera del modo de
mercado, normalizada para tener varianza muestral 1. La correlación muestral
entre el activo $i$ y esa cartera es

$$
\mathrm{corr}(z_i, f) = \sqrt{\lambda_1}\thinspace v_{1i} .
$$

*Demostración.* La covarianza muestral entre la columna $i$ de $Z$ y $Zv_1$
es la componente $i$ de $C v_1 = \lambda_1 v_1$, es decir, $\lambda_1 v_{1i}$.
La varianza de $Zv_1$ es $v_1^\top C v_1 = \lambda_1$. Dividiendo entre
$\sqrt{\lambda_1}$, que es la desviación estándar de $Zv_1$, queda
$\sqrt{\lambda_1}\thinspace v_{1i}$. $\square$

Es decir, el peso de cada acción en $v_1$ es proporcional a su $\beta$
respecto al mercado. En los datos, esas correlaciones van de 0.41 (Eli Lilly,
Netflix) a 0.81 (BlackRock, Berkshire Hathaway). Por sector, las finanzas son
las más ligadas al mercado (0.76 en promedio) y los servicios públicos los
menos (0.48).

### 2.3 Sectores

Si además del mercado hay factores de sector, aparecen más eigenvalores
grandes. Como los eigenvectores son ortogonales a $v_1$, que tiene todas sus
componentes positivas, $v_2, v_3, \dots$ deben tener componentes de **los
dos signos**: son **contrastes** entre grupos de acciones (unas suben cuando
otras bajan, una vez descontado el mercado). Por eso el cuaderno encuentra
$v_2$ = defensivas contra tecnología y $v_4$ = semiconductores contra
software, en lugar de un eigenvector por sector.

---

## 3. Ruido de estimación

### 3.1 La prueba por permutación

Aunque los activos fueran independientes, la $C$ muestral no sería la
identidad. Para cada par independiente, $\rho_{ij}$ tiene media
aproximadamente 0 y desviación estándar aproximadamente $1/\sqrt{T}$. Con
$T = 2952$ eso es 0.018, y es lo que se obtiene con los datos barajados.

**Hipótesis nula.** Los $N$ activos son independientes entre sí (cada uno
con su propia distribución, sin importar cuál sea).

**Estadístico.** El espectro de $C$; en particular, su mayor eigenvalor.

**Distribución nula por permutación.** Barajar los $T$ rendimientos de cada
activo con una permutación **distinta** conserva exactamente la distribución
de cada activo, incluidas sus colas pesadas, y destruye la alineación de los
días entre activos. Bajo la hipótesis nula, todas esas permutaciones son
igual de probables ([S.B5](repaso/estadistica.md#sb5)).

**Resultado.** Con los datos barajados, los eigenvalores quedan
aproximadamente entre 0.7 y 1.35 (los valores exactos cambian un poco con
cada permutación). Con los datos reales, 8 eigenvalores quedan por encima de
ese máximo y explican juntos el 58% de la varianza.

Dos observaciones:

1. Barajar con la **misma** permutación para todos los activos no cambia $C$:
   la correlación no depende del orden de los días. Por eso las
   permutaciones deben ser independientes entre activos.
2. La prueba destruye también la dependencia temporal de cada serie (la
   volatilidad agrupada de la sesión 2). Para la correlación contemporánea
   eso casi no importa, pero es otra diferencia entre datos y nula.

**Muchas pruebas a la vez.** Probar si cada uno de los 4,465 pares tiene
$\rho_{ij} \ne 0$ es un problema de pruebas múltiples
([S.I2](repaso/estadistica.md#si2)). Mirar el espectro completo evita ese
problema: es una sola prueba sobre la estructura de toda la matriz.

### 3.2 El ancho del espectro de ruido

Un cálculo elemental explica por qué los eigenvalores del ruido se dispersan
alrededor de 1.

**Proposición.** Para cualquier matriz de correlación,

$$
\frac{1}{N}\sum_{k=1}^{N} \lambda_k = 1,
\qquad
\frac{1}{N}\sum_{k=1}^{N} \lambda_k^2 = \frac{1}{N}\operatorname{tr} C^2 = 1 + \frac{1}{N}\sum_{i \ne j} \rho_{ij}^2 .
$$

Si los activos son independientes, $\mathbb{E}[\rho_{ij}^2] \approx 1/T$, así que la
varianza de los eigenvalores es aproximadamente

$$
\boxed{\ \frac{1}{N}\sum_{k} \lambda_k^2 - 1 \approx \frac{N-1}{T} \approx \frac{1}{Q},
\qquad Q = \frac{T}{N}\ }
$$

*Demostración.* $\operatorname{tr} C^2 = \sum_k \lambda_k^2$ porque los
eigenvalores de $C^2$ son $\lambda_k^2$, y también
$\operatorname{tr} C^2 = \sum_{i,j} \rho_{ij}^2 = N + \sum_{i \ne j} \rho_{ij}^2$. $\square$

Aunque cada correlación sea pequeña, hay $N(N-1)$ de ellas: el ruido se
acumula y abre el espectro con desviación estándar $1/\sqrt{Q}$. En el
cuaderno, $Q = 31.1$ y la varianza de los eigenvalores barajados es 0.032,
igual a $1/Q$.

### 3.3 [Avanzado] La ley de Marchenko-Pastur

El desarrollo de esta subsección está en [Para profundizar, P.1](#pp1).


---

## 4. Del mercado a una red

### 4.1 Distancia entre activos

**Definición.** $d_{ij} = \sqrt{2(1 - \rho_{ij})}$.

**Proposición.** $d$ es una distancia: $d_{ij} \ge 0$, con $d_{ij} = 0$ sólo
si $\rho_{ij} = 1$; $d_{ij} = d_{ji}$; y cumple la desigualdad del triángulo.

*Demostración.* Sea $\hat z_i = z_i / \sqrt{T-1}$, de modo que
$\Vert\hat z_i\Vert = 1$ y $\hat z_i \cdot \hat z_j = \rho_{ij}$. Entonces
$\Vert\hat z_i - \hat z_j\Vert^2 = 2 - 2\rho_{ij} = d_{ij}^2$: $d_{ij}$ es la
distancia euclidiana entre los vectores $\hat z_i$ y $\hat z_j$, que hereda
todas las propiedades ([L.I5](repaso/algebra_lineal.md#li5)). $\square$

Va de 0 ($\rho = 1$) a 2 ($\rho = -1$); dos activos no correlacionados están a
distancia $\sqrt{2}$.

### 4.2 El árbol de expansión mínima

Consideramos la red completa con $N$ nodos (los activos) y peso $d_{ij}$ en
cada enlace.

**Definición.** Un **árbol de expansión** es un subconjunto de $N-1$
enlaces que conecta a los $N$ nodos sin formar ciclos. El **árbol de
expansión mínima** (MST, por sus siglas en inglés) es el de menor suma de
pesos.

**Algoritmo de Kruskal (1956).**

1. Ordena los $N(N-1)/2$ enlaces de menor a mayor distancia.
2. Recórrelos en ese orden. Agrega cada enlace si une dos componentes
   distintas; descártalo si cerraría un ciclo.
3. Termina al tener $N-1$ enlaces.

**[Avanzado]** Proposición: en [Para profundizar, P.2](#pp2).

**[Avanzado]** Proposición (ultramétrica subdominante): en [Para profundizar, P.3](#pp3).

### 4.3 ¿Respeta el árbol los sectores?

Sea $X$ el número de enlaces del MST que unen dos acciones del mismo sector.

**Valor esperado sin relación con el sector.** Si las etiquetas de sector se
reparten al azar entre las acciones, cada enlace une dos acciones del mismo
sector con probabilidad

$$
p = \sum_{s} \frac{n_s(n_s - 1)}{N(N-1)},
$$

donde $n_s$ es el número de acciones del sector $s$. Con los 11 sectores del
S&P 100, $p = 0.129$ y $\mathbb{E}[X] = 94p = 12.1$, lo mismo que dan las
999 permutaciones del cuaderno.

**Observado.** $X = 65$. Ninguna permutación se acerca, así que el p-valor
por permutación es $1/1000$, el mínimo posible con 999 permutaciones
([S.B5](repaso/estadistica.md#sb5)). El árbol se construyó sólo con
correlaciones, sin usar los sectores: la clasificación sectorial **emerge**
de los datos.

---

## 5. Resumen formal

| Objeto | Definición | Resultado clave |
|---|---|---|
| Matriz de correlación | $C = Z^\top Z/(T-1)$ | Simétrica, semidefinida positiva, traza $N$ |
| Espectro | $C v_k = \lambda_k v_k$ | $\lambda_k$ = varianza de la cartera $v_k$; $\lambda_1$ es la máxima |
| Modelo de un factor | $z_{ti} = \beta_i f_t + \sqrt{1-\beta_i^2}\thinspace\varepsilon_{ti}$ | $\lambda_1 \approx 1 + (N-1)\bar\rho$; $v_{1i} \propto \beta_i$ |
| Ruido | Activos independientes | Varianza del espectro $\approx 1/Q$ |
| Marchenko-Pastur | $N, T \to \infty$, $Q = T/N$ | Eigenvalores en $[\lambda_-, \lambda_+]$, $\lambda_\pm = (1 \pm 1/\sqrt{Q})^2$ |
| Distancia | $d_{ij} = \sqrt{2(1-\rho_{ij})}$ | Euclidiana; de 0 a 2 |
| MST | $N-1$ enlaces de suma mínima | Kruskal; ultramétrica subdominante |

---

## 6. Ejercicios formales

**Ruta básica**

- **F1.** Demuestra que si $C$ es una matriz de correlación de $N$ activos,
  su eigenvalor más grande cumple $\lambda_1 \ge 1 + (N-1)\bar\rho$, donde
  $\bar\rho$ es la correlación promedio entre pares. (Sugerencia: evalúa
  $w^\top C w$ en $w = \mathbf{1}/\sqrt{N}$.)
- **F2.** En el modelo de un factor con $\beta$ igual para todos, ¿qué
  fracción de la varianza total explica el mercado? ¿Qué le pasa a la
  varianza de una cartera con el mismo peso en las $N$ acciones cuando
  $N \to \infty$? Interpreta: ¿qué riesgo no se puede diversificar?
- **F3.** Calcula $\lambda_\pm$ para $N = 95$ y $T = 250$ (un año de datos).
  ¿Cuántos de los eigenvalores del cuaderno (34.5, 6.45, 5.01, 2.38, 2.22,
  1.68, 1.48, 1.38) superarían ese umbral?
- **F4.** Dos acciones tienen $\rho = 0.88$ y otras dos $\rho = 0.10$. Calcula
  sus distancias. ¿Qué correlación corresponde a $d = 1$?

**Ruta avanzada**

- **F5.** Para el modelo de un factor con $\beta_i$ distintos, demuestra que
  $C = D + \beta\beta^\top$ con $D$ diagonal, y que si $D \approx (1-\bar\beta^2)I$
  el primer eigenvector es aproximadamente $\beta/\Vert\beta\Vert$.
- **F6.** Verifica que la densidad de Marchenko-Pastur integra 1 y tiene
  media 1 (puedes hacerlo numéricamente para $Q = 2$ y $Q = 10$).
- **F7.** Demuestra que todo árbol con $N$ nodos tiene exactamente $N-1$
  enlaces.
- **F8.** Demuestra que la ultramétrica subdominante cumple
  $u_{ij} \le d_{ij}$, y que es la mayor ultramétrica con esa propiedad.

---

## Para profundizar

Material opcional de la ruta avanzada: demostraciones y temas que no se ven
en clase. Cada apartado indica de qué sección de las notas viene.

<a id="pp1"></a>

### P.1 La ley de Marchenko-Pastur (sección 3.3)


**Teorema (Marchenko y Pastur, 1967).** Sean $z_{ti}$ variables i.i.d. de
media 0 y varianza 1, y $C = Z^\top Z / T$. Si $N, T \to \infty$ con
$Q = T/N \ge 1$ fijo, la distribución de los eigenvalores de $C$ converge a la
densidad

$$
\rho(\lambda) = \frac{Q}{2\pi}\thinspace\frac{\sqrt{(\lambda_+ - \lambda)(\lambda - \lambda_-)}}{\lambda},
\qquad
\lambda_\pm = \left(1 \pm \frac{1}{\sqrt{Q}}\right)^2,
$$

para $\lambda_- \le \lambda \le \lambda_+$, y a cero fuera de ese intervalo.

No damos la demostración (ver Requião da Cunha, sec. 5.2). Es consistente con
la sección 3.2: esta densidad tiene media 1 y varianza $1/Q$.

**Uso (Laloux et al., 1999; Plerou et al., 1999).** Los eigenvalores
empíricos dentro de $[\lambda_-, \lambda_+]$ son indistinguibles del ruido.
Para el S&P 100, $\lambda_+ = 1.39$ y 7 eigenvalores quedan por encima; el
octavo, 1.38, queda en el borde. Los artículos originales, con cientos de
acciones, encontraron lo mismo: sólo unos cuantos eigenvalores (alrededor del
6%) llevan información, y el resto de la matriz es ruido de estimación.

**Tres cuidados.**

1. El teorema supone varianza finita e independencia en el tiempo. Los
   rendimientos tienen colas pesadas y volatilidad agrupada; la prueba por
   permutación (sección 3.1) conserva las colas y por eso es la comparación
   más directa.
2. Los eigenvalores grandes se llevan parte de la traza $N$, así que el
   grueso del espectro real queda a la izquierda del de ruido. Una corrección
   habitual compara con Marchenko-Pastur con varianza $1 - \lambda_1/N$ en
   lugar de 1.
3. Es un resultado asintótico: con $N = 95$ los bordes no son exactos, y unos
   pocos eigenvalores de ruido pueden asomarse fuera del intervalo.

**Por qué importa.** Una cartera de mínima varianza usa $C^{-1}$, y los
eigenvalores pequeños del ruido se vuelven enormes al invertirlos: la cartera
apuesta en direcciones que son puro ruido. "Limpiar" $C$ reemplazando los
eigenvalores dentro del intervalo de ruido por su promedio da carteras con
riesgo realizado menor (Laloux et al., 1999).

<a id="pp2"></a>

### P.2 Proposición (sección 4.2)

**Proposición.** El algoritmo de Kruskal produce un árbol de
expansión mínima. Si todas las distancias son distintas, ese árbol es único.

*Demostración (propiedad del corte).* Supongamos que las distancias son
distintas; con empates el argumento se adapta. Sea $S$ un conjunto de nodos y $e$ el
enlace más corto entre $S$ y su complemento. Todo MST contiene a $e$: si un
árbol $A$ no lo contuviera, agregarle $e$ formaría un ciclo, ese ciclo cruza
de $S$ a su complemento por otro enlace $e'$ con $d(e') > d(e)$, y cambiar
$e'$ por $e$ daría un árbol más corto. Cada enlace que Kruskal agrega es el
más corto entre la componente de uno de sus extremos y el resto: cualquier
enlace más corto que saliera de esa componente ya se habría examinado y, al
no cerrar ciclo, agregado. Así, todos los enlaces que agrega pertenecen a
todo MST; como son $N-1$, el MST es único y es el que construye Kruskal. $\square$

<a id="pp3"></a>

### P.3 Proposición (ultramétrica subdominante) (sección 4.2)

**Proposición (ultramétrica subdominante).** Para dos activos
$i$ y $j$, sea $u_{ij}$ la mayor distancia entre los enlaces del camino que
los une en el MST. Entonces

$$
u_{ij} \le \max\left(u_{ik},\ u_{kj}\right) \quad \text{para todo } k,
$$

una desigualdad más fuerte que la del triángulo (es una **ultramétrica**).
Mantegna (1999) la usó para ordenar las acciones en una **jerarquía**: es el
mismo árbol que produce el agrupamiento jerárquico de enlace simple.

*Idea de la demostración.* El camino de $i$ a $j$ en el árbol está contenido
en la unión de los caminos de $i$ a $k$ y de $k$ a $j$, así que su enlace más
largo no puede superar al más largo de esos dos caminos. $\square$

---

## Referencias

- Kruskal, J. B. (1956). On the shortest spanning subtree of a graph and the
  traveling salesman problem. *Proceedings of the American Mathematical
  Society*, 7(1), 48–50. https://doi.org/10.1090/S0002-9939-1956-0078686-7
- Laloux, L., Cizeau, P., Bouchaud, J.-P., & Potters, M. (1999). Noise
  dressing of financial correlation matrices. *Physical Review Letters*,
  83(7), 1467–1470. https://doi.org/10.1103/PhysRevLett.83.1467
- Mantegna, R. N. (1999). Hierarchical structure in financial markets. *The
  European Physical Journal B*, 11(1), 193–197.
  https://doi.org/10.1007/s100510050929
- Mantegna, R. N., & Stanley, H. E. (1999). *An Introduction to Econophysics:
  Correlations and Complexity in Finance*. Cambridge University Press.
  https://doi.org/10.1017/CBO9780511755767
- Marčenko, V. A., & Pastur, L. A. (1967). Distribution of eigenvalues for
  some sets of random matrices. *Mathematics of the USSR-Sbornik*, 1(4),
  457–483. https://doi.org/10.1070/SM1967v001n04ABEH001994
- Plerou, V., Gopikrishnan, P., Rosenow, B., Amaral, L. A. N., & Stanley,
  H. E. (1999). Universal and nonuniversal properties of cross correlations
  in financial time series. *Physical Review Letters*, 83(7), 1471–1474.
  https://doi.org/10.1103/PhysRevLett.83.1471
- Requião da Cunha, C. (2021). *Introduction to Econophysics: Contemporary
  Approaches with Python Simulations*. CRC Press.
  https://doi.org/10.1201/9781003127956
- Sharpe, W. F. (1963). A simplified model for portfolio analysis.
  *Management Science*, 9(2), 277–293. https://doi.org/10.1287/mnsc.9.2.277
- Slanina, F. (2014). *Essentials of Econophysics Modelling*. Oxford
  University Press.
