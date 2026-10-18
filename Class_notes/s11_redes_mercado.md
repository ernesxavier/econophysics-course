# Sesión 11: Redes: el mercado y la robustez. Notas formales

> **Cómo leer estas notas.** Siguen el orden del programa: primero se observa
> en el cuaderno (`Code/s11_redes_mercado.jl`), después se formaliza aquí. El
> cuerpo de las notas es la **ruta básica**. Las demostraciones y los temas de
> la ruta avanzada, opcionales, están al final, en **Para profundizar**; donde
> estaban queda una línea que los enlaza. Las definiciones de matriz de
> adyacencia, agrupamiento y de los tres modelos de redes aleatorias
> (Erdős-Rényi, Watts-Strogatz y Barabási-Albert) están en las [notas de la
> sesión 10](s10_redes_modelos.md).
>
> Repaso recomendado: [F.B4](repaso/fisica_estadistica.md#fb4) (transiciones
> de fase), [L.I2](repaso/algebra_lineal.md#li2) (matrices simétricas y
> teorema espectral), [L.I5](repaso/algebra_lineal.md#li5) (distancia entre
> activos); para la ruta avanzada, [L.A1](repaso/algebra_lineal.md#la1) (el
> Laplaciano de una red) y [F.A5](repaso/fisica_estadistica.md#fa5)
> (percolación).
>
> Lecturas: Slanina (2014), cap. 6. **[Avanzado]** Newman (2018), caps. 12–13
> y 15; Jackson (2008), caps. 2–5 (consulta).

---

## 1. Una red de correlaciones del mercado

A partir de la matriz de correlación $C$ de los rendimientos
([sesión 7](s07_correlaciones.md)) se define la **red de umbral**

$$
A_{ij} = \mathbb 1\left[\rho_{ij} > u\right], \qquad i \ne j .
$$

A diferencia del árbol de expansión mínima, que siempre tiene $n - 1$
enlaces y ningún ciclo, la red de umbral conserva los triángulos. Con las
95 acciones del S&P 100 (2015–2026) y $u = 0.5$:

| Red | Enlaces | $\langle k\rangle$ | $k_{\max}$ | Agrupamiento | Distancia promedio | Componente gigante |
|---|---|---|---|---|---|---|
| Mercado | 455 | 9.6 | 37 | 0.59 | 2.84 | 86% |
| Erdős-Rényi | 455 | 9.6 | 19 | 0.09 | 2.24 | 100% |
| Barabási-Albert | 450 | 9.5 | 30 | 0.16 | 2.23 | 100% |

**Por qué el agrupamiento es alto.** La correlación tiende a ser
"transitiva": si los vectores de rendimientos estandarizados $\hat z_i$,
$\hat z_j$ y $\hat z_l$ tienen norma 1, entonces
$\Vert \hat z_j - \hat z_l\Vert \le \Vert \hat z_j - \hat z_i\Vert + \Vert \hat z_i - \hat z_l\Vert$
([L.I5](repaso/algebra_lineal.md#li5)), y con $d = \sqrt{2(1-\rho)}$ eso
acota $\rho_{jl}$ por abajo cuando $\rho_{ij}$ y $\rho_{il}$ son altas. Por
ejemplo, si $\rho_{ij} = \rho_{il} = 0.75$, entonces $d_{jl} \le 2\sqrt{0.5}$
y $\rho_{jl} \ge 0$. La geometría no garantiza triángulos, pero los
favorece; los sectores hacen el resto.

Las redes de umbral fueron estudiadas por Boginski, Butenko y Pardalos
(2005) y por Onnela, Kaski y Kertész (2004), que compararon su agrupamiento
con el del árbol de expansión mínima.

---

## 2. Robustez y percolación

### 2.1 El problema

Se quita una fracción $f$ de los nodos, con todos sus enlaces, y se mide la
fracción $S(f)$ de los $n$ nodos originales que sigue en el componente
gigante ([F.B4](repaso/fisica_estadistica.md#fb4)). En **fallas aleatorias**
los nodos se eligen al azar. En un **ataque dirigido** se quitan en orden
decreciente de grado.

En el cuaderno, con $n = 10^4$ y $\langle k\rangle = 4$:

| Red | Fallas: $S$ baja de 1% en | Ataque: se rompe en |
|---|---|---|
| Erdős-Rényi | $f \approx 0.77$ | $f \approx 0.37$ |
| Barabási-Albert | $f \approx 0.83$ | $f \approx 0.15$ |

Albert, Jeong y Barabási (2000) observaron este contraste en internet y en
la web: robustas ante fallas, frágiles ante ataques.

### 2.2 [Avanzado] El criterio de Molloy y Reed

El desarrollo de esta subsección está en [Para profundizar, P.1](#pp1).


### 2.3 [Avanzado] Umbral de fallas aleatorias

El desarrollo de esta subsección está en [Para profundizar, P.2](#pp2).


---

## 3. [Avanzado] El Laplaciano

El desarrollo de esta sección está en [Para profundizar, P.3](#pp3).


---

## 4. [Avanzado] Modelo de configuración y asortatividad

El desarrollo de esta sección está en [Para profundizar, P.4](#pp4).


---

## 5. Resumen formal

| Concepto | Resultado |
|---|---|
| Red de umbral del mercado | $A_{ij} = \mathbb{1}[\rho_{ij} > u]$; agrupamiento alto por estructura sectorial |
| Robustez empírica | Libres de escala: resistentes a fallas al azar, frágiles ante ataques dirigidos |
| Molloy-Reed | Gigante si $\kappa = \langle k^2\rangle/\langle k\rangle > 2$ |
| Fallas aleatorias | $f_c = 1 - 1/(\kappa - 1)$; libre de escala: $f_c \to 1$ |
| Laplaciano | Ceros = componentes; vector de Fiedler = partición |
| Modelo de configuración | Conserva grados al azar; aísla cuánto agrupamiento viene sólo de la secuencia $k_i$ |

---

## 6. Ejercicios formales

**Ruta básica**

- **F1.** En la red del mercado, ¿por qué usamos la red de umbral
  $\rho_{ij} > u$ en lugar del árbol de expansión mínima de la sesión 7 para
  comparar el agrupamiento con los modelos aleatorios?
- **F2.** Las cadenas de suministro globales tienen algunos proveedores de los
  que dependen muchas empresas (por ejemplo, fabricantes de semiconductores).
  Según la sección 2, ¿ante qué tipo de choque son robustas y ante cuál son
  vulnerables?

**Ruta avanzada**

- **F3.** Usa el criterio de Molloy y Reed para una red con
  $p_k \propto k^{-\gamma}$ y $k \ge 1$. ¿Para qué valores de $\gamma$ hay
  componente gigante sin importar cuántos nodos se quiten al azar (en el
  límite $n \to \infty$)?
- **F4.** Demuestra que $\lambda_2 > 0$ si y sólo si la red es conexa, y
  calcula el espectro del Laplaciano de un anillo de $n$ nodos.

---

## Para profundizar

Material opcional de la ruta avanzada: demostraciones y temas que no se ven
en clase. Cada apartado indica de qué sección de las notas viene.

<a id="pp1"></a>

### P.1 El criterio de Molloy y Reed (sección 2.2)


Considera una red aleatoria con distribución de grados $p_k$ y por lo demás
al azar (el modelo de configuración de la sección 4).

**Proposición.** Hay componente gigante si y sólo si

$$
\kappa \equiv \frac{\langle k^2\rangle}{\langle k\rangle} > 2 .
$$

*Argumento.* Exploramos la red desde un nodo, siguiendo enlaces. Al llegar a
un vecino por un enlace, la probabilidad de llegar a un nodo de grado $k$ es
proporcional a $k p_k$, porque un nodo con $k$ enlaces tiene $k$ maneras de
ser alcanzado. Desde ahí puede seguir por $k - 1$ enlaces nuevos (el
**grado excedente**). El número promedio de enlaces nuevos es

$$
\sum_k (k - 1)\frac{k p_k}{\langle k\rangle}
= \frac{\langle k^2\rangle - \langle k\rangle}{\langle k\rangle}
= \kappa - 1 .
$$

La exploración es un proceso de ramificación: crece sin límite (componente
gigante) si cada rama genera en promedio más de una rama nueva, es decir, si
$\kappa - 1 > 1$. $\square$

Para Erdős-Rényi, $\langle k^2\rangle = c^2 + c$, así que $\kappa = c + 1$ y
el criterio da $c > 1$, como vimos en la [sesión 10](s10_redes_modelos.md).

<a id="pp2"></a>

### P.2 Umbral de fallas aleatorias (sección 2.3)


**Proposición** (Cohen et al., 2000). Si se quita al azar una fracción $f$
de los nodos de una red con parámetro $\kappa$, el componente gigante
desaparece en ([F.A5](repaso/fisica_estadistica.md#fa5))

$$
f_c = 1 - \frac{1}{\kappa - 1} .
$$

*Demostración.* Sea $q = 1 - f$. Cada vecino de un nodo que sobrevive
sobrevive con probabilidad $q$, así que el nuevo grado $K'$ de ese nodo es
binomial con parámetros $K$ y $q$:
$\langle k'\rangle = q\langle k\rangle$ y
$\langle k'^2\rangle = q^2\langle k^2\rangle + q(1-q)\langle k\rangle$.
Entonces

$$
\kappa' = \frac{\langle k'^2\rangle}{\langle k'\rangle} = q\kappa + 1 - q .
$$

La condición $\kappa' = 2$ da $q(\kappa - 1) = 1$. $\square$

**Consecuencias.**

- Erdős-Rényi: $f_c = 1 - 1/c$. Con $c = 4$, $f_c = 0.75$; la simulación da
  0.77.
- Barabási-Albert: $\kappa \approx \frac{m+1}{2}\ln n$ crece sin límite, así
  que $f_c \to 1$: en una red infinita **no hay umbral**. Con $n = 10^3$,
  $10^4$ y $10^5$ (y $m = 2$), $\kappa$ vale 11.1, 14.9 y 19.0, y $f_c$ vale
  0.90, 0.93 y 0.94. Cerca de $f_c$ el componente gigante es diminuto, por
  eso la simulación lo "pierde" antes.
- Un **ataque** quita los nodos de mayor grado: corta la cola de $p_k$, y
  $\langle k^2\rangle$ cae muy rápido. Para Barabási-Albert basta quitar
  una fracción pequeña de los nodos (Callaway et al., 2000; Cohen et al.,
  2001).

<a id="pp3"></a>

### P.3 El Laplaciano (sección 3)


El Laplaciano es $L = D - A$, con $D = \mathrm{diag}(k_1, \dots, k_n)$
([L.I2](repaso/algebra_lineal.md#li2), [L.A1](repaso/algebra_lineal.md#la1)).

**Proposición.**

1. $x^\top L x = \sum_{\text{enlaces } (i,j)} (x_i - x_j)^2 \ge 0$, así que
   los eigenvalores son reales y no negativos:
   $0 = \lambda_1 \le \lambda_2 \le \dots \le \lambda_n$.
2. La multiplicidad del eigenvalor 0 es el número de componentes conexas.

*Demostración.* (1) $x^\top L x = \sum_i k_i x_i^2 - \sum_{i,j} A_{ij}x_ix_j$.
Cada enlace $(i, j)$ aporta $x_i^2 + x_j^2$ al primer término y $2x_ix_j$ al
segundo, es decir, $(x_i - x_j)^2$. (2) $Lx = 0$ implica $x^\top L x = 0$,
así que $x_i = x_j$ en cada enlace: $x$ es constante en cada componente.
Los vectores indicadores de las componentes son una base de ese espacio. $\square$

**Conectividad algebraica.** Por (1), $\lambda_2 = \min x^\top L x$ sobre los
$x$ de norma 1 ortogonales a $\mathbf 1$. Si la red se puede partir en dos
grupos con pocos enlaces entre ellos, el vector que vale $+1$ en un grupo y
$-1$ en el otro (centrado y normalizado) da un valor pequeño. Por eso el
eigenvector de $\lambda_2$ (**vector de Fiedler**) sugiere una partición:
los nodos con componente positiva de un lado, los negativos del otro.

En el cuaderno, la red de umbral del mercado con $u = 0.5$ tiene 11
componentes (10 son acciones aisladas o pares) y su Laplaciano tiene 11
ceros. En el componente gigante, el vector de Fiedler separa a las
acciones defensivas (consumo básico, servicios públicos, buena parte de
salud) de las cíclicas, igual que el segundo eigenvector de la matriz de
correlación en la [sesión 7](s07_correlaciones.md).

<a id="pp4"></a>

### P.4 Modelo de configuración y asortatividad (sección 4)


#### 4.1 El modelo de configuración

**Modelo.** Se fija una sucesión de grados $k_1, \dots, k_n$. Cada nodo
recibe $k_i$ "medios enlaces" y se emparejan todos al azar. Es la red más
aleatoria posible con esos grados; Erdős-Rényi es el caso de grados Poisson.

La probabilidad de un enlace entre $i$ y $j$ es $p_{ij} \approx k_ik_j/(2m)$,
y el agrupamiento es

$$
C = \frac{1}{n}\frac{\left(\langle k^2\rangle - \langle k\rangle\right)^2}{\langle k\rangle^3} .
$$

Tiende a cero si $\langle k^2\rangle$ es finito, pero puede ser apreciable
con grados muy desiguales y $n$ pequeño. Con los grados del mercado
($u = 0.5$) la fórmula da 0.35 y la red recableada del cuaderno tiene 0.31,
contra 0.59 de la red real. El resto del agrupamiento es estructura de
sectores, no consecuencia de los grados.

En la práctica, la red se muestrea recableando la real: se eligen dos
enlaces $(a, b)$ y $(c, d)$ y se cambian por $(a, d)$ y $(c, b)$, si eso no
crea enlaces repetidos ni lazos. Cada cambio conserva todos los grados.

#### 4.2 Asortatividad

La **asortatividad** $r$ es el coeficiente de correlación de Pearson entre
los grados de los dos extremos de un enlace, calculado sobre todos los
enlaces en los dos sentidos (Newman, 2002). Va de $-1$ a 1.

- $r > 0$: los nodos de grado alto se conectan entre sí. Típico de redes
  sociales. El mercado tiene $r = 0.12$: un núcleo financiero denso.
- $r < 0$: los concentradores se conectan con nodos de grado bajo. Típico
  de redes tecnológicas y biológicas, y de redes simples con grados muy
  desiguales: si sólo puede haber un enlace por par, los concentradores no
  encuentran suficientes concentradores. Por eso la red recableada tiene
  $r = -0.23$.

---

## Referencias

- Albert, R., Jeong, H., & Barabási, A.-L. (2000). Error and attack
  tolerance of complex networks. *Nature*, 406(6794), 378–382.
  https://doi.org/10.1038/35019019
- Boginski, V., Butenko, S., & Pardalos, P. M. (2005). Statistical analysis
  of financial networks. *Computational Statistics & Data Analysis*, 48(2),
  431–443. https://doi.org/10.1016/j.csda.2004.02.004
- Callaway, D. S., Newman, M. E. J., Strogatz, S. H., & Watts, D. J. (2000).
  Network robustness and fragility: Percolation on random graphs. *Physical
  Review Letters*, 85(25), 5468–5471.
  https://doi.org/10.1103/PhysRevLett.85.5468
- Cohen, R., Erez, K., ben-Avraham, D., & Havlin, S. (2000). Resilience of
  the Internet to random breakdowns. *Physical Review Letters*, 85(21),
  4626–4628. https://doi.org/10.1103/PhysRevLett.85.4626
- Cohen, R., Erez, K., ben-Avraham, D., & Havlin, S. (2001). Breakdown of
  the Internet under intentional attack. *Physical Review Letters*, 86(16),
  3682–3685. https://doi.org/10.1103/PhysRevLett.86.3682
- Jackson, M. O. (2008). *Social and Economic Networks*. Princeton
  University Press.
- Molloy, M., & Reed, B. (1995). A critical point for random graphs with a
  given degree sequence. *Random Structures & Algorithms*, 6(2–3), 161–180.
  https://doi.org/10.1002/rsa.3240060204
- Newman, M. E. J. (2002). Assortative mixing in networks. *Physical Review
  Letters*, 89(20), 208701. https://doi.org/10.1103/PhysRevLett.89.208701
- Newman, M. (2018). *Networks* (2.ª ed.). Oxford University Press.
  https://doi.org/10.1093/oso/9780198805090.001.0001
- Onnela, J.-P., Kaski, K., & Kertész, J. (2004). Clustering and
  information in correlation based financial networks. *The European
  Physical Journal B*, 38(2), 353–362.
  https://doi.org/10.1140/epjb/e2004-00128-7
- Slanina, F. (2014). *Essentials of Econophysics Modelling*. Oxford
  University Press.
