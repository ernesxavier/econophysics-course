# Sesión 12: Contagio y riesgo sistémico. Notas formales

> **Cómo leer estas notas.** Siguen el orden del programa: primero se observa
> en el cuaderno (`Code/s12_contagio.jl`), después se formaliza aquí. El
> cuerpo de las notas es la **ruta básica**. Las demostraciones y los temas de
> la ruta avanzada, opcionales, están al final, en **Para profundizar**; donde
> estaban queda una línea que los enlaza.
>
> Repaso recomendado: [L.B4](repaso/algebra_lineal.md#lb4) (matriz de
> adyacencia), [C.B5](repaso/calculo.md#cb5) (serie geométrica),
> [L.I4](repaso/algebra_lineal.md#li4) (potencias de matrices),
> [C.I5](repaso/calculo.md#ci5) (puntos fijos); para la ruta avanzada,
> [L.A4](repaso/algebra_lineal.md#la4) (series de matrices y puntos fijos),
> [L.A5](repaso/algebra_lineal.md#la5) (Perron-Frobenius) y
> [F.A5](repaso/fisica_estadistica.md#fa5) (percolación).
>
> Lecturas: Battiston et al. (2016). Complementaria: Watts (2002). Ruta
> avanzada: Eisenberg & Noe (2001); Battiston et al. (2012); de consulta,
> Acemoglu, Ozdaglar & Tahbaz-Salehi (2015) y Glasserman & Young (2016).

---

## 1. Cascadas de umbral

### 1.1 El modelo de Watts

**Modelo** (Watts, 2002). En una red no dirigida, cada nodo $i$ está en uno
de dos estados, 0 (en pie) o 1 (caído), y tiene un umbral $\varphi_i$. Un
nodo en pie cae si la fracción de sus vecinos caídos es al menos su umbral:

$$
\frac{1}{k_i}\sum_j A_{ij}\thinspace s_j \ge \varphi_i .
$$

Los nodos caídos no se levantan. Al principio cae un nodo al azar y se
aplica la regla hasta que nada cambia. Como el conjunto de nodos caídos sólo
crece y la red es finita, el proceso siempre termina.

En el cuaderno, todos los nodos tienen el mismo umbral $\varphi$.

**Nodos vulnerables.** Un nodo de grado $k$ cae con un solo vecino caído si
$1/k \ge \varphi$, es decir, si $k \le 1/\varphi$. Los demás necesitan al
menos dos vecinos caídos.

**Cascada global.** En una red de $n$ nodos, una cascada es global si
alcanza una fracción finita de la red cuando $n \to \infty$. En una
simulación se usa un corte: en el cuaderno, más del 5% de los nodos.

### 1.2 La condición de cascada

**Proposición** (Watts, 2002). En una red aleatoria con distribución de
grados $p_k$, grado promedio $z$ y umbral común $\varphi$, una cascada global
es posible si

$$
\sum_{k=1}^{\lfloor 1/\varphi\rfloor} k(k-1)\thinspace p_k > z .
$$

*Argumento.* En una red grande y dispersa, al principio la cascada sólo
avanza a través de nodos vulnerables, porque es muy poco probable que un
nodo no vulnerable tenga dos vecinos caídos antes de que la cascada sea
grande. Basta, entonces, con que los nodos vulnerables formen un grupo
conectado gigante: es percolación restringida a ellos
([F.A5](repaso/fisica_estadistica.md#fa5)).

Como en la sesión 11 (criterio de Molloy y Reed), al seguir un enlace se
llega a un nodo de grado $k$ con probabilidad $kp_k/z$. Si ese nodo es
vulnerable, abre $k - 1$ caminos nuevos. El número promedio de caminos
nuevos hacia nodos vulnerables es

$$
\sum_{k \le 1/\varphi} (k-1)\thinspace\frac{k\thinspace p_k}{z} .
$$

El grupo crece sin límite si este número es mayor que 1. $\square$

Si la cascada arranca, puede alcanzar también a nodos no vulnerables que
tengan varios vecinos en el grupo vulnerable. Por eso las cascadas globales
suelen abarcar mucho más que los nodos vulnerables.

### 1.3 La ventana de cascadas

Para Erdős-Rényi, $p_k = e^{-z}z^k/k!$ y la condición queda

$$
F(z) \equiv \sum_{k=1}^{K} k(k-1)\thinspace e^{-z}\frac{z^k}{k!} - z > 0,
\qquad K = \lfloor 1/\varphi\rfloor .
$$

- Para $z$ pequeño, $F(z) \approx 2p_2 - z \approx z^2 - z < 0$ si $z < 1$:
  no hay cascadas porque no hay componente gigante.
- Para $z$ grande, los términos $e^{-z}z^k$ con $k \le K$ se vuelven
  despreciables y $F(z) \approx -z < 0$: casi no quedan nodos vulnerables.
- En medio, $F > 0$ si $K \ge 4$ ($\varphi \le 1/4$). Si $K \le 3$, es
  decir, $\varphi > 1/4$, la ventana desaparece.

Los bordes se encuentran numéricamente. Con $\varphi = 0.18$ ($K = 5$) son
$z \approx 1.02$ y $z \approx 5.76$; con $\varphi = 0.22$ ($K = 4$), 1.11 y
3.86. Como $K$ es entero, la ventana sólo cambia cuando $1/\varphi$ cruza un
entero.

**[Avanzado]** El borde inferior: en [Para profundizar, P.1](#pp1).

### 1.4 Robusto pero frágil

Cerca del borde superior, la probabilidad de que un nodo inicial dispare
una cascada global es pequeña: el grupo de nodos vulnerables apenas
percola. Pero si la dispara, la cascada encuentra una red muy conectada y
alcanza a casi todos. En el cuaderno, con $\varphi = 0.18$ y $z = 6$, la
probabilidad es alrededor de 0.14 y el tamaño de las cascadas globales es
casi 1.

Esta combinación tiene una consecuencia práctica: la historia de un sistema
así casi no tiene crisis pequeñas que avisen de las grandes. Los choques
locales se absorben una y otra vez, hasta que uno no se absorbe.

---

## 2. Contagio interbancario

### 2.1 El modelo de Gai y Kapadia

Cada banco $i$ tiene, como fracción de sus activos totales:

- préstamos interbancarios $\theta$, repartidos en partes iguales entre sus
  $k_i$ deudores;
- capital $K$: la pérdida que puede absorber sin quebrar.

Si $f_i$ de sus deudores quiebran sin pagar nada, el banco pierde
$\theta f_i/k_i$ y quiebra si

$$
\frac{f_i}{k_i} > \frac{K}{\theta} .
$$

**Proposición.** El modelo es una cascada de umbral con $\varphi = K/\theta$,
en una red **dirigida**: la quiebra de $j$ afecta a sus acreedores. Un banco
es vulnerable si $k_i < \theta/K$.

La versión dirigida de la condición de cascada usa la distribución conjunta
del número de deudores y de acreedores de cada banco (Gai y Kapadia, 2010).
Para una red dirigida al azar con $z$ deudores en promedio, el resultado es
cualitativamente el mismo: una ventana de contagio. En el cuaderno, con
$\theta = 0.2$ y $K = 0.04$ ($\varphi = 0.2$), la ventana va de $z \approx 1$
a $z \approx 7$.

### 2.2 Diversificación contra contagio

Aumentar el número de contrapartes $k$ tiene dos efectos:

1. **Diversificación**: la exposición a cada deudor, $\theta/k$, baja. Menos
   bancos son vulnerables.
2. **Conectividad**: hay más caminos para el contagio.

El primer efecto explica el borde superior de la ventana, y el segundo el
inferior. El argumento clásico (Allen y Gale, 2000) de que más conexiones
dan más estabilidad sólo considera el primero.

Las palancas regulatorias del modelo son subir $K$ (más capital) o bajar
$\theta$ (menos exposición interbancaria). Ambas suben $\varphi = K/\theta$.

---

## 3. DebtRank

### 3.1 Definición

Sea $A_{ij} \ge 0$ el préstamo del banco $i$ al banco $j$ y $E_i > 0$ el
capital de $i$. Sea $h_i \in [0, 1]$ la fracción del capital de $i$ que se
ha perdido. La **matriz de impacto** es

$$
\Lambda_{ij} = \frac{A_{ij}}{E_i} .
$$

**Dinámica** (DebtRank lineal; Bardoscia et al., 2015). Con $h(0) = 0$ y
$h(1)$ dado por el choque inicial,

$$
h_i(t+1) = \min\left\lbrace 1,\ h_i(t) + \sum_j \Lambda_{ij}\left[h_j(t) - h_j(t-1)\right]\right\rbrace .
$$

**Interpretación.** El préstamo de $i$ a $j$ se valora en
$A_{ij}(1 - h_j)$: pierde valor en la misma proporción que el capital del
deudor. Si $h_j$ aumenta en $\Delta h_j$, el capital de $i$ baja en
$A_{ij}\Delta h_j$, es decir, $h_i$ sube en $\Lambda_{ij}\Delta h_j$. El
mínimo con 1 dice que un banco no puede perder más que todo su capital.

**Pérdida del sistema.** $H = \sum_i E_i h_i / \sum_i E_i$. El **DebtRank del
banco** $j$ es la pérdida de los demás cuando $h_j(1) = 1$ y $h_i(1) = 0$
para $i \ne j$:

$$
R_j = \frac{\sum_{i} E_i h_i(\infty) - E_j}{\sum_i E_i} .
$$

**Sobre las versiones.** El DebtRank original (Battiston et al., 2012)
deja que cada banco transmita su deterioro una sola vez, lo que evita que
las pérdidas den vueltas en ciclos. La versión lineal de estas notas, que
es la del cuaderno, deja que reverberen.

### 3.2 Solución sin saturación

**Proposición.** Mientras ningún banco llegue a $h_i = 1$, los incrementos
$\Delta(t) = h(t) - h(t-1)$ cumplen $\Delta(t+1) = \Lambda\thinspace\Delta(t)$,
y

$$
h(t) = \sum_{s=0}^{t-1} \Lambda^s\thinspace h(1) .
$$

Si el radio espectral $r(\Lambda) < 1$, la serie converge a

$$
h(\infty) = (I - \Lambda)^{-1} h(1)
$$

([L.A4](repaso/algebra_lineal.md#la4)).

*Demostración.* Sin el mínimo, $h(t+1) - h(t) = \Lambda\left[h(t) - h(t-1)\right]$,
así que $\Delta(t) = \Lambda^{t-1}\Delta(1)$ con $\Delta(1) = h(1)$, y $h(t)$
es la suma de los incrementos. Si $r(\Lambda) < 1$, la serie de Neumann
converge. $\square$

**Corolario (apalancamiento uniforme).** Si cada banco presta $\ell$ veces
su capital, cada fila de $\Lambda$ suma $\ell$: $\Lambda\mathbf 1 = \ell\mathbf 1$.
Un choque uniforme $h(1) = s\mathbf 1$ da $\Lambda^t h(1) = \ell^t s\mathbf 1$,
y

$$
h(\infty) = \frac{s}{1 - \ell}\thinspace\mathbf 1 \quad (\ell < 1).
$$

La pérdida se amplifica por $1/(1-\ell)$: la serie geométrica de
[C.B5](repaso/calculo.md#cb5). En el cuaderno, con $\ell = 0.75$ y un choque
del 1% del capital, la pérdida final es del 4%.

### 3.3 [Avanzado] Estabilidad e inestabilidad

El desarrollo de esta subsección está en [Para profundizar, P.2](#pp2).


### 3.4 Demasiado central para caer

$R_j$ depende de cuánto le deben a $j$ (su tamaño como deudor), del capital
de sus acreedores y de lo que éstos a su vez transmiten. En el cuaderno,
con 50 bancos de tamaños Pareto y $\ell = 0.75$, la correlación entre
tamaño y $R_j$ es positiva pero lejos de 1: algunos bancos medianos tienen
un impacto mayor que el de los más grandes. Battiston et al. (2012)
encontraron, con los datos de los préstamos de emergencia de la Reserva
Federal, un grupo de unas 20 instituciones tan interconectadas que todas
eran sistémicamente importantes en el pico de la crisis.

---

## 4. [Avanzado] Eisenberg-Noe

El desarrollo de esta sección está en [Para profundizar, P.3](#pp3).


---

## 5. Resumen formal

| Concepto | Resultado |
|---|---|
| Cascada de Watts | Vulnerables: $k \le 1/\varphi$; cascada global si $\sum_{k \le 1/\varphi} k(k-1)p_k > z$ |
| Ventana en Erdős-Rényi | Desde $z \approx 1$ hasta un borde superior que depende de $\lfloor 1/\varphi\rfloor$; no hay ventana si $\varphi > 1/4$ |
| Gai-Kapadia | Cascada de umbral dirigida con $\varphi = K/\theta$ |
| DebtRank lineal | $h(\infty) = (I - \Lambda)^{-1}h(1)$ si $r(\Lambda) < 1$; con filas iguales a $\ell$, amplificación $1/(1-\ell)$ |
| Inestabilidad | $r(\Lambda) \ge 1$ y $\Lambda$ irreducible: cualquier choque produce quiebras |
| Eisenberg-Noe | $p = \Phi(p)$ con $\Phi$ creciente; iterar desde $\bar p$ da el mayor punto fijo |
| Anillo contra completa | Completa: 1 quiebra si $\varepsilon \le n\delta$, todas si $\varepsilon > n\delta$; anillo: quiebras $\approx \varepsilon/\delta$ |

---

## 6. Ejercicios formales

**Ruta básica**

- **F1.** Con $\varphi = 0.15$, ¿qué grados son vulnerables? Calcula el lado
  izquierdo de la condición de cascada para Erdős-Rényi con $z = 4$ y decide
  si hay cascadas globales.
- **F2.** Un banco tiene 15% de sus activos en préstamos a 5 bancos, en
  partes iguales, y 3% de capital. ¿Cuántos deudores tienen que quebrar
  para que quiebre? ¿Y si presta a 10 bancos?
- **F3.** Con $\ell = 0.8$, un choque uniforme del 2% del capital, ¿qué
  pérdida final produce según DebtRank lineal? ¿Cuántas rondas hacen falta
  para alcanzar el 90% de esa pérdida?
- **F4.** Explica por qué un sistema bancario puede pasar muchos años sin
  crisis y aun así ser frágil.

**Ruta avanzada**

- **F5.** Para Erdős-Rényi y $K = \lfloor 1/\varphi\rfloor = 3$, demuestra que
  $F(z) < 0$ para todo $z > 0$: no hay ventana.
- **F6.** Demuestra las cotas de la sección 3.3 para el radio espectral de
  una matriz no negativa. Pista: aplica $\Lambda$ al eigenvector de
  Perron-Frobenius y mira su componente más grande.
- **F7.** En el modelo de Eisenberg-Noe sin deuda externa ($v = 0$), con
  $n = 3$ bancos en anillo, $L_{12} = L_{23} = L_{31} = 1$ y
  $e = (0, 0, 0)$, encuentra todos los vectores de pagos. ¿Cuál da la
  iteración desde $\bar p$? ¿Se cumple la condición de regularidad?
- **F8.** En el anillo de la sección 4.3, ¿qué pasa si el choque supera la
  deuda interbancaria total del banco 1 ($\varepsilon - \delta > 1$)? ¿Quién
  absorbe el exceso?

---

## Para profundizar

Material opcional de la ruta avanzada: demostraciones y temas que no se ven
en clase. Cada apartado indica de qué sección de las notas viene.

<a id="pp1"></a>

### P.1 El borde inferior (sección 1.3)

**El borde inferior.** Desarrollando para $z$ pequeño,
$F(z) = z^2 - z + O(z^3)$, así que el borde inferior es $z \approx 1$: la
transición de percolación de Erdős-Rényi.

<a id="pp2"></a>

### P.2 Estabilidad e inestabilidad (sección 3.3)


**Proposición.** Si $r(\Lambda) \ge 1$ y $\Lambda$ es irreducible (cualquier
banco está expuesto, directa o indirectamente, a cualquier otro), entonces
para todo choque $h(1) \ge 0$ no nulo, por pequeño que sea, algún banco
llega a $h_i = 1$.

*Idea de la demostración.* Por Perron-Frobenius
([L.A5](repaso/algebra_lineal.md#la5)), $\Lambda$ tiene un eigenvector
izquierdo $u > 0$ con eigenvalor $r = r(\Lambda)$. Mientras no haya
saturación, $u^\top\Delta(t+1) = u^\top\Lambda\Delta(t) = r\thinspace u^\top\Delta(t)$,
así que $u^\top\Delta(t) = r^{t-1}u^\top h(1)$, y la suma
$u^\top h(t) = \sum_{s < t} r^s\thinspace u^\top h(1)$ crece sin límite (de forma
geométrica si $r > 1$, lineal si $r = 1$). Como
$u^\top h(t) \le u^\top\mathbf 1$ mientras $h \le \mathbf 1$, el crecimiento sin
saturación no puede durar. $\square$

**Cotas útiles.** Para una matriz no negativa, el radio espectral está entre
la menor y la mayor suma por filas:

$$
\min_i \sum_j \Lambda_{ij} \le r(\Lambda) \le \max_i \sum_j \Lambda_{ij} .
$$

La suma de la fila $i$ es el apalancamiento interbancario de $i$: sus
préstamos a otros bancos entre su capital. Si **todos** los bancos prestan
menos que su capital, $r(\Lambda) < 1$ y el sistema es estable en este
sentido; si todos prestan más, $r(\Lambda) > 1$, sin importar la forma de la
red. Con
apalancamientos distintos, la forma de la red decide.

<a id="pp3"></a>

### P.3 Eisenberg-Noe (sección 4)


#### 4.1 El vector de pagos

Sea $L_{ij} \ge 0$ la deuda de $i$ con $j$ (sin deudas consigo mismo),
$\bar p_i = \sum_j L_{ij}$ su deuda interbancaria total y
$\Pi_{ij} = L_{ij}/\bar p_i$ la proporción que le corresponde a cada
acreedor. Cada banco tiene activos externos $e_i \ge 0$ y una deuda externa
$v_i \ge 0$ con prioridad sobre la interbancaria.

Un **vector de pagos** es un $p \in [0, \bar p]$ con $p = \Phi(p)$, donde

$$
\Phi_i(p) = \max\left\lbrace 0,\ \min\left\lbrace \bar p_i,\ e_i - v_i + \sum_j \Pi_{ji}\thinspace p_j\right\rbrace\right\rbrace .
$$

Las reglas que codifica son: **responsabilidad limitada** (nadie paga más
de lo que tiene), **prioridad absoluta** (se paga toda la deuda antes de
que quede algo para los accionistas) y **proporcionalidad** (un banco en
quiebra paga a sus acreedores en proporción a lo que les debe).

#### 4.2 Existencia y cálculo

**Proposición** (Eisenberg y Noe, 2001). Existe un vector de pagos máximo
$p^+$ (y uno mínimo $p^-$). La sucesión $p^{(0)} = \bar p$,
$p^{(t+1)} = \Phi(p^{(t)})$ decrece y converge a $p^+$.

*Demostración.* $\Phi$ es **creciente**: si $p \le q$ componente a
componente, entonces $\Phi(p) \le \Phi(q)$, porque $\Pi_{ji} \ge 0$ y el
máximo y el mínimo conservan el orden. Además $\Phi$ manda el intervalo
$[0, \bar p]$ en sí mismo. Como $p^{(1)} = \Phi(\bar p) \le \bar p = p^{(0)}$,
aplicando $\Phi$ repetidamente $p^{(t+1)} \le p^{(t)}$. Una sucesión
decreciente y acotada por abajo converge; por continuidad de $\Phi$ el
límite es un punto fijo. Si $q$ es cualquier otro punto fijo, $q \le \bar p$
implica $q = \Phi^t(q) \le \Phi^t(\bar p) = p^{(t)}$ para todo $t$, así que
el límite es el mayor. (El teorema de Tarski da la existencia de $p^-$.) $\square$

Bajo una condición de regularidad (que de cada grupo de bancos que sólo se
deben entre sí salga algún activo externo positivo), el punto fijo es
único.

#### 4.3 Anillo contra red completa

Considera $n$ bancos, cada uno con deuda interbancaria total 1 y colchón
$\delta = e_i - v_i > 0$. Un choque reduce $e_1$ en $\varepsilon$.

**Red completa**, con $L_{ij} = 1/(n-1)$. Si sólo quiebra el banco 1, éste deja
de pagar $\varepsilon - \delta$ (si es menor que 1), y cada uno de los demás
pierde $(\varepsilon - \delta)/(n-1)$. Los demás sobreviven si esa pérdida
no supera su colchón:

$$
\frac{\varepsilon - \delta}{n - 1} \le \delta
\quad\Longleftrightarrow\quad
\varepsilon \le n\delta .
$$

Si $\varepsilon > n\delta$, el capital conjunto del sistema, $n\delta - \varepsilon$,
es negativo, así que no puede ser que sólo el banco 1 quiebre; por simetría,
los otros $n - 1$ quiebran juntos. Con $n = 10$ y $\delta = 0.08$, el salto
de 1 a 10 quiebras ocurre en $\varepsilon = 0.8$.

**Anillo**, con $L_{i,i+1} = 1$. El banco 1 deja de pagar $\varepsilon - \delta$
al banco 2, que absorbe $\delta$ y deja de pagar $\varepsilon - 2\delta$ al 3,
y así. El banco $j$ quiebra si $\varepsilon > j\delta$: el número de quiebras
crece de uno en uno con $\varepsilon$.

**Conclusión** (Acemoglu, Ozdaglar y Tahbaz-Salehi, 2015). Para choques
pequeños, la red completa es la más estable y el anillo la menos estable.
Para choques mayores que el colchón total, la red completa se vuelve tan
frágil como el anillo: la conexión que reparte los choques pequeños también
reparte los grandes. Los autores muestran que, en ese régimen, las redes
débilmente conectadas son más estables que ambas.

---

## Referencias

- Acemoglu, D., Ozdaglar, A., & Tahbaz-Salehi, A. (2015). Systemic risk and
  stability in financial networks. *American Economic Review*, 105(2),
  564–608. https://doi.org/10.1257/aer.20130456
- Allen, F., & Gale, D. (2000). Financial contagion. *Journal of Political
  Economy*, 108(1), 1–33. https://doi.org/10.1086/262109
- Bardoscia, M., Battiston, S., Caccioli, F., & Caldarelli, G. (2015).
  DebtRank: A microscopic foundation for shock propagation. *PLoS ONE*,
  10(6), e0130406. https://doi.org/10.1371/journal.pone.0130406
- Bardoscia, M., Battiston, S., Caccioli, F., & Caldarelli, G. (2017).
  Pathways towards instability in financial networks. *Nature
  Communications*, 8, 14416. https://doi.org/10.1038/ncomms14416
- Battiston, S., Puliga, M., Kaushik, R., Tasca, P., & Caldarelli, G.
  (2012). DebtRank: Too central to fail? Financial networks, the FED and
  systemic risk. *Scientific Reports*, 2, 541.
  https://doi.org/10.1038/srep00541
- Battiston, S., Farmer, J. D., Flache, A., Garlaschelli, D., Haldane, A.
  G., Heesterbeek, H., Hommes, C., Jaeger, C., May, R., & Scheffer, M.
  (2016). Complexity theory and financial regulation. *Science*, 351(6275),
  818–819. https://doi.org/10.1126/science.aad0299
- Eisenberg, L., & Noe, T. H. (2001). Systemic risk in financial systems.
  *Management Science*, 47(2), 236–249.
  https://doi.org/10.1287/mnsc.47.2.236.9835
- Gai, P., & Kapadia, S. (2010). Contagion in financial networks.
  *Proceedings of the Royal Society A*, 466(2120), 2401–2423.
  https://doi.org/10.1098/rspa.2009.0410
- Glasserman, P., & Young, H. P. (2016). Contagion in financial networks.
  *Journal of Economic Literature*, 54(3), 779–831.
  https://doi.org/10.1257/jel.20151228
- Watts, D. J. (2002). A simple model of global cascades on random
  networks. *Proceedings of the National Academy of Sciences*, 99(9),
  5766–5771. https://doi.org/10.1073/pnas.082090499
