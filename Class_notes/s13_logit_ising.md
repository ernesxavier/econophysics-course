# Sesiones 13 y 14: Ising y elección discreta. Notas formales

> **Cómo leer estas notas.** Siguen el orden del programa: primero se observa
> en el cuaderno (`Code/s13_logit_ising.jl`), después se formaliza aquí. Lo que no
> tiene marca es **ruta básica**. Las secciones y demostraciones marcadas
> **[Avanzado]** son opcionales.
>
> Repaso recomendado: [C.B4](repaso/calculo.md#cb4) (función logística),
> [F.B3](repaso/fisica_estadistica.md#fb3) (temperatura),
> [F.B4](repaso/fisica_estadistica.md#fb4) (transiciones de fase),
> [F.I1](repaso/fisica_estadistica.md#fi1) (Boltzmann),
> [F.I2](repaso/fisica_estadistica.md#fi2) (logit = Boltzmann),
> [F.I3](repaso/fisica_estadistica.md#fi3) (Ising),
> [F.I4](repaso/fisica_estadistica.md#fi4) (campo medio),
> [F.I5](repaso/fisica_estadistica.md#fi5) (Metropolis),
> [C.I1](repaso/calculo.md#ci1) (Taylor), [C.I5](repaso/calculo.md#ci5)
> (puntos fijos); para la ruta avanzada,
> [E.A3](repaso/procesos_estocasticos.md#ea3) (ecuación maestra),
> [C.A4](repaso/calculo.md#ca4) (Stirling),
> [F.A1](repaso/fisica_estadistica.md#fa1) (Landau) y
> [F.A3](repaso/fisica_estadistica.md#fa3) (campo aleatorio).
>
> Lecturas: Sen & Chakrabarti (2014), caps. 2–3; Requião da Cunha (2021),
> cap. 6. Ruta avanzada: Brock & Durlauf (2001); Bouchaud, Marsili & Nadal
> (2023).

---

## 1. El modelo logit

### 1.1 Utilidad aleatoria

Un agente elige una opción de un conjunto finito. La opción $i$ le da
utilidad $U_i = u_i + \varepsilon_i$, donde $u_i$ es la parte que observa el
analista y $\varepsilon_i$ un choque que el analista no observa. El agente
elige la opción con mayor $U_i$.

**Distribución de Gumbel.** Una variable tiene distribución de Gumbel con
escala $1/\beta$ (y posición 0) si su CDF es

$$
F(x) = \exp\left(-e^{-\beta x}\right), \qquad x \in \mathbb R .
$$

Es la distribución límite del máximo de muchas variables con colas de tipo
exponencial (teoría de valores extremos).

**Proposición** (McFadden). Si los $\varepsilon_i$ son independientes con
distribución de Gumbel de escala $1/\beta$,

$$
\mathbb P(\text{elige } i) = \frac{e^{\beta u_i}}{\sum_j e^{\beta u_j}} .
$$

*Demostración.* La opción $i$ gana si $\varepsilon_j < u_i - u_j + \varepsilon_i$
para todo $j \ne i$. Condicionando en $\varepsilon_i = x$, con densidad
$f(x) = \beta e^{-\beta x}F(x)$,

$$
\mathbb P(i) = \int_{-\infty}^{\infty} \beta e^{-\beta x} e^{-e^{-\beta x}}
\prod_{j \ne i} \exp\left(-e^{-\beta(x + u_i - u_j)}\right) dx
= \int_{-\infty}^{\infty} \beta e^{-\beta x}
\exp\left(-e^{-\beta x}\thinspace\frac{S}{e^{\beta u_i}}\right) dx ,
$$

con $S = \sum_j e^{\beta u_j}$ (el término $j = i$ es el primer factor). Con
el cambio $y = e^{-\beta x}$, $dy = -\beta y\thinspace dx$:

$$
\mathbb P(i) = \int_0^\infty \exp\left(-y\thinspace\frac{S}{e^{\beta u_i}}\right) dy
= \frac{e^{\beta u_i}}{S} . \qquad \square
$$

**Caso binario.** Con dos opciones,
$\mathbb P(A) = \sigma\left(\beta(u_A - u_B)\right)$, con
$\sigma(x) = 1/(1 + e^{-x})$ ([C.B4](repaso/calculo.md#cb4)).

### 1.2 Logit = Boltzmann

La distribución de Boltzmann a temperatura $T = 1/\beta$
([F.I1](repaso/fisica_estadistica.md#fi1)) asigna al estado $i$, de energía
$E_i$, probabilidad $e^{-\beta E_i}/Z$. Es la fórmula logit con
$E_i = -u_i$. El diccionario:

| Física estadística | Elección discreta |
|---|---|
| Estado $i$ | Opción $i$ |
| Energía $E_i$ | Menos la utilidad, $-u_i$ |
| Temperatura $T = 1/\beta$ | Escala de los choques no observados |
| Función de partición $Z$ | Denominador $\sum_j e^{\beta u_j}$ |
| $\ln Z/\beta$ | Excedente esperado del consumidor (fórmula *log-sum*), salvo una constante |

La última línea dice que $\mathbb E\left[\max_i U_i\right] = \beta^{-1}\ln\sum_j e^{\beta u_j}$
más la media de una Gumbel de posición 0 ($\gamma_E/\beta$, con
$\gamma_E \approx 0.577$ la constante de Euler-Mascheroni). Es la energía
libre con el signo cambiado ([F.A1](repaso/fisica_estadistica.md#fa1)).

---

## 2. El modelo de Ising

### 2.1 Definición

$N$ agentes en los nodos de una red eligen $s_i \in \lbrace -1, +1\rbrace$. La
energía de una configuración $s$ es

$$
H(s) = -J\sum_{\langle i,j\rangle} s_is_j - h\sum_i s_i ,
$$

donde $\langle i,j\rangle$ recorre los pares de vecinos, $J > 0$ premia
coincidir y $h$ es un incentivo común. En equilibrio, la configuración
aparece con probabilidad de Boltzmann
$\mathbb P(s) = e^{-\beta H(s)}/Z$. La **magnetización**
$m = \frac1N\sum_i s_i$ es la adopción neta.

**Simetría.** Con $h = 0$, $H(s) = H(-s)$: las dos opciones son
equivalentes y $\mathbb E[m] = 0$ para cualquier $N$ finito. La **ruptura
espontánea de simetría** sólo se ve en el límite $N \to \infty$: debajo de
$T_c$, el sistema se queda en uno de los dos consensos y nunca pasa al otro.
En una simulación finita, el sistema salta entre ambos de vez en cuando;
por eso se mide $|m|$.

### 2.2 La solución exacta en la red cuadrada

En la red cuadrada infinita con $h = 0$ (Onsager, 1944; Yang, 1952):

$$
T_c = \frac{2J}{\ln\left(1 + \sqrt 2\right)} \approx 2.269\thinspace J,
\qquad
m(T) = \left[1 - \sinh^{-4}(2J/T)\right]^{1/8} \quad (T < T_c),
$$

y $m = 0$ para $T \ge T_c$. Cerca de $T_c$, $m \sim (T_c - T)^{1/8}$: el
exponente 1/8 no es el del campo medio, 1/2
([F.A2](repaso/fisica_estadistica.md#fa2)).

### 2.3 Dinámicas: Metropolis y baño térmico

Ambas son cadenas de Markov que cambian un agente a la vez. Sea
$h_i^{\text{loc}} = J\sum_{j \sim i} s_j + h$ el campo que siente $i$ y
$\Delta E = 2s_ih_i^{\text{loc}}$ el cambio de energía si $i$ cambia de opción.

- **Metropolis**: cambiar con probabilidad $\min\lbrace 1, e^{-\beta\Delta E}\rbrace$.
- **Baño térmico** (Glauber, 1963): elegir $s_i = +1$ con probabilidad
  $\sigma\left(2\beta h_i^{\text{loc}}\right)$, sin importar el valor anterior.
  Es exactamente la regla **logit**: la ventaja en utilidad de $+1$ sobre
  $-1$ es $2h_i^{\text{loc}}$.

**Proposición.** Las dos dinámicas cumplen **balance detallado** respecto a
la distribución de Boltzmann: si $s'$ difiere de $s$ sólo en el agente $i$,

$$
\mathbb P(s)\thinspace W(s \to s') = \mathbb P(s')\thinspace W(s' \to s) .
$$

Por lo tanto la distribución de Boltzmann es estacionaria
([E.A3](repaso/procesos_estocasticos.md#ea3)).

*Demostración.* Sea $\Delta E = H(s') - H(s)$, así que
$\mathbb P(s')/\mathbb P(s) = e^{-\beta\Delta E}$. Metropolis: si
$\Delta E > 0$, $W(s \to s') = e^{-\beta\Delta E}$ y $W(s' \to s) = 1$; si
$\Delta E \le 0$, al revés. En ambos casos
$W(s \to s')/W(s' \to s) = e^{-\beta\Delta E}$. Baño térmico: si $s'_i = +1$,
$W(s \to s') = \sigma(2\beta h)$ y $W(s' \to s) = \sigma(-2\beta h)$, con
$h = h_i^{\text{loc}}$, y $\sigma(x)/\sigma(-x) = e^x$ da
$e^{2\beta h} = e^{-\beta\Delta E}$. (Los factores $1/N$ por elegir al agente
$i$ se cancelan.) $\square$

**Interpretación económica.** Una población de agentes que revisan su
elección de vez en cuando con una regla logit, mirando lo que hacen sus
vecinos, termina distribuida según el modelo de Ising. Blume (1993) lo
estudió como dinámica de juegos de coordinación.

---

## 3. Campo medio

### 3.1 Interacción con el promedio

Supón que cada agente interactúa con el promedio de toda la población:

$$
H(s) = -\frac{J}{2N}\sum_{i \ne j} s_is_j - h\sum_i s_i
\approx -N\left(\frac{J}{2}m^2 + hm\right) .
$$

Cada agente siente el mismo campo, $h + Jm$ (salvo un término de orden
$1/N$). Si $m$ es fijo, cada elección es un logit independiente:

$$
\mathbb E[s_i] = \tanh\left(\beta(h + Jm)\right) .
$$

**Autoconsistencia.** Si $N$ es grande, el promedio de las elecciones es
$m$ (ley de los grandes números), así que

$$
m = \tanh\left(\beta(h + Jm)\right) .
$$

### 3.2 [Avanzado] Derivación con la energía libre

El número de configuraciones con $N_+ = N(1+m)/2$ agentes en $+1$ es
$\binom{N}{N_+}$. Con Stirling ([C.A4](repaso/calculo.md#ca4)),
$\ln\binom{N}{N_+} \approx N\thinspace\mathcal S(m)$ con

$$
\mathcal S(m) = -\frac{1+m}{2}\ln\frac{1+m}{2} - \frac{1-m}{2}\ln\frac{1-m}{2} .
$$

La probabilidad de que la adopción neta valga $m$ es entonces

$$
\mathbb P(m) \propto \binom{N}{N_+}e^{-\beta H} \approx e^{-N\beta f(m)},
\qquad
f(m) = -\frac{J}{2}m^2 - hm - T\thinspace\mathcal S(m) .
$$

Para $N$ grande, $\mathbb P(m)$ se concentra en los mínimos de la **energía
libre** por agente $f$. Como $\mathcal S'(m) = -\operatorname{atanh}(m)$,

$$
f'(m) = -Jm - h + T\operatorname{atanh}(m) = 0
\quad\Longleftrightarrow\quad
m = \tanh\left(\beta(h + Jm)\right) .
$$

Además, $f''(m) = -J + T/(1 - m^2)$, así que una solución es un **mínimo**
(equilibrio estable) si $\beta J(1 - m^2) < 1$, y un máximo si
$\beta J(1 - m^2) > 1$. Es la condición gráfica: la curva
$\tanh(\beta(h + Jm))$ cruza la recta $y = m$ con pendiente menor que 1.

### 3.3 Análisis de la ecuación

**Proposición.** Con $h = 0$:

1. Si $\beta J \le 1$, la única solución es $m = 0$.
2. Si $\beta J > 1$, hay tres soluciones: $m = 0$, inestable, y $\pm m^\ast$,
   estables, con $m^\ast > 0$.

*Demostración.* $g(m) = \tanh(\beta Jm)$ es impar, creciente, cóncava para
$m > 0$ y acotada por 1. Su pendiente en 0 es $\beta J$. Si $\beta J \le 1$,
la concavidad da $g(m) < \beta Jm \le m$ para $m > 0$: no hay cruce
positivo. Si $\beta J > 1$, $g(m) > m$ cerca de 0 y $g(1) < 1$: por el
teorema del valor intermedio hay un cruce en $(0, 1)$, único por la
concavidad. Por simetría hay otro en $-m^\ast$. $\square$

**Cerca de la transición.** Con $\tanh x \approx x - x^3/3$
([C.I1](repaso/calculo.md#ci1)):

$$
m \approx \beta Jm - \frac{(\beta Jm)^3}{3}
\quad\Longrightarrow\quad
m^\ast \approx \sqrt{\frac{3(\beta J - 1)}{(\beta J)^3}} \approx \sqrt{3(\beta J - 1)} .
$$

Así $m^\ast \sim (\beta J - 1)^{1/2}$: el exponente de campo medio, que
coincide con la teoría de Landau ([F.A1](repaso/fisica_estadistica.md#fa1)).

**En una red.** Si cada agente tiene $z$ vecinos y se reemplazan por su
promedio, el campo es $h + Jzm$ y la transición queda en $T_c^{\text{cm}} = Jz$.
En la red cuadrada, $Jz = 4J$ contra el exacto $2.269J$: el campo medio
sobreestima el orden porque ignora que los vecinos de un agente están
correlacionados entre sí. Cuando todos interactúan con todos, el campo
medio es exacto para $N \to \infty$, como en el cuaderno.

---

## 4. Brock y Durlauf

### 4.1 El modelo

Brock y Durlauf (2001) escriben la utilidad del agente $i$ por la elección
$\omega_i \in \lbrace -1, +1\rbrace$ como

$$
V(\omega_i) = h\omega_i + J\omega_i\thinspace m^e_i + \varepsilon(\omega_i),
$$

con $m^e_i$ la adopción neta que el agente **espera** de los demás y
$\varepsilon(\pm 1)$ choques independientes de valores extremos tales que
$\varepsilon(-1) - \varepsilon(+1)$ es logístico con parámetro $\beta$. La
ventaja de adoptar es $2(h + Jm^e_i)$, y

$$
\mathbb P(\omega_i = +1) = \sigma\left(2\beta(h + Jm^e_i)\right),
\qquad
\mathbb E[\omega_i] = \tanh\left(\beta(h + Jm^e_i)\right) .
$$

**Expectativas racionales.** Si todos esperan lo mismo y aciertan,
$m^e_i = m$ y $m = \tanh(\beta h + \beta Jm)$: la ecuación de campo medio. Los
equilibrios del modelo económico son las soluciones de la sección 3.

### 4.2 Resultados

**Proposición** (Brock y Durlauf, 2001).

1. Si $\beta J < 1$, hay un único equilibrio.
2. Si $\beta J > 1$ y $h = 0$, hay tres equilibrios.
3. Si $\beta J > 1$, existe $H > 0$ tal que hay tres equilibrios si
   $|\beta h| < H$ y uno solo si $|\beta h| > H$.

Con tres equilibrios y $h > 0$, uno de los equilibrios estables tiene
$m < 0$: la mayoría elige contra su incentivo privado porque los demás lo
hacen. Es una **trampa de coordinación**.

**Multiplicador social.** Derivando la ecuación en un equilibrio estable,

$$
\frac{dm}{dh} = \frac{\beta(1 - m^2)}{1 - \beta J(1 - m^2)} .
$$

El numerador es el efecto de $h$ sin interacciones; el denominador, menor
que 1, lo amplifica. Un subsidio mueve la adopción más de lo que movería
a agentes aislados, porque cada adopción hace más atractiva la siguiente.
Cerca del punto donde un equilibrio desaparece, el denominador tiende a 0
y el multiplicador diverge.

**Identificación.** Si sólo se observa que agentes del mismo grupo hacen lo
mismo, no se puede distinguir la interacción social de características
comunes del grupo (el "problema de reflexión" de Manski, 1993). La no
linealidad del logit es la que, en el modelo de Brock y Durlauf, permite
identificar $J$ en principio.

---

## 5. Histéresis

Fija $\beta = 1$ y $J > 1$. La rama inferior (con $m < 0$) existe mientras la
curva $\tanh(h + Jm)$ toque a la recta $y = m$ por abajo. Desaparece cuando
se vuelven tangentes: $J(1 - m^2) = 1$, es decir, en
$m_0 = -\sqrt{1 - 1/J}$. Sustituyendo en $h = \operatorname{atanh}(m) - Jm$,

$$
h_s = J\sqrt{1 - 1/J} - \operatorname{atanh}\left(\sqrt{1 - 1/J}\right) .
$$

Con $J = 1.5$, $h_s \approx 0.21$. Al subir $h$ lentamente, el sistema salta
de la rama inferior a la superior en $h = h_s$; al bajarlo, salta de vuelta
en $h = -h_s$. Entre $-h_s$ y $h_s$, el estado depende de la historia.
Si $J \le 1$ no hay ciclo: la respuesta es continua y reversible.

En una simulación con $N$ finito y $h$ cambiando a velocidad finita, los
saltos se retrasan: cerca de $h_s$ la rama es casi plana y el sistema tarda
mucho en dejarla (frenado crítico).

---

## 6. [Avanzado] Campo aleatorio y avalanchas

**Modelo.** A temperatura cero, el agente $i$ adopta si
$h + h_i + Jm > 0$, con $h_i$ independientes con distribución
$\mathcal N(0, R^2)$. En campo medio y con $N$ grande, la fracción que
adopta es $\mathbb P(h_i > -h - Jm)$, así que

$$
m = 2\Phi\left(\frac{h + Jm}{R}\right) - 1 ,
$$

con $\Phi$ la CDF normal estándar.

**Proposición.** La curva $m(h)$ tiene un salto si y sólo si
$R < R_c = J\sqrt{2/\pi} \approx 0.80\thinspace J$.

*Demostración.* El lado derecho, como función de $m$, tiene pendiente
$2J\varphi\left((h + Jm)/R\right)/R$, que es máxima en el argumento 0 y vale
$2J/(R\sqrt{2\pi})$. Si es mayor que 1 para algún $h$, hay un rango de $h$ con
tres soluciones y, al subir $h$, la rama inferior desaparece: salto. Si es
menor o igual que 1, la solución es única y continua en $h$. La condición
$2J/(R\sqrt{2\pi}) > 1$ es $R < J\sqrt{2/\pi}$. $\square$

**Avalanchas.** Con $N$ finito, al subir $h$ cada adopción sube el campo de
todos en $2J/N$ y puede provocar otras. En el punto crítico $R = R_c$, el
tamaño de las avalanchas sigue una ley de potencias; en campo medio su
densidad decae como $s^{-3/2}$ (Sethna et al., 1993). Lejos de $R_c$, las
avalanchas tienen un tamaño típico, salvo la avalancha macroscópica del
salto cuando $R < R_c$.

**Lectura económica** (Bouchaud, 2013). $R$ mide la diversidad de
opiniones o situaciones privadas. Una población diversa responde de forma
suave a los cambios de incentivos; una homogénea puede cambiar de golpe:
pánicos, cambios de normas, caídas del consumo.

---

## 7. Resumen formal

| Concepto | Resultado |
|---|---|
| Logit | Choques de Gumbel con escala $1/\beta$ dan $\mathbb P(i) = e^{\beta u_i}/\sum_j e^{\beta u_j}$ |
| Logit = Boltzmann | $E_i = -u_i$, $T = 1/\beta$, $Z$ = denominador del logit |
| Ising en red cuadrada | $T_c = 2J/\ln(1+\sqrt2)$; $m \sim (T_c - T)^{1/8}$ |
| Dinámicas | Metropolis y baño térmico (logit) cumplen balance detallado |
| Campo medio | $m = \tanh(\beta(h + Jm))$; transición en $\beta J = 1$; $m^\ast \approx \sqrt{3(\beta J - 1)}$ |
| Estabilidad | Equilibrio estable si $\beta J(1 - m^2) < 1$ |
| Multiplicador social | $dm/dh = \beta(1-m^2)/(1 - \beta J(1-m^2))$ |
| Histéresis ($\beta = 1$) | Saltos en $\pm h_s$, con $h_s = J\sqrt{1-1/J} - \operatorname{atanh}\sqrt{1-1/J}$ |
| Campo aleatorio | Salto macroscópico si $R < J\sqrt{2/\pi}$ |

---

## 8. Ejercicios formales

**Ruta básica**

- **F1.** Un consumidor elige entre tres marcas con utilidades 2, 1 y 0 y
  $\beta = 1$. Calcula las probabilidades logit. ¿Qué pasa si se agrega una
  cuarta marca idéntica a la segunda?
- **F2.** Con $h = 0$ y $\beta J = 1.5$, resuelve $m = \tanh(1.5m)$ iterando
  desde $m = 1$. Compara con la aproximación $\sqrt{3(\beta J - 1)}$.
- **F3.** Calcula el multiplicador social con $\beta = 1$, $J = 0.8$ y
  $h = 0$. ¿Cuánto más responde la adopción a un subsidio que sin
  interacciones?
- **F4.** Explica por qué en la red cuadrada el campo medio predice
  $T_c = 4J$ y la solución exacta es menor.

**Ruta avanzada**

- **F5.** Demuestra que si $\varepsilon_A$ y $\varepsilon_B$ son Gumbel
  independientes de escala $1/\beta$, entonces $\varepsilon_A - \varepsilon_B$
  tiene CDF $\sigma(\beta x)$.
- **F6.** Demuestra que $\ln Z/\beta$ es, salvo una constante, la utilidad
  esperada del agente: $\mathbb E\left[\max_i U_i\right]$.
- **F7.** Encuentra $H$ de la proposición de la sección 4.2 como función de
  $\beta J$.
- **F8.** Para el modelo de la sección 6, encuentra el valor de $h$ donde
  ocurre el salto cuando $R < R_c$, al menos numéricamente para $J = 1$ y
  $R = 0.6$, y compáralo con el cuaderno.
- **F9.** Verifica el balance detallado del baño térmico para un modelo de
  Ising en una red de 3 agentes en triángulo, escribiendo las 8
  configuraciones y sus probabilidades.

---

## Referencias

- Blume, L. E. (1993). The statistical mechanics of strategic interaction.
  *Games and Economic Behavior*, 5(3), 387–424.
  https://doi.org/10.1006/game.1993.1023
- Bouchaud, J.-P. (2013). Crises and collective socio-economic phenomena:
  Simple models and challenges. *Journal of Statistical Physics*, 151(3–4),
  567–606. https://doi.org/10.1007/s10955-012-0687-3
- Bouchaud, J.-P., Marsili, M., & Nadal, J.-P. (2023). Application of spin
  glass ideas in social sciences, economics and finance. En *Spin Glass
  Theory and Far Beyond* (pp. 561–579). World Scientific.
  https://doi.org/10.1142/9789811273926_0027
- Brock, W. A., & Durlauf, S. N. (2001). Discrete choice with social
  interactions. *Review of Economic Studies*, 68(2), 235–260.
  https://doi.org/10.1111/1467-937X.00168
- Glaeser, E. L., Sacerdote, B., & Scheinkman, J. A. (1996). Crime and
  social interactions. *Quarterly Journal of Economics*, 111(2), 507–548.
  https://doi.org/10.2307/2946686
- Glauber, R. J. (1963). Time-dependent statistics of the Ising model.
  *Journal of Mathematical Physics*, 4(2), 294–307.
  https://doi.org/10.1063/1.1703954
- Manski, C. F. (1993). Identification of endogenous social effects: The
  reflection problem. *Review of Economic Studies*, 60(3), 531–542.
  https://doi.org/10.2307/2298123
- McFadden, D. (1974). Conditional logit analysis of qualitative choice
  behavior. En P. Zarembka (Ed.), *Frontiers in Econometrics*
  (pp. 105–142). Academic Press.
- Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., &
  Teller, E. (1953). Equation of state calculations by fast computing
  machines. *Journal of Chemical Physics*, 21(6), 1087–1092.
  https://doi.org/10.1063/1.1699114
- Onsager, L. (1944). Crystal statistics. I. A two-dimensional model with
  an order-disorder transition. *Physical Review*, 65(3–4), 117–149.
  https://doi.org/10.1103/PhysRev.65.117
- Sen, P., & Chakrabarti, B. K. (2014). *Sociophysics: An Introduction*.
  Oxford University Press.
- Sethna, J. P., Dahmen, K., Kartha, S., Krumhansl, J. A., Roberts, B. W.,
  & Shore, J. D. (1993). Hysteresis and hierarchies: Dynamics of
  disorder-driven first-order phase transformations. *Physical Review
  Letters*, 70(21), 3347–3350. https://doi.org/10.1103/PhysRevLett.70.3347
- Yang, C. N. (1952). The spontaneous magnetization of a two-dimensional
  Ising model. *Physical Review*, 85(5), 808–816.
  https://doi.org/10.1103/PhysRev.85.808
