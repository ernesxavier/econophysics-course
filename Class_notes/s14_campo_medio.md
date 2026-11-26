# Sesión 14: Campo medio, histéresis y crisis. Notas formales

> **Cómo usar estas notas.** La clase combina el cuaderno
> (`Code/s14_campo_medio.jl`) y el pizarrón. El **plan de clase** de abajo
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
> Las secciones marcadas con ★ son el **núcleo** que se ve en clase; las
> demás se leen por cuenta propia y se evalúan igual.
>
> Repaso recomendado: [P.B4](repaso/probabilidad.md#pb4) (distribución normal),
> [C.I1](repaso/calculo.md#ci1) (Taylor),
> [C.I5](repaso/calculo.md#ci5) (puntos fijos),
> [F.I4](repaso/fisica_estadistica.md#fi4) (campo medio); para la ruta
> avanzada, [C.A4](repaso/calculo.md#ca4) (Stirling),
> [F.A1](repaso/fisica_estadistica.md#fa1) (Landau),
> [F.A2](repaso/fisica_estadistica.md#fa2) (exponentes críticos) y
> [F.A3](repaso/fisica_estadistica.md#fa3) (campo aleatorio).
>
> Lecturas: Requião da Cunha (2021), cap. 6. Ruta avanzada: Brock & Durlauf
> (2001); Bouchaud, Marsili & Nadal (2023).

**Plan de clase** (120 minutos: 105 de núcleo y 15 de margen). El núcleo:
la ecuación de campo medio y su solución gráfica, el modelo de Brock y
Durlauf y la histéresis.

| Minutos | Dónde | Qué | Secciones |
|---|---|---|---|
| 0–8 | Proyector | El fenómeno: cambios abruptos que no se revierten al deshacer la causa | Cuaderno, sección 0 |
| 8–22 | Pizarra | Campo medio: interacción con el promedio y autoconsistencia; la solución gráfica de $m = \tanh(\beta(h + Jm))$; la proposición, con la demostración; cerca de la transición, $m^\ast \sim (\beta J - 1)^{1/2}$ (la cuenta es el ejercicio 6 del cuaderno) | 1.1, 1.3 |
| 22–27 | Pizarra | En una red, el campo medio sobreestima el orden: $T_c^{\text{cm}} = Jz$ contra el exacto | 1.3 |
| 27–38 | Pizarra | Brock y Durlauf: la utilidad, las expectativas racionales y la ecuación de campo medio; equilibrios múltiples y la trampa de coordinación; el multiplicador social | 2.1–2.2 |
| 38–46 | Pizarra | Histéresis: la tangencia, el salto en $\pm h_s$ y la dependencia de la historia | 3 |
| 46–90 | Parejas | Laboratorio: la simulación frente a la solución gráfica y el ciclo de histéresis; ejercicio de código 1 | Cuaderno, secciones 1 y 2 y ejercicio 1 |
| 90–105 | Pizarra | Cierre: puesta en común de las preguntas 3 a 5 y la tabla del resumen | Cuaderno, preguntas 3–5; notas, 5 |
| 105–120 | Margen | Si alcanza: el ejercicio de código 2 (el multiplicador social) y cómo cambia cerca del salto | |
| Después | Lectura | Lo que no entró del margen; la identificación (2.2), el frenado crítico (3) y Para profundizar | |

Lo del **margen** se ve sólo si sobra tiempo, en el orden indicado; lo
que no se vea pasa a la lectura y se evalúa igual.

**Pregunta de la sesión.** ¿Por qué un cambio gradual en los incentivos puede
producir un salto de comportamiento que no se revierte al deshacer la causa?

---

## 1. Campo medio

### 1.1 Interacción con el promedio ★

Supón que cada agente interactúa con el promedio de toda la población
([F.I4](repaso/fisica_estadistica.md#fi4)):

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
\boxed{\ m = \tanh\left(\beta(h + Jm)\right)\ }
$$

Es una ecuación de punto fijo ([C.I5](repaso/calculo.md#ci5)).

### 1.2 [Avanzado] Derivación con la energía libre

El desarrollo de esta subsección está en [Para profundizar, P.1](#pp1).


### 1.3 Análisis de la ecuación ★

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
coincide con la teoría de Landau ([F.A1](repaso/fisica_estadistica.md#fa1))
y difiere del exponente exacto $1/8$ de la red cuadrada
([F.A2](repaso/fisica_estadistica.md#fa2)).

**En una red.** Si cada agente tiene $z$ vecinos y se reemplazan por su
promedio, el campo es $h + Jzm$ y la transición queda en $T_c^{\text{cm}} = Jz$.
En la red cuadrada, $Jz = 4J$ contra el exacto $2.269J$: el campo medio
sobreestima el orden porque ignora que los vecinos de un agente están
correlacionados entre sí. Cuando todos interactúan con todos, el campo
medio es exacto para $N \to \infty$, como en el cuaderno.

---

## 2. Brock y Durlauf ★

### 2.1 El modelo ★

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
equilibrios del modelo económico son las soluciones de la sección 1.

### 2.2 Resultados ★

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
\boxed{\ \frac{dm}{dh} = \frac{\beta(1 - m^2)}{1 - \beta J(1 - m^2)}\ }
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

## 3. Histéresis ★

Fija $\beta = 1$ y $J > 1$. La rama inferior (con $m < 0$) existe mientras la
curva $\tanh(h + Jm)$ toque a la recta $y = m$ por abajo. Desaparece cuando
se vuelven tangentes: $J(1 - m^2) = 1$, es decir, en
$m_0 = -\sqrt{1 - 1/J}$. Sustituyendo en $h = \operatorname{atanh}(m) - Jm$,

$$
\boxed{\ h_s = J\sqrt{1 - 1/J} - \operatorname{atanh}\left(\sqrt{1 - 1/J}\right)\ }
$$

Con $J = 1.5$, $h_s \approx 0.21$. Al subir $h$ lentamente, el sistema salta
de la rama inferior a la superior en $h = h_s$; al bajarlo, salta de vuelta
en $h = -h_s$. Entre $-h_s$ y $h_s$, el estado depende de la historia.
Si $J \le 1$ no hay ciclo: la respuesta es continua y reversible.

En una simulación con $N$ finito y $h$ cambiando a velocidad finita, los
saltos se retrasan: cerca de $h_s$ la rama es casi plana y el sistema tarda
mucho en dejarla (frenado crítico).

---

## 4. [Avanzado] Campo aleatorio y avalanchas

El desarrollo de esta sección está en [Para profundizar, P.2](#pp2).


---

## 5. Resumen formal

| Concepto | Resultado |
|---|---|
| Campo medio | $m = \tanh(\beta(h + Jm))$; transición en $\beta J = 1$; $m^\ast \approx \sqrt{3(\beta J - 1)}$ |
| Estabilidad | Equilibrio estable si $\beta J(1 - m^2) < 1$ |
| Multiplicador social | $dm/dh = \beta(1-m^2)/(1 - \beta J(1-m^2))$ |
| Histéresis ($\beta = 1$) | Saltos en $\pm h_s$, con $h_s = J\sqrt{1-1/J} - \operatorname{atanh}\sqrt{1-1/J}$ |
| Campo aleatorio | Salto macroscópico si $R < J\sqrt{2/\pi}$ |

---

## 6. Ejercicios formales

**Ruta básica**

- **F1.** Con $h = 0$ y $\beta J = 1.5$, resuelve $m = \tanh(1.5m)$ iterando
  desde $m = 1$. Compara con la aproximación $\sqrt{3(\beta J - 1)}$.
- **F2.** Calcula el multiplicador social con $\beta = 1$, $J = 0.8$ y
  $h = 0$. ¿Cuánto más responde la adopción a un subsidio que sin
  interacciones?
- **F3.** Explica por qué en la red cuadrada el campo medio predice
  $T_c = 4J$ y la solución exacta de la sesión 13 es menor.

**Ruta avanzada**

- **F4.** Encuentra $H$ de la proposición de la sección 2.2 como función de
  $\beta J$.
- **F5.** Para el modelo de la sección 4, encuentra el valor de $h$ donde
  ocurre el salto cuando $R < R_c$, al menos numéricamente para $J = 1$ y
  $R = 0.6$, y compáralo con el cuaderno.
- **F6.** Desarrolla la energía libre $f(m)$ de la sección 1.2 en serie de
  Taylor alrededor de $m = 0$ hasta orden $m^4$ (para $h = 0$) y verifica
  que tiene la forma de Landau ([F.A1](repaso/fisica_estadistica.md#fa1)).

---

## Para profundizar

Material opcional de la ruta avanzada: demostraciones y temas que no se ven
en clase. Cada apartado indica de qué sección de las notas viene.

<a id="pp1"></a>

### P.1 Derivación con la energía libre (sección 1.2)


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
libre** por agente $f$ ([F.A1](repaso/fisica_estadistica.md#fa1)). Como
$\mathcal S'(m) = -\operatorname{atanh}(m)$,

$$
f'(m) = -Jm - h + T\operatorname{atanh}(m) = 0
\quad\Longleftrightarrow\quad
m = \tanh\left(\beta(h + Jm)\right) .
$$

Además, $f''(m) = -J + T/(1 - m^2)$, así que una solución es un **mínimo**
(equilibrio estable) si $\beta J(1 - m^2) < 1$, y un máximo si
$\beta J(1 - m^2) > 1$. Es la condición gráfica: la curva
$\tanh(\beta(h + Jm))$ cruza la recta $y = m$ con pendiente menor que 1.

<a id="pp2"></a>

### P.2 Campo aleatorio y avalanchas (sección 4)


**Modelo.** A temperatura cero, el agente $i$ adopta si
$h + h_i + Jm > 0$, con $h_i$ independientes con distribución normal
$\mathcal N(0, R^2)$ ([P.B4](repaso/probabilidad.md#pb4),
[F.A3](repaso/fisica_estadistica.md#fa3)). En campo medio y con $N$ grande, la fracción que
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

## Referencias

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
- Manski, C. F. (1993). Identification of endogenous social effects: The
  reflection problem. *Review of Economic Studies*, 60(3), 531–542.
  https://doi.org/10.2307/2298123
- Requião da Cunha, C. (2021). *Introduction to Econophysics: Contemporary
  Approaches with Python Simulations*. CRC Press.
  https://doi.org/10.1201/9781003127956
- Sethna, J. P., Dahmen, K., Kartha, S., Krumhansl, J. A., Roberts, B. W.,
  & Shore, J. D. (1993). Hysteresis and hierarchies: Dynamics of
  disorder-driven first-order phase transformations. *Physical Review
  Letters*, 70(21), 3347–3350. https://doi.org/10.1103/PhysRevLett.70.3347
