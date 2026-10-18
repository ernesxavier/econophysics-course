# Sesión 4: Caminata aleatoria y cadenas de Markov. Notas formales

> **Cómo leer estas notas.** Siguen el orden del programa: primero se observa
> en el cuaderno (`Code/s04_caminata_markov.jl`), después se formaliza aquí.
> El cuerpo de las notas es la **ruta básica**. Las demostraciones y los temas
> de la ruta avanzada, opcionales, están al final, en **Para profundizar**;
> donde estaban queda una línea que los enlaza.
>
> Repaso recomendado: [E.B2–B3](repaso/procesos_estocasticos.md#eb2)
> (caminata aleatoria y promediar simulaciones),
> [P.B2](repaso/probabilidad.md#pb2) (varianza de una suma),
> [E.I3](repaso/procesos_estocasticos.md#ei3) (cadenas de Markov),
> [L.I1](repaso/algebra_lineal.md#li1) (eigenvalores),
> [S.I4](repaso/estadistica.md#si4) (datos dependientes); para la ruta
> avanzada, [E.A1](repaso/procesos_estocasticos.md#ea1) (movimiento
> browniano), [E.A3](repaso/procesos_estocasticos.md#ea3) (ecuación maestra),
> [C.A3](repaso/calculo.md#ca3) (difusión) y
> [L.A5](repaso/algebra_lineal.md#la5) (Perron-Frobenius).
>
> Lecturas: Slanina (2014), sec. 2.1. Complementaria: Requião da Cunha
> (2021), sec. 3.1.

---

## 1. La caminata aleatoria

**Definición.** Sean $\xi_1, \xi_2, \dots$ i.i.d. con $\mathbb{P}(\xi = 1) = p$ y
$\mathbb{P}(\xi = -1) = 1 - p$. La **caminata aleatoria** es
$S_n = \xi_1 + \dots + \xi_n$, con $S_0 = 0$.

**Proposición.**

$$
\mathbb{E}[S_n] = n(2p - 1), \qquad \mathrm{Var}[S_n] = 4np(1 - p).
$$

*Demostración.* $\mathbb{E}[\xi] = p - (1 - p) = 2p - 1$ y
$\mathbb{E}[\xi^2] = 1$, así que $\mathrm{Var}[\xi] = 1 - (2p-1)^2 = 4p(1-p)$.
La media de una suma es la suma de las medias y, por independencia, lo mismo
vale para la varianza ([P.B2](repaso/probabilidad.md#pb2)). $\square$

**Difusión.** La deriva crece como $n$ pero la dispersión como $\sqrt{n}$.
Con $p = 1/2$ no hay deriva y $S_n/\sqrt{n} \to \mathcal{N}(0, 1)$ por el TLC.
Duplicar el horizonte multiplica la dispersión por $\sqrt{2}$, no por 2. Es el
modelo de Bachelier (1900) para los precios y el origen de la regla
$\sigma\sqrt{\Delta t}$ de la sesión 1.

**[Avanzado]** Límite continuo: en [Para profundizar, P.1](#pp1).

---

## 2. Cadenas de Markov

### 2.1 Definición y distribución estacionaria

**Definición.** Un proceso $\lbrace X_t\rbrace$ con estados $\lbrace 1, \dots, K\rbrace$ es una
**cadena de Markov** si
$\mathbb{P}(X_{t+1} = j \mid X_t = i, X_{t-1}, \dots) = P_{ij}$. La matriz $P$
tiene entradas no negativas y cada fila suma 1.

Si $\pi_t$ es el vector fila con las probabilidades de cada estado en el
tiempo $t$, entonces $\pi_{t+1} = \pi_t P$ y $\pi_t = \pi_0 P^t$. Una
**distribución estacionaria** cumple $\pi = \pi P$: es un eigenvector
izquierdo de $P$ con eigenvalor 1 ([L.I1](repaso/algebra_lineal.md#li1)). Si
la cadena es irreducible y aperiódica, es única y $\pi_t \to \pi$ desde
cualquier $\pi_0$ ([L.A5](repaso/algebra_lineal.md#la5)).

### 2.2 El caso de dos estados

Con estados calma (C) y turbulencia (T),

$$
P = \begin{pmatrix} 1 - p & p \cr q & 1 - q \end{pmatrix}, \qquad 0 < p, q < 1.
$$

**Proposición.**

1. $\pi = \left(\frac{q}{p+q},\ \frac{p}{p+q}\right)$.
2. Con $\lambda = 1 - p - q$, $\pi_t - \pi = \lambda^t(\pi_0 - \pi)$: la
   convergencia es geométrica.
3. La duración de un episodio en T es geométrica con media $1/q$.

*Demostración.* (1) La ecuación $\pi = \pi P$ en la primera coordenada da
$\pi_C = \pi_C(1-p) + \pi_T q$, es decir, $\pi_C p = \pi_T q$; junto con
$\pi_C + \pi_T = 1$ resulta (1). (2) Los eigenvalores de $P$ son 1 y
$\lambda$ (su traza es $2 - p - q$ y uno de ellos es 1). Como
$\pi_0 - \pi$ suma cero, es múltiplo del eigenvector izquierdo de $\lambda$,
que es $(1, -1)$, así que $(\pi_0 - \pi)P^t = \lambda^t(\pi_0 - \pi)$.
(3) Cada día en T se sale con probabilidad $q$, independientemente del
pasado. $\square$

### 2.3 Un mercado con dos regímenes de volatilidad

Sea $r_t = \sigma_{X_t}\varepsilon_t$, con $\varepsilon_t$ i.i.d.
$\mathcal{N}(0, 1)$ independientes de la cadena, y volatilidades
$\sigma_C < \sigma_T$.

**Proposición.** En el estado estacionario,

1. $\mathrm{Cov}(r_t, r_{t+\tau}) = 0$ para $\tau \ge 1$;
2. $\mathrm{Cov}(\lvert r_t\rvert, \lvert r_{t+\tau}\rvert) = \frac{2}{\pi}(\sigma_T - \sigma_C)^2\thinspace\pi_C\pi_T\thinspace\lambda^\tau$.

*Demostración.* (1) Los $\varepsilon$ son independientes y de media cero.
(2) Sea $I_t$ el indicador de turbulencia, de modo que
$\sigma_{X_t} = \sigma_C + (\sigma_T - \sigma_C) I_t$. Como
$\mathbb{E}\lvert\varepsilon\rvert = \sqrt{2/\pi}$ y los $\varepsilon$ son
independientes de la cadena,
$\mathrm{Cov}(\lvert r_t\rvert, \lvert r_{t+\tau}\rvert) = \frac{2}{\pi}(\sigma_T - \sigma_C)^2\thinspace\mathrm{Cov}(I_t, I_{t+\tau})$.
Por la proposición anterior,
$\mathbb{P}(I_{t+\tau} = 1 \mid I_t = 1) = \pi_T + \lambda^\tau \pi_C$, de donde
$\mathrm{Cov}(I_t, I_{t+\tau}) = \pi_T(\pi_T + \lambda^\tau\pi_C) - \pi_T^2 = \pi_C\pi_T\lambda^\tau$. $\square$

El modelo reproduce el agrupamiento de volatilidad y, por ser una mezcla de
normales, colas más pesadas que la normal ([sesión 2](s02_volatilidad.md)). Pero la
autocorrelación de $\lvert r\rvert$ decae como $\lambda^\tau$, exponencialmente,
mientras que en los datos decae como una ley de potencias.

---

## 3. Promediar simulaciones: el error de una corrida

Sea $I_t$ el indicador de turbulencia y $\bar I_n = \frac1n\sum_{t=1}^n I_t$ la
fracción de días turbulentos medida en una corrida de $n$ días, empezando en
el estado estacionario. Su media es $\pi_T$, pero su varianza no es la de $n$
observaciones independientes ([E.B3](repaso/procesos_estocasticos.md#eb3),
[S.I4](repaso/estadistica.md#si4)).

**Proposición.** Para $n$ grande,

$$
\mathrm{Var}[\bar I_n] \approx \frac{\pi_C\pi_T}{n}\thinspace\frac{1 + \lambda}{1 - \lambda},
\qquad \lambda = 1 - p - q .
$$

*Demostración.* Por la proposición de 2.3, $\mathrm{Cov}(I_s, I_t) = \pi_C\pi_T\lambda^{|t-s|}$.
Entonces

$$
\mathrm{Var}[\bar I_n] = \frac{1}{n^2}\sum_{s,t}\mathrm{Cov}(I_s, I_t)
\approx \frac{\pi_C\pi_T}{n}\Big(1 + 2\sum_{\tau \ge 1}\lambda^\tau\Big)
= \frac{\pi_C\pi_T}{n}\thinspace\frac{1+\lambda}{1-\lambda} . \quad\square
$$

El factor $(1+\lambda)/(1-\lambda)$ dice cuántas veces "menos datos" tiene una
corrida persistente: $n$ días equivalen a
$n(1-\lambda)/(1+\lambda)$ días independientes (el **tamaño efectivo de
muestra**). Con $p = 0.01$ y $q = 0.05$, $\lambda = 0.94$ y el factor vale
$1.94/0.06 \approx 32$: 3000 días valen como unos 93. Promediar $M$ corridas
independientes divide la varianza entre $M$.

---

## 4. [Avanzado] De la ecuación maestra a la difusión

El desarrollo de esta sección está en [Para profundizar, P.2](#pp2).


---

## 5. Resumen formal

| Proceso | Ecuación | Resultado clave |
|---|---|---|
| Caminata aleatoria | $S_n = \sum \xi_i$ | $\mathrm{Var}[S_n] \propto n$: difusión |
| Cadena de Markov (2 estados) | $\pi_{t+1} = \pi_t P$ | $\pi_T = p/(p+q)$; convergencia como $(1-p-q)^t$ |
| Dos regímenes de volatilidad | $r_t = \sigma_{X_t}\varepsilon_t$ | ACF de $\lvert r\rvert \propto (1-p-q)^\tau$ |
| Una corrida de $n$ días | $\bar I_n$ | Varianza inflada por $(1+\lambda)/(1-\lambda)$ |

---

## 6. Ejercicios formales

**Ruta básica**

- **F1.** En una caminata aleatoria con $p = 0.52$, ¿después de cuántos pasos
  la deriva esperada iguala una desviación estándar? Interprétalo: ¿cuánto
  tiempo hace falta para distinguir una ventaja pequeña del azar?
- **F2.** Para la cadena de dos estados con $p = 0.01$ y $q = 0.05$, calcula
  $\pi$, la duración media de la turbulencia y cuántos días tarda
  $\lvert\pi_t - \pi\rvert$ en reducirse a la mitad.
- **F3.** Con la misma cadena, ¿qué desviación estándar tiene la fracción de
  días turbulentos medida en una corrida de 3000 días? ¿Cuántas corridas hay
  que promediar para reducirla a la mitad?

**Ruta avanzada**

- **F4.** Muestra que el modelo de dos regímenes con $\sigma_C = \sigma$ y
  $\sigma_T = 3\sigma$ tiene curtosis $3(\pi_C + 81\pi_T)/(\pi_C + 9\pi_T)^2$.
  ¿Para qué $\pi_T$ es máxima?
- **F5.** Demuestra que la ecuación de difusión de la sección 4 tiene como
  solución, desde un punto, una normal con varianza $2Dt$.

---

## Para profundizar

Material opcional de la ruta avanzada: demostraciones y temas que no se ven
en clase. Cada apartado indica de qué sección de las notas viene.

<a id="pp1"></a>

### P.1 Límite continuo (sección 1)

**Límite continuo.** Si los pasos ocurren cada $\Delta t$ y miden
$\sqrt{\Delta t}$, el proceso $W^{\Delta t}(t) = \sqrt{\Delta t}\thinspace S_{\lfloor t/\Delta t\rfloor}$
converge, cuando $\Delta t \to 0$, al **movimiento browniano** $W_t$ (teorema de
Donsker): un proceso continuo con incrementos independientes y
$W_t - W_s \sim \mathcal{N}(0, t - s)$ ([E.A1](repaso/procesos_estocasticos.md#ea1)).

<a id="pp2"></a>

### P.2 De la ecuación maestra a la difusión (sección 4)


Una caminata con pasos $\pm\delta$ cada $\tau$ unidades de tiempo y
probabilidad $1/2$ cumple la **ecuación maestra**

$$
P(x, t + \tau) = \tfrac{1}{2}P(x - \delta, t) + \tfrac{1}{2}P(x + \delta, t).
$$

Desarrollando por Taylor a primer orden en $\tau$ y a segundo en $\delta$,
$\tau\thinspace\partial_t P = \frac{\delta^2}{2}\partial_x^2 P$, es decir,

$$
\frac{\partial P}{\partial t} = D\frac{\partial^2 P}{\partial x^2}, \qquad D = \frac{\delta^2}{2\tau},
$$

la ecuación de difusión ([C.A3](repaso/calculo.md#ca3)), cuya solución desde
un punto es una normal con varianza $2Dt$.

---

## Referencias

- Bachelier, L. (1900). Théorie de la spéculation. *Annales scientifiques de
  l'École normale supérieure*, 17, 21–86. https://doi.org/10.24033/asens.476
- Requião da Cunha, C. (2021). *Introduction to Econophysics: Contemporary
  Approaches with Python Simulations*. CRC Press.
  https://doi.org/10.1201/9781003127956
- Slanina, F. (2014). *Essentials of Econophysics Modelling*. Oxford
  University Press.
