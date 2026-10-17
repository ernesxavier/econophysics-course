# Sesión 3: Colas pesadas y leyes de potencias. Notas formales

> **Cómo leer estas notas.** Siguen el orden del programa: primero se observa
> en el cuaderno (`Code/s03_colas_pesadas.jl`), después se formaliza aquí. Lo
> que no tiene marca es **ruta básica**. Las secciones y demostraciones
> marcadas **[Avanzado]** son opcionales.
>
> Repaso recomendado: [P.I5](repaso/probabilidad.md#pi5) (leyes de potencias),
> [C.B2](repaso/calculo.md#cb2) (escalas logarítmicas),
> [S.B4](repaso/estadistica.md#sb4) (mínimos cuadrados); para la ruta avanzada,
> [P.A2](repaso/probabilidad.md#pa2), [P.A3](repaso/probabilidad.md#pa3) y
> [P.A4](repaso/probabilidad.md#pa4).
>
> Lecturas: Gabaix (2009); Mantegna & Stanley, sec. 4.4; Requião da Cunha,
> sec. 2.2.7; Thurner, Klimek & Hanel, cap. 3. Ruta avanzada: Clauset, Shalizi
> & Newman (2009); Mantegna & Stanley, caps. 3–4.

---

## 1. Leyes de potencias

### 1.1 Definición

**Definición (cola de ley de potencias).** Una variable aleatoria positiva $X$
tiene una cola de ley de potencias con exponente $\alpha > 0$ si

$$
\bar F(x) = \mathbb{P}(X > x) \sim C\thinspace x^{-\alpha} \qquad (x \to \infty).
$$

El caso exacto es la **distribución de Pareto** con parámetros $\alpha$ y
$x_{\min} > 0$:

$$
\bar F(x) = \left(\frac{x}{x_{\min}}\right)^{-\alpha},
\qquad
f(x) = \frac{\alpha\thinspace x_{\min}^{\alpha}}{x^{\alpha+1}},
\qquad x \ge x_{\min}.
$$

Dos convenciones conviven en la literatura: el exponente de la CCDF ($\alpha$)
y el de la densidad ($\alpha + 1$). En el curso usamos siempre el de la CCDF.

**En escala log-log.** De la definición,
$\ln \bar F(x) = \alpha \ln x_{\min} - \alpha \ln x$: una recta de pendiente
$-\alpha$. Una exponencial, $\ln \bar F(x) = -\lambda x$, es recta en escala
semilog y en log-log se curva hacia abajo.

**[Avanzado] La versión general.** En rigor, se pide
$\bar F(x) = L(x)\thinspace x^{-\alpha}$ con $L$ de **variación lenta**:
$L(tx)/L(x) \to 1$ cuando $x \to \infty$, para todo $t > 0$ (por ejemplo,
una constante o $\ln x$). Se dice que $\bar F$ es de **variación regular**
con índice $-\alpha$. Todas las propiedades de esta sección valen en esa
generalidad.

### 1.2 Invariancia de escala

**Proposición.** Para la Pareto, $f(\lambda x) = \lambda^{-(\alpha+1)} f(x)$ para
todo $\lambda > 0$ y $x \ge x_{\min}$ (con $\lambda x \ge x_{\min}$).

*Demostración.* Sustituir en la fórmula de la densidad. $\square$

Cambiar la unidad de medida sólo multiplica la densidad por una constante: la
forma de la distribución es la misma a cualquier escala. Por eso **no hay un
tamaño típico**. Una exponencial, en cambio, tiene la escala característica
$1/\lambda$, y una normal, la escala $\sigma$.

### 1.3 Momentos

**Proposición.** Si $X$ es Pareto con exponente $\alpha$, entonces
$\mathbb{E}[X^n] < \infty$ si y sólo si $n < \alpha$, y en ese caso

$$
\mathbb{E}[X^n] = \frac{\alpha}{\alpha - n}\thinspace x_{\min}^{n}.
$$

*Demostración.*
$\mathbb{E}[X^n] = \alpha x_{\min}^{\alpha}\int_{x_{\min}}^{\infty} x^{n-\alpha-1}\thinspace dx$.
La integral de $x^{-a}$ en $[x_{\min}, \infty)$ converge si y sólo si $a > 1$
([C.I2](repaso/calculo.md#ci2)); aquí $a = \alpha + 1 - n$, así que converge
si y sólo si $n < \alpha$, y vale $x_{\min}^{n-\alpha}/(\alpha - n)$. $\square$

| Exponente | Media | Varianza | Curtosis | Ejemplo |
|---|---|---|---|---|
| $\alpha \le 1$ | infinita | infinita | infinita | tamaño de ciudades |
| $1 < \alpha \le 2$ | finita | infinita | infinita | riqueza de los más ricos |
| $2 < \alpha \le 4$ | finita | finita | infinita | retornos financieros |
| $\alpha > 4$ | finita | finita | finita | |

Con un momento infinito, su estimador muestral no converge al crecer $N$:
crece o salta cada vez que aparece un dato extremo. Es lo que vimos en la
sesión 1 con la curtosis de los retornos.

### 1.4 Concentración

**Proposición.** Sea $X$ Pareto con $\alpha > 1$. La fracción del total que
concentra la fracción $p$ de observaciones más grandes es

$$
S(p) = p^{(\alpha - 1)/\alpha}.
$$

*Demostración.* El umbral $x_p$ que deja por encima una fracción $p$ cumple
$(x_p/x_{\min})^{-\alpha} = p$, es decir, $x_p = x_{\min}\thinspace p^{-1/\alpha}$. Entonces

$$
S(p) = \frac{\int_{x_p}^{\infty} x f(x)\thinspace dx}{\mathbb{E}[X]}
= \frac{\frac{\alpha}{\alpha-1} x_{\min}^{\alpha} x_p^{1-\alpha}}{\frac{\alpha}{\alpha-1} x_{\min}}
= \left(\frac{x_p}{x_{\min}}\right)^{1-\alpha}
= p^{(\alpha-1)/\alpha}. \quad\square
$$

**Ejemplo.** La "regla 80/20" (el 20% concentra el 80%) corresponde a
$0.2^{(\alpha-1)/\alpha} = 0.8$, es decir, $\alpha \approx 1.16$. Cuando
$\alpha \to 1$, $S(p) \to 1$ para cualquier $p$: casi todo el total está en
unas pocas observaciones. En las localidades de México, el 1% más grande
concentra el 62% de la población de las localidades con al menos 100
habitantes.

---

## 2. La ley de Zipf

### 2.1 Rango y CCDF

Con $N$ observaciones ordenadas de mayor a menor, $x_{(1)} \ge x_{(2)} \ge \dots \ge x_{(N)}$,
el **rango** de $x_{(R)}$ es $R$. Como exactamente $R$ observaciones son
mayores o iguales que $x_{(R)}$,

$$
R = N\thinspace \hat F_{\ge}(x_{(R)}),
$$

donde $\hat F_{\ge}(x)$ es la fracción de datos mayores o iguales que $x$
(la CCDF empírica). La gráfica de $\ln R$ contra $\ln x$ es la CCDF empírica
en log-log desplazada por $\ln N$.

**Proposición.** Si la cola es Pareto con exponente $\alpha$, la observación de
rango $R$ cumple aproximadamente

$$
x_{(R)} \approx x_{\min}\left(\frac{R}{N_{\min}}\right)^{-1/\alpha},
$$

donde $N_{\min}$ es el número de observaciones mayores que $x_{\min}$.

*Demostración.* $R \approx N_{\min}\thinspace (x_{(R)}/x_{\min})^{-\alpha}$; despejar
$x_{(R)}$. $\square$

**Definición (ley de Zipf).** Es el caso $\alpha = 1$: $x_{(R)} \propto 1/R$.
La segunda observación es la mitad de la primera, la tercera un tercio, y
el producto rango × tamaño es aproximadamente constante.

### 2.2 Medición en México

Con las localidades del Censo 2020 (INEGI), el método de la sección 4.4 elige
$x_{\min} \approx 665$ habitantes y da $\hat\alpha = 0.99 \pm 0.02$ (intervalo
de 95%), con 15 829 localidades en la cola: la ley de Zipf casi exacta.

Dos advertencias. Primero, una localidad no es una ciudad: la Ciudad de México
está dividida en alcaldías y las grandes zonas metropolitanas en varios
municipios, lo que deforma la parte alta. Segundo, la regularidad empírica no
explica su origen; la sesión 6 estudia un mecanismo (crecimiento
multiplicativo con una barrera) que la genera.

---

## 3. Ley de potencias contra lognormal

**Definición.** $X$ es **lognormal** con parámetros $\mu, \sigma$ si
$\ln X \sim \mathcal{N}(\mu, \sigma^2)$. Todos sus momentos son finitos:
$\mathbb{E}[X^n] = e^{n\mu + n^2\sigma^2/2}$. Por lo tanto **no** tiene cola de
ley de potencias.

Sin embargo, en una gráfica log-log puede parecerlo. Definimos el
**exponente local** como menos la pendiente de la CCDF en log-log:

$$
\alpha_{\text{loc}}(x) = -\frac{d \ln \bar F(x)}{d \ln x}.
$$

Para una ley de potencias exacta, $\alpha_{\text{loc}}(x) = \alpha$ en todo $x$.

**Proposición.** Para la lognormal, con $u = (\ln x - \mu)/\sigma$,

$$
\alpha_{\text{loc}}(x) = \frac{1}{\sigma}\thinspace\frac{\varphi(u)}{1 - \Phi(u)}
\approx \frac{\ln x - \mu}{\sigma^2} \quad \text{para } u \text{ grande},
$$

donde $\varphi$ y $\Phi$ son la densidad y la CDF de la normal estándar.

*Demostración.* $\bar F(x) = 1 - \Phi(u)$ y $du/d\ln x = 1/\sigma$. Por la
regla de la cadena,
$d\ln\bar F/d\ln x = -\varphi(u)/\big(\sigma(1 - \Phi(u))\big)$. La
aproximación usa que $\varphi(u)/(1 - \Phi(u)) \approx u$ cuando
$u \to \infty$ (la razón de Mills). $\square$

**Consecuencia.** El exponente local crece con $x$, pero despacio: al avanzar
un orden de magnitud ($\Delta\ln x = \ln 10 \approx 2.3$), cambia en
$2.3/\sigma^2$. Con $\sigma = 2$ cambia apenas 0.58 por década, y en dos o
tres décadas los datos parecen una recta. Con $\sigma = 0.5$ cambia más de 9
por década: la curvatura es evidente.

Por eso, afirmar que hay una ley de potencias exige (Clauset et al., 2009):
varios órdenes de magnitud en la cola, un estimador correcto de $\alpha$
(sección 4) y una comparación formal con alternativas (sección 4.5).

---

## 4. Estimación del exponente

### 4.1 Máxima verosimilitud

**Proposición.** Sea $x_1, \dots, x_n$ una muestra i.i.d. de una Pareto con
$x_{\min}$ conocido. El estimador de máxima verosimilitud de $\alpha$ es

$$
\hat\alpha = \frac{n}{\sum_{i=1}^{n} \ln(x_i/x_{\min})},
\qquad
\mathrm{SE}(\hat\alpha) \approx \frac{\hat\alpha}{\sqrt{n}}.
$$

*Demostración.* La log-verosimilitud es
$\ell(\alpha) = n\ln\alpha + n\alpha\ln x_{\min} - (\alpha+1)\sum_i \ln x_i$.
Igualar $\ell'(\alpha) = n/\alpha - \sum_i \ln(x_i/x_{\min})$ a cero da
$\hat\alpha$. La curvatura es $\ell''(\alpha) = -n/\alpha^2$, de donde el error
estándar ([S.A1](repaso/estadistica.md#sa1)). $\square$

**[Avanzado] Sesgo.** Si $X$ es Pareto, $\ln(X/x_{\min})$ es exponencial con
tasa $\alpha$ ([P.I6](repaso/probabilidad.md#pi6)), así que la suma del
denominador es una gamma con forma $n$ y tasa $\alpha$. De ahí
$\mathbb{E}[\hat\alpha] = \frac{n}{n-1}\thinspace\alpha$: el estimador sobreestima
ligeramente, y $\frac{n-1}{n}\hat\alpha$ es insesgado.

### 4.2 El estimador de Hill

En los datos reales sólo la cola sigue una ley de potencias, y $x_{\min}$ es
desconocido.

**Definición (Hill, 1975).** Con los datos ordenados de mayor a menor,

$$
\hat\alpha_k = \frac{k}{\sum_{i=1}^{k} \ln\left(x_{(i)}/x_{(k+1)}\right)}.
$$

Es el estimador de máxima verosimilitud de la sección 4.1 aplicado a los $k$
datos más grandes, con $x_{\min} = x_{(k+1)}$.

**Elección de $k$.** Hay una disyuntiva entre sesgo y varianza:

- $k$ pequeño: pocos datos, error estándar $\hat\alpha_k/\sqrt{k}$ grande.
- $k$ grande: se incluye el cuerpo de la distribución, que no es ley de
  potencias, y el estimador se sesga (en los retornos, hacia abajo).

**[Avanzado]** Para distribuciones con cola de variación regular,
$\hat\alpha_k$ es consistente si $k \to \infty$ y $k/N \to 0$ cuando
$N \to \infty$. La gráfica de Hill ($\hat\alpha_k$ contra $k$) ayuda a
encontrar una zona estable.

### 4.3 Por qué no ajustar una recta

Estimar $\alpha$ como menos la pendiente de mínimos cuadrados de
$\ln \hat F_{\ge}$ contra $\ln x$ viola los supuestos de la regresión
([S.I5](repaso/estadistica.md#si5)):

1. Los puntos de la CCDF empírica **no son independientes**: cada uno cuenta
   todos los datos a su derecha.
2. Los errores **no tienen varianza constante**: los puntos de la cola, con
   muy pocos datos, son mucho más ruidosos y pesan igual.
3. Los errores estándar de la regresión no son válidos.

En el experimento del cuaderno (1000 muestras de 200 datos de una Pareto con
$\alpha = 2$), la recta dio $\hat\alpha = 1.93$ en promedio con desviación
estándar 0.20, y la máxima verosimilitud, $2.01$ con desviación 0.14.

### 4.4 Elección de $x_{\min}$: el método de Clauset, Shalizi y Newman

**Distancia de Kolmogorov-Smirnov.** Para un candidato $x_{\min}$, sean
$S(x)$ la CDF empírica de los datos mayores o iguales que $x_{\min}$ y
$P(x) = 1 - (x/x_{\min})^{-\hat\alpha}$ la CDF de la Pareto ajustada. Su
distancia es

$$
D(x_{\min}) = \max_{x \ge x_{\min}} \lvert S(x) - P(x)\rvert.
$$

**Receta.** Calcular $\hat\alpha$ y $D$ para cada candidato y elegir el
$x_{\min}$ que minimiza $D$. Con $x_{\min}$ muy bajo, el cuerpo de la
distribución aleja a los datos de la ley de potencias y $D$ es grande; con
$x_{\min}$ muy alto, quedan pocos datos y $D$ vuelve a crecer por ruido.

**[Avanzado] Bondad de ajuste.** El p-valor no puede tomarse de las tablas de
Kolmogorov-Smirnov, porque $\hat\alpha$ y $x_{\min}$ se estimaron con los
mismos datos. Clauset et al. proponen un bootstrap semiparamétrico: generar
muchas muestras sintéticas del mismo tamaño (los datos reales por debajo de
$x_{\min}$, una Pareto ajustada por encima), aplicar a cada una la receta
completa y calcular la fracción con $D$ mayor que la observada. Si es menor
que 0.1, se rechaza la ley de potencias.

### 4.5 [Avanzado] Comparación con alternativas

Que una ley de potencias no se rechace no significa que sea el mejor modelo.
Para compararla con una alternativa no anidada, como la lognormal truncada
en $x_{\min}$, se usa la razón de log-verosimilitudes

$$
\mathcal{R} = \sum_{i=1}^{n}\left[\ln f_{\text{pot}}(x_i) - \ln f_{\text{ln}}(x_i)\right].
$$

Su signo indica qué modelo ajusta mejor y la prueba de Vuong (1989) dice si la
diferencia es significativa: bajo la hipótesis de que ambos ajustan igual,
$\mathcal{R}/(\hat s\sqrt{n})$ es aproximadamente normal estándar, con
$\hat s$ la desviación estándar de los sumandos
([S.A2](repaso/estadistica.md#sa2)). Con datos económicos es frecuente que la
prueba no distinga entre las dos.

---

## 5. La ley cúbica inversa de los retornos

**Hecho empírico** (Gopikrishnan et al., 1998). Para los retornos
estandarizados $z$ de acciones, índices y divisas,

$$
\mathbb{P}(\lvert z\rvert > x) \sim x^{-\alpha}, \qquad \alpha \approx 3.
$$

**Medición en los datos del curso** (método de la sección 4.4; intervalo de
95% para el exponente de $\lvert z\rvert$):

| Activo | $\alpha$ caídas | $\alpha$ subidas | $\alpha$ de $\lvert z\rvert$ | $x_{\min}$ (en $\sigma$) | Días en la cola |
|---|---|---|---|---|---|
| IPC (BMV) | 2.80 | 2.83 | 2.82 ± 0.22 | 1.56 | 638 |
| USD/MXN | 3.48 | 2.83 | 3.10 ± 0.31 | 1.79 | 387 |
| S&P 500 | 3.14 | 2.79 | 3.12 ± 0.47 | 2.59 | 169 |
| Bitcoin | 2.70 | 4.65 | 4.04 ± 0.76 | 2.76 | 110 |

**Consecuencias.**

- Con $\alpha \approx 3$, la varianza es finita y la curtosis infinita
  (sección 1.3).
- **Regla de duplicación.** $\mathbb{P}(\lvert z\rvert > 2x)/\mathbb{P}(\lvert z\rvert > x) \to 2^{-\alpha} = 1/8$.
  Para la normal, el mismo cociente tiende a cero más rápido que cualquier
  potencia: con $x = 3$ ya vale menos de $10^{-5}$.
- Como $\alpha > 2$, el TLC clásico se aplica a los retornos agregados, que
  tienden a la normal, aunque lentamente (sección 6).

---

## 6. [Avanzado] Distribuciones estables y TLC generalizado

**Definición.** Una variable $X$ es **estable** si para $X_1, X_2$ copias
independientes de $X$ y cualesquiera $a, b > 0$ existen $c > 0$ y $d$ tales
que $aX_1 + bX_2$ tiene la misma distribución que $cX + d$.

**Caso simétrico.** Las estables simétricas centradas tienen función
característica

$$
\varphi(s) = \mathbb{E}[e^{isX}] = e^{-\gamma\lvert s\rvert^{\alpha}},
\qquad 0 < \alpha \le 2,\ \gamma > 0.
$$

Con $\alpha = 2$ es la normal; con $\alpha = 1$, la de Cauchy. No hay fórmula
cerrada de la densidad en general.

**Proposición (escalamiento).** Si $X_1, \dots, X_n$ son i.i.d. estables
simétricas de índice $\alpha$, entonces $S_n = X_1 + \dots + X_n$ tiene la
misma distribución que $n^{1/\alpha} X_1$.

*Demostración.* Por independencia,
$\varphi_{S_n}(s) = \varphi(s)^n = e^{-n\gamma\lvert s\rvert^{\alpha}} = e^{-\gamma\lvert n^{1/\alpha}s\rvert^{\alpha}} = \varphi(n^{1/\alpha}s)$,
que es la función característica de $n^{1/\alpha}X_1$. $\square$

**Colas.** Para $\alpha < 2$, $\mathbb{P}(\lvert X\rvert > x) \sim C x^{-\alpha}$:
las estables no gaussianas tienen colas de ley de potencias y varianza
infinita.

**Teorema (TLC generalizado, enunciado).** Sean $X_i$ i.i.d. simétricas con
$\mathbb{P}(\lvert X\rvert > x) \sim C x^{-\alpha}$.

- Si $0 < \alpha < 2$, entonces $S_n/n^{1/\alpha}$ converge en distribución a
  una estable de índice $\alpha$.
- Si $\alpha > 2$, la varianza es finita y $S_n/\sqrt{n}$ converge a una
  normal (TLC clásico).

En el caso $\alpha > 2$, la convergencia ocurre en el centro de la
distribución; la cola de ley de potencias de $S_n$ persiste más allá de una
región que crece con $n$. Es la razón formal de que la normalidad
agregacional de los retornos llegue lentamente. Mandelbrot (1963) propuso
estables con $\alpha \approx 1.7$ para los precios; las mediciones de la
sección 5, con $\alpha \approx 3$, descartan esa hipótesis para los retornos.

---

## 7. Resumen formal

| Concepto | Enunciado | Dónde |
|---|---|---|
| Ley de potencias | $\bar F(x) \sim C x^{-\alpha}$; recta de pendiente $-\alpha$ en log-log | 1.1 |
| Momentos | $\mathbb{E}[X^n] < \infty$ si y sólo si $n < \alpha$ | 1.3 |
| Concentración | $S(p) = p^{(\alpha-1)/\alpha}$ | 1.4 |
| Ley de Zipf | $\alpha = 1$: $x_{(R)} \propto 1/R$; México: $\hat\alpha = 0.99 \pm 0.02$ | 2 |
| Lognormal | Exponente local $\approx (\ln x - \mu)/\sigma^2$: casi constante si $\sigma$ es grande | 3 |
| Estimación | Máxima verosimilitud o Hill, con $x_{\min}$ por KS; no usar la recta | 4 |
| Ley cúbica inversa | $\alpha \approx 3$ en retornos: varianza finita, curtosis infinita | 5 |
| Estables | $S_n \overset{d}{=} n^{1/\alpha}X_1$; atraen a las sumas con $\alpha < 2$ | 6 |

---

## 8. Ejercicios formales

**Ruta básica**

- **F1.** Si la localidad más grande de un país tiene 2 millones de habitantes
  y se cumple la ley de Zipf, ¿cuántos habitantes tienen las localidades de
  rango 10 y 100? ¿Cuánto vale el producto rango × tamaño?
- **F2.** Con la fórmula de concentración, ¿qué fracción de la riqueza tiene
  el 1% más rico si $\alpha = 1.5$? ¿Y si $\alpha = 3$?
- **F3.** En una muestra de retornos con $\alpha = 3$ hay 40 días al año con
  $\lvert z\rvert > 2$. ¿Cuántos días al año esperas con $\lvert z\rvert > 4$ y
  con $\lvert z\rvert > 8$?

**Ruta avanzada**

- **F4.** Demuestra que si $X$ es Pareto con $x_{\min}$, entonces
  $\ln(X/x_{\min})$ es exponencial con tasa $\alpha$, y usa ese resultado para
  probar que $\mathbb{E}[\hat\alpha] = \frac{n}{n-1}\alpha$.
- **F5.** Calcula el exponente local de una CCDF de la forma
  $\bar F(x) = C x^{-\alpha} e^{-x/x_c}$ (ley de potencias con corte
  exponencial). ¿Cómo se ve en log-log para $x \ll x_c$ y para $x \gg x_c$?
- **F6.** Demuestra la proposición de escalamiento de la sección 6 para
  sumas ponderadas: si $X_1, X_2$ son estables simétricas i.i.d. de índice
  $\alpha$, entonces $aX_1 + bX_2 \overset{d}{=} (a^\alpha + b^\alpha)^{1/\alpha} X_1$.
- **F7.** Muestra que la Pareto con $\alpha > 2$ y la lognormal con el mismo
  promedio y la misma varianza pueden tener valores de curtosis muy
  distintos. ¿Cuál puede tener curtosis infinita?

---

## Referencias

- Clauset, A., Shalizi, C. R., & Newman, M. E. J. (2009). Power-law
  distributions in empirical data. *SIAM Review*, 51(4), 661–703.
  https://doi.org/10.1137/070710111
- Gabaix, X. (2009). Power laws in economics and finance. *Annual Review of
  Economics*, 1, 255–294. https://doi.org/10.1146/annurev.economics.050708.142940
- Gopikrishnan, P., Meyer, M., Amaral, L. A. N., & Stanley, H. E. (1998).
  Inverse cubic law for the distribution of stock price variations.
  *The European Physical Journal B*, 3(2), 139–140.
  https://doi.org/10.1007/s100510050292
- Hill, B. M. (1975). A simple general approach to inference about the tail
  of a distribution. *The Annals of Statistics*, 3(5), 1163–1174.
  https://doi.org/10.1214/aos/1176343247
- INEGI (2021). *Censo de Población y Vivienda 2020: Principales resultados
  por localidad (ITER)*. Instituto Nacional de Estadística y Geografía.
- Mandelbrot, B. (1963). The variation of certain speculative prices. *The
  Journal of Business*, 36(4), 394–419. https://doi.org/10.1086/294632
- Mantegna, R. N., & Stanley, H. E. (1999). *An Introduction to Econophysics:
  Correlations and Complexity in Finance*. Cambridge University Press.
  https://doi.org/10.1017/CBO9780511755767
- Newman, M. E. J. (2005). Power laws, Pareto distributions and Zipf's law.
  *Contemporary Physics*, 46(5), 323–351.
  https://doi.org/10.1080/00107510500052444
- Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
  non-nested hypotheses. *Econometrica*, 57(2), 307–333.
  https://doi.org/10.2307/1912557
