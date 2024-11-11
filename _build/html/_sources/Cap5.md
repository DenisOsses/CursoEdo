---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

<!--######################################################################################################################################################################################################################################################################################################################################################
-->

+++

(U5)=
# Unidad 5

+++

## Soluciones de EDO mediante Series de Potencias

+++

### Introducción. Sistemas Masa-Resorte con Resorte Variable

Podemos suponer que en un sistema masa-resorte modelado como [](eqMAS), a saber $mx''+kx=0$, que está en movimiento durante un largo tiempo, el resorte se debilite; es decir, varía la *constante del resorte* y ésta decae con el tiempo. En un modelo para el **resorte cada vez más viejo** la constante de resorte $k$ se reemplaza con la función decreciente $K(t)=ke^{-\alpha t}$, $k,\alpha>0$, obteniendo la ecuación diferencial lineal 

$$
mx''+ke^{-\alpha t}x=0
$$ (ResorteViejo)

También, cuando un sistema masa-resorte se somete a un ambiente en el que la temperatura disminuye con rapidez, podría tener sentido reemplazar la constante $k$ con $K(t)=kt$, $k>0$, una función que se incrementa con el tiempo ya que es razonable pensar que el frío aumenta la resistencia o endurece de los resortes. El modelo resultante se conoce como **ecuación diferencial de Airy**

$$
mx''+ktx=0
$$ (ResorteAiry)

¿Es posible resolver estas ecuaciones diferenciales? No podemos hacerlo con los métodos vistos hasta ahora.

+++

### Soluciones Respecto a Puntos Ordinarios

Consideremos la EDO lineal de segundo orden homogénea 

$$
a_2(x)y'''+a_1(x)y'+a_0(x)y=0
$$ (Serie1)

que en forma estándar es 

$$
y''+P(x)y'+Q(x)y=0
$$ (Serie2)

**Definición**: Decimos que $x_0$ es un **punto ordinario** de [](Serie1) si tanto $P(x)$ como $Q(x)$ en la forma estándar [](Serie2) son **analíticas** en $x_0$, es decir $P(x)$ y $Q(x)$ se pueden desarrollar en serie de potencias en torno a $x_0$:

$$
P(x)=\sum_{n=0}^\infty a_n(x-x_0)^n~~~~,~~~~Q(x)=\sum_{n=0}^\infty b_n(x-x_0)^n
$$ 

Un punto que no es ordinario se denomina **punto singular** de la Ecuación Diferencial.

**Teorema**: Si $x_0$ es un punto ordinario de [](Serie1) entonces siempre es posible encontrar 2 soluciones L.I. en forma de series de potencias en torno a $x_0$. Una solución en serie converge, por lo menos, en un intervalo $]x_0-R,x_0+R[$ donde $R$ es la distancia de $x_0$ al punto singular más cercano.

```{admonition} Ejercicio Teórico
Determine los puntos ordinarios y singulares de 

$$
(x^2-x)y''-3xy'+y=0.
$$
```

**Nota**: Habitualmente analizaremos ecuaciones con $x_0=0$ como punto ordinario.

### Método de Solución Mediante Series 

Si $x_0=0$ es un punto ordinario de la ecuación [](Serie1), reemplazamos la solución 

$$
y=\sum_{n=0}^\infty c_nx^n
$$ 

en [](Serie1). Esto nos lleva a una **ecuación de recurrencia** para los coeficientes $c_n$ de la serie de potencias solución.

```{admonition} Ejercicio Teórico
Resuelva la ecuación diferencial de Airy (con $k=1$): $y''+xy=0$ en torno al punto ordinario $x_0=0$. 
```

