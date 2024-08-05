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

(Glosario)=
# Glosario

+++

(EDO)= 
<u>**Ecuación Diferencial Ordinaria (EDO)**</u>: Una Ecuación Diferencial es una ecuación que tiene como incógnita a una función (usualmente $y$) de una variable independiente (habitualmente $x$ o $t$) junto a sus derivadas. Se puede representar de modo general como 

$$
F(x,y,y',y'',\ldots,y^{(n)})=0
$$

(PrimerOrden)=
<u>**EDO de Primer Orden**</u>: Una EDO se llama de Primer Orden si sólo contiene hasta la primera derivada de una función, digamos $F(x,y,y')=0$.

(SegundoOrden)=
<u>**EDO de Segundo Orden**</u>: Una EDO se llama de Segundo Orden si contiene hasta la segunda derivada de una función, digamos $F(x,y,y',y'')=0$.

(DefPVI)=
<u>**Problema de Valor Inicial (PVI)**</u>: Un PVI consiste de una EDO $F(x,y,y',y'',\ldots,y^{(n)})=0$ junto a las $n$ condiciones iniciales $y^{(i)}(x_0)=y_i$, para $i=0,1,2,\ldots,n-1$.

(SepVar)=
<u>**Método de separación de variables**</u>: Una EDO está en variables separables si es de la forma 

$$
\frac{dy}{dx}=g(x)h(y).
$$ 

Se resuelve dejando todas las variables de un tipo a un lado de la ED0 e integrando: 

$$
\int\frac{1}{h(y)}dy=\int g(x)~dx.
$$

(Explicita)=
<u>**Solución explícita**</u>: Es la solución de un PVI que se obtiene mediante un método analítico y tiene la forma $y=f(x)$, para $x\in~I$, donde $I$ es un intervalo.

(Intervalo)=
<u>**Intervalo de definición**</u>: Es el conjunto $I$ donde la solución de una EDO existe.

(Lineal1er)=
<u>**Lineal de Primer Orden**</u>: Una EDO es lineal de primer orden si tiene la forma 

$$
a_1(x)\frac{dy}{dx}+a_0(x)y=g(x).
$$ 

Para resolver esta EDO, la escribimos en su **forma estándar** 

$$
\frac{dy}{dx}+P(x)y=f(x),
$$ 

donde $P(x)$ y $f(x)$ son funciones continuas en un intervalo $I$. Luego identificamos la función $P(x)$ y calculamos el **factor integrante** 

$$
e^{\int P(x)~dx}
$$

Posteriormente, multiplicamos la EDO (en su forma estándar) por este factor integrante. Notamos que obtenemos la derivada de un producto: 

$$
\frac{d}{dx}\left(e^{\int P(x)~dx}\cdot y\right)=e^{\int P(x)~dx}\cdot f(x).
$$ 

Finalmente, integramos esta igualdad para obtener la solución $y(x)$ de la EDO:

$$y(x)=\dfrac{1}{e^{\int P(x)~dx}}\left(\int e^{\int P(x)~dx}\cdot f(x)~dx\right)+C$$
