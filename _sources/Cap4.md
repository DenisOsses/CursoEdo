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

(U4)=
# Unidad 4

+++

## Transformada de Laplace

En este capítulo estudiaremos una herramienta muy útil para resolver EDO Lineales: La [**transformada integral**](https://en.wikipedia.org/wiki/Integral_transform) de Laplace, que mapea o transforma una función $f(t)$ en otra función $F(s)$ mediante una integración:

$$
T\left(f(t)\right)=\int_{t_1}^{t_2}K(s,t)f(t)~dt=F(s)
$$ 

donde $K$ es conocida como el **núcleo** o **kernel** de la transformación y $t_1,t_2$ son valores que dependen de la definición particular que se utilice, pudiendo ser $\pm\infty$.

La esencia de toda transformada integral es *transformar* un problema complicado -que involucre a $f(t)$- en uno nuevo para $F(s)$ que simplifique el original, para luego retornar a éste mediante una *transformada inversa* y resolverlo.

+++

### Definición y Existencia

**Definición**: Si $f(t)$ es una función definida $\forall~t>0$, entonces se dice que la integral 

$$
\mathscr{L}\{f(t)\}=\int_0^\infty e^{-st}f(t)~dt=F(s)
$$ (Laplace)

es la **Transformada de Laplace** de $f(t)$, siempre que la integral impropia converja.

```{admonition} Ejercicio Teórico
Calcule la Transformada de Laplace de $f(t)=e^{at}$, donde $a\in\mathbb{R}$.
```

En Python:

```{code-cell}
:tags: [Laplace1]
:tags: [hide-input]
:mystnb:
:  code_prompt_show: "Mostrar el código fuente"
:  code_prompt_hide: "Ocultar el código"
import sympy as sp

# Definir la variable simbólica 't' (tiempo) y 's' (frecuencia en el dominio de Laplace), a constante
t, s, a= sp.symbols('t s a')

# Definir algunas funciones simbólicas en el tiempo
f1 = sp.exp(a*t)  # Función exponencial

# Calcular las transformadas de Laplace
laplace_f1 = sp.laplace_transform(f1, t, s)

# Mostrar los resultados
print(f"Transformada de Laplace de exp(at): {laplace_f1[0]}")
```

**Propiedad**: La transformada de Laplace es lineal: 

$$
\mathscr{L}\{\alpha f(t)+\beta g(t)\}=\alpha\mathscr{L}\{f(t)\}+\beta\mathscr{L}\{g(t)\}
$$ 

para todo $\alpha,\beta\in\mathbb{R}$.

**Condiciones Suficientes para la Existencia de la Transformada de Laplace.** Las condiciones suficientes que garantizan la existencia de $\mathscr{L}\{f(t)\}$ son que:

1. $f(t)$ sea continua por tramos en $[0,\infty[~$.
2. $f(t)$ sea de orden exponencial $c$, para todo $t>T$.

```{figure} Laplace1.png
---
height: 150px
name: TEU
---
Continuidad por Tramos
```

```{figure} Laplace2.png
---
height: 150px
name: TEU
---
Orden Exponencial
```
Que $f(t)$ sea de orden exponencial $c$ significa que existen constantes $c,M>0$ tales que 

$$
|f(t)|\leq Me^{ct}
$$ 

es decir, $f(t)$ no crece más rápido que $e^{ct}$.

**Teorema}**: Si $f(t)$ es continua por tramos en $[0,\infty[$ y de orden exponencial $c$ entonces $\mathscr{L}\{f(t)\}$ existe para $s>c$.

**Teorema**: Si $f(t)$ es continua por tramos en $[0,\infty[$ y de orden exponencial $c$, con $\mathscr{L}\{f(t)\}=F(s)$ entonces 

$$
\lim_{s\to\infty}F(s)=0.
$$ 

```{admonition} Ejercicio Teórico
Justifique que $\mathscr{L}\{\frac{1}{t}\}$, $\mathscr{L}\{e^{t^2}\}$ no existen.
```

```{admonition} Ejercicio Teórico
¿Existe una $f(t)$ continua por tramos en $[0,\infty[$ y de orden exponencial $c$ tal que $\mathscr{L}\{f(t)\}=s^2$?
```

+++

### Transformadas Básicas

1. $\mathscr{L}\{1\}=\frac{1}{s}$.
2. $\mathscr{L}\{t\}=\frac{1}{s^2}$.
3. $\mathscr{L}\{t^n\}=\frac{n!}{s^{n+1}}$, $n\in\mathbb{N}_0$.
4. $\mathscr{L}\{e^{at}\}=\frac{1}{s-a}$, $s>a$.
5. $\mathscr{L}\{\sin(kt)\}=\frac{k}{s^2+k^2}$.
6. $\mathscr{L}\{\cos(kt)\}=\frac{s}{s^2+k^2}$.

```{admonition} Ejercicio Teórico
Calcule $\mathscr{L}\{\cos^2(t)\}$.
```

+++
