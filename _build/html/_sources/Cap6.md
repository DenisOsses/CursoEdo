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

(U6)=
# Unidad 6

+++

## Series de Fourier

+++

### Introducción

Las [series de Fourier](https://en.wikipedia.org/wiki/Fourier_series) resultan útiles en incontables aplicaciones actuales, desde la solución de la [ecuación del calor](https://en.wikipedia.org/wiki/Heat_equation), el [procesamiento de señales](https://en.wikipedia.org/wiki/Signal_processing) hasta el [procesamiento digital de imágenes](https://en.wikipedia.org/wiki/Digital_image_processing), pasando por casi cualquier campo en el que el análisis de frecuencias tenga alguna importancia.

Por ejemplo, pueden ayudar a caracterizar y comprender mejor la composición química de las estrellas, o el modo en que el aparato fonador produce el sonido.

Por analogía, podemos entenderlas como una serie de Taylor generalizada, que ya no trabaja con potencias de $x$ sino que con funciones trigonométricas seno y coseno.

### Definición

La idea esencial es que el conjunto 

$$
\left\{1,\cos\left(\dfrac{\pi x}{p}\right),\cos\left(\dfrac{2\pi x}{p}\right),\cos\left(\dfrac{3\pi x}{p}\right),\dots, \sin\left(\dfrac{\pi x}{p}\right), \sin\left(\dfrac{2\pi x}{p}\right), \sin\left(\dfrac{3\pi x}{p}\right),\dots\right\}
$$

forma una base ortogonal para el espacio vectorial de las funciones definidas en el intervalo $]-p,p[$. Luego, $f$ se puede escribir como una combinación lineal *infinita* de estas funciones 

$$
f(x)=c_0f_0(x)+c_1f_1(x)+\cdots+c_nf_n(x)+\cdots
$$

donde la igualdad entre función y serie se debe a la [convergencia uniforme](https://en.wikipedia.org/wiki/Uniform_convergence) de la sucesión de funciones del conjunto ortogonal a $f$, hecho que asumiremos en este curso. 

¿Cómo calcular los coeficientes $c_n$ de dicha serie?

**Definición**: La **serie de Fourier** de una función $f$ definida en el intervalo $]-p,p[$, está dada por: 

$$
f(x)=\dfrac{a_{0}}{2}+\sum_{n=1}^{\infty}\left(a_{n}\cos\left(\dfrac{n\pi}{p}x\right)+b_{n}\sin \left(\dfrac{n\pi}{p}x\right)\right)
$$

donde

\begin{eqnarray*}
a_{0} & = & \dfrac{1}{p}\int_{-p}^{p}f(x)dx\\
a_{n} & = & \dfrac{1}{p}\int_{-p}^{p}f(x)\cos\left(\dfrac{n\pi}{p}x\right)dx\\
b_{n} & = & \dfrac{1}{p}\int_{-p}^{p}f(x)\sin \left(\dfrac{n\pi}{p}x\right)dx
\end{eqnarray*} (Fourier1)

```{admonition} Ejercicio Teórico
Determine la serie de Fourier de $f(x)=x$ en $]-\pi,\pi[$.
```

En Python, 

```{code-cell}
:tags: [Fourier1]
:tags: [hide-input]
:mystnb:
:  code_prompt_show: "Mostrar el código fuente"
import sympy as sp

def fourier_series(f, p, n_terms=10):
    """
    Calcula el desarrollo en serie de Fourier de una función `f(x)` definida en el intervalo ]-p, p[.

    Args:
        f (sympy expression): La función a desarrollar.
        p (float): El intervalo del desarrollo ]-p, p[.
        n_terms (int): Número de términos de la serie de Fourier.

    Returns:
        sympy expression: La serie de Fourier como una expresión simbólica.
    """
    # Definir la variable simbólica
    x = sp.Symbol('x', real=True)
    
    # Coeficiente a_0
    a0 = (1 / (2 * p)) * sp.integrate(f, (x, -p, p))

    # Coeficientes a_n y b_n
    a_n = lambda n: (1 / p) * sp.integrate(f * sp.cos(n * sp.pi * x / p), (x, -p, p))
    b_n = lambda n: (1 / p) * sp.integrate(f * sp.sin(n * sp.pi * x / p), (x, -p, p))

    # Construir la serie de Fourier
    series = a0 / 2
    for n in range(1, n_terms + 1):
        series += a_n(n) * sp.cos(n * sp.pi * x / p) + b_n(n) * sp.sin(n * sp.pi * x / p)

    # Simplificar la serie
    series = sp.simplify(series)
    return series

# Definir la función f(x) y el intervalo
x = sp.Symbol('x', real=True)
f = x # Función f(x) a desarrollar
p = sp.pi  # Intervalo ]-p, p[

# Llamar a la función para calcular la serie de Fourier
n_terms = 3  # Número de términos de la serie
series = fourier_series(f, p, n_terms)

# Mostrar la serie
display(series)
```

### Convergencia de una serie de Fourier

**Teorema: Condiciones para la convergencia**: Sean $f$ y $f'$ continuas por tramos en $]-p, p[$ excepto en un número finito de puntos en el intervalo y con discontinuidades finitas sólo en esos puntos. Entonces, la serie de Fourier de $f$ en el intervalo converge a $f(x)$ en un punto de continuidad. En un punto de discontinuidad, la serie de Fourier converge hacia el promedio 

$$
\dfrac{f(x^{+})+f(x^{-})}{2}
$$

donde $f(x^{+})$ y $f(x^{-})$ denotan el límite de $f$ en $x$, por la derecha y por la izquierda.

```{admonition} Ejercicio Teórico
Analice la convergencia de la serie de Fourier del ejemplo anterior en $]-\pi,\pi[$.
```

### Extensión periódica

Una serie de Fourier no solo representa la función $f$ en el intervalo $]-p, p[$, sino que también a la extensión periódica de $f$ fuera de este intervalo, debido a la periodicidad (con periodo $2p$) de las funciones seno y coseno involucradas en su desarrollo en serie de Fourier.

```{admonition} Ejercicio Teórico
Determine la serie de Fourier de $f(x)=x$ en $]-\pi,\pi[$. Extienda esta serie periódicamente a todo $\mathbb{R}$ y analice la convergencia de la serie para $f$ en todos sus puntos.
```

```{code-cell}
:tags: [Fourier2]
:tags: [hide-input]
:mystnb:
:  code_prompt_show: "Mostrar el código fuente"
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

def plot_fourier(f, p, n_terms=10, num_points=1000):
    """
    Grafica la función original `f(x)` y su desarrollo en serie de Fourier en el intervalo ]-p, p[.

    Args:
        f (sympy expression): La función a desarrollar.
        p (float or sympy expression): El intervalo del desarrollo ]-p, p[.
        n_terms (int): Número de términos de la serie de Fourier.
        num_points (int): Número de puntos para la gráfica.
    """
    # Variables simbólicas
    x = sp.Symbol('x', real=True)
    
    # Calcular la serie de Fourier
    series = fourier_series(f, p, n_terms)
    
    # Convertir las expresiones simbólicas a funciones numéricas
    f_numeric = sp.lambdify(x, f, modules=["numpy"])
    series_numeric = sp.lambdify(x, series, modules=["numpy"])
    
    # Convertir p a un número si es simbólico
    p_numeric = float(p) if isinstance(p, sp.Basic) else p
    
    # Crear puntos para graficar
    x_vals = np.linspace(-p_numeric, p_numeric, num_points)
    f_vals = f_numeric(x_vals)
    series_vals = series_numeric(x_vals)
    
    # Graficar
    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, f_vals, label='f(x) (Original)', color='blue')
    plt.plot(x_vals, series_vals, label=f'Serie de Fourier (n={n_terms})', color='red', linestyle='--')
    plt.title(f'Serie de Fourier de f(x)={f} en ]-{p_numeric}, {p_numeric}[')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
    plt.axvline(0, color='black', linewidth=0.8, linestyle='--')
    plt.legend()
    plt.grid()
    plt.show()

# Ejemplo de uso
x = sp.Symbol('x', real=True)
f = x  # Función original
p = np.pi  # Intervalo ]-p, p[
n_terms = 3  # Número de términos en la serie de Fourier

plot_fourier(f, p, n_terms)
```

Note que $a_n=0$ para todo $n=0,1,2,\ldots$ y los únicos coeficientes que aparecen en la serie son los coeficientes $b_n$. ¿Tiene esto alguna razón?