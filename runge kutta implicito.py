from numpy import zeros, linspace, identity, ones, stack, array, sqrt
from numpy.linalg import solve
from matplotlib.pyplot import plot, show


def runge_kutta_implicito_sistema(f, fy, t0, tf, N, y0, tol, itermax, A, b):
    # Creamos el paso de malla
    h = (tf - t0)/float(N)

    # Dimensión de y y de z
    n = len(y0)

    # Definimos las matrices de la solucion
    # numerica e incluimos las condiciones iniciales
    y = zeros([n, N+1])
    y[:, 0] = y0

    y_iters = zeros(N)

    # Creamos el vector de tiempos
    t = linspace(t0, tf, N+1)

    # Calculamos el número de etapas del método RK
    s = len(b)

    # Vamos utilizar mucho lo siguiente
    ns = n*s
    I = identity(s)
    E = -h*A

    for k in range(0, N):
        cont = 0
        d = 1.0 + tol
        # Creamos los y_k^(1),...,y_k^(s),
        # z_k^(1),...,z_k^(s) iniciales
        yk_iter = y[:, k]*ones((s, n))

        # Creamos un único vector con todas las variables
        uk_iter = yk_iter.flatten(order='F')

        # Inciamos el método de Newton
        while (d >= tol) and (cont < itermax):
            F = zeros((s, n))
            Fy = zeros(s, dtype='O')
            # Evaluamos f, g, fy, fz, gy, gz
            for i in range(0, s):
                F[i] = f(yk_iter[i])
                Fy[i] = fy(yk_iter[i])

            # Evaluacion de cada etapa
            Fy = stack(Fy, axis=-1)

            # Creamos y rellenamos la matriz Jacobiana
            DH = zeros((ns, ns))
            for i in range(0, n):
                for j in range(0, n):
                    if i == j:
                        DH[i*s:(i+1)*s, j*s:(j+1)*s] = I + E*Fy[i, j]
                    else:
                        DH[i*s:(i+1)*s, j*s:(j+1)*s] = E*Fy[i, j]

            # Creamos H
            H = zeros(ns)
            for i in range(n):
                H[i*s: (i+1)*s] = yk_iter[:, i] - y[i, k]*ones(s) + E@F[:, i]

            uk_new = uk_iter - solve(DH, H)
            d = max(abs(uk_iter - uk_new))
            uk_iter = uk_new
            yk_iter = uk_iter[:ns].reshape((n, s)).transpose()

            cont += 1

        if cont == itermax:
            print('El metodo no va bien (y_{k+1})')
        else:
            # Calculamos y_{k+1}
            for i in range(0, n):
                y[i, k+1] = y[i, k] + h*b@F[:, i]

            # Iteraciones del método Newton para y_{k+1}
            y_iters[k] = cont

    return (t, y)


A = array([
    [5/36, 2/9 - sqrt(15)/15, 5/36 - sqrt(15)/30],
    [5/36 + sqrt(15)/24, 2/9, 5/36 - sqrt(15)/24],
    [5/36 + sqrt(15)/30, 2/9 + sqrt(15)/15, 5/36]
])
b = array([5/18, 4/9, 5/18])

itermax = 200
y0 = array([1, 2])
t0 = 0.0
tf = 2
tol = 1.e-12
N = 50


def f(y):
    f1, f2 = y
    return array([f1, f2])


def fy(y):
<<<<<<< HEAD
  return array([[1,0], [0,2]])


(t,y) = runge_kutta_implicito_sistema(f,fy,t0,tf,N,y0,tol,itermax,A,b,c)

plot(t,y[0],t,y[1])
=======
    return array([[1, 0], [0, 2]])
>>>>>>> 1935629d00010d08af39658eeeb08df96e076683


t, y = runge_kutta_implicito_sistema(
    f, fy, t0, tf,
    N, y0, tol, itermax, A, b
)

plot(t, y[0], t, y[1])
show()
