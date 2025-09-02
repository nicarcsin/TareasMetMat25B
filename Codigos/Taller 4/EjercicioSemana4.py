import sympy
from sympy import *

# Variables y límites

x, y = symbols('x y')
a1, a2 = -1, 1

# Bases monomiales

base_P = [1, x, x**2]
base_G = [1, y, y**2]

# Defino productos internos

def producto_x(f, g):
    
    return simplify(integrate(f * g, (x, a1, a2)))

def producto_y(f, g):
    
    return simplify(integrate(f * g, (y, a1, a2)))

# Ortogonalización para encontrar Legendre 

def gram_schmidt(base, producto):
    ort = []
    for v in base:
        u = sympy.simplify(v)  # empieza con v
        for o in ort:
            num = producto(v, o)
            den = producto(o, o)
            u = sympy.simplify(u - (num/den) * o)
        ort.append(sympy.simplify(u))
    return ort

# Bases para x y para y:

ortogonales1 = gram_schmidt(base_P, producto_x)  # para variable x
ortogonales2 = gram_schmidt(base_G, producto_y)  # para variable y

#Función para obtener la matriz de coeficientes en la base monomial tensorial

def tensor_coeff_matrix_monomial(P_poly, G_poly):
    poly = expand(P_poly * G_poly)  # polinomio en x,y

    # Base tensorial monomial (orden: [1, y, y^2] como columnas por cada fila de x^i)
    base_T = [p * g for p in base_P for g in base_G]
    coeffs = {mon: Poly(poly, x, y).coeff_monomial(mon) for mon in base_T}
    M = Matrix([[coeffs[p * g] for g in base_G] for p in base_P])
    return M

# Función para obtener la matriz de coeficientes en la base ortogonal (producto tensorial de ortogonales) ---
def tensor_coeff_matrix_ortogonal(P_poly, G_poly):
    poly = expand(P_poly * G_poly)  # polinomio en x,y
    M = Matrix([[0 for _ in ortogonales2] for __ in ortogonales1])
    for i, pi in enumerate(ortogonales1):
        for j, gj in enumerate(ortogonales2):

            # Numerador: doble integral de poly(x,y) * pi(x) * gj(y)

            numer = integrate(poly * pi * gj, (x, a1, a2), (y, a1, a2))

            # Denominador: producto de normas 

            denom = producto_x(pi, pi) * producto_y(gj, gj)
            
            M[i, j] = simplify(numer / denom)
    return M

# Ejemplo general

a, b, c, d, e, f = symbols('a b c d e f')
P = a * x**2 + b * x + c
G = d * y**2 + e * y + f

Cmat = tensor_coeff_matrix_monomial(P, G)
print("a) Componentes del tensor (base monomial):")
print(Cmat)

# Ejemplo concreto

pol1 = x**2 + x + 3
pol2 = y + 1

tensor1 = expand(pol1 * pol2)
print("\nb) Polinomio tensor (expandido):")
print(tensor1)

comptensor1 = tensor_coeff_matrix_monomial(pol1, pol2)
print("\nComponentes del tensor en base monomial:")
print(comptensor1)

# Bases de Legendre

print("\nBase ortogonal en x (ortogonales1):")
for i, p in enumerate(ortogonales1):
    print(f"p_{i}(x) =", sympy.simplify(p))

print("\nBase ortogonal en y (ortogonales2):")
for j, g in enumerate(ortogonales2):
    print(f"g_{j}(y) =", sympy.simplify(g))

# Componentes en la base ortogonal tensorial

Cmat_ort = tensor_coeff_matrix_ortogonal(pol1, pol2)
print("\nd) Componentes en la base ortogonal tensorial:")
print(Cmat_ort)

