import sympy as sp

# Definición de matrices
A = sp.Matrix([
    [1, 3, -1],
    [3, 4, -2],
    [-1, -2, 2]
])

B = sp.Matrix([
    [1, 0, 0],
    [-3, 1, 0],
    [4, -7, 1]
])

C = sp.Matrix([
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [1, 0, 0, 0],
    [0, 1, 0, 0]
])

D = sp.Matrix([
    [0, 1, 0, 0],
    [1, 0, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0]
])

E = sp.diag(1, 1, -1, -1)

# Colección de matrices
matrices = {"A": A, "B": B, "C": C, "D": D, "E": E}

# Función para verificar ortogonalidad (también para vectores complejos)
def verificar_ortogonalidad(vectores):
    n = len(vectores)
    resultados = []
    for i in range(n):
        for j in range(i + 1, n):
            prod = (vectores[i].conjugate().T @ vectores[j])[0]
            resultados.append((i + 1, j + 1, sp.simplify(prod)))
    return resultados

# Cálculo principal
for nombre, M in matrices.items():
    print(f"\n{'='*50}")
    print(f"Matriz {nombre}:")
    sp.pprint(M)
    
    # Cálculo de autovalores y autovectores
    autovectores = M.eigenvects()
    print("\nAutovalores y autovectores:")
    for val, mult, vects in autovectores:
        print(f"\nλ = {val}  (multiplicidad {mult})")
        for v in vects:
            sp.pprint(v)
    
    # Normalizar todos los autovectores
    vectores_normalizados = []
    for _, _, vects in autovectores:
        for v in vects:
            norm = sp.sqrt((v.conjugate().T @ v)[0])
            if norm != 0:
                vectores_normalizados.append(v / norm)
    
    # Verificar ortogonalidad
    if len(vectores_normalizados) > 1:
        print("\nVerificación de ortogonalidad ⟨vᵢ,vⱼ⟩:")
        ortogonales = verificar_ortogonalidad(vectores_normalizados)
        for i, j, prod in ortogonales:
            print(f"v{i}·v{j} = {sp.simplify(prod)}")
    else:
        print("\nNo hay suficientes autovectores para verificar ortogonalidad.")
