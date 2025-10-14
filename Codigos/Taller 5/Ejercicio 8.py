import sympy as sp

# Datos del sistema
m1, m2, m3 = 1, 2, 1
r1 = sp.Matrix([1, 1, -2])
r2 = sp.Matrix([-1, -1, 0])
r3 = sp.Matrix([1, 1, 2])

# Lista de masas y posiciones
m = [m1, m2, m3]
r = [r1, r2, r3]

# Cálculo del tensor de inercia (definición general)
# I = sum_i m_i * [(r_i·r_i) * I3 - r_i ⊗ r_i]
I_total = sp.zeros(3)
I3 = sp.eye(3)

for mi, ri in zip(m, r):
    I_total += mi * ((ri.dot(ri)) * I3 - ri * ri.T)

print("Tensor de inercia I =")
sp.pprint(sp.simplify(I_total))

# Diagonalización del tensor de inercia
autovectores = I_total.eigenvects()

print("\nAutovalores y autovectores (Ejes principales):")
for val, mult, vects in autovectores:
    print(f"\nλ = {val} (Momento principal)")
    for v in vects:
        print("Eje principal →")
        sp.pprint(v.normalized())

# Verificación de ortogonalidad de los ejes principales
vects = [v.normalized() for _, _, V in autovectores for v in V]
print("\nVerificación de ortogonalidad ⟨vᵢ,vⱼ⟩:")
for i in range(len(vects)):
    for j in range(i+1, len(vects)):
        prod = (vects[i].dot(vects[j])).simplify()
        print(f"v{i+1}·v{j+1} = {prod}")
