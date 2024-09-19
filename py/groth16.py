import numpy as np
import galois
from functools import reduce
from py_ecc.bn128 import curve_order, G1, G2, multiply, add, neg
import random

def evaluate_poly(coeffs, powers_of_tau):
    final_point = None
    for i in range(len(coeffs)):
        if coeffs[i] == 0:
            continue
        if final_point == None:
            final_point = multiply(powers_of_tau[i], int(coeffs[i]))
            continue
        final_point = add(final_point, multiply(powers_of_tau[i], int(coeffs[i])))
    return final_point

def evaluate_witness(witness, powers_of_tau):
    witness = [int(x) for x in witness]

    final_point = None
    for i in range(len(witness)):
        if powers_of_tau[i] == 0:
            continue
        if final_point == None:
            final_point = multiply(powers_of_tau[i], witness[i])
            continue
        final_point = add(final_point, multiply(powers_of_tau[i], witness[i]))
    return final_point


# # Original Equation:
# # out = x**4 - 5*y**2*x**2
# #
# # Constraints:
# # v1 = x * x
# # v2 = v1 * v1
# # v3 = -5y * y
# # -v2 + out = v3*v1
# #
# # Witness:
# # [1, out, x, y, v1, v2, v3]


# -----
# Define Matrices L, R, and O
# -----

L = np.array(
    [
        [0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, -5, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1],
    ]
)

R = np.array(
    [
        [0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0],
    ]
)

O = np.array(
    [
        [0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 1],
        [0, 1, 0, 0, 0, -1, 0],
    ]
)


# -----
# Verify R1CS
# -----

print("Verifying R1CS...")

x = 4
y = -2
v1 = x * x
v2 = v1 * v1  # x^4
v3 = -5 * y * y
out = v3 * v1 + v2  # -5y^2 * x^2

witness = np.array([1, out, x, y, v1, v2, v3])

assert all(
    np.equal(np.matmul(L, witness) * np.matmul(R, witness), np.matmul(O, witness))
), "not equal"
print("R1CS verified!")


# -----
# Redefine Matrices to Eliminate Negative Values
# -----

def remove_negatives(matrix):
    return np.where(matrix < 0, matrix + curve_order, matrix)


L = remove_negatives(L)
R = remove_negatives(R)
O = remove_negatives(O)


# -----
# Recalculate Witness and Verify Galois Field
# -----

print("Verifying Galois Field...")

GF = galois.GF(curve_order)

L_galois = GF(L)
R_galois = GF(R)
O_galois = GF(O)

x = GF(4)
y = GF(curve_order - 2)  # we are using 79 as the field size, so 79 - 2 is -2
v1 = x * x
v2 = v1 * v1  # x^4
v3 = GF(curve_order - 5) * y * y
out = v3 * v1 + v2  # -5y^2 * x^2

witness = GF(np.array([1, out, x, y, v1, v2, v3]))

assert all(
    np.equal(
        np.matmul(L_galois, witness) * np.matmul(R_galois, witness),
        np.matmul(O_galois, witness),
    )
), "not equal"
print("Galois Field verified!")


# -----
# Interpolate Polynomial
# -----

def interpolate_column(col):
    xs = GF(np.array([1, 2, 3, 4]))
    return galois.lagrange_poly(xs, col)


# axis 0 is the columns. apply_along_axis is the same as doing a for loop over the columns and collecting the results in an array
U_polys = np.apply_along_axis(interpolate_column, 0, L_galois)
V_polys = np.apply_along_axis(interpolate_column, 0, R_galois)
W_polys = np.apply_along_axis(interpolate_column, 0, O_galois)


# -----
# Calculate h(x)
# -----


def inner_product_polynomials_with_witness(polys, witness):
    mul_ = lambda x, y: x * y
    sum_ = lambda x, y: x + y
    return reduce(sum_, map(mul_, polys, witness))


term_1 = inner_product_polynomials_with_witness(U_polys, witness)
term_2 = inner_product_polynomials_with_witness(V_polys, witness)
term_3 = inner_product_polynomials_with_witness(W_polys, witness)
term_3_pub = inner_product_polynomials_with_witness(W_polys, witness[:2])
term_3_priv = inner_product_polynomials_with_witness(W_polys, witness[2:])

# t = (x - 1)(x - 2)(x - 3)(x - 4)
t = (
    galois.Poly([1, (curve_order - 1)], field=GF)
    * galois.Poly([1, (curve_order - 2)], field=GF)
    * galois.Poly([1, (curve_order - 3)], field=GF)
    * galois.Poly([1, (curve_order - 4)], field=GF)
)

h = (term_1 * term_2 - term_3) // t


# -----
# Verify QAP
# -----

print("Verifying QAP...")

assert term_1 * term_2 == term_3 + h * t, "division has a remainder"
print("QAP verified!")


# -----
# Trusted Setup
# -----

print("Performing trusted setup...")

alpha = 2
beta = 3
gamma = 4
delta = 5
tau = 6

alpha_G1 = multiply(G1, alpha)  # Random shift for A
beta_G1 = multiply(G1, beta)
beta_G2 = multiply(G2, beta)  # Random shift for B
gamma_G2 = multiply(G2, gamma)  # Gamma EC point
delta_G2 = multiply(G2, delta)  # Delta EC point

powers_of_tau_G1 = [
    multiply(G1, tau**4),
    multiply(G1, tau**3),
    multiply(G1, tau**2),
    multiply(G1, tau),
    G1,
]  # Powers of tau for A
powers_of_tau_G2 = [
    multiply(G2, tau**4),
    multiply(G2, tau**3),
    multiply(G2, tau**2),
    multiply(G2, tau),
    G2,
]  # Powers of tau for B

sum_polys = U_polys * beta + V_polys * alpha + W_polys

def pub_tau(n):
    if n == 0:
        return 0
    return multiply(
        evaluate_poly(n.coefficients(5)[:2], powers_of_tau_G1[:2]),
        (curve_order - gamma),
    )

def priv_tau(n):
    if n == 0:
        return 0
    return multiply(
        evaluate_poly(n.coefficients(5)[2:], powers_of_tau_G1[2:]),
        (curve_order - delta),
    )

p_tau_pub = list(map(pub_tau, sum_polys))  # Powers of tau for public inputs
p_tau_priv = list(map(priv_tau, sum_polys))  # Powers of tau for private inputs

powers_of_tau_TG1 = evaluate_poly(
    t.coefficients(5), powers_of_tau_G1
)  # Powers of tau for h(tau)t(tau)


# -----
# Evaluating Polynomials
# -----

print("Evaluating Polynomials...")


def evaluate_HT(coeffs, T):
    final_point = None
    for i in range(len(coeffs)):
        if coeffs[i] == 0:
            continue
        if final_point == None:
            final_point = multiply(T, int(coeffs[i]))
            continue
        final_point = add(final_point, multiply(T, int(coeffs[i])))
    return final_point


r = random.randint(1, 100)  # max lower than curve_order just to speed up computation
s = random.randint(1, 100)  # max lower than curve_order just to speed up computation

delta_G1 = multiply(G1, delta)

print("Calculating A, B...")

A = add(
    add(evaluate_poly(term_1.coefficients(5), powers_of_tau_G1), alpha_G1),
    multiply(delta_G1, r),
)
B2 = add(
    add(evaluate_poly(term_2.coefficients(5), powers_of_tau_G2), beta_G2),
    multiply(delta_G2, s),
)
B1 = add(
    add(evaluate_poly(term_2.coefficients(5), powers_of_tau_G1), beta_G1),
    multiply(delta_G1, s),
)

print("Calculating C...")

# Evaluate h(tau) to get coefficients to use a scalars to mulitply with G points for t
HT = evaluate_HT(h.coefficients(5), powers_of_tau_TG1)

C_prime = evaluate_witness(witness[2:], p_tau_priv)

C = add(
    add(add(add(C_prime, HT), multiply(A, s)), multiply(B1, r)),
    neg(multiply(delta_G1, r * s)),
)

def strG(G):
  return (
    repr(G)
    .replace("(", "[")
    .replace(")", "]")
    .replace("[[", "[")
    .replace("]]", "]")
    .replace("], ", "],\n\t")
  )

def print_G1_point(point, name):
    print("Verify.G1Point memory %s = Verify.G1Point(%s, %s);" % (name, point[0], point[1]))

def print_G2_point(point, name):
    print("Verify.G2Point memory %s = Verify.G2Point(\n\t%s\n);" % (name, strG(point)))

def print_GF_vector(vector, name):
    vector = [int(x) for x in vector]
    print("uint256[%s] memory %s = [" % (len(vector), name))
    for i in range(len(vector)):
        print("%s, " % (vector[i]))
    print("];")

def print_G1_vector(vector, name):
    print("Verify.G1Point[] memory %s = new Verify.G1Point[](%s);" % (name, len(vector)))
    for i in vector:
        if (i == 0):
            print("%s[%s] = Verify.G1Point(%s, %s);" % (name, vector.index(i), "0", "0"))
            continue
        print("%s[%s] = Verify.G1Point(%s, %s);" % (name, vector.index(i), i[0], i[1]))

def print_G2_vector(vector, name):
    print("Verify.G2Point[] memory %s = new Verify.G2Point[](%s);" % (name, len(vector)))
    for i in vector:
        print("%s[%s] = Verify.G2Point(\n\t%s\n);" % (name, vector.index(i), strG(i)))

print("\nGive to Verifier\n----------------")
print_G1_point(A, "a")
print_G2_point(B2, "b")
print_G1_point(alpha_G1, "alpha")
print_G2_point(beta_G2, "beta")
print_G1_vector(p_tau_pub[:2], "pTauPub")
print_GF_vector(witness[:2], "pubInputs")
print_G2_point(gamma_G2, "gamma")
print_G1_point(C, "c")
print_G2_point(delta_G2, "delta")
