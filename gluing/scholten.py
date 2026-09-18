from dim_one.add_components import add_components
from dim_two.theta import ThetaPointDimTwo


def fp2_norm(a):
    return a[0] ** 2 + a[1] ** 2


def to_x_only(P):
    """
    This wrapper is mainly here to catch bugs, as we're using SageMath
    elliptic curve arithmetic at the moment, so Z = 1 always currently
    """
    assert P[2] == 1, "expected an affine point"
    return P[0], P[2]


def symmetric_action_change_of_basis(P8):
    """
    Computes the matrix coefficients `mi` for the Scholten change of basis
    given a point `P8` of order 8
    """
    alpha, beta = to_x_only(2 * P8)

    n_alpha = fp2_norm(alpha)
    n_beta = fp2_norm(beta)

    m0 = n_alpha + n_beta
    m2 = n_beta - n_alpha

    t0 = 2 * (alpha * beta.conjugate())
    m1, m3 = t0[0], t0[1]

    return [m0, m1, m2, m3]


def point_to_squared_theta_point(mi, P):
    """
    Compute the theta point from (P, P^sigma) on E x E^sigma
    """
    x, z = to_x_only(P)
    m0, m1, m2, m3 = mi

    n_x = fp2_norm(x)
    n_z = fp2_norm(z)

    t0 = n_x + n_z
    t1 = n_x - n_z

    tmp = 2 * (x * z.conjugate())
    t2, t3 = tmp[0], tmp[1]

    X = m0 * t0 - m1 * t2
    Y = m0 * t2 - m1 * t0
    Z = m2 * t1 - m3 * t3
    T = m2 * t3 + m3 * t1

    return ThetaPointDimTwo([X**2, Y**2, Z**2, -(T**2)])


def j_precomputation(mi, P8, Q8):
    """
    The theta point J = (y : y : x : x) used for computing images
    """
    T1 = point_to_squared_theta_point(mi, P8).hadamard()
    T2 = point_to_squared_theta_point(mi, Q8).hadamard()

    xA, _, yC, yD = T1.coords
    zA, _, zY, tD = T2.coords
    assert yD == 0, "H(S(T1))[3] should vanish"
    assert tD == 0, "H(S(T2))[3] should vanish"

    x = xA * zY
    y = yC * zA

    return ThetaPointDimTwo([y, y, x, x])


def superglue_precomputation(mi):
    """
    The five superglue constants (A, B, C, D, E).
    """
    m0, m1, m2, m3 = mi

    m0_sqr = m0**2
    m1_sqr = m1**2
    m2_sqr = m2**2
    m3_sqr = m3**2

    assert m0_sqr - m3_sqr == m1_sqr + m2_sqr

    a = m1_sqr + m2_sqr
    b = m0_sqr - m1_sqr
    c = m0_sqr - m2_sqr
    d = 4 * (m0 * m1)
    e = 4 * (m2 * m3)

    return [a, b, c, d, e]


def scholten_image(E, P, shift, precomp, J, change_of_basis=True):
    """
    Image of a single curve point under the Scholten isogeny, eventually
    we push everything for the `change_of_basis`, but for some points we
    set this bool to false to postpone this.
    """
    u, v, w = add_components(E, P, shift)

    u_abs = fp2_norm(u)
    v_abs = fp2_norm(v)
    w_abs = fp2_norm(w)

    u_sqr = u**2
    v_sqr = v**2
    w_sqr = w**2

    t0 = u_sqr - v_sqr
    tp = t0 + w_sqr
    tm = t0 - w_sqr

    tp_abs = fp2_norm(tp)
    tm_abs = fp2_norm(tm)

    # U + iV = (u^2 - v^2 + w^2) * conj(u*w)
    t0 = tp * (u * w).conjugate()
    U, V = t0[0], t0[1]

    a, b, c, d, e = precomp
    t0 = 4 * (u_abs * w_abs)
    t1 = d * U
    t2 = e * V

    alpha = a * tp_abs + c * t0 - t1 - t2
    gamma = c * tp_abs + a * t0 - t1 + t2
    beta = b * tm_abs
    delta = -b * 4 * (v_abs * w_abs)

    image = ThetaPointDimTwo([alpha, beta, gamma, delta])
    image = image.pointwise_mul(J)
    image = image.hadamard()

    if change_of_basis:
        image = image.scholten_change_of_basis()

    return image


def U_from_V(V):
    """
    Rather than push points through the isogeny to compute U values,
    we can infer them directly from the V images
    """
    x, y, z, t = V.coords
    t0, t1 = x + t, x - t
    t2, t3 = y + z, y - z
    return ThetaPointDimTwo([t1 + t2, t0 + t3, t0 - t3, t2 - t1])


def scholten_isogeny(E, scholten_basis, diag_basis, gluing_basis, chain_kernel):
    """
    Computes the image of diag_basis, gluing_basis, chain_kernel through the scholten
    isogeny from E x E^sigma to A2
    """
    P8, Q8 = scholten_basis
    assert (4 * Q8).x() == 0

    # Precompute change of basis matrix, image computation and superglue constants
    mi = symmetric_action_change_of_basis(P8)
    J = j_precomputation(mi, P8, Q8)
    image_precomp = superglue_precomputation(mi)

    # Wrapper to make the below less verbose
    def image(P, change_of_basis=True):
        return scholten_image(E, P, P8, image_precomp, J, change_of_basis)

    # Diagonal kernel: the Ui can be computed directly from Vi providing we
    # apply the `scholten_change_of_basis` later.
    V1, V2 = [image(P, False) for P in diag_basis]
    U1, U2 = U_from_V(V2), U_from_V(V1)
    diagonal_kernel = [P.scholten_change_of_basis() for P in (U1, U2, V1, V2)]

    # Gluing kernel is directly pushed through
    gluing_kernel = [[image(P) for P in pair] for pair in gluing_basis]

    # For each chain kernel pair (P, Q) push (P + T, Q + T) and (P - T, Q - T)
    # which will be combined in the later dim four gluing.
    T = gluing_basis[0]
    image_plus = [[image(P + T[0]), image(Q + T[1])] for (P, Q) in chain_kernel]
    image_minus = [[image(P - T[0]), image(Q - T[1])] for (P, Q) in chain_kernel]

    return diagonal_kernel, gluing_kernel, image_plus, image_minus
