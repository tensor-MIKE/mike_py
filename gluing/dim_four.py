from dim_four.theta import ThetaPoint, ThetaPointCompressed, ThetaStructure


def proj_batch_pseudo_inversion(v):
    n = len(v)
    prefix = [v[0]]
    for i in range(1, n):
        prefix.append(prefix[-1] * v[i])
    suffix = [1]
    for i in range(1, n):
        suffix.append(suffix[-1] * v[n - i])
    out = [suffix[-1]]
    for i in range(1, n):
        out.append(suffix[n - 1 - i] * prefix[i - 1])
    return out


def hadamard_like(a, b, c, d):
    """
    Hadamard like function used for apply_gluing_basis_change
    """
    t1, t2 = a + b, a - b
    t3, t4 = c + d, c - d
    return t1 + t3, t2 + t4, t1 - t3, t2 - t4


def apply_gluing_basis_change(coords):
    """
    Apply the appropriate change of basis to the dim four point,
    which is essentially four independent dim two hadamrd like things.
    """
    c = coords
    o = [None] * 16
    o[0], o[1], o[2], o[3] = hadamard_like(c[0], c[5], c[10], c[15])
    o[4], o[5], o[6], o[7] = hadamard_like(c[4], c[1], c[14], c[11])
    o[8], o[9], o[10], o[11] = hadamard_like(c[2], c[7], c[8], c[13])
    o[12], o[13], o[14], o[15] = hadamard_like(c[6], c[3], c[12], c[9])
    return o


def dim_two_to_dim_four(P, Q):
    """
    Product of two dimension-two theta points, with the gluing basis change
    applied.
    """
    px, py, pz, pt = P.coords
    qx, qy, qz, qt = Q.coords

    c = [0] * 16
    c[0], c[8], c[4], c[12] = px * qx, px * qz, px * qy, px * qt
    c[2], c[10], c[6], c[14] = pz * qx, pz * qz, pz * qy, pz * qt
    c[1], c[9], c[5], c[13] = py * qx, py * qz, py * qy, py * qt
    c[3], c[11], c[7], c[15] = pt * qx, pt * qz, pt * qy, pt * qt

    return ThetaPoint(apply_gluing_basis_change(c))


def solve_HIIP_gluing(kernel):
    """Inverse dual null point of the codomain, plus the inverse of T3."""
    # 0--(4,2)--6
    leg_4 = [
        kernel[0][4] * kernel[0][6],
        kernel[0][0] * kernel[0][6],
        kernel[0][0] * kernel[0][2],
    ]
    # 0--(8,1)--9
    leg_8 = [
        kernel[1][8] * kernel[1][9],
        kernel[1][0] * kernel[1][9],
        kernel[1][0] * kernel[1][1],
    ]
    # 0--(12,3)--15
    leg_12 = [
        kernel[2][12] * kernel[2][15],
        kernel[2][0] * kernel[2][15],
        kernel[2][0] * kernel[2][3],
    ]

    inv_null_comp = ThetaPointCompressed([0] * 10)

    tmp_8_12 = leg_8[0] * leg_12[0]
    inv_null_comp[0] = tmp_8_12 * leg_4[0]
    inv_null_comp[2] = tmp_8_12 * leg_4[1]
    inv_null_comp[5] = tmp_8_12 * leg_4[2]

    tmp_4_12 = leg_4[0] * leg_12[0]
    inv_null_comp[1] = tmp_4_12 * leg_8[1]
    inv_null_comp[7] = tmp_4_12 * leg_8[2]

    tmp_4_8 = leg_4[0] * leg_8[0]
    inv_null_comp[3] = tmp_4_8 * leg_12[1]
    inv_null_comp[9] = tmp_4_8 * leg_12[2]

    inv_theta_null = inv_null_comp.decompress()

    # t3[6] comes out zero here, because inv_theta_null[10] = inv_null_comp[4]
    # = 0; corrected just below.
    t3 = [kernel[0][i] * inv_theta_null[i] for i in (0, 1, 2, 3, 8, 9, 10, 15)]

    # The true value is
    #   t3[6] = t3[0] * kernel[3][2] * inv_null[2] / (kernel[3][8] * inv_null[8])
    # so clear the denominator across the other seven entries: 7M, not 1I.
    t3[6] = t3[0] * kernel[3][2] * inv_theta_null[2]
    correction_den = kernel[3][8] * inv_theta_null[8]
    for i in (0, 1, 2, 3, 4, 5, 7):
        t3[i] *= correction_den

    t3 = proj_batch_pseudo_inversion(t3)

    # Expand 8 -> 16 with the duplication inv_T3[8i+j] = inv_T3[8i+j+4]
    inv_t3_coords = [0] * 16
    for i in range(2):
        for j in range(4):
            inv_t3_coords[8 * i + j] = t3[4 * i + j]
            inv_t3_coords[8 * i + j + 4] = t3[4 * i + j]

    return inv_theta_null, ThetaPoint(inv_t3_coords)


NON_ZERO = [0, 1, 2, 3, 4, 6, 8, 9, 12, 15]


def dim_four_gluing_isogeny(gluing_kernel, image_plus, image_minus):
    """Returns (codomain, chain_kernel)."""
    kernel = [dim_two_to_dim_four(u, v).square().hadamard() for (u, v) in gluing_kernel]

    inv_null_point_dual, inv_T3 = solve_HIIP_gluing(kernel)

    non_zero = proj_batch_pseudo_inversion([inv_null_point_dual[i] for i in NON_ZERO])
    codomain_null_point = ThetaPoint([0] * 16)
    for value, idx in zip(non_zero, NON_ZERO):
        codomain_null_point[idx] = value
    codomain_null_point = codomain_null_point.hadamard()
    codomain = ThetaStructure.from_null_point(codomain_null_point)

    # t1 from image_plus, t2 from image_minus, image = H(inv_T3 * H(t1 * t2))
    chain_kernel = []
    for (pu, pv), (mu, mv) in zip(image_plus, image_minus):
        t1 = dim_two_to_dim_four(pu, pv)
        t2 = dim_two_to_dim_four(mu, mv)
        t3 = t1.pointwise_mul(t2).hadamard()
        t3 = t3.pointwise_mul(inv_T3).hadamard()
        chain_kernel.append(t3)

    return codomain, chain_kernel
