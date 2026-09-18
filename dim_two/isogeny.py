from dim_two.theta import ThetaPointDimTwo, ThetaStructureDimTwo


def two_isogeny(T1, T2, image_points, hadamard):
    """
    (2,2)-isogeny from the 8-torsion (T1, T2) above the kernel.

    `hadamard` is a pair of bools: the first moves the input to the dual, the
    second moves the output there.

    Returns (codomain_null_point, images), where `images` is a new list in the
    same order as `image_points`.
    """
    if hadamard[0]:
        T1, T2 = T1.hadamard(), T2.hadamard()

    xA, xB, _, _ = T1.square().hadamard().coords
    zA, tB, zC, tD = T2.square().hadamard().coords

    xAtB = xA * tB
    zAxB = zA * xB
    zCtD = zC * tD

    A = zA * xAtB
    B = tB * zAxB
    C = zC * xAtB
    D = tD * zAxB
    null_point = ThetaPointDimTwo([A, B, C, D])

    # C_inv, D_inv are deliberately D, C
    inverses = ThetaPointDimTwo([xB * zCtD, xA * zCtD, D, C])

    if hadamard[1]:
        null_point = null_point.hadamard()

    codomain = ThetaStructureDimTwo(null_point)

    images = []
    for P in image_points:
        if hadamard[0]:
            P = P.hadamard()
        P = P.square()
        P = P.hadamard()
        P = P.pointwise_mul(inverses)
        if hadamard[1]:
            P = P.hadamard()
        images.append(P)

    return codomain, images
