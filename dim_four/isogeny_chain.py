from .theta import ThetaPoint, ThetaStructure


def compute_inv_dual_codomain(P8: ThetaPoint, Q8: ThetaPoint) -> ThetaPoint:
    r1 = P8[12]
    r2 = r1 * P8[5]
    r3 = r2 * Q8[1]
    r4 = r3 * P8[9]
    r5 = r4 * P8[11]
    r6 = r5 * Q8[15]
    r7 = r6 * Q8[14]
    r8 = r7 * P8[6]

    # Precompute the suffix chain
    # 7M
    b8 = P8[0]
    b7 = b8 * P8[2]
    b6 = b7 * Q8[6]
    b5 = b6 * Q8[7]
    b4 = b5 * P8[15]
    b3 = b4 * P8[13]
    b2 = b3 * Q8[9]
    b1 = b2 * P8[1]

    # default point, all zeros, but need to make the ThetaPoint type...
    coords = [0] * 16

    # First element
    # 1M
    coords[0] = r8 * P8[4]

    # Last element required an extra mul
    # 2M
    tmp = P8[8] * P8[5]
    coords[12] = tmp * b2

    # Remaining coordinates
    # 8M
    coords[5] = r1 * b1
    coords[1] = r2 * b2
    coords[9] = r3 * b3
    coords[11] = r4 * b4
    coords[15] = r5 * b5
    coords[14] = r6 * b6
    coords[6] = r7 * b7
    coords[4] = r8 * b8

    # Direct copies from symmetry
    coords[3] = coords[12]
    coords[10] = coords[5]
    coords[8] = coords[1]
    coords[13] = coords[11]
    coords[7] = coords[14]
    coords[2] = coords[4]

    return ThetaPoint(coords)


def isogeny_codomain(
    K1: ThetaPoint, K2: ThetaPoint
) -> tuple[ThetaStructure, ThetaPoint]:
    """
    Given two points of order eight, compute the codomain of the isogeny
    with kernel <[4]K1, [4]K2> together with a constant used for images.
    """
    P8 = K1.square()
    P8 = P8.hadamard()

    Q8 = K2.square()
    Q8 = Q8.hadamard()

    inv_dual_null = compute_inv_dual_codomain(P8, Q8)
    codomain = ThetaStructure.from_inv_dual(inv_dual_null)

    return codomain, inv_dual_null


def isogeny_eval(point, inv_dual_null):
    """
    Compute the image of P under the isogeny f : A -> B with
    inv_dual_null being the inverse of the dual of the null
    point of B
    """
    P = point.square()
    P = P.hadamard()
    P = P.pointwise_mul(inv_dual_null)
    P = P.hadamard()
    return P


def dim_four_isogeny_chain(domain: ThetaStructure, kernel: tuple[ThetaPoint], e: int):
    """
    Compute a 2^e-isogeny from `domain` as a chain of (2,2,2,2)-isogenies,
    using a balanced strategy, and return the codomain of the isogeny
    """
    OA = domain
    T1, T2 = kernel

    # Stack of (point, order) pairs, order tracked in a parallel list.
    # The top of the stack is always the most recently doubled point.
    strategy = [(T1, T2)]
    orders = [e]

    for _ in range(e):
        while orders[-1] != 1:
            K1, K2 = strategy[-1]
            m = orders[-1] // 2
            K1 = OA.double_iter(K1, m)
            K2 = OA.double_iter(K2, m)
            strategy.append((K1, K2))
            orders.append(orders[-1] - m)

        # Pop the order-4 point and use it as the kernel of a 4-isogeny.
        K1, K2 = strategy.pop()
        orders.pop()
        OA, inv_dual_null = isogeny_codomain(K1, K2)

        # Push every remaining point through the isogeny.
        strategy = [[isogeny_eval(P, inv_dual_null) for P in K] for K in strategy]
        orders = [order - 1 for order in orders]

    return OA
