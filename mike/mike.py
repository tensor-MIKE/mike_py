from sage.all import GF, EllipticCurve, proof
import os
from hashlib import sha3_256

from dim_one.torsion_basis import torsion_basis
from dim_one.four_isogeny_chain import four_isogeny_chain
from gluing.mike_gluing import gluing_isogeny
from dim_four.isogeny_chain import dim_four_isogeny_chain
from utilities.aes256_ctr_drbg import AES256_CTR_DRBG


# Parameters for MIKE for prime p = c * 2**f - 1 with isogeny chains of length e
PARAMS_I = {
    "f": 308,
    "c": 633,
    "e": 228,
}

PARAMS_III = {
    "f": 474,
    "c": 593,
    "e": 340,
}

PARAMS_V = {
    "f": 628,
    "c": 317,
    "e": 452,
}


class Mike:
    def __init__(self, parameters, seed=None):
        """
        Given the three constants f, c, e precompute all public parameters
        of the MIKE key exchange
        """
        proof.all(False)

        self.f = parameters["f"]
        self.c = parameters["c"]
        self.e = parameters["e"]

        p = self.c * 2**self.f - 1
        self.Fp = GF(p)
        self.Fp2 = GF((p, 2), name="i", modulus=[1, 0, 1])
        self.encoding_length = (int(p).bit_length() + 7) // 8

        self._precompute_starting_curve_data()

        if seed is None:
            self.random_bytes = os.urandom
        else:
            drbg = AES256_CTR_DRBG(seed)
            self.random_bytes = drbg.random_bytes

    # ============================================================
    # Public API
    # ============================================================

    def keygen(self) -> tuple[bytes]:
        sk = self._secret_key()
        pk = self._public_key(sk)

        pk_bytes = self._encode_pk(pk)
        sk_bytes = self._encode_sk(sk)
        return pk_bytes, sk_bytes

    def shared_secret(self, sk, pk_other) -> bytes:
        EB = self._decode_pk(pk_other)
        x = self._decode_sk(sk)

        # Compute the gluing and chain kernel using ec arithmetic
        gluing_kernel, chain_kernel = self._generate_kernel_data(EB, x)

        # Compute the image of the chain kernel through the gluing isogeny
        (domain, chain_kernel) = gluing_isogeny(EB, gluing_kernel, chain_kernel)

        # Compute the codomain of the dim 4 isogeny chain after gluing
        codomain = dim_four_isogeny_chain(domain, chain_kernel, self.e - 3)

        absolute_invariants = codomain.absolute_invariants()
        return self._hash_moduli_invariants(absolute_invariants)

    # ============================================================
    # Precomputation Methods
    # ============================================================

    def _precompute_E0_basis(self):
        """
        Compute the basis <P0, Q0> = E0[2^f] such that [2]P has Fp rational
        coordinates and [2]Q = (-x, i*y)
        """
        # The curve E0 is used for key generation, but we only store the
        # isomorphic curve in the Montgomery model for efficient isogeny
        # computation during the public key computation.
        E0 = EllipticCurve(self.Fp2, [-1, 0])

        # Compute E0[2^(f - 1)] = <P, Q>
        x = self.Fp.zero()
        while True:
            # x must be an element of Fp such that (1 - x) is a NQR
            x += self.Fp.one()
            if (1 - x).is_square():
                continue
            # P must be a point with Fp rational y in E[2^f-1]
            P = self.c * E0.lift_x(x)
            if P.y() not in self.Fp:
                continue
            T = 2 ** (self.f - 2) * P
            if not T:
                continue
            break
        Q = E0(-P.x(), self.Fp2.gen() * P.y())

        # Now we lift the points to have full even torsion 2^f, which allows
        # us to to use the four torsion above the kernel of order 2^e providing
        # e <= f - 2
        P0 = min([X for X in P.division_points(2) if X.order() == 2**self.f])
        Q0 = min([X for X in Q.division_points(2) if X.order() == 2**self.f])
        P0 = 2 ** (self.f - self.e - 2) * P0
        Q0 = 2 ** (self.f - self.e - 2) * Q0

        return P0, Q0

    def _precompute_starting_curve_data(self):
        """
        For key generation, we want to work with the following Elliptic curve:

        E0 : y^2 = x(x + 1)(x - 1), <P, Q> = E[2^(f - 1)]

        Where in particular, we pick P above (1, 0) and Q above (-1, 0). This is
        what allows for a clear description for secret key generation.

        However, for efficiency, we want to work with the Montgomery model, so
        once we have E0 and the basis, we compute another isomorphism to the curve
        E0_M : y^2 = x^3 + Ax^2 + x and push the basis through this isomorphism.

        Additionally, we do not want the kernel ker(phi) = P0 + [x]Q0 to be above
        the point (0 : 0), and so during the isomorphism we make the specific
        choice to map Q above (0 : 0).

        To do this, we pick alpha, a root of  x(x + 1)(x - 1) to be 0, 1, or -1.
        We then compute s = 1 / (3*alpha^2 - 1) such that
        E0_M : B * y^2 = x^3 + Ax^2 + x; A = 3 alpha s, B = s
        (x, y) -> (s(x - alpha), sy)

        We see that to ensure that Q is mapped above (0, 0) we must select
        alpha = 1. We also want B = 1, so we need to have another factor of
        u = s.sqrt() when we compute y.
        """
        P0, Q0 = self._precompute_E0_basis()

        alpha = 1  # Pick this root to send Q to be above zero

        s = ~self.Fp2(3 * alpha**2 - 1).sqrt()
        u = s.sqrt()

        # E0_M : y^2 = x^3 + Ax^2 + x
        A = -3 * alpha * s
        self.E0 = EllipticCurve(self.Fp2, [0, A, 0, 1, 0])

        def iso(P):
            x, y = P.xy()
            x_new = (x + alpha) * s
            y_new = y * (s * u)  # Ensure B = 1
            return self.E0(x_new, y_new)

        # Compute the new basis from the special E0 basis
        (self.P0, self.Q0) = [iso(P) for P in (P0, Q0)]

    # ============================================================
    # Utility Methods
    # ============================================================

    def _lexicographically_largest(self, Ai):
        """
        Rust encodes a Fp2 value to a || b, then interprets this as
        a single little endian integer. This is the same as sorting
        by (A_im, A_re) in python.
        """
        b = []
        for A in Ai:
            b.append([A[1], A[0]])
        b_max = max(b)
        return self.Fp2([b_max[1], b_max[0]])

    def _normalize_montgomery(self, four_torsion):
        """
        Given a point P of order four on E, select the largest of the 12 possible
        public keys which can be computed from the various symmetries of the MIKE
        key exchange
        """
        i = self.Fp2.gen()
        assert i**2 == -1

        a0, b0 = four_torsion
        a1, b1 = a0 + b0, a0 - b0
        a2, b2 = a1 * i + b1, a1 + i * b1

        A0 = -2 * (a0**4 + b0**4) / (a0**4 - b0**4)
        A1 = -2 * (a1**4 + b1**4) / (a1**4 - b1**4)
        A2 = -2 * (a2**4 + b2**4) / (a2**4 - b2**4)
        A3 = -A0
        A4 = -A1
        A5 = -A2
        A6, A7, A8, A9, A10, A11 = (x.conjugate() for x in (A0, A1, A2, A3, A4, A5))

        A = self._lexicographically_largest(
            [A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11]
        )
        return EllipticCurve(self.Fp2, [0, A, 0, 1, 0])

    # ============================================================
    # Keygen Methods
    # ============================================================

    def _secret_key(self):
        # Sample N bytes of randomness for the scalar
        # and represent it as an integer (little endian)
        N = (self.e + 7) // 8
        x_bytes = self.random_bytes(N)
        x = int.from_bytes(x_bytes, "little")

        # Ensure scalar % 8 = 3
        x = ((x >> 3) << 3) | 3

        # Ensure scalar has at most e bits of data
        x %= 2**self.e

        return x

    def _public_key(self, sk):
        # Compute the point K such that [4]K = ker(phi) using the
        # precomputed torsion E0[2^(e + 2)] = <P0, Q0>
        k = self.P0 + sk * self.Q0

        # Push k through phi : E0 -> EA to recover a point of order
        # four on EA. NOTE: this can be optimised by using a bespoke
        # isogeny chain
        four_torsion = four_isogeny_chain(self.E0, k, self.e)

        # Normalize the codomain to the lexicographically largest of
        # the 6 possible isomorphic Montgomery coefficients A and their
        # complex conjugates
        return self._normalize_montgomery(four_torsion)

    def _encode_sk(self, sk) -> bytes:
        return int(sk).to_bytes((self.e + 7) // 8, "little")

    def _encode_pk(self, pk) -> bytes:
        A = pk.a2()
        A_re = int(A[0]).to_bytes(self.encoding_length, "little")
        A_im = int(A[1]).to_bytes(self.encoding_length, "little")
        return A_re + A_im

    # ============================================================
    # Shared Secret Methods
    # ============================================================

    def _decode_pk(self, pk):
        N = self.encoding_length
        A_re = int.from_bytes(pk[:N], "little")
        A_im = int.from_bytes(pk[N:], "little")
        A = self.Fp2([A_re, A_im])
        return EllipticCurve(self.Fp2, [0, A, 0, 1, 0])

    def _decode_sk(self, sk):
        return int.from_bytes(sk, "little")

    def _generate_kernel_data(self, EB, x):
        s = (x - 1) // 2

        # <P, Q> = E[2^(e + 2)]
        P, Q = torsion_basis(EB, self.e + 2)

        sP = s * P
        sQ = s * Q
        chain_kernel = [[sQ + Q, sQ], [-sP - P, -sP]]

        P_32 = 2 ** (self.e - 3) * P
        Q_32 = 2 ** (self.e - 3) * Q
        sP_32 = 2 ** (self.e - 3) * sP
        sQ_32 = 2 ** (self.e - 3) * sQ
        gluing_basis = [
            [sQ_32 + Q_32, sQ_32],
            [-sP_32 - P_32, -sP_32],
            [sQ_32 + Q_32 - sP_32 - P_32, sQ_32 - sP_32],
        ]

        P_16 = 2 * P_32
        Q_16 = 2 * Q_32
        diag_basis = [Q_16, P_16]

        P_8 = 2 * P_16
        Q_8 = 2 * Q_16

        scholten_basis = [P_8, Q_8]

        return (scholten_basis, diag_basis, gluing_basis), chain_kernel

    def _hash_moduli_invariants(self, invariants):
        hasher = sha3_256()
        for a in invariants:
            hasher.update(int(a).to_bytes(self.encoding_length, "little"))
        return hasher.digest()
