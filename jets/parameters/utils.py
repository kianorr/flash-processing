def fiuza_mfp(A_alpha, A_beta, Z_alpha, Z_beta, n_beta, v_flow, coulomb_log):
    n_beta /= 1e19
    v_flow /= 1e8
    numerator = 670 * A_alpha ** 2 * A_beta ** 2 * v_flow ** 4
    denominator = Z_alpha ** 2 * Z_beta ** 2 * (A_alpha + A_beta) ** 2 * n_beta * coulomb_log
    return numerator / denominator