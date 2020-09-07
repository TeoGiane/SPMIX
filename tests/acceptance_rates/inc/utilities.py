import numpy as __np__
import jax.numpy as __jnp__
import scipy.stats as __st__
import jax.scipy.stats as __jst__
from jax import jit as __jit__
import math as __op__
import autodiff as __ad__

# Data
# y = __np__.array(
#     [[-1.32143353, -0.51898207, 0.66863585, 0.99202206],
#      [-0.90519613, 0.09574643, -0.52083444, 0.16385932],
#      [1.41940461, 1.62953713, 0.87751597, 0.99753334],
#      [0.5519889, -1.39318992, 0.83221352, -0.35200125],
#      [-0.12552541, 0.73103625, 0.65526245, 0.31633918],
#      [0.56361516, -2.14527027, 0.74484321, 0.10014209],
#      [1.16161642, 0.02174823, 1.95924569, 0.05048404],
#      [0.45726673, -1.35790161, -0.50491895, -0.05047161],
#      [2.25243613, 2.19197456, 0.15546908, 0.43300146],
#      [-1.12732651, -1.49262252, 0.29936592, 0.16520678],
#      [0.67141161, -0.24103818, 0.86089341, -0.82464932],
#      [-0.08816647, 0.24937073, -1.48903397, 0.38687034],
#      [-1.23876655, 0.03168559, -0.60180304, 0.87221273],
#      [-0.0558993, -0.97508989, 0.51154181, 0.09699324],
#      [0.63613972, -0.53903145, 0.22885281, 0.13799762],
#      [-2.48422091, -1.22388189, 1.28810692, -1.39470603],
#      [-0.65191361, -0.66526812, 0.0123634, 1.54042766],
#      [-1.52636157, 0.07602121, 0.94604507, 0.15137594],
#      [0.95777342, -1.61036781, 0.54730304, -1.2029418],
#      [-1.09180203, -1.15764184, 1.9398773, 0.47318077],
#      [-0.08194096, 0.20736639, -0.42074774, 0.21242351],
#      [-0.432337, -0.2813933, 0.17934206, 1.92860458],
#      [-0.54598378, 2.38994893, 0.14226166, 0.41960939],
#      [0.61810484, -0.53631034, 1.02996508, 0.13660732],
#      [-0.52962221, 0.30847022, -0.29510575, 0.58173469],
#      [0.03081507, 0.13645585, -0.51989438, 0.3429146],
#      [0.01236229, 0.75303065, -1.23016396, 0.4018715],
#      [-0.17429395, 1.38806564, -1.07198918, 1.11207754],
#      [-0.89915718, -1.45906498, -0.50948914, 0.79603691],
#      [-1.4188268, 0.59514439, 0.85593439, 1.61872801],
#      [-0.85043408, -0.97563748, 0.31297229, -0.75339741],
#      [0.06544567, -0.66339345, 0.4549622, -0.76045047],
#      [0.49893906, -0.11097365, -1.99760268, 1.09979979],
#      [-0.49696856, -0.76949064, -0.71931509, -0.67897545],
#      [1.2876968, -0.95219247, 2.67363253, -0.54602438],
#      [-1.98866913, -1.07478678, 0.04918448, 1.82712027],
#      [0.41596888, -0.11332518, -1.57783639, 0.39609482],
#      [-1.79613997, -0.27333736, -0.38155095, -0.03961926],
#      [0.10329596, 1.17571236, 1.74861323, -0.33948211],
#      [0.29232079, 0.71129325, -1.25567901, -1.22551903],
#      [-1.52033959, 0.02919915, -0.01564822, -0.90030934],
#      [-0.81112123, 0.7146676, 0.64493158, 1.82376779],
#      [1.26317911, 2.54995835, -0.53625716, 0.08307352],
#      [-0.02008845, 0.4304167, -0.8701899, 1.63649928],
#      [-1.44623605, -1.90751326, -0.56870643, 0.48349046],
#      [-0.70271125, 0.0816952, 1.10873575, 0.22286869],
#      [1.24027146, -1.49552159, -0.45168953, 0.26004218],
#      [1.40013403, -0.37023397, 0.98589442, 0.26646571],
#      [-1.84670783, -0.10267385, 2.19575321, -0.73559475],
#      [0.15096012, -0.50679119, 0.47337793, -0.25901381],
#      [-2.63628005, 0.32094278, 1.58620186, 1.16375435],
#      [-1.99147369, -1.47707985, 0.98054915, 0.18129843],
#      [0.74549788, -1.02526845, 0.48561926, 0.96623012],
#      [0.5820761, -1.62963135, 0.40388384, 0.55448351],
#      [0.85868362, -0.36041486, 0.48118614, 0.98695599],
#      [-0.7154938, -0.00316957, -0.87437363, 0.81118121],
#      [1.96028303, 0.89207452, -0.45251474, 0.90035339],
#      [-0.68476444, 1.23781123, 0.75242273, 0.36847808],
#      [-1.26493479, -0.26543445, 0.62213317, -0.03104386],
#      [-0.3901302, -0.06320059, 1.2760787, 0.67179998],
#      [-0.68469318, 0.84231257, 1.10925926, -0.59721375],
#      [-0.80063364, 0.17746275, -0.25931016, 1.04055739],
#      [-0.315639, -1.49413778, -1.45095521, 1.3055353],
#      [1.00468791, 0.1121166, 0.49022903, 0.12839647],
#      [-1.69432244, -0.56161006, -1.63340365, -0.37358948],
#      [-0.60211922, -1.77913617, -0.5006201, 0.36960813],
#      [0.56673238, -0.8684426, -0.11516759, -0.56254088],
#      [-0.7928956, -1.6111134, -0.64748268, -0.44068287],
#      [-0.75379291, 0.43310738, 0.54810526, -2.23053716],
#      [-0.69121434, 0.44700934, 1.02705095, -0.557429],
#      [-1.2049535, -0.08999076, -0.7270101, -1.66804638],
#      [-1.17433708, -0.61626874, -1.32008912, -0.46573664],
#      [0.54785172, -1.31376003, -0.47457154, 0.56947721],
#      [-0.04305617, -2.20978121, -1.11381415, -0.50932077],
#      [1.80426537, 0.28506718, -0.56588916, 1.5078392],
#      [-0.83675791, -1.21766359, 0.3995987, 0.96728738],
#      [0.71163282, -3.15357423, -1.35626651, -0.15758902],
#      [0.57992903, -0.11050793, 2.81862017, 0.50705471],
#      [1.08554373, 0.26468976, 0.65571645, 0.55532055],
#      [1.04116905, 0.58257291, -1.01427902, 0.44804831],
#      [0.56418175, 0.34977266, 1.3734518, 0.50180224],
#      [-0.78434856, -0.07055733, -0.25003196, 0.84257752],
#      [-0.86829478, -1.8296362, 0.70495677, 1.66719186],
#      [-0.82532208, -0.96614542, 1.17861216, 1.02199428],
#      [0.21867554, -2.00757565, 0.25067124, 0.53315692],
#      [-1.33311388, -1.03078579, 0.9630271, 0.69477063],
#      [-0.65515633, 0.72098059, 1.67882848, -1.50787376],
#      [-0.12253014, -1.52211323, 1.04205519, 0.84597039],
#      [-3.0654022, -0.2915516, 0.60928342, 0.5219904],
#      [-1.01886448, -1.30669153, 0.4700761, -0.29178501],
#      [0.17080785, 1.28297049, 0.10577506, 1.43150115],
#      [-1.67268638, -0.60710573, 1.7800286, 1.04677723],
#      [0.04104859, 0.62886498, 0.97947866, 0.57515702],
#      [-1.89736055, -0.74058903, 1.31843233, -0.48945251],
#      [0.38912019, 0.02869802, -1.40379928, 0.56420767],
#      [-0.09977907, -1.18402051, 0.63898124, 0.88325489],
#      [0.3263557, 1.38730562, 0.07353949, 0.41653438],
#      [1.0391162, 0.64606349, -0.04503566, -0.67482545],
#      [-1.00350591, -1.65081151, -0.18870751, 1.08081531],
#      [-1.63610675, -0.93623312, -0.66301958, 1.01796523]]
# )
#
# w_tilde = __np__.array([[-2.04413099, 3.61422721, -3.36737531, 2.31553919,
#                          -1.88079731, 1.94087401, -0.73716732, 0.93167129]])
#
# tau = __np__.array(
#     [[3.11629214, 1.63714994],
#      [-0.30394855, 1.22396014],
#      [0.60742063, 1.72996584]]
# )
#
# h_dim: int = 4
#
# input_list = {
#     "data": y,
#     "w_tilde": w_tilde,
#     "tau": tau,
#     "h_dim": h_dim
# }
#
# a = 2
# b = 2
# mu0 = 0
# l = 0.1
# m_tilde = __np__.array([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]])
# S_tilde = __jnp__.array(
#     [[2.97736114, 0., 0., 2.01071063, 0., 0.,
#       1.45940086, 0., 0., 1.20955099, 0., 0.],
#      [0., 2.97736114, 0., 0., 2.01071063, 0.,
#       0., 1.45940086, 0., 0., 1.20955099, 0.],
#      [0., 0., 2.97736114, 0., 0., 2.01071063,
#       0., 0., 1.45940086, 0., 0., 1.20955099],
#      [2.01071063, 0., 0., 2.42605137, 0., 0.,
#       1.76086076, 0., 0., 1.45940086, 0., 0.],
#      [0., 2.01071063, 0., 0., 2.42605137, 0.,
#       0., 1.76086076, 0., 0., 1.45940086, 0.],
#      [0., 0., 2.01071063, 0., 0., 2.42605137,
#       0., 0., 1.76086076, 0., 0., 1.45940086],
#      [1.45940086, 0., 0., 1.76086076, 0., 0.,
#       2.42605137, 0., 0., 2.01071063, 0., 0.],
#      [0., 1.45940086, 0., 0., 1.76086076, 0.,
#       0., 2.42605137, 0., 0., 2.01071063, 0.],
#      [0., 0., 1.45940086, 0., 0., 1.76086076,
#       0., 0., 2.42605137, 0., 0., 2.01071063],
#      [1.20955099, 0., 0., 1.45940086, 0., 0.,
#       2.01071063, 0., 0., 2.97736114, 0., 0.],
#      [0., 1.20955099, 0., 0., 1.45940086, 0.,
#       0., 2.01071063, 0., 0., 2.97736114, 0.],
#      [0., 0., 1.20955099, 0., 0., 1.45940086,
#       0., 0., 2.01071063, 0., 0., 2.97736114]]
# )
#
# param_list = {
#     "a": a,
#     "b": b,
#     "mu0": mu0,
#     "lambda": l,
#     "m_tilde": m_tilde,
#     "S_tilde": S_tilde
# }
#
# opt_options = {
#   "max_iter": 1,
#   "tol": 1e-8,
#   "w_0": __np__.array([[0., 0., 0., 0.]]),
#   "tau_0": __np__.array([[0., 1.]])
# }


# Functions definitions
def __alr_row__(w: __np__.array):
    h = w.shape[0]
    return __jnp__.log(w[range(h - 1)] / w[h - 1])


def alr(w: __jnp__.array, track: bool = False):
    if w.ndim == 1:
        output = __alr_row__(w)
        return output if track is True else __np__.array(output)
    elif w.ndim == 2:
        for elem in range(w.shape[0]):
            if elem == 0:
                output = __alr_row__(w[elem, :])
            else:
                output = __jnp__.concatenate((output, __alr_row__(w[elem, :])))
        output = output.reshape((w.shape[0], w.shape[1] - 1))
        return output if track is True else __np__.array(output)
    else:
        return "Wrong number of dimensions in input."


def __inv_alr_row__(w_tilde: __jnp__.ndarray):
    w_sum = __jnp__.sum(__jnp__.exp(w_tilde))
    return __jnp__.append(__jnp__.exp(w_tilde) / (1. + w_sum), __jnp__.array([1. / (1. + w_sum)]))


def inv_alr(w_tilde: __jnp__.ndarray, track: bool = False):
    if w_tilde.ndim == 1:
        output = __inv_alr_row__(w_tilde)
        return output if track is True else __np__.array(output)
    elif w_tilde.ndim == 2:
        for elem in range(w_tilde.shape[0]):
            if elem == 0:
                output = __inv_alr_row__(w_tilde[elem, :])
            else:
                output = __jnp__.concatenate((output, __inv_alr_row__(w_tilde[elem, :])))
        output = output.reshape((w_tilde.shape[0], w_tilde.shape[1] + 1))
        return output if track is True else __np__.array(output)
    else:
        return "Wrong number of dimensions in input."


def spmix_loglikelihood(data: __np__.array, w_tilde: __np__.array, tau: __np__.array, h_dim: int,
                        param: dict, w_tilde_new: __np__.array = None, tau_new: __np__.array = None,
                        include_h: bool = True):
    # Eliciting problem dimensions
    n_dim, i_dim = data.shape
    # TODO: Estendere al caso di N_i diversi
    # Initialize output
    output = 0.
    # Extending w_tilde vector and kernel matrix
    # print("w_tilde\n", w_tilde)
    # print("w_tilde_new\n", w_tilde_new)
    w_tilde_ext = __np__.append(w_tilde, w_tilde_new, axis=0) if w_tilde_new is not None else w_tilde
    # print("w_tilde_ext\n", w_tilde_ext, "\n")
    tau_ext = __np__.append(tau, tau_new, axis=0) if tau_new is not None else tau
    # print("tau_ext\n", tau_ext, "\n")
    # Computing contribution from data
    for i in range(i_dim):
        for j in range(n_dim):
            tmp = 0.
            for h in range(h_dim):
                tmp += inv_alr(w_tilde_ext.reshape((i_dim, h_dim - 1)))[i, h] * \
                       __st__.norm.pdf(data[j, i], loc=tau_ext[h, 0], scale=__op__.sqrt(tau_ext[h, 1]))
            output += __op__.log(tmp)
    # Computing contribution from kernels
    tmp = 0.
    for h in range(h_dim):
        tmp += __st__.gamma.logpdf(1 / tau_ext[h, 1], param['a'], scale=1 / param['b']) + \
               __st__.norm.logpdf(tau_ext[h, 0], param['mu0'], __op__.sqrt(tau_ext[h, 1] / param['lambda']))
    output += tmp
    # Computing contribution from weights
    output += __st__.multivariate_normal.logpdf(w_tilde_ext, mean=param['m_tilde'], cov=param['S_tilde'])
    # Computing contribution from H (if required)
    if include_h is True:
        output += __st__.poisson.logpmf(h_dim, mu=1)
    # Return output
    return output


def __spmix_loglikelihood_jax__(data: __np__.array, w_tilde: __np__.array, tau: __np__.array, h_dim: int,
                                param: dict, w_tilde_new: __np__.array = None, tau_new: __np__.array = None,
                                include_h: bool = True):
    # Eliciting problem dimensions
    n_dim, i_dim = data.shape
    # TODO: Estendere al caso di N_i diversi

    # Initialize output
    output = 0.

    # Conversions
    data = __jnp__.array(data)
    w_tilde = __jnp__.array(w_tilde)
    tau = __jnp__.array(tau)
    w_tilde_new = __jnp__.array(w_tilde_new)
    tau_new = __jnp__.array(tau_new)

    # Extending w_tilde vector and kernel matrix
    w_tilde_ext = __jnp__.append(w_tilde, w_tilde_new, axis=1) if w_tilde_new is not None else w_tilde
    tau_ext = __jnp__.append(tau, tau_new, axis=0) if tau_new is not None else tau

    # Computing contribution from data
    w_mat = inv_alr(w_tilde_ext.reshape((i_dim, h_dim - 1)), track=True)
    for i in range(i_dim):
        for j in range(n_dim):
            tmp = 0.
            for h in range(h_dim):
                tmp += w_mat[i, h] * \
                       __jst__.norm.pdf(data[j, i], loc=tau_ext[h, 0], scale=__jnp__.sqrt(tau_ext[h, 1]))
            output += __jnp__.log(tmp)

    # Computing contribution from kernels
    tmp = 0.
    for h in range(h_dim):
        tmp += __jst__.gamma.logpdf(1 / tau_ext[h, 1], param['a'], scale=1 / param['b']) + \
               __jst__.norm.logpdf(tau_ext[h, 0], param['mu0'], __jnp__.sqrt(tau_ext[h, 1] / param['lambda']))
    output += tmp

    # Computing contribution from weights (not using jax since a bug occurs)
    output += __jst__.multivariate_normal.logpdf(w_tilde_ext, mean=param['m_tilde'], cov=param['S_tilde']).squeeze()

    # Computing contribution from H (if required)
    if include_h is True:
        output += __jst__.poisson.logpmf(h_dim, mu=1)

    # Return output
    return output

# def norets_loss_function(w_toadd, tau_toadd):
#     return -__spmix_loglikelihood_jax__(input_list['data'], input_list['w_tilde'], input_list['tau'],
#                                         input_list['h_dim'], param_list, w_toadd, tau_toadd, False)

# TODO: Fare il test per passi
def norets_optimal_proposal(input_list: dict, param_list: dict, opt_options: dict):
    # Defining the loss funciton and jitting

    def norets_loss_function(w_toadd, tau_toadd):
        return -__spmix_loglikelihood_jax__(input_list['data'], input_list['w_tilde'], input_list['tau'],
                                            input_list['h_dim'], param_list, w_toadd, tau_toadd, False)

    # Newton method for optimization
    w_input = __np__.array(opt_options['w_0'])
    tau_input = __np__.array(opt_options['tau_0'])
    for iter in range(opt_options['max_iter']):
        hess = __ad__.hess_recast(__ad__.hessian(norets_loss_function)(w_input, tau_input))
        grad = __ad__.grad_recast(__ad__.gradient(norets_loss_function)(w_input, tau_input)).transpose()
        h = __jnp__.dot(__np__.linalg.inv(hess), grad)
        w_h, tau_h = __np__.split(h.transpose().squeeze(), [input_list['h_dim']])
        w_h = __np__.array([w_h])
        tau_h = __np__.array([tau_h])
        w_input = w_input - w_h
        tau_input = tau_input - tau_h
        # Convergence check
        grad_norm = __np__.linalg.norm(__ad__.grad_recast(__ad__.gradient(norets_loss_function)(w_input, tau_input)))
        if grad_norm < opt_options['tol']:
            break

    # Return output
    output = {
        "post_mode": __np__.append(w_input, tau_input, axis=1),
        "post_variance": __ad__.hess_recast(__ad__.hessian(norets_loss_function)(w_input, tau_input)),
        "norm": grad_norm
    }
    return output

# Testing
# __spmix_loglikelihood_jax__(input_list['data'], input_list['w_tilde'], input_list['tau'],
#                             input_list['h_dim'], param_list, opt_options['w_0'], opt_options['tau_0'])


# tmp = norets_optimal_proposal(input_list, param_list, opt_options)
# print(tmp['post_mode'])
# print(tmp['post_variance'])
