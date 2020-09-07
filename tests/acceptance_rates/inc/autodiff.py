import numpy as __np__
from jax import grad as __grad__, jacfwd as __jacfwd__, jacrev as __jacrev__
from inspect import signature as __signature__
import warnings as __w__

# Avoid GPU UserWarning
__w__.filterwarnings("ignore")

def gradient(f):
    nargs = len(__signature__(f).parameters)
    return __grad__(f, argnums=range(nargs))

def grad_recast(grad_output):
    return __np__.append(*grad_output, axis=1)

def hessian(f):
    nargs = len(__signature__(f).parameters)
    return __jacfwd__(__jacrev__(f, argnums=range(nargs)), argnums=range(nargs))

def hess_recast(hess_output):
    for i in range(len(hess_output)):
        row_block = __np__.block([elem for elem in hess_output[i]])
        if i == 0:
            res = row_block
        else:
            res = __np__.concatenate([res, row_block], axis=1)
    return res.squeeze()

# def fun(x):
#     return x ** x + x ** 2.0
#
#
# fun_fast = __grad__(fun)
# print(fun_fast(2.0))
