
from abc import ABC, abstractmethod
import netket as nk
from netket.operator.spin import sigmax, sigmaz, identity
import numpy as np
from typing import List
from scipy import sparse

from netket.hilbert import AbstractHilbert
from netket.utils.types import DType

from netket.operator._local_operator.compile_helpers import _append_matrix, _append_matrix_sparse, max_nonzero_per_row

LocalOperator = nk.operator._local_operator.LocalOperator



class ConstraintLocalOperator(LocalOperator):
    def _setup(self, force: bool = False):
        """Analyze the operator strings and precompute arrays for get_conn inference"""
        if force or not self._initialized:
            data = pack_internals(
                self.hilbert,
                self._operators_dict,
                self.constant,
                self.dtype,
                self.mel_cutoff,
            )

            self._acting_on = data["acting_on"]
            self._acting_size = data["acting_size"]
            self._diag_mels = data["diag_mels"]
            self._mels = data["mels"]
            self._x_prime = data["x_prime"]
            self._n_conns = data["n_conns"]
            self._local_states = data["local_states"]
            self._basis = data["basis"]
            self._nonzero_diagonal = data["nonzero_diagonal"]
            self._max_conn_size = data["max_conn_size"]
            self._initialized = True


def pack_internals(
    hilbert: AbstractHilbert,
    operators_dict: dict,
    constant,
    dtype: DType,
    mel_cutoff: float,
):
    """
    Take the internal lazy representation of a local operator and returns the arrays
    needed for the numba implementation.

    This takes as input a dictionary with Tuples as keys, the `acting_on` and matrices as values.
    The keys represent the sites upon which the matrix acts.
    It is assumed that the integer in the tuples are sorted.

    Returns a dictionary with all the data fields
    """
    op_acting_on = list(operators_dict.keys())
    operators = list(operators_dict.values())
    n_operators = len(operators_dict)

    """Analyze the operator strings and precompute arrays for get_conn inference"""

    # how many sites each operator is acting on
    acting_size = np.array([len(aon) for aon in op_acting_on], dtype=np.intp)

    # compute the maximum number of off-diagonal nonzeros (over all rows) of each operator
    op_n_conns_offdiag = max_nonzero_per_row(operators, mel_cutoff)

    # Support empty LocalOperators such as the identity.
    if len(acting_size) > 0:
        # maximum number of sites any operator is acting on
        max_acting_on_sz = np.max(acting_size)

        # max local hilbert size of all sites acted on by any operator
        max_local_hilbert_size = max(
            [max(map(hilbert.size_at_index, aon)) for aon in op_acting_on]
        )
        # maximum size of any operator (maximum size of the matrix / prod of local hilbert spaces)
        max_op_size = max(map(lambda x: x.shape[0], operators))
        # maximum number of off-diagonal nonzeros of any operator
        max_op_size_offdiag = np.max(op_n_conns_offdiag)
    else:
        max_acting_on_sz = 0
        max_local_hilbert_size = 0
        max_op_size = 0
        max_op_size_offdiag = 0

    # matrix storing which sites each operator acts on, padded with -1
    acting_on = np.full((n_operators, max_acting_on_sz), -1, dtype=np.intp)
    for i, aon in enumerate(op_acting_on):
        acting_on[i][: len(aon)] = aon

    ###
    # allocate empty arrays which are filled below

    # array which will be storing the local states
    # of each site each operator is acting on
    local_states = np.full(
        (n_operators, max_acting_on_sz, max_local_hilbert_size), np.nan
    )

    # array storing the basis for each site each operator is acting on
    # The basis is an integer used to map (indices of) local states on sites
    # to (indices of) states in the space spanned by all sites the op is acting on
    # (it is simply the product of number of local states of the sites before)
    basis = np.full((n_operators, max_acting_on_sz), 0x7FFFFFFF, dtype=np.int64)

    diag_mels = np.full((n_operators, max_op_size), np.nan, dtype=dtype)

    mels = np.full(
        (n_operators, max_op_size, max_op_size_offdiag),
        np.nan,
        dtype=dtype,
    )
    # x_prime contains the local state after the operator has been applied
    # for the sites the operator is acting on, for each row
    # (each row also corresponds to a list of local states)
    x_prime = np.full(
        (n_operators, max_op_size, max_op_size_offdiag, max_acting_on_sz),
        -1,
        dtype=np.float64,
    )
    # store the number of off-diagonal nonzeros per row of each operator
    n_conns = np.full((n_operators, max_op_size), 0, dtype=np.intp)

    ###
    # iterate over all operators
    for i, (aon, op) in enumerate(operators_dict.items()):
        # how many sites this operator is acting on
        aon_size = len(aon)

        n_local_states_per_site = np.asarray([hilbert.size_at_index(i) for i in aon])

        ## add an operator to local_states
        for j, site in enumerate(aon):
            local_states[i, j, : hilbert.shape[site]] = np.asarray(
                hilbert.states_at_index(site)
            )

        # compute the basis of each site of this operator
        # i.e. the product of the number of local states of all sites before it
        ba = 1
        for s in range(aon_size):
            basis[i, s] = ba
            ba *= hilbert.shape[aon[aon_size - s - 1]]

        if sparse.issparse(op):
            if not isinstance(op, sparse.csr_matrix):
                op = op.tocsr()
            # Extract the sparse matrix representation to numpy arrays
            data = np.array(op.data, copy=False)
            indices = np.array(op.indices, copy=False)
            indptr = np.array(op.indptr, copy=False)

            _append_matrix_sparse(
                data,
                indices,
                indptr,
                aon_size,
                local_states[i],
                n_local_states_per_site,
                mel_cutoff,
                diag_mels[i],
                mels[i],
                x_prime[i],
                n_conns[i],
            )

        else:
            _append_matrix(
                op,
                aon_size,
                local_states[i],
                n_local_states_per_site,
                mel_cutoff,
                diag_mels[i],
                mels[i],
                x_prime[i],
                n_conns[i],
            )

    nonzero_diagonal = (
        np.any(np.abs(diag_mels) >= mel_cutoff) or np.abs(constant) >= mel_cutoff
    )

    max_conn_size = 1 if nonzero_diagonal else 0
    # estimate max_conn_size with the sum of the
    # maximum number of off-diagonal nonzeros of all operators
    max_conn_size = max_conn_size + np.sum(op_n_conns_offdiag)

    return {
        "acting_on": acting_on,
        "acting_size": acting_size,
        "diag_mels": diag_mels,
        "mels": mels,
        "x_prime": x_prime,
        "n_conns": n_conns,
        "local_states": local_states,
        "basis": basis,
        "nonzero_diagonal": nonzero_diagonal,
        "max_conn_size": max_conn_size,
    }