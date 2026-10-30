%% Records of FS 4.4 and 4.5.

-record(circuit, {
    n_modes  = 0         :: non_neg_integer(),
    input    = vacuum    :: term(),                  % input_spec()
    hbar     = undefined :: undefined | number(),    % undefined: the configured default (FS 3.1)
    instrs   = []        :: [tuple()],               % instruction()
    params   = #{}       :: #{atom() => term()},     % default bindings (sym_expr)
    symbols  = #{}       :: #{atom() => real | positive | nonnegative | complex},
    meta     = #{}       :: map()
}).

-record(gaussian_state, {
    n_modes :: non_neg_integer(),
    hbar    = 2 :: number(),
    mean    :: term(),                               % packed_rvector() | sym_expr(), length 2m, xxpp
    cov     :: term()                                % packed_rmatrix() | sym_expr(), 2m x 2m, xxpp
}).

-record(fock_state, {
    n_modes :: non_neg_integer(),
    cutoff  :: non_neg_integer(),
    amps    :: term()                                % packed_cvector() | sym_expr(), length D (FS 3.4)
}).

-record(fock_dm, {
    n_modes :: non_neg_integer(),
    cutoff  :: non_neg_integer(),
    rho     :: term()                                % packed_cmatrix(), D x D
}).
