NEURON {
    SUFFIX decay_inj
    NONSPECIFIC_CURRENT i
}

PARAMETER {
    iinj0 = 0
    tau = 1
}

ASSIGNED {
}

STATE {
    iinj
}

INITIAL {
    iinj = iinj0
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    i = iinj
}

DERIVATIVE state {
    iinj' = -iinj/tau
}
