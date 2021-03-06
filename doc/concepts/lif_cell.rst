.. _lifcell:

LIF cells
===========

The description of a LIF cell is used to control the leaky integrate-and-fire dynamics:

    * Resting potential.
    * Reset potential.
    * Initial value of membrane potential.
    * Membrane potential decaying constant.
    * Membrane capacitance.
    * Firing threshold.
    * Refractory period.

The morphology of a LIF cell is automatically modelled as a single :term:`compartment <control volume>`; each cell has one built-in
**source** and one built-in **target** which do not need to be explicitly added in the cell description.
LIF cells do not support adding additional **sources** or **targets** to the description. They do not support
**gap junctions**. They do not support adding density or point mechanisms.

API
---

* :ref:`Python <pylifcell>`
* :ref:`C++ <cpplifcell>`
