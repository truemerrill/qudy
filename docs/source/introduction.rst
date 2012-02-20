Introduction
============

``qudy`` is an open-source quantum dynamics package implemented in
Python.  Although many routines provided in the module are quite
general, the program is specifically designed to solve problems in
quantum control theory.  In recent years several control theoretic
methods have found wide application across physics.  In particular,
certain problems relevant to coherent spectroscopy (e.g. NMR, ESR,
laser spectroscopy), quantum computation, and in AMO physics may be
studied using quantum control theory.

``qudy`` is intended to be an aid for scientific research by
meeting the following design goals:

1. **Solve control theory problems:** The program should solve several
common problems in quantum control theory using a variety of
techniques, such as:
   
   * Trotter formulas
   * Time-dependent perturbation theory (Dyson series)
   * Magnus expansions
   * Markovian master equations
   * Open quantum systems
   
The program should also solve several families of optimal control
problems, including:
   
   * Optimal fidelity
   * Time optimal control
   * Constrained optimal control, for example constraints applied by
     experimental hardware
   * Optimal control problems of the Mayer, Lagrange and Bolza type

The program should also include several subroutines for studying the
dynamical Lie algebra and for establishing the controllability of the
dynamical system.

2. **Open source:** The program should have an nonrestrictive
open-source license, and encourage uses to actively participate in
code development.

3. **Modular:** The program design should be modular, and should take
advantage of the Python OOP paradigm.  The source code should be well
organized, clean, easy to read, and easily extendable.

4. **Distributable:** The program should be distributable to many
different architectures without modification.  All build dependencies
should also be open-source and distributable.  When possible, the
program should detect and make use of system resources, including
multiple cores and scratch space.
