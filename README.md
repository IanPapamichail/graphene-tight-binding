![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)
![NumPy](https://img.shields.io/badge/NumPy-Scientific%20Computing-orange.svg)
![Status](https://img.shields.io/badge/Status-Active%20Development-brightgreen)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)
![Domain](https://img.shields.io/badge/Domain-Computational%20Physics-purple)
![Model](https://img.shields.io/badge/Model-Tight--Binding-red)
![Last Commit](https://img.shields.io/github/last-commit/IanPapamichail/graphene-tight-binding)
![Repo Size](https://img.shields.io/github/repo-size/IanPapamichail/graphene-tight-binding)

Graphene Tight-Binding Simulation (Hexagonal Fragment)

This project implements a tight-binding model of a graphene hexagonal fragment (6-site ring) using Python and NumPy.

Overview

Graphene’s electronic properties can be approximated using the tight-binding formalism, where electrons hop between nearest-neighbor carbon atoms.

The Hamiltonian is constructed as:

H[i,j] = -t for nearest neighbors
H[i,i] = εᵢ (on-site energy)

Eigenvalue decomposition provides the allowed energy spectrum and corresponding eigenstates.

Features
	•	Real-space graphene fragment modeling
	•	On-site energy support
	•	Nearest-neighbor hopping
	•	Eigenvalue computation
	•	Eigenstate visualization on lattice geometry

Physical Interpretation

For t = 1 and εᵢ = 0:

Energy spectrum:
[-2, -1, -1, 1, 1, 2]

The degeneracies arise from symmetry in the hexagonal lattice.

Future Extensions
	•	Hubbard interaction term (electron correlation)
	•	Fermionic operator formulation
	•	Jordan–Wigner mapping
	•	Qubit Hamiltonian construction
	•	Variational Quantum Eigensolver (VQE)

Technologies Used
	•	Python
	•	NumPy
	•	Matplotlib
