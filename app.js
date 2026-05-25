// Graphene Tight-Binding Simulator Core Logic
// Author: Antigravity AI Code Assistant

// App State
const state = {
    t: 1.0,                 // Hopping parameter
    eps: [0, 0, 0, 0, 0, 0], // On-site energies for 6 atoms
    selectedStateIndex: 0,  // Currently active eigenstate (0-5)
    eigenstates: [],        // Output of diagonalization: array of { eigenvalue, eigenvector }
    animationEnabled: true, // Flag to animate quantum phase factor
    animationTime: 0       // Time tracker for phase animation
};

// Lattice definition matching structure.py & Jupyter Notebook
const N_SITES = 6;
const EDGES = [
    [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]
];

// Scale factor and coordinates for 400x400 SVG representation
// Cartesian x is mapped to SVG x, positive Cartesian y is mapped to negative SVG y (upwards)
const SCALE = 135;
const COORDS = {
    0: { x: 1.0 * SCALE, y: 0.0 },
    1: { x: 0.5 * SCALE, y: -Math.sqrt(3)/2 * SCALE },
    2: { x: -0.5 * SCALE, y: -Math.sqrt(3)/2 * SCALE },
    3: { x: -1.0 * SCALE, y: 0.0 },
    4: { x: -0.5 * SCALE, y: Math.sqrt(3)/2 * SCALE },
    5: { x: 0.5 * SCALE, y: Math.sqrt(3)/2 * SCALE }
};

// Jacobi Eigenvalue Algorithm for real symmetric 6x6 matrices
// Solves H|psi> = E|psi> exactly on the client-side
function diagonalize(matrix) {
    let A = matrix.map(row => [...row]); // Deep copy
    let V = Array.from({ length: N_SITES }, (_, i) => 
        Array.from({ length: N_SITES }, (_, j) => (i === j ? 1.0 : 0.0))
    ); // Identity matrix to accumulate eigenvectors

    const maxIterations = 150;
    const eps = 1e-15;

    for (let iter = 0; iter < maxIterations; iter++) {
        // 1) Find the largest off-diagonal element
        let maxVal = 0.0;
        let p = 0, q = 0;
        for (let i = 0; i < N_SITES - 1; i++) {
            for (let j = i + 1; j < N_SITES; j++) {
                if (Math.abs(A[i][j]) > maxVal) {
                    maxVal = Math.abs(A[i][j]);
                    p = i;
                    q = j;
                }
            }
        }

        // Converged! Off-diagonals are zero to machine precision
        if (maxVal < eps) break;

        // 2) Compute Jacobi rotation components
        let diff = A[q][q] - A[p][p];
        let t;
        if (Math.abs(A[p][q]) < Math.abs(diff) * 1e-15) {
            t = A[p][q] / diff;
        } else {
            let phi = diff / (2.0 * A[p][q]);
            t = 1.0 / (Math.abs(phi) + Math.sqrt(1.0 + phi * phi));
            if (phi < 0.0) t = -t;
        }

        let c = 1.0 / Math.sqrt(1.0 + t * t);
        let s = t * c;
        let tau = s / (1.0 + c);

        // 3) Update matrix A elements (Givens rotation transformation)
        let tempApq = A[p][q];
        A[p][q] = 0.0;
        A[p][p] -= t * tempApq;
        A[q][q] += t * tempApq;

        for (let i = 0; i < N_SITES; i++) {
            if (i !== p && i !== q) {
                let tempAip = A[i][p];
                let tempAiq = A[i][q];
                A[i][p] = c * tempAip - s * tempAiq;
                A[q][i] = A[i][q] = s * tempAip + c * tempAiq;
                A[p][i] = A[i][p]; // Force symmetry
            }
        }

        // 4) Accumulate eigenvectors in matrix V
        for (let i = 0; i < N_SITES; i++) {
            let tempVip = V[i][p];
            let tempViq = V[i][q];
            V[i][p] = c * tempVip - s * tempViq;
            V[i][q] = s * tempVip + c * tempViq;
        }
    }

    // Extract eigenvalues (diagonals) & pack eigenvectors
    let results = [];
    for (let i = 0; i < N_SITES; i++) {
        results.push({
            eigenvalue: A[i][i],
            eigenvector: V.map(row => row[i])
        });
    }

    // Sort states ascending by eigenvalue (Energy spectrum)
    results.sort((a, b) => a.eigenvalue - b.eigenvalue);

    // Enforce standardized gauge / sign phase (e.g. let the first non-zero element be positive)
    results.forEach(res => {
        let firstSignificant = res.eigenvector.find(val => Math.abs(val) > 1e-4);
        if (firstSignificant && firstSignificant < 0) {
            res.eigenvector = res.eigenvector.map(val => -val);
        }
    });

    return results;
}

// Build Hamiltonian matrix from current state parameters
function buildHamiltonian() {
    let H = Array.from({ length: N_SITES }, () => Array(N_SITES).fill(0.0));
    
    // Set diagonal elements (on-site energies)
    for (let i = 0; i < N_SITES; i++) {
        H[i][i] = state.eps[i];
    }
    
    // Set off-diagonal elements (nearest-neighbor hopping)
    EDGES.forEach(([i, j]) => {
        H[i][j] = -state.t;
        H[j][i] = -state.t;
    });
    
    return H;
}

// Solve simulation and update UI components
function runSimulation() {
    const H = buildHamiltonian();
    state.eigenstates = diagonalize(H);
    
    // Render and synchronize everything
    updateMatrixUI(H);
    updateSpectrumLadderUI();
    updateLatticeUI();
    updateStatsTableUI();
}

// Presets loader
function loadPreset(presetName) {
    // Reset buttons
    document.querySelectorAll('.presets-grid .btn').forEach(btn => btn.classList.remove('active'));
    
    if (event && event.target) {
        event.target.classList.add('active');
    }

    switch(presetName) {
        case 'ideal':
            state.t = 1.0;
            state.eps = [0, 0, 0, 0, 0, 0];
            break;
        case 'alternating':
            state.t = 1.0;
            state.eps = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0];
            break;
        case 'doped':
            state.t = 1.0;
            state.eps = [2.0, 0, 0, 0, 0, 0];
            break;
        case 'dimers':
            // Custom decoupling is represented by zero-hopping, but since we are keeping ring geometry constant,
            // we can simulate localized on-site isolation: e.g. very strong potential boundaries.
            // Alternatively, setting high potentials decoupling dimers:
            state.t = 1.0;
            state.eps = [2.0, 2.0, -2.0, -2.0, 0.0, 0.0];
            break;
    }
    
    // Update HTML sliders to match loaded preset
    document.getElementById('t-slider').value = state.t;
    document.getElementById('t-value').innerText = state.t.toFixed(1);
    
    for (let i = 0; i < N_SITES; i++) {
        document.getElementById(`eps${i}-slider`).value = state.eps[i];
        document.getElementById(`eps${i}-value`).innerText = state.eps[i].toFixed(1);
    }
    
    // Re-run simulator
    runSimulation();
}

// Sliders and interaction updates
function updateParam(type, value, index) {
    const val = parseFloat(value);
    
    if (type === 't') {
        state.t = val;
        document.getElementById('t-value').innerText = val.toFixed(1);
    } else if (type === 'eps') {
        state.eps[index] = val;
        document.getElementById(`eps${index}-value`).innerText = val.toFixed(1);
    }
    
    // Remove active state from preset buttons since settings are custom now
    document.querySelectorAll('.presets-grid .btn').forEach(btn => btn.classList.remove('active'));
    
    runSimulation();
}

// Select state on energy ladder click
function selectState(index) {
    state.selectedStateIndex = index;
    
    // Update active visual tags
    document.querySelectorAll('.ladder-state-bar').forEach((bar, idx) => {
        if (idx === index) {
            bar.classList.add('active');
        } else {
            bar.classList.remove('active');
        }
    });
    
    // Update wavefunction analytics
    updateLatticeUI();
    updateStatsTableUI();
}

// ----------------- UI RENDERERS -----------------

// Populate 6x6 Matrix Elements grid
function updateMatrixUI(H) {
    const container = document.getElementById('matrix-elements');
    container.innerHTML = '';
    
    for (let i = 0; i < N_SITES; i++) {
        for (let j = 0; j < N_SITES; j++) {
            const cell = document.createElement('div');
            cell.classList.add('matrix-cell');
            
            const value = H[i][j];
            cell.innerText = value.toFixed(1);
            
            if (i === j) {
                cell.classList.add('diagonal');
                if (Math.abs(value) > 0.01) {
                    cell.classList.add('non-zero');
                }
            } else if (Math.abs(value) > 0.01) {
                cell.classList.add('hopping');
            } else {
                cell.style.opacity = '0.25';
            }
            
            container.appendChild(cell);
        }
    }
}

// Render vertical spectrum ladder
function updateSpectrumLadderUI() {
    const container = document.getElementById('states-ladder');
    container.innerHTML = '';
    
    const energies = state.eigenstates.map(s => s.eigenvalue);
    const minE = Math.min(...energies, -2.2);
    const maxE = Math.max(...energies, 2.2);
    const span = maxE - minE;
    
    // Check for double degeneracies to offset bars horizontally (so they don't overlap)
    const offsets = Array(N_SITES).fill(0);
    const tolerance = 0.02; // grouping degenerate states closer than 0.02t
    
    for (let i = 0; i < N_SITES; i++) {
        let degenerateCount = 0;
        let degenerateIndex = 0;
        
        for (let j = 0; j < N_SITES; j++) {
            if (Math.abs(energies[i] - energies[j]) < tolerance) {
                degenerateCount++;
                if (j < i) degenerateIndex++;
            }
        }
        
        if (degenerateCount > 1) {
            // Apply horizontal offset for degenerate splitting
            offsets[i] = (degenerateIndex - (degenerateCount - 1) / 2) * 55;
        }
    }

    state.eigenstates.forEach((est, idx) => {
        const val = est.eigenvalue;
        // Map eigenvalue position to vertical percentage (20% to 80%)
        const yPct = ((val - minE) / span) * 75 + 12.5; 
        
        const stateBar = document.createElement('div');
        stateBar.classList.add('ladder-state-bar');
        if (idx === state.selectedStateIndex) {
            stateBar.classList.add('active');
        }
        
        stateBar.style.bottom = `${yPct}%`;
        
        // Offset left & right limits to shift degenerate nodes
        if (offsets[idx] !== 0) {
            stateBar.style.left = `calc(15% + ${offsets[idx]}px)`;
            stateBar.style.right = `calc(15% - ${offsets[idx]}px)`;
        }
        
        stateBar.onclick = () => selectState(idx);
        
        stateBar.innerHTML = `
            <span class="ladder-state-badge">E<sub>${idx}</sub></span>
            <span class="ladder-state-energy-label">${val.toFixed(3)} eV</span>
        `;
        
        container.appendChild(stateBar);
    });
}

// Render dynamic elements inside Real-Space lattice
function updateLatticeUI(t = 0) {
    const selectedState = state.eigenstates[state.selectedStateIndex];
    if (!selectedState) return;
    
    const energy = selectedState.eigenvalue;
    const psi = selectedState.eigenvector;
    
    // Draw edges / bonds
    const bondsGroup = document.getElementById('bonds-group');
    bondsGroup.innerHTML = '';
    
    EDGES.forEach(([i, j]) => {
        const start = COORDS[i];
        const end = COORDS[j];
        
        const line = document.createElementNS("http://www.w3.org/2000/svg", "line");
        line.setAttribute("x1", start.x);
        line.setAttribute("y1", start.y);
        line.setAttribute("x2", end.x);
        line.setAttribute("y2", end.y);
        line.classList.add("lattice-bond");
        
        // Glow effect for bonds based on parameter t
        const opacity = Math.min(Math.abs(state.t) * 0.4 + 0.1, 0.9);
        const color = state.t >= 0 ? "var(--accent-teal)" : "var(--accent-purple)";
        line.setAttribute("stroke", color);
        line.setAttribute("stroke-opacity", opacity);
        
        bondsGroup.appendChild(line);
    });
    
    // Draw nodes / atoms
    const atomsGroup = document.getElementById('atoms-group');
    atomsGroup.innerHTML = '';
    
    for (let i = 0; i < N_SITES; i++) {
        const coord = COORDS[i];
        const amp = psi[i];
        const prob = Math.pow(amp, 2);
        
        // Base size represents standard atom radius (e.g. 15px). Max expansion scales with probability.
        // If animated, we add an oscillatory factor proportional to phase and amplitude:
        let dynamicRadius = 16 + prob * 35;
        let pulseFactor = 1.0;
        
        if (state.animationEnabled) {
            // Quantum time-dependent phase factor: exp(-i E t) -> Re = cos(E * t)
            // Frequency of oscillation corresponds to energy level
            const freq = energy * 1.5;
            pulseFactor = Math.cos(freq * t);
            
            // Modulate node radius: expands when dynamic amplitude goes positive, contracts on negative phase
            dynamicRadius += amp * pulseFactor * 16;
            // Prevent zero or negative radius values
            if (dynamicRadius < 4) dynamicRadius = 4;
        }
        
        // Glow backdrop group
        const glow = document.createElementNS("http://www.w3.org/2000/svg", "circle");
        glow.setAttribute("cx", coord.x);
        glow.setAttribute("cy", coord.y);
        glow.setAttribute("r", dynamicRadius * 1.8);
        glow.setAttribute("fill", "url(#node-glow)");
        
        const colorVar = amp * pulseFactor >= 0 ? "var(--accent-teal)" : "var(--accent-gold)";
        glow.style.setProperty("--glow-color", colorVar);
        atomsGroup.appendChild(glow);
        
        // Solid core circle
        const circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
        circle.setAttribute("cx", coord.x);
        circle.setAttribute("cy", coord.y);
        circle.setAttribute("r", dynamicRadius);
        circle.classList.add("lattice-node");
        
        const nodeColor = amp * pulseFactor >= 0 ? "var(--accent-teal)" : "var(--accent-gold)";
        circle.setAttribute("fill", nodeColor);
        circle.setAttribute("fill-opacity", "0.85");
        
        // Border color highlights individual atom potentials
        const borderGlow = state.eps[i] > 0.05 ? "var(--node-0)" : "rgba(255,255,255,0.4)";
        circle.setAttribute("stroke", borderGlow);
        circle.setAttribute("stroke-width", "2");
        
        // Clicking atom increases its potential directly (+0.5 eV per click, loops at +3.0)
        circle.onclick = () => {
            let nextEps = state.eps[i] + 0.5;
            if (nextEps > 3.0) nextEps = -3.0;
            updateParam('eps', nextEps, i);
            document.getElementById(`eps${i}-slider`).value = nextEps;
        };
        
        atomsGroup.appendChild(circle);
        
        // Node indexes text
        const text = document.createElementNS("http://www.w3.org/2000/svg", "text");
        text.setAttribute("x", coord.x);
        text.setAttribute("y", coord.y + 5);
        text.setAttribute("text-anchor", "middle");
        text.classList.add("node-text");
        text.textContent = i;
        
        atomsGroup.appendChild(text);
    }
}

// Update Analysis Stats Table
function updateStatsTableUI() {
    const selectedState = state.eigenstates[state.selectedStateIndex];
    if (!selectedState) return;
    
    // Update summary headers
    document.getElementById('state-index').innerText = `State ${state.selectedStateIndex}`;
    document.getElementById('state-energy-val').innerText = selectedState.eigenvalue.toFixed(4);
    
    const tbody = document.getElementById('eigenstate-table-body');
    tbody.innerHTML = '';
    
    const psi = selectedState.eigenvector;
    
    for (let i = 0; i < N_SITES; i++) {
        const tr = document.createElement('tr');
        const amp = psi[i];
        const prob = Math.pow(amp, 2);
        
        // Phase class
        let phaseTagHTML = '';
        if (amp > 0.005) {
            phaseTagHTML = `<span class="phase-tag positive">+</span>`;
        } else if (amp < -0.005) {
            phaseTagHTML = `<span class="phase-tag negative">-</span>`;
        } else {
            phaseTagHTML = `<span class="phase-tag zero">0</span>`;
        }
        
        tr.innerHTML = `
            <td><strong style="color: var(--node-${i})">Atom ${i}</strong></td>
            <td class="energy-cell-value">${state.eps[i].toFixed(1)} eV</td>
            <td class="amp-value" style="color: ${amp >= 0 ? 'var(--accent-teal)' : 'var(--accent-gold)'}">${amp.toFixed(6)}</td>
            <td class="prob-value">${prob.toFixed(6)}</td>
            <td>${phaseTagHTML}</td>
        `;
        
        tbody.appendChild(tr);
    }
}

// ----------------- ANIMATION CONTROL -----------------

function toggleAnimation(enabled) {
    state.animationEnabled = enabled;
}

// Phase Animation Loop
let lastTime = 0;
function quantumPhaseLoop(timestamp) {
    if (!lastTime) lastTime = timestamp;
    const elapsed = (timestamp - lastTime) * 0.001; // elapsed seconds
    
    if (state.animationEnabled) {
        state.animationTime += elapsed;
        updateLatticeUI(state.animationTime);
    }
    
    lastTime = timestamp;
    requestAnimationFrame(quantumPhaseLoop);
}

// Initialize Application
window.addEventListener('DOMContentLoaded', () => {
    runSimulation();
    requestAnimationFrame(quantumPhaseLoop);
});
