export const GRID   = 8;
export const NCELLS = GRID * GRID;

const MAX_HIST = 300;

export const state = {
  nRabbits: 0,
  nFoxes:   0,
  stepN:    0,
  history:  [],   // [{s, r, f}]
};

function randCell() { return Math.floor(Math.random() * NCELLS); }

/**
 * Advance one step. Returns the per-cell grid (post-interaction) for rendering.
 * @param {number} pDeath  Per-step probability that each fox dies.
 */
export function tick(pDeath) {
  // 1. Randomly scatter all agents onto the grid
  const grid = Array.from({ length: NCELLS }, () => ({ r: 0, f: 0, event: null }));
  for (let i = 0; i < state.nRabbits; i++) grid[randCell()].r++;
  for (let i = 0; i < state.nFoxes;   i++) grid[randCell()].f++;

  // 2. Per-cell interactions
  let totalR = 0, totalF = 0;
  for (const c of grid) {
    if (c.f > 0 && c.r > 0) {
      // Predation: each fox eats one rabbit (up to min(f,r)), spawns one new fox
      const eaten = Math.min(c.f, c.r);
      c.r -= eaten;
      c.f += eaten;
      c.event = 'predation';
    } else if (c.f === 0 && c.r >= 2) {
      // Breeding: every pair of rabbits produces one offspring
      c.r += Math.floor(c.r / 2);
      c.event = 'breeding';
    }
    totalR += c.r;
    totalF += c.f;
  }

  // 3. Stochastic fox death
  let survived = 0;
  for (let i = 0; i < totalF; i++) if (Math.random() >= pDeath) survived++;

  state.nRabbits = totalR;
  state.nFoxes   = survived;
  state.stepN++;

  state.history.push({ s: state.stepN, r: state.nRabbits, f: state.nFoxes });
  if (state.history.length > MAX_HIST) state.history.shift();

  return grid;
}

export function resetState(r0, f0) {
  state.nRabbits = r0;
  state.nFoxes   = f0;
  state.stepN    = 0;
  state.history  = [{ s: 0, r: r0, f: f0 }];
}
