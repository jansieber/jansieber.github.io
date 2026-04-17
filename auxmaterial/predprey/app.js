import { state, tick, resetState } from './simulation.js';
import { buildBoard, renderBoard, renderStats, resizeChart, drawChart } from './renderer.js';

// ── DOM refs ─────────────────────────────────────────────────────────────────

const boardEl  = document.getElementById('board');
const chartEl  = document.getElementById('chart');
const statusEl = document.getElementById('status');
const speedEl  = document.getElementById('speed');
const speedLbl = document.getElementById('speed-label');

const statEls = {
  step:    document.getElementById('v-step'),
  rabbits: document.getElementById('v-rabbits'),
  foxes:   document.getElementById('v-foxes'),
};

// ── init board ───────────────────────────────────────────────────────────────

const cellEls = buildBoard(boardEl);

resizeChart(chartEl);
window.addEventListener('resize', () => { resizeChart(chartEl); drawChart(state.history, chartEl); });

// ── params ───────────────────────────────────────────────────────────────────

function getParams() {
  return {
    r0:     Math.max(1, parseInt(document.getElementById('init-rabbits').value) || 30),
    f0:     Math.max(1, parseInt(document.getElementById('init-foxes').value)   || 8),
    pDeath: Math.min(1, Math.max(0, parseFloat(document.getElementById('fox-death').value) || 0.05)),
    speed:  parseInt(speedEl.value) || 600,
  };
}

// ── simulation loop ───────────────────────────────────────────────────────────

let running = false, timerId = null;

function step() {
  const { pDeath } = getParams();
  const grid = tick(pDeath);

  renderBoard(cellEls, grid);
  renderStats(state.stepN, state.nRabbits, state.nFoxes, statEls);
  drawChart(state.history, chartEl);

  if (state.nRabbits === 0 && state.nFoxes === 0) {
    stop(); statusEl.textContent = 'All agents extinct.';
  } else if (state.nFoxes === 0) {
    stop(); statusEl.textContent = 'Foxes extinct — rabbits will grow unchecked.';
  } else if (state.nRabbits === 0) {
    statusEl.textContent = 'Rabbits extinct — foxes starving…';
  }
}

function start() {
  if (running) return;
  if (state.nRabbits === 0 && state.nFoxes === 0) reset();
  running = true;
  timerId = setInterval(step, getParams().speed);
  statusEl.textContent = 'Running…';
}

function stop() {
  running = false;
  if (timerId) { clearInterval(timerId); timerId = null; }
}

function reset() {
  stop();
  const { r0, f0 } = getParams();
  resetState(r0, f0);
  renderBoard(cellEls, null);
  renderStats(state.stepN, state.nRabbits, state.nFoxes, statEls);
  drawChart(state.history, chartEl);
  statusEl.textContent = 'Ready — press Start.';
}

// ── event listeners ───────────────────────────────────────────────────────────

speedEl.addEventListener('input', () => {
  speedLbl.textContent = speedEl.value + ' ms';
  if (running) { stop(); start(); }
});

document.getElementById('btn-start').addEventListener('click', start);
document.getElementById('btn-stop').addEventListener('click', () => { stop(); statusEl.textContent = 'Paused.'; });
document.getElementById('btn-reset').addEventListener('click', reset);

// ── boot ──────────────────────────────────────────────────────────────────────

reset();
