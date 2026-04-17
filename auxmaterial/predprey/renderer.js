import { NCELLS } from './simulation.js';

// ── board ────────────────────────────────────────────────────────────────────

function cellBase(i) {
  return ((i >> 3) + (i & 7)) % 2 === 0 ? 'empty' : 'alt';
}

export function buildBoard(boardEl) {
  const cellEls = [];
  for (let i = 0; i < NCELLS; i++) {
    const d = document.createElement('div');
    d.className = 'cell ' + cellBase(i);
    boardEl.appendChild(d);
    cellEls.push(d);
  }
  return cellEls;
}

export function renderBoard(cellEls, grid) {
  for (let i = 0; i < NCELLS; i++) {
    const el = cellEls[i];

    if (!grid) {
      el.className = 'cell ' + cellBase(i);
      el.innerHTML = '';
      continue;
    }

    const { r, f, event } = grid[i];

    if (r === 0 && f === 0) {
      el.className = 'cell ' + cellBase(i);
      el.innerHTML = '';
      continue;
    }

    let cls, icon, cnt;
    if (event === 'predation') {
      cls  = 'predation';
      icon = f > 0 && r > 0 ? '🦊🐇' : (f > 0 ? '🦊' : '🐇');
      cnt  = `${f}🦊 ${r > 0 ? r + '🐇' : ''}`.trim();
    } else if (event === 'breeding') {
      cls  = 'breeding';
      icon = '🐇';
      cnt  = `×${r} ✦`;
    } else if (f > 0) {
      cls  = 'fox';
      icon = '🦊';
      cnt  = f > 1 ? `×${f}` : '';
    } else {
      cls  = 'rabbit';
      icon = '🐇';
      cnt  = r > 1 ? `×${r}` : '';
    }

    el.className = `cell ${cls}`;
    el.innerHTML = `<span class="icon">${icon}</span><span class="cnt">${cnt}</span>`;
  }
}

export function renderStats(stepN, nRabbits, nFoxes, els) {
  els.step.textContent    = stepN;
  els.rabbits.textContent = nRabbits;
  els.foxes.textContent   = nFoxes;
}

// ── chart ────────────────────────────────────────────────────────────────────

export function resizeChart(chartEl) {
  chartEl.width = chartEl.parentElement.clientWidth - 32;
}

export function drawChart(history, chartEl) {
  const ctx = chartEl.getContext('2d');
  const W = chartEl.width, H = chartEl.height;
  ctx.clearRect(0, 0, W, H);

  ctx.fillStyle = '#252535';
  ctx.fillRect(0, 0, W, H);

  if (history.length < 2) return;

  const pad    = { l: 44, r: 14, t: 28, b: 28 };
  const w      = W - pad.l - pad.r;
  const h      = H - pad.t - pad.b;
  const maxPop = Math.max(10, ...history.map(e => Math.max(e.r, e.f)));

  // Grid lines + y-axis labels
  ctx.strokeStyle = '#3a3a50';
  ctx.lineWidth   = 1;
  ctx.font        = '10px system-ui';
  ctx.textAlign   = 'right';
  ctx.fillStyle   = '#666688';
  const ticks = 4;
  for (let i = 0; i <= ticks; i++) {
    const y = pad.t + (h * i / ticks);
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(W - pad.r, y); ctx.stroke();
    ctx.fillText(Math.round(maxPop * (ticks - i) / ticks), pad.l - 4, y + 3);
  }

  function plotLine(values, color) {
    if (values.length < 2) return;
    ctx.strokeStyle = color;
    ctx.lineWidth   = 2;
    ctx.lineJoin    = 'round';
    ctx.beginPath();
    values.forEach((v, i) => {
      const x = pad.l + (i / (values.length - 1)) * w;
      const y = pad.t + h - Math.min(v / maxPop, 1) * h;
      i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
    });
    ctx.stroke();
  }

  plotLine(history.map(e => e.r), '#a6e3a1');
  plotLine(history.map(e => e.f), '#fab387');

  // Legend
  ctx.font      = 'bold 11px system-ui';
  ctx.textAlign = 'left';
  ctx.fillStyle = '#a6e3a1'; ctx.fillText('● Rabbits', pad.l + 6,  pad.t - 8);
  ctx.fillStyle = '#fab387'; ctx.fillText('● Foxes',   pad.l + 90, pad.t - 8);

  // Step counter
  ctx.fillStyle = '#555577';
  ctx.textAlign = 'right';
  ctx.font      = '10px system-ui';
  ctx.fillText(`step ${history[history.length - 1].s}`, W - pad.r, pad.t - 8);
}
