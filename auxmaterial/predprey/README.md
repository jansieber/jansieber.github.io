# Predator–Prey Simulation

A browser-based agent simulation of rabbits and foxes on an 8×8 grid. Each step, all agents are randomly scattered across the board and interactions are resolved per cell:

- **Two rabbits** on the same cell → they breed (one new rabbit)
- **A rabbit and a fox** on the same cell → the fox eats the rabbit and reproduces (one new fox)
- **Each fox** has a fixed probability of dying every step (starvation)

## Files

```
index.html       — markup and styles
app.js           — controls, event listeners, main loop
simulation.js    — pure simulation logic
renderer.js      — board and chart rendering
```

## Running locally

The simulation uses ES modules, so it must be served over HTTP — opening `index.html` directly as a file will not work.

### Option 1: Python (no install required)

```bash
python3 -m http.server 8080
```

### Option 2: Node.js (`npx`)

```bash
npx serve .
```

### Option 3: Node.js (`http-server`)

```bash
npm install -g http-server
http-server -p 8080
```

Once the server is running, open your browser and go to:

```
http://localhost:8080
```

## Controls

| Control | Description |
|---|---|
| Rabbits | Initial rabbit population |
| Foxes | Initial fox population |
| Fox death / step | Probability (0–1) each fox dies per step |
| Speed | Milliseconds between steps |
| ▶ Start | Begin the simulation |
| ■ Stop | Pause |
| ↺ Reset | Reset to initial populations |
