# Phase Portrait Sketcher

A browser tool for sketching qualitative phase portraits of planar flows, the kind drawn in a dynamical systems course. You draw trajectories, equilibria, periodic orbits and inflow/outflow boundaries by hand, as Bézier curves. The tool keeps the picture topologically consistent: trajectories must not cross, attachments must respect the direction of the flow, and so on.

* **One file, no installation:** open `index.html` in Chrome, Edge, Firefox or Safari. It works offline from disk.
* **Native format:** `*.phase.json`, a plain JSON description of the drawing.
* **Export:** SVG, PNG (3×), TikZ (`.tex`), and PDF via the browser's print dialog.

## Drawing tools

| Tool | Key | What it does |
|---|---|---|
| Select / edit | `V` | Click or drag a rectangle to select; `Shift` adds to the selection. Drag to move. Corner handles scale uniformly (`Shift`: free, `Alt`: about the centre), edge handles stretch one direction, the round handle on top rotates (`Shift`: 15° steps). |
| Pan | `H` | Drag to pan. `Space`+drag or the middle mouse button also pans. |
| Trajectory (freehand) | `T` | Draw in the direction of the flow. The stroke is smoothed and converted to cubic Béziers. Ends snap to targets (see *Magnet*). If an end only fits the other way round, the trajectory is reversed automatically. |
| Bézier pen | `P` | Click for corner points; click-drag to pull out smooth handles. `Enter` or double-click finishes, `Esc` cancels. |
| Sink | `1` | Filled disc. Can be shown as a focus (cw/ccw spiral). |
| Source | `2` | Open disc. Can be shown as a focus. |
| Saddle | `3` | Has 4 terminals: two unstable separatrices (pointing out) and two stable ones (pointing in). Drag the diamond handles to rotate it or to change the angle between its eigendirections. |
| Saddle-node | `4` | Codimension 1, only when *Structurally stable* is off. Transversally stable or unstable. Half-filled disc on the node side, with three separatrix terminals plus the node sector. |
| Periodic orbit | `C` | Click for an ellipse, or draw a closed loop. Stable = solid, unstable = dashed. Semi-stable (codim 1) = solid with a dashed shadow on the repelling side. Arrows show the orientation. |
| Boundary | `B` | A transversal the flow crosses in one direction only, shown by small crossing arrows. A closed loop is created as an *inflow* boundary (trapping region); flip it in the panel. |
| Label | `L` | Click to place text. Inside `$…$`: `_` gives a subscript, `^` a superscript, and `\alpha`, `\Gamma`, … give Greek letters. TikZ receives the text verbatim. |

## Editing

* **Nodes:** a single selected curve shows its Bézier nodes (squares) and handles (dots). Drag them; `Alt`-drag a handle to break smoothness.
* **Adding and removing nodes:** double-click a curve to insert a node, double-click a node to toggle smooth/corner, or select a node and press `Del` to remove it.
* **Detaching ends:** dragging an attached end of a trajectory detaches it; drop it near a target to re-attach. An end attached to a sink, source or node sector has no handle; use *Detach end* in the panel.
* **Joining:** drop the free end of a trajectory on the free start of another (or draw a new one starting or ending there) to merge them into one trajectory.
* **Simplify** (`S` or panel button): smooths the selected curve(s) and refits them with about 30% fewer nodes. Press repeatedly to straighten wobbly hand-drawn curves. The smoothing applied to freehand strokes (low/medium/high) is set in the right panel when nothing is selected.
* **Other keys:**
  * `Del` deletes.
  * `Ctrl+D` duplicates. A duplicate is never refused, so you can drag the copy into place.
  * `R` reverses a trajectory or orbit, or flips a boundary.
  * `Ctrl+A` selects all; arrow keys nudge.
  * `Ctrl+Z` / `Ctrl+Shift+Z` undo and redo.
  * `Esc` goes back to Select.
* **View:** the mouse wheel zooms around the cursor; `+`/`−` zoom and `0` fits the drawing.
* **Panels:** the right panel edits the selected element (type, rotation, stability, arrows, line weight, …). With nothing selected it shows drawing-wide settings: arrow size and freehand smoothing.

## Magnet (snapping)

A trajectory end within about 14 screen pixels of one of these targets snaps to it:

* the centre of a sink or source;
* a saddle or saddle-node terminal (the curve's tangent is aligned with the separatrix);
* a periodic orbit (the trajectory runs into or out of it tangentially, in the orbit's direction, just inside or outside it);
* a boundary;
* the free end or start of another trajectory (join).

The ring is green for a valid target and red for an invalid one; the status bar says why.

## Consistency checks

The *Consistency* panel lists:

* **Problems (red):**
  * curves that cross: trajectories, separatrices, periodic orbits, or a curve with itself;
  * a curve running through an equilibrium;
  * a crossing of a boundary against its direction;
  * a periodic orbit crossing a boundary;
  * a trajectory attached the wrong way, e.g. flowing out of a sink or into an unstable separatrix;
  * two trajectories on one separatrix terminal.
* **Warnings (amber):**
  * the index sum inside a periodic orbit or closed boundary is not 1 (nodes and foci count +1, saddles −1, saddle-nodes 0);
  * non-generic elements present while *Structurally stable* is on.

**Refuse crossings** (off by default): when on, an edit that introduces a new problem is undone, the offending spot flashes red and the status bar gives the reason. When off, problems are only flagged.

**Structurally stable** (on by default):

* **On:** only hyperbolic elements are allowed, and new saddle connections are refused.
* **Off:** saddle-nodes, semi-stable periodic orbits and saddle connections / homoclinic loops are allowed.

## Files and export

* **Saving:** *Save*, *Save as…* and *Open…* use a real file dialog in Chrome and Edge, so you choose the folder and name. After that, *Save* (`Ctrl+S`) keeps writing to the same file. Safari and Firefox do not give web pages this ability, so there files go to the Downloads folder.
* **Autosave:** the drawing is also autosaved in the browser (localStorage) and restored when the page is reopened.
* **Export menu:**
  * SVG
  * PNG (3× resolution)
  * TikZ code (shown in a box to copy, or saved as `.tex`)
  * Print / save as PDF (prints just the drawing)
* **TikZ:** needs `\usepackage{tikz}` and `\usetikzlibrary{arrows.meta,decorations.markings}`. 50 px of drawing = 1 cm; change the `x`/`y` options of the `tikzpicture` to rescale.

## Files in this folder

* `index.html` – the application (HTML + CSS + JavaScript, no dependencies).
* `README.md` – this help.
* `PLAN.md` – specification, design decisions and implementation notes, sufficient to regenerate or extend the tool.
