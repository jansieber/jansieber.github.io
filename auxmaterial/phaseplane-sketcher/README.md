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
| Region | `G` | Drag a rectangle. A small bar offers **Copy image** (exactly that rectangle, as a PNG for other apps; curves are clipped at its edge) and **Copy sketch** (the elements lying completely inside, for pasting into this tool). `Esc` or ✕ cancels. |
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

## Copy and paste

There are two kinds of copy: a native sketch (for pasting back into this tool) and an image (for other apps).

* **Copy** (button or `Ctrl/⌘+C`) puts the selected elements on the clipboard as a sketch. `Ctrl/⌘+X` cuts.
* **Paste** (button or `Ctrl/⌘+V`) inserts the sketch.
  * The pasted elements are selected. They are offset by 20 px if they would sit exactly on the originals; if the copied area is off-screen, they go to the centre of the view.
  * Attachments among copied elements are kept; attachments to elements that were not copied are dropped.
  * This works between drawings and between tabs.
  * A copy is also kept in the browser, so Paste still works where the browser refuses to read the clipboard (e.g. some Firefox versions). It also works when the clipboard meanwhile holds something else, such as an image.
  * On iPad, the Paste button may show a small system "Paste" bubble that has to be tapped once.
* **Copy image** (button or `I`) puts a PNG (3×, white background) on the clipboard, ready to paste into note-taking apps, documents or chats.
  * It copies the selection, or the whole drawing if nothing is selected. Use the **Region** tool for an exact rectangle.
  * If the browser refuses image clipboard access, iPad and Android offer the share sheet instead (which also has *Copy*), and desktops save a PNG file.

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
* **iPad / iPhone (any browser):** web pages cannot open a folder dialog there. *Save*, *Save as…* and the exports open a small box instead:
  * **Share… / Save to Files** opens the iOS share sheet; choose *Save to Files* and pick a folder (iCloud Drive, On My iPad, …). *Open…* can read the file back from there.
  * **Download** puts the file into the browser's own downloads area instead.
* **Android (Chrome, Firefox, …):** web pages cannot open a folder dialog there either. *Save*, *Save as…* and the exports open the same small box:
  * **Download** saves into the phone's or tablet's *Downloads* folder; *Open…* opens the system file picker, where Downloads (and Google Drive etc.) can be browsed to load it again.
  * **Share…** sends the file to another app (Google Drive, Files, e-mail, …). It is only shown when the browser supports sharing files.
* **Library…:** keeps named drawings inside the browser (localStorage) with no files involved, and lists them for opening or deleting. Handy on iPad, but tied to that browser on that device.
  * **Export library…** saves all stored drawings as one backup file (`phase-library-DATE.json`).
  * **Import…** adds drawings from such a backup, or from one or more `.phase.json` files. Identical entries are skipped, and name clashes get a suffix such as `name (2)`.
  * Use export and import to move a library between browsers, devices, or copies of the tool loaded from different places.
  * **Which copies share a library:** browser storage belongs to an *origin* (protocol + host + port), not to a folder or file. All tabs showing the tool from the same website share one library; copies on different websites each have their own. Copies opened from disk share one library in Chrome, while Firefox may keep local files separate. Private windows get a temporary library that is deleted when the window closes. Settings (grid, snap, magnet, refuse crossings, smoothing) are shared per origin in the same way.
* **Drag and drop:** a `.phase.json` file can be dropped onto the drawing area to open it.
* **Autosave:** each browser tab autosaves its drawing in its own slot, so two tabs never overwrite each other's restore point.
  * Reloading a tab restores that tab's drawing.
  * A newly opened tab starts from the most recently autosaved drawing of any tab, then keeps its own slot.
  * A duplicated tab is detected and gets its own slot.
  * The 12 most recent slots are kept.
* **Export menu:**
  * **PDF file, cropped to drawing (for LaTeX):** a vector PDF whose page is exactly the drawing plus a small margin. It is made directly, without the print dialog, and can go straight into `\includegraphics`. 50 px of drawing = 1 cm, as in TikZ. Fonts are the standard PDF Times and Symbol, so Greek letters and the label symbols work.
  * **Print / save as PDF – page = drawing…:** the browser's print dialog with the page size set to the drawing. This works in Chrome, Edge and Firefox; other browsers may ignore the page size.
  * **Print / save as PDF – full page…:** the drawing scaled onto a normal paper page.
  * SVG
  * PNG (3× resolution)
  * TikZ code (shown in a box to copy, or saved as `.tex`)
* **TikZ:** needs `\usepackage{tikz}` and `\usetikzlibrary{arrows.meta,decorations.markings}`. 50 px of drawing = 1 cm; change the `x`/`y` options of the `tikzpicture` to rescale.

## Files in this folder

* `index.html` – the application (HTML + CSS + JavaScript, no dependencies).
* `README.md` – this help.
* `PLAN.md` – specification, design decisions and implementation notes, sufficient to regenerate or extend the tool.
