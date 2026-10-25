# Day 23 — 2R Robot Integration Project

Build a modular NumPy 2R robot simulation. Reuse your Day 22 formulas but implement module interfaces.

## Tasks

- TODO 1 in `controllers.py`: four controllers.
- TODO 2 in `controllers.py`: total torque clipping + **boolean** saturation mask.
- TODO 3 in `simulator.py`: semi-implicit Euler integration, logging, and end-effector paths.

The true plant uses `robot.forward_dynamics`; the controller nominal model uses `robot.model_terms`.

## Run

```bash
pip install numpy matplotlib
python test_robot.py
python main.py
```

Tests check model consistency and regression against Day 22's ideal/saturation CTC results. `main.py` plots histories and animates the 2R arm.

## Stretch goals

- Save the animation to GIF or MP4.
- Plot both motor torque commands and actual torques.
- Compare all four controllers on one plot.
- Explore integration step sensitivity.
