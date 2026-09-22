# first-hires-exoplanet

Reducing and analyzing the first Keck/HIRES spectra taken for exoplanet
detection, using [PypeIt](https://github.com/pypeit/PypeIt).

The goal is to assess how well a modern reduction of that original data recovers
the published radial-velocity result, and to share the reduction recipe and
diagnostics with the public.

## Installation

```
conda activate pypeit14
pip install -r requirements.txt
pip install -e .
```

## Layout

- `first_hires_exoplanet/` — the Python package (`tests/` beneath it)
- `claude_prompts/` — task prompts and a dated log of the work
- `CLAUDE.md` — conventions for working in this repository

## Authors

- J. Xavier Prochaska
- Claude
