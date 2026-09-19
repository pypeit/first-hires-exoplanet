# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Project

This repository is used to reduce and analyze the first Keck/HIRES spectra taken
for exoplanet detection, using PypeIt. We assess how well a modern reduction of
that original data recovers the published radial-velocity result, and share the
reduction recipe and diagnostics with the public.

The repository is new: at present it holds only `README.md`, `LICENSE` (BSD
3-Clause) and `claude_prompts/`. Structure will be added by the numbered tasks
in `claude_prompts/start_up.md`.

## Workflow conventions

- **Git is handled by the user.** I (the user) will perform all git commands —
  staging, committing, branching, pushing, tagging, etc. Do not run `git add`,
  `git commit`, `git push`, `git reset`, `git rebase` or any other
  state-changing git command unless I explicitly ask. Read-only git inspection
  (`git status`, `git diff`, `git log`, `git show`, `git branch`) is fine.
- **Calculations become scripts on disk.** If you do any calculation, generate
  it as a Python script and write it to disk so that I can add it to the
  repository. Do not leave results as ephemeral inline snippets — they should be
  committable and re-runnable.
- **`claude_prompts/` is the source of instruction.** Read the relevant prompt
  doc before acting, and do only the numbered task you were pointed at — not the
  whole file.
- **Log your work.** After finishing a task, append a dated entry under the
  `## Logs` section of the prompt doc you were working from, using the format
  given in that doc's `## Logging` section: `### YYYY-MM-DD (short summary)`
  followed by what was done *and what you learned about the repository*.

## Environment

- Python code runs in the `pypeit14` conda environment. Never the system Python.
  Invoke it as `conda run -n pypeit14 python ...`.

## Related repositories on disk

Consult these for HIRES reduction behavior rather than guessing at it:

- `/Users/xavier/Projects/PypeIt/PypeIt` — the PypeIt source. The HIRES
  spectrograph class is `pypeit/spectrographs/keck_hires.py` (PypeIt name
  `keck_hires`).
- `/Users/xavier/Projects/PypeIt/PypeIt-development-suite` — the development
  suite, for worked reduction examples and test data conventions.
