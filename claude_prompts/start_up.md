# Getting started

## Goals

This repository will be used to reduce and analyze the first Keck/HIRES spectra
taken for exoplanet detection, using PypeIt.  We will assess how well a modern
reduction of that original data recovers the published radial-velocity result,
and share the reduction recipe and diagnostics with the public.

## Prompts

1. Read this file.  Execute the 1st task under "Claude/CLAUDE.md file"
2. Read this file.  Execute the 1st task under "Claude/Skills"
3. Read this file.  Execute the 1st task under "Claude/Settings"
4. Read this file.  Execute the 1st task under "Basic start up"

## Claude

### CLAUDE.md file

1. Examine all of the current files in the repository.  Then generate a basic
   CLAUDE.md file for this project.  Have it indicate:

    - I will perform git commands.  Read-only git (status, diff, log, show,
      branch) is fine; anything that changes repository state is mine to run.
    - If you do any calculation, generate it as a python script and write it to
      disk so that I can add it to the Repository.
    - If you need to run Python, use the "pypeit14" conda environment.
    - The PypeIt source itself lives at
      `/Users/xavier/Projects/PypeIt/PypeIt` and the development suite at
      `/Users/xavier/Projects/PypeIt/PypeIt-development-suite`; consult them for
      HIRES reduction behavior rather than guessing.
    - Read `claude_prompts/` before acting, and do only the numbered task you
      were pointed at.
    - Log your work in the Logs section of the relevant prompt doc.

### Skills

1. Copy over the skills/ files from the IOPtics repository
   (`/Users/xavier/Oceanography/python/IOPtics/.claude/skills`) into
   `.claude/skills/` here.  That is `critical-partner` and `grill-me`.

### Settings

1. Generate a `.claude/settings.json` for the project that allows you to run
   most bash commands and run Python (use the "pypeit14" conda environment).
   Start from the IOPtics settings.json, but copy the policy, not the
   accumulated path-specific allow entries — those point at IOPtics' own
   scratchpad and test files.  Keep the deny list, in particular
   `git push` / `git commit` / `git reset` / `git rebase`, so the permission
   system enforces the git rule above.  Change the conda environment from
   `ocean14` to `pypeit14`.

## Basic start up

1. Generate the basic files that one needs for a Python GitHub repository, e.g.
   a file for dependencies.  Examine the other Repositories in
   `Projects/PypeIt` (and `Oceanography/python` if useful) to see how I tend to
   organize things.  Note the package directory must be `first_hires_exoplanet`
   or similar — the repo name has hyphens, which are invalid in module names.
   Log your work in the Logs section.

### Report

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs

### 2026-09-19 (Generated CLAUDE.md)

Examined every file in the repository — it is essentially empty: `README.md` (a
single `# first-hires-exoplanet` title line), `LICENSE` (BSD 3-Clause,
"Copyright (c) 2026, PypeIt"), and `claude_prompts/start_up.md`. Single commit
`bb8c4c4 "Initial commit"` touching only LICENSE and README. So the generated
CLAUDE.md is intentionally minimal and will grow as the codebase develops.

Created `CLAUDE.md` at the repo root with five sections: Project (goal text from
`## Goals` above, plus a note on the current empty state), Workflow conventions
(git is the user's — read-only git permitted, state-changing git not;
calculations written to disk as Python scripts; read `claude_prompts/` and do
only the task pointed at; log dated entries under `## Logs`), Environment
(`pypeit14` conda env, invoked as `conda run -n pypeit14 python ...`), and
Related repositories on disk (PypeIt source and development suite).

Modeled the file's shape on `shane-telescope/CLAUDE.md`, the closest PypeIt
sibling.

What I learned:
- The conda environment on this machine is **`pypeit14`**, not `pypeit14b`.
  `shane-telescope/claude_prompts/start_up.md` and its CLAUDE.md both specify
  `pypeit14b`, which does not exist in `conda env list` here — presumably a
  different machine. Used `pypeit14` throughout.
- PypeIt's HIRES support is `pypeit/spectrographs/keck_hires.py`, spectrograph
  name `keck_hires`. Recorded the path in CLAUDE.md.
- `PypeIt-development-suite/RAW_DATA/` has no HIRES directory, so the dev suite
  will not supply example raw HIRES frames — real data will have to come from
  KOA or elsewhere. Flagging for the data-acquisition step later.
- Caution for the remaining prompts: `shane-telescope/CLAUDE.md`'s "## Project"
  paragraph is stale — it still describes IOPtics' IOP algorithms rather than
  the Shane telescope. These docs get cloned between repos, so sibling files are
  a style template, not a content source.

No git commands were run other than read-only `git log`.

### 2026-09-19 (Copied skills from IOPtics)

Copied the skills from the IOPtics repository
(`/Users/xavier/Oceanography/python/IOPtics/.claude/skills`) into this repo at
`.claude/skills/`. Two skills came over, one `SKILL.md` each:

- `critical-partner` — constructive-disagreement mode: surface an untested
  assumption before agreeing, argue the strongest opposing case, retreat only on
  new evidence rather than on pushback, review weaknesses first, name emotional
  investment, say plainly when no flaw can be found rather than inventing one,
  end with a question instead of a summary. No praise or flattery openers.
- `grill-me` — interview the user one question at a time down each branch of the
  design tree until shared understanding is reached; explore the codebase
  instead of asking anything the codebase can answer; offer a recommended answer
  with every question.

Created the `.claude/` tree, which did not previously exist. Verified with
`diff -r` that the copy is byte-identical to the IOPtics source.

What I learned:
- The skills live under `.claude/skills/`, not a top-level `skills/` directory —
  `.claude/skills/` is the standard Claude Code location, so the phrase "the
  skills/ files" in the task refers there.
- Both SKILL.md files are entirely repo-agnostic: neither mentions IOPtics,
  ocean optics, or the `ocean14` environment. So no adaptation was needed for
  this repo, unlike settings.json (next task), which the prompt doc already
  warns carries IOPtics-specific cruft.
- `.claude/settings.json` does not exist yet — prompt #3 will create it. Note
  `.claude/settings.local.json` is machine-local and should not be committed;
  worth checking the `.gitignore` handles it during prompt #4.

No git commands were run.

### 2026-09-19 (Generated .claude/settings.json)

Created `.claude/settings.json`, starting from IOPtics'
(`/Users/xavier/Oceanography/python/IOPtics/.claude/settings.json`) and keeping
the policy rather than the accumulated entries. Validated as JSON: 42 allow,
7 deny, 1 ask; no remaining `ocean14` references.

Kept unchanged: the read-only shell allow-list (`ls`, `cat`, `cd`, `pwd`,
`echo`, `find`, `grep`, `rg`, `head`, `tail`, `wc`, `sort`, `which`, `env`), the
benign file ops (`mkdir`, `touch`, `cp`, `mv`), read-only git (`status`, `diff`,
`log`, `show`, `branch`), the python/pip/pytest/jupyter/conda entries, the full
deny list (`sudo`, `rm -rf /`, `rm -rf ~`, and `git push`/`commit`/`reset`/
`rebase`), and `ask` on `Bash(rm:*)`.

Changed for this repo:
- `conda run -n ocean14:*` → `conda run -n pypeit14:*`.
- Dropped the six trailing path-specific entries: two Crossref `curl` one-offs
  and their two matching `python3 -c` JSON-parsing incantations, two
  hard-coded-interpreter `pytest` invocations (one of which still points at a
  Linux `/home/xavier/miniconda3` path that does not exist on this machine), a
  `python -c "import bing..."` probe, and a `grep` against an IOPtics report
  file. All were session cruft, per the prompt doc.
- Swapped the WebFetch publisher domains. IOPtics allows Optica, Wiley (ASLO),
  and Taylor & Francis — ocean-optics journals, useless here. Substituted the
  astronomy equivalents: ADS, arXiv, IOP Science, plus the Keck Observatory
  Archive (KOA) and the PypeIt docs. Kept `doi.org`. This is a judgment call
  beyond the literal instruction: I read "copy the policy" to mean literature
  checking stays first-class, with the domains translated to this field.
- Added `Bash(pypeit_setup:*)` and `Bash(run_pypeit:*)` — the two core PypeIt
  CLI entry points this project will lean on. Also added `Bash(diff:*)`, which
  I wanted during the skills task.

What I learned:
- IOPtics' allow-list is a good illustration of the accretion problem the
  start-up notes warn about: roughly a sixth of its entries are one-off escaped
  shell incantations that can never match again, and one references a machine
  the user no longer works on. Worth re-pruning this file periodically rather
  than only at creation.
- `Bash(cd:*)` is in the inherited allow-list, but note that `cd` inside a
  compound Bash command can still trigger a prompt; absolute paths are the safer
  habit in this repo.
- Still outstanding for prompt #4: there is no `.gitignore` yet, and
  `.claude/settings.local.json` (machine-local, uncommitted) needs to be listed
  in it. `.claude/settings.json` itself *should* be committed — it is the
  committed policy.

No git commands were run.
