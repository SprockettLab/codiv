# codiv 0.1.0 — CRAN readiness report

_Generated 2026-09-03. This file lives in `dev/` and is `.Rbuildignore`d, so it
does not ship with the package._

## Current status

`R CMD check --as-cran` on the built tarball (`codiv_0.1.0.tar.gz`):

```
Status: 0 ERRORs | 0 WARNINGs | 3 NOTEs
```

Test suite (run inside the check): **FAIL 0 | WARN 29 | SKIP 0 | PASS 252**.

The three NOTEs:

| # | NOTE | Real? | Action |
|---|------|-------|--------|
| 1 | "New submission" + 3 invalid URLs (HTTP 404) | **Partly real** | "New submission" is expected. The 404 URLs must be fixed — see Blocker 1. |
| 2 | "unable to verify current time" / future file timestamps | No | Local machine had no network time source. Will not appear on CRAN / win-builder. |
| 3 | "'tidy' doesn't look like recent enough HTML Tidy" | No | Local HTML Tidy is old. Environmental; CRAN's is current. |

So the package is **one blocker away** from a clean submission.

---

## Blockers (must fix before submitting)

### 1. DESCRIPTION URLs return 404

`R CMD check` flagged all three URLs:

- `https://github.com/SprockettLab/codiv` — 404 (repo private or not created)
- `https://github.com/SprockettLab/codiv/issues` — 404
- `https://sprockettlab.github.io/codiv/` — 404 ("moved to https://www.sprockettlab.com/codiv/", which also 404s)

CRAN will reject or hand-review a submission with dead URLs in DESCRIPTION.
Pick one:

- **Preferred:** make the GitHub repo `SprockettLab/codiv` public and deploy the
  pkgdown site (the `pkgdown.yaml` workflow is already in `.github/workflows/`)
  *before* submitting, so both URLs resolve.
- Or edit `DESCRIPTION` (`URL:` and `BugReports:`) and `_pkgdown.yml` (`url:`) to
  whatever location is actually live, then rerun `roxygen2::roxygenise()` and
  rebuild. The vignette also links to the issues URL (via `index.md` / articles),
  so update those too.

### 2. Confirm the maintainer email

DESCRIPTION lists `Daniel.Sprockett@wfusm.edu`. CRAN sends the submission
confirmation link there and all future correspondence. Make sure that inbox is
monitored (your account/memory also has `daniel.sprockett@gmail.com`). If you
want the gmail address as maintainer, change it in `DESCRIPTION` `Authors@R`
before submitting — it cannot be changed later without a manual CRAN request.

---

## What was changed to get here

Source / documentation edits (all committed on the working branch as
uncommitted changes right now — review and commit):

1. **`vignettes/codiv.Rmd`** — the tutorial was broken: it called `codiv()` on
   the bundled (unrooted) symbiont tree, but `check_inputs()` now requires a
   rooted symbiont tree (commit `c3a30d2`). Added
   `Symbiont_Tree <- castor::root_at_midpoint(Symbiont_Tree)` with an
   explanation. Also fixed stale `Collapsed_*` column names in the "Interpreting
   Results" section (the `Collapsed_` prefix was dropped in commit `26eed2a`) and
   rewrote the "full vs. collapsed" section to describe current behaviour
   (all methods run on the collapsed subtree; `uncollapsed_hommola = TRUE` is
   opt-in and adds `Uncollapsed_Hommola_*`).
2. **`R/filter_results.R`** — the `@examples` block used
   `Uncollapsed_Hommola_pvalue`, a column that does not exist unless
   `uncollapsed_hommola = TRUE`. The donttest example would have errored on CRAN.
   Switched examples and param docs to `Hommola_pvalue` / `Hommola_r`.
3. **`R/codiv.R`, `R/codiv_class.R`** — fixed the misspelled output column
   `MatchingSplitDistanc` → `MatchingSplitDistance` (in both the column list and
   the assignment). Expanded the `@return` column list to include `Hommola_r`,
   `Topology_*`, and the opt-in `Uncollapsed_Hommola_*`.
4. **`R/codiv_null_scans.R`, `R/host_codiv_summary.R`, `R/molecular_clock.R`** —
   the `statistic`/`pvalue` param docs claimed the `NULL` default "resolves to
   `Uncollapsed_Hommola_r` when present"; the code
   (`.default_hommola_cols()`) always returns `Hommola_r`. Docs corrected.
5. **`R/zzz.R`** — deleted. Its `.onAttach` set an option
   (`codiv.suppress_phylo_warnings`) that nothing reads, and setting options in
   `.onAttach` is discouraged by CRAN policy.
6. **`README.Rmd` / `README.md`** — rewrote for accuracy (mentions all four
   methods, the bundled dataset, correct install command). Regenerated
   `README.md`.
7. **`cran-comments.md`** — rewritten to describe the actual 0.1.0 state and the
   test environments still to be filled in.
8. **`inst/guides/*.md`** — moved to `dev/guides/`. These two beginner guides
   were unreferenced and badly out of date (every column name still had the
   `Collapsed_` prefix). They no longer ship with the package. Either update and
   reinstate them, or leave them as dev notes; the vignette + reference docs
   already cover the same ground.
9. **`vignettes/codiv_bac120_v1.Rmd`** — moved to `dev/`. It was an untracked
   scratch draft with hard-coded local paths (`/Users/danielsprockett/...`) and
   `output: html_document` (no vignette engine), which broke `R CMD build`.
10. **`.Rbuildignore` / `.gitignore`** — added `dev/`, `..Rcheck`, `*.tar.gz`,
    `README.html`, `.Rhistory`, etc.
11. **`_pkgdown.yml`** — added `filter_results` to the reference index (it was
    exported and documented but missing from the site).

## Non-blocking issues / recommendations

- **`..Rcheck/` is committed to git** (`..Rcheck/00check.log`,
  `..Rcheck/R_check_bin/*`). It is `.Rbuildignore`d so it does not affect the
  CRAN tarball, but it should not be in the repo. Run:
  `git rm -r --cached ..Rcheck && git commit -m "Drop committed check artifacts"`.
- **29 incidental test warnings.** All are `codiv()`'s own advisory warnings
  (`permutations < 99`, "not fully bifurcating") surfacing in tests that use tiny
  trees and don't assert on them. Harmless (`R CMD check` reports tests OK).
  Optional: bump `permutations` to 99 or wrap with `suppressWarnings()` in those
  tests for a silent run.
- **`data-raw/make_logo.R`** uses `hexSticker` and `ggplot2`, neither declared;
  fine because `data-raw/` is `.Rbuildignore`d and the script is dev-only.
- **`tests/performance/benchmark.R`** is not a testthat test and not run by
  `test_check()`; it is inside `tests/` but harmless. Consider moving to `dev/`.
- **`inst/extdata/Symbiont_Tree.treefile` is unrooted.** Left as-is (real
  inference output) and the vignette now roots it explicitly, which is a better
  teaching example than shipping a pre-rooted tree. If you would rather ship a
  rooted tree, re-root it and update the vignette back.
- Consider adding a `NEWS.md` (`* Initial CRAN release.`) — CRAN likes it and
  pkgdown renders a changelog from it.
- Consider `_R_CHECK_CRAN_INCOMING_REMOTE_` will re-test URLs on CRAN's side;
  make sure they are live at submission time, not just "soon".

---

## Step-by-step: submitting codiv to CRAN

### 0. Prerequisites (one time)

- A LaTeX toolchain for the PDF manual. This machine now has TinyTeX
  (`~/Library/TinyTeX`); `tinytex::install_tinytex()` if you move machines.
- Pandoc (RStudio bundles it; the standalone check here used
  `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64`).
- Packages: `install.packages(c("devtools", "rhub", "urlchecker", "spelling"))`.

### 1. Fix the blockers above

- Resolve the URLs (Blocker 1). Then:
  ```r
  urlchecker::url_check()      # should report no broken URLs
  ```
- Confirm the maintainer email (Blocker 2).

### 2. Final local checks

```r
# from the package root
devtools::document()                     # regenerate man/ + NAMESPACE
devtools::check(remote = TRUE, manual = TRUE)   # aim for 0/0/0 (1 note: "New submission")
```

Or from the shell (what was used here):

```sh
export RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64"
export PATH="$HOME/Library/TinyTeX/bin/universal-darwin:$RSTUDIO_PANDOC:$PATH"
R CMD build codiv
R CMD check --as-cran codiv_0.1.0.tar.gz
```

Optional but recommended:

```r
devtools::spell_check()
goodpractice::gp()        # style / complexity hints, not required
```

### 3. Check on other platforms

CRAN requires the package to pass on Windows and a second OS. Use:

```r
# Windows (R-devel and R-release), results emailed to the maintainer address
devtools::check_win_devel()
devtools::check_win_release()

# macOS builder
devtools::check_mac_release()

# R-hub v2 (runs in GitHub Actions on your repo; repo must be a git remote)
rhub::rhub_setup()      # one time
rhub::rhub_check(platforms = c("linux", "windows", "macos"))
```

Wait for the win-builder emails (usually < 1 hour). Fix anything that appears
only there (common: non-ASCII, `\donttest` timing, missing Suggests).

### 4. Fill in cran-comments.md

Replace the `<fill in>` lines in `cran-comments.md` with the actual results,
e.g.:

```
## Test environments
- local macOS 26 (aarch64), R 4.5.1
- win-builder devel + release (2026-09-03) — OK
- macOS builder — OK
- R-hub linux/windows/macos (R-devel) — OK

## R CMD check results
0 errors | 0 warnings | 1 note
* New submission.
```

### 5. Submit

```r
devtools::submit_cran()
```

This uploads the tarball to <https://cran.r-project.org/submit.html> and
pre-fills the form from `cran-comments.md`. Alternatively upload
`codiv_0.1.0.tar.gz` manually at that URL.

Then:

1. CRAN emails the maintainer a confirmation link — **click it** (submission is
   not queued until you do).
2. Automated incoming checks run (minutes to a few hours). If they pass, a
   human reviewer may still look, since it is a new submission.
3. Respond to any reviewer email promptly and in-thread. Typical asks for a
   first submission: single-quote package names in DESCRIPTION, expand acronyms,
   add `\value` to any `.Rd` missing it (currently all have it), reduce example
   runtime. If asked for changes, bump to 0.1.1 (or keep 0.1.0 if they say so),
   fix, re-run checks, `submit_cran()` again, reply to the thread.
4. On acceptance the package appears on CRAN within a day and builds start for
   all platforms. Watch the CRAN check results page for a week or two.

### 6. After acceptance

```sh
git tag v0.1.0 && git push --tags
```

- Create a GitHub release for the tag.
- Add the CRAN badge to the README:
  `[![CRAN status](https://www.r-pkg.org/badges/version/codiv)](https://CRAN.R-project.org/package=codiv)`
- Start a `codiv 0.1.1` section in `NEWS.md` for the next cycle.
