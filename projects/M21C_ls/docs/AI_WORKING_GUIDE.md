# AI Working Guide

Work from the local repository:

`/Users/amfox/Desktop/geosldas-analysis`

Branch:

`paper-august-rebuild`

The manuscript project is under:

`projects/M21C_ls/docs/`

## Authority order

Use the following files as the authoritative working context, in this order:

1. **Current accepted science**

   `projects/M21C_ls/docs/notes/current_scientific_state.md`

2. **Manuscript structure and editorial strategy**

   `projects/M21C_ls/docs/notes/manuscript_rebuild_plan.md`

3. **Coauthor comments and their current disposition**

   `projects/M21C_ls/docs/notes/coauthor_feedback_ledger.md`

4. **Current manuscript being rebuilt**

   `projects/M21C_ls/docs/manuscript/manuscript.md`

5. **Current figure set**

   `projects/M21C_ls/docs/paper_figures/`

6. **Current analysis reports**

   `projects/M21C_ls/docs/`

   In particular, use the relevant snow, trend, breakpoint, and figure-generation reports there for numerical details and provenance.

## Historical source material

`projects/M21C_ls/docs/source_material/draft_052426.pdf`

This May 24, 2026 PDF is available locally but is **not authoritative**. Use it only for:

- older prose worth salvaging
- references
- established methods wording
- dataset descriptions
- historical context

If the May draft conflicts with:

- `current_scientific_state.md`
- current analysis reports
- the current figure set
- `manuscript_rebuild_plan.md`
- later documented decisions

the newer material wins.

## Manuscript rules

- Do not rebuild the paper by simply editing the May draft.
- Treat this as a destructive rebuild around the current science.
- Use P1–P9 from Fig. 1 as the authoritative observing-system periods.
- Treat April 2015 / P6 / SMAP as the dominant observing-system transition.
- Treat the common-window snow-DA water budget as the main hydrologic process result.
- Main figure numbering is currently fixed through Fig. 17.
- Use Fox et al. (2025) and Fox et al. (2026) heavily in Methods to avoid repeating established GEOS LDAS soil-moisture DA details.
- Do not revive analyses marked superseded or unsupported in `current_scientific_state.md`.
- If something remains a factual check in the feedback ledger, do not guess it. Flag it as `[CHECK]`.
- Keep the manuscript concise and restrained; do not expand the old draft just because more analyses now exist.

## Before drafting or editing

Before drafting or editing any manuscript section:

1. Read `current_scientific_state.md`.
2. Read `manuscript_rebuild_plan.md`.
3. Check `coauthor_feedback_ledger.md` for relevant unresolved items.
4. Inspect the current manuscript section.
5. Inspect the relevant analysis reports and figures.
6. Only then draft.

When uncertain, prefer the newest accepted GitHub analysis product over older manuscript text.
