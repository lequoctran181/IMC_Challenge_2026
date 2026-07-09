# Exact23 Bumpy Proxy And Generic GN Notes

- New exact23 diagnostics show the hidden `23201/46398` case is regular-valence but not smooth-sphere-like:
  `radcv > .25`, `smooth <= .70`, `coarse > .93`, `sharp > .10`, `d6 > .50`,
  `d5+d6+d7 > .80`, `d4+d8 <= .12`.
- Generated regular-valence bumpy sphere proxies with exact `23201/46398`.
  Best fingerprint matches are `synthetic_exact23_20260709b/blob_lat112_lon209_a80_f11.obj`
  and `a95_f11`: `radcv ~= .25-.27`, `smooth ~= .62`, `coarse ~= .93-.96`,
  `sharp ~= .11-.12`, valence-6 dominated.
- GJ-order versions of the same proxies prove the existing `FL` structured-grid route is not a viable breakthrough:
  forced trial 1/2/3/4 outputs `340/614/904/1406` vertices but only reaches
  `vps_eval512 ~= .61/.65/.67/.70` on bumpy proxies, far below the real visual target.
- W2C literal edits (`AI/BB/angle/lim/rm/w`) rerun isolated are output-identical to active base on the proxy portfolio.
- The generic smooth block in `GN()` is not the useful control point for exact23:
  `skip_generic`, `loose`, `no_tight`, and `range23_loose` are output-identical or timing-only on bumpy exact23.
  Broad `pth90` changes Lucy to a deeper output family (`13067/26130`) and is unsafe without a much stronger guard.
- No candidate from these batches is evidence-backed enough to submit. Next exact23 direction should either improve the existing
  `~11k` QEM output quality/vertex count with a truly exact guard, or pivot to another hidden case with more headroom.
