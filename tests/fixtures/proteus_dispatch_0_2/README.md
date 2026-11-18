Specs written against ymf 0.2 / main before 2026-11-17, kept as they were, so that every ymf release still loads them
(tests/test_stabilization_compat.py). Copied byte for byte from proteus `cekees/main` b9d98f5d, `test/dispatch/`:

| file | git blob |
|---|---|
| `two_phase_level_set.ymf` | 59cb5110b87d3557718f98014f2cf6f44b54e4eb |
| `dambreak.ymf` | 88aa84631839ec89bbdbda6bd190240da22bec80 |

Do not edit them or bring them up to date: their old spellings (`lag`, `sc_beta`, `sc_uref`) and RDLS's `stabilization`
block, which states `shock_capturing` and no `method`, are what they test.
