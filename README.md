# 6N Cousin & Sexy Two-Centre Mechanism (Part X)

A quantitative two-centre CRT mechanism for the three omega-distortions of Part IX
— with no fitted parameters.

**Background.** Part IX found that the conditional omega-dependence of prime pairs
on the 6N skeleton splits into three geometry-determined shapes: the twin rises,
cousin and sexy-A fall and coincide, sexy-B is non-monotone. It left a
quantitative mechanism as the open problem. This is it.

**The model.** Extend the Part VII closed-form CRT survival from one centre to a
pair straddling two consecutive centres (N, N+1). Write each pair as two members,
each a wing 6(N+off)+s with offset off in {0,1} and wing s in {±1}:

```
  twin   = ((0,-1),(0,+1))    cousin = ((0,+1),(1,-1))
  sexy-A = ((0,-1),(1,-1))    sexy-B = ((0,+1),(1,+1))
```

A member 6(N+off)+s is divisible by q when N+off ≡ -s·6⁻¹ (mod q). The two-centre
CRT factor f_q is, per prime q>3:
- **q | N** (N≡0): deterministic — each member q-safe unless off ≡ -s·6⁻¹ (mod q);
- **q ∤ N**: averaged over admissible nonzero residues r of N.

Pair-rate model = K·∏_q f_q, evaluated per centre by its true small-prime
divisibility, averaged within each left-centre omega stratum, normalised to omega=1.
**No fitted parameters** — only each pair's wing geometry enters.

**Result (S₁₀, 1.5×10⁹ centres), all three shapes from one formula:**

| shape | max err (omega ≤ 6) |
|-------|--------------------:|
| twin (rise) | 1.5% |
| cousin (fall) | 2.3% |
| sexy-A (fall) | 2.0% |
| sexy-B (non-monotone peak) | 1.9% |

**Why the signs.** For q | N, f_q checks whether the partner's wing avoids the
locked residue. The twin's members both sit on N (offset 0) — rarely locked, so
q|N adds a factor ≈1, and more such primes (higher omega) lift the rate. Cousin
and sexy-A have their partner on N+1 (offset 1) — the offset lands in the locked
residue for more primes, so q|N gives factors <1 and the rate falls. Sexy-B
(two right wings) interpolates, giving the non-monotone peak. The sign is set by
whether the partner lives on the conditioned centre N (rise) or the neighbour N+1
(fall).

**The cousin ≡ sexy-A coincidence.** They share the right member 6(N+1)−1; the
model gives identical curves. Confirmed in data through omega=6 (measured ratio
1.00). At omega=7 they diverge (sexy-A residual 14.4%), but that stratum holds
only ~50 pairs of each; the difference is 0.9σ — small-sample noise, not a model
defect.

> **Scope.** Experimental / computational number theory. Main results for
> omega ≤ 6 on S₁₀; omega=7 is small-sample. The q∤N branch uses a uniform-residue
> approximation (source of the few-percent residuals, largest for the non-monotone
> sexy-B). No claim about the infinitude of any prime constellation.

Part I: doi:10.5281/zenodo.20470367 · V: doi:10.5281/zenodo.20510700 ·
VII: doi:10.5281/zenodo.20518470 · IX: doi:10.5281/zenodo.20520492

---

## Layout

```
.
├── README.md
├── LICENSE                 (MIT)
├── CITATION.cff
├── data/
│   └── cs_mechanism_S10_data.csv   pair, omega, measured, model, err_pct  (S10)
├── code/
│   ├── cs_mechanism.py     two-centre CRT model vs measured, per pair & omega;
│   │                       emits cs_mechanism_S{K}_data.csv  (streaming, low memory)
│   └── make_mech_fig.py    builds the 2-panel figure from ../data
├── figures/                fig_paper10_mechanism.{pdf,png}
└── paper/                  Chen_6N_Paper10.{tex,pdf} + figure
```

## Reproducing

Requirements: Python 3.8+, `numpy`, `matplotlib`.

```bash
pip install numpy matplotlib

# 1. Model vs measured. Default S10 (~18 min). Memory-light streaming.
#    Emits cs_mechanism_S{K}_data.csv.
python code/cs_mechanism.py            # S10
MAXK=9 python code/cs_mechanism.py     # S9 (faster, for validation)

# 2. Figure (reads ../data/cs_mechanism_S10_data.csv).
cd code && python make_mech_fig.py
```

### Conventions (same as Parts I–IX)

- Twin centre N: 6N−1, 6N+1 both prime. omega₍>3₎(N) = #distinct prime factors >3.
- Pairs attributed to the LEFT centre's omega; rates normalised to omega=1.
- dead(q,s,off) = (-s·6⁻¹ − off) mod q : the residue of N making member q-divisible.
- Working primes q ∈ {5,…,47}; tail absorbed in the overall constant (divides out
  under normalisation).
- Engine: complete segmented-sieve factorisation + deterministic interval-sieve
  primality; S₁₀ twin count 23,988,173 matches Part I.

## License

MIT — see `LICENSE`.
