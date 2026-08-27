# Physicians' Online Consultation Engagement — De-identified Data & Analysis Code

This repository accompanies the manuscript:

> Wang Y. *Attitude as a modifiable target in the knowledge–practice pathway of physicians' online consultation engagement.* Submitted to *PLOS Digital Health*.

It provides the **de-identified underlying data** and the **analysis code** required to reproduce the main quantitative results reported in the paper, in line with PLOS's data-availability policy.

---

## Repository contents

| Path | Description |
|---|---|
| `data/deidentified_dataset.csv` | De-identified analytic dataset (n = 68 physicians): demographics + Knowledge/Attitude/Practice (KAP) dimension total scores + participation status. |
| `codebook.csv` | Variable-level data dictionary (English + Chinese labels, types, value ranges). |
| `analysis/mediation_analysis.py` | Hayes PROCESS Model 4 mediation (K→A→P), standardized, Bootstrap 5000, seed = 42. |
| `analysis/power_analysis.py` | Post-hoc Monte Carlo power analysis (reproduces the manuscript's reported power). |
| `analysis/generate_figures.py` | Regenerates Fig 1–3 (PNG/SVG/TIFF 300 dpi) from `figures_source/`. |
| `analysis/requirements.txt` | Python dependencies. |
| `figures_source/*.json` | Source data for the three figures. |
| `figures/` | Rendered figure outputs (PNG/SVG/TIFF). |

---

## Ethics & de-identification

- All records were desensitized before analysis. The dataset contains only a pseudonymized ID (`DOC_001`–`DOC_068`), age, sex, professional title, online-consultation participation status, and three KAP dimension total scores.
- **No** real names, contact information, department names, or free-text responses are included.
- The study was approved by the Ethics Committee of the First Affiliated Hospital of Guangxi Medical University (Approval No. 2025-E0657) and conducted in accordance with the Declaration of Helsinki; written informed consent was obtained from all participants.

---

## How to reproduce

```bash
cd analysis
pip install -r requirements.txt
python mediation_analysis.py     # prints indirect effects + 95% CIs
python power_analysis.py         # prints post-hoc power
python generate_figures.py       # writes ../figures/Fig1-3.{png,svg,tif}
```

### Reproduced results (should match the manuscript)

**Mediation (indirect effect a×b, 95% CI):**
| Group | Indirect effect | 95% CI | Significant |
|---|---|---|---|
| Full sample (n=68) | 0.179 | [−0.029, 0.383] | No |
| Participants (n=39) | 0.160 | [0.048, 0.294] | Yes |
| Non-participants (n=29) | 0.069 | [−0.094, 0.262] | No |

**Post-hoc power (Monte Carlo):** Full sample 52.5% · Participants 48.5% · Non-participants 12.0%.

> Methodological note: path coefficients use intercept-including OLS; K/A/P are standardized on the full sample before subgroup slicing. The bootstrap seed is fixed (42) so results are exactly reproducible. The subgroup non-significant mediation should be read as inconclusive (low power), not as evidence of an absent pathway.

---

## License

Data and code are released under **CC-BY 4.0**. You are free to share and adapt with attribution.

## Citation

Wang Y. Physicians' online consultation engagement — de-identified KAP data and analysis code. GitHub repository. [URL to be inserted after repository creation]
