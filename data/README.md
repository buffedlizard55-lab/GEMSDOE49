# `data/` — the official competition rasters (not committed)

| filename here | official artifact | bytes | sha256 |
| --- | --- | --- | --- |
| `training_features.tif` | `gems-geodawn-numerical-features.tif` | 418,912,844 | `4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5` |
| `existing_faults.tif` | `existing_faults.tif` (the catalogue / labels) | 425,830 | `7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093` |
| `example_submission.tif` | `example_submission.tif` (a.k.a. `sample_submission.tif`) | 1,599,597 | `2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc` |

They are excluded by `.gitignore` because the largest is 419 MB. Fetch and verify them with:

```bash
bash scripts/download_competition_data.sh          # official mirrors, sha256-pinned, fail-closed
python3 scripts/prepare_data.py                    # audits grid, CRS, dtype, footprint, band tags
```

`prepare_data.py` exits non-zero if any pin in `src/gems49/spec.py` does not match, so a green run is
the evidence that the data in this directory is the data these results were computed from.

`example_submission.tif` is **not** all zeros: it carries 1.0 at exactly the 60,988 catalogue pixels.
See `docs/research/knowledge-base.md` §8.

---

## Provenance and the portal caveat (kept from the earlier merged session)


Large raster files are ignored under `.cache/gems_data/`; none are checked in.

```bash
bash scripts/download_competition_data.sh
```

The official [DrivenData data tab](https://www.drivendata.org/competitions/306/competition-doe-gems/data/) requires a registered login and redirects unauthenticated access to the login page. This sandbox cannot use that portal account. For this research run, the script restores the same named inputs from public, owner-hosted GitHub mirrors at pinned commits and verifies exact byte counts and SHA-256 digests before accepting them. This verifies mirror consistency, **not** organizer authentication. The provenance caveat is repeated in `manifest.json` and every experiment receipt.

The official competition overview and problem page remain authoritative for the challenge feature/label descriptions and output format. The local `labels.tif` is the existing-fault catalogue proxy only; it is not the hidden expert-labelled test set. The `sample_submission.tif` is used strictly as a grid/footprint template and is never treated as model output.

The data-fetch script requires the GitHub CLI (`gh`) to be installed and authenticated in this environment. It does not ask for or store credentials.
