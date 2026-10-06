# Data placement and provenance

Large raster files are ignored under `.cache/gems_data/`; none are checked in.

```bash
bash scripts/download_competition_data.sh
```

The official [DrivenData data tab](https://www.drivendata.org/competitions/306/competition-doe-gems/data/) requires a registered login and redirects unauthenticated access to the login page. This sandbox cannot use that portal account. For this research run, the script restores the same named inputs from public, owner-hosted GitHub mirrors at pinned commits and verifies exact byte counts and SHA-256 digests before accepting them. This verifies mirror consistency, **not** organizer authentication. The provenance caveat is repeated in `manifest.json` and every experiment receipt.

The official competition overview and problem page remain authoritative for the challenge feature/label descriptions and output format. The local `labels.tif` is the existing-fault catalogue proxy only; it is not the hidden expert-labelled test set. The `sample_submission.tif` is used strictly as a grid/footprint template and is never treated as model output.

The data-fetch script requires the GitHub CLI (`gh`) to be installed and authenticated in this environment. It does not ask for or store credentials.
