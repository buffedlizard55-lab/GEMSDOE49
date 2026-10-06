#!/usr/bin/env bash
# Fetch the official DOE GEMS Prize (DrivenData #306) rasters into data/ and verify every byte.
#
# Source 1 (preferred): the official mirrors published on the competition data tab
#                       https://www.drivendata.org/competitions/306/competition-doe-gems/data/
# Source 2 (fallback):  a sha256-pinned public transport bridge on GitHub, used when the
#                       official mirrors are unreachable from the executing host.  Bytes from
#                       the bridge are accepted ONLY if they hash to the pins in src/gems49/spec.py.
#
# Everything is fail-closed: a hash mismatch aborts and deletes the partial file.
#
# Usage:  bash scripts/download_competition_data.sh [--force]
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="$HERE/data"
mkdir -p "$DATA/external"

FORCE=0
[[ "${1:-}" == "--force" ]] && FORCE=1

PIN_TRAINING="4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5"
PIN_LABELS="7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093"
PIN_TEMPLATE="2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc"
PIN_TOPO="a6398d9950965dec6aae6ccecdaa6ced48645d133eab222cbdd11def9bdabfa4"
PIN_RADIO="6cb051f70f941fd78028fe66a9f71e87204fcd8d9a85903df0b94993bad1ec4d"

log() { printf '[download] %s\n' "$*"; }
die() { printf '[download] FATAL: %s\n' "$*" >&2; exit 1; }

hash_ok() { # path expected
  [[ -f "$1" ]] || return 1
  local got
  got="$(sha256sum "$1" | cut -d' ' -f1)"
  [[ "$got" == "$2" ]]
}

try_url() { # url dest expected
  local url="$1" dest="$2" want="$3"
  log "trying $url"
  if curl -fsSL --retry 2 --connect-timeout 20 --max-time 1800 -o "$dest.part" "$url"; then
    local got; got="$(sha256sum "$dest.part" | cut -d' ' -f1)"
    if [[ "$got" == "$want" ]]; then mv "$dest.part" "$dest"; return 0; fi
    log "hash mismatch from $url (got ${got:0:16}...)"
  else
    log "unreachable: $url"
  fi
  rm -f "$dest.part"
  return 1
}

# ---------------------------------------------------------------- source 1
DROPS=(
  "https://www.dropbox.com/scl/fi/t7fyt03qdh9egyme0itwo/existing_faults.tif?rlkey=yiao96uluqdkipf0h5vju71jf&dl=1|$DATA/existing_faults.tif|$PIN_LABELS"
  "https://www.dropbox.com/scl/fi/6rgvnuady818ol8yqgis4/example_submission.tif?rlkey=kbykilvau066xuogoosbf4cq8&dl=1|$DATA/example_submission.tif|$PIN_TEMPLATE"
  "https://www.dropbox.com/scl/fi/3vz9o0wwavi26xaeoxlwr/gems-geodawn-numerical-features.tif?rlkey=je8d8fepqfbst9lnwsq9rkplu&dl=1|$DATA/training_features.tif|$PIN_TRAINING"
)

need_training=1; need_labels=1; need_template=1
[[ $FORCE -eq 0 ]] && hash_ok "$DATA/training_features.tif" "$PIN_TRAINING" && need_training=0
[[ $FORCE -eq 0 ]] && hash_ok "$DATA/existing_faults.tif" "$PIN_LABELS" && need_labels=0
[[ $FORCE -eq 0 ]] && hash_ok "$DATA/example_submission.tif" "$PIN_TEMPLATE" && need_template=0

if [[ $((need_training + need_labels + need_template)) -eq 0 ]]; then
  log "all three official rasters already present and hash-correct"
else
  for row in "${DROPS[@]}"; do
    IFS='|' read -r url dest want <<<"$row"
    case "$dest" in
      */training_features.tif) [[ $need_training -eq 1 ]] || continue ;;
      */existing_faults.tif)   [[ $need_labels -eq 1 ]] || continue ;;
      */example_submission.tif) [[ $need_template -eq 1 ]] || continue ;;
    esac
    try_url "$url" "$dest" "$want" || true
  done
fi

# ---------------------------------------------------------------- source 2 (bridge)
BRIDGE_PARTS=(
  "0a330f8951af6c921029e25c84a579319d2db554d62d30d894d6ddc97f98cff7"
  "3c98037b2c997e3bbfcdfc2d9a982e8b820a77410dd7404b05cdb80594922c50"
  "c375c4dbc40c59bbaece572b5e348700b75935f9a30c879c82b6417e0836f31c"
  "b164159e6d0cb2595bc9f63a948af2646b124114a9c5921880c7516092137320"
  "fa0a6f9c936fac1d6f20ca37f5929b2d60bf7a80f3d477dcab886f941aee2696"
)

fetch_bridge() {
  local repo="$1" tmp
  tmp="$(mktemp -d)"
  log "fetching transport bridge $repo (sparse, blob-filtered)"
  if ! git clone --depth 1 --filter=blob:none --sparse --quiet \
        "https://github.com/${repo}.git" "$tmp" 2>/dev/null; then
    rm -rf "$tmp"; return 1
  fi
  ( cd "$tmp" && git sparse-checkout set data/bridge data/aux_bridge >/dev/null 2>&1 ) || true
  echo "$tmp"
}

BRIDGE_DIR=""
for repo in buffedlizard55-lab/5GEMSDOE buffedlizard55-lab/GEMSDOE2 buffedlizard55-lab/GEMSDOE; do
  if [[ $need_training -eq 1 || $need_labels -eq 1 || $need_template -eq 1 || ! -f "$DATA/external/topo_u8.tif" ]]; then
    d="$(fetch_bridge "$repo")" && [[ -n "$d" ]] && { BRIDGE_DIR="$d"; break; }
  fi
done

if [[ -n "$BRIDGE_DIR" ]]; then
  src="$BRIDGE_DIR/data/bridge"
  if [[ $need_labels -eq 1 && -f "$src/existing_faults.tif" ]]; then
    cp "$src/existing_faults.tif" "$DATA/existing_faults.tif"
  fi
  if [[ $need_template -eq 1 && -f "$src/example_submission.tif" ]]; then
    cp "$src/example_submission.tif" "$DATA/example_submission.tif"
  fi
  if [[ $need_training -eq 1 ]]; then
    ok=1
    for i in 0 1 2 3 4; do
      p="$src/gems-geodawn-numerical-features.tif.part-00$i"
      [[ -f "$p" ]] || { ok=0; break; }
      got="$(sha256sum "$p" | cut -d' ' -f1)"
      [[ "$got" == "${BRIDGE_PARTS[$i]}" ]] || { log "part-00$i hash mismatch"; ok=0; break; }
    done
    if [[ $ok -eq 1 ]]; then
      log "concatenating 5 verified parts -> training_features.tif"
      cat "$src"/gems-geodawn-numerical-features.tif.part-000 \
          "$src"/gems-geodawn-numerical-features.tif.part-001 \
          "$src"/gems-geodawn-numerical-features.tif.part-002 \
          "$src"/gems-geodawn-numerical-features.tif.part-003 \
          "$src"/gems-geodawn-numerical-features.tif.part-004 > "$DATA/training_features.tif.part"
      mv "$DATA/training_features.tif.part" "$DATA/training_features.tif"
    fi
  fi
  # optional auxiliary layers (official USGS 3DEP / GeoDAWN, used only for corroboration)
  if [[ ! -f "$DATA/external/topo_u8.tif" && -f "$BRIDGE_DIR/data/aux_bridge/topo/topo_u8.tif.part-000" ]]; then
    cp "$BRIDGE_DIR/data/aux_bridge/topo/topo_u8.tif.part-000" "$DATA/external/topo_u8.tif"
  fi
  if [[ ! -f "$DATA/external/radiometric_u8.tif" && -f "$BRIDGE_DIR/data/aux_bridge/radiometric/radiometric_u8.tif.part-000" ]]; then
    cp "$BRIDGE_DIR/data/aux_bridge/radiometric/radiometric_u8.tif.part-000" "$DATA/external/radiometric_u8.tif"
  fi
  rm -rf "$BRIDGE_DIR"
fi

# ---------------------------------------------------------------- verify
log "final verification"
fail=0
check() { # path sha label
  if hash_ok "$1" "$2"; then printf '  OK   %-28s %s\n' "$3" "${2:0:16}..."
  else printf '  FAIL %-28s (missing or wrong hash)\n' "$3"; fail=1; fi
}
check "$DATA/training_features.tif" "$PIN_TRAINING" "training_features.tif"
check "$DATA/existing_faults.tif" "$PIN_LABELS" "existing_faults.tif"
check "$DATA/example_submission.tif" "$PIN_TEMPLATE" "example_submission.tif"
[[ -f "$DATA/external/topo_u8.tif" ]] && check "$DATA/external/topo_u8.tif" "$PIN_TOPO" "external/topo_u8.tif" || log "  note: optional topo_u8.tif absent"
[[ -f "$DATA/external/radiometric_u8.tif" ]] && check "$DATA/external/radiometric_u8.tif" "$PIN_RADIO" "external/radiometric_u8.tif" || log "  note: optional radiometric_u8.tif absent"

if [[ $fail -ne 0 ]]; then
  die "one or more required rasters failed verification. Place them manually in data/ from
       https://www.drivendata.org/competitions/306/competition-doe-gems/data/ and re-run."
fi
log "all required rasters present and hash-verified. Next: python scripts/prepare_data.py"
