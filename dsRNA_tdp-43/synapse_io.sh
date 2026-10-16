# ============================================================================
# synapse_io.sh -- thin wrappers over the `synapse` CLI (ships with synapseclient).
# Matches the lab convention in job_scripts/download_synapse.sh.
#
# Source this from the numbered scripts:  source synapse_io.sh
# Every call is a no-op unless USE_SYNAPSE=1 (set in params.sh).
#
# AUTH: the CLI reads SYNAPSE_AUTH_TOKEN or ~/.synapseConfig. Run `synapse config`
# (or export SYNAPSE_AUTH_TOKEN=...) once; no token is stored in the code.
# ============================================================================

syn_enabled() { [ "${USE_SYNAPSE:-0}" = "1" ]; }

syn_require() {
  command -v synapse >/dev/null 2>&1 || {
    echo "ERROR: 'synapse' CLI not found. Install with: pip install synapseclient" >&2
    exit 1; }
}

# syn_get <synID> <dest_dir>
#   Download an entity into dest_dir under its Synapse filename (no path parsing;
#   the caller references dest_dir/<known_filename>). No-op unless USE_SYNAPSE=1.
syn_get() {
  syn_enabled || return 0
  syn_require
  local id="$1" dest="${2:-.}"
  mkdir -p "$dest"
  echo "[synapse] get $id -> $dest/" >&2
  synapse get "$id" --downloadLocation "$dest"
}

# syn_store <file> <parentFolderId> [used] [executed]
#   Upload (or create a new version of) a file in a Synapse folder, with provenance.
#   `used` and `executed` are optional (synID/URL lists and the producing script).
#   No-op unless USE_SYNAPSE=1.
syn_store() {
  syn_enabled || return 0
  syn_require
  local file="$1" parent="$2" used="${3:-}" exec="${4:-}"
  [ -s "$file" ] || { echo "ERROR: syn_store: file not found: $file" >&2; exit 1; }
  local args=(store "$file" --parentId "$parent")
  [ -n "$used" ] && args+=(--used "$used")
  [ -n "$exec" ] && args+=(--executed "$exec")
  echo "[synapse] store $(basename "$file") -> $parent" >&2
  synapse "${args[@]}"
}
