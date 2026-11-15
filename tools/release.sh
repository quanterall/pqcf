#!/bin/sh
# The files of a GitHub release for a tag: the source archive, its SHA-256 and the release notes (the tag's section
# of CHANGELOG.md).
#   sh tools/release.sh v1.4.0                  writes dist/v1.4.0/
#   gh release create v1.4.0 dist/v1.4.0/pqcf-v1.4.0.tar.gz dist/v1.4.0/SHA256SUMS --notes-file dist/v1.4.0/NOTES.md
# The archive is git archive of the tag, so anyone can rebuild it and compare the checksum (same git and gzip):
#   git archive --format=tar --prefix=pqcf-v1.4.0/ v1.4.0 | gzip -n | sha256sum
set -eu
tag=${1:?usage: sh tools/release.sh <tag>}
git rev-parse -q --verify "refs/tags/$tag" >/dev/null || { echo "no tag $tag" >&2; exit 1; }
version=${tag#v}
out=dist/$tag
mkdir -p "$out"
git archive --format=tar --prefix="pqcf-$tag/" "$tag" | gzip -n > "$out/pqcf-$tag.tar.gz"
(cd "$out" && sha256sum "pqcf-$tag.tar.gz" > SHA256SUMS)
# the section "## <version> (<date>)" of the tag's own CHANGELOG.md (the current one if the tag has none), up to the
# next "## "
{ git show "$tag:CHANGELOG.md" 2>/dev/null || cat CHANGELOG.md; } | awk -v v="$version" '
  index($0, "## " v " ") == 1 { on = 1; next }
  on && /^## / { exit }
  on { print }' > "$out/notes.tmp"
[ -s "$out/notes.tmp" ] || { echo "no section for $version in CHANGELOG.md" >&2; exit 1; }
{ cat "$out/notes.tmp"
  printf '\n## Checksum\n\n```\n%s\n```\n\n' "$(cat "$out/SHA256SUMS")"
  printf 'The archive is `git archive` of the tag: `git archive --format=tar --prefix=pqcf-%s/ %s | gzip -n | sha256sum`\n' "$tag" "$tag"
  printf 'gives the same checksum (with the same versions of git and gzip).\n'; } > "$out/NOTES.md"
rm "$out/notes.tmp"
echo "written $out: $(ls "$out" | tr '\n' ' ')"
