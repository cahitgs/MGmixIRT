#!/usr/bin/env bash
# Pre-push content check: the repository must not contain tool-generated
# attribution strings. Run from the package root; exits non-zero on any hit.
set -u
PATTERN='claude|anthropic|copilot|generated with|co-authored-by'
HITS=$(grep -rniE "$PATTERN" --exclude-dir=.git --exclude="hygiene_check.sh" . || true)
if [ -n "$HITS" ]; then
  echo "Hygiene check FAILED:"
  echo "$HITS"
  exit 1
fi
LOGHITS=$(git log --all --format='%H %an %ae %B' 2>/dev/null | grep -niE "$PATTERN" || true)
if [ -n "$LOGHITS" ]; then
  echo "Hygiene check FAILED in git history:"
  echo "$LOGHITS"
  exit 1
fi
echo "Hygiene check passed."
