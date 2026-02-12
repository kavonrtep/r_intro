#!/bin/bash
# Update local repo and sync working copy to ~/r_intro
# Usage: run from the cloned repository directory
#   ./update_r_intro.sh          - pull and copy, skip existing student files
#   ./update_r_intro.sh --force  - pull and overwrite everything

set -e

DEST="$HOME/r_intro"
SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"

# Pull latest changes
echo "Pulling latest changes..."
git -C "$SCRIPT_DIR" pull

# Create destination if needed
mkdir -p "$DEST"

if [ "$1" = "--force" ]; then
    echo "Syncing to $DEST (overwriting all files)..."
    rsync -a --delete --exclude='.git' --exclude='work_dir/' "$SCRIPT_DIR/" "$DEST/"
else
    echo "Syncing to $DEST (keeping existing .R and .qmd files)..."
    rsync -a --exclude='.git' --exclude='work_dir/' \
        --ignore-existing --include='*.R' --include='*.qmd' \
        "$SCRIPT_DIR/" "$DEST/"
    # Now sync everything except .R/.qmd (data, notes, etc.) with overwrite
    rsync -a --exclude='.git' --exclude='work_dir/' \
        --exclude='*.R' --exclude='*.qmd' \
        "$SCRIPT_DIR/" "$DEST/"
fi

echo "Done. Working copy is at $DEST"
