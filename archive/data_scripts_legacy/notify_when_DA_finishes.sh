#!/usr/bin/env bash
# Wait for the MSstatsTMT DA R script to finish, then show a macOS notification.
# Run from repo root: ./data/scripts/notify_when_DA_finishes.sh

echo "Watching for R process (DA_subtype_MSstatsTMT_PDC000120.R)..."
while pgrep -f "DA_subtype_MSstatsTMT_PDC000120.R" > /dev/null 2>&1; do
  sleep 15
done
echo "DA finished at $(date)."
osascript -e 'display notification "MSstatsTMT subtype DA finished." with title "Proteomics"'
exit 0
