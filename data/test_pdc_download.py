"""
Test PDC GraphQL API: list PSM files for a study and download the first file.

Requires: pip install requests  (or use a venv: python3 -m venv .venv && .venv/bin/pip install requests)

Note: As of 2025, pdc.cancer.gov GraphQL may use "uiFile" instead of "files"
and different field names. If you get 400 Bad Request, the API schema may have
changed; check https://pdc.cancer.gov/ or use the manifest-based downloader
(pdc_manifest_downloader.py) with a CSV manifest from the PDC portal.
"""
import os
import sys

try:
    import requests
except ImportError:
    print("This script requires the 'requests' library. Install with:")
    print("  pip install requests")
    print("  or: python3 -m venv .venv && .venv/bin/pip install requests && .venv/bin/python test_pdc_download.py")
    sys.exit(1)

# Change this to the study you want
STUDY_ID = "PDC000402"

GRAPHQL_URL = "https://pdc.cancer.gov/graphql"

query = """
{
  files(filter: {study_id: "%s", data_category: "PSM"}) {
    file_id
    file_name
    file_size
    download_link
  }
}
""" % STUDY_ID

print("Querying PDC API...")

response = requests.post(
    GRAPHQL_URL,
    json={"query": query},
    timeout=30,
)

data = response.json()

if response.status_code != 200:
    print("HTTP error:", response.status_code, response.text[:500])
    sys.exit(1)

if "errors" in data:
    print("GraphQL errors:", data["errors"])
    sys.exit(1)

files = data.get("data", {}).get("files")

if files is None:
    print("Unexpected response (no data.files):", list(data.get("data", {}).keys()))
    sys.exit(1)

print(f"\nFound {len(files)} PSM files\n")

if len(files) == 0:
    print("No files returned — check study ID")
    sys.exit(1)

# Show first few files
for f in files[:5]:
    print(f["file_name"], f["file_size"])

# Download first file as test
test_file = files[0]

url = test_file["download_link"]
filename = test_file["file_name"]

print("\nDownloading test file:", filename)

r = requests.get(url, stream=True, timeout=60)
r.raise_for_status()

with open(filename, "wb") as f:
    for chunk in r.iter_content(chunk_size=8192):
        f.write(chunk)

print("Download finished.")

# Verify file
print("\nChecking downloaded file...")

size = os.path.getsize(filename)
print("File size:", size, "bytes")

print("\nFirst 5 lines of file:\n")

with open(filename) as f:
    for i in range(5):
        print(f.readline().strip())
