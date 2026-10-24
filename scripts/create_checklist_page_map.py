#!/usr/bin/env python3
"""Extract existing RECORD-PE and TARGET working crosswalk items into a page-map CSV.

Standard library only; does not infer Word page numbers or determine compliance.
Usage:
    python scripts/create_checklist_page_map.py RECORD_PE_Checklist_FINAL.docx TARGET_Checklist_FINAL.docx output.csv
"""
import csv
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
ITEM = re.compile(r"^((?:\d+[a-z]?|[A-Za-z]+\d+(?:\.\d+)?))\s*\|\s*(.*?)\s*\|\s*([A-Z_]+)\s*$")


def docx_paragraphs(path: Path):
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
    for p in root.iter(W + "p"):
        txt = "".join(node.text or "" for node in p.iter(W + "t")).strip()
        if txt:
            yield txt


def parse(path: Path, guideline: str):
    rows = []
    for paragraph in docx_paragraphs(path):
        match = ITEM.match(paragraph)
        if match:
            rows.append({
                "guideline": guideline,
                "item": match.group(1),
                "item_topic": match.group(2),
                "working_status": match.group(3),
                "reported_location": "",
                "page_numbers": "",
                "verified_by": "",
            })
    if not rows:
        raise ValueError(f"No item headings found in {path}; confirm source version")
    if len({row["item"] for row in rows}) != len(rows):
        raise ValueError(f"Duplicate checklist item IDs in {path}")
    return rows


def main(args):
    if len(args) != 3:
        raise SystemExit("Expected RECORD_PE.docx TARGET.docx output.csv")
    record, target, dest = map(Path, args)
    if not record.exists() or not target.exists():
        raise SystemExit("Both source Word checklists must exist")
    rows = parse(record, "RECORD-PE/STROBE/RECORD") + parse(target, "TARGET")
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        raise SystemExit(f"Output exists; will not overwrite completed page mapping: {dest}")
    with dest.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"Created {dest} with {len(rows)} source-derived checklist items. Final pagination and verification remain manual.")


if __name__ == "__main__":
    main(sys.argv[1:])
