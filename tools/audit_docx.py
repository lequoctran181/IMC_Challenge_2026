#!/usr/bin/env python3
"""Fail-closed structural audit for the released Word article.

This deliberately uses only the Python standard library so GitHub Actions can
verify the committed DOCX without rebuilding it or installing authoring tools.
"""

from __future__ import annotations

import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
ARTICLE = ROOT / "paper" / "IMC_Challenge_Round2_NEU_AddictedTribes.docx"
NS = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    "m": "http://schemas.openxmlformats.org/officeDocument/2006/math",
    "wp": "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing",
    "cp": "http://schemas.openxmlformats.org/package/2006/metadata/core-properties",
    "dc": "http://purl.org/dc/elements/1.1/",
}


def main() -> int:
    errors: list[str] = []
    if not ARTICLE.is_file():
        print(f"missing article: {ARTICLE}", file=sys.stderr)
        return 1
    try:
        with zipfile.ZipFile(ARTICLE) as archive:
            document = ET.fromstring(archive.read("word/document.xml"))
            core = ET.fromstring(archive.read("docProps/core.xml"))
            media = [name for name in archive.namelist() if name.startswith("word/media/")]
    except (zipfile.BadZipFile, KeyError, ET.ParseError) as error:
        print(f"invalid DOCX package: {error}", file=sys.stderr)
        return 1

    headings = []
    for paragraph in document.findall(".//w:p", NS):
        style = paragraph.find("./w:pPr/w:pStyle", NS)
        if style is not None:
            value = style.get(f"{{{NS['w']}}}val", "")
            # The organizer template uses numeric style IDs 1/2/3 for the
            # built-in Heading 1/2/3 names; Word documents created from other
            # templates commonly use Heading1/Heading2/Heading3.
            if value in {"1", "2", "3", "Heading1", "Heading2", "Heading3"}:
                headings.append(value)
    levels = ({"1", "Heading1"}, {"2", "Heading2"}, {"3", "Heading3"})
    if len(headings) < 50 or not all(any(value in headings for value in level) for level in levels):
        errors.append(f"native heading hierarchy is incomplete: {len(headings)} heading paragraphs")

    drawings = document.findall(".//wp:docPr", NS)
    described = [node for node in drawings if node.get("descr", "").strip()]
    if len(drawings) != 7 or len(described) != 7 or len(media) < 7:
        errors.append(f"expected 7 described figures; drawings={len(drawings)}, alt={len(described)}, media={len(media)}")

    equations = document.findall(".//m:oMath", NS)
    if len(equations) < 40:
        errors.append(f"expected native Word equations; found {len(equations)} OMML nodes")
    math_text = "".join(node.text or "" for node in document.findall(".//m:t", NS))
    if "_" in math_text:
        errors.append("a literal underscore is visible inside an OMML equation")

    visible_text = "".join(node.text or "" for node in document.findall(".//w:t", NS))
    for forbidden in ("<!--", "second-place finish", "final Round 2 placement"):
        if forbidden.lower() in visible_text.lower():
            errors.append(f"forbidden/stale visible text: {forbidden}")
    for required in ("Quoc Tran Anh Le", "Quang Minh Ha", "93.830074", "1,498,780"):
        if required not in visible_text:
            errors.append(f"missing required article text: {required}")

    title = core.find("dc:title", NS)
    creator = core.find("dc:creator", NS)
    if title is None or not (title.text or "").startswith("Perception-Aware Mesh Simplification"):
        errors.append("unexpected DOCX title metadata")
    if creator is None or "Quoc Tran Anh Le" not in (creator.text or "") or "Quang Minh Ha" not in (creator.text or ""):
        errors.append("DOCX creator metadata does not name both authors")

    if errors:
        print("DOCX audit failed:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        return 1
    print(f"DOCX audit passed: {len(headings)} native headings, {len(drawings)} described figures, {len(equations)} OMML nodes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
