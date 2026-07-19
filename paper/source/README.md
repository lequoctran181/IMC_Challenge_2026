# Article source

This directory contains the authoring pipeline for the Round 2 article.

- `build_article.py` replaces the organizer-template placeholders with `../manuscript.md`, creates native Word tables and hyperlinks, converts LaTeX to editable Office Math (OMML), embeds inline figures, and applies the template-derived page/typography system.
- `generate_figures.py` regenerates the six schematic/chart figures in both 300-DPI PNG and editable SVG/PDF form.
- `generate_qualitative_figure.py` composes a seventh, unretouched normal/depth/error montage from evaluator PPM dumps using one common crop.
- `data/accepted_milestones.csv` is the evidence table behind the score-progression chart.
- `requirements.txt` pins the Python authoring dependencies at compatible major versions.

The organizer's blank DOCX is intentionally absent. Download it from the official Round 2 documentation and supply it with `--template`. Generated documents must be rendered to PDF and visually checked page by page; a successful script exit is not sufficient layout validation.

The released DOCX also preserves template package parts such as its theme, numbering, font table, and opaque organizer custom properties where feasible.
