# Round 2 article

**Certified Perceptual Mesh Simplification under a 21-Second and 128-KiB Budget**

*Hybrid QEM, cluster-normal memory, renderer-aware topology replay, and fail-closed optimization*

- [`IMC_Challenge_Round2_NEU_AddictedTribes.pdf`](IMC_Challenge_Round2_NEU_AddictedTribes.pdf): final submission-ready article, 20 A4 pages.
- [`IMC_Challenge_Round2_NEU_AddictedTribes.docx`](IMC_Challenge_Round2_NEU_AddictedTribes.docx): editable organizer-template version.
- [`figures/`](figures/): publication-resolution diagrams used in the article.
- [`manuscript.md`](manuscript.md): complete text source.
- [`source/`](source/): deterministic figure and DOCX builders.

The article follows the organizer's template geometry and evaluates the method through the requested abstract, assumptions and symbols, problem analysis, model building, model solving, model summary, test description, references, and appendices. Its result claims trace to the immutable release record in [`../release/final/`](../release/final/).

Artifact integrity:

```text
PDF  052ba91e04601eab90635b24b6e317b591e51d01d634e074a7875dd2296a5014
DOCX 2e3b5f4c094d5acae21a78ae74dcfbd1acc0ea565f039d12b2761a1f684c3238
```

The original blank organizer template is not redistributed. The DOCX retains its A4 page system, typography hierarchy, and opaque custom properties while replacing all placeholders with the team's work.

## Rebuild from the organizer template

Create an isolated Python environment and install the small authoring toolchain:

```bash
python3 -m pip install -r paper/source/requirements.txt
python3 paper/source/generate_figures.py
python3 paper/source/build_article.py \
  --template "/path/to/Article template for IMC Challenge.docx" \
  --output /tmp/IMC_Challenge_Round2_NEU_AddictedTribes.docx
```

The generated DOCX should then be rendered and visually inspected before release; DOCX layout depends on the office renderer and installed fonts. The committed PDF/DOCX remain the signed publication artifacts.
