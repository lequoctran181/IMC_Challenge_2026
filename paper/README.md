# Round 2 article

**Certified Perceptual Mesh Simplification under a 21-Second and 128-KiB Budget**

*Hybrid QEM, cluster-normal memory, renderer-aware topology replay, and fail-closed optimization*

- [`IMC_Challenge_Round2_NEU_AddictedTribes.pdf`](IMC_Challenge_Round2_NEU_AddictedTribes.pdf): final submission-ready article, 28 A4 pages.
- [`IMC_Challenge_Round2_NEU_AddictedTribes.docx`](IMC_Challenge_Round2_NEU_AddictedTribes.docx): editable organizer-template version with native Word Equation objects.
- [`figures/`](figures/): publication-resolution diagrams used in the article.
- [`manuscript.md`](manuscript.md): complete text source.
- [`source/`](source/): deterministic figure and DOCX builders.

The article follows the organizer's template geometry while using the requested five-part scientific structure: Introduction, Related Literature, Methodology, Results and Discussions, and Conclusion. Within those parts it explicitly covers every organizer rubric item—abstract, assumptions and symbols, problem analysis, model building, model solving, model summary, test description, references, and appendices. It also documents controlled hidden-constraint discovery, isolated score decoding, invariant fingerprints, acceptance-frontier reconstruction, and the fail-closed deployment workflow. Its result claims trace to the immutable release record in [`../release/final/`](../release/final/).

All mathematical displays and inline symbols in the released DOCX are editable native Office Math (OMML), with true subscripts and superscripts rather than visible underscore notation. Public sources and the official GitHub artifact are linked directly from the article.

Artifact integrity:

```text
PDF  b473094cf814930f0fef166b58574aad5a299e377fc5b3d4f8e6337d0c766b36
DOCX 0636fd9e81eff3d72eb269e032e3925ba216c6867e71e7bd73833dd69ab65ada
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
