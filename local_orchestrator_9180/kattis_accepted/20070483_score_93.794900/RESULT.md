# Kattis submission 20070483

- Judgement: Accepted
- Tests: 7/7
- Score: 93.794900
- Output vertex counts: `[25, 4340, 2848, 3040, 7800, 17100]`
- Authoritative source bytes: 130396
- SHA-256: `7cafab01b3e266989fb485097e485b707ad1397134f86cbe31228265a257a5b0`

This submission lowers the accepted Nefertiti output from 17660 to 17100
vertices while preserving the other five outputs.  The count-only score is
`93.794900599247`, matching the rounded Kattis result.

Local Nefertiti QA before submission:

- 17100 vertices and 34196 faces; validator `VALID`.
- Symmetric Hausdorff ratio: `0.436365967` of the allowed tolerance.
- VPS1024: normal `0.809053016364`, depth `0.997212769120`, combined
  `0.903132892742`.
- The heap pruning cap was reduced to `max(1000000, 6*L+100000)`, allowing
  test 7 to finish while keeping the generated mesh valid.

The 7/7 status and score were verified with `kattis_manager.py list`.  The
archived source was fetched back with `kattis_manager.py source 20070483`.
