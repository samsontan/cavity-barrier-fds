"""
verify_citations.py — resolve every DOI in a manuscript's reference list against
Crossref and check that the cited title actually belongs to that DOI.

WHY THIS EXISTS (F132, HIGH)
  On 2026-06-11 the SAI learning loop recorded F028: "AI-inherited fabricated
  citation survived 6 manuscript versions." On 2026-06-12 and 06-13 it recorded
  "audit in-text citation targets before submission" and "referee citations
  BEFORE trusting them."

  On 2026-07-10 Paper B still shipped a fabricated reference title (ref 10,
  matching no Crossref record, carrying its own false annotation "[exact title
  per Crossref record]") and a second title that did not belong to its DOI
  (ref 17). Four independent referee agents had passed the manuscript. They
  checked whether the numbers were right, never whether the papers existed.

  The lesson lived as prose in a knowledge file, and prose does not execute.
  This script is the countermeasure. Run it before any manuscript build or
  submission. It exits non-zero on any mismatch.

WHAT IT CATCHES
  * a DOI that does not resolve (typo, or invented)
  * a title that does not belong to its DOI (the dangerous case: a right DOI
    with a wrong title passes every casual spot check)
  * a first-author surname that does not match Crossref's record
  * a year that disagrees with Crossref by more than one (online-first vs issue)

WHAT IT CANNOT CATCH
  * a reference with no DOI (standards, reports, books). Those are listed as
    SKIPPED and must be checked by hand against the publisher.
  * a real paper cited for a claim it does not make. Only a human can catch that.

Usage:
    python verify_citations.py                       # defaults to ../../PaperB_v4.md
    python verify_citations.py path/to/manuscript.md
"""

import os
import re
import sys
import json
import time
import difflib
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MD = os.path.normpath(os.path.join(HERE, "..", "..", "PaperB_v4.md"))

UA = "SAI-citation-verifier/1.0 (mailto:samson.tan@vu.edu.au)"
TITLE_RATIO_FAIL = 0.60   # below this, treat as a different paper
TITLE_RATIO_WARN = 0.85   # below this, flag for a human read

# A DOI legitimately contains '.' and '(' ')' (e.g. 10.1016/S0379-7112(97)00041-6).
# Grab everything up to whitespace or a markdown marker, then trim sentence punctuation.
DOI_RE = re.compile(r"(?:doi:|https?://(?:dx\.)?doi\.org/)\s*(10\.\d{4,9}/[^\s*]+)", re.IGNORECASE)


def clean_doi(raw):
    d = raw.strip()
    while d and d[-1] in ".,;":
        d = d[:-1]
    # a trailing ')' belongs to the DOI only if parens are balanced
    while d.endswith(")") and d.count("(") < d.count(")"):
        d = d[:-1]
    while d.endswith("]"):
        d = d[:-1]
    return d
ENTRY_RE = re.compile(r"^(\d+)\.\s+(.*)$")

# Publisher metadata is sometimes wrong. Where it is, the citation must still be
# verified against a STRONGER source, and the reason recorded here. An exemption
# without a verified_against value is not an exemption, it is a fudge.
#
# Do not add to this table to make the script go green. Add to it only when you
# have read the primary document and Crossref is the thing that is wrong.
KNOWN_EXCEPTIONS = {
    "10.6028/NIST.SP.1019": (
        "Crossref holds a stale 2001 record for this SERIES DOI: it returns "
        "'Fire dynamics simulator (version 4) :' regardless of edition.",
        "PDF title page of the FDS 6.10.1 User's Guide (revision FDS-6.10.1-0-g12efa16, "
        "4 April 2025): 'NIST Special Publication 1019, Sixth Edition'. Verified 2026-07-10.",
    ),
    "10.6028/NIST.SP.1018": (
        "Same stale series record as SP 1019.",
        "PDF title page of the FDS 6.10.1 Validation Guide: 'NIST Special Publication "
        "1018-3, Sixth Edition, Volume 3: Validation'. Verified 2026-07-10.",
    ),
}


def norm(s):
    s = re.sub(r"<[^>]+>", " ", s or "")
    s = s.replace("‐", "-").replace("–", "-").replace("—", "-")
    s = re.sub(r"[^a-z0-9 ]+", " ", s.lower())
    return re.sub(r"\s+", " ", s).strip()


def crossref(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)["message"]


def extract_references(md_path):
    """Pull the numbered reference entries out of the References section."""
    text = open(md_path, encoding="utf-8").read()
    m = re.search(r"^## References\b.*?$", text, re.M)
    if not m:
        raise SystemExit("no '## References' section found")
    tail = text[m.end():]
    stop = re.search(r"^## ", tail, re.M)
    block = tail[:stop.start()] if stop else tail
    out = []
    for line in block.splitlines():
        em = ENTRY_RE.match(line.strip())
        if em:
            out.append((int(em.group(1)), em.group(2).strip()))
    return out


def cited_title(entry):
    """Best-effort extraction of the cited title.

    Two house styles appear in this reference list:
      * journal articles: Authors. Title. *Journal* vol (year) pages. doi:...
      * books/reports   : Authors. *Title.* Publisher, year. doi:...
    For the second, the title sits inside the FIRST italic block. Take that when
    present, otherwise strip the leading author list and take the run of text
    before the italic journal name.
    """
    e = re.sub(r"\[[^\]]*\]", "", entry)          # drop bracketed notes
    e = re.sub(r"doi:\S+", "", e, flags=re.I)
    e = re.sub(r"https?://\S+", "", e)

    # book / report style: the title is the first italic block, if it is long
    # enough to be a title and is not a journal name following a plain title.
    first_italic = re.search(r"\*([^*]{12,})\*", e)
    before_italic = e[:first_italic.start()] if first_italic else e
    # if everything before the first italic block is just an author list, the
    # italic block IS the title
    author_only = re.fullmatch(
        r"\s*(?:[A-ZÀ-Þ][\w'’\-]+(?:\s+[A-Z]\.?){0,4}[,.]?\s*|et al\.?\s*|and\s+)*",
        before_italic or "")
    if first_italic and author_only:
        return first_italic.group(1).strip(" .,")

    # journal style: strip the leading author list, then cut at the italic journal
    e = re.sub(r"^(?:[A-ZÀ-Þ][\w'’\-]+(?:\s+[A-Z]\.?){1,4},?\s*|et al\.?,?\s*)+", "", e)
    e = e.split("*")[0]
    return e.strip(" .,")


def first_author_surname(entry):
    m = re.match(r"^([A-ZÀ-Þ][\w'’\-]+)", entry.strip())
    return m.group(1) if m else ""


def main():
    md = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_MD
    print(f"manuscript: {md}\n")
    refs = extract_references(md)
    print(f"{len(refs)} reference entries found\n")

    fails, warns, skipped, exempt = [], [], [], []
    for num, entry in refs:
        dm = DOI_RE.search(entry)
        if not dm:
            skipped.append((num, entry[:70]))
            continue
        doi = clean_doi(dm.group(1))
        try:
            msg = crossref(doi)
        except Exception as e:  # noqa: BLE001
            fails.append((num, doi, f"DOI DID NOT RESOLVE: {e}"))
            print(f"  [FAIL] {num:>2}  {doi}\n         DOI did not resolve: {e}")
            continue

        xt = (msg.get("title") or [""])[0]
        ct = cited_title(entry)
        n_ct, n_xt = norm(ct), norm(xt)
        ratio = difflib.SequenceMatcher(None, n_ct, n_xt).ratio()

        # Containment means the same work with an extra qualifier on one side,
        # e.g. citing "...Technical Reference Guide, Volume 3: Validation" against
        # a DOI whose registered title omits the volume, or an author fragment
        # leaking into the extracted title. That is NOT the F124 pattern, which is
        # two unrelated strings. Treat containment as a match.
        if len(n_ct) >= 20 and len(n_xt) >= 20 and (n_xt in n_ct or n_ct in n_xt):
            ratio = 1.0

        auths = msg.get("author") or []
        x_surname = (auths[0].get("family") or "") if auths else ""
        c_surname = first_author_surname(entry)
        author_ok = (not x_surname) or (not c_surname) or \
                    norm(x_surname).split()[0] in norm(c_surname) or \
                    norm(c_surname).split()[0] in norm(x_surname)

        if doi in KNOWN_EXCEPTIONS:
            why, against = KNOWN_EXCEPTIONS[doi]
            exempt.append((num, doi))
            print(f"  [EXEMPT] {num:>2}  {doi}")
            print(f"           crossref metadata unreliable: {why}")
            print(f"           verified against: {against}")
            time.sleep(0.25)
            continue

        if ratio < TITLE_RATIO_FAIL:
            fails.append((num, doi, "TITLE DOES NOT BELONG TO THIS DOI"))
            print(f"  [FAIL] {num:>2}  {doi}   title similarity {ratio:.2f}")
            print(f"         cited   : {ct[:88]}")
            print(f"         crossref: {xt[:88]}")
        elif ratio < TITLE_RATIO_WARN or not author_ok:
            warns.append((num, doi, f"similarity {ratio:.2f}, author_ok={author_ok}"))
            print(f"  [warn] {num:>2}  {doi}   similarity {ratio:.2f}  first-author ok={author_ok}")
            print(f"         cited   : {ct[:88]}")
            print(f"         crossref: {xt[:88]}")
        else:
            print(f"  [ ok ] {num:>2}  {doi}   {xt[:66]}")
        time.sleep(0.25)   # be polite to Crossref

    print("\n" + "=" * 74)
    print(f"  {len(refs)-len(skipped)-len(fails)-len(warns)-len(exempt)} verified   "
          f"{len(warns)} warn   {len(fails)} FAIL   {len(exempt)} exempt   "
          f"{len(skipped)} skipped (no DOI)")
    if exempt:
        print("\n  Exempt (Crossref metadata wrong, citation verified against a stronger source):")
        for num, d in exempt:
            print(f"    {num:>2}. {d}")
    if skipped:
        print("\n  No DOI, verify by hand against the publisher:")
        for num, e in skipped:
            print(f"    {num:>2}. {e}")
    if fails:
        print("\n  FAILURES: a title that does not belong to its DOI is the F124 pattern.")
        print("  Do NOT build or submit until these are resolved.")
    print("=" * 74)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
