#!/usr/bin/env python3
"""Check the statements this repository and the paper make about this deposit.

Known limitation: a claim spread across two sentences -- one naming a file,
the next asserting it is deposited -- is not detected.  Widening the regular
expressions to reach across sentence boundaries costs more in false positives
than it recovers, so it is recorded here rather than papered over.

Scope: four checks over the repository's own .md files, extended to the
manuscript, Supplementary Information and response letter when those are
given on the command line.  It does not check the science, the numbers in
the tables, or prose anywhere else.

Statements of the form "the repository contains X", "all values derive from Y",
"twenty-two of the fifty-five scripts" are true when written and can stop being
true when the deposit changes under them.  During revision this class produced
a false sentence six times running, twice in the very edit that fixed the
previous one, so it is checked here mechanically rather than by reading.

Three checks:

  A. Path references.  Every path-like token in the repository's own .md files
     -- and in the manuscript, Supplementary Information and response letter,
     when those are given on the command line -- must name something that
     exists here.  Paths written inside a README are resolved relative to that
     README's own folder, which is how they are meant to be read.

  B. Self-counts.  The sentences in README.md that count this repository's own
     scripts are parsed out of the prose and recomputed from the files.  The
     expected values are not hard-coded: they are read from the README, so the
     check compares the documentation against the deposit rather than against
     a constant that can itself go stale.

  D. Cross-references.  Every "Table N", "Fig. N" and "Section N.N" reference in
     this documentation must name something that exists in the manuscript or the
     Supplementary Information.  This is the form the rot actually took during
     revision -- the structural material moved into the Supplementary Information
     and several folder READMEs went on citing the submitted version's table
     numbers -- and check A does not see it, because a table number is not a path.
     D needs the documents, so it is skipped when they are not given.

  C. Universal and exclusive claims.  Sentences carrying "all ...", "only",
     "nothing here", "derive from" and the like are listed in full.  Whether
     each is true depends on what it means, not on what it names, so no
     program can settle them.  The value of the list is that it is exhaustive:
     changing the deposit means re-reading these sentences rather than
     remembering which ones existed.

Usage
-----
    python3 verification/check_deposit_claims.py
    python3 verification/check_deposit_claims.py manuscript.docx si.docx rl.docx

With no arguments it checks the repository's own documentation, which is
everything a clone contains.  The three manuscript files are not deposited --
they are the files under review -- so pass them explicitly to extend checks A
and C over the paper itself.  The run says which sources it used.

Exit status is 0 when A and B pass and non-zero otherwise.  C never fails a
run; it prints a list for a person to read.

Requires only the Python standard library.
"""
import hashlib
import io
import os
import re
import subprocess
import sys
import zipfile

# --------------------------------------------------------------------------
# Reading .docx without a dependency: paragraph text, with deleted runs from
# tracked changes removed so that what is checked is what the reader sees.
# --------------------------------------------------------------------------
_WP = re.compile(r'<w:p(?:\s[^>]*)?>.*?</w:p>|<w:p(?:\s[^>]*)?/>', re.S)
_WT = re.compile(r'<w:t(?:\s[^>]*)?>(.*?)</w:t>', re.S)
_WDEL = re.compile(r'<w:delText(?:\s[^>]*)?>.*?</w:delText>', re.S)


def docx_paragraphs(path):
    with zipfile.ZipFile(path) as z:
        xml = z.read('word/document.xml').decode('utf-8')
    out = []
    for para in _WP.findall(xml):
        text = ''.join(_WT.findall(_WDEL.sub('', para)))
        out.append(text.replace('&amp;', '&').replace('&lt;', '<')
                       .replace('&gt;', '>').replace('&quot;', '"')
                       .replace('&apos;', "'"))
    return out


# --------------------------------------------------------------------------
# What counts as a path.  Built from the repository itself, so the matcher
# cannot drift from what is actually deposited: a token qualifies if it
# carries one of the deposit's file extensions, or is rooted at one of its
# top-level directories.  This keeps gene pairs (YAP/TAZ), units (mg/kg) and
# ratios (1/1,001) out, which a bare "contains a slash" rule does not.
# --------------------------------------------------------------------------
EXTS = ('csv tsv py md txt toml yml json ipynb smi sdf pdb pdbqt xlsx docx '
        'xml gz jpg png ccp4 h5ad prior').split()
SKIP = re.compile(r'^(https?|www|doi|10\.\d)', re.I)
RECORD = None

# A path with a directory in it is a claim about where something lives here, so
# it is always resolved.  A bare filename is ambiguous -- most name an input
# the reader supplies or an output a script writes at run time -- so it is
# resolved only when the sentence around it says the file is held here, and not
# when the sentence marks it as an input or a run-time product.
CLAIMS_HELD = re.compile(
    r'\b(?:are|is|were|was) (?:recorded|deposited|held|stored|given|listed|'
    r'tabulated|contained|provided|supplied|included) (?:in|with|as)\b'
    r'|\bdeposited (?:as|in|here|with)\b'
    r'|\baccompan(?:y|ies|ied) this (?:revision|submission)\b'
    r'|\b(?:ships|comes) with\b'
    r'|\b(?:this |the )?repository (?:holds|contains)\b', re.I)
NOT_A_CLAIM = re.compile(
    r'\b(?:Outputs?|Inputs?|Requires?|Usage)\s*:|available from|downloaded'
    r'|\bwrites\b|\bproduced by\b|\bgenerated by\b|at run time'
    r'|\bplaced in\b|\bregenerated by\b|\bexports?\b', re.I)


def claims_deposited(sent, pos):
    """Does `sent` say the name at offset `pos` is held in this deposit?

    Both patterns can match the same sentence -- "is deposited in X and the
    tool exports it" is a true deposit claim with an output verb trailing it --
    so an unconditional veto would discard true claims.  Whichever pattern
    matches nearer the name wins, which is the clause the name belongs to.
    """
    def nearest(rx):
        best = None
        for m in rx.finditer(sent):
            d = min(abs(m.start() - pos), abs(m.end() - pos))
            if best is None or d < best:
                best = d
        return best

    claim, veto = nearest(CLAIMS_HELD), nearest(NOT_A_CLAIM)
    if claim is None:
        return False
    return veto is None or claim <= veto

# Tokens that cannot resolve for reasons that are not defects, each with the
# reason.  A checker that keeps reporting known-good findings stops being read.
KNOWN = {
    'reproducibility/ProTox-3.0':
        'filename contains spaces; the table row names a glob matching the '
        'six ProTox workbooks, and a whitespace-delimited token cannot span it',
    'out/obs.csv':
        'usage example in phase1_singlecell/README.md; out/ is created by the '
        'script at run time and is not a deposited path',
}

QUANT = re.compile(
    r'\b(all (?:values|of the|ten|five|six|four|three|two|the)'
    r'|every (?:value|file|script|residue|one|compound|number)'
    r'|the only\b|nothing (?:here|in this)|not used for any'
    r'|no (?:value|file|script) |are not the basis of|derive[sd]? from'
    r'|is the source of|sole\b|solely\b|entirely from|exclusively'
    # superlatives: "the largest by volume", "the poorest profile",
    # "second-highest" -- claims of rank, true only of the current set
    r'|the (?:largest|smallest|deepest|highest|lowest|poorest|best|worst|'
    r'most|least|first|second|third)\b|(?:second|third)-(?:highest|lowest|best)'
    # negatives with no quantifier: "is not redistributed", "are not deposited"
    r'|(?:is|are|was|were) not (?:redistributed|deposited|included|retained|used)'
    # reporting verbs: where a value is said to live
    r'|(?:tabulated|listed|given|reported|recorded|shown) (?:in|as|here|above|below))\b',
    re.I)

DEPOSIT_WORD = re.compile(
    r'\b(deposit|repositor|archiv|accompan|available at|provided (?:in|with)|'
    r'contains?|holds?|carr(?:y|ies)|records?|stored|supplied)\b', re.I)

WORDNUM = {w: n for n, w in enumerate(
    'zero one two three four five six seven eight nine ten eleven twelve '
    'thirteen fourteen fifteen sixteen seventeen eighteen nineteen'.split())}
WORDNUM.update({'twenty': 20, 'thirty': 30, 'forty': 40, 'fifty': 50,
                'sixty': 60, 'seventy': 70, 'eighty': 80, 'ninety': 90})


def wordnum(s):
    """'fifty-five' -> 55, 'twelve' -> 12, '22' -> 22; None if not a number."""
    s = s.strip().lower().replace('–', '-')
    if s.replace(',', '').isdigit():
        return int(s.replace(',', ''))
    total = 0
    for part in s.split('-'):
        if part not in WORDNUM:
            return None
        total += WORDNUM[part]
    return total


def sentence_at(text, pos):
    """The sentence of `text` containing character offset `pos`."""
    start = 0
    for m in re.finditer(r'(?<=[.!?;])\s+', text):
        if m.end() > pos:
            break
        start = m.end()
    m = re.search(r'[.!?;]\s', text[pos:])
    end = pos + m.start() + 1 if m else len(text)
    return text[start:end]


def sentences(text):
    return [s.strip() for s in re.split(r'(?<=[.!?])\s+', text) if s.strip()]


def git_root(start):
    r = subprocess.run(['git', '-C', start, 'rev-parse', '--show-toplevel'],
                       capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else start


def read(path):
    return io.open(path, encoding='utf-8', errors='replace').read()


def main(argv):
    global RECORD
    docs = [a for a in argv if a.lower().endswith('.docx')]
    repo = git_root(os.path.dirname(os.path.abspath(__file__)))
    RECORD = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'claims_adjudicated.tsv')

    files = set(f for f in subprocess.run(
        ['git', '-C', repo, 'ls-files'], capture_output=True,
        text=True).stdout.split('\n') if f)
    if not files:
        print('No git-tracked files found under %s' % repo)
        return 2
    dirs = set()
    for f in files:
        parts = f.split('/')
        for i in range(1, len(parts)):
            dirs.add('/'.join(parts[:i]))

    tops = sorted({d.split('/')[0] for d in dirs}, key=len, reverse=True)
    path_re = re.compile(
        r'(?<![\w/.-])((?:[A-Za-z0-9_.*-]+/)*[A-Za-z0-9_.*-]+\.(?:%s)\b'
        r'|(?:%s)(?:/[A-Za-z0-9_.*-]+)*/?)'
        % ('|'.join(EXTS), '|'.join(re.escape(t) for t in tops)))

    basenames = {os.path.basename(f) for f in files}

    def hit(t):
        if t in files or t in dirs:
            return True
        # A bare filename with no directory resolves if the deposit holds a
        # file of that name anywhere.  Skipping these is how a README came to
        # name a run-time output as though it were deposited.
        if '/' not in t:
            if t in basenames:
                return True
            if '*' in t:
                pre, suf = t.split('*')[0], t.rsplit('*', 1)[1]
                return any(b.startswith(pre) and b.endswith(suf)
                           for b in basenames)
        if t.endswith('/'):
            return t.rstrip('/') in dirs
        if '*' in t:
            pre, suf = t.split('*')[0], t.rsplit('*', 1)[1]
            return any(f.startswith(pre) and f.endswith(suf) for f in files)
        return False

    def resolves(tok, base):
        t = tok.rstrip('.,;:)')
        if hit(t):
            return True
        return bool(base) and hit(
            os.path.normpath(os.path.join(base, t)).replace('\\', '/'))

    sources, doc_paras = [], []
    for d in docs:
        if not os.path.exists(d):
            print('  (skipped, not found: %s)' % d)
            continue
        tag = os.path.basename(d)[:28]
        doc_paras += [(tag, i, p) for i, p in enumerate(docx_paragraphs(d))]
    sources += doc_paras
    md = sorted(f for f in files if f.endswith('.md'))
    for f in md:
        sources += [('repo:' + f, i, line) for i, line
                    in enumerate(read(os.path.join(repo, f)).split('\n'), 1)]

    print('Repository : %s' % repo)
    print('Sources    : %d repository .md files%s' %
          (len(md), ', plus %d document(s) given on the command line' % len(docs)
           if docs else ' (no manuscript files given; A and C cover the '
                        'repository only)'))
    print()

    # ---------------- A ----------------
    print('=' * 74)
    print('A. Path references')
    print('=' * 74)
    seen, bad, known, total = set(), [], [], 0
    for tag, i, text in sources:
        for m in path_re.finditer(text):
            tok = m.group(1)
            if SKIP.match(tok):
                continue
            if '/' not in tok:
                # Judge on the sentence the name sits in, not the paragraph:
                # a Methods paragraph can deposit five things in its closing
                # sentence while naming an undeposited input earlier on.
                sent = sentence_at(text, m.start())
                off = m.start() - text.find(sent) if sent in text else 0
                if not claims_deposited(sent, max(off, 0)):
                    continue
            total += 1
            if (tag, tok) in seen:
                continue
            seen.add((tag, tok))
            base = os.path.dirname(tag[5:]) if tag.startswith('repo:') else ''
            if not resolves(tok, base):
                entry = (tag, i, tok, text[max(0, m.start() - 55):m.start() + 55])
                (known if tok in KNOWN else bad).append(entry)
    print('  %d tokens, %d distinct (source, path) pairs' % (total, len(seen)))
    print('  unresolved %d   known-and-explained %d' % (len(bad), len(known)))
    for tag, i, tok, ctx in bad:
        print('    FAIL  %-26s %-5s %-44s | %s'
              % (tag, i, tok, ctx.replace('\n', ' ')))
    for tag, i, tok, _ in known:
        print('    known %-26s %-5s %-30s %s' % (tag, i, tok, KNOWN[tok]))

    # ---------------- B ----------------
    print()
    print('=' * 74)
    print('B. Counts README.md states about this repository')
    print('=' * 74)
    me = os.path.relpath(os.path.abspath(__file__), repo).replace(os.sep, '/')
    py = sorted(f for f in files if f.endswith('.py') and f != me)
    excluded = me if me in files else None
    ABS = re.compile(r'["\']([A-Za-z]:[/\\]|/home/|/mnt/)|r["\'][A-Za-z]:')
    # Built from code points so that this file stays pure ASCII: a
    # literal Hangul range here would make the checker match itself.
    HANGUL = re.compile('[' + chr(0xAC00) + '-' + chr(0xD7A3) + ']')
    absn = [f for f in py if ABS.search(read(os.path.join(repo, f)))]
    kor = [f for f in py if HANGUL.search(read(os.path.join(repo, f)))]
    measured = {
        'scripts': len(py),
        'absolute paths': len(absn),
        'of those in deprecated/': len([f for f in absn
                                        if f.startswith('deprecated/')]),
        'Korean comments': len(kor),
        'of those in revision_outputs/scripts/':
            len([f for f in kor if f.startswith('revision_outputs/scripts/')]),
    }

    rm = read(os.path.join(repo, 'README.md'))
    NUM = r'([A-Za-z]+(?:-[A-Za-z]+)?|\d+)'
    stated = {}
    m = re.search(NUM + r' of the ' + NUM + r' scripts carry absolute paths', rm)
    if m:
        stated['absolute paths'] = wordnum(m.group(1))
        stated['scripts'] = wordnum(m.group(2))
    m = re.search(NUM + r' of those are in `deprecated/`', rm)
    if m:
        stated['of those in deprecated/'] = wordnum(m.group(1))
    m = re.search(NUM + r' scripts carry comments in Korean', rm)
    if m:
        stated['Korean comments'] = wordnum(m.group(1))
    m = re.search(r'comments in Korean, most of them\s+in `revision_outputs/scripts/`', rm)
    if m:
        # "most of them" is a claim about a majority, not an exact count
        stated['of those in revision_outputs/scripts/'] = 'majority'

    failed_b = False
    if not stated:
        print('  No countable sentences found in README.md -- check the '
              'patterns in this script against the current wording.')
        failed_b = True
    for key, want in sorted(stated.items()):
        got = measured[key]
        if want == 'majority':
            ok = got * 2 > measured['Korean comments']
            print('  %-38s stated "most"  measured %-4s %s'
                  % (key, '%d of %d' % (got, measured['Korean comments']),
                     'OK' if ok else '<<< NOT A MAJORITY'))
        else:
            ok = (got == want)
            print('  %-38s stated %-6s measured %-6s %s'
                  % (key, want, got, 'OK' if ok else '<<< MISMATCH'))
        failed_b = failed_b or not ok
    print('  (for reference: %d git-tracked files in this deposit)' % len(files))
    if excluded:
        print('  (%s is excluded from the script counts: it produced no\n   deposited output, so README.md does not count it)' % excluded)

    # ---------------- D ----------------
    print()
    print('=' * 74)
    print('D. Cross-references to tables, figures and sections')
    print('=' * 74)
    failed_d = False
    if not doc_paras:
        print('  Skipped: this check needs the manuscript and Supplementary')
        print('  Information, which are not deposited here. Pass them on the')
        print('  command line to run it.')
    else:
        # The valid targets are whatever the documents actually define, read
        # from their own numbered captions -- not a list kept in this script,
        # which could go stale in exactly the way the check exists to catch.
        defined = set()
        for tag, _i, para in doc_paras:
            m = re.match(r'\s*Table (S?\d+)[.．]', para)
            if m:
                defined.add('Table ' + m.group(1))
            m = re.match(r'\s*(?:Figure|Fig\.) (S?\d+)[.．]', para)
            if m:
                defined.add('Fig. ' + m.group(1))
            m = re.match(r'\s*(\d+\.\d+)\s+\S', para)
            if m:
                defined.add('Section ' + m.group(1))
        # The keyword, its first number, and any list or range that follows.
        # Without the fourth group "Tables 3 and 4" is checked as Table 3
        # alone and "Sections 3.4-3.7" as Section 3.4 alone, which would let
        # through half of the fourteenth-round finding this check exists for.
        REF = re.compile(r'\b(?:(Supplementary|Supp\.)\s+)?'
                         r'(Tables?|Figs?\.|Figures?|Sections?)\s+'
                         r'(S?\d+(?:\.\d+)?)'
                         r'((?:\s*(?:,|and|to|\u2013|\u2014|-)\s*'
                         r'S?\d+(?:\.\d+)?)*)')
        DASH = re.compile(r'[\u2013\u2014-]')

        def numbers_in(m):
            """Every number the reference names, ranges expanded."""
            nums = [m.group(3)] + re.findall(r'S?\d+(?:\.\d+)?', m.group(4) or '')
            # A dash between two plain integers is a range and every number in
            # it is referenced.  A dash between dotted section numbers is left
            # as its endpoints: 3.4-3.7 does not require 3.5 and 3.6 to be
            # captions in their own right.
            if (len(nums) == 2 and DASH.search(m.group(4) or '')
                    and all('.' not in n for n in nums)
                    and all(n.lstrip('S').isdigit() for n in nums)):
                a, b = nums
                pre = 'S' if a.startswith('S') and b.startswith('S') else ''
                lo, hi = int(a.lstrip('S')), int(b.lstrip('S'))
                if 0 <= hi - lo <= 30:
                    nums = ['%s%d' % (pre, v) for v in range(lo, hi + 1)]
            return nums

        bad_d, known_d, nref = [], [], 0
        for tag, i, text in sources:
            for m in REF.finditer(text):
                kind = m.group(2)
                k = ('Table ' if kind.lower().startswith('table')
                     else 'Fig. ' if kind.lower().startswith('fig')
                     else 'Section ')
                for num in numbers_in(m):
                    # "Supplementary Table 8" means S8; section numbering is
                    # not prefixed that way.
                    if m.group(1) and not num.startswith('S') and k != 'Section ':
                        num = 'S' + num
                    nref += 1
                    if k + num in defined:
                        continue
                    # Three kinds of reference are meant not to resolve against
                    # the revised documents, and each says so where it stands.
                    why = None
                    if re.search(r'\bsubmitted\b', text, re.I):
                        why = ('names the submitted version explicitly, which '
                               'is what the sentence is about')
                    elif tag.startswith('repo:verification/'):
                        why = ('verification/ preserves the July-August 2026 '
                               'audit of the submitted manuscript; '
                               'README_verification_ledger_scope.md declares it')
                    elif re.match(r'\s*Sections?\s+\d+\s*$', text):
                        why = ("a heading of the response letter's own "
                               'numbering, not a reference into the manuscript')
                    (known_d if why else bad_d).append(
                        (tag, i, (k + num).strip(),
                         why or text[max(0, m.start() - 55):m.start() + 55]))
        print('  %d reference(s) checked against %d numbered captions and '
              'headings' % (nref, len(defined)))
        print('  unresolved %d   known-and-explained %d'
              % (len(bad_d), len(known_d)))
        for tag, i, ref, ctx in bad_d:
            print('    FAIL  %-26s %-5s %-24s | %s'
                  % (tag, i, ref, ctx.replace('\n', ' ')))
        for tag, i, ref, why in known_d:
            print('    known %-26s %-5s %-24s %s' % (tag, i, ref, why))
        failed_d = bool(bad_d)

    # ---------------- C ----------------
    print()
    print('=' * 74)
    print('C. Universal and exclusive claims about the deposit  (read these)')
    print('=' * 74)
    found = []
    for tag, i, text in sources:
        for sent in sentences(text):
            if QUANT.search(sent) and (DEPOSIT_WORD.search(sent)
                                       or path_re.search(sent)):
                found.append((hashlib.sha1(
                    ' '.join(sent.split()).encode('utf-8')).hexdigest()[:12],
                    tag, i, sent))

    # A record of which sentences have been read, so that a later run shows the
    # ones that changed rather than reprinting the whole list.  Without it the
    # check degrades into a wall of text that nobody reads, which is how this
    # class of error survived six rounds in the first place.
    seen_before = {}
    if os.path.exists(RECORD):
        for line in io.open(RECORD, encoding='utf-8'):
            if line.strip() and not line.startswith('#'):
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    seen_before[parts[0]] = parts[1]

    fresh = [f for f in found if f[0] not in seen_before]
    old_n = len(found) - len(fresh)
    print('  %d sentence(s) in this class; %d already recorded as read, '
          '%d new or changed.' % (len(found), old_n, len(fresh)))
    print()
    for h, tag, i, sent in (fresh if seen_before else found):
        print('  [%s %s] %s' % (tag, i, sent[:280]))
    if seen_before and not fresh:
        print('  (nothing new since the recorded run)')
    print()
    print('  Each is true only of the deposit as it now stands. Changing the')
    print('  deposit means re-reading these, not remembering them.')
    print('  Re-run with --record after reading to write %s.'
          % os.path.relpath(RECORD, repo).replace(os.sep, '/'))
    if '--record' in argv:
        sha = subprocess.run(['git', '-C', repo, 'rev-parse', '--short', 'HEAD'],
                             capture_output=True, text=True).stdout.strip()
        with io.open(RECORD, 'w', encoding='utf-8') as fh:
            fh.write('# Sentences making a universal, exclusive or superlative '
                     'claim about this deposit,\n# recorded as read at the '
                     'commit named on each line.  A sentence that is later\n'
                     '# reworded gets a new hash and reappears as new.\n')
            for h, tag, i, sent in found:
                fh.write('%s\t%s\t%s %s\t%s\n'
                         % (h, sha, tag, i, ' '.join(sent.split())[:200]))
        print('  recorded %d sentence(s) at %s' % (len(found), sha))

    return 1 if (bad or failed_b or failed_d) else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
