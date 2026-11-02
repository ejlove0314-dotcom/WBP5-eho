#!/usr/bin/env python3
"""Check every claim this repository and the paper make about this repository.

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
    r'|is the source of|sole\b|solely\b|entirely from|exclusively)\b', re.I)

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


def sentences(text):
    return [s.strip() for s in re.split(r'(?<=[.!?])\s+', text) if s.strip()]


def git_root(start):
    r = subprocess.run(['git', '-C', start, 'rev-parse', '--show-toplevel'],
                       capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else start


def read(path):
    return io.open(path, encoding='utf-8', errors='replace').read()


def main(argv):
    docs = [a for a in argv if a.lower().endswith('.docx')]
    repo = git_root(os.path.dirname(os.path.abspath(__file__)))

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

    def hit(t):
        if t in files or t in dirs:
            return True
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

    sources = []
    for d in docs:
        if not os.path.exists(d):
            print('  (skipped, not found: %s)' % d)
            continue
        tag = os.path.basename(d)[:28]
        sources += [(tag, i, p) for i, p in enumerate(docx_paragraphs(d))]
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
            if SKIP.match(tok) or '/' not in tok:
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

    # ---------------- C ----------------
    print()
    print('=' * 74)
    print('C. Universal and exclusive claims about the deposit  (read these)')
    print('=' * 74)
    n = 0
    for tag, i, text in sources:
        for s in sentences(text):
            if QUANT.search(s) and (DEPOSIT_WORD.search(s) or path_re.search(s)):
                n += 1
                print('  [%s %s] %s' % (tag, i, s[:280]))
    print()
    print('  %d sentence(s) in this class. Each is true only of the deposit as'
          % n)
    print('  it now stands. Changing the deposit means re-reading this list.')

    return 1 if (bad or failed_b) else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
