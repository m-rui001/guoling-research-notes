import glob, sys
LQ = '\u201c'
RQ = '\u201d'


def is_cjkish(ch):
    if not ch:
        return False
    o = ord(ch)
    return (0x4E00 <= o <= 0x9FFF) or (0x3400 <= o <= 0x4DBF) or (0x3000 <= o <= 0x303F) or (0xFF00 <= o <= 0xFFEF)


def repair(s):
    out = []
    for i, ch in enumerate(s):
        prev = s[i - 1] if i > 0 else ''
        nxt = s[i + 1] if i + 1 < len(s) else ''
        if ch == LQ and prev in '([':
            out.append('"')
            continue
        if s.startswith('(""', i-1):
            out.append('"')
            continue
            out.append('"')
            continue
        if ch == RQ and nxt in '),]':
            out.append('"')
            continue
        out.append(ch)
    return ''.join(out)


def fix(s):
    out = []
    in_math = False
    flip = True
    n = len(s)
    for i, ch in enumerate(s):
        if ch == '$':
            in_math = not in_math
            out.append(ch)
            continue
        if in_math:
            out.append(ch)
            continue
        if ch == '"':
            prev = s[i - 1] if i > 0 else ''
            j = i + 1
            while j < n and s[j] == ' ':
                j += 1
            nn = s[j] if j < n else ''
            if (prev in '([,') or (nn in '),]:'):
                out.append('"')
                flip = True
                continue
            if is_cjkish(prev) or prev in (LQ, RQ) or (prev == '"' and is_cjkish(nn)):
                out.append(LQ if flip else RQ)
                flip = not flip
                continue
        out.append(ch)
    return ''.join(out)


targets = sys.argv[1:] or sorted(glob.glob('parts/*.typ'))
for f in targets:
    s = open(f, encoding='utf-8').read()
    t = fix(repair(s))
    if t != s:
        open(f, 'w', encoding='utf-8').write(t)
        print('fixed', f)
