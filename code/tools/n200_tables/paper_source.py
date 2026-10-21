"""Text read from the paper's LaTeX source (tables that build_tables.py and check_tables.py read, never write)."""

# The paper source marks edited text with one-argument macros that print their argument.
MARKS = ("ownerfix", "newfix", "orderfix", "releaseedit", "reviewedit", "redfix")


def unmark(s):
    """The text as printed: every \\name{X} of MARKS replaced by X (balanced braces, nested marks too)."""
    for name in MARKS:
        tok = "\\" + name + "{"
        while tok in s:
            i = s.index(tok)
            j, depth = i + len(tok), 1
            while depth:
                depth += {"{": 1, "}": -1}.get(s[j], 0)
                j += 1
            s = s[:i] + s[i + len(tok):j - 1] + s[j:]
    return s
