"""Staged outputs and parity reporting shared by the notebooks.

A notebook writes its results into a Stage; they reach their final paths only if every parity
check passes, so a failed parity check leaves no results behind.
"""
import os
import shutil
from pathlib import Path


class Stage:
    """Outputs written under a staging folder, moved into place by commit()."""

    def __init__(self, staging_dir):
        self.dir = Path(staging_dir)
        self._files = {}  # final path -> staged path

    def path(self, final):
        """Staged path for a final output file (same file name, so the suffix still picks the format)."""
        final = Path(final)
        if final not in self._files:
            self._files[final] = self.dir / str(len(self._files)) / final.name
            self._files[final].parent.mkdir(parents=True, exist_ok=True)
        return self._files[final]

    def stem(self, final_stem):
        """Staged stem for a figure written as <stem>.png and <stem>.pdf (see viz.save_fig)."""
        final_stem = Path(final_stem)
        png = final_stem.with_name(final_stem.name + ".png")
        if png not in self._files:
            folder = self.dir / str(len(self._files))
            folder.mkdir(parents=True, exist_ok=True)
            for ext in (".png", ".pdf"):
                self._files[final_stem.with_name(final_stem.name + ext)] = folder / (final_stem.name + ext)
        return str(self._files[png].parent / final_stem.name)

    def commit(self):
        """Move every staged file to its final path and remove the staging folder."""
        missing = [str(f) for f, s in self._files.items() if not s.exists()]
        if missing:
            self.discard()
            raise FileNotFoundError("staged outputs were never written: " + ", ".join(sorted(missing)))
        for final, staged in sorted(self._files.items()):
            final.parent.mkdir(parents=True, exist_ok=True)
            os.replace(str(staged), str(final))
        self.discard()

    def discard(self):
        shutil.rmtree(str(self.dir), ignore_errors=True)
        self._files = {}


class Parity:
    """Collects named checks against legacy values; finish() writes the report and commits or raises."""

    def __init__(self):
        self.checks = []
        self.skipped = []  # (name, reason): checks that could not run; they do not count as passed

    def skip(self, name, reason):
        self.skipped.append((name, str(reason)))

    def check(self, name, ok, detail=""):
        self.checks.append((name, bool(ok), str(detail)))

    @property
    def ok(self):
        return all(ok for _, ok, _ in self.checks)

    def report(self):
        skipped = ", {} skipped".format(len(self.skipped)) if self.skipped else ""
        head = "parity: {} ({} checks{})".format("PASS" if self.ok else "FAIL", len(self.checks), skipped)
        lines = [head]
        for name, ok, detail in self.checks:
            lines.append("{}  {}{}".format("PASS" if ok else "FAIL", name, "  [" + detail + "]" if detail else ""))
        for name, reason in self.skipped:
            lines.append("SKIP  {}  [{}]".format(name, reason))
        return "\n".join(lines) + "\n"

    def finish(self, path, stage=None):
        """Write the report to path; commit the stage if all checks passed, else discard it and raise."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.report())
        if not self.ok:
            if stage is not None:
                stage.discard()
            raise AssertionError("parity failed, no results written:\n" + self.report())
        if stage is not None:
            stage.commit()
        print(self.report())
