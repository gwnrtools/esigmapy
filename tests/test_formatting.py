"""Enforce black formatting on the entire source tree.

Running `black --check` here means any unformatted commit will be caught by
`pytest` (and by CI) before it is merged.  To fix a failure, run:

    black esigmapy/ tests/
"""

import subprocess
import sys
from pathlib import Path

_ROOT = Path(__file__).parent.parent


def test_black_formatting():
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "black",
            "--check",
            str(_ROOT / "esigmapy"),
            str(_ROOT / "tests"),
        ],
        capture_output=True,
    )
    assert result.returncode == 0, (
        "black found unformatted files:\n" + result.stderr.decode()
    )
