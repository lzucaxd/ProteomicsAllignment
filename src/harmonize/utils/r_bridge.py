"""Subprocess bridge for calling R scripts from Python."""

from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def run_r_script(
    script: str | Path,
    args: list[str] | None = None,
    env: dict[str, str] | None = None,
    cwd: str | Path | None = None,
    timeout: int | None = 3600,
    vanilla: bool = True,
) -> subprocess.CompletedProcess:
    """
    Run an R script via Rscript subprocess.

    Parameters
    ----------
    script : path to .R file
    args : CLI arguments passed after the script path
    env : extra environment variables (merged with current env)
    cwd : working directory
    timeout : seconds before killing the process (None = no limit)
    vanilla : if True, pass --vanilla to skip .Rprofile / renv
    """
    script = Path(script).resolve()
    if not script.exists():
        raise FileNotFoundError(f"R script not found: {script}")

    cmd = ["Rscript"]
    if vanilla:
        cmd.append("--vanilla")
    cmd.append(str(script))
    if args:
        cmd.extend(args)

    run_env = {**os.environ}
    run_env["OMP_NUM_THREADS"] = run_env.get("OMP_NUM_THREADS", "1")
    if env:
        run_env.update(env)

    logger.info("Running: %s", " ".join(cmd))

    result = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        env=run_env,
        capture_output=True,
        text=True,
        timeout=timeout,
    )

    if result.stdout:
        for line in result.stdout.strip().split("\n"):
            logger.info("[R stdout] %s", line)
    if result.stderr:
        for line in result.stderr.strip().split("\n"):
            logger.warning("[R stderr] %s", line)

    if result.returncode != 0:
        raise RuntimeError(
            f"R script failed (exit {result.returncode}): {script}\n"
            f"stderr:\n{result.stderr[-2000:]}"
        )

    return result


def run_python_script(
    script: str | Path,
    args: list[str] | None = None,
    python_cmd: str = "python3",
    env: dict[str, str] | None = None,
    cwd: str | Path | None = None,
    timeout: int | None = 3600,
) -> subprocess.CompletedProcess:
    """Run a Python script via subprocess (for Celligner conda env etc.)."""
    script = Path(script).resolve()
    if not script.exists():
        raise FileNotFoundError(f"Python script not found: {script}")

    cmd = [python_cmd, str(script)]
    if args:
        cmd.extend(args)

    run_env = {**os.environ}
    if env:
        run_env.update(env)

    logger.info("Running: %s", " ".join(cmd))

    result = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        env=run_env,
        capture_output=True,
        text=True,
        timeout=timeout,
    )

    if result.stdout:
        for line in result.stdout.strip().split("\n"):
            logger.info("[py stdout] %s", line)
    if result.stderr:
        for line in result.stderr.strip().split("\n"):
            logger.warning("[py stderr] %s", line)

    if result.returncode != 0:
        raise RuntimeError(
            f"Python script failed (exit {result.returncode}): {script}\n"
            f"stderr:\n{result.stderr[-2000:]}"
        )

    return result
