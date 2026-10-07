"""Concatenate explicitly tracked inputs without a shell argument-length limit."""
import gzip
from pathlib import Path

with open(snakemake.output.fq, "wb") as out, \
     open(snakemake.output.manifest, "w") as manifest, \
     open(snakemake.log[0], "w") as log:
    manifest.write("path\tsize_bytes\tmtime_ns\n")
    for name in snakemake.input:
        path = Path(name)
        stat = path.stat()
        manifest.write(f"{path}\t{stat.st_size}\t{stat.st_mtime_ns}\n")
        log.write(f"Reading {path}\n")
        opener = gzip.open if path.name.endswith(".gz") else open
        last = b""
        with opener(path, "rb") as source:
            while chunk := source.read(1024 * 1024):
                out.write(chunk)
                last = chunk[-1:]
        # Prevent the last quality line joining the next file's record header.
        if last and last != b"\n":
            out.write(b"\n")
    log.write(f"Concatenated {len(snakemake.input)} input files\n")
