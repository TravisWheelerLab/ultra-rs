#!/usr/bin/env python3
"""Write the body of a GitHub release: a table of the downloads.

usage: release_table.py <owner/repo> <tag> <archive>...

Each archive is named <name>-<platform>.tar.gz or .zip. The table goes to
standard output, one row per archive, in the order of PLATFORMS.
"""

import os
import sys

# platform in the archive name -> (OS, architecture), in table order
PLATFORMS = {
    "macos-universal": ("macOS", "Universal (Intel and Apple Silicon)"),
    "x86_64-apple-darwin": ("macOS", "Intel 64 bit"),
    "aarch64-apple-darwin": ("macOS", "ARM 64 bit (Apple Silicon)"),
    "x86_64-unknown-linux-gnu": ("Linux", "Intel/AMD 64 bit"),
    "aarch64-unknown-linux-gnu": ("Linux", "ARM 64 bit"),
    "x86_64-pc-windows-msvc": ("Windows", "Intel/AMD 64 bit"),
}


def platform_of(archive: str) -> str:
    """'ultra-x86_64-apple-darwin.tar.gz' -> 'x86_64-apple-darwin'"""
    for ext in (".tar.gz", ".zip"):
        if archive.endswith(ext):
            return archive[: -len(ext)].split("-", 1)[1]
    sys.exit(f"{archive}: not a .tar.gz or .zip archive")


def main() -> None:
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    repo, tag = sys.argv[1], sys.argv[2]
    archives = {platform_of(a): a for a in map(os.path.basename, sys.argv[3:])}

    unknown = sorted(set(archives) - set(PLATFORMS))
    missing = [p for p in PLATFORMS if p not in archives]
    if unknown or missing:
        sys.exit(f"unexpected platforms: {unknown}; no archive for: {missing}")

    print("### Release Assets")
    print("| OS | Architecture | Link |")
    print("|---------|----------|-------------|")
    for platform, (system, arch) in PLATFORMS.items():
        url = f"https://github.com/{repo}/releases/download/{tag}/{archives[platform]}"
        print(f"| {system} | {arch} | [Download]({url}) |")

    print(
        """
Each download unpacks to a folder that holds two programs, `ultra` and `ultrafast`. Copy them to a folder on your `PATH`. The Linux programs need glibc 2.17 or newer.

On macOS, you may see the following message when you try to run `ultra` or `ultrafast`:

```
"ultra" cannot be opened because the developer cannot be verified.
```

To fix this, run the following command:

```
xattr -dr com.apple.quarantine <path to the unpacked folder>
```"""
    )


if __name__ == "__main__":
    main()
