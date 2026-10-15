from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PROTOCOL = ROOT / "research/ale_external_screen_protocol_20261015.json"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout", type=float, default=240)
    args = parser.parse_args()
    if args.timeout <= 0:
        raise ValueError("下载时限必须为正")
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    manifest_path = ROOT / protocol["source_manifest"]
    if hashlib.sha256(manifest_path.read_bytes()).hexdigest() != protocol["input_sha256"][protocol["source_manifest"]]:
        raise ValueError("来源清单与冻结 SHA 不一致")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    urls = {row["path"]: row["source_url"] for row in manifest["files"]}
    for pmc, metadata in manifest["metadata"].items():
        urls[metadata["path"]] = f"https://pmc-oa-opendata.s3.amazonaws.com/metadata/{pmc}.1.json"
        urls[f"data/raw/ale_external_sources_20261013/{pmc}.xml"] = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML"
    downloaded, existing = 0, 0
    for rel, url in urls.items():
        path = ROOT / rel
        expected = protocol["input_sha256"][rel]
        if path.exists():
            if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError(f"现有文件与冻结来源不一致，拒绝覆盖：{rel}")
            existing += 1
            continue
        print(f"下载：{rel}", flush=True)
        started = time.monotonic()
        blocks = []
        with urllib.request.urlopen(url, timeout=min(20, args.timeout)) as response:
            while True:
                block = response.read(65536)
                if not block:
                    break
                blocks.append(block)
                if time.monotonic() - started > args.timeout:
                    raise TimeoutError(f"下载超时，未写入原件：{rel}")
        content = b"".join(blocks)
        if hashlib.sha256(content).hexdigest() != expected:
            raise ValueError(f"服务内容与冻结 SHA 不一致，未写入：{rel}")
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(content)
        downloaded += 1
    print(f"来源核查通过：已有 {existing} 个，下载 {downloaded} 个；服务变化时保留原件并拒绝自动更新协议")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
