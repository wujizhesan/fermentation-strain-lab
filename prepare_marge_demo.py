import argparse
from http.client import HTTPException
import json
from pathlib import Path
import time
from urllib.request import Request, urlopen
import zipfile

from marge_sources import ARCHIVE_PATH, archive_manifest, digest, validate_archive_members, validate_marge_inputs


ROOT = Path(__file__).resolve().parent
APPENDIX_PATH = 'data/raw/ale_marge_20261012/MSB-17-e10253-s002.pdf'
APPENDIX_MEMBER = 'MSB-17-e10253-s002.pdf'


def download_archive(manifest, attempts=3, timeout=60):
    request = Request(manifest['source_url'], headers={'User-Agent': 'ale-marge-demo/1.1.0'})
    for attempt in range(1, attempts + 1):
        print(f'官方补材下载：第{attempt}/{attempts}次，每次读取超时{timeout}秒。', flush=True)
        try:
            chunks = []
            size = 0
            next_progress = 1024 * 1024
            with urlopen(request, timeout=timeout) as response:
                while True:
                    chunk = response.read1(1024 * 1024)
                    if not chunk:
                        break
                    chunks.append(chunk)
                    size += len(chunk)
                    if size >= next_progress:
                        print(f'已收到补材：{size / 1024 / 1024:.2f} MiB。', flush=True)
                        next_progress = size + 1024 * 1024
            body = b''.join(chunks)
        except (OSError, HTTPException) as error:
            if attempt == attempts:
                raise RuntimeError(f'官方补材下载连续{attempts}次失败；未写入任何输入。') from error
            print(f'下载中断：{type(error).__name__}；正在重新请求。', flush=True)
            time.sleep(attempt)
            continue
        validate_archive_members(body, manifest)
        return body
    raise ValueError('Download attempts must be positive')


def store_input(root, name, body):
    raw_root = (root / 'data/raw').resolve()
    destination = (root / name).resolve()
    if not destination.is_relative_to(raw_root):
        raise ValueError('Input target must remain in the raw-data directory')
    if destination.exists():
        raise ValueError('Refuse to replace an existing input: ' + name)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + '.download')
    if temporary.exists():
        raise ValueError('Download staging file already exists')
    temporary.write_bytes(body)
    temporary.replace(destination)


def prepare(root=ROOT, check_only=False):
    root = Path(root)
    protocol = json.loads((root / 'research/ale_marge_protocol.json').read_text(encoding='utf-8'))
    downloadable = {ARCHIVE_PATH, APPENDIX_PATH}
    missing_snapshots = [name for name in protocol['input_sha256'] if name not in downloadable and not (root / name).is_file()]
    if missing_snapshots:
        raise FileNotFoundError('Release snapshot incomplete: ' + '; '.join(missing_snapshots))
    if check_only:
        return validate_marge_inputs(root, protocol)
    manifest = archive_manifest(root)
    archive = root / ARCHIVE_PATH
    if not archive.exists():
        print('正在下载官方补材；随后核对27个冻结成员。', flush=True)
        body = download_archive(manifest)
        store_input(root, ARCHIVE_PATH, body)
    else:
        validate_archive_members(archive.read_bytes(), manifest)
    if not (root / APPENDIX_PATH).exists():
        with zipfile.ZipFile(archive) as bundle:
            body = bundle.read(APPENDIX_MEMBER)
        if digest(body) != protocol['input_sha256'][APPENDIX_PATH]:
            raise ValueError('Appendix content changed')
        store_input(root, APPENDIX_PATH, body)
    return validate_marge_inputs(root, protocol)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args()
    checks = prepare(check_only=args.check_only)
    print(json.dumps({'passed': True, 'input_count': len(checks), 'input_validation': checks}, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main()
