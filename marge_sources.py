import hashlib
import io
import json
from pathlib import Path
import zipfile


ARCHIVE_PATH = 'data/raw/ale_chassis_2021/all_supplementary.zip'
MEMBER_MANIFEST = 'research/ale_marge_archive_members.json'
MEMBER_MANIFEST_SHA256 = '517455cb82b002f26eaf72fc9060e564064cbc8ce5c4f877f681991afda5de3c'


def digest(body):
    return hashlib.sha256(body).hexdigest()


def archive_manifest(root):
    body = (Path(root) / MEMBER_MANIFEST).read_bytes()
    if digest(body) != MEMBER_MANIFEST_SHA256:
        raise ValueError('Frozen supplement member manifest changed')
    return json.loads(body)


def validate_archive_members(body, manifest):
    expected = manifest['members']
    with zipfile.ZipFile(io.BytesIO(body)) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError('Duplicate supplement member')
        if set(names) != set(expected):
            raise ValueError('Supplement member coverage changed')
        for name, reference in expected.items():
            member = archive.read(name)
            if len(member) != reference['bytes'] or digest(member) != reference['sha256']:
                raise ValueError('Supplement member content changed: ' + name)
    return len(expected)


def validate_marge_inputs(root, protocol):
    root = Path(root)
    checks = []
    for name, reference in protocol['input_sha256'].items():
        path = root / name
        if not path.is_file():
            raise FileNotFoundError('Missing frozen input: ' + name)
        body = path.read_bytes()
        actual = digest(body)
        if name == ARCHIVE_PATH:
            manifest = archive_manifest(root)
            if manifest['reference_archive_sha256'] != reference:
                raise ValueError('Supplement reference archive changed')
            count = validate_archive_members(body, manifest)
            checks.append({'path': name, 'reference_sha256': reference, 'actual_sha256': actual,
                'validation': 'archive_member_sha256', 'member_count': count,
                'member_manifest_sha256': MEMBER_MANIFEST_SHA256})
        else:
            if actual != reference:
                raise ValueError('Input fingerprint changed: ' + name)
            checks.append({'path': name, 'reference_sha256': reference, 'actual_sha256': actual,
                'validation': 'sha256'})
    return checks
