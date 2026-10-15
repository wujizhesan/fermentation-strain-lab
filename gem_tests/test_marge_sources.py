import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import warnings
import zipfile

from marge_sources import ARCHIVE_PATH, MEMBER_MANIFEST, archive_manifest, validate_archive_members, validate_marge_inputs


ROOT = Path(__file__).resolve().parents[1]


def repack(records):
    body = io.BytesIO()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        with zipfile.ZipFile(body, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            for name, value in records:
                entry = zipfile.ZipInfo(name, date_time=(2026, 10, 13, 0, 0, 0))
                entry.compress_type = zipfile.ZIP_DEFLATED
                archive.writestr(entry, value)
    return body.getvalue()


class MargeSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = archive_manifest(ROOT)
        with zipfile.ZipFile(ROOT / ARCHIVE_PATH) as archive:
            cls.records = [(name, archive.read(name)) for name in archive.namelist()]

    def test_repackaged_archive_passes_the_actual_input_gate(self):
        body = repack(self.records)
        reference = self.manifest['reference_archive_sha256']
        self.assertNotEqual(hashlib.sha256(body).hexdigest(), reference)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = root / ARCHIVE_PATH; archive.parent.mkdir(parents=True); archive.write_bytes(body)
            manifest = root / MEMBER_MANIFEST; manifest.parent.mkdir(parents=True); manifest.write_bytes((ROOT / MEMBER_MANIFEST).read_bytes())
            receipt = validate_marge_inputs(root, {'input_sha256': {ARCHIVE_PATH: reference}})
        self.assertEqual(receipt[0]['member_count'], 27)
        self.assertEqual(receipt[0]['validation'], 'archive_member_sha256')
        self.assertNotEqual(receipt[0]['actual_sha256'], receipt[0]['reference_sha256'])

    def test_changed_member_is_rejected(self):
        records = list(self.records); name, value = records[0]; records[0] = (name, value + b'changed')
        with self.assertRaisesRegex(ValueError, 'member content changed'):
            validate_archive_members(repack(records), self.manifest)

    def test_missing_member_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'member coverage changed'):
            validate_archive_members(repack(self.records[:-1]), self.manifest)

    def test_duplicate_member_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate supplement member'):
            validate_archive_members(repack(self.records + [self.records[0]]), self.manifest)

    def test_member_manifest_cannot_be_resigned_by_editing_it(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); path = root / MEMBER_MANIFEST; path.parent.mkdir(parents=True)
            body = json.loads((ROOT / MEMBER_MANIFEST).read_text(encoding='utf-8'))
            body['members'][next(iter(body['members']))]['sha256'] = 'changed'
            path.write_text(json.dumps(body), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'member manifest changed'):
                archive_manifest(root)

    def test_nonarchive_inputs_keep_exact_byte_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); path = root / 'table.csv'; path.write_bytes(b'original')
            protocol = {'input_sha256': {'table.csv': hashlib.sha256(b'original').hexdigest()}}
            self.assertEqual(validate_marge_inputs(root, protocol)[0]['validation'], 'sha256')
            path.write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'Input fingerprint changed'):
                validate_marge_inputs(root, protocol)


if __name__ == '__main__':
    unittest.main()
