import hashlib
from http.client import IncompleteRead
import io
import unittest
from unittest.mock import patch
import zipfile

from prepare_marge_demo import download_archive


def supplement(content=b'reference'):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w') as archive:
        archive.writestr('table.csv', content)
    return buffer.getvalue()


def manifest():
    return {'source_url': 'https://example.invalid/supplementaryFiles',
        'members': {'table.csv': {'bytes': 9, 'sha256': hashlib.sha256(b'reference').hexdigest()}}}


class PrepareMargeTests(unittest.TestCase):
    def test_stream_timeout_retries_from_the_start_and_validates(self):
        class Interrupted(io.BytesIO):
            def read1(self, size=-1):
                raise TimeoutError('Interrupted stream')
        body = supplement()
        with patch('prepare_marge_demo.urlopen', side_effect=[Interrupted(b'partial'), io.BytesIO(body)]) as request:
            with patch('prepare_marge_demo.time.sleep') as sleep:
                self.assertEqual(download_archive(manifest()), body)
        self.assertEqual(request.call_count, 2)
        sleep.assert_called_once_with(1)

    def test_chunked_disconnect_is_retried(self):
        body = supplement()
        with patch('prepare_marge_demo.urlopen', side_effect=[IncompleteRead(b'partial'), io.BytesIO(body)]):
            with patch('prepare_marge_demo.time.sleep'):
                self.assertEqual(download_archive(manifest()), body)

    def test_network_failure_has_a_bounded_attempt_count(self):
        with patch('prepare_marge_demo.urlopen', side_effect=TimeoutError('Unavailable')) as request:
            with patch('prepare_marge_demo.time.sleep'):
                with self.assertRaisesRegex(RuntimeError, '连续3次失败'):
                    download_archive(manifest())
        self.assertEqual(request.call_count, 3)

    def test_changed_content_is_rejected_without_network_retry(self):
        with patch('prepare_marge_demo.urlopen', return_value=io.BytesIO(supplement(b'changed'))) as request:
            with self.assertRaisesRegex(ValueError, 'member content changed'):
                download_archive(manifest())
        self.assertEqual(request.call_count, 1)


if __name__ == '__main__':
    unittest.main()
