"""Exercise every client against HTTP fixtures, without starting Jenkins."""
from http.server import BaseHTTPRequestHandler, HTTPServer
import json
import os
from pathlib import Path
import subprocess
import tempfile
import threading
import unittest

CI = Path(__file__).resolve().parents[1]
CLIENTS = ('personal', 'team-workstation', 'perlmutter_gpu', 'stellar_cpu')


class SharedClientTests(unittest.TestCase):
    def test_queue_transitions_and_credentials(self):
        queue = {'why': 'Waiting for an executor | still queued'}
        requests = []

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def do_GET(self):
                requests.append(self.path)
                self.send_response(200)
                self.end_headers()
                if '/queue/item/' in self.path:
                    payload = queue
                elif '/api/json' in self.path:
                    payload = {'building': False, 'result': 'SUCCESS'}
                else:
                    payload = {}
                self.wfile.write(json.dumps(payload).encode())

        with tempfile.TemporaryDirectory() as directory, HTTPServer(('127.0.0.1', 0), Handler) as server:
            directory = Path(directory)
            auth = directory / 'auth'
            auth.write_text('test:fixture-token\n')
            auth.chmod(0o600)
            # Cluster clients retain their controller startup behavior; fake only
            # the tmux session and Java presence, and answer readiness via HTTP.
            (directory / 'tmux').write_text('#!/bin/sh\nexit 0\n')
            (directory / 'tmux').chmod(0o755)
            (directory / 'bin').mkdir()
            (directory / 'bin/java').write_text('#!/bin/sh\nexit 0\n')
            (directory / 'bin/java').chmod(0o755)
            (directory / 'jenkins.war').touch()
            env = dict(os.environ, JENKINS_URL=f'http://127.0.0.1:{server.server_port}',
                       JENKINS_CLI_AUTH_FILE=str(auth), GKEYLL_CI_ROOT=str(directory),
                       JAVA_HOME=str(directory), PATH=str(directory) + ':' + os.environ['PATH'],
                       JENKINS_JOB='fixture', TMPDIR=str(directory))
            worker = threading.Thread(target=server.serve_forever, daemon=True)
            worker.start()
            try:
                for client in CLIENTS:
                    command = [str(CI / f'jenkins-{client}.sh'), 'status', '--queue', '42']
                    with self.subTest(client=client):
                        queue.clear()
                        queue['why'] = 'Waiting for an executor | still queued'
                        result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=10)
                        self.assertEqual(result.returncode, 0, result.stderr)
                        self.assertIn('Queue item 42: Waiting for an executor | still queued', result.stdout)
                        queue.update(cancelled=True)
                        result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=10)
                        self.assertNotEqual(result.returncode, 0)
                        self.assertIn('CANCELLED', result.stdout)
                        queue.clear()
                        queue['executable'] = {'number': 7}
                        result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=10)
                        self.assertEqual(result.returncode, 0, result.stderr)
                        self.assertIn('SUCCESS', result.stdout)
                        auth.chmod(0o644)
                        result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=10)
                        self.assertNotEqual(result.returncode, 0)
                        self.assertIn('mode 600', result.stderr)
                        auth.chmod(0o600)
                self.assertFalse(list(directory.glob('gkeyll-jenkins.??????')))
            finally:
                server.shutdown()
                worker.join(timeout=5)

    def test_help_is_available_without_credentials(self):
        for client in CLIENTS:
            for args in (['--help'], ['run', '--help'], ['artifact', '--help']):
                with self.subTest(client=client, args=args):
                    result = subprocess.run([str(CI / f'jenkins-{client}.sh')] + args,
                                            capture_output=True, text=True, timeout=10)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertIn(f'jenkins-{client}.sh', result.stdout)
