"""
Minimal HTTP server for the PRESTUS web alignment viewer.

Serves static files from --workdir and handles two JSON endpoints:
  GET  /placement_init.json   — initial transducer/focus positions
  POST /save_placement        — write placement_result.json and signal MATLAB
  POST /cancel_placement      — write cancelled result and signal MATLAB

Usage:
  python3 server.py --port 8742 --workdir /path/to/work_dir
"""

import argparse
import json
import os
import sys
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path


def make_handler(workdir: Path):

    class Handler(BaseHTTPRequestHandler):

        def log_message(self, fmt, *args):  # suppress access log noise
            pass

        def _send_json(self, data: dict, status: int = 200):
            body = json.dumps(data).encode()
            self.send_response(status)
            self.send_header('Content-Type', 'application/json')
            self.send_header('Content-Length', str(len(body)))
            self.send_header('Access-Control-Allow-Origin', '*')
            self.end_headers()
            self.wfile.write(body)

        def _send_file(self, path: Path):
            ext = path.suffix.lower()
            mime = {
                '.html': 'text/html',
                '.js':   'application/javascript',
                '.css':  'text/css',
                '.json': 'application/json',
                '.gz':   'application/gzip',
                '.nii':  'application/octet-stream',
            }.get(ext, 'application/octet-stream')

            data = path.read_bytes()
            self.send_response(200)
            self.send_header('Content-Type', mime)
            self.send_header('Content-Length', str(len(data)))
            self.send_header('Access-Control-Allow-Origin', '*')
            # Allow NiiVue to read the file
            self.send_header('Cross-Origin-Opener-Policy', 'same-origin')
            self.send_header('Cross-Origin-Embedder-Policy', 'require-corp')
            self.end_headers()
            self.wfile.write(data)

        def do_OPTIONS(self):
            self.send_response(204)
            self.send_header('Access-Control-Allow-Origin', '*')
            self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
            self.send_header('Access-Control-Allow-Headers', 'Content-Type')
            self.end_headers()

        def do_GET(self):
            path = self.path.split('?')[0].lstrip('/')
            if not path:
                path = 'index.html'

            target = (workdir / path).resolve()

            # Restrict to workdir only
            if not str(target).startswith(str(workdir.resolve())):
                self.send_response(403)
                self.end_headers()
                return

            if target.exists() and target.is_file():
                self._send_file(target)
            else:
                self.send_response(404)
                self.end_headers()

        def do_POST(self):
            length = int(self.headers.get('Content-Length', 0))
            body = self.rfile.read(length).decode() if length else '{}'

            if self.path == '/save_placement':
                result_path = workdir / 'placement_result.json'
                result_path.write_text(body)
                self._send_json({'status': 'ok'})

            elif self.path == '/cancel_placement':
                result_path = workdir / 'placement_result.json'
                result_path.write_text(json.dumps({'cancelled': True}))
                self._send_json({'status': 'cancelled'})

            else:
                self._send_json({'error': 'unknown endpoint'}, status=404)

    return Handler


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--port',    type=int,  default=8742)
    parser.add_argument('--workdir', type=str,  required=True)
    args = parser.parse_args()

    workdir = Path(args.workdir).resolve()
    if not workdir.is_dir():
        print(f'workdir does not exist: {workdir}', file=sys.stderr)
        sys.exit(1)

    server = HTTPServer(('localhost', args.port), make_handler(workdir))
    print(f'PRESTUS alignment server listening on http://localhost:{args.port}')
    print(f'  workdir: {workdir}')
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == '__main__':
    main()
