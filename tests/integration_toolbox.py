#!/usr/bin/env python

"""Integration tests for toolbox

Copyright 2019 SpacePy contributors
"""

try:
    import http.server
    requesthandlerclass = http.server.SimpleHTTPRequestHandler
except ImportError: #py2k
    import SimpleHTTPServer
    requesthandlerclass = SimpleHTTPServer.SimpleHTTPRequestHandler
    import SocketServer
import os
import shutil
import socket
import tempfile
import threading
import unittest

import spacepy_testing
import spacepy.toolbox


class SilentLoggingHTTPRequestHandler(requesthandlerclass):
    """HTTP request handler that doesn't log to stderr"""
    def log_message(self, *args):
        """Suppress the message"""
        pass


class WebGettingIntegration(unittest.TestCase):
    """Tests of functions that get data from the web"""

    def setUp(self):
        """Set up an HTTP server on a test directory"""
        self.td = tempfile.mkdtemp()
        """Root directory for HTTP server"""
        self.oldwd = os.getcwd()
        """Working directory before running this test"""
        #The server handles requests from the current directory
        try: #py3k
            serverclass = http.server.HTTPServer
        except NameError: #Py2k
            serverclass = SocketServer.TCPServer
        os.chdir(self.td)
        for port in range(8080, 8089):
            try:
                self.server = serverclass(
                    ('localhost', port), SilentLoggingHTTPRequestHandler)
            except (OSError, socket.error): #Port in use, try another
                continue
            self.port = port
            break
        else:
            raise RuntimeError('Unable to find a free port for http server.')
        self.server_thread = threading.Thread(target=self.server.serve_forever)
        self.server_thread.setDaemon(True)
        self.server_thread.start()

    def tearDown(self):
        """Quit HTTP server and remove test directory"""
        self.server.shutdown()
        self.server_thread.join()
        self.server.server_close()
        #This is a bit redundant, and doesn't necessarily do the job:
        #the socket will stay in WAIT for awhile regardless, seems to be
        #worse on Python 2
        self.server.socket.close()
        os.chdir(self.oldwd) #Back to old working directory
        shutil.rmtree(self.td)

    def testGetUrlReturn(self):
        """Call get_url, return data directly"""
        with open(os.path.join(self.td, 'foo.txt'), 'wb') as f:
            f.write(b'This is a test\n')
        data = spacepy.toolbox.get_url(
            'http://localhost:{}/foo.txt'.format(self.port))
        self.assertEqual(
            b"This is a test\n", data)

    def testGetUrlKeepaliveReturn(self):
        """Call get_url, return data directly, keep connection"""
        # This isn't the greatest test because the Python http server
        # doesn't support keepalive....
        with open(os.path.join(self.td, 'foo.txt'), 'wb') as f:
            f.write(b'This is a test\n')
        data, conn = spacepy.toolbox.get_url(
            'http://localhost:{}/foo.txt'.format(self.port),
            keepalive=True)
        try:
            self.assertEqual(
                b"This is a test\n", data)
            data, conn = spacepy.toolbox.get_url(
                'http://localhost:{}/foo.txt'.format(self.port),
                conn=conn, keepalive=True)
        finally:
            conn.close()
        self.assertEqual(
            b"This is a test\n", data)

    def testGetUrlToFile(self):
        """Call get_url, write to file"""
        with open(os.path.join(self.td, 'foo.txt'), 'wb') as f:
            f.write(b'This is a test\n')
        #Set its modification time to the distant past
        os.utime(os.path.join(self.td, 'foo.txt'),
                 (0, 0))
        data = spacepy.toolbox.get_url(
            'http://localhost:{}/foo.txt'.format(self.port),
            outfile=os.path.join(self.td, 'output.txt'))
        self.assertEqual(
            b"This is a test\n", data)
        with open(os.path.join(self.td, 'output.txt'), 'rb') as f:
            data = f.read()
        self.assertEqual(
            b"This is a test\n", data)
        self.assertEqual(
            0, os.path.getmtime(os.path.join(self.td, 'output.txt')))

    def testGetUrlToFileCached(self):
        """Call get_url, write to file with cache"""
        with open(os.path.join(self.td, 'foo.txt'), 'wb') as f:
            f.write(b'This is a test\n')
        outfile = os.path.join(self.td, 'output.txt')
        data = spacepy.toolbox.get_url(
            'http://localhost:{}/foo.txt'.format(self.port),
            outfile=outfile)
        data = spacepy.toolbox.get_url(
            'http://localhost:{}/foo.txt'.format(self.port),
            outfile=outfile, cached=True)
        self.assertTrue(data is None)

    def testGetUrlToFileCachedKeepalive(self):
        """Call get_url keepalive, write to file with cache"""
        with open(os.path.join(self.td, 'foo.txt'), 'wb') as f:
            f.write(b'This is a test\n')
        outfile = os.path.join(self.td, 'output.txt')
        data, conn = spacepy.toolbox.get_url(
            'http://localhost:{}/foo.txt'.format(self.port),
            outfile=outfile, keepalive=True)
        try:
            data, conn = spacepy.toolbox.get_url(
                'http://localhost:{}/foo.txt'.format(self.port),
                outfile=outfile, cached=True, conn=conn, keepalive=True)
        finally:
            conn.close()
        self.assertTrue(data is None)

    def testGetUrlToFileCachedKeepaliveChanged(self):
        """Call get_url keepalive, write to file with cache"""
        with open(os.path.join(self.td, 'foo.txt'), 'wb') as f:
            f.write(b'This is a test\n')
        outfile = os.path.join(self.td, 'output.txt')
        data, conn = spacepy.toolbox.get_url(
            'http://localhost:{}/foo.txt'.format(self.port),
            outfile=outfile, keepalive=True)
        try:
            os.utime(outfile, (0., 0.))
            data, conn = spacepy.toolbox.get_url(
                'http://localhost:{}/foo.txt'.format(self.port),
                outfile=outfile, cached=True, conn=conn, keepalive=True)
        finally:
            conn.close()
        self.assertEqual(b'This is a test\n', data)


if __name__ == "__main__":
    unittest.main()
