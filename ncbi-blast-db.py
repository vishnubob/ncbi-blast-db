#!/usr/bin/env python

import os
import re
import ftplib
import pprint
import fnmatch
import subprocess
import threading
import md5
import logging
import argparse
import ConfigParser

logger = logging.basicConfig(level=logging.DEBUG)

class ThreadQueue(object):
    def __init__(self):
        self.que = []
        self.que_cv = threading.Condition()

    def deque(self):
        self.que_cv.acquire()
        if len(self.que) == 0:
            self.que_cv.wait()
        item = self.que.pop()
        self.que_cv.release()
        return item

    def enque(self, work):
        self.que_cv.acquire()
        self.que.insert(0, work)
        self.que_cv.notify()
        self.que_cv.release()

    def __iter__(self):
        while 1:
            item = self.deque()
            if item == None:
                break
            yield item

class DownloadThread(threading.Thread):
    def __init__(self, connection):
        self.connection = connection
        self.work_que = ThreadQueue()
        self.finished_que = ThreadQueue()
        super(DownloadThread, self).__init__()
        self.daemon = True

    def run(self):
        try:
            for work in self.work_que:
                if self.download(*work):
                    self.finished_que.enque(work)
        finally:
            self.finished_que.enque(None)

    def download(self, target, filename_hash):
        md5hash = md5.new()
        temp_filename = "%s_download" % target
        filename = os.path.split(target)[-1]
        msg = "Downloading %s" % filename
        logging.info(msg)
        with open(temp_filename, 'w') as fh:
            def callback(block):
                md5hash.update(block)
                fh.write(block)
            cmd = "RETR %s" % filename
            self.connection.retrbinary(cmd, callback)
        if filename_hash and md5hash.hexdigest() != filename_hash:
            msg = "%s does not match hash" % filename
            logging.error(msg)
            os.unlink(temp_filename)
            return False
        os.rename(temp_filename, target)
        return True

class NCBI_BlastMirror(object):
    NCBI_BLASTDB_URL = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/"
    NCBI_HOST = "ftp.ncbi.nlm.nih.gov"
    NCBI_BLAST_DIR = "/blast/db"
    HASH_FILENAME = "blastdb.md5"

    def __init__(self, blastdb_dir=None, archive_dir=None, config_path=None, hashfile=None, include=None, exclude=None):
        self.blastdb_dir = blastdb_dir if blastdb_dir != None else os.getcwd()
        self.archive_dir = archive_dir if archive_dir != None else os.path.join(self.blastdb_dir, "archives")
        self.config_path = config_path if config_path != None else os.path.join(self.blastdb_dir, ".config.ini")
        self.hashfile = hashfile if hashfile != None else os.path.join(self.archive_dir, self.HASH_FILENAME)
        self.include = include if include != None else ["*"]
        self.exclude = exclude if exclude != None else []
        self.re_whitespace = re.compile("\s+")

    def initialize(self):
        if not os.path.isdir(self.blastdb_dir):
            os.makedirs(self.blastdb_dir)
        if not os.path.isdir(self.archive_dir):
            os.makedirs(self.archive_dir)
        config = ConfigParser.RawConfigParser()
        config.add_section("BLAST_MirrorConfig")
        config.set("BLAST_MirrorConfig", "blastdb_dir", self.blastdb_dir)
        config.set("BLAST_MirrorConfig", "archive_dir", self.archive_dir)
        config.set("BLAST_MirrorConfig", "hashfile", self.hashfile)
        config.set("BLAST_MirrorConfig", "include", ["*"])
        config.set("BLAST_MirrorConfig", "exclude", [])
        with open(self.config_path, 'wb') as configfile:
            config.write(configfile)
        msg = "Initialized configuration file at %s" % self.config_path
        print(msg)

    def load_local_hashmap(self):
        hashmap = {}
        if not os.path.exists(self.hashfile):
            return hashmap
        with open(self.hashfile, 'r') as fh:
            for line in fh:
                (md5hash, filename) = self.re_whitespace.split(line.strip())
                hashmap[filename] = md5hash
        return hashmap

    def write_local_hashmap(self, hashmap):
        fnlist = hashmap.keys()
        fnlist.sort()
        content = ''
        for fn in fnlist:
            md5hash = hashmap[fn]
            content += "%s  %s\n" % (fn, md5hash)
        with open(self.hashfile, 'w') as fh:
            fh.write(content)

    def load_remote_hashmap(self):
        files = self.retrlines("NLST")
        hashes = [fn for fn in files if fn.endswith(".md5")]
        dblist = set([fn.split('.')[0] for fn in hashes])
        dblist = set(self.filter_dblist(dblist))
        hashmap = {}
        for fn in hashes:
            dbname = fn.split('.')[0]
            if dbname not in dblist:
                continue
            (filename, md5hash) = self.download_hash(fn)
            hashmap[filename] = md5hash
        return hashmap

    def manifest(self):
        self.local_hashmap = self.load_local_hashmap()
        self.remote_hashmap = self.load_remote_hashmap()
        for fn in self.remote_hashmap:
            if fn in self.local_hashmap and self.local_hashmap[fn] == self.remote_hashmap[fn]:
                continue
            yield fn

    def sync(self):
        self.connect()
        download_thread = DownloadThread(self.ncbi)
        download_thread.start()
        for fn in self.manifest():
            target = os.path.join(self.archive_dir, fn)
            work = (target, self.remote_hashmap[fn])
            download_thread.work_que.enque(work)
        download_thread.work_que.enque(None)
        for (target, local_hash) in download_thread.finished_que:
            filename = os.path.split(target)[-1]
            if not self.install(target):
                msg = "%s failed to uncompress" % target
                logging.error(msg)
                continue
            self.local_hashmap[fn] = self.remote_hashmap[fn]
            self.write_local_hashmap(self.local_hashmap)

    def filter_dblist(self, dblist):
        def match(dbname, pattern_list):
            return any([fnmatch.fnmatch(dbname, pattern) for pattern in pattern_list])
        def filter_rule(dbname):
            return match(dbname, self.include) and not match(dbname, self.exclude)
        return filter(filter_rule, dblist)

    def connect(self):
        self.ncbi = ftplib.FTP(self.NCBI_HOST)
        self.ncbi.set_debuglevel(0)
        self.ncbi.login()
        self.ncbi.cwd(self.NCBI_BLAST_DIR)

    def download_hash(self, fn):
        content = []
        def callback(block):
            content.append(block)
        cmd = "RETR %s" % fn
        self.ncbi.retrbinary(cmd, callback)
        content = str.join('', content)
        content = content.strip()
        (md5hash, filename) = self.re_whitespace.split(content)
        return (filename, md5hash)
    
    def install(self, fn):
        fn = os.path.join(self.archive_dir, fn)
        cmd = ["tar", "xzf", fn]
        cmd = str.join(' ', cmd)
        msg = "Installing %s into %s" % (fn, self.blastdb_dir)
        logging.info(msg)
        proc = subprocess.Popen(cmd, cwd=self.blastdb_dir, shell=True)
        return_code = proc.wait()
        return return_code == 0
    
    def retrlines(self, cmd):
        buf = []
        def callback(line):
            buf.append(line)
        self.ncbi.retrlines(cmd, callback)
        return buf

def get_cli():
    parser = argparse.ArgumentParser(description="blastdb command")
    parser.add_argument("-c", "--config_path", help="path to configuration file")
    parser.add_argument("-d", "--blastdb_dir", help="path to BLAST database")
    parser.add_argument("-a", "--archive_dir", help="path to archive directory")
    parser.add_argument("-H", "--hashfile", help="path to local hash file")
    parser.add_argument("-i", "--init", action="store_true", help="initialize BLAST database")
    defaults = {
        "blastdb_dir": os.environ.get("BLASTDB", None),
        "archive_dir": None,
        "config_path": ".config.ini",
        "hashfile": None,
        "init": False,
    }
    parser.set_defaults(**defaults)
    args = parser.parse_args()
    return args

def load_config_path(config_path):
    config = ConfigParser.ConfigParser()
    config.read(config_path)
    conf = dict(config.items("BLAST_MirrorConfig"))
    conf["include"] = eval(conf.get("include", "[*]"))
    conf["exclude"] = eval(conf.get("exclude", "[]"))
    return conf

if __name__ == "__main__":
    args = get_cli()
    ns = args.__dict__.copy()
    init = ns.pop("init", False)
    if init:
        mirror = NCBI_BlastMirror(**ns)
        mirror.initialize()
    else:
        ns = load_config_path(args.config_path)
        mirror = NCBI_BlastMirror(**ns)
        mirror.sync()
