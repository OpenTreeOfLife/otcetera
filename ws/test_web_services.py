#!/usr/bin/env python
import subprocess
import requests
import json
import time
import logging
try:
    from Queue import Queue
except:
    from queue import Queue
from threading import Thread, RLock
_LOG = logging.getLogger(__name__)
_LOG.setLevel(logging.DEBUG)
_lh = logging.StreamHandler()
_lh.setFormatter(logging.Formatter("[%(asctime)s] %(filename)s (%(lineno)3d): %(levelname) 8s: %(message)s"))
_LOG.addHandler(_lh)

NUM_TESTS = 0
FAILED_TESTS = []
FRAC_FLOAT_DIFF_TOL = 0.001

def _extend_diff_list(diff_list, r):
    if r:
        if isinstance(r, list):
            diff_list.extend(r)
        else:
            diff_list.append(r)

def gen_dict_diff_str(expected, observed, ex_pref, obs_pref):
    if expected == observed:
        return None
    diff_list = []
    if isinstance(expected, dict):
        if not isinstance(observed, dict):
            return '{} is a dict, but {} is a {}'.format(ex_pref, obs_pref, type(observed))
        for ek, ev in expected.items():
            if ek in observed:
                ov = observed[ek]
                if ov != ev:
                    matched_as_floats = False
                    if isinstance(ev, float):
                        try:
                            if abs(ev - ov)/abs(ev) <= FRAC_FLOAT_DIFF_TOL:
                                matched_as_floats = True
                        except:
                            pass
                    if not matched_as_floats:
                        r = gen_dict_diff_str(ev, ov, '{}["{}"]'.format(ex_pref, ek), '{}["{}"]'.format(obs_pref, ek))
                        _extend_diff_list(diff_list, r)

            else:
                diff_list.append('{}["{}"] is absent'.format(obs_pref, ek))
        for k in observed.keys():
            if k not in expected:
                diff_list.append('{}["{}"] was present, but not an expected key'.format(obs_pref, k))
    elif isinstance(expected, list) or isinstance(expected, tuple):
        if not isinstance(observed, list) or isinstance(observed, tuple):
            return '{} is a list, but {} is a {}'.format(ex_pref, obs_pref, type(observed))
        if len(expected) != len(observed):
            diff_list.append('{} had {} elments but {} has {}'.format(ex_pref, len(expected), obs_pref, len(observed)))
        else:
            ml = len(expected)
            for ind in range(ml):
                eel, oel = expected[ind], observed[ind]
                if eel != oel:
                    r = gen_dict_diff_str(eel, oel, '{}[{}]'.format(ex_pref, ind), '{}[{}]'.format(obs_pref, ind))
                    _extend_diff_list(diff_list, r)
    elif type(expected) == type(observed):
        return ['{} = {}, but {} = {}'.format(ex_pref, repr(expected), obs_pref, repr(observed))]
    else:
        return ['{} is the {} equal to {}, but {} is a {}'.format(ex_pref, type(expected), repr(expected), obs_pref, type(observed))]
    return diff_list

def gen_expected_obs_diff(expected, observed, tag):
    return gen_dict_diff_str(expected, observed, 'Expected {}'.format(tag), 'Observed {}'.format(tag))

#########################################################################################
# The following code for execution in a non-blocking thread is from pyraphyletic. If
#    we moved it to peyotl, we could import it from there (at the cost of making)
#    otcetera depend on peyot.
class JobQueue(Queue):
    """Thread-safe Queue that logs the addition of a job to debug"""
    def put(self, item, block=None, timeout=None):
        """Logs `item` at the debug level then calls base-class put"""
        _LOG.debug("%s queued" % str(item))
        Queue.put(self, item, block=block, timeout=timeout)


_jobq = JobQueue()


def worker():
    """Infinite loop of getting jobs off of _jobq and performing them."""
    while True:
        job = _jobq.get()
        _LOG.debug('"{}" started"'.format(job))
        try:
            job.start()
        except:
            _LOG.exception("Worker dying.")
        else:
            try:
                job.get_results()
            except:
                _LOG.exception("Worker exception.  Error in job.get_results")
        _LOG.debug('"{}" completed'.format(job))
        _jobq.task_done()


_WORKER_THREADS = []


def start_worker(num_workers):
    """Spawns worker threads such that at least `num_workers` threads will be
    launched for processing jobs in the jobq.

    The only way that you can get more than `num_workers` threads is if you
    have previously called the function with a number > `num_workers`.
    (worker threads are never killed).
    """
    assert num_workers > 0, "A positive number must be passed as the number of worker threads"
    num_currently_running = len(_WORKER_THREADS)
    for i in range(num_currently_running, num_workers):
        _LOG.debug("Launching Worker thread #%d" % i)
        t = Thread(target=worker)
        _WORKER_THREADS.append(t)
        t.setDaemon(True)
        t.start()

#########################################################################################
_verb_name_to_req_method = {"GET": requests.get,
                            "PUT": requests.put,
                            "POST": requests.post,
                            "DELETE": requests.delete,
                            "HEAD": requests.head,
                            "OPTIONS": requests.options,
                            }
API_HEADERS = {'content-type' : 'application/json',
               'accept' : 'application/json',
              }
class WebServiceTestJob(object):
    def __init__(self, test_par, test_description, service_prefix):
        self.url_fragment = test_description["url_fragment"]
        self.arguments = test_description["arguments"]
        v = test_description.get("verb", "GET").upper()
        self.requests_method = _verb_name_to_req_method[v]
        self.service_prefix = service_prefix
        self.url = service_prefix + self.url_fragment
        self.expected = test_description.get('expected_response_payload')
        self.expected_status = test_description.get('expected_status_code')
        self._status_str = None
        self.passed = False
        self.failed = False
        self.erred = False
        self.test_par = test_par
        self.test_dir = test_description.get("test_dir")
        self.test_subdir = os.path.relpath(self.test_dir, self.test_par)
        self.name = test_description.get("name", self.test_subdir or self.url_fragment)
        self.stat_lock = RLock()
    
    @property
    def status_str(self):
        with self.stat_lock:
            if self._status_str is None:
                x = None
            else:
                x = str(self._status_str)
            
        return x
    @status_str.setter
    def status_str(self, value):
        with self.stat_lock:
            self._status_str = value
        

    def __str__(self):
        return 'WebServiceTestJob {}'.format(self.name)

    def run_ws_test(self):
        self.status_str = ''
        try:
            # 1. Make the call
            if self.arguments:
                _LOG.debug("{} arguments = {}".format(self.name, repr(self.arguments)))
                response = self.requests_method(self.url, headers=API_HEADERS, data=json.dumps(self.arguments))
            else:
                response = self.requests_method(self.url)

            # 2.A Raise exception if we expected status 200 and didn't get it.
            if self.expected_status == 200:
                try:
                    response.raise_for_status()
                except Exception as sce:
                    _LOG.exception('exception url: {}'.format(self.url))
                    try:
                        self.status_str = "Non-200 response body = {}\n".format(response.text)
                    except:
                        pass
                    raise sce
            # 2.B. Return if we expected an error code and didn't get it
            elif response.status_code != self.expected_status:
                self.failed = True
                try:
                    self.status_str = "Expected status {} but got {}.  response body = {}\n".format(self.expected_status, response.status_code, response.text)
                except:
                    pass
                return

            # 3. Check JSON body
            _LOG.debug('name: {}  Expected: {}'.format(self.name, self.expected))
            if self.expected is not None:
                try:
                    j = response.json()
                except:
                    _LOG.error("{} no JSON in response: {}".format(self.name, response.text))
                    raise
                _LOG.debug('name: {} Observed: {}'.format(self.name, j))

                if j != self.expected:
                    dd = gen_expected_obs_diff(self.expected, j, 'x')
                    if dd:
                        self.failed = True
                        if self.test_dir:
                            dbout_observed = os.path.join(self.test_dir, "observed.json")
                            with codecs.open(dbout_observed, 'w', encoding="utf-8") as obsfo:
                                json.dump(j, obsfo, sort_keys=True, indent=2, separators=(',', ': '))
                            m = 'Response written to {}'.format(dbout_observed)
                        else:
                            m = ''
                        self.status_str = "Wrong response:\n{}\n{}".format('\n'.join(dd), m)
                        return
            self.passed = True
            self.status_str = "Completed"
        except Exception as x:
            self.erred = True
            _LOG.exception('writing exception to status string')
            self.status_str += "Exception: {}".format(x)

    def start(self):
        """Trigger to start push - blocking"""
        self.run_ws_test()

    def get_results(self):
        """:return self.status_str"""
        return self.status_str
#########################################################################################

PIDFILE_NAME = "pidfile.txt"
RUNNING_SERVER = None
SERVER_PORT = 1985 # global, set by CLI. Needed by server launch and threads
SERVER_OUT_ERR_FN = "test-server-stdouterr.txt"

def launch_server(exe_dir, taxonomy_dir, synth_par, server_threads=4):
    global RUNNING_SERVER
    exe_path = os.path.join(exe_dir, 'otc-tol-ws')
    pidfile_path = os.path.join(exe_dir, PIDFILE_NAME)
    server_std_out = os.path.join(exe_dir, SERVER_OUT_ERR_FN)
    invocation = [exe_path,
                  taxonomy_dir,
                  "-D" + synth_par,
                  "-p{}".format(pidfile_path),
                  "-P{}".format(SERVER_PORT), 
                  "--num-threads={}".format(server_threads)]
    _LOG.debug('Launching with: "{}"'.format('" "'.join(invocation)))
    with open(server_std_out, 'w') as sstdoe:
        RUNNING_SERVER = subprocess.Popen(invocation,
                                          stdout=sstdoe,
                                          stderr=subprocess.STDOUT)
        wc = 0
        while (RUNNING_SERVER.poll() is None) and (not os.path.exists(pidfile_path)):
            time.sleep(0.1)
            if wc > 100:
                raise RuntimeError("Assuming that the server has hung after waiting for pidfile")
            wc += 1
    return (RUNNING_SERVER.poll() is None) and (os.path.exists(pidfile_path))

def kill_server(exe_dir):
    if os.environ.get("PROMPT_BEFORE_KILLING_SERVER"):
        c = raw_input("type any key to kill the server... ")
    pidfile_path = os.path.join(exe_dir, PIDFILE_NAME)
    if RUNNING_SERVER.poll() is None:
        RUNNING_SERVER.terminate()
        if RUNNING_SERVER.poll() is None:
            RUNNING_SERVER.kill()
        wc = 0
        while RUNNING_SERVER.poll() is None:
            time.sleep(0.1)
            wc += 1
            if wc > 100:
                break
    if RUNNING_SERVER.poll() is None:
        sys.stderr.write("Could not kill server! Kill it then remove the pidfile.txt\n")
    else:
        if os.path.exists(pidfile_path):
            os.remove(pidfile_path)
        sys.stderr.write("Server no longer running and pidfile removed.\n")

FAILED_TESTS, ERRORED_TESTS = [], []

def run_tests(test_par, dirs_to_run, test_threads):
    assert test_threads > 0
    td_list = []
    for test_dir in dirs_to_run:
        with codecs.open(os.path.join(test_dir, "method.json")) as inp:
            td = json.load(inp)
        if os.path.exists(os.path.join(test_dir, "expected.json")):
            with codecs.open(os.path.join(test_dir, "expected.json")) as inp:
                td["expected_response_payload"] = json.load(inp)

        if os.path.exists(os.path.join(test_dir, "expected_code.txt")):
            with codecs.open(os.path.join(test_dir, "expected_code.txt")) as inp:
                lines = inp.readlines()
            td["expected_status_code"] = int(lines[0])
        else:
            td["expected_status_code"] = 200

        td["test_dir"] = test_dir
        td_list.append(td)

    start_worker(test_threads)
    service_prefix = "http://127.0.0.1:{}/".format(SERVER_PORT)
    all_jobs = [WebServiceTestJob(test_par=test_par, test_description=td, service_prefix=service_prefix) for td in td_list]
    running_jobs = list(all_jobs)
    for j in all_jobs:
        _jobq.put(j)

    # now we block until all jobs have a status_str
    num_passed = 0 
    num_failed = 0
    num_errors = 0
    while True:
        srj = []
        for j in running_jobs:
            jss = j.status_str
            if not jss:
                # _LOG.debug('putting {} back in queue'.format(j.name))
                srj.append(j)
                continue
            _LOG.debug('test {} status_str = {} resolved'.format(j.name, repr(jss)))
            if j.erred or j.failed:
                if j.failed:
                    fc = "FAILURE"
                    num_failed += 1
                    FAILED_TESTS.append(j.name)
                else:
                    fc = "ERROR"
                    num_errors += 1
                    ERRORED_TESTS.append(j.name)
                _LOG.error("{} {}: {}".format(j.name, fc, j.status_str))
            else:
                num_passed += 1
                _LOG.debug("{} passed.".format(j.name))
        if not srj:
            break
        running_jobs = srj
        time.sleep(0.1)
    return num_passed, num_failed, num_errors


def get_test_dirs_under(top_test_dir):
    test_dirs = []
    for root, dirs, files in os.walk(top_test_dir):
        if "method.json" in files:
            path = os.path.relpath(root, top_test_dir)
            test_dirs.insert(0,path)
    return test_dirs


if __name__ == '__main__':
    import argparse
    import codecs
    import sys
    import os
    parser = argparse.ArgumentParser(description="Runs the otc-tol-ws and tests described in method.json files")
    parser.add_argument('--taxonomy-dir', required=True, help='Directory that is the parent of the taxonomy files')
    parser.add_argument('--synthesis-parent', required=True, help='Directory that is the parent of synthesis directories (if there is more than one subdirectory, then there will be multiple trees served - that option is not well tested).')
    parser.add_argument('--exe-dir', required=True, help='Directory that holds the otc-tol-ws executable and which will be the working directory of the server.')
    parser.add_argument('--tests-parent', required=True, help='Directory. Each subdir that holds a "method.json" file will be interpreted as a test.')
    parser.add_argument('--test-name', default=None, required=False, help='Name of a subdir of the tests-parent dir. If provided only that test will be run; otherwise all of the tests will be run.')
    parser.add_argument('--server-port', default=1985, type=int, required=False, help='Port number for the server')
    parser.add_argument('--server-threads', default=4, type=int, required=False, help='Number of threads for the server')
    parser.add_argument('--test-threads', default=8, type=int, required=False, help='Number of threads launched for running tests.')
    parser.add_argument('--secs-to-recheck-pid-file', default=0, type=int, required=False, help='If the pid file exists, the process will enter a loop sleeping and rechecking for this number of seconds.')
    
    args = parser.parse_args()
    if args.server_threads < 1 or args.test_threads < 1:
        sys.exit("The number of threads must be positive.")
    taxonomy_dir = args.taxonomy_dir
    if not os.path.isdir(taxonomy_dir):
        sys.exit('Taxonomy directory "{}" does not exist.\n'.format(taxonomy_dir))
    synth_par_path = args.synthesis_parent
    if not os.path.isdir(synth_par_path):
        sys.exit('Synthetic tree parent directory "{}" does not exist.\n'.format(synth_par_path))
    exe_dir = args.exe_dir
    if not os.path.isdir(exe_dir):
        sys.exit('Executable directory "{}" does not exist.\n'.format(exe_dir))
    test_par = args.tests_parent
    if not os.path.isdir(test_par):
        sys.exit('Tests parent directory "{}" does not exist.\n'.format(test_par))
    if args.test_name is not None:
        e_dir_list = [args.test_name]
    else:
        e_dir_list = get_test_dirs_under(test_par)
        e_dir_list.sort()
    SERVER_PORT = args.server_port

    # Get test paths
    to_run = []
    for e_subdir_name in e_dir_list:
        e_path = os.path.join(test_par, e_subdir_name)
        if not os.path.isdir(e_path):
            sys.stderr.write("Skipping test {} due to missing dir {} \n".format(e_subdir_name, e_path))
            continue
        mfile = os.path.join(e_path, "method.json")
        if not os.path.isfile(mfile):
            sys.stderr.write("Skipping test {} due to missing file {}\n".format(e_subdir_name, mfile))
            continue
        to_run.append(e_path)
    if not to_run:
        sys.exit("No test were found!")

    # Check that there are no PIDfiles in the way
    pidfile_path = os.path.join(exe_dir, PIDFILE_NAME)
    if os.path.exists(pidfile_path):
        recheck = 0
        checks_per_sec = 3
        while recheck < checks_per_sec*args.secs_to_recheck_pid_file:
            time.sleep(1.0/checks_per_sec)
            if not os.path.exists(pidfile_path):
                break
        if os.path.exists(pidfile_path):
            sys.exit("{} is in the way!\n".format(pidfile_path))

    # try launching otc-tol-ws and running the tests against it.
    for i in range(2):
        if launch_server(exe_dir=exe_dir,
                        taxonomy_dir=taxonomy_dir,
                        synth_par=synth_par_path,
                        server_threads=args.server_threads):
            try:
                num_passed, nf, ne = run_tests(test_par, to_run, args.test_threads)
            finally:
                kill_server(exe_dir)
            NUM_TESTS = nf + ne + num_passed
            assert nf == len(FAILED_TESTS)
            assert ne == len(ERRORED_TESTS)
            sys.stderr.write('Passed {p:d}/{t:d} tests.'.format(p=num_passed, t=NUM_TESTS))
            if FAILED_TESTS:
                sys.stderr.write(' Failed:\n    {}\n'.format('\n    '.join(FAILED_TESTS)))
            if ERRORED_TESTS:
                sys.stderr.write(' Errors in:\n    {}\n'.format('\n    '.join(ERRORED_TESTS)))
            if nf + ne > 0:
                sys.exit(nf + ne)
            sys.stderr.write('SUCCESS\n')
            sys.exit(0)
        else:
            time.sleep(1) # relaunch (most likely cause is the port not being freed from previous test) 
    _LOG.error("Server launch failed: ")
    with open(os.path.join(exe_dir, SERVER_OUT_ERR_FN), 'r') as seo:
        sys.stderr.write(seo.read())
    sys.exit(-1)

