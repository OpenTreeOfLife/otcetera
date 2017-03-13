#!/usr/bin/env python
import subprocess
import time

NUM_TESTS = 0
FAILED_TESTS = []

def debug(m):
    sys.stderr.write('DEBUG: ')
    sys.stderr.write(m)
    sys.stderr.write('\n')
def error(m):
    sys.stderr.write('ERROR: ')
    sys.stderr.write(m)
    sys.stderr.write('\n')

def test_invoc(tag, invocation, result_dir):
    global NUM_TESTS, FAILED_TESTS
    obt_outf = os.path.join(result_dir, 'obtained-output')
    obt_errf = os.path.join(result_dir, 'obtained-error')
    obt_exitf = os.path.join(result_dir, 'obtained-exit')
    obt_extraf = os.path.join(result_dir, 'obtained-extra')
    NUM_TESTS += 1
    with codecs.open(obt_outf, 'w', encoding='utf-8') as ob_o:
        with codecs.open(obt_errf, 'w', encoding='utf-8') as ob_e:
            invoc = '"{}"'.format('" "'.join(invocation))
            debug('Running:\n    ' + invoc + ' >"' + obt_outf + '" 2>"' + obt_errf + '" ; echo $? >"' + obt_exitf + '"')
            p = subprocess.Popen(invocation, cwd=result_dir, stdout=ob_o, stderr=ob_e)
            exit_code = p.wait()
            with codecs.open(obt_exitf, 'w', encoding='utf-8') as ob_ex:
                ob_ex.write('{e:d}\n'.format(e=exit_code))
    exp_outf = os.path.join(result_dir, 'output')
    exp_extraf = os.path.join(result_dir, 'extra')
    if os.path.exists(exp_outf):
        expected_output = codecs.open(exp_outf, 'r', encoding='utf-8').read()
        obtained_output = codecs.open(obt_outf, 'r', encoding='utf-8').read()
        if obtained_output != expected_output:
            FAILED_TESTS.append(tag)
            error('OUTPUT differed for {}:\n'.format(tag))
            subprocess.call(['diff', exp_outf, obt_outf])
    if os.path.exists(exp_extraf):
        expected_output = codecs.open(exp_extraf, 'r', encoding='utf-8').read()
        obtained_output = codecs.open(obt_extraf, 'r', encoding='utf-8').read()
        if obtained_output != expected_output:
            FAILED_TESTS.append(tag)
            error('EXTRA file output differed for {}:\n'.format(tag))
            subprocess.call(['diff', exp_extraf, obt_extraf])
            
    exp_exitf = os.path.join(result_dir, 'exit')
    if os.path.exists(exp_exitf):
        expected_exit = int(codecs.open(exp_exitf, 'r', encoding='utf-8').read())
    else:
        expected_exit = 0
    if expected_exit != exit_code:
        m = 'EXIT CODE differed for {t}: Expected {e:d} Obtained {o:d}\n'
        error(m.format(t=tag, e=expected_exit, o=exit_code))
        if tag not in FAILED_TESTS:
            FAILED_TESTS.append(tag)



def run_tests_for_invocation(descrip, data_dir, expected_par, executable_par):
    invoc = descrip['invocation']
    # make abs path to exe
    invoc[0] = os.path.join(executable_par, invoc[0])
    test_tag = os.path.split(expected_par)[1]
    expected_results_path = os.path.join(expected_par, descrip['expected'])
    if 'infile_list' in descrip:
        # tests that need multiple input files
        assert '<INFILELIST>' in invoc
        infile_list = [os.path.join(data_dir, i) for i in descrip['infile_list']]
        ind = invoc.index('<INFILELIST>')
        p = invoc[:ind]
        if len(invoc) > ind + 1:
            s = invoc[(1 + ind):]
        else:
            s = []
        invoc = p + infile_list + s
        tag = test_tag + '/' + descrip['expected']
        test_invoc(tag, invoc, expected_results_path)
    elif '<INFILE>' in invoc:
        ind = invoc.index('<INFILE>')
        inp_file_names = os.listdir(expected_results_path)
        inp_file_names.sort()
        for inp_file in inp_file_names:
            expected_for_full_invoc = os.path.join(expected_results_path, inp_file)
            inp_path = os.path.join(data_dir, inp_file)
            sub_invoc = list(invoc)
            sub_invoc[ind] = inp_path
            tag = test_tag + '/' + os.path.join(descrip['expected'], inp_file)
            test_invoc(tag, sub_invoc, expected_for_full_invoc)
    else:
        assert False # no inputs needed

PIDFILE_NAME = "pidfile.txt"
RUNNING_SERVER = None

def launch_server(ws_build_dir, data_dir):
    global RUNNING_SERVER
    exe_path = os.path.join(ws_build_dir, 'otc-tol-ws')
    pidfile_path = os.path.join(ws_build_dir, PIDFILE_NAME)
    taxonomy_path = os.path.join(data_dir, "ex-tax-1")
    summary_tree_path = os.path.join(data_dir, "ex-synth-par")
    server_std_out = os.path.join(ws_build_dir, "test-server-stdouterr.txt")
    invocation = [exe_path, taxonomy_path, "-D" + summary_tree_path, "-p" + pidfile_path]
    print(invocation)
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
    return True

def kill_server(ws_build_dir):
    pidfile_path = os.path.join(ws_build_dir, PIDFILE_NAME)
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
        os.remove(pidfile_path)
        sys.stderr.write("Server no longer running and pidfile removed.\n")

def run_tests(dirs_to_run):
    pass #c = raw_input("enter something to kill this  ")

if __name__ == '__main__':
    import codecs
    import json
    import sys
    import os
    SCRIPT_DIR = os.path.split(sys.argv[0])[0]
    data_dir, expected_dir, ws_build_dir = [os.path.abspath(i) for i in sys.argv[1:4]]
    if len(sys.argv) > 4:
        test_sub_dir_to_run = sys.argv[4]
        e_dir_list = [test_sub_dir_to_run]
    else:
        e_dir_list = os.listdir(expected_dir)
        e_dir_list.sort()

    to_run = []
    for e_subdir_name in e_dir_list:
        e_path = os.path.join(expected_dir, e_subdir_name)
        if not os.path.isdir(e_path):
            sys.stderr.write("Skipping test {} due to missing dir\n".format(e_path))
            continue
        to_run.append(e_path)
    
    pidfile_path = os.path.join(ws_build_dir, PIDFILE_NAME)
    if os.path.exists(pidfile_path):
        sys.exit("{} is in the way!\n".format(pidfile_path))
    if launch_server(ws_build_dir, data_dir):
        try:
            run_tests(to_run)
        finally:
            kill_server(ws_build_dir)
        passed = NUM_TESTS - len(FAILED_TESTS)
        sys.stderr.write('Passed {p:d}/{t:d} tests.'.format(p=passed, t=NUM_TESTS))
        if FAILED_TESTS:
            sys.exit(' Failed:\n    {}\n'.format('\n    '.join(FAILED_TESTS)))
        sys.stderr.write('SUCCESS\n')

