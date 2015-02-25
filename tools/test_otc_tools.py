#!/usr/bin/env python
import subprocess

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
    if os.path.exists(exp_outf):
        expected_output = codecs.open(exp_outf, 'r', encoding='utf-8').read()
        obtained_output = codecs.open(obt_outf, 'r', encoding='utf-8').read()
        if obtained_output != expected_output:
            if expected_output.strip() != '' or obtained_output.strip() != '':
                FAILED_TESTS.append(tag)
                error('OUTPUT differed for {}:\n'.format(tag))
                subprocess.call(['diff', exp_outf, obt_outf])
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
        tag = descrip['expected']
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
            tag = os.path.join(descrip['expected'], inp_file)
            test_invoc(tag, sub_invoc, expected_for_full_invoc)
    else:
        assert False # no inputs needed
if __name__ == '__main__':
    import codecs
    import json
    import sys
    import os
    dat_dir, expected_dir, tools_build_dir = [os.path.abspath(i) for i in sys.argv[1:]]
    e_dir_list = os.listdir(expected_dir)
    e_dir_list.sort()
    for e_subdir_name in e_dir_list:
        e_path = os.path.join(expected_dir, e_subdir_name)
        t_path = os.path.join(tools_build_dir, e_subdir_name)
        if os.path.isdir(t_path):
            th = [i for i in os.listdir(e_path) if i.endswith('.json')]
            th.sort()
            for jfn in th:
                test_json_path = os.path.join(e_path, jfn)
                with codecs.open(test_json_path, 'r', encoding='utf-8') as jinf:
                    test_list = json.load(jinf)
                assert isinstance(test_list, list)
                for test_descrip in test_list:
                    run_tests_for_invocation(test_descrip, dat_dir, e_path, t_path)
    passed = NUM_TESTS - len(FAILED_TESTS)
    sys.stderr.write('Passed {p:d}/{t:d} tests.'.format(p=passed, t=NUM_TESTS))
    if FAILED_TESTS:
        sys.exit(' Failed:\n    {}\n'.format('\n    '.join(FAILED_TESTS)))
    sys.stderr.write('SUCCESS\n')

