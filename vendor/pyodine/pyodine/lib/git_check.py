# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 10:42:24 2021

Some functions to return information about the current git branch and
commit hash.

@author: pheeren
"""


import subprocess, os
import logging


def get_git_revision_hash(pathname=None) -> str:
    if not isinstance(pathname, str):
        dirname = os.path.dirname(os.path.abspath(__file__))
    else:
        dirname = pathname
    try:
        return subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=dirname).decode('ascii').strip()
    except Exception as e:
        logging.error('git revision hash failed:', exc_info=True)
        return 'git revision hash failed:\n' + str(e.args)

def get_git_revision_short_hash(pathname=None) -> str:
    if not isinstance(pathname, str):
        dirname = os.path.dirname(os.path.abspath(__file__))
    else:
        dirname = pathname
    try:
        return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=dirname).decode('ascii').strip()
    except Exception as e:
        logging.error('git revision short hash failed:', exc_info=True)
        return 'git revision short hash failed:\n' + str(e.args)
        

def get_git_branch_name(pathname=None) -> str:
    if not isinstance(pathname, str):
        dirname = os.path.dirname(os.path.abspath(__file__))
    else:
        dirname = pathname
    try:
        return subprocess.check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD'], cwd=dirname, universal_newlines=True).strip()
    except Exception as e:
        logging.error('git branch name failed:', exc_info=True)
        return 'git branch name failed:\n' + str(e.args[0])