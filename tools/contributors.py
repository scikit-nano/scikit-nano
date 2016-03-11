#!/usr/bin/env python
from __future__ import unicode_literals
import subprocess
import sys
import string
import shlex

if len(sys.argv) != 2:
    print("Usage: ./contributors.py tag-of-previous-release")
    sys.exit(-1)

tag = sys.argv[1]


if not sys.version_info[:2] == (2, 6):

    def call(cmd):
        return subprocess.check_output(shlex.split(cmd),
                                       universal_newlines=True).split('\n')

    tag_date = call("git log --pretty='format:%ci' {}".format(tag))[0]
    # print('tag_date: {}'.format(tag_date))

    print("Release {} was on {}\n".format(tag, tag_date))

    merges = \
        call("git log --since='{}' --merges --format='>>>%%B' --reverse".format(tag_date))
    merges = [m for m in merges if m.strip()]
    merges = '\n'.join(merges).split('>>>')
    merges = [m.split('\n')[:2] for m in merges]
    merges = [m for m in merges if len(m) == 2 and m[1].strip()]

    num_commits = call("git rev-list {}..HEAD --count".format(tag))[0]
    print("A total of {} changes have been committed.\n".format(num_commits))

    print("It contained the following {:d} merges:\n".format(len(merges)))
    for (merge, message) in merges:
        if merge.startswith('Merge pull request #'):
            PR = ' ({})'.format(merge.split())[3]
        else:
            PR = ''

        print('- ' + message + PR)

    print("\nMade by the following committers [alphabetical by last name]:\n")

    authors = call("git log --since='{}' --format=%aN".format(tag_date))
    authors = [a.strip() for a in authors if a.strip()]

    def key(author):
        author = [v for v in author.split() if v[0] in string.ascii_letters]
        if len(author) > 0:
            return author[-1]

    authors = sorted(set(authors), key=key)

    for a in authors:
        print('- ' + a)
