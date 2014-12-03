# Copyright (c) 2011-2014 Mathias Laurin, 3-clause BSD License
"""Coroutine library."""
from __future__ import print_function
import sys
from functools import wraps


def coroutine(func):
    @wraps(func)
    def start(*args, **kwargs):
        cr = func(*args, **kwargs)
        next(cr)
        return cr
    return start


@coroutine
def broadcast(targets):
    while True:
        line = (yield)
        for target in targets:
            target.send(line)


@coroutine
def parse_section(section, intarget, offtarget):
    target = offtarget
    while True:
        line = (yield)
        if line.startswith(section):
            offtarget.send(line)  # section header always off target
            target = intarget
            continue
        elif line.startswith("[") or line.startswith(section):
            target = offtarget
        target.send(line)


@coroutine
def periodic_split(range, period, intarget, offtarget):
    counter = 0
    while True:
        counter %= period
        obj = (yield)
        (intarget if counter in range else offtarget).send(obj)
        counter += 1


@coroutine
def convert(type_, target):
    while True:
        obj = (yield)
        target.send(type_(obj))


@coroutine
def append(lst):
    while True:
        obj = (yield)
        lst.append(obj)


@coroutine
def extend(dct):
    while True:
        key, value = (yield)
        dct.setdefault(key, []).append(value)


@coroutine
def dct_set(dct):
    while True:
        key, value = (yield)
        dct[key] = value


@coroutine
def printer(file=sys.stdout):
    printer.close = file.close
    while True:
        obj = (yield)
        print(obj, file=file, end="")
