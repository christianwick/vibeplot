from functools import wraps


def coroutine(func):
    @wraps(func)
    def start(*args, **kwargs):
        cr = func(*args, **kwargs)
        cr.next()
        return cr
    return start

@coroutine
def broadcast(targets):
    while True:
        line = (yield)
        for target in targets:
            target.send(line)

@coroutine
def parse_section(section_name, target):
    in_section = False
    while True:
        line = (yield)
        if line.startswith(section_name):
            in_section = True
            continue
        elif line.startswith("["):
            in_section = False
        if in_section:
            target.send(line)

@coroutine
def to_type(type_, target):
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
def printer():
    while True:
        obj = (yield)
        print obj


