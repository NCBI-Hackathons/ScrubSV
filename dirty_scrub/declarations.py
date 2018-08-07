#!/usr/bin/env python
from __future__ import print_function, unicode_literals


class YomoDict(dict):
    """You Only Map Once.

    Dictionary that only allows keys to be added once

    TODO: There's probably a built-in for this, use that"""

    def __init__(self, *args, **kwargs):
        super(YomoDict, self).__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        if key in self:
            raise KeyError
        super(YomoDict, self).__setitem__(key, value)
