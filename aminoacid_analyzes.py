#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2015.04.17.
"""


import wx
from front_end import FrontEnd


# Run the program
if __name__ == "__main__":
    app = wx.App(False)
    frame = FrontEnd()
    frame.Show()
    app.MainLoop()
