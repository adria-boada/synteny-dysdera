#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# missatges.py
#
# 10 d’abr. 2024  <adria@molevol-OptiPlex-9020>

# Measure the time it takes to run this script
from datetime import datetime
import locale
locale.setlocale(locale.LC_TIME, '') # sets locale to the one used by user

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Missatge:
    """
    Tipificació del missatge. Més endavant embolcallarem la classe pare amb
    diferents colors.
    """
    def __init__(self, missatge):
        self.data = datetime.now().strftime("%d %b, %H:%M")
        self.colors = {
            "predet": "\033[0;0m",
            "vermell": "\033[0;31m",
            "verd": "\033[0;32m",
            "groc": "\033[0;33m" }

class Estat(Missatge):
    def __init__(self, missatge):
        super().__init__(missatge)
        missatge = "STATUS:" + str(self.data) + ": " + str(missatge)
        # Canvia el color d'impressió a verd, a través dels colors ANSI.
        missatge = self.colors["verd"] + missatge + self.colors["predet"]

        print(missatge)

class Avis(Missatge):
    def __init__(self, missatge):
        super().__init__(missatge)
        missatge = "WARNING:" + str(self.data) + ": " + str(missatge)
        # Canvia el color d'impressió a groc, a través dels colors ANSI.
        missatge = self.colors["groc"] + missatge + self.colors["predet"]

        print(missatge)

class Error(Missatge):
    def __init__(self, missatge):
        super().__init__(missatge)
        missatge = "ERROR:" + str(self.data) + ": " + str(missatge)
        # Canvia el color d'impressió a vermell, a través dels colors ANSI.
        missatge = self.colors["vermell"] + missatge + self.colors["predet"]

        print(missatge)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

