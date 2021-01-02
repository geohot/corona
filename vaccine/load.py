#!/usr/bin/env python3

import pathlib
import os

# code from world health organization
# https://mednet-communities.net/inn/db/media/docs/11889.doc

with open(os.path.join(pathlib.Path(__file__).parent.absolute(), "code")) as f:
  dat = f.read().strip().replace("\n", "").replace(" ", "")

if __name__ == "__main__":
  print(len(dat))
  print(dat)



