"""Misc utilities."""

import collections
import csv

import numpy as np

def read_csv(filename, header=True):

  with open(filename, 'r') as f:
    reader = csv.reader(f)

    if header is True:
      row = reader.next()
      record = collections.namedtuple('record', row)
    elif header:
      record = collections.namedtuple('record', header)
    else:  
      record = collections.namedtuple('record', map(str, range(len(row))))

    rows = []
    for row in reader:
      rows.append(record._make(row))
      
  return rows
