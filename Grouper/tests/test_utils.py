from Grouper import GroupGraph, AtomGraph
from Grouper.io import has_mbuild, has_torch
from Grouper.tests.base_test import BaseTest
import networkx as nx
import pytest
import numpy as np

class TestGroupGraph(BaseTest):
    def convert_to_nx(self):
       pass