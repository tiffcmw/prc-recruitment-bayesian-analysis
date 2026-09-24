import numpy as np
from sortedcontainers import SortedSet

class Lattice:
  def __init__(self, microtubule_length, microtubule_width, microtubule_spacing, microtubule_separation, microtubule_offset):
    self.microtubule_length = microtubule_length
    self.microtubule_width = microtubule_width
    self.microtubule_spacing = microtubule_spacing
    self.microtubule_separation = microtubule_separation
    self.microtubule_offset = microtubule_offset

    self.num_rows = int(microtubule_width / microtubule_spacing)
    self.num_cols = int(microtubule_length / microtubule_spacing) + 1
    self.num_sites = self.num_rows * self.num_cols

    self.top_taken = SortedSet()
    self.bottom_taken = SortedSet()

    self.top_single_adj = SortedSet()
    self.top_double_adj = SortedSet()
    self.top_triple_adj = SortedSet()
    self.top_quad_adj = SortedSet()
    self.bottom_single_adj = SortedSet()
    self.bottom_double_adj = SortedSet()
    self.bottom_triple_adj = SortedSet()
    self.bottom_quad_adj = SortedSet()

  def to_index(self, row, col):
    return row * self.num_cols + col

  def to_coord(self, index):
    row = index // self.num_cols
    col = index % self.num_cols
    return row, col

  def get_neighbors(self, index):
    row, col = self.to_coord(index)
    nbrs = []
    for dr, dc in [(-1, 0), (+1, 0), (0, -1), (0, +1)]:
      nbr_row, nbr_col = row + dr, col + dc
      if 0 <= nbr_row < self.num_rows and 0 <= nbr_col < self.num_cols:
        nbrs.append(self.to_index(nbr_row, nbr_col))
    return nbrs

  def get_distance(self, bottom_index, top_index):
    row_bottom, col_bottom = self.to_coord(bottom_index)
    row_top, col_top = self.to_coord(top_index)
    col_dist = (col_top + self.offset - col_bottom) * self.spacing
    row_dist = (row_top - row_bottom) * self.spacing
    z_dist = self.separation
    return np.sqrt(col_dist**2 + row_dist**2 + z_dist**2)