import numpy as np

def get_public_attr(obj):
    # obtain the attr to be explored
    attr = dir(obj)

    # filtrate the ones that are public
    attr = [atr for atr in attr if atr[:2] != '__' and atr[-2:] != '__']

    return attr


class ColumnFunctionsMixin:
    def add_angle(self, axes = ['px', 'py', 'pz'], angle_axis = 'pz'):
      """
      Calculate the angle of the particles starting from a list of the form [px, py, pz].
      This with respect to the "axis" element.
      """
      momentumTotal = np.linalg.norm(np.array([self.data[axes[0]].to_numpy(), 
                                              self.data[axes[1]].to_numpy(), 
                                              self.data[axes[2]].to_numpy()]), axis=0)
      self.data['angle'] = np.arccos(self.data[angle_axis] /  momentumTotal)