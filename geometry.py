import numpy as np

class Geometry:
    def __init__(self, **kwargs):
        self.c_a = kwargs.get("c_a")
        self.l_a = kwargs.get("l_a")

        # self.x_1 = kwargs["x_1"]
        # self.x_2 = kwargs["x_2"]
        # self.x_3 = kwargs["x_3"]

        #TODO: add the rest of the geometry inputs
        pass

    @property
    def MMoI(self):
        return f"Here is the length of the aileron : {self.l_a}"

    @property
    def centroid(self):
        return

    @property
    def shearcenter(self):
        return


if __name__ == "__main__": # is called when you run the script
    params = {"l_a" : 1.2,
              "c_a" : 0.4}

    geo = Geometry(**params)

    print(geo.MMoI)