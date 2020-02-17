from data.consts import parameters_geometry

class Geometry:
    def __init__(self, **kwargs):
        self.c_a = kwargs.get("c_a")
        self.l_a = kwargs.get("l_a")

        self.x_1 = kwargs.get("x_1")
        self.x_2 = kwargs.get("x_2")
        self.x_3 = kwargs.get("x_3")

        self.x_a = kwargs.get("x_a")

        self.h = kwargs.get("h")
        self.t_sk = kwargs.get("t_sk")
        self.t_sp = kwargs.get("t_sp")

        self.t_st = kwargs.get("t_st")
        self.h_st = kwargs.get("h_st")
        self.w_st = kwargs.get("w_st")
        self.n_st = kwargs.get("n_st")
        pass

    @property
    def MMoI(self):
        print(self.c_a) # we can access things now
        print(self.x_a)

        return self.c_a + self.x_a

    @property
    def centroid(self):
        return

    @property
    def shearcenter(self):
        return


if __name__ == "__main__": # is called when you run the script
    # call an instance of the class
    geo = Geometry(**parameters_geometry) 

    print(geo.MMoI)