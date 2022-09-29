import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RectBivariateSpline
import pathlib as pl
import time
import progressbar

class ASC_surface:
    # constructor with an asc file as input
    def __init__(self, asc_file):
        # storing the file name so the user knows from what file the surface is made
        self.file = asc_file

        # reading point cloud data
        points = np.loadtxt(asc_file)

        # sorting the points
        points = points[points[:, 1].argsort()]
        points = points[points[:, 0].argsort(kind='mergesort')]

        # making necessary variables for the coordinates
        self.x = np.unique(points[:, 0])
        self.y = np.unique(points[:, 1])
        self.z = np.reshape(points[:, 2], (self.x.shape[0], self.y.shape[0]))

        # translating the surface to origin
        self.x = np.arange(0.0, self.x.shape[0] * self.dx, self.dx)
        self.y = np.arange(0.0, self.y.shape[0] * self.dy, self.dy)

        # TODO: not understandable behaviour 51 - 50, temporary fix
        if self.z.shape[0] != self.x.shape[0]:
            self.x = self.x[0:self.z.shape[0]]
        if self.z.shape[1] != self.y.shape[0]:
            self.y = self.y[0:self.z.shape[1]]

    @property
    def dx(self):
        return abs(np.average(np.diff(self.x)))

    @property
    def dy(self):
        return abs(np.average(np.diff(self.y)))

    @property
    def size(self):
        return [self.x[-1], self.y[-1]]

    @property
    def shape(self):
        return (len(self.x), len(self.y))

    @property
    def mean(self):
        return np.mean(self.z)

    @property
    def min(self):
        return np.amin(self.z)

    @property
    def max(self):
        return np.amax(self.z)

    def __getitem__(self, key):
        return self.z[key]

    # print function for debugging
    def __str__(self):
        # string = "#" * 35 + "\n"
        string = "The surface information is: \n"
        x_range = [np.amin(self.x), np.amax(self.x)]
        string += ("xmin = " + str(x_range[0]) + "\t xmax = " + str(x_range[1]) + "\t Nx = " + str(len(self.x)) + "\t dx = " + str(
            abs(np.average(np.diff(self.x)))) + "\n")
        y_range = [np.amin(self.y), np.amax(self.y)]
        string += ("ymin = " + str(y_range[0]) + "\t ymax = " + str(y_range[1]) + "\t Ny = " + str(len(self.y)) + "\t dy = " + str(
            abs(np.average(np.diff(self.y)))) + "\n")
        string += ("z shape = " + str(self.z.shape) + "\n")
        string += "-" * 35
        return string

    def contour(self):
        """ plot a contour plot of the asc surface """
        X, Y = np.meshgrid(self.x, self.y)
        Z = self.z.transpose()
        fig, ax = plt.subplots()
        ax.contourf(X, Y, Z)
        plt.contourf(X, Y, Z)
        plt.colorbar()
        plt.show()

    def plot3D(self):
        """ plots a 3D plot of ASC_surface """

        X, Y = np.meshgrid(self.x, self.y)
        Z = self.z.transpose()

        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()

    # interpolates as per new x and y range
    def interpolate(self, new_x, new_y):
        """
        interpolates the z values to the new x and y grid values
        """

        # interpolate z values to the new axis values
        interpolation_function = RectBivariateSpline(self.y, self.x, self.z.transpose(), ky=2)
        self.z = np.array(interpolation_function(new_y, new_x)).transpose()

        # set the axis to new axis
        self.x = new_x
        self.y = new_y

    def id_maker(self, i, j, k):
        return 1 + i + j * self.shape[0] + k * self.shape[0] * self.shape[1]


class Node:
    def __init__(self, ID, x, y, z):
        self.ID = ID
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        string = str(self.ID) + "\t" + str(self.x) + "\t" + str(self.y) + "\t" + str(self.z)
        return string

    def __getitem__(self):
        return self


class Element:
    def __init__(self, ID, A, B, C, D, E, F, G, H):
        """
        Element is a list of node ID's in a specific order
        """
        self.ID = ID
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E
        self.F = F
        self.G = G
        self.H = H

    @property
    def node_IDs(self):
        return [self.A, self.B, self.C, self.D, self.E, self.F, self.G, self.H]


class Material:
    def __init__(self, density=8050e-18, E=1000, mu=0.3):
        """
        defines a linear elastic material
        CAUTION : PLEASE TAKE CARE TO ASSIGN VALUES IN THE RIGHT UNITS ; Default values are in kg, micrometers
        :param density: density of the surface part
        :param E: Young's modulus
        :param mu: Poisson's ratio
        """
        self.density = density
        self.E = E
        self.mu = mu


class VolumeMesh:
    def __init__(self, surface, layers=5, direction=-1, material=Material()):
        """
        constructor for VolumeMesh class
        :param surface: ASC_surface from which volume would be made
        :param layers: number of layers in Z direction
        :param direction: direction in which the layers should be added.  Accepted values -1 or 1 representing +Z and  -Z directions.
        """
        # ensuring the right direction value is given
        while True:
            if not (direction in [-1, 1]):
                direction = int(
                    input("Invalid direction value.\ndirection can only be +1 or -1 indicating +Z or -Z direction respectively.  Please choose a direction."))
            else:
                break

        # storing the file name so the user knows from what file the surface is made
        self.file = surface.file

        # storing material of the surface
        self.material = material

        # surface dimensions setting
        self.x = surface.x
        self.y = surface.y
        self.dx = surface.dx
        self.dy = surface.dy
        self.min = surface.min
        self.max = surface.max
        self.mean = surface.mean
        start = surface.min if direction == -1 else surface.max
        self.z = np.arange(start, start + direction * layers * self.dz, direction * self.dz)

        # making nodes
        bar = progressbar.ProgressBar(max_value=np.prod(self.shape))
        self.nodes = np.array([], dtype=Node)
        ID = 1
        for z in self.z:
            for y in self.y:
                for x in self.x:
                    if z == start:
                        i = np.where(self.x == x)[0][0]
                        j = np.where(self.y == y)[0][0]

                        node = Node(ID, x, y, surface.z[i, j])
                    else:
                        node = Node(ID, x, y, z)
                    self.nodes = np.append(self.nodes, node)
                    # print(str(ID) + "/" + str(np.prod(self.shape)) + " nodes created")
                    bar.update(ID)
                    ID += 1

        # making elements
        n_elements = (self.shape[0] - 1) * (self.shape[1] - 1) * layers
        bar = progressbar.ProgressBar(max_value=n_elements)
        ID = 1
        self.elements = np.array([], dtype=Element)
        for node in self.nodes:

            # finding i, j, k using the formula node ID = (1 + i)+ (j * Nx) + (k * Nx * Ny)
            k = int((node.ID - 1) / (self.shape[0] * self.shape[1]))
            j = int(((node.ID - 1) / self.shape[0]) % self.shape[1])
            i = int((node.ID - 1) % self.shape[0])

            # if node is boundary node
            if (i == self.shape[0] - 1) or (j == self.shape[1] - 1) or (k == self.shape[2] - 1):
                continue

            # if not make an element
            A = node.ID
            B = self.get_node(node.ID + self.shape[0]).ID
            C = self.get_node(node.ID + self.shape[0] + 1).ID
            D = self.get_node(node.ID + 1).ID
            # print(A.ID, B.ID, C.ID, D.ID)

            next_layer_ID = int(1 + i + j * self.shape[0] + (k + 1) * self.shape[0] * self.shape[1])
            E = self.get_node(next_layer_ID).ID
            F = self.get_node(next_layer_ID + self.shape[0]).ID
            G = self.get_node(next_layer_ID + self.shape[0] + 1).ID
            H = self.get_node(next_layer_ID + 1).ID

            self.elements = np.append(self.elements, Element(ID, A, B, C, D, E, F, G, H))
            bar.update(ID)
            ID += 1

    @property
    def dz(self):
        return min(self.dx, self.dy)

    @property
    def shape(self):
        return tuple([len(self.x), len(self.y), len(self.z)])

    @property
    def node_sets(self):
        '''
        makes necessary node sets for simulation purposes
        :return:
        '''
        node_sets = {'rough': [], 'smooth': [], 'Ymin': [], 'Ymax': [], 'Xmin': [], 'Xmax': []}
        for node in self.nodes:
            # finding i, j, k using the formula node ID = (1 + i)+ (j * Nx) + (k * Nx * Ny)
            k = int((node.ID - 1) / (self.shape[0] * self.shape[1]))
            j = int(((node.ID - 1) / self.shape[0]) % self.shape[1])
            i = int((node.ID - 1) % self.shape[0])

            if self.z[-1] - self.z[0] < 0:  # if direction==-1
                if k == 0:
                    node_sets["rough"].append(node.ID)

                if k == self.shape[2] - 1:
                    node_sets["smooth"].append(node.ID)
            else:  # if direction==+1
                if k == self.shape[2] - 1:
                    node_sets["smooth"].append(node.ID)

                if k == 0:
                    node_sets["rough"].append(node.ID)

            if j == 0:
                node_sets["Ymin"].append(node.ID)

            if j == self.shape[1] - 1:
                node_sets["Ymax"].append(node.ID)

            if i == 0:
                node_sets["Xmin"].append(node.ID)

            if i == self.shape[0] - 1:
                node_sets["Xmax"].append(node.ID)

        return node_sets

    @property
    def surface(self):
        surface = []
        for node in self.nodes:
            # finding i, j, k using the formula node ID = (1 + i)+ (j * Nx) + (k * Nx * Ny)
            k = int((node.ID - 1) / (self.shape[0] * self.shape[1]))

            if k == 0:
                for element in self.elements:
                    if node.ID in element.node_IDs:
                        surface.append(element.ID)
        surface = list(dict.fromkeys(surface))
        return surface

    def plot(self):
        x = []
        y = []
        z = []
        for node in self.nodes:
            x.append(node.x)
            y.append(node.y)
            z.append(node.z)
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(x, y, z, c='r', marker='o')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()

    def get_node(self, ID):
        """
        :param ID:  node ID
        :return: node object fron self.nodes list
        """
        for node in self.nodes:
            if node.ID == ID:
                return node
        print("Node with ID = ", ID, " does not exist")
        return []

    def write_inp_file(self, file_path):
        '''
        writes input file for VolumeMesh object
        :param file_path
        :return: input file written to the file path
        '''

        # writing node data
        with file_path.open('w') as f:
            # writing the heading
            f.write("** File manually created by a python program **\n")
            now = time.strftime("%c")
            f.write("** Created on : " + str(now) + " **\n")
            f.write("** max_height : " + str(self.max) + "\n")
            f.write("*HEADING\n*NODE\n")

            # print("Writing node data ...")
            for node in self.nodes:
                f.write(str(node.ID) + ",\t" + str(node.x) + ",\t" + str(node.y) + ",\t" + str(node.z) + "\n")

        # writing elements data
        # print("Writing elements data...")
        with file_path.open('a') as f:
            f.write("**\n** SOLID ELEMENTS\n**\n")
            f.write("*ELEMENT, TYPE=C3D8, ELSET=" + file_path.stem + "\n")
            for element in self.elements:
                f.write(str(element.ID) + ",\t" + str(element.A) + ",\t" + str(element.B) + ",\t" + str(element.C) + ",\t" + str(element.D) + ",\t" + str(
                    element.E) + ",\t" + str(element.F) + ",\t" + str(element.G) + ",\t" + str(element.H) + "\n")

        # writng node sets data
        # print("Writing node sets data...")
        with file_path.open('a') as f:
            f.write("**\n** NODE SETS\n**\n")
            for name, nodes in self.node_sets.items():
                f.write("*NSET, NSET=" + file_path.stem + "_" + name + "\n")
                for i in range(0, len(nodes)):
                    # ABAQUS only reads 9 values in one line
                    if i % 9 == 0 and i != 0:
                        f.write(str(nodes[i]) + ",\n")
                    elif i == len(nodes) - 1:
                        f.write(str(nodes[i]) + "\n")
                    else:
                        f.write(str(nodes[i]) + ",")

        # writing surfaces
        with file_path.open('a') as f:
            f.write("**\n** SURFACE DEFINITIONS\n**\n")
            f.write("*SURFACE, NAME=" + file_path.stem + "_roughFacets, TYPE=ELEMENT\n")
            for elementID in self.surface:
                f.write(str(elementID) + ", " + str("S1\n"))

        # writing the section data
        # print("Writing section data...")
        with file_path.open('a') as f:
            f.write("**\n** SECTION DATA\n**\n")
            f.write("*SOLID SECTION, ELSET=" + file_path.stem + ", MATERIAL=" + file_path.stem + "_mat\n")

        # writing material data
        # print("Writing material data...")
        with file_path.open('a') as f:
            f.write("**\n** MATERIAL DATA\n**\n")
            f.write("*MATERIAL, NAME=" + file_path.stem + "_mat\n")
            f.write("*DENSITY\n" + str(self.material.density) + ",\n")
            f.write("*ELASTIC\n")
            f.write(str(self.material.E) + ", " + str(self.material.mu))

def surfaces_reinterpolation(surface1, surface2):
    """
    Make the two surfaces of equal size and interval
    :param surface1: Class asc_surface object
    :param surface2: Class asc_surface object
    :return: None
    """
    # finding a common surface sizes
    ideal_x_max = min(surface1.size[0], surface2.size[0])
    ideal_y_max = min(surface1.size[1], surface2.size[1])

    # cutting excess surface areas
    surface1.x = surface1.x[surface1.x <= ideal_x_max]
    surface1.y = surface1.y[surface1.y <= ideal_y_max]
    surface1.z = surface1.z[0:surface1.shape[0], 0:surface1.shape[1]]

    surface2.x = surface2.x[surface2.x <= ideal_x_max]
    surface2.y = surface2.y[surface2.y <= ideal_y_max]
    surface2.z = surface2.z[0:surface2.shape[0], 0:surface2.shape[1]]

    # finding ideal surface intervals
    ideal_dx = min(surface1.dx, surface2.dx)
    ideal_dy = min(surface1.dy, surface2.dy)

    # making the ideal axis
    ideal_x = np.arange(0.0, ideal_x_max + 0.1 * ideal_dx, ideal_dx)
    ideal_y = np.arange(0.0, ideal_y_max + 0.1 * ideal_dy, ideal_dy)

    # interpolating as per new axis
    surface1.interpolate(ideal_x, ideal_y)
    surface2.interpolate(ideal_x, ideal_y)


def appendFiles(file1, file2):
    """
    Appends file2 line into file1
    :param file1: pathlib path of file1
    :param file2: pathlib path of file2
    :return: None
    """
    with file2.open('r') as f:
        lines = f.readlines()

    with file1.open('a') as f:
        for line in lines:
            f.write(line)