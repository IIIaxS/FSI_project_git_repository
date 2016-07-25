class Node:
    def __init__(self, id_num, x, y):
        self.id = int(id_num)
        self.x = float(x)
        self.y = float(y)

    def print_id_val(self, val, output):
        print(str(self.id).rjust(4), "0".rjust(8), str(val).rjust(8),
              file = output)

class Element:
    def __init__(self, id_num, n1, n2, n3):
        self.id = int(id_num)
        self.nodes_id = [int(n1),
                         int(n2),
                         int(n3)]

class Cond:
    def __init__ (self, id_num, n1, n2):
        self.id = int(id_num)
        self.nodes_id = [int(n1),
                         int(n2)]

    def print_cond(self, output):
        print(str(self.id).rjust(4),
              "0".rjust(8),
              str(self.nodes_id[0]).rjust(8),
              str(self.nodes_id[1]).rjust(8),
              file = output)

def read_words(line):
    i  = line.find("//")
    if i != -1:
        return line[:i].split() 
    return line.split() 

class Model:
    def __init__(self, filename):
        self.nodes = [ ]
        self.elements = [ ]
        self.conditions = [ ]
        self.read_input(filename)

        self.output = open(filename, 'a')


    def read_input(self, inputfilename):
        input_file = open(inputfilename)
        if input_file:
            for line in input_file:
                
                if line.find("Begin Nodes") != -1:
                    words = read_words(line)
                    if words[1] == "Nodes":
                        for line in input_file:
                            if line.find("End Nodes") == -1:
                                words = read_words(line)
                                new_node = Node(words[0], words[1], words[2])
                                self.nodes.append(new_node)
                            else:
                                break

                if line.find("Begin Elements") != -1:
                    words = read_words(line)
                    if words[1] == "Elements":
                        for line in input_file:
                            if line.find("End Elements") == -1:
                                words = read_words(line)
                                new_el = Element(words[0], words[2], words[3], words[4])
                                self.elements.append(new_el)
                            else:
                                break
        input_file.close()
                            
    def set_initial_vel(self, value):
        # set initial value of all nodes
        print("Begin NodalData VELOCITY_X Init", file = self.output)
        for node in self.nodes:
            node.print_id_val(value, self.output)
        print("End NodalData", end='\n\n', file = self.output)

        print("Begin NodalData VELOCITY_Y Init", file = self.output)
        for node in self.nodes:
            node.print_id_val(value, self.output)
        print("End NodalData", end='\n\n', file = self.output)

    def set_wall_vertical(self, x, tol = 1e-4):
        print("Begin NodalData VELOCITY_X Wall", file = self.output)
        for node in self.nodes:
            if abs(node.x - x) < tol:
                print(node.id,  file = self.output)
        print("End NodalData", end='\n\n', file = self.output)
        
        print("Begin NodalData VELOCITY_Y Wall", file = self.output)
        for node in self.nodes:
            if abs(node.x - x) < tol:
                print(node.id, file = self.output)
        print("End NodalData", end='\n\n', file = self.output)

    def set_wall_horizontal(self, y, tol = 1e-4):
        print("Begin NodalData VELOCITY_X Wall", file = self.output)
        for node in self.nodes:
            if abs(node.y - y) < tol:
                print(node.id, file = self.output)
        print("End NodalData", end='\n\n', file = self.output)
        
        print("Begin NodalData VELOCITY_Y Wall", file = self.output)
        for node in self.nodes:
            if abs(node.y - y) < tol:
                print(node.id, file = self.output)
        print("End NodalData", end='\n\n', file = self.output)

    def set_structure(self, xmin, xmax, ymin, ymax, value = 0, tol = 1e-4):
        print("Begin NodalData VELOCITY_X IS_STRUCTURE", file = self.output)
        for node in self.nodes:
            if xmin - tol < node.x < xmax + tol\
            and ymin - tol < node.y < ymax + tol:
                node.print_id_val(value, self.output)
        print("End NodalData", end='\n\n', file = self.output)

        print("Begin NodalData VELOCITY_Y IS_STRUCTURE", file = self.output)
        for node in self.nodes:
            if xmin - tol < node.x < xmax + tol\
            and ymin - tol < node.y < ymax + tol:
                node.print_id_val(value, self.output)
        print("End NodalData", end='\n\n', file = self.output)

    def set_inlet_vertical(self, x, vx, tol = 1e-4):
        print("Begin NodalData VELOCITY_X Inlet", file = self.output)
        for node in self.nodes:
            if abs(node.x - x) < tol:
                node.print_id_val(vx, self.output)
        print("End NodalData", end='\n\n', file = self.output)
        
        print("Begin NodalData VELOCITY_Y Inlet", file = self.output)
        for node in self.nodes:
            if abs(node.x - x) < tol:
                node.print_id_val(0, self.output)
        print("End NodalData", end='\n\n', file = self.output)
        
    def set_inlet_horizontal(self, y, vy, tol = 1e-4):
        print("Begin NodalData VELOCITY_X Inlet", file = self.output)
        for node in self.nodes:
            if abs(node.y - y) < tol:
                node.print_id_val(0, self.output)
        print("End NodalData", end='\n\n', file = self.output)
        
        print("Begin NodalData VELOCITY_Y Inlet", file = self.output)
        for node in self.nodes:
            if abs(node.y - y) < tol:
                node.print_id_val(vy, self.output)
        print("End NodalData", end='\n\n', file = self.output)

    def set_condition_vertical(self, x,
                               condition = "robust_outflow_condition",
                               tol = 1e-4):
        print("Begin Conditions", condition, file = self.output)

        node_list = [ ]
        for node in self.nodes:
            if abs(node.x - x) < tol:
                node_list.append(node)
        node_list.sort(key = lambda n: n.x)
        
        counter = 1
        node_iter = iter(node_list)
        i = next(node_iter)
        for j in node_iter:
            new_cond = Cond(counter, i.id, j.id)
            self.conditions.append(new_cond)
            i = j
            counter += 1

        for cond in self.conditions:
            cond.print_cond(self.output)

        print("End Conditions", end = '\n\n', file = self.output)

model = Model("/Users/massimosferza/LRZ Sync+Share/TUM/TUM SoSe16/Courses/Fluid Structure Interaction/FSI_project_git_repository/our test/Test-robust - write s0 CLEAN/test2.mdpa")
model.set_condition_vertical(1.5)
model.set_initial_vel(0)
model.set_wall_horizontal(-10)
model.set_wall_horizontal(10)
model.set_structure(-0.5, 0.5,-0.5, 0.5, value = 0)
model.set_inlet_vertical(-10, 30.0)

model.output.close()
