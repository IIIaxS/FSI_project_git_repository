from numpy import *
from .variables import *
from .model_part import *

# ----------------------------------------------------------------------
#author  : martin.ruchti@tum.de
#
# comment : modifications for the existing gid_io.py have been done
#           to be able to read in GiD geometry .mdpa files 
#
#notes: this changes have to be done currently manually to the .mdpa file
#inlet velocity has to be written as:
#            Begin NodalData VELOCITY_X Inlet
#            Begin NodalData VELOCITY_Y Inlet
#wall velocity has to be written as:
#            Begin NodalData VELOCITY_X Wall
#            Begin NodalData VELOCITY_Y Wall
#             
#the movable boundary has to be written as: 
#            Begin NodalData VELOCITY_X IS_BOUNDARY
#            Begin NodalData VELOCITY_Y IS_BOUNDARY
#the structure boundary for drag force computation has to be written as:
#            Begin NodalData VELOCITY_X IS_STRUCTURE
#            Begin NodalData VELOCITY_Y IS_STRUCTURE
#
#additionally, the following changes have to be done in this script
#    line: 207 element type
#    
# ----------------------------------------------------------------------

class GidIONS:
    def __init__(self,inputfilename, outputfilename,zero_based_indices_for_nodes=False):
        try:
            self.input_file = open( inputfilename)
            print('tried input file')
        except:
            self.input_file = 0
            print('failed')
        self.mesh_file = open( outputfilename + ".post.msh", 'w')
        self.result_file = open( outputfilename + ".post.res", 'w')
        self.result_file.write("GiD Post Results File 1.0\n")
        self.zero_based_indices_for_nodes = zero_based_indices_for_nodes
        
    def ReadModelPart(self, model_part):
        if self.input_file:
            for line in self.input_file:
                if line.find("Begin Nodes") != -1:
                    #print('found in line1:',line.find("Begin Nodes"))
                    words = self.ReadWords(line)
                    #print(words)  
                    if words[1] == "Nodes":
                        print('in coord')
                        self.ReadNodes(model_part)
                        
                elif line.find("Begin Elements") != -1:
                    #print('found in line1:',line.find("Begin Elements"))
                #if line.find("MESH") != -1:
                    words = self.ReadWords(line)
                    #print(words)
                    if words[1] == "Elements":
                        print('in elements')
                        self.ReadElements(words[2],model_part)
                        
                elif line.find("Begin Conditions") != -1:
                    #print('found in line1:',line.find("Begin Conditions"))
                    words = self.ReadWords(line)
                    #print(words)
                    if words[1] == "Conditions":
                        print('in conditions')
                        self.ReadConditions(words[2],model_part)
                        
                # fixing nodal values on wall
                elif line.find("Begin NodalData VELOCITY_X Init") != -1:
                    #print('found in line1:',line.find("Begin Conditions"))
                    words = self.ReadWords(line)
                    #print(words)
                    if ((words[2] == "VELOCITY_X") and (words[3] =="Init")):
                        print('in vel_x Init')
                        self.ReadNodalDataInit(words[2],model_part)
                        
                # fixing nodal values on wall
                elif line.find("Begin NodalData VELOCITY_Y Init") != -1:
                    #print('found in line1:',line.find("Begin Conditions"))
                    words = self.ReadWords(line)
                    #print(words)
                    if ((words[2] == "VELOCITY_Y") and (words[3] =="Init")):
                        print('in vel_y Init')
                        self.ReadNodalDataInit(words[2],model_part)
                        
                # fixing nodal values on wall
                elif line.find("Begin NodalData VELOCITY_X Wall") != -1:
                    #print('found in line1:',line.find("Begin Conditions"))
                    words = self.ReadWords(line)
                    #print(words)
                    if ((words[2] == "VELOCITY_X") and (words[3] =="Wall")):
                        print('in vel_x wall')
                        self.ReadNodalDataWall(words[2],model_part)
                        
                elif line.find("Begin NodalData VELOCITY_Y Wall") != -1:
                    #print('found in line1:',line.find("Begin Conditions"))
                    words = self.ReadWords(line)
                    #print(words)
                    if ((words[2] == "VELOCITY_Y") and (words[3] =="Wall")):
                        print('in vel_y wall')
                        self.ReadNodalDataWall(words[2],model_part)
                        
                elif line.find("Begin NodalData VELOCITY_X IS_BOUNDARY") != -1:
                    #print('found in line1:',line.find("Begin Conditions"))
                    words = self.ReadWords(line)
                    #print(words)
                    if ((words[2] == "VELOCITY_X") and (words[3] =="IS_BOUNDARY")):
                        print('boundary x')
                        self.ReadNodalDataBoundary(words[2],model_part)
                        
                elif line.find("Begin NodalData VELOCITY_Y IS_BOUNDARY") != -1:
                    #print('found in line1:',line.find("Begin Conditions"))
                    words = self.ReadWords(line)
                    #print(words)
                    if ((words[2] == "VELOCITY_Y") and (words[3] =="IS_BOUNDARY")):
                        print('boundary y')
                        self.ReadNodalDataBoundary(words[2],model_part)
                        
                elif line.find("Begin NodalData VELOCITY_X IS_STRUCTURE") != -1:
                    #print('found in line1:',line.find("Begin Conditions"))
                    words = self.ReadWords(line)
                    #print(words)
                    if ((words[2] == "VELOCITY_X") and (words[3] =="IS_STRUCTURE")):
                        print('structure x')
                        self.ReadNodalDataStructure(words[2],model_part)
                        
                elif line.find("Begin NodalData VELOCITY_Y IS_STRUCTURE") != -1:
                    #print('found in line1:',line.find("Begin Conditions"))
                    words = self.ReadWords(line)
                    #print(words)
                    if ((words[2] == "VELOCITY_Y") and (words[3] =="IS_STRUCTURE")):
                        print('structure y')
                        self.ReadNodalDataStructure(words[2],model_part)
                        
                #on inlet   
                elif line.find("Begin NodalData VELOCITY_X Inlet") != -1:
                    #print('found in line1:',line.find("Begin Conditions"))
                    words = self.ReadWords(line)
                    #print(words)
                    if ((words[2] == "VELOCITY_X") and (words[3] =="Inlet")):
                        print('in vel_x inlet')
                        self.ReadNodalDataInlet(words[2],model_part)
                        
                elif line.find("Begin NodalData VELOCITY_Y Inlet") != -1:
                    #print('found in line1:',line.find("Begin Conditions"))
                    words = self.ReadWords(line)
                    #print(words)
                    if ((words[2] == "VELOCITY_Y") and (words[3] =="Inlet")):
                        print('in vel_y inlet')
                        self.ReadNodalDataInlet(words[2],model_part)
                        
                #on outlet
                elif line.find("Begin NodalData PRESSURE Outlet") != -1:
                    #print('found in line1:',line.find("Begin Conditions"))
                    words = self.ReadWords(line)
                    #print(words)
                    if ((words[2] == "PRESSURE") and (words[3] =="Outlet")):
                        print('in pres outlet')
                        self.ReadNodalDataOutlet(words[2],model_part)
          
    def ReadWords(self,line):
        i  = line.find("//")
        if i != -1:
            return line[:i].split() 
        return line.split() 
         
    def ReadProperties(self,properties):
        for line in self.input_file:
            if line.find("End") == -1:
                words = self.ReadWords(line)
                if words[1].find("[") ==-1: #scalar case
                    variable,value = words[0],float(words[1])
                    if variable not in properties:
                        properties.update({variable:value})
                    else:
                        properties[variable] = value
            else:
                break
                
    def ReadNodes(self,model_part):
        #print('in READ nodes')
        counter = 0
        for line in self.input_file:
            if line.find("End Nodes") == -1:
                counter += 1
                words = self.ReadWords(line)
                
                #-1 for being consequent - GiD generates numbering from 1, solving needs it from 0
                #format: node_id-1, array[node_coord_x,node_coord_y] ->-1 at ID to be as the algorithm
                model_part.AddNode(int(words[0])-1,array([float(words[1]), float(words[2])]))

            else:
                break
        print(counter)
                   
    def ReadElements(self,element_name, model_part):
        #print('in READ elements')
        counter = 0
        for line in self.input_file:
            if line.find("End Elements") == -1:
                counter += 1
                words = self.ReadWords(line)
                
                element_name="navier_stokes_element_2d"
                
                #-1 for being consequent - GiD generates numbering from 1, solving needs it from 0
                #format: elem_id, element = [0, [node1, node1, node3]] for 2 d 3 nodes for triangular elem -> -1 for being consequent
                model_part.AddElements(element_name,{int(words[0])-1:[int(words[1]), [int(words[2])-1, int(words[3])-1, int(words[4])-1]]})
            else:
                break
        print(counter)
                   
    def ReadConditions(self,condition_type,model_part):
        #print('in READ conditions')
        counter = 0
        for line in self.input_file:
            if line.find("End Conditions") == -1:
                counter += 1
                words = self.ReadWords(line)
                
                #-1 for being consequent - GiD generates numbering from 1, solving needs it from 0
                model_part.AddConditions(condition_type, {int(words[0])-1:[int(words[1]), [int(words[2])-1, int(words[3])-1]]})
            else:
                break
        print(counter)
                   
    def ReadNodalDataWall(self,variable,model_part):
        print(variable)
        counter = 0
        for line in self.input_file:
            if line.find("End NodalData") == -1:
                counter += 1
                words = self.ReadWords(line)
                
                #-1 for being consequent - GiD generates numbering from 1, solving needs it from 0         
                if (variable == "VELOCITY_Y"):
                    model_part.Nodes[int(words[0])-1].Fix(variable)
                if (variable == "VELOCITY_X"):
                    model_part.Nodes[int(words[0])-1].Fix(variable)
                    model_part.Nodes[int(words[0])-1].SetSolutionStepValue(variable,0,0.0)
            else:
                break
        print(counter)
        
    '''
        this function fixes the boundary and sets the sensitivity boundary
    '''
    def ReadNodalDataBoundary(self,variable,model_part):
        print(variable)
        counter = 0
        for line in self.input_file:
            if line.find("End NodalData") == -1:
                counter += 1
                words = self.ReadWords(line)
                
                #-1 for being consequent - GiD generates numbering from 1, solving needs it from 0         
                if (variable == "VELOCITY_Y"):
                    model_part.Nodes[int(words[0])-1].Fix(variable)
                    model_part.Nodes[int(words[0])-1].SetSolutionStepValue(IS_BOUNDARY, 0, True)
                if (variable == "VELOCITY_X"):
                    model_part.Nodes[int(words[0])-1].Fix(variable)
                    model_part.Nodes[int(words[0])-1].SetSolutionStepValue(IS_BOUNDARY, 0, True)
                    model_part.Nodes[int(words[0])-1].SetSolutionStepValue(variable,0,0.0)
            else:
                break
        print(counter)
        
    '''
        this function sets an unfixed initial condition on the velocity
    '''
    def ReadNodalDataInit(self,variable,model_part):
        print(variable)
        counter = 0
        for line in self.input_file:
            if line.find("End NodalData") == -1:
                counter += 1
                words = self.ReadWords(line)
                #read value from .mdpa file
                val = float(words[2])
                model_part.Nodes[int(words[0])-1].SetSolutionStepValue(variable,0,val)
            else:
                break
        print(counter)
        
    '''
        this function fixes the boundary and sets the flag for the drag force
    '''
    def ReadNodalDataStructure(self,variable,model_part):
        print(variable)
        counter = 0
        for line in self.input_file:
            if line.find("End NodalData") == -1:
                counter += 1
                words = self.ReadWords(line)
                val = float(words[2])
                
                #-1 for being consequent - GiD generates numbering from 1, solving needs it from 0         
                if (variable == "VELOCITY_Y"):
                    model_part.Nodes[int(words[0])-1].Fix(variable)
                    model_part.Nodes[int(words[0])-1].SetSolutionStepValue(IS_STRUCTURE, 0, True)
                    model_part.Nodes[int(words[0])-1].SetSolutionStepValue(variable,0,val)
                if (variable == "VELOCITY_X"):
                    model_part.Nodes[int(words[0])-1].Fix(variable)
                    model_part.Nodes[int(words[0])-1].SetSolutionStepValue(IS_STRUCTURE, 0, True)
                    model_part.Nodes[int(words[0])-1].SetSolutionStepValue(variable,0,val)
            else:
                break
        print(counter)
        
    def ReadNodalDataInlet(self,variable,model_part):
        print(variable)
        counter = 0
        for line in self.input_file:
            if line.find("End NodalData") == -1:
                counter += 1
                words = self.ReadWords(line)
                
                #-1 for being consequent - GiD generates numbering from 1, solving needs it from 0
                if model_part.Nodes[int(words[0])-1].Fix(variable) != True:
                    model_part.Nodes[int(words[0])-1].Fix(variable)
                    #read value from .mdpa file
                    val = float(words[2])
                    model_part.Nodes[int(words[0])-1].SetSolutionStepValue(variable,0,val)

            else:
                break
        print(counter)
        
    def ReadNodalDataOutlet(self,variable,model_part):
        print(variable)
        counter = 0
        for line in self.input_file:
            if line.find("End NodalData") == -1:
                counter += 1
                words = self.ReadWords(line)
                #read value from .mdpa file
                val = float(words[2])
                #fix pressure on outlet for channel flow out BC
                model_part.Nodes[int(words[0])-1].Fix(PRESSURE)
                model_part.Nodes[int(words[0])-1].SetSolutionStepValue(PRESSURE,0,val)

            else:
                break
        print(counter)
        
    def gid_id( self, kratos_id ):
        if self.zero_based_indices_for_nodes:
            return kratos_id+1
        else:
            return kratos_id
        
    def WriteMesh(self,model_part, mesh_name):
            self.mesh_file.write("MESH \"")
            self.mesh_file.write(mesh_name)
            self.mesh_file.write("\" dimension 3 ElemType Triangle Nnode 3\n")
            self.mesh_file.write("Coordinates\n")
            if 0 in model_part.Nodes:
                self.zero_based_indices_for_nodes = True
            for node in model_part.NodeIterators():
                if(len(node.coordinates) == 3):
                    self.mesh_file.write("{} {} {} {}\n".format(self.gid_id(node.Id), node.coordinates[0], node.coordinates[1], node.coordinates[2]))
                else:
                    self.mesh_file.write("{} {} {} {}\n".format( self.gid_id(node.Id) , node.coordinates[0], node.coordinates[1], 0.0))
            self.mesh_file.write("end coordinates\n")
            self.mesh_file.write("Elements\n")   
            for element in model_part.ElementIterators():
                self.mesh_file.write("{} {} {} {} 1\n".format(self.gid_id(element.Id), self.gid_id(element.geometry[0].Id), self.gid_id(element.geometry[1].Id), self.gid_id(element.geometry[2].Id) ))
            self.mesh_file.write("end elements\n")
            self.mesh_file.flush()

    def WriteNodalResults(self, variable, nodes, time):
        if isinstance(variable, list):
            self.result_file.write("Result \"")
            self.result_file.write(variable[0])
            self.result_file.write('" "pyKratos" {} vector OnNodes\n'.format(time))
            self.result_file.write("values\n")
            for node in nodes:
                self.result_file.write("{} {} {}\n".format(self.gid_id(node.Id), node.GetSolutionStepValue(variable[1],0), node.GetSolutionStepValue(variable[2],0)))
        else:
            self.result_file.write("Result \"")
            self.result_file.write(variable)
            self.result_file.write('" "pyKratos" {} scalar OnNodes\n'.format(time))
            self.result_file.write("values\n")
            for node in nodes:
                self.result_file.write("{} {}\n".format(self.gid_id(node.Id), node.GetSolutionStepValue(variable,0)))
        self.result_file.write("end values\n")
        self.result_file.flush()
        
    def WriteNodalResultsInTimeInterval(self, variableDictIn, nodes, timeDictIn):
        #write results for every timestep in the interval
        
        #wrap single variables and single times to a list
        variableDict    = []
        timeDict        = []
        #this is to access the right buffer positions
        timeCounter     = 1
        
        if type(timeDictIn) is list:
            #timeDictIn is a list
            for times in timeDictIn:
                timeDict.append(times)
                timeCounter += 1
        else:
            #timeDictIn is a single value
            timeDict.append(timeDictIn)
            #timeCounter+= 1
        
        if type(variableDictIn) is list:
            #variableDictIn is a list
            for vars in variableDictIn:
                variableDict.append(vars)
        else:
            #variableDictIn is a single value
            variableDict.append(variableDictIn)
            
        
        for time in timeDict:
            #decrement buffer counter for each timestep accessed
            #also decrement in the first access, since 0 is the buffer position for the actual value
            #but timeCounter = 1 if just one time step is to be outputted
            timeCounter -= 1
            for variable in variableDict:
                if isinstance(variable, list):
                    self.result_file.write("Result \"")
                    self.result_file.write(variable[0])
                    self.result_file.write('" "pyKratos" {} vector OnNodes\n'.format(time))
                    self.result_file.write("values\n")
                    for node in nodes:
                        self.result_file.write("{} {} {}\n".format(self.gid_id(node.Id), node.GetSolutionStepValue(variable[1],timeCounter), node.GetSolutionStepValue(variable[2],timeCounter)))
                else:
                    self.result_file.write("Result \"")
                    self.result_file.write(variable)
                    self.result_file.write('" "pyKratos" {} scalar OnNodes\n'.format(time))
                    self.result_file.write("values\n")
                    for node in nodes:
                        self.result_file.write("{} {}\n".format(self.gid_id(node.Id), node.GetSolutionStepValue(variable,timeCounter)))
                self.result_file.write("end values\n")
            self.result_file.flush()