import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops
import time
from math import atan2,pi,pow,sqrt
from mpl_toolkits.mplot3d import Axes3D
from PyPonding import PondingLoadManager3d

def camber(xi,L,c):
    if c == 0:
        z = 0.
    else:
        r = (c**2 + (L/2)**2)/(2*c)
        z = sqrt(r**2 - (xi*L - L/2)**2) + (c-r)
    return z

class AnalysisResults:
    print_each_analysis_time_increment = True
    total_analysis_time = 0.
    
    def __init__(self):
        pass
        
    def add_to_analysis_time(self,tic,toc):
        self.total_analysis_time += toc - tic
        if self.print_each_analysis_time_increment:
            print(f"Adding {toc - tic:0.4f} seconds to total analysis run time.")

    def print_total_analysis_time(self):
        print(f"Total analysis time: {self.total_analysis_time:0.4f} seconds")

class ExampleRoof:

    # Bay Widths
    L_AB = 480
    L_BC = 480
    L_CD = 480
    L_12 = 480
    L_23 = 480

    # Roof Elevations at Columns
    z_A1 =   0
    z_A2 =  -8
    z_A3 =   0
    z_B1 =   0
    z_B2 = -10
    z_B3 =   0
    z_C1 =   0
    z_C2 = -10
    z_C3 =   0
    z_D1 =   0
    z_D2 =  -8
    z_D3 =   0

    # Joist and Joist Girder Section Properties
    E = 29000

    A_J   = 100.
    Iz_J  = 215.
    Iy_J  = 100.
    GJ_J  = 100.

    A_JG  = 100.
    Iz_JG = 2029.
    Iy_JG = 100.
    GJ_JG = 100.

    # Other Properties
    c_J     = 1.
    c_JG    = 1.
    nspaces = 8

    # Loading Parameters
    wdJG    = 50./1000/12 # Joist Girder self-weight (forcer per unit length, downward positive)
    qD      = 10./1000/12**2 # Uniform dead load (forcer per unit area, downward positive)
    gamma   = 62.4/1000/12**3
    zw      = -6

    # Analysis Options
    use_CBDI = False
    include_ponding_effect = True
    ndiv_J  = 10
    na      = 4
    nb      = 4
    element_type = 'dispBeamColumn'
    num_steps_zw = 1000
    # @todo add analysis options
    #   1. path analysis, ramping up volume and simple step incremental
    #   2. lumped analysis, going directly to zw and iterating

    # Analysis Tags
    girder_transf_tag   = 1
    girder_section_tag  = 1
    girder_beamint_tag  = 1

    joist_transf_tag    = 2
    joist_section_tag   = 2
    joist_beamint_tag   = 2

    dead_load_pattern_tag = 1
    dead_load_ts_tag      = 1

    ponding_load_pattern_tag_start = 2
    ponding_load_ts_tag = 2

    # Run Options
    plot_load_cells = False

    def __init__(self):
        pass

    def lowest_point(self):
        zo = float('inf')
        for i in ops.getNodeTags():
            iz = ops.nodeCoord(i, 3) + ops.nodeDisp(i, 3)
            if iz < zo:
                zo = iz
        return zo

    @property
    def nele_J(self):
        if self.use_CBDI:
            return 1
        else:
            return self.ndiv_J

    def BuildModel(self):

        z_position = dict()


        ops.wipe()
        ops.model('basic', '-ndm', 3, '-ndf', 6)

        # Define objects for girders
        girder_transf  = ops.geomTransf('Linear', self.girder_transf_tag, 0.0, -1.0, 0.0)
        girder_section = ops.section('Elastic', self.girder_section_tag, self.E, self.A_JG, self.Iz_JG, self.Iy_JG, 1, self.GJ_JG)
        girder_beamint = ops.beamIntegration('Lobatto', self.girder_beamint_tag, self.girder_section_tag, 4)

        joist_transf   = ops.geomTransf('Linear', self.joist_transf_tag, 1.0, 0.0, 0.0)
        joist_section  = ops.section('Elastic', self.joist_section_tag, self.E, self.A_J, self.Iz_J, self.Iy_J, 1, self.GJ_J)
        joist_beamint  = ops.beamIntegration('Lobatto', self.joist_beamint_tag, self.joist_section_tag, 4)


        ###########################################################
        # Define Joist Girders
        ###########################################################

        # Define Joist Girder B12
        for i in range(self.nspaces+1):
            n = 110101+i

            xi = (i/self.nspaces)
            L = self.L_12
            x = xi*L
            y = self.L_AB
            z = self.z_B1 + xi*(self.z_B2-self.z_B1) + camber(xi,L,self.c_JG)
            z_position[n] = z

            ops.node(n,x,y,z)
            if i == 0:
                ops.fix(n,1,1,1,1,0,1)
            elif i == self.nspaces:
                ops.fix(n,0,1,1,1,0,1)
            else:
                ops.fix(n,0,1,0,1,0,1)

        for i in range(self.nspaces):
            ops.element(self.element_type, 110101+i, 110101+i, 110102+i, self.girder_transf_tag, self.girder_beamint_tag)

        # Define Joist Girder B23
        for i in range(self.nspaces+1):
            n = 110201+i

            xi = (i/self.nspaces)
            L = self.L_23
            x = self.L_12 + xi*L
            y = self.L_AB
            z = self.z_B2 + xi*(self.z_B3-self.z_B2) + camber(xi,L,self.c_JG)
            z_position[n] = z

            ops.node(n,x,y,z)
            if i == 0:
                ops.fix(n,1,1,1,1,0,1)
            elif i == self.nspaces:
                ops.fix(n,0,1,1,1,0,1)
            else:
                ops.fix(n,0,1,0,1,0,1)

        for i in range(self.nspaces):
            ops.element(self.element_type, 110201+i, 110201+i, 110202+i, self.girder_transf_tag, self.girder_beamint_tag)

        # Define Joist Girder C12
        for i in range(self.nspaces+1):
            n = 120101+i

            xi = (i/self.nspaces)
            L = self.L_12
            x = xi*L
            y = self.L_AB + self.L_BC
            z = self.z_C1 + xi*(self.z_C2-self.z_C1) + camber(xi,L,self.c_JG)
            z_position[n] = z

            ops.node(n,x,y,z)
            if i == 0:
                ops.fix(n,1,1,1,1,0,1)
            elif i == self.nspaces:
                ops.fix(n,0,1,1,1,0,1)
            else:
                ops.fix(n,0,1,0,1,0,1)

        for i in range(self.nspaces):
            ops.element(self.element_type, 120101+i, 120101+i, 120102+i, self.girder_transf_tag, self.girder_beamint_tag)

        # Define Joist Girder C23
        for i in range(self.nspaces+1):
            n = 120201+i

            xi = (i/self.nspaces)
            L = self.L_23
            x = self.L_12 + xi*L
            y = self.L_AB + self.L_BC
            z = self.z_C2 + xi*(self.z_C3-self.z_C2) + camber(xi,L,self.c_JG)
            z_position[n] = z

            ops.node(n,x,y,z)
            if i == 0:
                ops.fix(n,1,1,1,1,0,1)
            elif i == self.nspaces:
                ops.fix(n,0,1,1,1,0,1)
            else:
                ops.fix(n,0,1,0,1,0,1)

        for i in range(self.nspaces):
            ops.element(self.element_type, 120201+i, 120201+i, 120202+i, self.girder_transf_tag, self.girder_beamint_tag)


        ###########################################################
        # Define Joists
        ##########################################################

        # Define joists between grid lines A and B
        for i in range(1,2*self.nspaces):

            if i < self.nspaces:
                nj = 110101+i
                zi = self.z_A1 + (i/self.nspaces)*(self.z_A2-self.z_A1)
                zj = z_position[nj]
                x = (i/self.nspaces)*self.L_12
            elif i == self.nspaces:
                zi = self.z_A2
                zj = self.z_B2
                x = self.L_12
            else:
                nj = 110201+(i-self.nspaces)
                zi = self.z_A2 + ((i-self.nspaces)/self.nspaces)*(self.z_A3-self.z_A2)
                zj = z_position[nj]
                x = self.L_12 + ((i-self.nspaces)/self.nspaces)*self.L_23


            for j in range(self.nele_J+1):
                n = 210001+100*i+j

                xi = j/self.nele_J
                L = self.L_AB
                y = xi*L
                z = zi + xi*(zj-zi) + camber(xi,L,self.c_J)

                ops.node(n,x,y,z)
                if j == 0:
                    ops.fix(n,1,1,1,0,1,1)
                elif j == self.nele_J:
                    if i == self.nspaces:
                        ops.fix(n,1,0,1,0,1,1)
                    else:
                        ops.fix(n,1,0,0,0,1,1)
                        ops.equalDOF(nj,n,3)
                else:
                    ops.fix(n,1,0,0,0,1,1)

            for j in range(self.nele_J):
                ops.element(self.element_type, 210001+100*i+j, 210001+100*i+j, 210002+100*i+j, self.joist_transf_tag, self.joist_beamint_tag)


        # Define joists between grid lines B and C
        for i in range(1,2*self.nspaces):

            if i < self.nspaces:
                ni = 110101+i
                nj = 120101+i
                zi = z_position[ni]
                zj = z_position[nj]
                x = (i/self.nspaces)*self.L_12
            elif i == self.nspaces:
                zi = self.z_B2
                zj = self.z_C2
                x = self.L_12
            else:
                ni = 110201+(i-self.nspaces)
                nj = 120201+(i-self.nspaces)
                zi = z_position[ni]
                zj = z_position[nj]
                x = self.L_12 + ((i-self.nspaces)/self.nspaces)*self.L_23


            for j in range(self.nele_J+1):
                n = 220001+100*i+j

                xi = j/self.nele_J
                L = self.L_BC
                y = self.L_AB + xi*L
                z = zi + xi*(zj-zi) + camber(xi,L,self.c_J)

                ops.node(n,x,y,z)
                if j == 0:
                    if i == self.nspaces:
                        ops.fix(n,1,1,1,0,1,1)
                    else:
                        ops.fix(n,1,1,0,0,1,1)
                        ops.equalDOF(ni,n,3)
                elif j == self.nele_J:
                    if i == self.nspaces:
                        ops.fix(n,1,0,1,0,1,1)
                    else:
                        ops.fix(n,1,0,0,0,1,1)
                        ops.equalDOF(nj,n,3)
                else:
                    ops.fix(n,1,0,0,0,1,1)

            for j in range(self.nele_J):
                ops.element(self.element_type, 220001+100*i+j, 220001+100*i+j, 220002+100*i+j, self.joist_transf_tag, self.joist_beamint_tag)


        # Define joists between grid lines C and D
        for i in range(1,2*self.nspaces):

            if i < self.nspaces:
                ni = 120101+i
                zi = z_position[ni]
                zj = self.z_D1 + (i/self.nspaces)*(self.z_D2-self.z_D1)
                x = (i/self.nspaces)*self.L_12
            elif i == self.nspaces:
                zi = self.z_C2
                zj = self.z_D2
                x = self.L_12
            else:
                ni = 120201+(i-self.nspaces)
                zi = z_position[ni]
                zj = self.z_D2 + ((i-self.nspaces)/self.nspaces)*(self.z_D3-self.z_D2)
                x = self.L_12 + ((i-self.nspaces)/self.nspaces)*self.L_23


            for j in range(self.nele_J+1):
                n = 230001+100*i+j

                xi = j/self.nele_J
                L = self.L_CD
                y = self.L_AB + self.L_BC + xi*L
                z = zi + xi*(zj-zi) + camber(xi,L,self.c_J)

                ops.node(n,x,y,z)
                if j == 0:
                    if i == self.nspaces:
                        ops.fix(n,1,1,1,0,1,1)
                    else:
                        ops.fix(n,1,1,0,0,1,1)
                        ops.equalDOF(ni,n,3)
                elif j == self.nele_J:
                    ops.fix(n,1,0,1,0,1,1)
                else:
                    ops.fix(n,1,0,0,0,1,1)

            for j in range(self.nele_J):
                ops.element(self.element_type, 230001+100*i+j, 230001+100*i+j, 230002+100*i+j, self.joist_transf_tag, self.joist_beamint_tag)


        ###########################################################
        # Define Ponding Load Cells
        ###########################################################

        PondingLoadManager = PondingLoadManager3d()

        # Define ponding load cells between A and B
        for i in range(1,2*self.nspaces+1):
            for j in range(1,self.ndiv_J+1):
                z_offsetI = 0.
                z_offsetJ = 0.
                z_offsetK = 0.
                z_offsetL = 0.

                # Vertex I
                if i == 1:
                    x = 0.
                    y = ((j-1)/self.ndiv_J)*self.L_AB
                    z = self.z_A1 + ((j-1)/self.ndiv_J)*(self.z_B1-self.z_A1)
                    vertexI = ('fixed',x,y,z)
                else:
                    if j == 1:
                        y = 0.
                        if i <= self.nspaces:
                            x = ((i-1)/self.nspaces)*self.L_12
                            z = self.z_A1 + ((i-1)/self.nspaces)*(self.z_A2-self.z_A1)
                        else:
                            x = self.L_12 + ((i-1-self.nspaces)/self.nspaces)*self.L_23
                            z = self.z_A2 + ((i-1-self.nspaces)/self.nspaces)*(self.z_A3-self.z_A2)
                        vertexI = ('fixed',x,y,z)
                    else:
                        if self.use_CBDI:
                            n = 210001 + 100*(i-1)
                            xi = (j-1)/self.ndiv_J
                            z_offsetI = camber(xi,self.L_AB,self.c_J)
                            vertexI = ('element',n,xi)
                        else:
                            n = 210001 + 100*(i-1) + (j-1)
                            vertexI = ('node',n)

                # Vertex J
                if i == 2*self.nspaces:
                    x = self.L_12 + self.L_23
                    y = ((j-1)/self.ndiv_J)*self.L_AB
                    z = self.z_A3 + ((j-1)/self.ndiv_J)*(self.z_B3-self.z_A3)
                    vertexJ = ('fixed',x,y,z)
                else:
                    if j == 1:
                        y = 0.
                        if i <= self.nspaces:
                            x = (i/self.nspaces)*self.L_12
                            z = self.z_A1 + (i/self.nspaces)*(self.z_A2-self.z_A1)
                        else:
                            x = self.L_12 + ((i-self.nspaces)/self.nspaces)*self.L_23
                            z = self.z_A2 + ((i-self.nspaces)/self.nspaces)*(self.z_A3-self.z_A2)
                        vertexJ = ('fixed',x,y,z)
                    else:
                        if self.use_CBDI:
                            n = 210001 + 100*i
                            xi = (j-1)/self.ndiv_J
                            z_offsetJ = camber(xi,self.L_AB,self.c_J)
                            vertexJ = ('element',n,xi)
                        else:
                            n = 210001 + 100*i + (j-1)
                            vertexJ = ('node',n)

                # Vertex K
                if i == 2*self.nspaces:
                    x = self.L_12 + self.L_23
                    y = (j/self.ndiv_J)*self.L_AB
                    z = self.z_A3 + (j/self.ndiv_J)*(self.z_B3-self.z_A3)
                    vertexK = ('fixed',x,y,z)
                else:
                    if j == self.ndiv_J:
                        if i <= self.nspaces:
                            n = 110101 + i
                        elif i == self.nspaces:
                            n = 210001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 110201 + (i-self.nspaces)
                        vertexK = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 210001 + 100*i
                            xi = j/self.ndiv_J
                            z_offsetK = camber(xi,self.L_AB,self.c_J)
                            vertexK = ('element',n,xi)
                        else:
                            n = 210001 + 100*i + j
                            vertexK = ('node',n)

                # Vertex L
                if i == 1:
                    x = 0.
                    y = (j/self.ndiv_J)*self.L_AB
                    z = self.z_A1 + (j/self.ndiv_J)*(self.z_B1-self.z_A1)
                    vertexL = ('fixed',x,y,z)
                else:
                    if j == self.ndiv_J:
                        if i <= self.nspaces:
                            n = 110101 + (i-1)
                        elif i == self.nspaces+1:
                            n = 210001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 110201 + (i-1-self.nspaces)
                        vertexL = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 210001 + 100*(i-1)
                            xi = j/self.ndiv_J
                            z_offsetL = camber(xi,self.L_AB,self.c_J)
                            vertexL = ('element',n,xi)
                        else:
                            n = 210001 + 100*(i-1) + j
                            vertexL = ('node',n)

                # add load cell
                id = 'AB_%03ix_%03iy' % (i,j)
                PondingLoadManager.add_cell(id,vertexI,vertexJ,vertexK,vertexL,self.gamma,self.na,self.nb)
                PondingLoadManager.cells[id].vertexI.z_offset = z_offsetI
                PondingLoadManager.cells[id].vertexJ.z_offset = z_offsetJ
                PondingLoadManager.cells[id].vertexK.z_offset = z_offsetK
                PondingLoadManager.cells[id].vertexL.z_offset = z_offsetL
                PondingLoadManager.cells[id].update_coord()


        # Define ponding load cells between B and C
        for i in range(1,2*self.nspaces+1):
            for j in range(1,self.ndiv_J+1):
                z_offsetI = 0.
                z_offsetJ = 0.
                z_offsetK = 0.
                z_offsetL = 0.

                # Vertex I
                if i == 1:
                    x = 0.
                    y = self.L_AB + ((j-1)/self.ndiv_J)*self.L_BC
                    z = self.z_B1 + ((j-1)/self.ndiv_J)*(self.z_C1-self.z_B1)
                    vertexI = ('fixed',x,y,z)
                else:
                    if j == 1:
                        if i <= self.nspaces:
                            n = 110101 + (i-1)
                        elif i == self.nspaces+1:
                            n = 210001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 110201 + (i-1-self.nspaces)
                        vertexI = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 220001 + 100*(i-1)
                            xi = (j-1)/self.ndiv_J
                            z_offsetI = camber(xi,self.L_BC,self.c_J)
                            vertexI = ('element',n,xi)
                        else:
                            n = 220001 + 100*(i-1) + (j-1)
                            vertexI = ('node',n)

                # Vertex J
                if i == 2*self.nspaces:
                    x = self.L_12 + self.L_23
                    y = self.L_AB + ((j-1)/self.ndiv_J)*self.L_BC
                    z = self.z_B3 + ((j-1)/self.ndiv_J)*(self.z_C3-self.z_B3)
                    vertexJ = ('fixed',x,y,z)
                else:
                    if j == 1:
                        if i <= self.nspaces:
                            n = 110101 + i
                        elif i == self.nspaces:
                            n = 210001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 110201 + (i-self.nspaces)
                        vertexJ = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 220001 + 100*i
                            xi = (j-1)/self.ndiv_J
                            z_offsetJ = camber(xi,self.L_BC,self.c_J)
                            vertexJ = ('element',n,xi)
                        else:
                            n = 220001 + 100*i + (j-1)
                            vertexJ = ('node',n)

                # Vertex K
                if i == 2*self.nspaces:
                    x = self.L_12 + self.L_23
                    y = self.L_AB + (j/self.ndiv_J)*self.L_BC
                    z = self.z_B3 + (j/self.ndiv_J)*(self.z_C3-self.z_B3)
                    vertexK = ('fixed',x,y,z)
                else:
                    if j == self.ndiv_J:
                        if i <= self.nspaces:
                            n = 120101 + i
                        elif i == self.nspaces:
                            n = 220001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 120201 + (i-self.nspaces)
                        vertexK = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 220001 + 100*i
                            xi = j/self.ndiv_J
                            z_offsetK = camber(xi,self.L_BC,self.c_J)
                            vertexK = ('element',n,xi)
                        else:
                            n = 220001 + 100*i + j
                            vertexK = ('node',n)

                # Vertex L
                if i == 1:
                    x = 0.
                    y = self.L_AB + (j/self.ndiv_J)*self.L_BC
                    z = self.z_B1 + (j/self.ndiv_J)*(self.z_C1-self.z_B1)
                    vertexL = ('fixed',x,y,z)
                else:
                    if j == self.ndiv_J:
                        if i <= self.nspaces:
                            n = 120101 + (i-1)
                        elif i == self.nspaces+1:
                            n = 220001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 120201 + (i-1-self.nspaces)
                        vertexL = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 220001 + 100*(i-1)
                            xi = j/self.ndiv_J
                            z_offsetL = camber(xi,self.L_BC,self.c_J)
                            vertexL = ('element',n,xi)
                        else:
                            n = 220001 + 100*(i-1) + j
                            vertexL = ('node',n)

                # add load cell
                id = 'BC_%03ix_%03iy' % (i,j)
                PondingLoadManager.add_cell(id,vertexI,vertexJ,vertexK,vertexL,self.gamma,self.na,self.nb)
                PondingLoadManager.cells[id].vertexI.z_offset = z_offsetI
                PondingLoadManager.cells[id].vertexJ.z_offset = z_offsetJ
                PondingLoadManager.cells[id].vertexK.z_offset = z_offsetK
                PondingLoadManager.cells[id].vertexL.z_offset = z_offsetL
                PondingLoadManager.cells[id].update_coord()


        # Define ponding load cells between C and D
        for i in range(1,2*self.nspaces+1):
            for j in range(1,self.ndiv_J+1):
                z_offsetI = 0.
                z_offsetJ = 0.
                z_offsetK = 0.
                z_offsetL = 0.

                # Vertex I
                if i == 1:
                    x = 0.
                    y = self.L_AB + self.L_BC + ((j-1)/self.ndiv_J)*self.L_CD
                    z = self.z_C1 + ((j-1)/self.ndiv_J)*(self.z_D1-self.z_C1)
                    vertexI = ('fixed',x,y,z)
                else:
                    if j == 1:
                        if i <= self.nspaces:
                            n = 120101 + (i-1)
                        elif i == self.nspaces+1:
                            n = 220001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 120201 + (i-1-self.nspaces)
                        vertexI = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 230001 + 100*(i-1)
                            xi = (j-1)/self.ndiv_J
                            z_offsetI = camber(xi,self.L_CD,self.c_J)
                            vertexI = ('element',n,xi)
                        else:
                            n = 230001 + 100*(i-1) + (j-1)
                            vertexI = ('node',n)

                # Vertex J
                if i == 2*self.nspaces:
                    x = self.L_12 + self.L_23
                    y = self.L_AB + self.L_BC + ((j-1)/self.ndiv_J)*self.L_CD
                    z = self.z_C3 + ((j-1)/self.ndiv_J)*(self.z_D3-self.z_C3)
                    vertexJ = ('fixed',x,y,z)
                else:
                    if j == 1:
                        if i <= self.nspaces:
                            n = 120101 + i
                        elif i == self.nspaces:
                            n = 220001 + 100*self.nspaces + self.nele_J
                        else:
                            n = 120201 + (i-self.nspaces)
                        vertexJ = ('node',n)
                    else:
                        if self.use_CBDI:
                            n = 230001 + 100*i
                            xi = (j-1)/self.ndiv_J
                            z_offsetJ = camber(xi,self.L_CD,self.c_J)
                            vertexJ = ('element',n,xi)
                        else:
                            n = 230001 + 100*i + (j-1)
                            vertexJ = ('node',n)

                # Vertex K
                if i == 2*self.nspaces:
                    x = self.L_12 + self.L_23
                    y = self.L_AB + self.L_BC + (j/self.ndiv_J)*self.L_CD
                    z = self.z_C3 + (j/self.ndiv_J)*(self.z_D3-self.z_C3)
                    vertexK = ('fixed',x,y,z)
                else:
                    if j == self.ndiv_J:
                        y = self.L_AB + self.L_BC + self.L_CD
                        if i <= self.nspaces:
                            x = (i/self.nspaces)*self.L_12
                            z = self.z_D1 + (i/self.nspaces)*(self.z_D2-self.z_D1)
                        else:
                            x = self.L_12 + ((i-self.nspaces)/self.nspaces)*self.L_23
                            z = self.z_D2 + ((i-self.nspaces)/self.nspaces)*(self.z_D3-self.z_D2)
                        vertexK = ('fixed',x,y,z)
                    else:
                        if self.use_CBDI:
                            n = 230001 + 100*i
                            xi = j/self.ndiv_J
                            z_offsetK = camber(xi,self.L_CD,self.c_J)
                            vertexK = ('element',n,xi)
                        else:
                            n = 230001 + 100*i + j
                            vertexK = ('node',n)

                # Vertex L
                if i == 1:
                    x = 0.
                    y = self.L_AB + self.L_BC + (j/self.ndiv_J)*self.L_CD
                    z = self.z_C1 + (j/self.ndiv_J)*(self.z_D1-self.z_C1)
                    vertexL = ('fixed',x,y,z)
                else:
                    if j == self.ndiv_J:
                        y = self.L_AB + self.L_BC + self.L_CD
                        if i <= self.nspaces:
                            x = ((i-1)/self.nspaces)*self.L_12
                            z = self.z_D1 + ((i-1)/self.nspaces)*(self.z_D2-self.z_D1)
                        else:
                            x = self.L_12 + ((i-1-self.nspaces)/self.nspaces)*self.L_23
                            z = self.z_D2 + ((i-1-self.nspaces)/self.nspaces)*(self.z_D3-self.z_D2)
                        vertexL = ('fixed',x,y,z)
                    else:
                        if self.use_CBDI:
                            n = 230001 + 100*(i-1)
                            xi = j/self.ndiv_J
                            z_offsetL = camber(xi,self.L_CD,self.c_J)
                            vertexL = ('element',n,xi)
                        else:
                            n = 230001 + 100*(i-1) + j
                            vertexL = ('node',n)

                # add load cell
                id = 'CD_%03ix_%03iy' % (i,j)
                PondingLoadManager.add_cell(id,vertexI,vertexJ,vertexK,vertexL,self.gamma,self.na,self.nb)
                PondingLoadManager.cells[id].vertexI.z_offset = z_offsetI
                PondingLoadManager.cells[id].vertexJ.z_offset = z_offsetJ
                PondingLoadManager.cells[id].vertexK.z_offset = z_offsetK
                PondingLoadManager.cells[id].vertexL.z_offset = z_offsetL
                PondingLoadManager.cells[id].update_coord()


        ###########################################################
        # Define Dead Load
        ###########################################################

        ops.timeSeries("Constant", self.dead_load_ts_tag)
        ops.pattern("Plain", self.dead_load_pattern_tag, self.dead_load_ts_tag)

        # Define uniform dead load on joists
        for i in range(1,2*self.nspaces):

            if i < self.nspaces:
                wD = self.qD*self.L_12/self.nspaces
            elif i == self.nspaces:
                wD = self.qD*(0.5*self.L_12/self.nspaces+0.5*self.L_23/self.nspaces)
            else:
                wD = self.qD*self.L_23/self.nspaces

            for j in range(self.nele_J):
                n = 210001+100*i+j
                nodes = ops.eleNodes(n)
                coordi = ops.nodeCoord(nodes[0])
                coordj = ops.nodeCoord(nodes[1])
                Ly = coordj[1]-coordi[1]
                Lz = coordj[2]-coordi[2]
                L = sqrt(Ly**2 + Lz**2)
                Wy = -wD*(Ly/L)*(Ly/L)
                Wx = -wD*(Ly/L)*(Lz/L)
                ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

                n = 220001+100*i+j
                nodes = ops.eleNodes(n)
                coordi = ops.nodeCoord(nodes[0])
                coordj = ops.nodeCoord(nodes[1])
                Ly = coordj[1]-coordi[1]
                Lz = coordj[2]-coordi[2]
                L = sqrt(Ly**2 + Lz**2)
                Wy = -wD*(Ly/L)*(Ly/L)
                Wx = -wD*(Ly/L)*(Lz/L)
                ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

                n = 230001+100*i+j
                nodes = ops.eleNodes(n)
                coordi = ops.nodeCoord(nodes[0])
                coordj = ops.nodeCoord(nodes[1])
                Ly = coordj[1]-coordi[1]
                Lz = coordj[2]-coordi[2]
                L = sqrt(Ly**2 + Lz**2)
                Wy = -wD*(Ly/L)*(Ly/L)
                Wx = -wD*(Ly/L)*(Lz/L)
                ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

        # Define Self Weight on Joist Girder
        for i in range(self.nspaces):
            n = 110101+i
            nodes = ops.eleNodes(n)
            coordi = ops.nodeCoord(nodes[0])
            coordj = ops.nodeCoord(nodes[1])
            Lx = coordj[0]-coordi[0]
            Lz = coordj[2]-coordi[2]
            L = sqrt(Lx**2 + Lz**2)
            Wy = -self.wdJG*(Lx/L)*(Lx/L)
            Wx = -self.wdJG*(Lx/L)*(Lz/L)
            ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

            n = 110201+i
            nodes = ops.eleNodes(n)
            coordi = ops.nodeCoord(nodes[0])
            coordj = ops.nodeCoord(nodes[1])
            Lx = coordj[0]-coordi[0]
            Lz = coordj[2]-coordi[2]
            L = sqrt(Lx**2 + Lz**2)
            Wy = -self.wdJG*(Lx/L)*(Lx/L)
            Wx = -self.wdJG*(Lx/L)*(Lz/L)
            ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

            n = 120101+i
            nodes = ops.eleNodes(n)
            coordi = ops.nodeCoord(nodes[0])
            coordj = ops.nodeCoord(nodes[1])
            Lx = coordj[0]-coordi[0]
            Lz = coordj[2]-coordi[2]
            L = sqrt(Lx**2 + Lz**2)
            Wy = -self.wdJG*(Lx/L)*(Lx/L)
            Wx = -self.wdJG*(Lx/L)*(Lz/L)
            ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

            n = 120201+i
            nodes = ops.eleNodes(n)
            coordi = ops.nodeCoord(nodes[0])
            coordj = ops.nodeCoord(nodes[1])
            Lx = coordj[0]-coordi[0]
            Lz = coordj[2]-coordi[2]
            L = sqrt(Lx**2 + Lz**2)
            Wy = -self.wdJG*(Lx/L)*(Lx/L)
            Wx = -self.wdJG*(Lx/L)*(Lz/L)
            ops.eleLoad('-ele', n, '-type', '-beamUniform', Wy, 0.0, Wx)

        ###########################################################
        # Run Analysis
        ###########################################################

        # Initilize data
        results = AnalysisResults()
        results.water_volume = np.zeros((self.num_steps_zw+1,1))
        results.water_level  = np.zeros((self.num_steps_zw+1,1))
        results.col_react_B2 = np.zeros((self.num_steps_zw+1,1))
        
        # Set up analysis
        ops.numberer("RCM")
        ops.constraints("Plain")
        ops.system("BandSPD")
        ops.test("NormUnbalance", 1.0e-4, 25, 1)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        tic = time.perf_counter()
        ops.analyze(1)
        toc = time.perf_counter()
        results.add_to_analysis_time(tic,toc)
        ops.reactions()

        #print(ops.systemSize())

        # Find lowest point
        zo = self.lowest_point()
    
        # Store Reuslts
        results.water_volume[0] = 0.
        results.water_level[0]  = zo
        results.col_react_B2[0] = self.ColumnReaction('B2')

        # Run Ponding Analysis
        ops.timeSeries("Constant", self.ponding_load_ts_tag)
        for iStep in range(1,self.num_steps_zw+1):

            # Update ponding load cells
            if self.include_ponding_effect:
                PondingLoadManager.update()

            # Compute load vector
            izw = zo + (iStep/self.num_steps_zw)*(self.zw-zo)
            (iV,idVdz) = PondingLoadManager.get_volume(izw)
            PondingLoadManager.compute_current_load_vector(izw)

            # Apply difference to model
            ops.pattern("Plain", self.ponding_load_pattern_tag_start+iStep, self.ponding_load_ts_tag)
            PondingLoadManager.apply_load_increment()
            PondingLoadManager.commit_current_load_vector()

            # Run analysis
            tic = time.perf_counter()
            ops.analyze(1)
            toc = time.perf_counter()
            results.add_to_analysis_time(tic,toc)
            ops.reactions()

            # Store Reuslts
            results.water_volume[iStep] = iV
            results.water_level[iStep]  = izw
            results.col_react_B2[iStep] = self.ColumnReaction('B2')

        results.print_total_analysis_time()

        # Plot Ponding Load Cells
        if self.plot_load_cells:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            for i in PondingLoadManager.cells:
                coordI = PondingLoadManager.cells[i].vertexI.coord()
                coordJ = PondingLoadManager.cells[i].vertexJ.coord()
                coordK = PondingLoadManager.cells[i].vertexK.coord()
                coordL = PondingLoadManager.cells[i].vertexL.coord()
                #if i=='BC_001x_001y':
                #    print(coordI)
                #    print(coordJ)
                #    print(coordK)
                #    print(coordL)
                X = np.array([[coordJ[0], coordK[0]], [coordI[0], coordL[0]]])
                Y = np.array([[coordJ[1], coordK[1]], [coordI[1], coordL[1]]])
                Z = np.array([[coordJ[2], coordK[2]], [coordI[2], coordL[2]]])
                surf = ax.plot_surface(X, Y, Z, linewidth=1)
            plt.show()
        
        return results

    def PlotModel(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        NodeCoord = dict()
        NodeTags = ops.getNodeTags()
        for i in NodeTags:
            coord = ops.nodeCoord(i)
            NodeCoord[i] = coord
            #ax.scatter(coord[0], coord[1], coord[2], marker='o')

        EleTags = ops.getEleTags()
        for i in EleTags:
            nodes = ops.eleNodes(i)
            coordi = NodeCoord[nodes[0]]
            coordj = NodeCoord[nodes[1]]
            ax.plot([coordi[0],coordj[0]],[coordi[1],coordj[1]],[coordi[2],coordj[2]])

        plt.show()


    def ColumnReaction(self,column):
        if column == 'B2':
            #print(ops.nodeReaction(110101+self.nspaces))
            #print(ops.nodeReaction(110201))
            #print(ops.nodeReaction(210001+100*self.nspaces+self.nele_J))
            #print(ops.nodeReaction(220001+100*self.nspaces))

            R = ops.nodeReaction(110101+self.nspaces, 3) + \
                ops.nodeReaction(110201, 3) + \
                ops.nodeReaction(210001+100*self.nspaces+self.nele_J, 3) + \
                ops.nodeReaction(220001+100*self.nspaces, 3)
        elif column == 'C2':
            R = ops.nodeReaction(120101+self.nspaces, 3) + \
                ops.nodeReaction(120201, 3) + \
                ops.nodeReaction(220001+100*self.nspaces+self.nele_J, 3) + \
                ops.nodeReaction(230001+100*self.nspaces, 3)
        else:
            raise Exception('Unknown column')

        return R

