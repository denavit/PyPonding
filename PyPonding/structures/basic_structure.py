from PyPonding import FE

class basic_structure:
    # Loads
    LF_D    = 1.2 # Dead
    wd      = 10/1000/12**2
    LF_P    = 1.2 # Impounded Water
    gamma   = 62.4/1000/12**3        
    LF_S1   = 1.2 # Snow in Ponding Load Cell
    LF_S2   = 0.0 # Snow as Simple Load
    gammas  = 20/1000/12**3
    hs      = 12
    
    def __init__(self):
        pass    
    
    def Run_To_Strength_Limit(self):
        datam = self.lowest_point()
        level = datam + 1
        max_level = datam + 10
        
        self.BuildModel();
        self.model.use_sparse_matrix_solver = True
        
        PA = FE.PondingAnalysis(self.model,'Constant_Level')
        while (level <= max_level):
            res = PA.run({'DEAD':self.LF_D,'SNOW':self.LF_S2},level)
            if res != 0:
                print('Not converged')
            (SR,SR_note) = self.Strength_Ratio(PA)
            print('Level = %6.3f, Strength Ratio = %8.5f (%s)' % (level,SR,SR_note))
            if SR > 1:
                break
            else:
                level = level + 1
        else:
            print('Maximum water level reached')
    
        level += -0.9
        for i in range(10):
            res = PA.run({'DEAD':self.LF_D,'SNOW':self.LF_S2},level)
            if res != 0:
                print('Not converged')
            (SR,SR_note) = self.Strength_Ratio(PA)
            print('Level = %6.3f, Strength Ratio = %8.5f (%s)' % (level,SR,SR_note))
            if SR > 1:
                break
            else:
                level = level + 0.1

        level += -0.09
        for i in range(10):
            res = PA.run({'DEAD':self.LF_D,'SNOW':self.LF_S2},level)
            if res != 0:
                print('Not converged')
            (SR,SR_note) = self.Strength_Ratio(PA)
            print('Level = %6.3f, Strength Ratio = %8.5f (%s)' % (level,SR,SR_note))
            if SR > 1:
                break
            else:
                level = level + 0.01
        
        level += -0.009
        for i in range(10):
            res = PA.run({'DEAD':self.LF_D,'SNOW':self.LF_S2},level)
            if res != 0:
                print('Not converged')
            (SR,SR_note) = self.Strength_Ratio(PA)
            print('Level = %6.3f, Strength Ratio = %8.5f (%s)' % (level,SR,SR_note))
            if SR > 1:
                break
            else:
                level = level + 0.001
                
        #self.plot_results(PA)                
