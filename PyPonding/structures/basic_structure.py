from PyPonding import FE

class basic_structure:
    # Loads
    alpha   = 1
    LF_D    = 1.2 # Dead
    wd      = 10/1000/12**2
    LF_P    = 1.2 # Impounded Water
    gamma   = 62.4/1000/12**3        
    LF_S1   = 1.2 # Snow in Ponding Load Cell
    LF_S2   = 0.0 # Snow as Simple Load
    gammas  = 20/1000/12**3
    hs      = 12
    include_ponding_effect = True
    modified_rain_load = False
    
    def __init__(self):
        pass    
    
    def Run_To_Strength_Limit(self,start_level=None,max_level=None,min_level=None,incr=1,tol=0.0001,use_stored=True,use_sparse=False):
    
        if start_level is None:
            start_level = self.lowest_point()
        if max_level is None:
            max_level = start_level + 240
        if min_level is None:
            min_level = start_level - 48
        
        self.BuildModel();
        self.model.use_sparse_matrix_solver = use_sparse
        
        if self.include_ponding_effect:   
            if self.modified_rain_load:
                PA = FE.PondingAnalysis(self.model,'Modified_Rain_Load')
            else:
                PA = FE.PondingAnalysis(self.model,'Constant_Level')
        else:
            PA = FE.PondingAnalysis(self.model,'No_Ponding_Effect')
        
        PA.max_iterations_z = 100
        PA.use_stored_analysis = use_stored
        if use_stored:
            self.model.StoreAnalysis()
        
        # Run Initial Analysis
        level = start_level
        res = PA.run({'DEAD':self.alpha*self.LF_D,'SNOW':self.alpha*self.LF_S2},level)
        if res != 0:
            print('Not converged')
            return float('nan')
        (SR,SR_note) = self.Strength_Ratio(PA)
        print('Level = %7.4f, Strength Ratio = %10.7f (%s)' % (level,SR,SR_note))
        
        # First Level of Iteration
        if abs(SR-1) < tol:
            return level
        elif SR < 1:
            # If start level is okay then increase level
            while (level <= max_level):
                if SR > 1:
                    above_level   = level
                    above_SR      = SR
                    above_SR_note = SR_note
                    break
                else:
                    below_level   = level
                    below_SR      = SR
                    below_SR_note = SR_note
                    level = level + incr
                    
                res = PA.run({'DEAD':self.alpha*self.LF_D,'SNOW':self.alpha*self.LF_S2},level)
                if res != 0:
                    print('Not converged')
                    return float('nan')
                (SR,SR_note) = self.Strength_Ratio(PA)
                print('Level = %7.4f, Strength Ratio = %10.7f (%s)' % (level,SR,SR_note))                    
            else:
                print('Maximum water level reached')
                return float('nan')
        else:
            # If start level is overstressed then decrease level
            while (level >= min_level):
                if SR < 1:
                    below_level   = level
                    below_SR      = SR
                    below_SR_note = SR_note
                    break
                else:
                    above_level   = level
                    above_SR      = SR
                    above_SR_note = SR_note
                    level = level - incr
                    
                res = PA.run({'DEAD':self.alpha*self.LF_D,'SNOW':self.alpha*self.LF_S2},level)
                if res != 0:
                    print('Not converged')
                    return float('nan')
                (SR,SR_note) = self.Strength_Ratio(PA)
                print('Level = %7.4f, Strength Ratio = %10.7f (%s)' % (level,SR,SR_note))
            else:
                print('Maximum water level reached')
                return float('nan')        
                
        # Second Level of Iteration
        SR = 0
        while (abs(SR-1) > tol):
            level = below_level + (above_level-below_level)*(1-below_SR)/(above_SR-below_SR)
            res = PA.run({'DEAD':self.alpha*self.LF_D,'SNOW':self.alpha*self.LF_S2},level)
            if res != 0:
                print('Not converged')
                return float('nan')
            (SR,SR_note) = self.Strength_Ratio(PA)
            print('Level = %7.4f, Strength Ratio = %10.7f (%s)' % (level,SR,SR_note))
            
            if SR > 1:
                above_level   = level
                above_SR      = SR
                above_SR_note = SR_note
            else:
                below_level   = level
                below_SR      = SR
                below_SR_note = SR_note           
                
        return level
