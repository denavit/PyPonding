import numpy as np
import matplotlib.pyplot as plt
import PyPonding.Marino as Marino
from libdenavit import GoalSeekMonotonic
from PyPonding.structures import IdealizedBay
from math import pi

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

class IdealizedJoistBay:
    """
    Class to perform ponding analyses on an idealized rectangular bay with joists and joist girders. 
    
    
    """
    
    def __init__(self, **attrs):
    
        # Members
        self.joist          = attrs['joist']
        self.joist_girder   = attrs['joist_girder']

        # Loads and load factors
        self.dead_load_uniform         = attrs['dead_load_uniform'] 
        self.dead_load_primary_member  = attrs['dead_load_primary_member'] # Self-weight of joist girder
        self.water_density             = attrs['water_density'] 
        self.snow_density              = attrs['snow_density'] 
        self.snow_height               = attrs['snow_height'] 
        self.alpha                     = attrs['alpha']
        self.load_factor_dead          = attrs['load_factor_dead']
        self.load_factor_ponding       = attrs['load_factor_ponding']
        self.load_factor_snow          = attrs['load_factor_snow']
        self.consider_snow_and_water_overlap = attrs.get('consider_snow_and_water_overlap', True)

        self.run_factored_analysis_after_ponding_analysis = attrs.get('run_factored_analysis_after_ponding_analysis', False)
        self.additional_load_factor_dead    = attrs.get('additional_load_factor_dead', 1.0)
        self.additional_load_factor_ponding = attrs.get('additional_load_factor_ponding', 1.0)
        self.additional_load_factor_snow    = attrs.get('additional_load_factor_snow', 1.0)

        # Top of roof elevation
        self.z_TL = attrs['z_TL']
        self.z_TR = attrs['z_TR']
        self.z_BL = attrs['z_BL']
        self.z_BR = attrs['z_BR']

        # Camber
        self.primary_member_camber_T = attrs['primary_member_camber_T']
        self.primary_member_camber_B = attrs['primary_member_camber_B']
        self.secondary_member_camber = attrs['secondary_member_camber']

        # Edge conditions ('mirrored', 'rigid', or '')
        self.edge_condition_L = attrs['edge_condition_L']
        self.edge_condition_R = attrs['edge_condition_R']
        self.edge_condition_T = attrs['edge_condition_T']
        self.edge_condition_B = attrs['edge_condition_B']

        # Material and section properties
        self.E  = attrs['E']
        self.Ap = attrs['Ap']
        self.As = attrs['As']

        # Analysis options
        self.analsis_engine         = attrs.get('analsis_engine', 'FE')
        self.include_ponding_effect = attrs.get('include_ponding_effect', True)
        self.num_ele_secondary      = attrs.get('num_ele_secondary', 20)
        self.num_subcell_X          = attrs.get('num_subcell_X', 4) # Number of ponding sub-cells along joist direction
        self.num_subcell_Y          = attrs.get('num_subcell_Y', 4) # Number of ponding sub-cells along joist girder direction       
        self.MAX_ITER               = attrs.get('MAX_ITER', 50) # Maximum number of ponding analysis iterations
        self.tol                    = attrs.get('tol', 0.00001) # Ponding analysis tolerance
        
    def IdealizedBayObject(self):
    
        inch = 1.0
        ft  = 12.0*inch
      
        input = {
            'primary_member_span': self.joist_girder.span_ft*ft,
            'secondary_member_span': self.joist.span_ft*ft,
            'number_of_joist_spaces': self.joist_girder.num_spaces,
            'dead_load_uniform': self.dead_load_uniform,
            'dead_load_primary_member': self.dead_load_primary_member,
            'water_density': self.water_density,  
            'snow_density': self.snow_density,
            'snow_height': self.snow_height,
            'alpha': self.alpha,
            'load_factor_dead':    self.load_factor_dead,
            'load_factor_ponding': self.load_factor_ponding,
            'load_factor_snow':    self.load_factor_snow,
            'consider_snow_and_water_overlap': self.consider_snow_and_water_overlap,
            'run_factored_analysis_after_ponding_analysis': self.run_factored_analysis_after_ponding_analysis,
            'additional_load_factor_dead': self.additional_load_factor_dead,
            'additional_load_factor_ponding': self.additional_load_factor_ponding,
            'additional_load_factor_snow': self.additional_load_factor_snow,
            'z_TL': self.z_TL,
            'z_TR': self.z_TR,
            'z_BL': self.z_BL,
            'z_BR': self.z_BR,
            'secondary_member_camber': self.secondary_member_camber,
            'primary_member_camber_T': self.primary_member_camber_T,
            'primary_member_camber_B': self.primary_member_camber_B,
            'edge_condition_L': self.edge_condition_L,
            'edge_condition_R': self.edge_condition_R,
            'edge_condition_T': self.edge_condition_T,
            'edge_condition_B': self.edge_condition_B,
            'E': self.E,
            'As': self.As,
            'Ap': self.Ap,
            'Is': self.Is(),
            'Ip': self.Ip(),
            'analsis_engine': 'FE',
        }

        return IdealizedBay(**input)
        
    def RunPondingAnalysis(self,water_level):
        # Run analysis and spit out strength ratios like in the roof bay analysis tool.
        
        bay = self.IdealizedBayObject()
        results = bay.Run_Analysis(water_level)
        max_SR = 0.0
        
        print('\n==== Member Output ====') 
        print('  Member         Max Shear    Strength Ratio    Max Moment    Strength Ratio')
        print('                  (kips)          (---)          (kip-ft)         (---)     ')
        
        # Joists
        for i in range(bay.number_of_joist_spaces+1):
            if i == 0 and bay.edge_condition_L == 'rigid':
                continue
            if i == bay.number_of_joist_spaces and bay.edge_condition_R == 'rigid':
                continue
            Vmax = max(results.secondary_members_bot_reaction[i],results.secondary_members_top_reaction[i])
            Mmax = max(results.secondary_members_moment[i,:])/12
            
            # To match the SJI roof bay analysis tool. Shear should be checked at 
            # memeber ends with the reaction and the middle of each element.
            x_shear = np.zeros(self.num_ele_secondary+2)
            V_shear = np.zeros(self.num_ele_secondary+2)
            x_shear[0] = 0
            V_shear[0] = results.secondary_members_bot_reaction[i]
            for iele in range(1,self.num_ele_secondary+1):
                x_shear[iele] = 0.5*(results.secondary_members_position[2*iele-2]+results.secondary_members_position[2*iele-1])
                V_shear[iele] = results.secondary_members_shear[i,2*iele-2]
            x_shear[self.num_ele_secondary+1] = self.joist.span_ft*12
            V_shear[self.num_ele_secondary+1] = -results.secondary_members_top_reaction[i]
            
            SRV  = max(self.joist.shear_strength_ratio(x_shear/12,1000*V_shear))
            SRM  = max(self.joist.moment_strength_ratio(results.secondary_members_position/12,1000/12*results.secondary_members_moment[i,:]))
            max_SR = max(max_SR,SRV,SRM)
            if max(SRV,SRM) <= 1.0:
                strength_result = 'OKAY'
            else:
                strength_result = 'NO GOOD'
            print(f'Secondary {i+1:<2d}   {Vmax:>8.2f}        {SRV:>8.3f}         {Mmax:>8.2f}       {SRM:>8.3f}         {strength_result}')
            
        # Top Joist Girder
        Vmax = max(abs(results.top_primary_member_shear))
        Mmax = max(results.top_primary_member_moment)/12
        SRV  = max(self.joist_girder.shear_strength_ratio(results.primary_members_position/12,1000*results.top_primary_member_shear))
        SRM  = max(self.joist_girder.moment_strength_ratio(results.primary_members_position/12,1000/12*results.top_primary_member_moment))
        max_SR = max(max_SR,SRV,SRM)
        if max(SRV,SRM) <= 1.0:
            strength_result = 'OKAY'
        else:
            strength_result = 'NO GOOD'
        print(f'Primary Top    {Vmax:>8.2f}        {SRV:>8.3f}         {Mmax:>8.2f}       {SRM:>8.3f}         {strength_result}')
        
        # Bottom Joist Girder
        Vmax = max(abs(results.bot_primary_member_shear))
        Mmax = max(results.bot_primary_member_moment)/12
        SRV  = max(self.joist_girder.shear_strength_ratio(results.primary_members_position/12,1000*results.bot_primary_member_shear))
        SRM  = max(self.joist_girder.moment_strength_ratio(results.primary_members_position/12,1000/12*results.bot_primary_member_moment))
        max_SR = max(max_SR,SRV,SRM)
        if max(SRV,SRM) <= 1.0:
            strength_result = 'OKAY'
        else:
            strength_result = 'NO GOOD'
        print(f'Primary Bot    {Vmax:>8.2f}        {SRV:>8.3f}         {Mmax:>8.2f}       {SRM:>8.3f}         {strength_result}')
        
        return max_SR
    
    def LimitPointPodningAnalysis(self):
        goal_seek = GoalSeekMonotonic(1.000,0.001,tolerance_mode='UnderOnly')
   
        max_steps = 20
        for i in range(max_steps):
            water_level = goal_seek.suggest_input()
            print(f'\nRunning at {water_level = :.3f} in.')
            SRmax = joist_bay.RunPondingAnalysis(water_level)
            result = goal_seek.add_and_check(water_level,SRmax)        
            if result:
                print(f'\nThe water rising capacity is {water_level:.3f} in.')
                break
            
            if i == (max_steps-1):
                raise Exception(f'Goal seek did not find a solution in {max_steps} steps')
        
        return water_level
        
    def Is(self):
        return self.joist.moment_of_inertia()/1.15

    def Ip(self):
        return self.joist_girder.moment_of_inertia()/1.15
    
    def Cp(self):
        Ls = self.joist.span_ft*12.0
        Lp = self.joist_girder.span_ft*12.0        
        Cp = (self.water_density*Ls*Lp**4)/(pi**4*self.E*self.Ip())
        #print(f'{Ls = }, {Lp = }, {Cp = }')
        return Cp
    
    def Cs(self):
        S  = self.joist_girder.span_ft*12.0/self.joist_girder.num_spaces
        Ls = self.joist.span_ft*12.0
        Cs = (self.water_density*S*Ls**4)/(pi**4*self.E*self.Is())
        #print(f'{S = }, {Ls = }, {Cs = }')
        return Cs
    
    def Appendix2_Simplified(self):
        pass_tf = (self.Cp() + 0.9*self.Cs() <= 0.25)
        return pass_tf
        
    def Appendix2_Improved(self,water_level):
        S  = self.joist_girder.span_ft*12.0/self.joist_girder.num_spaces
        Ls = self.joist.span_ft*12.0 
        w  = self.dead_load_uniform + \
            max(self.snow_height*self.snow_density,water_level*self.water_density)

        xp = (w*S*Ls+self.dead_load_primary_member*S)/self.joist_girder.P_kips
        if self.joist_girder.strength_type == 'ASD':
            Up = (0.8-0.6*xp)/(0.6*xp)
        elif self.joist_girder.strength_type == 'LRFD':
            Up = (0.8-0.9*xp)/(0.9*xp)
        else:
            raise Exception(f'Unknown joist girder strength type: {self.joist_girder.strength_type}')
            
        xs = (w*S)/(self.joist.wTL_plf/12000)
        if self.joist_girder.strength_type == 'ASD':
            Us = (0.8-0.6*xs)/(0.6*xs)
        elif self.joist_girder.strength_type == 'LRFD':
            Us = (0.8-0.9*xs)/(0.9*xs)
        else:
            raise Exception(f'Unknown joist strength type: {self.joist.strength_type}')        
        
        #print(f'{Up = :.3f}')
        #print(f'{Us = :.3f}')

        Cp = self.Cp()
        Cs = self.Cs()
        Up_min = Marino.Up_min(Cp,Cs)
        Us_min = Marino.Us_min(Cp,Cs)

        #print(f'Up_min = {Up_min:.3f}')
        #print(f'Us_min = {Us_min:.3f}')

        pass_tf = (Up >= Up_min) and (Us >= Us_min)
        max_diff = max(Up_min-Up,Us_min-Us)
        
        return pass_tf, max_diff

    def LimitPointPondingAnalysis(self):
        tol = 0.001
        goal_seek = GoalSeekMonotonic(1.0-0.5*tol,0.5*tol)
   
        max_steps = 20
        for i in range(max_steps):
            water_level = goal_seek.suggest_input()
            print(f'\nRunning at {water_level = :.3f} in.')
            SRmax = self.RunPondingAnalysis(water_level)
            result = goal_seek.add_and_check(water_level,SRmax)        
            if result:
                print(f'\nThe water rising capacity is {water_level:.3f} in.')
                break
            
            if i == (max_steps-1):
                raise Exception(f'Goal seek did not find a solution in {max_steps} steps')
        
        return water_level
        
    def LimitPointAppendix2(self):
        tol = 0.001
        goal_seek = GoalSeekMonotonic(0.0-0.5*tol,0.5*tol)
   
        max_steps = 20
        for i in range(max_steps):
            water_level = goal_seek.suggest_input()
            print(f'\nRunning at {water_level = :.3f} in.')
            pass_tf, max_diff = self.Appendix2_Improved(water_level)
            print(f'{pass_tf} {max_diff:.3f}')
            result = goal_seek.add_and_check(water_level,max_diff)        
            if result:
                print(f'\nThe water rising capacity is {water_level:.3f} in.')
                break
            
            if i == (max_steps-1):
                raise Exception(f'Goal seek did not find a solution in {max_steps} steps')
        
        return water_level
    
    
def run_example():   
    from libdenavit.design import OpenWebSteelJoist,JoistGirder

    # Define units
    inch = 1.0
    kip = 1.0
    lb  = kip/1000.0
    ft  = 12.0*inch
    in_per_ft = inch/ft
    ksi = kip/inch**2
    plf = lb/ft
    psf = lb/ft**2
    pcf = psf/ft
    kipft = kip*ft

    strength_type = 'ASD'
    joist = OpenWebSteelJoist(strength_type,45,334,198)
    joist.minimum_shear_reversal_strength_ratio = 0.125
    joist_girder = JoistGirder(strength_type,40,48,5,14.3)

    input = {
        'joist': joist,
        'joist_girder': joist_girder,
        'dead_load_uniform': 13.73*psf,
        'dead_load_primary_member': 40*plf,
        'water_density': 62.4*pcf,  
        'snow_density': 15.30*pcf,
        'snow_height': 0.5*0.7*(25.00*psf)/(15.30*pcf),
        'alpha': 1.0,
        'load_factor_dead':    1.0,
        'load_factor_ponding': 1.0,
        'load_factor_snow':    1.0,
        'z_TL':  0*inch,
        'z_TR': 20*inch,
        'z_BL':  0*inch,
        'z_BR': 20*inch,
        'secondary_member_camber': 0.625*inch,
        'primary_member_camber_T': 0.625*inch,
        'primary_member_camber_B': 0.625*inch,
        'edge_condition_L': 'mirrored',
        'edge_condition_R': 'rigid',
        'edge_condition_T': 'mirrored',
        'edge_condition_B': 'mirrored',
        'E': 29000*ksi,
        'As': 100*inch**2,
        'Ap': 100*inch**2,
        'analsis_engine': 'FE',
    }

    joist_bay = IdealizedJoistBay(**input)
    
    zw = 4.0
    joist_bay.RunPondingAnalysis(zw)
    
    #from pprint import pprint
    #pprint(vars(joist_bay))
    #pprint(vars(bay))
    #pprint(vars(results))
    return


if __name__ == "__main__":
    run_example()

   