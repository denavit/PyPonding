from math import pi as π

def Up_min(Cp,Cs):
    αp = Cp/(1-Cp)
    αs = Cs/(1-Cs)
    ρ = Cs/Cp
    return αp*(1+0.25*π*αs+0.25*π*ρ*(1+αs))/(1-0.25*π*αp*αs)

def Us_min(Cp,Cs):
    αp = Cp/(1-Cp)
    αs = Cs/(1-Cs)
    ρ = Cs/Cp
    return αs*(1+0.03125*π**3*αp+0.125*π**2*(1+αp)/ρ+0.185*αp*αs)/(1-0.25*π*αp*αs)


def run_example():
    Cp = 0.3
    Cs = 0.3
    
    print(f'Up_min = {Up_min(Cp,Cs):.3f}')
    print(f'Us_min = {Us_min(Cp,Cs):.3f}')
    
if __name__ == "__main__":
    run_example()