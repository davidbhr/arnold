def linear_stiffness(stiffness):
    """
    Defines the material properties for a linear elastic material.
    
    Args:
        stiffness(float) : Bulk modulus of the material in Pascal for SAENO Simulation (see [Steinwachs,2015])
   
    """
    return {'K_0': stiffness, 'D_0': 1e30, 'L_S': 1e30, 'D_S': 1e30}

    

def youngs_modulus(youngs_modulus):
    """
    Defines the material properties for a linear-elastic hydrogel (such as Matrigel) hydrogel with a poission ratio of 0.25 (see [Steinwachs,2015], in this case Youngsmodulus equals K_0/6).
    
    Args:
        youngs_modulus(float) : Young's modulus of the material in Pascal for SAENO Simulation (see [Steinwachs,2015])
    """
    return {'K_0': youngs_modulus*6, 'D_0': 1e30, 'L_S': 1e30, 'D_S': 1e30}


def youngs_modulus_saenopy(youngs_modulus):
    """
    Defines the material properties for a linear-elastic hydrogel (such as Matrigel) hydrogel with a poission ratio of 0.25 (see [Steinwachs,2015], in this case Youngsmodulus equals K_0/6).
    
    Use None in saenopy
    
    Args:
        youngs_modulus(float) : Young's modulus of the material in Pascal for SAENO Simulation (see [Steinwachs,2015])
    """
    return {'K_0': youngs_modulus*6, 'D_0': None, 'L_S': None, 'D_S': None}


def custom(K_0, D_0, L_S, D_S):
    """
    Defines the material properties for a custom nonlinear material.
    
    Args:
        K_0(float) : Bulk modulus of the material in Pascal for SAENO Simulation (see [Steinwachs,2015])
        D_0(float) : Buckling coefficient of the fibers for SAENO Simulation (see [Steinwachs,2015])
        L_S(float) : Onset of strain stiffening of the fibers (see [Steinwachs,2015])
        D_S(float) : Strain stiffening coefficient of the fibers (see [Steinwachs,2015])   
    """
    return {'K_0': K_0, 'D_0': D_0, 'L_S': L_S, 'D_S': D_S}


collagen_06 = {'K_0': 447, 'D_0': 0.0008, 'L_S': 0.0075, 'D_S': 0.033}
collagen_12 = {'K_0': 1645, 'D_0': 0.0008, 'L_S': 0.0075, 'D_S': 0.033}
collagen_24 = {'K_0': 5208, 'D_0': 0.0008, 'L_S': 0.0075, 'D_S': 0.033}
matrigel_10mg_ml = {'K_0': 200*6, 'D_0': 1e30, 'L_S': 1e30, 'D_S': 1e30}
