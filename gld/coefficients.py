"""
The different coefficients defined in the paper by Scrymgeour et al.


Theta

# material-constants:
Epsilon11, Epsilon22

Q = ??? ??? ??? ???
    ??? ??? ??? ???
    Q31 ??? Q33 ???
    ??? Q42 ??? Q44

C = C11 C12 C13 C14
    ??? ??? ??? ???
    ??? C33 ??? ???
    ??? ??? ??? C44

"""

import numpy as np
from scipy.constants import epsilon_0
########### coefficients ###########



def coefs(Theta, Epsilon11, Epsilon33, Q, C, Ph, print_pars=False, return_all=False):
    '''
    Handles the coefficient-jungle.
    Needs: Angle, Epsilon11, Epsilon33, Q (4x4-array), C (4x4-array), Ph
    Returns: S1, S2, Xi, Phi
    '''

    #first layer:
    Alpha1 = 1/(2*Epsilon33*epsilon_0)
    Alpha2 = 0
    Alpha3 = 1/(Epsilon11*epsilon_0)

    Beta1 = (0.5 * C[2, 2]) # does paper use abs() ???
    Beta2 = 0.25 * (C[0, 0] + C[0, 1])
    Beta3 = 0.25 * (C[0, 0] - C[0, 1])
    Beta4 = C[0, 2]
    Beta5 = 0.5 * C[3, 3]
    Beta6 = C[0, 3]

    Gamma1 = 0.5 * (C[0, 0] + C[0, 1])*Q[2, 0] + 0.5 * C[0, 2]*Q[2, 2]
    Gamma2 = 0.5 * C[2, 2]*Q[2, 2] + 0.5 * C[0, 2]*Q[2, 0]
    Gamma3 = 2 * C[0, 3]*Q[3, 3] - 0.5 * (C[0, 0] - C[0, 1])*Q[3, 1]
    Gamma4 = C[3, 3] * Q[3, 3]

    # second layer:
    Psi1 = ( 2*Gamma1*Beta1 - Gamma2*Beta4 ) / ( 2*( Beta4**2 - 4*Beta1*Beta2 ) )
    Psi2 = ( 2*Gamma2*Beta2 - Gamma1*Beta4 ) / ( Beta4**2 - 4*Beta1*Beta2 )

    #Ph = ( Alpha1 / ( Alpha2 + 4*( Beta1*(Psi2**2) + 4*Beta2*(Psi1**2) \
    #                            + 2*Beta4*Psi1*Psi2 + 2*Gamma1*Psi1 + Gamma2*Psi2 ) ) )**0.5
    Alpha2 = Alpha1/(Ph**2) \
        - 4*( Beta1*(Psi2**2) + 4*Beta2*(Psi1**2) + 2*Beta4*Psi1*Psi2 + 2*Gamma1*Psi1 + Gamma2*Psi2 )


    Lambda1 = Psi1 * (Ph**2)
    Lambda2 = Psi2 * (Ph**2)
    


    # third layer:
    Aij = np.array([
        [ 2*(Beta2+Beta3), Beta6*np.sin(3*Theta), 0 ],
        [ Beta6*np.sin(3*Theta), 2*Beta5, Beta6*np.cos(3*Theta) ],
        [ 0, Beta6*np.cos(3*Theta), 2*Beta3 ]
    ])

    Bij = np.array([
        [ -Gamma1, -Gamma3*np.sin(3*Theta), -Gamma3*np.cos(3*Theta) ],
        [ 0, -Gamma4, 0 ], 
        [ 0, -Gamma3*np.cos(3*Theta), Gamma3*np.sin(3*Theta) ]
    ])

    Mij = np.matmul( np.linalg.inv(Aij), Bij )

    Nu1 = -Gamma3*Mij[2, 0]*np.cos(3*Theta) - Gamma3*Mij[0, 0]*np.sin(3*Theta) - Gamma4*Mij[1, 0]
    Nu2 = -Gamma3*Mij[0, 0]*np.cos(3*Theta) + Gamma3*Mij[2, 0]*np.sin(3*Theta)

    Mu11 = Gamma3*Mij[2, 1]*np.cos(3*Theta) + Gamma3*Mij[0, 1]*np.sin(3*Theta) + Gamma4*Mij[1, 1] 
    Mu12 = Gamma3*Mij[2, 2]*np.cos(3*Theta) + Gamma3*Mij[0, 2]*np.sin(3*Theta) + Gamma4*Mij[1, 2] 

    Mu21 = Gamma3*Mij[0, 1]*np.cos(3*Theta) - Gamma3*Mij[2, 1]*np.sin(3*Theta)
    Mu22 = Gamma3*Mij[0, 2]*np.cos(3*Theta) - Gamma3*Mij[2, 2]*np.sin(3*Theta)

    Xi = np.array([
        [ -Nu1*(Ph**2)/Alpha3 , Nu1/Alpha3 - (Nu1*Mu22-Nu2*Mu12)*(Ph**2)/(Alpha3**2) , (Nu1*Mu22 - Nu2*Mu12)/(Alpha3**2)],
        [ -Nu2*(Ph**2)/Alpha3 , Nu2/Alpha3 - (Nu1*Mu21-Nu2*Mu11)*(Ph**2)/(Alpha3**2) , (Nu2*Mu11 - Nu1*Mu21)/(Alpha3**2)],
        [1, 0, 0]
    ])

    Phi = np.zeros( (3, 4) )
    Phi[:, 0] = - Mij[:, 0] * (Ph**2)
    Phi[:, 1] = Mij[:, 0] + Mij[:, 1]*Xi[0, 0] + Mij[:, 2]*Xi[1, 0]
    Phi[:, 2] = Mij[:, 1]*Xi[0, 1] + Mij[:, 2]*Xi[1, 1]
    Phi[:, 3] = Mij[:, 1]*Xi[0, 2] + Mij[:, 2]*Xi[1, 2]


    # fourth layer
    S1 = Alpha1 - 4*Gamma1*Lambda1 - 2*Gamma1*Phi[0, 0] - 2*Gamma2*Lambda2 \
        - Gamma3*( Phi[0, 0]*Xi[1, 0] + Phi[2, 0]*Xi[0, 0] )*np.cos(3*Theta) \
        - Gamma3*( Phi[0, 0]*Xi[0, 0] - Phi[2, 0]*Xi[1, 0] )*np.sin(3*Theta) \
        - Gamma4*Phi[1, 0]*Xi[0, 0]
    
    S3 = Alpha2 + 2*Gamma1*Phi[0, 1] \
        + Gamma3*( Phi[0, 0]*Xi[1, 1] + Phi[0, 1]*Xi[1, 0] + Phi[2, 0]*Xi[0, 1] + Phi[2, 1]*Xi[0, 0] ) * np.cos(3*Theta)\
        + Gamma3*( Phi[0, 0]*Xi[0, 1] + Phi[0, 1]*Xi[0, 0] - Phi[2, 0]*Xi[1, 1] - Phi[2, 1]*Xi[1, 0] ) * np.sin(3*Theta)\
        + Gamma4*(Phi[1, 0]*Xi[0, 1] + Phi[1, 1]*Xi[0, 0])
    
    if print_pars:
        dummy = [Alpha1, Alpha2, Alpha3, Beta1, Beta2, Beta3, Beta4, Beta5, Beta6, Gamma1, Gamma2, Gamma3, Gamma4, Lambda1, Lambda2]
        for dum in dummy:
            print(f'{dum: .3E}')

    # can set flag, to return the paramters needed for calculation of full energy as well
    if return_all:
        return S1, S3, Xi, Phi, [Alpha1, Alpha2, Alpha3], [Beta1, Beta2, Beta3, Beta4, Beta5, Beta6], [Gamma1, Gamma2, Gamma3, Gamma4], [Lambda1, Lambda2]
    
    return S1, S3, Xi, Phi



