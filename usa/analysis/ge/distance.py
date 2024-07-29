import numpy as np

def get_angle(x, y, z, sweep_radius, minimise_memory = False, verbose = False):
    """
    given a 3D-coordinate return the radian angle of to origin
    pBAngle : in b-mode plane
    pVAngle : in b-mode plane

    x,y,z : physical position of the image co-ordinate
    sweep_radius : resolution * offset

    NOTE: This method seems to crash on large files. For now, I have implemented a crude workaround using minimise_memory.
    """
    if not minimise_memory:

        rDsq = sweep_radius ** 2

        xsq = x ** 2
        ysq = y ** 2
        zsq = z ** 2

        tworDz = 2 * sweep_radius * z
        tmpone = rDsq + ysq + zsq + tworDz
        tmptwo = xsq + ysq + zsq + tworDz + 2 * rDsq
        tmpthree = np.sqrt(tmpone) * 2 * sweep_radius

        pRB = np.sqrt(tmptwo - tmpthree)
        pBAngle = np.arcsin(x / pRB)
        pVAngle = np.arcsin(y / (sweep_radius + pRB * np.cos(pBAngle)))
    
    else: # minimise memory
        sh = np.shape(x)
        pRB = np.zeros(sh)

        for i in range(sh[0]):
            if verbose:
                print(str(int(100*i/sh[0])),end='\r')
            # for j in range(sh[1]):
            #     for k in range(sh[2]):

            #         # if not self.seg_array[k,i,j] == self.vol_id: # note transpose
            #         #     pRB[i,j,k] = 0
            #         #     continue

            #         rDsq = sweep_radius ** 2

            #         xsq = x[i,j,k] ** 2
            #         ysq = y[i,j,k] ** 2
            #         zsq = z[i,j,k] ** 2

            #         tworDz = 2 * sweep_radius * z[i,j,k]
            #         tmpone = rDsq + ysq + zsq + tworDz
            #         tmptwo = xsq + ysq + zsq + tworDz + 2 * rDsq
            #         tmpthree = np.sqrt(tmpone) * 2 * sweep_radius

            #         pRB[i,j,k] = np.sqrt(tmptwo - tmpthree)
                    # if not self.seg_array[k,i,j] == self.vol_id: # note transpose
                    #     pRB[i,j,k] = 0
                    #     continue

            rDsq = sweep_radius ** 2

            xsq = x[i,:,:] ** 2
            ysq = y[i,:,:] ** 2
            zsq = z[i,:,:] ** 2

            tworDz = 2 * sweep_radius * z[i,:,:]
            tmpone = rDsq + ysq + zsq + tworDz
            tmptwo = xsq + ysq + zsq + tworDz + 2 * rDsq
            tmpthree = np.sqrt(tmpone) * 2 * sweep_radius

            pRB[i,:,:] = np.sqrt(tmptwo - tmpthree)
            # print(str(int(100*i/sh[0])), end='\r')
        pBAngle = 0
        pVAngle = 0

    # self.vrbs("There are some non.inf's: "+str(np.sum(pRB<np.inf))) # can't remember what this is about?

    return pBAngle, pVAngle, pRB