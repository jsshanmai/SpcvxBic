
# coding: utf-8

# In[34]:


import numpy as np
from scipy import linalg


def elk_python(n,p):

    el1 = np.zeros([n,int(n*(n-1)/2)])
    el2 = np.zeros([n,int(n*(n-1)/2)])

    count = -1
    for i in list(range(n-1,0,-1)):
        temp = np.zeros([n,i])
        temp[n-i-1,:] = 1
        el1[:,(count+1):(count+i+1)] = temp
        el2[(n-i):n, (count + 1):(count + i + 1)] = np.eye(i)
        count = count+i

    ek1 = np.zeros([p,int(p*(p-1)/2)])
    ek2 = np.zeros([p,int(p*(p-1)/2)])

    count = -1
    for i in list(range(p-1,0,-1)):
        temp = np.zeros([p,i])
        temp[p-i-1,:] = 1
        ek1[:,(count+1):(count+i+1)] = temp
        ek2[(p-i):p, (count + 1):(count + i + 1)] = np.eye(i)
        count = count+i

    return(el1,el2,ek1,ek2)


def SCB_ADMM_python_WS(X, A, v, z, g, nu1, nu2, nu3, gamma_1, gamma_2, gamma_3,  w_l, u_k, feature_weight,  tol, niter=1,output = 1):

    n = X.shape[0]
    p = X.shape[1]

    n2 = int(n*(n-1)/2)
    p2 = int(p*(p-1)/2)

    elks = elk_python(n,p)
    el1 = elks[0]
    el2 = elks[1]
    ek1 = elks[2]
    ek2 = elks[3]

    En = np.diag(list(range(n))) + np.diag(list(range(n-1,-1,-1))) - np.ones([n,n]) + np.eye(n)
    Ep = np.diag(list(range(p))) + np.diag(list(range(p-1,-1,-1))) - np.ones([p,p]) + np.eye(p)

    M = np.eye(n) + nu1 * En
    N = nu2 * Ep + nu3 * np.eye(p)

    # 初始化
    
    if A is None : A = np.zeros(n*p).reshape(n,p)

    if v is None : v = np.zeros(p*n2).reshape(p,n2)
    
    if z is None : z = np.zeros(p2*n).reshape(n,p2)

    if g is None : g = np.zeros(n*p).reshape(n,p)
    

    lambda_1 = np.zeros(p*n2).reshape(p,n2)
    lambda_2 = np.zeros(p2*n).reshape(n,p2)
    lambda_3 = np.zeros(n*p).reshape(n,p)

    ## iterations
    for iters in range(int(niter)):

        A_old = A; v_old = v; z_old = z; g_old = g;
        lambda_1_old = lambda_1; lambda_2_old = lambda_2; lambda_3_old = lambda_3;

        # update A

        lv = lambda_1 + nu1 * v
        lz = lambda_2 + nu2 * z
        lg = lambda_3 + nu3 * g
        C2 = 0 -np.dot((el2-el1),lv.T)
        C3 = np.dot(lz,(ek1-ek2).T)
        C4 = np.dot(lg,np.eye(p))
        C = X +  C2 + C3 + C4

        A = linalg.solve_sylvester(M, N.T, C)

        al1 = np.dot(A.T,el1)
        al2 = np.dot(A.T,el2)
        ak1 = np.dot(A,ek1)
        ak2 = np.dot(A,ek2)

        # update v z g

        sigma_1 = gamma_1 * w_l/nu1
        sigma_1 = sigma_1.flatten()
        vtemp = al1 - al2 - 1/nu1 * lambda_1

        temp1 = np.where((1 - sigma_1/np.sqrt(np.sum(vtemp**2,axis=0)) < 0),0,1 - sigma_1/np.sqrt(np.sum(vtemp**2,axis=0)))
        temp2 = np.repeat(temp1,p).reshape(n2,p).T * vtemp

        v = temp2

        sigma_2 = gamma_2 * u_k/nu2
        sigma_2 = sigma_2.flatten()
        ztemp = ak1 - ak2 - 1/nu2 * lambda_2

        temp3 = np.where((1 - sigma_2/np.sqrt(np.sum(ztemp**2,axis=0)) < 0), 0 ,1 - sigma_2/np.sqrt(np.sum(ztemp**2,axis=0)))
        temp4 = np.repeat(temp3,n).reshape(p2,n).T * ztemp

        z = temp4
        
        sigma_3 = gamma_3 * feature_weight/nu3
        sigma_3 = sigma_3.flatten()
        gtemp = A - 1/nu3 * lambda_3
    
        temp5 = np.where((1 - sigma_3/np.sqrt(np.sum(gtemp**2,axis=0))) < 0, 0,1 - sigma_3/np.sqrt(np.sum(gtemp**2,axis=0)))
        temp6 = np.repeat(temp5,n).reshape(p,n).T * gtemp
    
        g = temp6

        # update lambda

        lambda_1 = lambda_1 + nu1 * (v - al1 + al2)
        lambda_2 = lambda_2 + nu2 * (z - ak1 + ak2)
        lambda_3 = lambda_3 + nu3 * (g - A)

        # update gamma 
        # gamma_1 = gamma_1 * t
        # gamma_2 = gamma_2 * t
        # gamma_3 = gamma_3 * t


        # output print

        if output == 1:
            print('iteration time: ' + str(iters))
            print('A_new - A_old: ' + str(round(np.mean(np.abs(A - A_old)),3)))
            print('v_new - v_old: ' + str(round(np.mean(np.abs(v - v_old)),3)))
            print('z_new - z_old: ' + str(round(np.mean(np.abs(z - z_old)),3)))
            print('g_new - g_old: ' + str(round(np.mean(np.abs(g - g_old)),3)))
            print('lambda1_new - lambda1_old: ' + str(round(np.mean(np.abs(lambda_1 - lambda_1_old)),3)))
            print('lambda2_new - lambda2_old: ' + str(round(np.mean(np.abs(lambda_2 - lambda_2_old)),3)))
            print('lambda3_new - lambda3_old: ' + str(round(np.mean(np.abs(lambda_3 - lambda_3_old)),3)))
            #print('gamma_3 =: ' + str(round(np.mean(np.abs(gamma_3)),3)))

        # return

        conditions = ((np.mean(np.abs(A - A_old)) < tol) &
                      (np.mean(np.abs(v - v_old)) < tol) &
                      (np.mean(np.abs(z - z_old)) < tol) &
                      (np.mean(np.abs(g - g_old)) < tol) &
                      (np.mean(np.abs(lambda_1 - lambda_1_old))< tol) &
                      (np.mean(np.abs(lambda_2 - lambda_2_old))< tol) &
                      (np.mean(np.abs(lambda_3 - lambda_3_old))< tol) )

        if conditions:
            return(A, v, z, g, lambda_1, lambda_2, lambda_3,gamma_1,gamma_2,gamma_3, iters)
            break


    return(A, v, z, g, lambda_1, lambda_2, lambda_3, gamma_1,gamma_2,gamma_3, iters)



