#Import modules
import numpy as np
from sklearn.linear_model import ElasticNetCV
import pandas as pd
import time
from joblib import Parallel, delayed
import multiprocessing

#Define functions
def make_weights(E_dist, theta):
    w = [np.exp(-1 * theta * E_dist[i] / np.mean(E_dist)) for i in range(E_dist.shape[0])]
    return w

def weight_data(sub_block, w):
    wp = np.empty(shape=sub_block.shape)
    for i in range(sub_block.shape[0]):
        for j in range(sub_block.shape[1]):
            wp[i, j] = sub_block[i, j] * w[i]
    return wp

def Regularized_Smap(abund, target_otu, theta, l_grid, iteration, cv, train_len):
    print('Process data for otu No. %s' % str(target_otu+1))
    # Make input for the elastic_net
    block = np.append(abund[1:, target_otu], abund[0:-1, ], axis=1)
    ##Delete the uncontinuous states
    block = np.delete(block, [abund.shape[0] / 3 - 1, abund.shape[0] / 3 * 2 - 1], axis=0)
    ##Scaling the input
    ##Each time series is normalized to have a mean of 0 and standard deviation of 1 before analysis with S-maps
    block = (block - np.average(block, axis=0)) / np.std(block, axis=0)

    ##Select data and fitting
    print('Start fitting.')
    lib = range(block.shape[0])
    coefs = np.empty(shape=(block.shape[0], block.shape[1] - 1))
    fit_results = np.empty(shape=(block.shape[0], 13))

    for ipred in lib:
        print('\r', 'Complete percentage: %.2f%%' % (ipred / len(lib) * 100), end="", flush=True)
        sub_block = np.delete(block, ipred, axis=0)
        q = block[lib[ipred], :]
        ###Calculate weights
        E_dist = np.sqrt(np.sum(np.array(sub_block[:, 1:] - q[:, 1:]) ** 2, axis=1))
        w = make_weights(E_dist, theta)
        ###Weighted predictors and responses
        X_wp = weight_data(sub_block[:, 1:], w)
        Y_wp = np.ravel(weight_data(sub_block[:, 0], w))
        X_target = block[ipred, 1:]
        Y_target = block[ipred, 0]

        ##Split training and test data
        pick_test = np.random.choice(range(X_wp.shape[0]), size=train_len, replace=False)
        X_train = np.append(np.delete(X_wp, pick_test, axis=0), X_target, axis=0)
        X_test = X_wp[pick_test, :]
        Y_train = np.append(np.delete(Y_wp, pick_test, axis=0), Y_target)
        Y_test = Y_wp[pick_test]

        ###Fit function
        regr = ElasticNetCV(cv=cv, random_state=0, max_iter=iteration,
                            l1_ratio=[(i + 1) * l_grid for i in range(int(1 / l_grid))])
        regr.fit(X_train, Y_train)
        rmse = np.sqrt(np.mean((regr.predict(X_train) - Y_train) ** 2))
        rmse_o = np.sqrt(np.mean((regr.predict(X_test) - Y_test) ** 2))
        coefs[ipred, :] = regr.coef_
        fit_results[ipred, :] = regr.intercept_, regr.alpha_, regr.l1_ratio_, rmse, np.std(Y_train), rmse_o, np.std(
            Y_test), regr.score(X_test, Y_test), regr.score(X_train, Y_train), max(Y_train), min(Y_train), max(
            Y_test), min(Y_test)
        print('\r', 'Complete percentage: %.2f%%' % ((ipred + 1) / len(lib) * 100), end="", flush=True)

        # Output results
    coefs = pd.DataFrame(data=coefs)
    coefs.to_csv('../Output/test/0/coefs/%s_%s_coefs.csv' % (target_otu, theta))
    fit_results = pd.DataFrame(
        columns=['Intercept', 'Best alpha', 'Best l1_ratio', 'RMSE', 'Std', 'RMSE_o', 'Std_o', 'Test set score',
                 'Test set score_train', 'ymax_train', 'ymin_train', 'ymax_test', 'ymin_test'],
        data=fit_results)
    fit_results.to_csv('../Output/test/0/fit_result/%s_%s_fit_results.csv' % (target_otu, theta))


if __name__ == '__main__':
    print('Start!')
    #Imnput data and parameter setting
    ts = pd.read_csv('../Input/intersect_OTU_0.csv',index_col=0) ##Triplicate data
    abund = np.mat(ts)
    l_grid = 0.01
    iteration = 1000000
    cv = 10
    train_len = 5
    t_grid = 0.01

    #Work in parallel
    num_cores = multiprocessing.cpu_count()
    for target_otu in range(abund.shape[1]):
        Parallel(n_jobs=num_cores, backend='multiprocessing')(
            delayed(Regularized_Smap)(abund, target_otu, theta, l_grid, iteration, cv, train_len) for theta in
            [0.1, 0.5, 1, 2, 5, 10])
    print('\nFinished!')