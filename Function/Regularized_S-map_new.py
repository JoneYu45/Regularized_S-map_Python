#Import modules
import numpy as np
from sklearn.linear_model import ElasticNetCV
import pandas as pd
import time
from joblib import Parallel, delayed
import multiprocessing
from optparse import OptionParser
import os

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

def select_frequent_dominate_genera(input, dominate_threshold, zero_frequency_threshold, select):
    # Process data
    rawdata = pd.read_csv(input).iloc[:, 1:]
    abundance = rawdata

    # Calculate the relative abundance profile
    read_num = np.sum(rawdata, axis=1)
    for i in range(rawdata.shape[0]):
        abundance.iloc[i, :] = rawdata.iloc[i, :] / read_num[i] * 100

    # Process or not
    if select == True:
        # Select the most frequent and dominate genera
        dominate = np.where(np.mean(abundance, axis=0) > dominate_threshold)
        wanted_abundance = abundance.iloc[:, dominate[0]]
        frequency = np.where(
            (wanted_abundance == 0).astype(int).sum(axis=0) / abundance.shape[0] * 100 < zero_frequency_threshold)
        wanted_abundance = wanted_abundance.iloc[:, frequency[0]]
    else:
        wanted_abundance = abundance

    #Output selection
    return wanted_abundance

def Regularized_Smap(abund, target_otu, theta, l_grid, iteration, cv, train_len, uncontinuous):
    print('Process data for otu No. %s' % str(target_otu+1))
    # Make input for the elastic_net
    block = np.append(abund[1:, target_otu], abund[0:-1, ], axis=1)
    ##Delete the uncontinuous states
    ##Commonly, we inferer the Jacobian matrices using the continuous time series. However, if you don't have enough time points but have the replicate time series from independet reactors.
    ##You can combine the replicate OTU tables as the input but delete uncontinuous states in the block.
    if uncontinuous == True:
        block = np.delete(block, [abund.shape[0] / 3 - 1, abund.shape[0] / 3 * 2 - 1], axis=0) 
        ##Triplicate time series are used as example, so we remove two uncontiunous states, i.e., [abund.shape[0] / 3 - 1, abund.shape[0] / 3 * 2 - 1], in the block. 
        ##You can also specify the list of uncontiunous states in the block using the following line. Remember to uncomment the the following line and delete line 55.
        # block = np.delete(block, [uncontiunous states], axis=0)
        ##The [uncontiunous states] can be like [34, 55] or [1, 2, 10] as you need it to be.
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
    coefs.to_csv('/'.join([output_dir,'coefs/%s_%s_coefs.csv' % (target_otu, theta)]))
    fit_results = pd.DataFrame(
        columns=['Intercept', 'Best alpha', 'Best l1_ratio', 'RMSE', 'Std', 'RMSE_o', 'Std_o', 'Test set score',
                 'Test set score_train', 'ymax_train', 'ymin_train', 'ymax_test', 'ymin_test'],
        data=fit_results)
    fit_results.to_csv('/'.join([output_dir,'fit_result/%s_%s_fit_results.csv' % (target_otu, theta)]))


if __name__ == '__main__':
    # Imnput data and setting
    parse = OptionParser()
    parse.add_option('-I', '--input', dest='input', default='../inputs/month_sample.csv')
    parse.add_option('-O', '--output', dest='output', default='../outputs')
    parse.add_option('-C', '--CV', dest='CV', default=10)
    parse.add_option('-T', '--train-length', dest='tl', default=5)

    (options, args) = parse.parse_args()
    input = options.input
    output_dir = options.output
    l_grid = 0.01
    iteration = 1000000
    cv = options.CV
    train_len = int(options.tl)
    uncontinuous = False
    dominate_threshold = 1  # OTUs whose abundances are more than the threshold will be included in the inference.
    zero_frequency_threshold = 20  # OTUs whose absence frequency are less than threshold will be included in the inference.
    select = True # If true, we will select the OTUs consistent with the aforementioned resitriction. Otherwise, all OTUs will be used for the inference.

    print('Start!')
    # Make ouput direction
    path_coefs = '/'.join([output_dir, 'coefs'])
    if not os.path.exists(path_coefs):
        os.makedirs('/'.join([output_dir, 'coefs']))
    path_fitresult = '/'.join([output_dir, 'fit_result'])
    if not os.path.exists(path_fitresult):
        os.makedirs('/'.join([output_dir, 'fit_result']))

    # Imnput data and parameter setting
    ts = select_frequent_dominate_genera(input, dominate_threshold, zero_frequency_threshold, select)
    print('Output analyzed OTUs')
    ts.to_csv('/'.join([output_dir, 'abund.csv']))
    abund = np.mat(ts)

    #Work in parallel
    num_cores = multiprocessing.cpu_count()
    for target_otu in range(abund.shape[1]):
        Parallel(n_jobs=num_cores, backend='multiprocessing')(
            delayed(Regularized_Smap)(abund, target_otu, theta, l_grid, iteration, cv, train_len, output_dir) for theta in
            [0.1, 0.5, 1, 2, 5, 10]) #You can specify the thetas (θ) list you want to try. Only the states closer to the target state will be used for regression when θ is large.
    print('\nFinished!')
