from math import sqrt

import numpy as np


def LSM(x, y, n):
    x_y = np.mean(x * y)
    x_ = np.mean(x)
    y_ = np.mean(y)

    x_2 = np.mean(x ** 2)
    y_2 = np.mean(y ** 2)

    b = (x_y - (x_ * y_))/\
        (x_2 - (x_ ** 2))
    sigma_b = 1/sqrt(n) * sqrt((y_2 - (y_ ** 2))/
                               (x_2 - (x_ ** 2)) - b ** 2)
    epsilon = sigma_b / b

    a = y_ - b * x_

    return b, a, sigma_b, epsilon


def mean_square_error(data):
    data_mean = np.mean(data)
    delta = data - data_mean
    sigma = sqrt(len(data) / (len(data) - 1) * sum(delta ** 2))
    return data_mean, sigma