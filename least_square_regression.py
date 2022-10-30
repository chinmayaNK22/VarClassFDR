import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

##plt.style.use('seaborn-poster')
##
### generate x and y
##x = np.linspace(0, 1, 101)
##y = 1 + x + x * np.random.random(len(x))
##
### assemble matrix A
##A = np.vstack([x, np.ones(len(x))]).T
##
### turn y into a column vector
##y = y[:, np.newaxis]
##
### Direct least square regression
##alpha = np.dot((np.dot(np.linalg.inv(np.dot(A.T,A)),A.T)),y)
##print(alpha)
##
### plot the results
##plt.figure(figsize = (10,8))
##plt.plot(x, y, 'b.')
##plt.plot(x, alpha[0]*x + alpha[1], 'r')
##plt.xlabel('x')
##plt.ylabel('y')
##plt.show()


from sklearn.datasets import make_regression
from sklearn.linear_model import LinearRegression

X, y, coefficients = make_regression(
    n_samples=50,
    n_features=1,
    n_informative=1,
    n_targets=1,
    noise=5,
    coef=True,
    random_state=1
)

print (X, y, coefficients)
