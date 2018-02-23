import numpy as np
import openturns as ot

########################

field = np.loadtxt('H_added_symmetric_points.dat', skiprows=1, delimiter=" ")
theta = field[:, 0]
phi = field[:, 1]
g = field[:, 2]
porosity = field[:, 3]
Re = field[:, 4]
k11 = field[:, 5]
k22 = field[:, 6]
k33 = field[:, 7]


x = np.column_stack((Re, theta, phi, porosity))
X = ot.Sample(x)

z1 = ot.Sample(np.reshape(k11, (len(g), 1)))
z2 = ot.Sample(np.reshape(k22, (len(g), 1)))
z3 = ot.Sample(np.reshape(k33, (len(g), 1)))

dimension = 4
# basis = ot.ConstantBasisFactory(dimension).build()
# basis = ot.LinearBasisFactory(dimension).build()
basis = ot.QuadraticBasisFactory(dimension).build()
# covarianceModel = ot.SquaredExponential([1.,1.3,1.3,2.], [3.1])
covarianceModel = ot.MaternModel([1.2, 1., 1., 1.], 1.9)


myLeastSquares1 = ot.KrigingAlgorithm(X, z1, covarianceModel, basis)
myLeastSquares2 = ot.KrigingAlgorithm(X, z2, covarianceModel, basis)
myLeastSquares3 = ot.KrigingAlgorithm(X, z3, covarianceModel, basis)
# myLeastSquares = ot.QuadraticLeastSquares(X, z)
# myLeastSquares = ot.LinearLeastSquares(X, z)
myLeastSquares1.run()
myLeastSquares2.run()
myLeastSquares3.run()

result1 = myLeastSquares1.getResult()
result2 = myLeastSquares2.getResult()
result3 = myLeastSquares3.getResult()

metamodel1 = result1.getMetaModel()
metamodel2 = result2.getMetaModel()
metamodel3 = result3.getMetaModel()


Re = np.linspace(0, 100, 41)
th = np.linspace(0, 45, 23)
ph = np.linspace(0, 90, 46)
por = np.linspace(0.4, 0.8, 41)

rev, thv, phv, porv = np.meshgrid(Re, th, ph, por)

tot_len = len(Re)*len(th)*len(ph)*len(por)

Re_ = rev.reshape(tot_len,)
theta_ = thv.reshape(tot_len,)
phi_ = phv.reshape(tot_len,)
porosity_ = porv.reshape(tot_len,)


K11 = np.column_stack((Re_, theta_,
                       phi_, porosity_))



mmk11 = np.array(metamodel1(K11)).reshape(tot_len,)
mmk22 = np.array(metamodel2(K11)).reshape(tot_len,)
mmk33 = np.array(metamodel3(K11)).reshape(tot_len,)


np.savetxt("kriging_K.dat", np.column_stack((Re_, theta_, phi_, porosity_,
                                          1./mmk11, 1./mmk22, 1./mmk33)), delimiter=" ")



