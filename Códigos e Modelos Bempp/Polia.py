from IPython import get_ipython
get_ipython().magic('reset -sf')

import time

import bempp.api
import numpy as np
bempp.api.global_parameters.assembly.boundary_operator_assembly_type="dense"
bempp.api.global_parameters.quadrature.medium.double_order = 4
bempp.api.global_parameters.quadrature.far.double_order = 4


import Mesh
file = 'Polia1';
dirichlet_segments = [32,3,7]
mymesh = Mesh.Mesh()

mymesh.read_msh(file + '.msh')
    

XYZ = mymesh.Verts
tri = mymesh.Elmts[2][1]-1
surf = mymesh.Elmts[2][0]


domain_index=surf
segments=np.unique(surf)


def dirichlet_data_fun(x,domain_index):
    if(domain_index==32):
        return 220
    else:
        return 0

def neumann_data_fun(x):
    return 0

order_neumann = 0
order_dirichlet = 1

neumann_segments=[iseg for iseg in segments if iseg not in dirichlet_segments]
grid= bempp.api.grid_from_element_data(XYZ.transpose(), tri.transpose(),domain_index)

global_neumann_space = bempp.api.function_space(grid, "DP", order_neumann)
global_dirichlet_space = bempp.api.function_space(grid, "P", order_dirichlet)



dirichlet_space_neumann_segment = bempp.api.function_space(
    grid, "P", order_dirichlet, domains=neumann_segments, closed=False)

dual_dirichlet_space = bempp.api.function_space(
    grid, "P", order_dirichlet, domains=dirichlet_segments,
    closed=True, strictly_on_segment=True)

neumann_space_dirichlet_segment = bempp.api.function_space(
    grid, "DP", order_neumann, domains=dirichlet_segments,
    closed=True, element_on_segment=True)

neumann_space_neumann_segment = bempp.api.function_space(
    grid, "DP", order_neumann, domains=neumann_segments,
    closed=False, element_on_segment=True, reference_point_on_segment=False)

dirichlet_space_dirichlet_segment = bempp.api.function_space(
    grid, "P", order_dirichlet, domains=dirichlet_segments, closed=True)

slp_DD = bempp.api.operators.boundary.laplace.single_layer(
    neumann_space_dirichlet_segment,
    dirichlet_space_dirichlet_segment,
    neumann_space_dirichlet_segment)  # used in "blocked"

dlp_DN = bempp.api.operators.boundary.laplace.double_layer(
    dirichlet_space_neumann_segment,
    dirichlet_space_dirichlet_segment,
    neumann_space_dirichlet_segment)  # used in "blocked"

adlp_ND = bempp.api.operators.boundary.laplace.adjoint_double_layer(
    neumann_space_dirichlet_segment,
    neumann_space_neumann_segment,
    dirichlet_space_neumann_segment)  # used in "blocked"

hyp_NN = bempp.api.operators.boundary.laplace.hypersingular(
    dirichlet_space_neumann_segment,
    neumann_space_neumann_segment,
    dirichlet_space_neumann_segment)  # used in "blocked"

slp_DN = bempp.api.operators.boundary.laplace.single_layer(
    neumann_space_neumann_segment,
    dirichlet_space_dirichlet_segment,
    neumann_space_dirichlet_segment)  # not used in "blocked"

dlp_DD = bempp.api.operators.boundary.laplace.double_layer(
    dirichlet_space_dirichlet_segment,
    dirichlet_space_dirichlet_segment,
    neumann_space_dirichlet_segment)  # not used in "blocked"

id_DD = bempp.api.operators.boundary.sparse.identity(
    dirichlet_space_dirichlet_segment,
    dirichlet_space_dirichlet_segment,
    neumann_space_dirichlet_segment)  # not used in "blocked"

adlp_NN = bempp.api.operators.boundary.laplace.adjoint_double_layer(
    neumann_space_neumann_segment,
    neumann_space_neumann_segment,
    dirichlet_space_neumann_segment)  # not used in "blocked"

id_NN = bempp.api.operators.boundary.sparse.identity(
    neumann_space_neumann_segment,
    neumann_space_neumann_segment,
    dirichlet_space_neumann_segment)  # not used in "blocked"

hyp_ND = bempp.api.operators.boundary.laplace.hypersingular(
    dirichlet_space_dirichlet_segment,
    neumann_space_neumann_segment,
    dirichlet_space_neumann_segment)  # not used in "blocked"

blocked = bempp.api.BlockedOperator(2, 2)

blocked[0, 0] = slp_DD
blocked[0, 1] = -dlp_DN
blocked[1, 0] = adlp_ND
blocked[1, 1] = hyp_NN


    
def dirichlet_data(x, n, domain_index, res):
    res[0] = dirichlet_data_fun(x,domain_index)

 
def neumann_data(x, n, domain_index, res):
    res[0] = neumann_data_fun(x)

dirichlet_grid_fun = bempp.api.GridFunction(
    dirichlet_space_dirichlet_segment,
    fun=dirichlet_data, dual_space=dual_dirichlet_space)

neumann_grid_fun = bempp.api.GridFunction(
    neumann_space_neumann_segment,
    fun=neumann_data, dual_space=dirichlet_space_neumann_segment)
start_time = time.time()
rhs_fun1 = (.5 * id_DD + dlp_DD) * dirichlet_grid_fun \
           - slp_DN * neumann_grid_fun

rhs_fun2 = - hyp_ND * dirichlet_grid_fun \
           + (.5 * id_NN - adlp_NN) * neumann_grid_fun

lhs = blocked.weak_form()
rhs = np.hstack([rhs_fun1.projections(neumann_space_dirichlet_segment), 
                 rhs_fun2.projections(dirichlet_space_neumann_segment)])
elapsed_time1= time.time() - start_time
from scipy.sparse.linalg import gmres
start_time = time.time()

x, info = gmres(lhs, rhs)
elapsed_time2= time.time() - start_time

nx0 = neumann_space_dirichlet_segment.global_dof_count

neumann_solution = bempp.api.GridFunction(
    neumann_space_dirichlet_segment, coefficients=x[:nx0])
dirichlet_solution = bempp.api.GridFunction(
    dirichlet_space_neumann_segment, coefficients=x[nx0:])

# bempp.api.PLOT_BACKEND = "ipython_notebook"

neumann_imbedding_dirichlet_segment = \
    bempp.api.operators.boundary.sparse.identity(
        neumann_space_dirichlet_segment,
        global_neumann_space,
        global_neumann_space)

neumann_imbedding_neumann_segment = \
    bempp.api.operators.boundary.sparse.identity(
        neumann_space_neumann_segment,
        global_neumann_space,
        global_neumann_space)

dirichlet_imbedding_dirichlet_segment = \
    bempp.api.operators.boundary.sparse.identity(
        dirichlet_space_dirichlet_segment,
        global_dirichlet_space,
        global_dirichlet_space)

dirichlet_imbedding_neumann_segment = \
    bempp.api.operators.boundary.sparse.identity(
        dirichlet_space_neumann_segment,
        global_dirichlet_space,
        global_dirichlet_space)

dirichlet = (dirichlet_imbedding_dirichlet_segment * dirichlet_grid_fun +
             dirichlet_imbedding_neumann_segment * dirichlet_solution)

neumann = (neumann_imbedding_neumann_segment * neumann_grid_fun +
           neumann_imbedding_dirichlet_segment * neumann_solution)
dirichlet.plot()

unknown_dofs=x.shape[0]
elements = grid.leaf_view.entity_count(0)
number_of_global_neumann_dofs = global_neumann_space.global_dof_count
number_of_global_dirichlet_dofs = global_dirichlet_space.global_dof_count
total_time=elapsed_time1+elapsed_time2

print("time to assembly matrices blocks: "+str(elapsed_time1))
print("time to solve linear system: "+str(elapsed_time2))
print("total time: "+str(total_time))
print("Number of elements "+str(elements))
print("Number of global Neumann dofs: "+str(number_of_global_neumann_dofs))
print("Number of global Dirichlet dofs: "+str(number_of_global_dirichlet_dofs))
print("Number of unknown dofs: "+str(unknown_dofs))